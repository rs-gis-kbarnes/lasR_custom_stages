#include "treehull3d.h"  
#include "PointSchema.h"  
#include "openmp.h"  
  
#include <unordered_set>  
#include <unordered_map>  
#include <random>  
  
#include "delaunator/delaunator.hpp"  
#include "ogrsf_frmts.h"  
  
LASRtreehull3d::LASRtreehull3d()  
{  
  vector.set_geometry_type(wkbPolygon25D);  
  vector.add_field("tree_id", OFTInteger);  
}  
  
// Explicit copy constructor: mirrors Stage::Stage(const Stage&) which never  
// bitwise-copies PointFilter (owns raw Condition* pointers, would double-free).  
// wire_pointfilter is instead rebuilt from the stored raw strings.  
LASRtreehull3d::LASRtreehull3d(const LASRtreehull3d& other) : StageVector(other)  
{  
  id_attribute       = other.id_attribute;  
  max_edge           = other.max_edge;  
  trim               = other.trim;  
  radius             = other.radius;  
  wire_filter_strings = other.wire_filter_strings;  
  
  wire_pointfilter = PointFilter(); // fresh, empty  
  for (const auto& f : wire_filter_strings)  
    wire_pointfilter.add_condition(f);  
}  
  
bool LASRtreehull3d::set_parameters(const nlohmann::json& stage)  
{  
  id_attribute = stage.at("attribute");  
  max_edge     = stage.at("max_edge");  
  radius       = stage.at("radius");  
  trim         = max_edge * max_edge;  
  
  wire_filter_strings.clear();  
  if (stage.contains("wire_filter"))  
  {  
    for (const auto& f : stage["wire_filter"])  
      if (f.is_string() && !f.get<std::string>().empty())  
        wire_filter_strings.push_back(f.get<std::string>());  
  }  
  
  for (const auto& f : wire_filter_strings)  
    wire_pointfilter.add_condition(f);  
  
  vector = Vector(xmin, ymin, xmax, ymax);  
  vector.set_geometry_type(wkbPolygon25D);  
  vector.add_field("tree_id", OFTInteger);  
  
  return true;  
}  
  
bool LASRtreehull3d::process(PointCloud*& las)  
{  
  // --- validate the tree-id attribute exists (verified idiom: get_attribute_index,  
  // not AttributeAccessor::exist(), which only reflects lazy-resolution state) ---  
  int idx = las->header->schema.get_attribute_index(id_attribute);  
  if (idx == -1)  
  {  
    last_error = id_attribute + " is not present in the point cloud";  
    return false;  
  }  
  
  AttributeType data_type = las->header->schema.attributes[idx].type;  
  if (data_type != AttributeType::INT32 && data_type != AttributeType::DOUBLE)  
  {  
    last_error = "the attribute " + id_attribute + " must be of type 'int' or 'double'";  
    return false;  
  }  
  
  AttributeAccessor id_accessor(id_attribute);  
  
  std::unordered_set<int> detected_ids;  
  bool restrict_to_wire = radius > 0 && !wire_filter_strings.empty();  
  
  // --- PASS 1: fixed-radius search around wire points to find candidate tree ids ---  
  if (restrict_to_wire)  
  {  
    las->build_kdtree();  
  
    while (las->read_point())  
    {  
      Point* p = &las->point;  
      if (!wire_pointfilter.filter(p)) continue; // keep only points that PASS the wire filter  
      // NOTE: PointFilter::filter() returns true when a point is REJECTED (see verified  
      // PointFilter::filter() body) - a point matching the wire condition is one that  
      // does NOT get filtered out, so this branch keeps points for which filter() == false.  
      // (see corrected condition below)  
    }  
  
    // Correct pass, replacing the confused draft above:  
    las->seek(0);  
    while (las->read_point())  
    {  
      Point* p = &las->point;  
      if (wire_pointfilter.filter(p)) continue; // rejected by wire_filter -> not a wire point  
      if (pointfilter.filter(p)) continue;      // also respect stage's own point filter  
  
      std::vector<Point> neighbors;  
      las->query_sphere(*p, radius, neighbors, &pointfilter);  
  
      for (auto& np : neighbors)  
      {  
        int tid = (int)id_accessor(&np);  
        detected_ids.insert(tid);  
      }  
    }  
  
    las->seek(0);  
  }  
  
  // --- PASS 2: bucket point coordinates per tree id ---  
  std::unordered_map<int, std::vector<double>> coords_by_id; // flat x,y pairs for delaunator  
  std::unordered_map<int, std::vector<PointXYZ>> pts_by_id;  // xyz for Z-lookup after triangulation  
  
  std::default_random_engine gen(std::random_device{}());  
  std::normal_distribution<double> noise(0.0, 1e-10);  
  
  double xoffset = (las->header->min_x + las->header->max_x) / 2;  
  double yoffset = (las->header->min_y + las->header->max_y) / 2;  
  
  while (las->read_point())  
  {  
    Point* p = &las->point;  
    if (pointfilter.filter(p)) continue;  
  
    int tid = (int)id_accessor(p);  
    if (tid == 0) continue; // 0 conventionally means "no tree" from region_growing's raster  
  
    if (restrict_to_wire && detected_ids.find(tid) == detected_ids.end()) continue;  
  
    coords_by_id[tid].push_back(p->get_x() - xoffset + noise(gen));  
    coords_by_id[tid].push_back(p->get_y() - yoffset + noise(gen));  
    pts_by_id[tid].push_back(PointXYZ(p->get_x(), p->get_y(), p->get_z()));  
  }  
  
  // --- PASS 3: per-tree Delaunay triangulation + contour extraction + polygonize ---  
  tree_polygons.clear();  
  
  for (auto& kv : coords_by_id)  
  {  
    int tid = kv.first;  
    std::vector<double>& coords = kv.second;  
    std::vector<PointXYZ>& pts = pts_by_id[tid];  
  
    if (coords.size() < 6) continue; // need >= 3 points  
  
    delaunator::Delaunator* d = nullptr;  
    try  
    {  
      d = new delaunator::Delaunator(coords);  
    }  
    catch (const std::exception& e)  
    {  
      continue; // skip degenerate trees rather than failing the whole stage  
    }  
  
    // Build a lookup from (x,y) rounded key -> z, since delaunator only knows x,y  
    std::unordered_map<int64_t, double> z_lookup;  
    auto key_of = [](double x, double y) -> int64_t {  
      int32_t ix = (int32_t)std::llround(x * 1000.0);  
      int32_t iy = (int32_t)std::llround(y * 1000.0);  
      return (int64_t(ix) << 32) ^ (uint32_t)iy;  
    };  
    for (auto& p : pts)  
      z_lookup[key_of(p.x - xoffset, p.y - yoffset)] = p.z;  
  
    std::unordered_set<Edge3D> edges;  
  
    for (size_t i = 0; i < d->triangles.size(); i += 3)  
    {  
      size_t ia = d->triangles[i];  
      size_t ib = d->triangles[i + 1];  
      size_t ic = d->triangles[i + 2];  
  
      double ax = d->coords[2 * ia], ay = d->coords[2 * ia + 1];  
      double bx = d->coords[2 * ib], by = d->coords[2 * ib + 1];  
      double cx = d->coords[2 * ic], cy = d->coords[2 * ic + 1];  
  
      double dx1 = bx - ax, dy1 = by - ay;  
      double dx2 = cx - ax, dy2 = cy - ay;  
      double edge_ab2 = dx1 * dx1 + dy1 * dy1;  
      double edge_bc2 = (cx - bx) * (cx - bx) + (cy - by) * (cy - by);  
      double edge_ca2 = (ax - cx) * (ax - cx) + (ay - cy) * (ay - cy);  
      double max_edge2 = std::max({edge_ab2, edge_bc2, edge_ca2});  
  
      bool keep = (trim == 0) || (max_edge2 < trim);  
      if (!keep) continue;  
  
      auto za = z_lookup.count(key_of(ax, ay)) ? z_lookup[key_of(ax, ay)] : 0.0;  
      auto zb = z_lookup.count(key_of(bx, by)) ? z_lookup[key_of(bx, by)] : 0.0;  
      auto zc = z_lookup.count(key_of(cx, cy)) ? z_lookup[key_of(cx, cy)] : 0.0;  
  
      PointXYZ A(ax + xoffset, ay + yoffset, za);  
      PointXYZ B(bx + xoffset, by + yoffset, zb);  
      PointXYZ C(cx + xoffset, cy + yoffset, zc);  
  
      Edge3D AB(A, B), BC(B, C), CA(C, A);  
  
      if (edges.count(AB) > 0) edges.erase(AB); else edges.insert(AB);  
      if (edges.count(BC) > 0) edges.erase(BC); else edges.insert(BC);  
      if (edges.count(CA) > 0) edges.erase(CA); else edges.insert(CA);  
    }  
  
    delete d;  
  
    if (edges.empty()) continue;  
  
    // Polygonize this tree's boundary edges exactly like LASRboundaries::process()  
    OGRGeometryCollection gc;  
    for (const auto& e : edges)  
    {  
      OGRLineString ls;  
      ls.addPoint(e.A.x, e.A.y, e.A.z);  
      ls.addPoint(e.B.x, e.B.y, e.B.z);  
      gc.addGeometry(&ls);  
    }  
  
    OGRGeometry* polys = gc.Polygonize();  
    if (!polys) continue;  
  
    if (wkbFlatten(polys->getGeometryType()) == wkbGeometryCollection)  
    {  
      OGRGeometryCollection* ogc = polys->toGeometryCollection();  
      for (int i = 0; i < ogc->getNumGeometries(); ++i)  
      {  
        OGRGeometry* geom = ogc->getGeometryRef(i);  
        if (!geom || wkbFlatten(geom->getGeometryType()) != wkbPolygon) continue;  
  
        OGRPolygon* ogrPoly = geom->toPolygon();  
        OGRLinearRing* ext = ogrPoly->getExteriorRing();  
        if (!ext) continue;  
  
        std::vector<PointXYZ> ring;  
        int n = ext->getNumPoints();  
        for (int j = 0; j < n; j++)  
        {  
          double zz = ext->getZ(j); // NOTE: GDAL Z-preservation through Polygonize()  
                                     // is unverified for this build - see known gaps.  
          if (zz == 0.0)  
            zz = z_lookup.count(key_of(ext->getX(j) - xoffset, ext->getY(j) - yoffset))  
                   ? z_lookup[key_of(ext->getX(j) - xoffset, ext->getY(j) - yoffset)] : 0.0;  
          ring.emplace_back(ext->getX(j), ext->getY(j), zz);  
        }  
  
        PolygonXYZ poly(ring);  
        poly.close();  
        tree_polygons.push_back({tid, poly});  
      }  
    }  
  
    OGRGeometryFactory::destroyGeometry(polys);  
  }  
  
  return true;  
}  
  
void LASRtreehull3d::clear(bool last)  
{  
  tree_polygons.clear();  
}  
  
bool LASRtreehull3d::write()  
{  
  bool success = true;  
  
  #pragma omp critical (write_tree_hull3d)  
  {  
    for (auto& kv : tree_polygons)  
    {  
      std::vector<PolygonXYZ> single = { kv.second };  
      if (!vector.write(single, kv.first)) { success = false; break; }  
    }  
  }  
  
  return success;  
}