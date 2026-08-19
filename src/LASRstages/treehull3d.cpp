// src/LASRstages/treehull3d.cpp  
#include "treehull3d.h"  
#include "PointSchema.h"  
#include "openmp.h"  

#include <unordered_set>  
#include <unordered_map>  
#include <random>  
#include <algorithm>  

#include "delaunator/delaunator.hpp"  
#include "ogrsf_frmts.h"  

// ---- helper: quantized (x,y) key for the Z-lookup fallback map ----  
static inline std::pair<long, long> key_of(double x, double y, double eps = 1e-6)
{
	return { std::llround(x / eps), std::llround(y / eps) };
}
struct PairHash
{
	std::size_t operator()(const std::pair<long, long>& p) const
	{
		return std::hash<long>{}(p.first) ^ (std::hash<long>{}(p.second) << 1);
	}
};

// -------------------- constructor (real definition, NOT '= default' in header) --------------------  
LASRtreehull3d::LASRtreehull3d()
{
	vector.set_geometry_type(wkbPolygon25D);
	vector.add_field("tree_id", OFTInteger);

	// NEW: mesh output geometry/schema  
	mesh_vector.set_geometry_type(wkbMultiPolygon25D);
	mesh_vector.add_field("tree_id", OFTInteger);
}

// -------------------- copy constructor: PointFilter is never bitwise-copied --------------------  
LASRtreehull3d::LASRtreehull3d(const LASRtreehull3d& other) : StageVector(other)
{
	id_attribute = other.id_attribute;
	max_edge = other.max_edge;
	trim = other.trim;
	radius = other.radius;
	wire_filter_strings = other.wire_filter_strings;

	wire_pointfilter = PointFilter(); // rebuilt fresh, never copied directly  
	for (const auto& f : wire_filter_strings)
		wire_pointfilter.add_condition(f);

	// NEW  
	mesh_vector = other.mesh_vector;
	mesh_ofile = other.mesh_ofile;
	mesh_template_filename = other.mesh_template_filename;
}

// -------------------- parameters --------------------  
bool LASRtreehull3d::set_parameters(const nlohmann::json& stage)
{
	id_attribute = stage.value("attribute", "tree_id");
	max_edge = stage.value("max_edge", 0.0);
	radius = stage.value("radius", 0.0);

	trim = max_edge * max_edge; // squared, matches LASRtriangulate's own 'trim' convention  

	if (stage.contains("wire_filter"))
		wire_filter_strings = stage.at("wire_filter").get<std::vector<std::string>>();

	wire_pointfilter = PointFilter();
	for (const auto& f : wire_filter_strings)
		wire_pointfilter.add_condition(f);

	return true;
}

// -------------------- output file wiring (NEW) --------------------  
bool LASRtreehull3d::set_output_file(const std::string& file)
{
	// Let the base class set up the primary hull output (vector, template_filename, written[])  
	if (!StageVector::set_output_file(file)) return false;
	if (file.empty()) return true; // no mesh output requested  

	// Derive the mesh path, preserving the '*' wildcard (if any) at the same  
	// relative position so per-chunk substitution still works in set_input_file_name().  
	size_t dot = file.find_last_of('.');
	if (dot == std::string::npos)
		mesh_template_filename = file + "_mesh";
	else
		mesh_template_filename = file.substr(0, dot) + "_mesh" + file.substr(dot);

	size_t pos = mesh_template_filename.find('*');
	if (pos == std::string::npos)
	{
		// merged/single-file case: create immediately, mirroring StageVector::set_output_file  
		mesh_ofile = mesh_template_filename;
		mesh_vector.set_file(mesh_ofile);
		if (!mesh_vector.create_file()) return false;
		written.push_back(mesh_ofile);
	}
	// else: wildcard case, defer creation to set_input_file_name()  

	return true;
}

bool LASRtreehull3d::set_input_file_name(const std::string& file)
{
	// Let the base class handle the primary hull output's per-chunk file first.  
	if (!StageVector::set_input_file_name(file)) return false;

	if (mesh_template_filename.empty()) return true;

	size_t pos = mesh_template_filename.find('*');
	if (pos != std::string::npos)
	{
		mesh_ofile = mesh_template_filename;
		mesh_ofile.replace(pos, 1, file);
		mesh_vector.set_file(mesh_ofile);
		if (!mesh_vector.create_file()) return false;
		written.push_back(mesh_ofile);
	}

	return true;
}

// -------------------- main processing --------------------  
bool LASRtreehull3d::process(PointCloud*& las)
{
	progress->reset();
	progress->set_prefix("Per-tree 3D hulls");

	// Resolve the tree-id attribute up front and fail early if missing/wrong type  
	int index = las->header->schema.get_attribute_index(id_attribute);
	if (index == -1)
	{
		last_error = "attribute " + id_attribute + " not found in point schema";

		progress->done(); //proper pregress reporting  

		return false;
	}

	AttributeType data_type = las->header->schema.attributes[index].type;
	if (data_type != AttributeType::INT32 && data_type != AttributeType::DOUBLE)
	{
		last_error = "the attribute " + id_attribute + " must be of type 'int' or 'double'";

		progress->done(); //proper pregress reporting  

		return false;
	}

	AttributeAccessor id_accessor(id_attribute);

	// -------------------- Pass 0 (optional): restrict to tree ids near "wire" points --------------------  
	std::unordered_set<int> detected_ids;
	bool restrict_ids = (radius > 0.0) && !wire_filter_strings.empty();

	if (restrict_ids)
	{
		las->build_kdtree();

		Point* p;
		while (las->read_point())
		{
			p = &las->point;
			if (pointfilter.filter(p)) continue;        // respect the stage's own global filter  
			if (wire_pointfilter.filter(p)) continue;    // filter() returns true = REJECTED; skip points that don't pass the wire filter  

			// 'p' passed the wire filter -> it IS a wire point; search its tree neighborhood  
			std::vector<Point> neighbors;
			las->query_sphere(*p, radius, neighbors, &pointfilter);

			for (auto& np : neighbors)
			{
				int id = (int)id_accessor(&np);
				if (id != 0) detected_ids.insert(id); // 0 = "no tree" sentinel from region_growing()  
			}
		}

		las->seek(0);
	}

	// -------------------- Pass 1: bucket points by tree id --------------------  
	std::unordered_map<int, std::vector<PointXYZ>> pts_by_id;

	double xoffset = (las->header->min_x + las->header->max_x) / 2;
	double yoffset = (las->header->min_y + las->header->max_y) / 2;

	std::default_random_engine gen(std::random_device{}());
	std::normal_distribution<double> noise(0.0, 1e-10);

	Point* p;
	while (las->read_point())
	{
		p = &las->point;
		if (pointfilter.filter(p)) continue;

		int id = (int)id_accessor(p);
		if (id == 0) continue;                                  // no tree assigned  
		if (restrict_ids && detected_ids.count(id) == 0) continue; // not near a wire  

		pts_by_id[id].emplace_back(p->get_x() - xoffset + noise(gen),
			p->get_y() - yoffset + noise(gen),
			p->get_z());
	}

	// -------------------- Pass 2: per-tree Delaunay triangulation + contour + polygonize --------------------  
	tree_polygons.clear();
	tree_triangles.clear(); // NEW  

	progress->set_total(pts_by_id.size()); //proper pregress reporting  
	for (auto& kv : pts_by_id)
	{
		int id = kv.first;
		std::vector<PointXYZ>& pts = kv.second;

		(*progress)++; //proper pregress reporting  
		progress->show(); //proper pregress reporting  

		if (pts.size() < 3) continue;

		std::vector<double> coords;
		coords.reserve(pts.size() * 2);
		std::unordered_map<std::pair<long, long>, double, PairHash> z_lookup;

		for (auto& pt : pts)
		{
			coords.push_back(pt.x);
			coords.push_back(pt.y);
			z_lookup[key_of(pt.x, pt.y)] = pt.z;
		}

		delaunator::Delaunator* d = nullptr;
		try
		{
			d = new delaunator::Delaunator(coords);
		}
		catch (const std::exception& e)
		{
			last_error = std::string("In delaunator (tree ") + std::to_string(id) + "): " + e.what();
			continue; // skip this tree, keep processing others  
		}

		// ---- extract boundary edges (same toggle logic as LASRtriangulate::contour) ----  
		// ---- and retain kept facets for the 3D mesh output (NEW) ----  
		std::unordered_set<Edge3D> edges;

		for (unsigned int i = 0; i < d->triangles.size(); i += 3)
		{
			PointXYZ A(coords[2 * d->triangles[i]], coords[2 * d->triangles[i] + 1], 0);
			PointXYZ B(coords[2 * d->triangles[i + 1]], coords[2 * d->triangles[i + 1] + 1], 0);
			PointXYZ C(coords[2 * d->triangles[i + 2]], coords[2 * d->triangles[i + 2] + 1], 0);

			A.z = z_lookup[key_of(A.x, A.y)];
			B.z = z_lookup[key_of(B.x, B.y)];
			C.z = z_lookup[key_of(C.x, C.y)];

			double dx1 = B.x - A.x, dy1 = B.y - A.y;
			double dx2 = C.x - B.x, dy2 = C.y - B.y;
			double dx3 = A.x - C.x, dy3 = A.y - C.y;
			double max_edge2 = std::max({ dx1 * dx1 + dy1 * dy1, dx2 * dx2 + dy2 * dy2, dx3 * dx3 + dy3 * dy3 });

			bool keep = (trim == 0) || (max_edge2 < trim);
			if (!keep) continue;

			// NEW: retain this facet for the 3D mesh output. This produces a 2.5D  
			// triangulated surface (one Z per unique x,y vertex), NOT a  
			// topologically watertight enclosed 3D solid -- it cannot represent  
			// overhangs/concave crowns with multiple Z at the same (x,y). A  
			// literal watertight volume would require replacing this 2D Delaunay  
			// approach (alpha shapes / ball-pivoting / Poisson reconstruction).  
			tree_triangles[id].emplace_back(A, B, C);

			Edge3D AB(A, B), BC(B, C), CA(C, A);
			if (edges.count(AB) > 0) edges.erase(AB); else edges.insert(AB);
			if (edges.count(BC) > 0) edges.erase(BC); else edges.insert(BC);
			if (edges.count(CA) > 0) edges.erase(CA); else edges.insert(CA);
		}

		delete d;

		if (edges.empty()) continue;

		// ---- polygonize (2D geometry collection, Z restored from z_lookup afterward) ----  
		OGRGeometryCollection gc;
		for (const auto& e : edges)
		{
			OGRLineString ls;
			ls.addPoint(e.A.x, e.A.y);
			ls.addPoint(e.B.x, e.B.y);
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
					double x = ext->getX(j);
					double y = ext->getY(j);
					double zz = ext->getZ(j); // NOTE: GDAL Z-preservation through Polygonize() is  
					// unverified for this build - falls back below if 0.  
					if (zz == 0.0)
					{
						auto it = z_lookup.find(key_of(x, y));
						zz = (it != z_lookup.end()) ? it->second : 0.0;
					}
					ring.emplace_back(x, y, zz);
				}

				PolygonXYZ poly(ring);
				poly.close();
				tree_polygons[id] = poly; // one exterior ring per tree; keeps the last if multiple  
			}
		}

		OGRGeometryFactory::destroyGeometry(polys);
	}

	progress->done(); //proper pregress reporting  

	return true;
}

// -------------------- clear / write --------------------  
void LASRtreehull3d::clear(bool last)
{
	tree_polygons.clear();
	tree_triangles.clear(); // NEW  

	if (last)
	{
		vector.finalize_extent();
		mesh_vector.finalize_extent();
	}
}

bool LASRtreehull3d::write()
{
	bool success = true;

#pragma omp critical (write_tree_hull3d)  
	{
		for (const auto& kv : tree_polygons)
		{
			int id = kv.first;
			const PolygonXYZ& poly = kv.second;
			success = success && vector.write(poly, id);
		}
	}

	// NEW: write the per-tree mesh facets  
#pragma omp critical (write_tree_hull3d_mesh)  
	{
		for (const auto& kv : tree_triangles)
		{
			int id = kv.first;
			const std::vector<TriangleXYZ>& tris = kv.second;
			if (!tris.empty())
				success = success && mesh_vector.write(tris, id);
		}
	}

	return success;
}