#include "treehull3d.h"  
#include "PointCloud.h"  
#include "Progress.h"  

#include "delaunator.hpp"  

#include <random>  

LASRtreehull3d::LASRtreehull3d()
{
	vector.set_geometry_type(wkbPolygon25D);
	vector.add_field("tree_id", OFTInteger);

	mesh_vector.set_geometry_type(wkbMultiPolygon25D);
	mesh_vector.add_field("tree_id", OFTInteger);
}

LASRtreehull3d::LASRtreehull3d(const LASRtreehull3d& other) : StageVector(other)
{
	id_attribute = other.id_attribute;
	max_edge = other.max_edge;
	trim = other.trim;
	radius = other.radius;
	wire_filter_strings = other.wire_filter_strings;
	wire_pointfilter = PointFilter(); // must be rebuilt, never bitwise-copied  
	for (const auto& f : wire_filter_strings) wire_pointfilter.add_condition(f);

	mesh_vector = other.mesh_vector;
	mesh_ofile = other.mesh_ofile;
	mesh_template_filename = other.mesh_template_filename;
}

bool LASRtreehull3d::set_parameters(const nlohmann::json& stage)
{
	id_attribute = stage.value("attribute", id_attribute);
	max_edge = stage.value("max_edge", max_edge);
	radius = stage.value("radius", radius);
	trim = max_edge * max_edge;

	if (stage.contains("wire_filter"))
	{
		wire_filter_strings = stage["wire_filter"].get<std::vector<std::string>>();
		for (const auto& f : wire_filter_strings) wire_pointfilter.add_condition(f);
	}

	return true;
}

void LASRtreehull3d::set_crs(const CRS& crs)
{
	StageVector::set_crs(crs); // sets Stage::crs and vector.set_crs(crs)  
	mesh_vector.set_crs(crs);  // NEW: mesh_vector was never getting the CRS before this fix  
}

bool LASRtreehull3d::set_output_file(const std::string& file)
{
	if (!StageVector::set_output_file(file)) return false;
	if (file.empty()) return true;

	size_t dot = file.find_last_of('.');
	if (dot == std::string::npos)
		mesh_template_filename = file + "_mesh";
	else
		mesh_template_filename = file.substr(0, dot) + "_mesh" + file.substr(dot);

	size_t pos = mesh_template_filename.find('*');
	if (pos == std::string::npos)
	{
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
	if (!StageVector::set_input_file_name(file)) return false;
	if (mesh_template_filename.empty()) return true;

	std::string ofile2 = mesh_template_filename;
	size_t pos = ofile2.find('*');
	if (pos != std::string::npos)
	{
		ofile2.replace(pos, 1, file);
		mesh_ofile = ofile2;
		mesh_vector.set_file(mesh_ofile);
		if (!mesh_vector.create_file()) return false;
		written.push_back(mesh_ofile);
	}

	return true;
}

bool LASRtreehull3d::process(PointCloud*& las)
{
	progress->reset();
	progress->set_prefix("Per-tree 3D hulls");

	int index = las->header->schema.get_attribute_index(id_attribute);
	if (index == -1)
	{
		last_error = "attribute " + id_attribute + " not found in point schema";
		progress->done();
		return false;
	}

	AttributeType data_type = las->header->schema.attributes[index].type;
	if (data_type != AttributeType::INT32 && data_type != AttributeType::DOUBLE)
	{
		last_error = "the attribute " + id_attribute + " must be of type 'int' or 'double'";
		progress->done();
		return false;
	}

	AttributeAccessor id_accessor(id_attribute);

	bool restrict_ids = !wire_filter_strings.empty();
	std::unordered_set<int> detected_ids;

	if (restrict_ids)
	{
		Point* p;
		while (las->read_point())
		{
			p = &las->point;
			if (pointfilter.filter(p)) continue;
			if (wire_pointfilter.filter(p)) continue; // filter() returns true = REJECTED; skip points that don't pass the wire filter  

			std::vector<Point> neighbors;
			las->query_sphere(*p, radius, neighbors, &pointfilter);
			for (auto& n : neighbors)
				detected_ids.insert((int)id_accessor(&n));
		}
		las->seek(0);
	}

	// -------------------- per-tree point collection, TRUE coordinates + index_map --------------------  
	// Mirrors LASRtriangulate::process(): shift-for-Delaunator-stability coords are kept separate  
	// from the index that maps back to the ORIGINAL point, so no offset ever needs to be re-applied  
	// and no z_lookup/quantization hack is needed to recover Z.  

	double xoffset = (las->header->min_x + las->header->max_x) / 2;
	double yoffset = (las->header->min_y + las->header->max_y) / 2;

	std::default_random_engine gen(std::random_device{}());
	std::normal_distribution<double> noise(0.0, 1e-10);

	std::unordered_map<int, std::vector<double>> coords_by_id; // shifted x,y pairs fed to Delaunator only  

	Point* p;
	while (las->read_point())
	{
		p = &las->point;
		if (pointfilter.filter(p)) continue;

		int id = (int)id_accessor(p);
		if (id == 0) continue;                                    // no tree assigned  
		if (restrict_ids && detected_ids.count(id) == 0) continue; // not near a wire  

		coords_by_id[id].push_back(p->get_x() - xoffset + noise(gen));
		coords_by_id[id].push_back(p->get_y() - yoffset + noise(gen));
		index_by_id[id].push_back(las->current_point); // original point index, for las->get_point() later  
	}

	progress->set_total(coords_by_id.size());

	for (auto& kv : coords_by_id)
	{
		int id = kv.first;
		std::vector<double>& coords = kv.second;
		std::vector<int>& index_map = index_by_id[id];

		if (coords.size() < 6) continue; // fewer than 3 points  

		delaunator::Delaunator d(coords);

		std::vector<TriangleXYZ> kept_triangles;
		std::unordered_map<Edge3D, int> edge_count;

		Point A, B, C;
		A.set_schema(&las->header->schema);
		B.set_schema(&las->header->schema);
		C.set_schema(&las->header->schema);

		for (std::size_t i = 0; i < d.triangles.size(); i += 3)
		{
			int ia = index_map[d.triangles[i]];
			int ib = index_map[d.triangles[i + 1]];
			int ic = index_map[d.triangles[i + 2]];

			las->get_point(ia, &A);
			las->get_point(ib, &B);
			las->get_point(ic, &C);

			TriangleXYZ tri(A, B, C);

			bool keep = (trim == 0) || (tri.square_max_edge_size() < trim);
			if (!keep) continue;

			tri.make_clock_wise();
			kept_triangles.push_back(tri);

			PointXYZ pa(A.get_x(), A.get_y(), A.get_z());
			PointXYZ pb(B.get_x(), B.get_y(), B.get_z());
			PointXYZ pc(C.get_x(), C.get_y(), C.get_z());

			Edge3D e1(pa, pb), e2(pb, pc), e3(pc, pa);
			edge_count[e1]++;
			edge_count[e2]++;
			edge_count[e3]++;
		}

		if (kept_triangles.empty()) continue;

		tree_triangles[id] = kept_triangles; // NEW: store mesh facets for this tree  

		// ---- boundary edges: appear exactly once (not shared by two kept triangles) ----  
		std::vector<Edge3D> edges;
		for (auto& ec : edge_count)
			if (ec.second == 1) edges.push_back(ec.first);

		if (edges.empty()) continue;

		OGRGeometryCollection gc;
		for (const auto& e : edges)
		{
			OGRLineString ls;
			ls.addPoint(e.A.x, e.A.y);
			ls.addPoint(e.B.x, e.B.y);
			gc.addGeometry(&ls);
		}

		OGRGeometry* merged = nullptr;
		OGRGeometry* polys = nullptr;
		merged = gc.Union(nullptr); // # nocov safe-guard if Union unsupported handled below  
		if (merged)
		{
			polys = OGRGeometryFactory::organizePolygons(&merged, 1, nullptr, nullptr);
		}

		if (polys && wkbFlatten(polys->getGeometryType()) == wkbPolygon)
		{
			OGRPolygon* p2 = polys->toPolygon();
			OGRLinearRing* ring = p2->getExteriorRing();

			std::vector<PointXYZ> coords_out;
			for (int i = 0; i < ring->getNumPoints(); i++)
				coords_out.emplace_back(ring->getX(i), ring->getY(i), 0.0);

			PolygonXYZ poly(coords_out);
			poly.close();
			tree_polygons[id] = poly;
		}

		if (polys) OGRGeometryFactory::destroyGeometry(polys);
	}

	progress->done();

	return true;
}

// -------------------- clear / write --------------------  

void LASRtreehull3d::clear(bool last)
{
	pts_by_id.clear();
	index_by_id.clear();
	tree_polygons.clear();
	tree_triangles.clear();

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
		TransactionGuard tg(vector); // one transaction for the whole chunk's trees  
		for (const auto& kv : tree_polygons)
		{
			int id = kv.first;
			const PolygonXYZ& poly = kv.second;
			success = success && vector.write(poly, id);
		}
	}

#pragma omp critical (write_tree_hull3d_mesh)  
	{
		TransactionGuard tg(mesh_vector); // one transaction for the whole chunk's mesh facets  
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