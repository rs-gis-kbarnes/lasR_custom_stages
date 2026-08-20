#include "treewireintersect.h"
#include "treeseedsphere.h"
#include "PointCloud.h"

#include <cmath>

bool LASRtreewireintersect::set_parameters(const nlohmann::json& stage)
{
  search_radius = stage.value("search_radius", search_radius);

  // Output: plain xyz points, no attribute fields registered. write(hits, false)
  // skips all field population (Vector.cpp write_attributes==false branch), so
  // no add_field() calls are needed here.
  vector = Vector(xmin, ymin, xmax, ymax);
  vector.set_geometry_type(wkbPoint25D);

  return true;
}

bool LASRtreewireintersect::connect(const std::list<std::unique_ptr<Stage>>& pipeline, const std::string& uid)
{
  Stage* s = search_connection(pipeline, uid);
  if (s == nullptr) return false;

  // Single upstream connection: must be a tree_seed_sphere stage.
  LASRtreeseedsphere* p = dynamic_cast<LASRtreeseedsphere*>(s);
  if (p)
  {
    set_connection(p);
  }
  else
  {
    last_error = "Incompatible stage combination for 'tree_wire_intersect': expected 'tree_seed_sphere'";
    return false;
  }

  return true;
}

bool LASRtreewireintersect::process(PointCloud*& las)
{
  auto it = connections.begin();
  LASRtreeseedsphere* seed_stage = dynamic_cast<LASRtreeseedsphere*>(it->second);
  if (seed_stage == nullptr)
  {
    last_error = "invalid dynamic cast. Expecting a pointer to LASRtreeseedsphere";
    return false;
  }

  // query_sphere() throws std::runtime_error if the KD-tree is not built yet
  // (PointCloud.cpp L600-604). Build it once per chunk before the loop,
  // mirroring LASRnnmetrics::process() (nnmetrics.cpp L58-59).
  las->build_kdtree();

  const auto& seeds = seed_stage->get_seeds();
  hits.clear();

  for (const auto& seed : seeds)
  {
    double radius = seed.get_extrabyte("radius");
    if (std::isnan(radius) || radius <= 0) continue;

    // Exact geometric test shape: true HAG-derived sphere centered at the
    // seed (already dropped to ground by tree_seed_sphere).
    Sphere sph(seed.x, seed.y, seed.z, radius);

    Point pp(&las->header->schema);
    pp.set_x(seed.x);
    pp.set_y(seed.y);
    pp.set_z(seed.z);

    // Coarse speed-up query: over-estimated search_radius, restricted by
    // this stage's own generic pointfilter (populated from the pipeline's
    // filter= argument -- e.g. filter = keep_class(14) selects "wire"
    // points as candidates). This is NOT a generic input-restriction filter
    // in the usual sense; for this stage it specifically plays the role of
    // the "wire/strike" class selector.
    std::vector<Point> candidates;
    las->query_sphere(pp, search_radius, candidates, &pointfilter);

    for (auto& c : candidates)
    {
      double cx = c.get_x();
      double cy = c.get_y();
      double cz = c.get_z();

      // Exact test against the true (HAG-derived) sphere, not search_radius.
      if (!sph.contains(cx, cy, cz)) continue;

      PointLAS hit;              // zero-initialized via memset (PointLAS.cpp L10-13)
      hit.x = cx;
      hit.y = cy;
      hit.z = cz;
      hit.FID = hit_counter++;   // unique FID per hit -- required, see header comment

      hits.push_back(std::move(hit));
    }
  }

  return true;
}

bool LASRtreewireintersect::write()
{
  if (ofile.empty()) return true;
  if (hits.empty()) return true;

  bool success;
  #pragma omp critical (write_tree_wire_intersect)
  {
    success = vector.write(hits, false); // false => plain xyz, no fields
  }
  return success;
}

void LASRtreewireintersect::clear(bool last)
{
  hits.clear();
  hit_counter = 0;

  if (last) vector.finalize_extent();
}