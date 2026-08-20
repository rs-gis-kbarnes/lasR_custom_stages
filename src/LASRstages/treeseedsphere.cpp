#include "treewireintersect.h"
#include "treeseedsphere.h"
#include "PointCloud.h"

#include <cmath>

LASRtreewireintersect::LASRtreewireintersect()
{
  // Allocated once per top-level instance; copied by shared_ptr (not
  // duplicated) into every clone(), so all clones/chunks share one counter.
  hit_counter = std::make_shared<unsigned int>(0);
}

bool LASRtreewireintersect::set_parameters(const nlohmann::json& stage)
{
  search_radius = stage.value("search_radius", search_radius);

  vector = Vector(xmin, ymin, xmax, ymax);
  vector.set_geometry_type(wkbPoint25D);

  return true;
}

bool LASRtreewireintersect::connect(const std::list<std::unique_ptr<Stage>>& pipeline, const std::string& uid)
{
  Stage* s = search_connection(pipeline, uid);
  if (s == nullptr) return false;

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

  las->build_kdtree();

  const auto& seeds = seed_stage->get_seeds();
  hits.clear();

  for (const auto& seed : seeds)
  {
    double radius = seed.get_extrabyte("radius");
    if (std::isnan(radius) || radius <= 0) continue;

    Sphere sph(seed.x, seed.y, seed.z, radius);

    Point pp(&las->header->schema);
    pp.set_x(seed.x);
    pp.set_y(seed.y);
    pp.set_z(seed.z);

    std::vector<Point> candidates;
    las->query_sphere(pp, search_radius, candidates, &pointfilter);

    for (auto& c : candidates)
    {
      double cx = c.get_x();
      double cy = c.get_y();
      double cz = c.get_z();

      if (!sph.contains(cx, cy, cz)) continue;

      PointLAS hit;
      hit.x = cx;
      hit.y = cy;
      hit.z = cz;

      // Serialize the increment: hit_counter is shared across clones/threads,
      // so a plain (*hit_counter)++ outside a critical section is a data race
      // on the underlying unsigned int, independent of the FID-uniqueness issue.
      #pragma omp critical (assign_wire_hit_ids)
      {
        hit.FID = (*hit_counter)++;
      }

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
  // Do NOT reset hit_counter here -- it must persist across chunks and stay
  // shared across clones. Resetting it reintroduces duplicate-FID collisions.

  if (last) vector.finalize_extent();
}