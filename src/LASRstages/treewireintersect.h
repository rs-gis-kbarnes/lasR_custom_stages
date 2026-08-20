#ifndef LASRTREEWIREINTERSECT_H
#define LASRTREEWIREINTERSECT_H

#include "Stage.h"
#include "Vector.h"
#include "Shape.h"
#include "PointLAS.h"

#include <vector>
#include <string>

class LASRtreeseedsphere; // fwd decl, full type only needed in the .cpp

// Stage 2 of the tree/wire workflow.
//
// Consumes sphere seeds produced by an upstream tree_seed_sphere() stage.
// For each seed:
//   1. Coarse speed-up: las->query_sphere(seed_center, search_radius, ...)
//      restricted by this stage's own generic filter= (e.g. keep_class(14)
//      to select "wire"/strike points) -- this is the "find trees within
//      arbitrary large distance of wire points" speed-up step.
//   2. Exact test: Sphere::contains(x,y,z) on each coarse candidate, using
//      the seed's real (HAG-derived) radius, not search_radius.
// Surviving points are written as plain xyz features (no attribute fields).
class LASRtreewireintersect : public StageVector
{
public:
  LASRtreewireintersect() = default;

  bool process(PointCloud*& las) override;
  bool write() override;
  void clear(bool last) override;

  bool need_points() const override { return true; }
  double need_buffer() const override { return search_radius; } // cross-tile candidates must be included
  bool is_streamable() const override { return false; }         // needs all seeds resolved via connect()

  bool set_parameters(const nlohmann::json&) override;
  bool connect(const std::list<std::unique_ptr<Stage>>&, const std::string& uid) override;

  std::string get_name() const override { return "tree_wire_intersect"; }

  LASRtreewireintersect* clone() const override { return new LASRtreewireintersect(*this); };

private:
  // --- parameters ---
  double search_radius = 250.0; // coarse over-estimate speed-up radius (map units)

  // --- per-chunk transient state ---
  std::vector<PointLAS> hits;

  // Monotonically incrementing counter used to assign a unique OGR FID to
  // every hit. NOTE: PointLAS's default ctor zero-initializes FID via memset,
  // so without this counter every hit would get FID == 0 and
  // Vector::write()'s duplicate-FID check (Vector.cpp L111-118) would
  // silently drop every hit after the first one written in the whole run.
  unsigned int hit_counter = 0;
};

#endif // LASRTREEWIREINTERSECT_H