#ifndef LASRTREEWIREINTERSECT_H
#define LASRTREEWIREINTERSECT_H

#include "Stage.h"
#include "Vector.h"
#include "Shape.h"
#include "PointLAS.h"

#include <memory>
#include <string>
#include <vector>

class LASRtreeseedsphere; // fwd decl, full type only needed in the .cpp

class LASRtreewireintersect : public StageVector
{
public:
  LASRtreewireintersect();
  LASRtreewireintersect(const LASRtreewireintersect& other) = default; // shared_ptr copies fine by default

  bool process(PointCloud*& las) override;
  bool write() override;
  void clear(bool last) override;

  bool need_points() const override { return true; }
  double need_buffer() const override { return search_radius; }
  bool is_parallelized() const override { return true; }   // NEW: explicit, matches localmaximum/nnmetrics precedent

  bool set_parameters(const nlohmann::json&) override;
  bool connect(const std::list<std::unique_ptr<Stage>>&, const std::string& uid) override;

  std::string get_name() const override { return "tree_wire_intersect"; }

  LASRtreewireintersect* clone() const override { return new LASRtreewireintersect(*this); }

private:
  double search_radius = 250.0;
  std::vector<PointLAS> hits;

  // Shared across every clone (multi-threaded/concurrent-files execution) so that
  // FIDs stay globally unique in the single shared GDALDataset/layer that all
  // clones write into (see GDALdataset::dataset being a std::shared_ptr, and
  // Engine's copy constructor comment: "shared resources such as ... gdal
  // datasets"). Mirrors LASRlocalmaximum::counter exactly. Must NOT be reset
  // in clear() -- doing so reintroduces duplicate-FID collisions across
  // chunks/clones, silently dropping hits via Vector::write's dupfid logic.
  std::shared_ptr<unsigned int> hit_counter;
};

#endif