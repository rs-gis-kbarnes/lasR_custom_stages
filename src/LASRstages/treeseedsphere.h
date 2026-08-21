#ifndef LASRTREESEEDSPHERE_H
#define LASRTREESEEDSPHERE_H

#include "Stage.h"
#include "Vector.h"
#include "PointLAS.h"

#include <vector>
#include <string>

class LASRlocalmaximum; // fwd decl, full type only needed in the .cpp (includes localmaximum.h)

// Stage 1 of the tree/wire workflow.
//
// Consumes ttops produced by an upstream local_maximum() stage and, for each
// ttop that has a valid HAG (height-above-ground) extrabyte value, derives a
// "seed" point representing a sphere:
//   center = (ttop.x, ttop.y, ground_z)   where ground_z = ttop.z - HAG   ("drop to ground")
//   radius = HAG * radius_multiplier                                      ("buffer to sphere by HAG")
//
// This stage performs pure geometry derivation from the upstream stage's
// in-memory results (LASRlocalmaximum::get_maxima()) and does not need to
// touch any LAS points itself -- hence need_points() == false.
class LASRtreeseedsphere : public StageVector
{
public:
  LASRtreeseedsphere() = default;

  // No per-point access needed; only reads upstream maxima via connections.
  bool process(PointCloud*& las) override;
  bool write() override;
  void clear(bool last) override;

  // Pure geometry derivation from ttops already resolved by local_maximum;
  // no LAS points are read by this stage.
  bool need_points() const override { return false; }

  bool LASRtreeseedsphere::set_parameters(const nlohmann::json& stage)
  {
    hag_attribute = stage.value("hag_attribute", hag_attribute);
    radius_multiplier = stage.value("radius_multiplier", radius_multiplier);

    vector = Vector(xmin, ymin, xmax, ymax);
    vector.set_geometry_type(wkbPoint25D);

    // Vector::write(seeds, true) unconditionally populates these 5 standard
    // fields (see Vector.cpp write_attributes==true branch) regardless of
    // whether the values are meaningful for a derived seed point (they are
    // all left at their zero-initialized default from PointLAS's memset
    // ctor here). They MUST be registered via add_field() or GDAL logs
    // "ERROR 1: Invalid index : -1" once per field per feature.
    vector.add_field("Intensity", OFTInteger);
    vector.add_field("gpstime", OFTReal);
    vector.add_field("ReturnNumber", OFTInteger);
    vector.add_field("Classification", OFTInteger);
    vector.add_field("ScanAngle", OFTReal);

    // Fields actually meaningful for this stage.
    vector.add_field("radius", OFTReal);
    vector.add_field(hag_attribute, OFTReal);

    return true;
  }


  bool connect(const std::list<std::unique_ptr<Stage>>&, const std::string& uid) override;

  std::string get_name() const override { return "tree_seed_sphere"; }

  // Accessor used by the downstream LASRtreewireintersect stage.
  std::vector<PointLAS>& get_seeds() { return seeds; };

  LASRtreeseedsphere* clone() const override { return new LASRtreeseedsphere(*this); };

private:
  // --- parameters (from JSON via set_parameters()) ---
  std::string hag_attribute = "HAG";   // name of the HAG extrabyte stored on each ttop by local_maximum()
  double radius_multiplier = 1.0;      // sphere radius = HAG * radius_multiplier

  // --- per-chunk transient output ---
  // One PointLAS per valid ttop: x/y = ttop x/y, z = ground z (ttop.z - HAG),
  // FID = ttop.FID (tree id), extrabytes["radius"] and extrabytes[hag_attribute]
  // carry the derived sphere radius and the HAG value used to compute it.
  std::vector<PointLAS> seeds;
};

#endif // LASRTREESEEDSPHERE_H