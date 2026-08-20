#include "treeseedsphere.h"
#include "localmaximum.h"

#include <cmath>
#include <unordered_map>

bool LASRtreeseedsphere::set_parameters(const nlohmann::json& stage)
{
  hag_attribute = stage.value("hag_attribute", hag_attribute);
  radius_multiplier = stage.value("radius_multiplier", radius_multiplier);

  // Output: one wkbPoint25D feature per seed, carrying the derived radius
  // and the HAG value it was derived from as real-valued fields. These are
  // populated automatically by Vector::write(seeds, true) because they are
  // registered here via add_field() and looked up via PointLAS::extrabytes
  // in Vector::write's field-population loop (same mechanism local_maximum
  // uses for its use_attribute field -- see Vector.cpp L129-151).
  vector = Vector(xmin, ymin, xmax, ymax);
  vector.set_geometry_type(wkbPoint25D);
  vector.add_field("radius", OFTReal);
  vector.add_field(hag_attribute, OFTReal);

  return true;
}

bool LASRtreeseedsphere::connect(const std::list<std::unique_ptr<Stage>>& pipeline, const std::string& uid)
{
  Stage* s = search_connection(pipeline, uid);
  if (s == nullptr) return false;

  // Single upstream connection: must be a local_maximum stage.
  // Mirrors LASRnnmetrics::connect() (nnmetrics.cpp L133-150).
  LASRlocalmaximum* p = dynamic_cast<LASRlocalmaximum*>(s);
  if (p)
  {
    set_connection(p);
  }
  else
  {
    last_error = "Incompatible stage combination for 'tree_seed_sphere': expected 'local_maximum'";
    return false;
  }

  return true;
}

bool LASRtreeseedsphere::process(PointCloud*& las)
{
  // 'las' is unused: need_points() == false means the engine may pass a
  // null/unpopulated point cloud here. All data comes from the upstream
  // LASRlocalmaximum connection instead (mirrors LASRnnmetrics::process()
  // pattern of pulling get_maxima() from the connected stage).
  auto it = connections.begin();
  LASRlocalmaximum* lmx = dynamic_cast<LASRlocalmaximum*>(it->second);
  if (lmx == nullptr)
  {
    last_error = "invalid dynamic cast. Expecting a pointer to LASRlocalmaximum";
    return false;
  }

  const auto& maxima = lmx->get_maxima();
  seeds.clear();
  seeds.reserve(maxima.size());

  for (const auto& ttop : maxima)
  {
    double hag = ttop.get_extrabyte(hag_attribute);

    // Skip ttops with no valid HAG (attribute missing or non-positive height).
    if (std::isnan(hag) || hag <= 0) continue;

    double ground_z = ttop.z - hag; // "drop to ground"
    double radius = hag * radius_multiplier; // "buffer to sphere by HAG"

    PointLAS seed;               // default ctor zero-initializes via memset (PointLAS.cpp L10-13)
    seed.x = ttop.x;
    seed.y = ttop.y;
    seed.z = ground_z;
    seed.FID = ttop.FID;         // preserves the tree id from local_maximum's unicity counter

    seed.extrabytes = new std::unordered_map<std::string, double>();
    (*seed.extrabytes)["radius"] = radius;
    (*seed.extrabytes)[hag_attribute] = hag;

    seeds.push_back(std::move(seed));
  }

  return true;
}

bool LASRtreeseedsphere::write()
{
  if (ofile.empty()) return true;
  if (seeds.empty()) return true;

  bool success;
  #pragma omp critical (write_tree_seed_sphere)
  {
    success = vector.write(seeds, true); // true => write standard fields + registered extra fields (radius, HAG)
  }
  return success;
}

void LASRtreeseedsphere::clear(bool last)
{
  seeds.clear();

  // Consistent with LASRtreehull3d::clear()'s established pattern of
  // finalizing the output vector's extent on the last chunk.
  if (last) vector.finalize_extent();
}