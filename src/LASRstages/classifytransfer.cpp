#include "classifytransfer.h"
#include "LASio.h"
#include "Header.h"
#include "PointSchema.h" // Point, AttributeAccessor

#include "nanoflann/nanoflann.h"

#include <cmath>

// Adaptor over plain std::vector<double> buffers for the reference points,
// mirroring PointCloudAdaptor's shape but for an external point set.
struct RefCloudAdaptor
{
  const std::vector<double>& xs;
  const std::vector<double>& ys;
  const std::vector<double>& zs;
  
  RefCloudAdaptor(const std::vector<double>& xs, const std::vector<double>& ys, const std::vector<double>& zs)
    : xs(xs), ys(ys), zs(zs) {}
  
  inline size_t kdtree_get_point_count() const { return xs.size(); }
  inline double kdtree_get_pt(const size_t idx, int dim) const
  {
    if (dim == 0) return xs[idx];
    if (dim == 1) return ys[idx];
    return zs[idx];
  }
  template<class BBOX> bool kdtree_get_bbox(BBOX&) const { return false; }
};

using RefKDTree = nanoflann::KDTreeSingleIndexAdaptor<
  nanoflann::L2_Simple_Adaptor<double, RefCloudAdaptor>, RefCloudAdaptor, 3>;

bool LASRclassifytransfer::process(PointCloud*& las)
{
  if (reference_files.empty())
  {
    last_error = "No reference files provided.";
    return false;
  }
  
  // 1. Read the reference points inside [xmin,ymin,xmax,ymax] + ref_buffer.
  //xmin/ymin/xmax/ymax/circular are inherited Stage members set by set_chunk().
  LASio refio;
  bool ok;
  
  try
  {
    ok = refio.query(reference_files, {}, xmin, ymin, xmax, ymax, ref_buffer, circular, filters);
  }
  catch (const std::exception& e)
  {
    last_error = e.what();
    return false;
  }
  
  if (!ok)
  {
    last_error = "Failed to query reference files.";
    return false;
  }
  
  Header ref_header;
  refio.populate_header(&ref_header);
  
  // NOTE: ref_header.number_of_point_records reflects the WHOLE reference
  // file's point count (raw LAS header field), not the number of points
  // actually inside the [xmin,ymin,xmax,ymax] + ref_buffer query region
  // applied above. Do NOT reserve() against it here: doing so caused every
  // chunk to over-allocate 4 vectors sized to the entire reference file,
  // which multiplies badly with chunk count/concurrency on large files.
  // Let the vectors grow dynamically instead.
  std::vector<double> ref_x, ref_y, ref_z;
  std::vector<int> ref_class;
  
  Point rp(&ref_header.schema);
  AttributeAccessor get_classification("Classification");
  
  while (refio.read_point(&rp))
  {
    ref_x.push_back(rp.get_x());
    ref_y.push_back(rp.get_y());
    ref_z.push_back(rp.get_z());
    ref_class.push_back((int)get_classification(&rp));
  }
  
  refio.close();
  
  if (ref_x.empty())
  {
    // Nothing to transfer for this chunk: not an error, just a no-op.
    return true;
  }
  
  // 2. Build a KD-tree over the reference subset, same params lasR uses
  //internally in PointCloud::build_kdtree() (leaf_max_size = 10).
  RefCloudAdaptor adaptor(ref_x, ref_y, ref_z);
  RefKDTree kdtree(3, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  kdtree.buildIndex();
  
  // 3. Stream the flight chunk and classify.
  AttributeAccessor set_classification("Classification");
  
  progress->reset();
  progress->set_total(las->npoints);
  progress->set_prefix("classify_transfer");
  
  while (las->read_point())
  {
    if (progress->interrupted()) break;
    
    double query_pt[3] = { las->point.get_x(), las->point.get_y(), las->point.get_z() };
    
    RefKDTree::IndexType nn_idx;
    double nn_dist_sq;
    kdtree.knnSearch(query_pt, 1, &nn_idx, &nn_dist_sq);
    
    if (std::sqrt(nn_dist_sq) < dist_threshold && source_classes.count(ref_class[nn_idx]) > 0)
    {
      set_classification(&las->point, (double)target_class);
    }
    
    (*progress)++;
    progress->show();
  }
  
  return true;
}

bool LASRclassifytransfer::set_parameters(const nlohmann::json& stage)
{
  reference_files = get_vector<std::string>(stage["reference_files"]);
  dist_threshold= stage.value("dist_threshold", 3.0);
  ref_buffer= stage.value("ref_buffer", 25.0);
  target_class= stage.value("class", 6);
  
  std::vector<int> classes = get_vector<int>(stage["source_classes"]);
  source_classes = std::unordered_set<int>(classes.begin(), classes.end());
  
  if (reference_files.empty())
  {
    last_error = "'reference_files' must not be empty.";
    return false;
  }
  
  if (source_classes.empty())
  {
    last_error = "'source_classes' must not be empty.";
    return false;
  }
  
  if (dist_threshold <= 0)
  {
    last_error = "'dist_threshold' must be positive.";
    return false;
  }
  
  return true;
}