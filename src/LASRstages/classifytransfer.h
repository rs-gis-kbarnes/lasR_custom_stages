#ifndef CLASSIFYTRANSFER_H
#define CLASSIFYTRANSFER_H

#include "Stage.h"

#include <unordered_set>
#include <vector>
#include <string>

class LASRclassifytransfer: public Stage
{
public:
  LASRclassifytransfer() = default;
  bool process(PointCloud*& las) override;
  double need_buffer() const override { return dist_threshold; };
  bool set_parameters(const nlohmann::json&) override;
  std::string get_name() const override { return "classify_transfer"; };

  // multi-threading: base Stage(const Stage&) copy ctor already copies
  // progress/ncpu/filters/crs; our own members are trivially copyable.
  LASRclassifytransfer* clone() const override { return new LASRclassifytransfer(*this); };

private:
  std::vector<std::string> reference_files;
  double dist_threshold = 3.0;
  double ref_buffer = 25.0;
  std::unordered_set<int> source_classes;
  int target_class = 6;
};
#endif