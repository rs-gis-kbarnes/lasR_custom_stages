#ifndef TREEHULL3D_H 
#define TREEHULL3D_H 
 
#include "Stage.h" 
#include "Vector.h" 
#include "Shape.h" 
#include <unordered_map> 
#include <unordered_set> 
 
namespace delaunator { class Delaunator; } 
 
class LASRtreehull3d : public StageVector 
{ 
public: 
  LASRtreehull3d(); //declare only, imp in cpp
  LASRtreehull3d(const LASRtreehull3d& other); // explicit: must rebuild wire_pointfilter, see .cpp 
 
  bool process(PointCloud*& las) override; 
  void clear(bool last) override; 
  bool write() override; 
  bool need_points() const override { return true; } 
  bool is_streamable() const override { return false; } // needs full per-tree point sets in one chunk 
  bool set_parameters(const nlohmann::json&) override; 
  std::string get_name() const override { return "tree_hull3d"; } 
  double need_buffer() const override { return radius > 0 ? radius : 0.0; } 
  bool is_parallelized() const override { return true; } 
  LASRtreehull3d* clone() const override { return new LASRtreehull3d(*this); } 
 
private: 
  std::string id_attribute = "tree_id"; 
  double max_edge = 0; 
  double trim = 0;
  double radius = 0; 
  std::vector<std::string> wire_filter_strings; // raw filter strings, source of truth for wire_pointfilter 
  PointFilter wire_pointfilter;                 // rebuilt from wire_filter_strings in every copy, never bitwise-copied 
 
  std::unordered_map<int, std::vector<double>> coords_by_id; 
  std::unordered_map<int, std::vector<PointXYZ>> pts_by_id; 
  std::unordered_map<int, PolygonXYZ> tree_polygons;
}; 
 
#endif