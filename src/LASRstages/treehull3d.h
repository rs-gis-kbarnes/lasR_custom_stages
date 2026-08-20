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
    LASRtreehull3d();
    LASRtreehull3d(const LASRtreehull3d& other);

    bool process(PointCloud*& las) override;
    void clear(bool last) override;
    bool write() override;
    bool need_points() const override { return true; }
    bool is_streamable() const override { return false; } // needs full per-tree point sets in one chunk  
    bool set_parameters(const nlohmann::json&) override;
    bool set_output_file(const std::string& file) override;      // derives mesh_ofile too  
    bool set_input_file_name(const std::string& file) override;  // per-chunk mesh file creation  
    void set_crs(const CRS& crs) override;                       // NEW: propagate CRS to mesh_vector too  
    LASRtreehull3d* clone() const override { return new LASRtreehull3d(*this); }

private:
    std::string id_attribute = "tree_id";
    double max_edge = 0;
    double trim = 0;
    double radius = 0;
    std::vector<std::string> wire_filter_strings; // raw filter strings, source of truth for wire_pointfilter  
    PointFilter wire_pointfilter;                 // rebuilt from wire_filter_strings in every copy, never bitwise-copied  

    std::unordered_map<int, std::vector<PointXYZ>> pts_by_id;   // kept for compatibility; positions are TRUE coords now  
    std::unordered_map<int, std::vector<int>> index_by_id;      // NEW: original las point index per tree, mirrors triangulate.cpp  
    std::unordered_map<int, PolygonXYZ> tree_polygons;
    std::unordered_map<int, std::vector<TriangleXYZ>> tree_triangles; // NEW: kept mesh facets per tree  

    Vector mesh_vector;                 // NEW: second output for the 3D mesh  
    std::string mesh_ofile;
    std::string mesh_template_filename; // NEW: mirrors StageVector's template_filename, for wildcard mesh paths  
};

#endif //TREEHULL3D_H