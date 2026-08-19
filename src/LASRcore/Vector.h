#ifndef VECTOR_H  
#define VECTOR_H  

#include "GDALdataset.h"  
#include "PointLAS.h"  
#include "Chunk.h"  

#include <limits>  

typedef std::pair<std::string, OGRFieldType> Field;

class Vector : public GDALdataset
{
public:
	Vector();
	Vector(double xmin, double ymin, double xmax, double ymax, int nattr = 1);
	Vector(const Vector& vector, const Chunk& chunk);
	bool create_file();
	bool write(const std::vector<PointLAS>& batch, bool write_attributes = false);
	bool write(const PointXYZAttrs& p);
	bool write(const std::vector<TriangleXYZ>& triangles);
	bool write(const std::vector<PolygonXY>& poly);
	bool write(const PolygonXYZ& poly, int tree_id);
	bool write(const std::vector<TriangleXYZ>& triangles, int tree_id);   // NEW: per-tree mesh  
	bool finalize_extent();                                               // NEW: push accumulated bbox into gpkg_contents  
	void add_field(const std::string& name, OGRFieldType type);
	void set_chunk(const Chunk& chunk);
	int get_dupfid() const { return dupfid; };
	//void set_fields_for(writable type) { writetype = type; };  

private:
	int nattr;
	int writetype;
	int dupfid;
	double extent[4];       // write-time bbox FILTER (from Chunk) -- unrelated to file metadata extent  
	double bbox[4] = { std::numeric_limits<double>::max(),
						std::numeric_limits<double>::max(),
					   -std::numeric_limits<double>::max(),
					   -std::numeric_limits<double>::max() }; // NEW: tracks the actual written-geometry extent  
	std::vector<Field> fields;

	void update_bbox(const Shape& s);        // NEW: for TriangleXYZ and anything else deriving from Shape  
	void update_bbox(const PolygonXYZ& p);   // NEW: PolygonXYZ isn't a Shape, needs its own overload  
};

#endif //VECTOR_H