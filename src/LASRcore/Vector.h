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
	bool write(const PolygonXYZ& poly, int tree_id);                    // per-tree hull  
	bool write(const std::vector<TriangleXYZ>& triangles, int tree_id); // per-tree mesh  
	bool finalize_extent();
	void add_field(const std::string& name, OGRFieldType type);
	void set_chunk(const Chunk& chunk);
	int get_dupfid() const { return dupfid; };
	//void set_fields_for(writable type) { writetype = type; };  

private:
	int nattr;
	int writetype;
	int dupfid;
	double extent[4];
	double bbox[4] = { std::numeric_limits<double>::max(),
						std::numeric_limits<double>::max(),
					   -std::numeric_limits<double>::max(),
					   -std::numeric_limits<double>::max() };
	std::vector<Field> fields;

	void update_bbox(const Shape& s);        // for TriangleXYZ and anything deriving from Shape  
	void update_bbox(const PolygonXYZ& p);   // PolygonXYZ is not a Shape, needs its own overload  
};

// Wraps a batch of writes in a single GDAL/OGR transaction so per-feature  
// CreateFeature() calls don't each trigger an implicit SQLite commit/fsync.  
// Construct ONE of these around an entire write loop (e.g. all trees in a  
// chunk), not once per feature/call — that's what actually amortizes the I/O.  
class TransactionGuard
{
public:
	TransactionGuard(Vector& v) : ds(v.dataset.get())
	{
		active = ds && (ds->StartTransaction() == OGRERR_NONE);
	}
	~TransactionGuard()
	{
		if (active) ds->CommitTransaction();
	}

	TransactionGuard(const TransactionGuard&) = delete;
	TransactionGuard& operator=(const TransactionGuard&) = delete;

private:
	GDALDataset* ds;
	bool active;
};

#endif //VECTOR_H