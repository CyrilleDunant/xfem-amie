//
// C++ Interface: triangle writer
//
// Description: writer for 2D triangle file
//
//
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __TRIANGLE_WRITER_H__
#define __TRIANGLE_WRITER_H__

#include "../../features/features.h"
#include "../../elements/integrable_entity.h"
#include "../itoa.h"

namespace Amie
{

typedef enum
{
	TWFT_COORDINATE,
	TWFT_DISPLACEMENTS,
	TWFT_SCALAR,
	TWFT_DOH,
	TWFT_PRINCIPAL_ANGLE,
    TWFT_TRANSFORM,
    TWFT_VOLUME,
	TWFT_TRIANGLE_ANGLE,
	TWFT_CRACK_ANGLE,
	TWFT_CRITERION,
	TWFT_STIFFNESS,
    TWFT_VISCOSITY,
    TWFT_STIFFNESS_X,
	TWFT_STIFFNESS_Y,
	TWFT_STIFFNESS_Z,
	TWFT_ENRICHMENT,
	TWFT_IMPOSED_STRESS_NORM,
	TWFT_PRINCIPAL_STRESS,
	TWFT_PRINCIPAL_STRAIN,
	TWFT_DAMAGE,
	TWFT_GRADIENT,
	TWFT_FLUX,
	TWFT_GRADIENT_AND_FLUX,
	TWFT_VON_MISES,
	TWFT_CRACKS,
	TWFT_INTERSECTION,
    TWFT_FIELD_TYPE,
    TWFT_INTERNAL_VARIABLE,
    TWFT_LARGE_DEFORMATION_TRANSFORM,
} TWFieldType ;

/** \brief utility class to write the values of various fields of a 2D sample in a text file. You need to use getField() to get your data before writing the file.*/
class TriangleWriter
{
protected:
	std::string filename ;
	std::string head ;
  std::string base ;
	FeatureTree * source ;
	std::vector<std::vector< std::vector<std::valarray<double> > > > values ;
	std::vector<int> nTriangles ;
	std::vector<int> timePlane ;
	std::vector<int> layers ;
	std::vector<TWFieldType> extraFields ;
	std::vector<TWFieldType> fields ;
	std::map<size_t, FieldType> fieldsOther ;
    std::map<size_t, std::string> fieldsInternal ;
	std::map<int, size_t> layerTranslator ;
	int counter ;
	const Geometry * intersection = nullptr;
public:
	
	int getCounter() { return counter ; }
	void nextCounter() 
	{
		counter++ ; 
		filename = base + std::string("_") + itoa(counter) ; 
	}
	
	void setGeometry (const Geometry * in) 
	{
		intersection = in ;
	}
	
	/** \brief simple constructor, get immediately the coordinates of the triangles contained in F */
	TriangleWriter(std::string f, FeatureTree * F = nullptr, int t = 0) ;

	/** \brief clears the data stored*/
	virtual void reset(FeatureTree * F, int t = 0) ;

	/** \brief write the data stored in the writer in a simple txt file, overwriting the file*/
	virtual void write() ;
	virtual void writeSvg(double factor = 1., bool incolor = false ) ;

	/** \brief write the data stored in the writer in a simple txt file*/
	virtual void append() ;

	/** \brief store the values of the field in the writer*/
    virtual void getField(TWFieldType field, bool extra = true, std::string type = std::string(), double offset = 0.) ;

    virtual void getField(std::string field, double offset = 0., bool extra = true ) { getField(TWFT_INTERNAL_VARIABLE, extra, field, offset) ; }

	/** \brief store the values of the field in the writer*/
	virtual void getField(FieldType field, bool extra = true) ;

	/** \brief get the raw values of the specified field from the source*/
    virtual std::vector<std::valarray<double> > getDoubleValues(TWFieldType field, size_t index, int layer, std::string name = std::string() , double offset = 0.) ;

	/** \brief get the raw values of the specified field from the source*/
	virtual std::vector<std::valarray<double> > getDoubleValues(FieldType field, int layer) ;
	
	/** \brief get the raw values of the specific field from a single DelaunayTriangle */
	virtual std::pair<bool, std::vector<double> > getDoubleValue(DelaunayTriangle * tri, TWFieldType field, size_t index) ;

protected:
	/** \brief write the header */
	void writeHeader(int layer, bool append = true) ;

} ;

class BinaryTriangleWriter : public TriangleWriter
{
public:
	BinaryTriangleWriter(std::string f, FeatureTree * F, int t = 0) ;

	virtual void write() ;
	virtual void append() ;

protected:
	virtual void writeHeader(int layer,bool append = true) ;

} ;

class MultiTriangleWriter : public TriangleWriter
{
protected:

public:
    MultiTriangleWriter(std::string head, std::string base, FeatureTree * F, int t = 0) ;

    virtual void append() ;

    void resetCounter() { counter = 0 ; }

protected:
    void writeHeader(int layer, bool firstlayer, bool append = true) ;
    virtual void writeIndexFile() ;
};


/** \brief indicates the number of "columns" needed by a specific field type*/
int numberOfFields(TWFieldType field, size_t index, std::map<size_t, FieldType> & fieldsOther, std::map<size_t, std::string> & external) ;
int numberOfFields(FieldType field) ;

std::string nameOfField(TWFieldType field, size_t index , std::map<size_t, FieldType> & fieldsOther, std::map<size_t, std::string> & external) ;
std::string nameOfField(FieldType field) ;



}

#endif
