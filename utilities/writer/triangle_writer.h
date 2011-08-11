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
#include "../itoa.h"

namespace Mu
{

typedef enum
{
	TWFT_COORDINATE,
	TWFT_DISPLACEMENTS,
	TWFT_PRINCIPAL_ANGLE,
	TWFT_CRACK_ANGLE,
	TWFT_CRITERION,
	TWFT_STIFFNESS,
	TWFT_STRAIN,
	TWFT_STRESS,
	TWFT_PRINCIPAL_STRESS,
	TWFT_PRINCIPAL_STRAIN,
	TWFT_STRAIN_AND_STRESS,
	TWFT_DAMAGE,
	TWFT_GRADIENT,
	TWFT_FLUX,
	TWFT_GRADIENT_AND_FLUX,
	TWFT_VON_MISES,
} TWFieldType ;

/** \brief utility class to write the values of various fields of a 2D sample in a text file. You need to use getField() to get your data before writing the file.*/
class TriangleWriter
{
protected:
	std::string filename ;
	FeatureTree * source ;
	std::vector< std::vector<std::valarray<double> > > values ;
	std::vector<int> nTriangles ;
	std::vector<int> timePlane ;
	std::vector<int> layers ;
	std::vector<TWFieldType> extraFields ;
	std::map<int, size_t> layerTranslator ;

public:
	/** \brief simple constructor, get immediately the coordinates of the triangles contained in F */
	TriangleWriter(std::string f, FeatureTree * F = NULL, int t = 0) ;

	/** \brief clears the data stored*/
	virtual void reset(FeatureTree * F, int t = 0) ;

	/** \brief write the data stored in the writer in a simple txt file, overwriting the file*/
	virtual void write() ;

	/** \brief write the data stored in the writer in a simple txt file*/
	virtual void append() ;

	/** \brief store the values of the field in the writer*/
	virtual void getField(TWFieldType field, bool extra = true) ;

	/** \brief get the raw values of the specified field from the source*/
	virtual std::vector<std::valarray<double> > getDoubleValues(TWFieldType field, int layer) ;

	/** \brief get the raw values of the specific field from a single DelaunayTriangle */
	virtual std::pair<bool, std::vector<double> > getDoubleValue(DelaunayTriangle * tri, TWFieldType field) ;

protected:
	/** \brief write the header */
	virtual void writeHeader(int layer, bool append = true) ;

} ;

class BinaryTriangleWriter : public TriangleWriter
{
public:
	BinaryTriangleWriter(std::string f, FeatureTree * F, int t = 0) ;

	virtual void write() ;
	virtual void append() ;

protected:
	virtual void writeHeader(bool append = true) ;

} ;

class MultiTriangleWriter : public TriangleWriter
{
protected:
    int counter ;
    std::string head ;
    std::string base ;

public:
    MultiTriangleWriter(std::string head, std::string base, FeatureTree * F, int t = 0) ;

    virtual void append() ;
    int getCounter() { return counter ; }
    void nextCounter() { counter++ ; filename = base + "_" + itoa(counter) ; }
    void resetCounter() { counter = 0 ; }

protected:
    virtual void writeHeader(int layer, bool firstlayer, bool append = true) ;
};


/** \brief indicates the number of "columns" needed by a specific field type*/
int numberOfFields(TWFieldType field) ;



}

#endif
