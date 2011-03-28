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

namespace Mu
{

typedef enum
{
	TWFT_COORDINATE,
	TWFT_DISPLACEMENTS,
	TWFT_PRINCIPAL_ANGLE,
	TWFT_STIFFNESS,
	TWFT_STRAIN,
	TWFT_STRESS,
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
	std::vector<std::valarray<double> > values ;
	int nTriangles ;
		
public:
	/** \brief simple constructor, get immediately the coordinates of the triangles contained in F */
	TriangleWriter(std::string f, FeatureTree * F) ;

	/** \brief write the data stored in the writer in a simple txt file*/
	void write() ;
	
	/** \brief store the values of the field in the writer*/
	void getField(TWFieldType field) ;	
	
	/** \brief get the raw values of the specified field from the source*/
	std::vector<std::valarray<double> > getDoubleValues(TWFieldType field) ;
	
	/** \brief get the raw values of the specific field from a single DelaunayTriangle */
	std::pair<bool, std::vector<double> > getDoubleValue(DelaunayTriangle * tri, TWFieldType field) ;

private:
	/** \brief write the header */
	void writeHeader() ;
	
} ;

/** \brief indicates the number of "columns" needed by a specific field type*/
int numberOfFields(TWFieldType field) ;



}

#endif
