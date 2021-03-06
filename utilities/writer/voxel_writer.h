//
// C++ Interface: voxel writer
//
// Description: writer for binary voxel file
//
//
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __VOXEL_WRITER_H__
#define __VOXEL_WRITER_H__

#include "../../features/features.h"

namespace Amie
{

typedef enum
{
	VWFT_PRINCIPAL_ANGLE,
	VWFT_STIFFNESS,
	VWFT_STRAIN,
	VWFT_STRESS,
        VWFT_PRINCIPAL_STRAIN,
        VWFT_PRINCIPAL_STRESS,
	VWFT_STRAIN_AND_STRESS,
	VWFT_CONCENTRATION,
	VWFT_GRADIENT,
	VWFT_FLUX,
	VWFT_GRADIENT_AND_FLUX,
	VWFT_VON_MISES,
	VWFT_ENRICHEMENT,
	VWFT_DAMAGE
} VWFieldType ;

/** \brief utility class to write various fields of a 3D sample in a binary file. You need to call getField() to get your data before writing the file.*/
class VoxelWriter
{
protected:
	std::string filename ;
	std::vector<std::valarray<unsigned char> > values ;
	std::valarray<bool> voids ; 
	int nVoxelX ;
	int nVoxelY ;
	int nVoxelZ ;
	Point bottom_left ;
	Point top_right ;
	bool fullSample ;
		
public:
	/** \brief simple constructor, with same number of sampling points on all axis*/
	VoxelWriter(std::string f, int n) ;
	
	/** \brief simple constructor, with different number of sampling points on all axis*/
	VoxelWriter(std::string f, int nx, int ny, int nz) ;
	
	/** \brief simple constructor, with same number of sampling points on all axis, and which does not take the full sample, but a cube defined by the bottom left and the top right corners*/
	VoxelWriter(std::string f, Point bl, Point tr, int n) ;
	
	/** \brief simple constructor, with same number of sampling points on all axis, and which does not take the full sample, but a cube defined by the bottom left and the top right corners*/
	VoxelWriter(std::string f, Point bl, Point tr, int nx, int ny, int nz) ;
	
	/** \brief write the data stored in the writer in a binary file*/
	void write() ;

	/** \brief get the values of the specified field in the FeatureTree, and store it in the writer.*/
	void getField(FeatureTree * F, VWFieldType field) ;	
	
	/** \brief get the values of the specified field in the FeatureTree*/
	std::vector<std::valarray<double> > getDoubleValues(FeatureTree * F, VWFieldType field) ;
	
	/** \brief get the values of a specific field at a specific point in a specific tetrahedron*/
	std::pair<bool, std::vector<double> > getDoubleValue(DelaunayTetrahedron * tet, const Point & p, VWFieldType field) ;

	/** \brief number of sampling points*/
	int nPoints() { return nVoxelX*nVoxelY*nVoxelZ ; }
	
	/** \brief utility static method to write a 2D slice of a 3D sample, according to a plane parallel to the specified axis.*/
	void writeMap(std::string filename, FeatureTree * F, Variable axis, double pos, int n, VWFieldType field, int k = 0, int min = 1, int max = 255) ;

private:
	/** \brief write the header of the binary file*/
	void writeHeader() ;
	
} ;

/** \brief normalize an array of double to an array of unsigned short integers*/
std::valarray<unsigned char> normalizeArray(const std::valarray< double >& val, const std::valarray< bool >& voids, short unsigned int min = 1, short unsigned int max = 255) ;

/** \brief number of different fields contained within a VoxelWriterFieldType, typically 1*/
int numberOfFields(VWFieldType field) ;



}

#endif
