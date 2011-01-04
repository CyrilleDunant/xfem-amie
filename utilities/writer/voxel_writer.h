//
// C++ Interface: voxel writer
//
// Description: writer for binary voxel file
//
//
// Author: Alain Giorla
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __VOXEL_WRITER_H__
#define __VOXEL_WRITER_H__

#include "../../features/features.h"

namespace Mu
{

typedef enum
{
	VWFT_PRINCIPAL_ANGLE,
	VWFT_STIFFNESS,
	VWFT_STRAIN,
	VWFT_STRESS,
	VWFT_STRAIN_AND_STRESS,
	VWFT_VON_MISES,
} VWFieldType ;

/** \brief utility class to write various fields of a 3D sample in a binary file. You need to call getField() to get your data before writing the file.*/
class VoxelWriter
{
protected:
	std::string filename ;
	std::vector<std::valarray<unsigned short int> > values ;
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
	static void writeMap(std::string filename, FeatureTree * F, Variable axis, double pos, int n, VWFieldType field, int k = 0, int min = 0, int max = 255) ;

private:
	/** \brief write the header of the binary file*/
	void writeHeader() ;
	
} ;

/** \brief normalize an array of double to an array of unsigned short integers*/
std::valarray<unsigned short int> normalizeArray(std::valarray<double> val, unsigned short int min = 0, unsigned short int max = 255) ;

/** \brief number of different fields contained within a VoxelWriterFieldType, typically 1*/
const int numberOfFields(VWFieldType field) ;



}

#endif
