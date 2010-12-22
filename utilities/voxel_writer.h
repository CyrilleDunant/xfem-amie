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

#include "../features/features.h"

namespace Mu
{

typedef enum
{
	F_STIFFNESS,
	F_STRAIN,
	F_STRESS,
} VWFieldType ;

class VoxelWriter
{
protected:
	std::string filename ;
	int nVoxelX ;
	int nVoxelY ;
	int nVoxelZ ;
	int nFields ;
	int cFields ;
	Point bottom_left ;
	Point top_right ;
	bool fullSample ;
		
public:
	VoxelWriter(std::string f, int n, int nf = 1) ;
	VoxelWriter(std::string f, int nx, int ny, int nz, int nf = 1) ;
	VoxelWriter(std::string f, Point bl, Point tr, int n, int nf = 1) ;
	VoxelWriter(std::string f, Point bl, Point tr, int nx, int ny, int nz, int nf = 1) ;
	
	void write(FeatureTree * F, VWFieldType field) ;
	std::valarray<double> getDoubleValues(FeatureTree * F, VWFieldType field) ;
	std::pair<bool, std::vector<double> > getDoubleValue(DelaunayTetrahedron * tet, const Point & p, VWFieldType field) ;
	bool remaining() { return nFields - cFields ; }

private:
	void writeHeader() ;
	
} ;

std::valarray<unsigned short int> normalizeArray(std::valarray<double> val, unsigned short int min = 0, unsigned short int max = 255) ;
int numberOfFields(VWFieldType field) ;



}

#endif
