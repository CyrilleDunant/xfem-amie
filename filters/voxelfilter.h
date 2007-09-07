//
// C++ Interface: voxelfilter
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MUVOXELFILTER_H
#define MUVOXELFILTER_H

#include <vector>
#include <map>

#include "../elements/elements.h"
#include "../delaunay_3d.h"

namespace Mu {

/**
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
class VoxelFilter{
protected:
	std::vector<Point *> points ;
	std::vector<DelaunayTetrahedron *> elems;
	
public:
	VoxelFilter();

	virtual ~VoxelFilter();
	
	void read(const char * filename) ;
	
	std::map<int,LinearForm *> behaviourMap ;
	
	std::vector<Point *> & getPoints() ;
	std::vector<DelaunayTetrahedron *> & getElements() ;
	
	const std::vector<Point *> & getPoints() const ;
	const std::vector<DelaunayTetrahedron *> & getElements() const ;

};

}

#endif //MUVOXELFILTER_H
