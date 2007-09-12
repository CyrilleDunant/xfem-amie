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
#ifndef MUVOXELPOREFILTER_H
#define MUVOXELPOREFILTER_H

#include <vector>
#include <map>

#include "../elements/elements.h"
#include "../delaunay_3d.h"

namespace Mu {

/**
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
class VoxelPoreFilter{
protected:
	std::vector<Point *> points ;
	std::vector<DelaunayTetrahedron *> elems;
	
public:
	VoxelPoreFilter();

	LinearForm * behaviour ;
	LinearForm * behaviourAlt ;
	virtual ~VoxelPoreFilter();
	
	void read(const char * filename) ;
	
	int poreIndex ;
	
	std::vector<Point *> & getPoints() ;
	std::vector<DelaunayTetrahedron *> & getElements() ;
	
	const std::vector<Point *> & getPoints() const ;
	const std::vector<DelaunayTetrahedron *> & getElements() const ;

};

}

#endif //MUVOXELFILTER_H
