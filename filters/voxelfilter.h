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

namespace Mu {

/**
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
class VoxelFilter{
protected:
	std::vector<Point *> points ;
	std::vector<HexahedralElement *> elems;
	
public:
	VoxelFilter();

	~VoxelFilter();
	
	void read(const char * filename) ;
	
	std::map<int,LinearForm *> behaviourMap ;
	
	std::vector<Point *> & getPoints() ;
	std::vector<HexahedralElement *> & getElements() ;
	
	const std::vector<Point *> & getPoints() const ;
	const std::vector<HexahedralElement *> & getElements() const ;

};

}

#endif
