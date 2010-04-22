
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef __GROWING_EXPANSIVE_ZONE_H__
#define __GROWING_EXPANSIVE_ZONE_H__

#include "expansiveZone.h"
#include "../polynomial/vm_function_base.h"

namespace Mu
{

class GrowingExpansiveZone :  public ExpansiveZone
{
	bool changed ;
	double time_pos ;
	Function growth ;
public:

	GrowingExpansiveZone(Feature *father, double r_init, double x, double y, const Matrix & cgTensor, Vector deformation, Function r) ;
	virtual ~GrowingExpansiveZone() ;
	
	virtual void print() const
	{
		std::cout << "I am a growing expansive zone" << std::endl ;
	}
	
	virtual void enrich(size_t & , Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree) ;

	virtual void step(double dt, Vector *, const Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree) ;

	virtual bool moved() const { return changed ; } ;
	
public:
	GEO_DERIVED_OBJECT(Circle) ;
	
} ;

}

#endif
