// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution

#include "growingExpansiveZone.h"

using namespace Mu ;

GrowingExpansiveZone::GrowingExpansiveZone(Feature *father, double r_init, double x, double y, const Matrix & tensor, Vector def, Function r) : ExpansiveZone(father,r_init,x,y,tensor,def), growth(r)
{
	changed = (r_init >= POINT_TOLERANCE) ;
	time_pos = 0. ;
}

GrowingExpansiveZone::~GrowingExpansiveZone() {} ;

void GrowingExpansiveZone::enrich(size_t & counter, Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree)
{
	if(this->getGeometry()->getRadius() >= POINT_TOLERANCE)
	{
		ExpansiveZone::enrich(counter,dtree) ;
	}
}


void GrowingExpansiveZone::step(double dt, Vector *, const Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree)
{
	double r = this->getGeometry()->getRadius() ;

	double previous_r = VirtualMachine().eval(growth, time_pos) ;
	double dr = VirtualMachine().eval(growth,time_pos+dt) - previous_r;

	if(std::abs(dr) < POINT_TOLERANCE)
		changed = false ;
	else
	{
		this->setRadius(r+dr) ;
		changed = (r+dr >= POINT_TOLERANCE) ;	
	}

	time_pos += dt ;

}
	


