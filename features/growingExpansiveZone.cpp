// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution

#include "growingExpansiveZone.h"

using namespace Mu ;

GrowingExpansiveZone::GrowingExpansiveZone(Feature* father, Function & g, double x, double y, const ViscoelasticityAndImposedDeformation* i) : ExpansiveZone(father,VirtualMachine().eval(g, Point(0,0,0,0)),x,y,i->getTensor(Point(0,0,0,0)),i->getImposedStrain(Point(0,0,0,0))), growth(g)
{
	changed = (radius >= POINT_TOLERANCE_2D) ;
}

GrowingExpansiveZone::~GrowingExpansiveZone() {} ;

void GrowingExpansiveZone::enrich(size_t & counter, Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree)
{
	if(getRadius() >= POINT_TOLERANCE_2D)
	{
		ExpansiveZone::enrich(counter,dtree) ;
	}
}


void GrowingExpansiveZone::step(double dt, Vector *, const Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree)
{
// 	double r = getRadius() ;
// 
// 	double previous_r = VirtualMachine().eval(growth, time_pos) ;
// 	double dr = VirtualMachine().eval(growth,time_pos+dt) - previous_r;
// 
// 	if(std::abs(dr) < POINT_TOLERANCE_2D)
// 		changed = false ;
// 	else
// 	{
// 		this->setRadius(r+dr) ;
// 		changed = (r+dr >= POINT_TOLERANCE_2D) ;	
// 	}
// 
// 	time_pos += dt ;

}
	


