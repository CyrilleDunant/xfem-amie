// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution

#include "vibratingcircularmembrane.h"
#include "../physics/wave.h"

using namespace Mu ;

VibratingMembrane::VibratingMembrane(Feature *father, double radius, double x, double y, Matrix tensor) : EnrichmentInclusion(father, radius, x, y), cgTensor(tensor)
{
	
}

VibratingMembrane::~VibratingMembrane() {}
	

void VibratingMembrane::reset() 
{
	cache.clear() ;
	updated = true ;
}

void VibratingMembrane::enrich(size_t & lastId, Mesh< DelaunayTriangle, DelaunayTreeItem >* dtree)
{
	/*
	if(cache.empty())
		cache = dtree->getConflictingElements(getPrimitive()) ;

	//we get a unique list of the points
	std::map<Point *, int> done ;

	for(size_t i = 0 ; i < cache.size() ; i++)
	{
		done[&cache[i]->getBoundingPoint(0)] = -1 ;
		done[&cache[i]->getBoundingPoint(1)] = -1 ;
		done[&cache[i]->getBoundingPoint(2)] = -1 ;
		done[&cache[i]->getBoundingPoint(3)] = -1 ;
		done[&cache[i]->getBoundingPoint(4)] = -1 ;
		done[&cache[i]->getBoundingPoint(5)] = -1 ;
	}
	
	std::valarray<Function> shapefunc = TriElement(LINEAR_TIME_LINEAR).getShapeFunctions() ;
	VirtualMachine vm ;
	std::vector<Point> hint ;
	hint.push_back(Point(1./3., 1./3.)) ;
	
	for(size_t i = 0 ; i < cache.size() ; i++)
	{
		Function x = cache[i]->getXTransform() ;
		Function y = cache[i]->getYTransform() ;
		Function x_ = x - getCenter().getX() ;
		Function y_ = y - getCenter().getY() ;
		Function theta = f_atan2 ( y_, x_ );
		Function r = f_sqrt ( ( x_^2 ) + ( y_^2 ) );
		
		double besselRoot01 = 2.40482555769577/getRadius();
		double besselRoot02 = 5.52007811028631/getRadius();
		double besselRoot11 = 3.83170597020751/getRadius();
		double besselRoot12 = 7.01558666981561/getRadius();
		
		Function f_01 = f_cyl_bessel_j(0, r*besselRoot01) ;
		Function f_02 = f_cyl_bessel_j(0, r*besselRoot02) ;
		Function f_11 = f_cyl_bessel_j(1, r*besselRoot11)*f_cos(theta) ;
		Function f_12 = f_cyl_bessel_j(1, r*besselRoot12)*f_cos(theta) ;
		Function f_11s = f_cyl_bessel_j(1, r*besselRoot11)*f_sin(theta) ;
		Function f_12s = f_cyl_bessel_j(1, r*besselRoot12)*f_sin(theta) ;
		for(size_t j = 0 ; j < 6 ; j++)
		{
			Point current = cache[i]->inLocalCoordinates(cache[i]->getBoundingPoint(j)) ; 
			Point * currentPointer = & cache[i]->getBoundingPoint(j) ;
			int id = cache[i]->getBoundingPoint(j).getId() ;
			Function f = shapefunc[j]* ( f_01 - vm.eval ( f_01, current) )  ;
			f.setIntegrationHint ( hint ) ;
			f.setPoint ( &cache[i]->getBoundingPoint(j) ) ;
			if(done[currentPointer] == -1)
			{
				done[currentPointer] = lastId ;
				 lastId+=6 ;
			}
			f.setDofID ( done[cache[i]->first] ) ;
			cache[i]->setEnrichment (  f , static_cast<Circle *>(this) ) ;
			
			f = shapefunc[j]* ( f_02 - vm.eval ( f_02, current ) ) ;
			f.setIntegrationHint ( hint ) ;
			f.setPoint ( &cache[i]->getBoundingPoint(j) ) ;
			f.setDofID ( done[currentPointer]+1 ) ;
			cache[i]->setEnrichment (  f , static_cast<Circle *>(this) ) ;
			
			f = shapefunc[j]* ( f_11 - vm.eval ( f_11, current ) ) ;
			f.setIntegrationHint ( hint ) ;
			f.setPoint ( &cache[i]->getBoundingPoint(j) ) ;
			f.setDofID ( done[currentPointer]+2 ) ;
			cache[i]->setEnrichment (  f , static_cast<Circle *>(this) ) ;
			
			f = shapefunc[j]* ( f_12 - vm.eval ( f_12,current ) ) ;
			f.setIntegrationHint ( hint ) ;
			f.setPoint ( &cache[i]->getBoundingPoint(j) ) ;
			f.setDofID ( done[currentPointer]+3 ) ;
			cache[i]->setEnrichment (  f , static_cast<Circle *>(this) ) ;
			
			f = shapefunc[j]* ( f_11s - vm.eval ( f_11s, current ) ) ;
			f.setIntegrationHint ( hint ) ;
			f.setPoint ( &cache[i]->getBoundingPoint(j) ) ;
			f.setDofID ( done[currentPointer]+4 ) ;
			cache[i]->setEnrichment (  f , static_cast<Circle *>(this) ) ;
			
			f = shapefunc[j]* ( f_12s - vm.eval ( f_12s, current ) ) ;
			f.setIntegrationHint ( hint ) ;
			f.setPoint ( &cache[i]->getBoundingPoint(j) ) ;
			f.setDofID ( done[currentPointer]+5 ) ;
			cache[i]->setEnrichment (  f , static_cast<Circle *>(this) ) ;
		}
		
		if(in(cache[i]->getCenter()))
			cache[i]->setBehaviour(new Wave(cgTensor)) ;
	}
	*/
}
	

