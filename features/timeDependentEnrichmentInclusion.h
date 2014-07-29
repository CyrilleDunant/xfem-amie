
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef __TIME_ENR_INCLUSION_H__
#define __TIME_ENR_INCLUSION_H__

#include "features.h"
#include "enrichmentInclusion.h"
#include "../physics/viscoelasticity_and_imposed_deformation.h"
#include "../geometry/space_time_geometry_2D.h"

namespace Amie
{

class TimeDependentEnrichmentInclusion : public EnrichmentInclusion, public TimeDependentCircle
{
protected:  
	std::map<const Point *, int> dofIdPrev ;
	std::map<const Point *, int> dofIdCurrent ;
	std::set<DelaunayTriangle *> enrichedElem;
  
public:
	TimeDependentEnrichmentInclusion(Feature * father, Function & r, double x, double y) ;
	TimeDependentEnrichmentInclusion( Function & r, double x, double y) ;
	virtual ~TimeDependentEnrichmentInclusion() ;
  	
	virtual void enrich(size_t & , Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree) ;
	
	virtual void print() const
	{
		std::cout << "I am a time-dependent enriched inclusion" << std::endl ;
	}
	
	virtual void step(double dt, std::valarray<double> *,  Amie::Mesh <Amie::DelaunayTriangle, Amie::DelaunayTreeItem > * dtree);

	void update(Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree) ;

public:
	GEO_DERIVED_OBJECT(TimeDependentCircle) ;
	
} ;

class TimeDependentHomogenisingInclusion : public EnrichmentInclusion, public TimeDependentCircle
{
protected:  
	std::map<DelaunayTriangle *, Form *> homogeneised;
	ViscoelasticityAndImposedDeformation * imposed ;
	std::set<DelaunayTriangle *> enrichedElem;
	std::set<DelaunayTriangle *> zoneElem;
  
public:
	TimeDependentHomogenisingInclusion(Feature * father, Function & r, double x, double y, ViscoelasticityAndImposedDeformation * imp) ;
	TimeDependentHomogenisingInclusion( Function & r, double x, double y, ViscoelasticityAndImposedDeformation * imp) ;
	virtual ~TimeDependentHomogenisingInclusion() ;
  	
	virtual void enrich(size_t & , Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree) ;
	
	virtual void print() const
	{
		std::cout << "I am a time-dependent enriched inclusion" << std::endl ;
	}
	
	virtual void step(double dt, std::valarray<double> *,  Amie::Mesh <Amie::DelaunayTriangle, Amie::DelaunayTreeItem > * dtree);

	void update(Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree) ;

public:
	GEO_DERIVED_OBJECT(TimeDependentCircle) ;
	

} ;


  

}

#endif
