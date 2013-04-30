
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef __TIME_ENR_INCLUSION_H__
#define __TIME_ENR_INCLUSION_H__

#include "features.h"
#include "enrichmentInclusion.h"

namespace Mu
{

class TimeDependentEnrichmentInclusion : public EnrichmentInclusion
{
protected:
	Function radius_t ;
	
public:
	TimeDependentEnrichmentInclusion(Feature * father, Function & r, double x, double y) ;
	TimeDependentEnrichmentInclusion( Function & r, double x, double y) ;
	virtual ~TimeDependentEnrichmentInclusion() ;
  	
	virtual void enrich(size_t & , Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree) ;
	
	virtual void print() const
	{
		std::cout << "I am a time-dependent enriched inclusion" << std::endl ;
	}
	
	virtual void step(double dt, std::valarray<double> *, const  Mu::Mesh <Mu::DelaunayTriangle, Mu::DelaunayTreeItem > * dtree);

	void update(Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree) ;

} ;



  

}

#endif
