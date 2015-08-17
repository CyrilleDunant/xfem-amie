
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef __GROWING_EXPANSIVE_ZONE_H__
#define __GROWING_EXPANSIVE_ZONE_H__

#include "timeDependentEnrichmentInclusion.h"
#include "../polynomial/vm_function_base.h"
#include "../physics/viscoelasticity_and_imposed_deformation.h"

namespace Amie
{

class GrowingExpansiveZone :  public TimeDependentEnrichmentInclusion
{
    std::set<DelaunayTriangle *> bimateralInterfaced ;
    std::set<DelaunayTriangle *> expansive ;
    Form * imp ;
    std::map<Point *, std::vector<double> > pointsAndValues ;

public:

    GrowingExpansiveZone(Feature *father, const Function & g, double x, double y, Form * i, double s = -10) ;
    GrowingExpansiveZone(Feature *father, const Function & g, double x, double y) ;
    virtual ~GrowingExpansiveZone() ;

    virtual void print() const
    {
        std::cout << "I am a growing expansive zone" << std::endl ;
    }

    virtual void enrich(size_t & , Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree) ;

    virtual bool moved() const {
        return changed ;
    } ;

    virtual void setBehaviour(Form * const b) { imp = b ; Feature::setBehaviour(b) ; }

    virtual void step(double dt, std::valarray<double> *, Mesh <Amie::DelaunayTriangle, Amie::DelaunayTreeItem > * dtree);

public:
    GEO_DERIVED_OBJECT(TimeDependentCircle) ;

} ;

}

#endif
