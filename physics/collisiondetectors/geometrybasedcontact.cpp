
//
// C++ Implementation: vonmises
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "geometrybasedcontact.h"
#include "../../mesher/delaunay.h"
#include "../damagemodels/damagemodel.h"
namespace Amie {

GeometryBasedContact::GeometryBasedContact(Geometry *geo) : geo(geo)
{

}


GeometryBasedContact::~GeometryBasedContact()
{
}

double GeometryBasedContact::grade(ElementState &s)
{
    double sc = -1 ;
    int dim = s.getParent()->spaceDimensions() ;
    VirtualMachine vm ;
    for(size_t i = 0 ; i < s.getParent()->getBoundingPoints().size() ; i++)
    {
        Vector disp(dim) ;
        s.getField(DISPLACEMENT_FIELD,s.getParent()->getBoundingPoint(i), disp,false,  &vm);
        Point test(s.getParent()->getBoundingPoint(i) + disp) ;
        Point base(test) ;
        geo->project(&test);
        double dx = test.x- base.x ;
        double dy = test.y- base.y ;
        double dz = test.z- base.z ;
        
        if(geo->in(base))
        {
            sc = std::max(sc, std::max(std::abs(dx), std::max(std::abs(dy),std::abs(dz)))) ;
        }
        else
        {
            sc = std::max(std::max(sc, -1.), -std::max(sc, std::max(std::abs(dx), std::max(std::abs(dy),std::abs(dz))))) ;
        }
    }
    
    return sc/geo->getRadius() ;

}

FractureCriterion * GeometryBasedContact::getCopy() const
{
    return new GeometryBasedContact(geo) ;
}


}
