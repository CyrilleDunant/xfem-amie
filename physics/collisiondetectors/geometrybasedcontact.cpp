
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
#include "../contactmodels/linearcontactforce.h"
namespace Amie {

GeometryBasedContact::GeometryBasedContact(Geometry *geo, double delta) : geo(geo), delta(delta)
{

}


GeometryBasedContact::~GeometryBasedContact()
{
}

double GeometryBasedContact::grade(ElementState &s)
{
    VirtualMachine vm ;
    int dim = s.getParent()->spaceDimensions() ;
//     Vector deltaForce(0., dim*s.getParent()->getBoundingPoints().size()) ;
    Vector disp(dim) ;
    
    double num = -1 ;
    int count = 0 ;

    for(size_t i = 0 ; i < s.getParent()->getBoundingPoints().size() ; i++)
    {
//         s.getField(DISPLACEMENT_FIELD,s.getParent()->getBoundingPoint(i), disp,false,  &vm);
        
        Point test(s.getParent()->getBoundingPoint(i) + Vector({s.getDisplacements()[i*dim], s.getDisplacements()[i*dim+1]})) ;
        Point base(test) ;
        geo->project(&test);

        double dx = test.x- base.x ;
        double dy = test.y- base.y ;
        double dz = test.z- base.z ;
        double nom = sqrt(dx*dx+dy*dy+dz*dz) ;
        
//         std::cout << base.x << " " << base.y << "  " << nom << std::endl ;
        
        if(geo->in(base) && nom < s.getParent()->getRadius())
        {
            num = std::max(nom, num) ;
            count++ ;
        }
        else
        { 
            num = std::max(-nom, num) ;
        }
    }
    
    if(count)
        num *= count ;
    
    return num/delta;
  
}

FractureCriterion * GeometryBasedContact::getCopy() const
{
    return new GeometryBasedContact(geo, delta) ;
}


}
