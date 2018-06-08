
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
    
//     s.getField(DISPLACEMENT_FIELD,s.getParent()->getCenter(), disp,false,  &vm);
//     Point test(s.getParent()->getCenter() + disp) ;
//     Point base(test) ;
//     geo->project(&test);
// 
//     double dx = test.x- base.x ;
//     double dy = test.y- base.y ;
//     double dz = test.z- base.z ;
//     double cnom = sqrt(dx*dx+dy*dy+dz*dz) ;
//     bool cin = geo->in(base) ;
    
//     double inpoints = 0 ;
//     double overlapPoints = 0 ;
//     for(double i = -s.getParent()->getRadius() ; i < s.getParent()->getRadius() ; i += 0.05*s.getParent()->getRadius())
//     {
//         for(double j = -s.getParent()->getRadius() ; j < s.getParent()->getRadius() ; j += 0.05*s.getParent()->getRadius())
//         {
//             Point pt(i+s.getParent()->getCenter().getX(), j+s.getParent()->getCenter().getY()) ;
//             if(s.getParent()->in(pt))
//             {
//                 inpoints++ ;
//                 if(geo->in(pt))
//                     overlapPoints++ ;
//             }
//         }
//     }
//     
//     if(overlapPoints > 1)
//         return overlapPoints/inpoints ;
//     return -1 ;
    
    for(size_t i = 0 ; i < s.getParent()->getBoundingPoints().size() ; i++)
    {
        s.getField(DISPLACEMENT_FIELD,s.getParent()->getBoundingPoint(i), disp,false,  &vm);
        Point test(s.getParent()->getBoundingPoint(i) + disp) ;
        Point base(test) ;
        geo->project(&test);

        double dx = test.x- base.x ;
        double dy = test.y- base.y ;
        double dz = test.z- base.z ;
        double nom = sqrt(dx*dx+dy*dy+dz*dz) ;
        
//         std::cout << base.x << " " << base.y << "  " << nom << std::endl ;
        
        if(geo->in(base))
        {
            num = std::max(nom, num) ;
            count++ ;
        }
        else
        { 
            num = std::max(-nom, num) ;
        }
    }
    
    if(count == 1)
    {
//         if(cin && cnom < num)
            num = -1 ;
    }
    
    
    
    LinearContactForce * cf = dynamic_cast<LinearContactForce *>(s.getParent()->getBehaviour()->getContactModel()) ;
    if(cf && num < POINT_TOLERANCE)
    {
        
        if(std::abs(cf->forces+cf->deltaForce*cf->getState()[0]).max() > POINT_TOLERANCE)
            return -std::abs(num/delta-POINT_TOLERANCE) ;
    }
    
    return num/delta -POINT_TOLERANCE;
  
}

FractureCriterion * GeometryBasedContact::getCopy() const
{
    return new GeometryBasedContact(geo, delta) ;
}


}
