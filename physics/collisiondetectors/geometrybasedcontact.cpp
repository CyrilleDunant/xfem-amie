
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
        VirtualMachine vm ;
        int dim = s.getParent()->spaceDimensions() ;
        Vector deltaForce(0.,s.getParent()->getBoundingPoints().size()*dim) ;
        Vector tDeltaForce(0.,s.getParent()->getBoundingPoints().size()*dim) ;
        
        Vector disp(dim) ;
        s.getField(DISPLACEMENT_FIELD,s.getParent()->getCenter(), disp,false,  &vm);
        Point test(s.getParent()->getCenter() + disp) ;
        Point base(test) ;
        geo->project(&test);
        double mdx = test.x- base.x ;
        double mdy = test.y- base.y ;
        double mdz = test.z- base.z ;
        int count = 0 ;
        for(size_t i = 0 ; i < s.getParent()->getBoundingPoints().size() ; i++)
        {
            s.getField(DISPLACEMENT_FIELD,s.getParent()->getBoundingPoint(i), disp,false,  &vm);
            Point test(s.getParent()->getBoundingPoint(i) + disp) ;
            Point base(test) ;
            geo->project(&test);
            
            double dx = test.x- base.x ;
            double dy = test.y- base.y ;
            double dz = test.z- base.z ;
            
            double num = sqrt((dx*dx+dy*dy+dz*dz)*(mdx*mdx+mdy*mdy+mdz*mdz)) ;
            if(num <  1e-12)
                continue ;
            double normalFactor = std::abs((dx*mdx+dy*mdy+dz*mdz)/num);
            double tangentFactor = 1.-normalFactor ;
            if(geo->in(base))
            {

                count++ ;

                deltaForce[i*dim] = std::abs(dx*normalFactor) ;
                deltaForce[i*dim+1] = std::abs(dy*normalFactor) ;
                if(dim == 3)
                {
                    deltaForce[i*dim+2] = std::abs(dz*normalFactor) ;
                }
                
                tDeltaForce[i*dim] = std::abs(disp[0]*tangentFactor) ;
                tDeltaForce[i*dim+1] = std::abs(disp[1]*tangentFactor) ;
                if(dim == 3)
                    tDeltaForce[i*dim+2] = std::abs(disp[2]*tangentFactor) ;
                
            }
            else
            { 
                deltaForce[i*dim] = -std::abs(dx*normalFactor) ;
                deltaForce[i*dim+1] = -std::abs(dy*normalFactor) ;
                if(dim == 3)
                {
                    deltaForce[i*dim+2] = -std::abs(dz*normalFactor) ;
                }
                
                tDeltaForce[i*dim] = -std::abs(disp[0]*tangentFactor) ;
                tDeltaForce[i*dim+1] = -std::abs(disp[1]*tangentFactor) ;
                if(dim == 3)
                    tDeltaForce[i*dim+2] = -std::abs(disp[2]*tangentFactor) ;

            }
        }
        
        if(count == 1)
        {
            deltaForce = -1 ;
            tDeltaForce = -1 ;
        }
//         std::cout << deltaForce.max() << std::endl ;
        return std::max(deltaForce.max(), tDeltaForce.max())*s.getParent()->getRadius()*1e3 ;
    
   /* 
    
    
    double sc = -1 ;
    int dim = s.getParent()->spaceDimensions() ;
    VirtualMachine vm ;
    
    Vector disp(dim) ;
    s.getField(DISPLACEMENT_FIELD,s.getParent()->getCenter(), disp,false,  &vm);
    Point test(s.getParent()->getCenter() + disp) ;
    Point base(test) ;
    geo->project(&test);
    double dx = test.x- base.x ;
    double dy = test.y- base.y ;
    double dz = test.z- base.z ;
    double dnorm = sqrt(dx*dx+dy*dy+dz*dz)/s.getParent()->getRadius() ;
    int count = 0 ;
    for(size_t i = 0 ; i < s.getParent()->getBoundingPoints().size() ; i++)
    {
        s.getField(DISPLACEMENT_FIELD,s.getParent()->getBoundingPoint(i), disp,false,  &vm);
        Point test(s.getParent()->getBoundingPoint(i) + disp) ;
        Point base(test) ;
        geo->project(&test);
        double dx = test.x- base.x ;
        double dy = test.y- base.y ;
        double dz = test.z- base.z ;
        
        if(geo->in(base))
        {
            count++ ;
            if(sc > 0)
                sc = std::min(sqrt(dx*dx+dy*dy+dz*dz), sc) ;
            else
                sc = sqrt(dx*dx+dy*dy+dz*dz) ;
        }
        else
        {
            sc = std::max(std::max(sc, -1.), -sqrt(dx*dx+dy*dy+dz*dz)) ;
        }
    }
    if(count == 1)
        return -1 ;
//     if(dnorm  > 1e-6)
//         return sco/dnorm ;
    return sc;*/

}

FractureCriterion * GeometryBasedContact::getCopy() const
{
    return new GeometryBasedContact(geo) ;
}


}
