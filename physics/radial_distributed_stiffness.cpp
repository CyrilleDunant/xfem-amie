//
// C++ Implementation: RadialDistributedStiffness
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
//

#include "radial_distributed_stiffness.h"
#include "physics_base.h"
#include "stiffness.h"


using namespace Amie ;

RadialDistributedStiffness::RadialDistributedStiffness(std::vector<std::pair<double, Matrix> > rig) : LinearForm(rig[0].second, false, false, rig[0].second.numRows()/3+1)
{
    v.push_back(XI);
    v.push_back(ETA);
    if(param.size() > 9)
        v.push_back(ZETA);

    stiff = rig ;
    angle = 0 ;
}


RadialDistributedStiffness::~RadialDistributedStiffness() { } 

void RadialDistributedStiffness::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{

    vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v, ret) ;
}

bool RadialDistributedStiffness::fractured() const
{
    return false ;
}

void RadialDistributedStiffness::setAngle(double a)
{
    angle = a ;
    while(angle > 2 * M_PI)
        angle = angle - 2*M_PI ;
}

Form * RadialDistributedStiffness::getCopy() const
{
    size_t i = 0 ;
    while(angle < stiff[i].first && i < stiff.size())
        i++ ;

    Stiffness* copy = new Stiffness( stiff[i].second ) ;

    return copy ;
}

RadialInclusion::RadialInclusion(Feature * f, double r, Point c) : Circle(r,c), Inclusion(f,r,c)
{
}

RadialInclusion::RadialInclusion(Feature * f, double r, double x, double y) : Circle(r,x,y), Inclusion(f,r,x,y)
{
}

RadialInclusion::RadialInclusion(double r, Point c) : Circle(r,c), Inclusion(r,c)
{
}

RadialInclusion::RadialInclusion(double r, double x, double y) : Circle(r,x,y), Inclusion(r,x,y)
{
}

Form * RadialInclusion::getBehaviour( const Point & p ) const 
{
    static_cast<RadialDistributedStiffness *>(behaviour)->setAngle(p.angle()) ;
    return behaviour->getCopy() ;
}


