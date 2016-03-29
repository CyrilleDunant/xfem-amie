
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
#include "boundedvonmises.h"
namespace Amie {

BoundedVonMises::BoundedVonMises(double thres, double damageThreshold) : threshold(thres), damageThreshold(damageThreshold), dmodel(nullptr)
{

}


BoundedVonMises::~BoundedVonMises()
{
}

double BoundedVonMises::grade(ElementState &s)
{
    dmodel = s.getParent()->getBehaviour()->getDamageModel() ;
    if(dmodel && dmodel->getState().max() > damageThreshold)
        return -1. ;

    Vector v(0.,1) ;
    s.getAverageField( VON_MISES_REAL_STRESS_FIELD, v, nullptr, 0) ;
    double maxStress = v[0] ;

    if(maxStress > threshold )
    {
        return 1. - std::abs(threshold/maxStress) ;
    }
    else
    {
        return -1.+ std::abs(maxStress/threshold);
    }
}

FractureCriterion * BoundedVonMises::getCopy() const
{
    BoundedVonMises * ret = new BoundedVonMises(threshold, damageThreshold) ;
    ret->copyEssentialParameters( this ) ;
    return ret ;
}



}
