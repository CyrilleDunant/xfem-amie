
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
#include "nonlocalvonmises.h"
#include "../../mesher/delaunay.h"
#include "../damagemodels/damagemodel.h"
namespace Amie {

NonLocalVonMises::NonLocalVonMises(double thresh, double radius) : threshold(std::abs(thresh))
{
    setMaterialCharacteristicRadius(radius);
    met = false ;
//     smoothingType = GAUSSIAN_NONCOMPACT ;
}


NonLocalVonMises::~NonLocalVonMises()
{
}

double NonLocalVonMises::grade(ElementState &s)
{
    met = false ;
//     Vector str(2) ;
//     s.getField(PRINCIPAL_REAL_STRESS_FIELD, Point(), str, true);
    Vector str( getSmoothedField( PRINCIPAL_REAL_STRESS_FIELD, s ) ) ;

    double vm = (s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)? sqrt(str[0]*str[0]+str[1]*str[1]-str[0]*str[1]) : sqrt(str[0]*str[0]+str[1]*str[1]+str[2]*str[2]-str[0]*str[1]-str[0]*str[2]-str[1]*str[2]);
//     std::cout << "\n" << vm << std::endl ;
    if( vm >= threshold )
        met = true ;


    return -1. + std::abs( vm / threshold );

}

FractureCriterion * NonLocalVonMises::getCopy() const
{
    NonLocalVonMises * ret = new NonLocalVonMises(threshold, getMaterialCharacteristicRadius()) ;
    ret->copyEssentialParameters( this ) ;
    return ret ;
}

}
