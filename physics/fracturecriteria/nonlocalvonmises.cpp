
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
    Vector str( getSmoothedField( REAL_STRESS_FIELD, s ) ) ;
    double tr = (str.size()==3)?(str[0]+str[1]):(str[0]+str[1]+str[2]) ;
    
    if(str.size()==3)
    {
      for(size_t i = 0 ; i < 2 ; i++)
	str[i] -= tr/3. ;
      
    }
    else
    {
      for(size_t i = 0 ; i < 3 ; i++)
	str[i] -= tr/3. ;

    }

//     double vm = (s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)? sqrt(str[0]*str[0]+str[1]*str[1]-str[0]*str[1]) : sqrt(str[0]*str[0]+str[1]*str[1]+str[2]*str[2]-str[0]*str[1]-str[0]*str[2]-str[1]*str[2]);
    double vm = sqrt(3./2.*(str*str).sum()) ;
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
