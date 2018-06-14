
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

NonLocalDeviatoricVonMises::NonLocalDeviatoricVonMises(double thresh, double radius) : threshold(std::abs(thresh))
{
    setMaterialCharacteristicRadius(radius);
//     smoothingType = GAUSSIAN_NONCOMPACT ;
//         smoothingType = GAUSSIAN_NONCOMPACT ;

}


NonLocalDeviatoricVonMises::~NonLocalDeviatoricVonMises()
{
}

double NonLocalDeviatoricVonMises::grade(ElementState &s)
{
    Vector str = getSmoothedField( PRINCIPAL_REAL_STRESS_FIELD, s )  ;
    
//     int dim = s.getParent()->spaceDimensions() ;
//     VirtualMachine vim ;
//     Vector disp(dim) ;
//     s.getField(PRINCIPAL_REAL_STRESS_FIELD,s.getParent()->getCenter(), disp,false,  &vim);
    
    double tr = str.sum() ;
    
    str -= tr/str.size() ;
    
    double vm = 0 ;
    if(str.size() == 2)
        vm = sqrt(str[0]*str[0]-str[0]*str[1]+str[1]*str[1]) ;
    else
        vm = sqrt(0.5*((str[0]-str[1])*(str[0]-str[1])+(str[1]-str[2])*(str[1]-str[2])+(str[2]-str[0])*(str[2]-str[0]))) ;

//     std::cout << disp[0] << "  " << disp[1] << "  "<< std::endl ;

    double v = std::abs( vm / threshold )-1;
    return v;

}

FractureCriterion * NonLocalDeviatoricVonMises::getCopy() const
{
    NonLocalDeviatoricVonMises * ret = new NonLocalDeviatoricVonMises(threshold, getMaterialCharacteristicRadius()) ;
    ret->copyEssentialParameters( this ) ;
    return ret ;
}

}
