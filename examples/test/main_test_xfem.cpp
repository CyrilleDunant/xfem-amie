
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../main.h"
#include "../../features/features.h"
#include "../../features/sample.h"
#include "../../features/expansiveZone.h"
#include "../../physics/stiffness.h"
#include "../../physics/stiffness_with_imposed_deformation.h"


#include <fstream>
#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h>


using namespace Amie ;

int main( int argc, char *argv[] )
{
    Matrix C = Tensor::cauchyGreen( 10e9, 0.2, true, SPACE_TWO_DIMENSIONAL, PLANE_STRAIN ) ;
    Vector v(3) ; v[0] = 0.01 ; v[1] = 0.01 ;

    Sample box(0.1,0.1,0,0) ;
    box.setBehaviour( new Stiffness(C) ) ;

    FeatureTree F(&box) ;
    ExpansiveZone * exp = new ExpansiveZone( &box, 0.02, 0,0, new StiffnessWithImposedDeformation( C*2, v ) ) ;
    F.addFeature(&box, exp) ;

    F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, LEFT ) ) ;
    F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM ) ) ;

    F.setSamplingNumber(4) ;
    F.setDeltaTime(1) ;

    std::ofstream out ;
    out.open("../examples/test/test_xfem_current", std::ios::out) ;

    F.step() ;
    Vector str = F.getAverageField( STRAIN_FIELD ) ;
    out << F.getCurrentTime() << "\t" << str[0] << "\t" << str[1] << "\t" << str[2] << std::endl ;

    exp->setRadius(0.03) ;
    F.step() ;
    str = F.getAverageField( STRAIN_FIELD ) ;
    out << F.getCurrentTime() << "\t" << str[0] << "\t" << str[1] << "\t" << str[2] << std::endl ;

    return 0 ;
}
