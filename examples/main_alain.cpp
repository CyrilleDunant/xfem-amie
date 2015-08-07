
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "./main.h"
#include "./../features/features.h"
#include "./../features/sample.h"
#include "./../features/growingExpansiveZone.h"
#include "./../physics/viscoelasticity.h"
#include "./../physics/viscoelasticity_and_imposed_deformation.h"
#include "./../utilities/writer/triangle_writer.h"


#include <fstream>
#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h>


using namespace Amie ;

int main( int argc, char *argv[] )
{
    Matrix C = Tensor::cauchyGreen( 10e9, 0.2, true, SPACE_TWO_DIMENSIONAL, PLANE_STRESS ) ;
    Vector v(3) ; v[0] = 0.01 ; v[1] = 0.01 ;

    Sample box(0.1,0.1,0,0) ;
    box.setBehaviour( new Viscoelasticity(PURE_ELASTICITY, C) ) ;

    FeatureTree F(&box) ;
    F.setOrder(LINEAR_TIME_LINEAR) ;
    Function radius(0.02) ; //radius += "t 0.01 *" ;
    GrowingExpansiveZone * exp = new GrowingExpansiveZone( &box, radius, 0,0, new Viscoelasticity( PURE_ELASTICITY, C ) ) ;
    F.addFeature(&box, exp) ;

    F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, LEFT_AFTER ) ) ;
    F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM_AFTER ) ) ;
    F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_STRESS_ETA, TOP_AFTER, 1e6 ) ) ;

    F.setSamplingNumber(4) ;
    F.setDeltaTime(1) ;

    std::ofstream out ;
    out.open("../examples/test/test_stxfem_base", std::ios::out) ;

    F.step() ;
    Vector str = F.getAverageField( STRAIN_FIELD, -1, 1 ) ;
    std::cout << F.getCurrentTime() << "\t" << exp->radiusAtTime(Point(0,0,0,F.getCurrentTime())) << "\t" << str[0] << "\t" << str[1] << "\t" << str[2] << std::endl ;

    F.step() ;
    str = F.getAverageField( STRAIN_FIELD, -1, 1 ) ;
    std::cout << F.getCurrentTime() << "\t" << exp->radiusAtTime(Point(0,0,0,F.getCurrentTime())) << "\t" << str[0] << "\t" << str[1] << "\t" << str[2] << std::endl ;

    F.step() ;
    str = F.getAverageField( STRAIN_FIELD, -1, 1 ) ;
    std::cout << F.getCurrentTime() << "\t" << exp->radiusAtTime(Point(0,0,0,F.getCurrentTime())) << "\t" << str[0] << "\t" << str[1] << "\t" << str[2] << std::endl ;
    
    TriangleWriter trg( "fields", &F, 1) ;
    trg.getField( TWFT_STIFFNESS ) ;
    trg.getField( STRAIN_FIELD ) ;
    trg.getField( REAL_STRESS_FIELD ) ;
    trg.write() ;

/*    F.getAssembly()->print() ;
    Vector d = F.getDisplacements() ;
    for(size_t i = 0 ; i < d.size() ; i++)
        std::cout << d[i] << std::endl ;*/

    return 0 ;
}
