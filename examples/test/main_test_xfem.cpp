
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
#include "../../utilities/parser/command_line_parser.h"

#include <fstream>
#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h>


using namespace Amie ;

int main( int argc, char *argv[] )
{
	CommandLineParser parser("Test a single XFEM inclusion") ;
	parser.addFlag("--renew-base", "renew the base of results") ;
	parser.addString("--output-directory","../examples/test/","directory where the results are stored", "-D") ;
	parser.parseCommandLine(argc, argv) ;
	bool renew = parser.getFlag("--renew-base") ;
	std::string outdir = parser.getString("--output-directory") ;


    Matrix C = Tensor::cauchyGreen( 10e9, 0.2, SPACE_TWO_DIMENSIONAL, PLANE_STRAIN, YOUNG_POISSON ) ;
    Vector v(3) ; v[0] = 0.01 ; v[1] = 0.01 ;

    Sample box(0.1,0.1,0,0) ;
    box.setBehaviour( new Stiffness(C) ) ;

    FeatureTree F(&box) ;
    ExpansiveZone * exp = new ExpansiveZone( &box, 0.02, 0,0, new StiffnessWithImposedDeformation( C*2, v ) ) ;
//    Inclusion * iexp = new Inclusion( &box, 0.02, 0,0) ; iexp->setBehaviour(new StiffnessWithImposedDeformation( C*2, v ) ) ;
    F.addFeature(&box, exp) ;

    F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, BOTTOM_LEFT ) ) ;
    F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM ) ) ;

    F.setSamplingNumber(4) ;
    F.setDeltaTime(1) ;

    std::ofstream out ;
    if(renew)
        out.open(outdir+"/test_xfem_base", std::ios::out) ;
    else
        out.open(outdir+"/test_xfem_current", std::ios::out) ;

    F.step() ;
    F.step() ;
    Vector str = F.getAverageField( TOTAL_STRAIN_FIELD )*1e3 ;
    Vector sig = F.getAverageField( REAL_STRESS_FIELD )/1e6 ;
    out << F.getCurrentTime() << "\t" << str[0] << "\t" << str[1] << "\t" << str[2] << "\t" << sig[0] << "\t" << sig[1] << "\t" << sig[2] << std::endl ;

    exp->setRadius(0.03) ;
    F.step() ;
    F.step() ;
    str = F.getAverageField( TOTAL_STRAIN_FIELD )*1e3 ;
    sig = F.getAverageField( REAL_STRESS_FIELD )/1e6 ;
    out << F.getCurrentTime() << "\t" << str[0] << "\t" << str[1] << "\t" << str[2] << "\t" << sig[0] << "\t" << sig[1] << "\t" << sig[2] << std::endl ;

    exp->setRadius(0.5) ;
    F.step() ;
    F.step() ;
    str = F.getAverageField( TOTAL_STRAIN_FIELD )*1e3 ;
    sig = F.getAverageField( REAL_STRESS_FIELD )/1e6 ;
    out << F.getCurrentTime() << "\t" << str[0] << "\t" << str[1] << "\t" << str[2] << "\t" << sig[0] << "\t" << sig[1] << "\t" << sig[2] << std::endl ;

    return 0 ;
}
