
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../main.h"
#include "../../features/features.h"
#include "../../features/sample.h"
#include "../../features/growingExpansiveZone.h"
#include "../../physics/viscoelasticity.h"
#include "../../physics/viscoelasticity_and_imposed_deformation.h"
#include "../../utilities/writer/triangle_writer.h"
#include "../../utilities/parser.h"


#include <fstream>
#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h>


using namespace Amie ;

int main( int argc, char *argv[] )
{
	CommandLineParser parser("Test a single space-time XFEM inclusion") ;
	parser.addFlag("--renew-base", "renew the base of results") ;
	parser.addString("--output-directory","../examples/test/","directory where the results are stored", "-D") ;
	parser.parseCommandLine(argc, argv) ;
	bool renew = parser.getFlag("--renew-base") ;
	std::string outdir = parser.getString("--output-directory") ;


    Matrix C = Tensor::cauchyGreen( 10e9, 0.2, true, SPACE_TWO_DIMENSIONAL, PLANE_STRESS ) ;
    Vector v(3) ; v[0] = 0.01 ; v[1] = 0.01 ;

    Sample box(0.1,0.1,0,0) ;
    box.setBehaviour( new Viscoelasticity(PURE_ELASTICITY, C) ) ;

    FeatureTree F(&box) ;
    F.setOrder(LINEAR_TIME_LINEAR) ;
    Function radius(0.02) ; radius += "t 0.01 *" ;
    GrowingExpansiveZone * exp = new GrowingExpansiveZone( &box, radius, 0,0, new ViscoelasticityAndImposedDeformation( PURE_ELASTICITY, C*2., v ) ) ;
    F.addFeature(&box, exp) ;

    F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, LEFT_AFTER ) ) ;
    F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM_AFTER ) ) ;
    F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_STRESS_ETA, TOP_AFTER, 1e6 ) ) ;

    F.setSamplingNumber(4) ;
    F.setDeltaTime(1) ;
    F.setSolverPrecision(1e-16) ;

    std::ofstream out ;
    if(renew)
        out.open(outdir+"/test_stxfem_base", std::ios::out) ;
    else
        out.open(outdir+"/test_stxfem_current", std::ios::out) ;

    F.step() ;
    Vector str = F.getAverageField( STRAIN_FIELD ) ;
    out << F.getCurrentTime() << "\t" << exp->radiusAtTime(Point(0,0,0,F.getCurrentTime())) << "\t" << str[0]*1e3 << "\t" << str[1]*1e3 << "\t" << str[2]*1e3 << std::endl ;

    F.step() ;
    str = F.getAverageField( STRAIN_FIELD ) ;
    out << F.getCurrentTime() << "\t" << exp->radiusAtTime(Point(0,0,0,F.getCurrentTime())) << "\t" << str[0]*1e3 << "\t" << str[1]*1e3 << "\t" << str[2]*1e3 << std::endl ;

    return 0 ;
}
