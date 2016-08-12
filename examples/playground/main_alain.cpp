// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../main.h"
#include "../../features/features.h"
#include "../../physics/viscoelasticity.h"
#include "../../physics/viscoelasticity_and_fracture.h"
#include "../../physics/fracturecriteria/spacetimemultilinearsofteningfracturecriterion.h"
#include "../../physics/fracturecriteria/spacetimemultisurfacefracturecriterion.h"
#include "../../physics/fracturecriteria/spacetimelimitsurfacefracturecriterion.h"
#include "../../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h"
#include "../../physics/damagemodels/spacetimefiberbasedfixedcrack.h"
#include "../../physics/damagemodels/spacetimebifurcation.h"
#include "../../physics/material_laws/mechanical_material_laws.h"
#include "../../physics/materials/paste_behaviour.h"
#include "../../utilities/writer/triangle_writer.h"
#include "../../features/sample.h"
#include "../../utilities/parser/command_line_parser.h"

#include <fstream>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include <cmath>
#include <typeinfo>
#include <limits>
#include <sys/time.h>
#define DEBUG


using namespace Amie ;

int main(int argc, char *argv[])
{
        CommandLineParser parser("Test a fixed crack damage behaviour on two elements") ;
        parser.addFlag("--renew-base", "renew the base of results") ;
        parser.addString("--output-directory","../examples/test/","directory where the results are stored", "-D") ;
        parser.parseCommandLine(argc, argv) ;
        bool renew = parser.getFlag("--renew-base") ;
        std::string outdir = parser.getString("--output-directory") ;

        std::ofstream out ;
        if(renew)
                out.open(outdir+"/test_spacetime_fixed_crack_base", std::ios::out) ;
        else
                out.open(outdir+"/test_spacetime_fixed_crack_current", std::ios::out) ;

        Sample rect(nullptr, 0.01,0.01,0,0) ;

        LogarithmicCreepWithExternalParameters * toto = new LogarithmicCreepWithExternalParameters("young_modulus = 10e9, poisson_ratio = 0.2, creep_modulus = 40e9, creep_characteristic_time = 1, imposed_deformation = 0") ;
        LogarithmicCreepWithExternalParameters * tata = new LogarithmicCreepWithExternalParameters("young_modulus = 60e9, poisson_ratio = 0.2, imposed_deformation = 0.001") ;

	Inclusion * inc = new Inclusion( 0.003,0,0 ) ;
	inc->setBehaviour(tata) ;

        rect.setBehaviour( toto ) ;

        FeatureTree f(&rect) ;
        f.setSamplingNumber(4) ;
        f.setDeltaTime(1) ;
        f.setMinDeltaTime(1e-9) ;
	f.addFeature( &rect, inc ) ;
        f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, LEFT_AFTER) ) ;
        f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM_AFTER ) ) ;

	for(size_t i = 0 ; i < 10 ; i++)
	{
		f.setDeltaTime( i+1) ;
	        f.step() ;
		Vector strain = f.getAverageField( STRAIN_FIELD, 1. ) ;
		std::cout << strain[0] << "\t" << strain[1] << "\t"  << strain[2] << std::endl ;
	}


        return 0 ;
}

