// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../main.h"
#include "../../features/features.h"
#include "../../features/sample.h"
#include "../../features/polygonSample.h"
#include "../../physics/stiffness.h"
#include "../../physics/logarithmic_creep_with_external_parameters.h"
#include "../../physics/material_laws/mechanical_material_laws.h"
#include "../../physics/material_laws/material_laws.h"
#include "../../utilities/parser.h"
#include "../../utilities/writer/triangle_writer.h"
#include "../../geometry/level_set.h" 
#include "../../utilities/mineral.h" 


#include <fstream>
#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h>
#include <sys/time.h>


using namespace Amie ;


int main( int argc, char *argv[] )
{

        Sample rect(nullptr, 0.1,0.1,0,0) ;
	rect.setBehaviour(new Stiffness( 10e9, 0.2 ) ) ;
        
	std::vector<VoronoiGrain> grains ;
	LogarithmicCreepWithExternalParameters test("young_modulus = 15e9, poisson_ratio = 0.2") ;
	test.plane = PLANE_STRAIN ;
        test.addMaterialLaw( new MineralMaterialLaw("/home/ag3/Code/amie-ornl/data_minerals/labradorite-elasticity.sci", "-", -1, 1e9 ) ) ;
	test.addMaterialLaw( new GetParticleOrientationMaterialLaw("angle_z") ) ;
	test.addMaterialLaw( new UniformDistributedPerParticleMaterialLaw("angle_x", -0.1*M_PI, 0.1*M_PI ) ) ;
	test.addMaterialLaw( new UniformDistributedPerParticleMaterialLaw("angle_y", -0.1*M_PI, 0.1*M_PI ) ) ;
	grains.push_back(VoronoiGrain(&test, 0.005, 1, 1) ) ;

	FeatureTree f(&rect) ;
	PSDGenerator::get2DVoronoiPolygons( &f, grains, 0, 0.005, 0.01, 32, true, 0) ;
	f.setSamplingNumber( 64 ) ;
	f.setSamplingRestriction(0.002) ;
		
	f.step() ;

	std::string name = "toto" ;
 
	TriangleWriter trg( name.c_str(), &f, 1. ) ;
	trg.getField( "angle_x" ) ;
	trg.getField( "angle_y" ) ;
	trg.getField( "angle_z" ) ;
	trg.getField( "E11" ) ;
	trg.getField( "E22" ) ;
	trg.getField( "E33" ) ;
	trg.write() ;

    return 0 ;
}
