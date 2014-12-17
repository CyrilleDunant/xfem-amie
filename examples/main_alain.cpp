// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../utilities/samplingcriterion.h"
#include "../features/features.h"
#include "../physics/physics_base.h"
#include "../physics/kelvinvoight.h"
#include "../physics/maxwell.h"
#include "../physics/stiffness.h"
#include "../physics/dual_behaviour.h"
#include "../physics/logarithmic_creep.h"
#include "../physics/logarithmic_creep_with_external_parameters.h"
#include "../physics/parallel_behaviour.h"
//#include "../physics/generalized_spacetime_viscoelasticity.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/spacetimemultilinearsofteningfracturecriterion.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../features/pore.h"
#include "../utilities/writer/triangle_writer.h"
#include "../physics/materials/gel_behaviour.h"
#include "../physics/material_laws/material_laws.h"
#include "../physics/material_laws/temperature_material_laws.h"
#include "../physics/orthotropicstiffness.h"
#include "../physics/materials/paste_behaviour.h"
#include "../physics/materials/aggregate_behaviour.h"
#include "../physics/homogenization/homogenization_base.h"
#include "../features/sample.h"
#include "../features/sample3d.h"
#include "../features/inclusion.h"
#include "../features/expansiveZone.h"
#include "../features/crack.h"
#include "../features/features.h"
#include "../features/enrichmentInclusion.h"
#include "../mesher/delaunay_3d.h"
#include "../solvers/assembly.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"
#include "../utilities/itoa.h"
#include "../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h"
#include "../physics/fracturecriteria/maxstrain.h"

#include <fstream>
#include <omp.h>
#include <cmath>
#include <typeinfo>
#include <limits>
#include <sys/time.h>
#define DEBUG


using namespace Amie ;

int main(int argc, char *argv[])
{
//    omp_set_num_threads(1) ;

    Sample box(nullptr, 0.08,0.08,0.,0.) ;


    FeatureTree F(&box) ;
    F.setSamplingNumber(32) ;
    F.setOrder(LINEAR_TIME_LINEAR) ;
    double time_step = 0.01 ;
    F.setDeltaTime(time_step) ;
    F.setMinDeltaTime(1e-9) ;
    F.setSamplingRestriction( SAMPLE_RESTRICT_4 ) ;

    std::vector<Point> p ;
    p.push_back( Point( 1e6/1e9, 1e6 ) ) ;
    p.push_back( Point( 2e6/1e9, 0e6 ) ) ;
    std::vector<Point> c ;
    AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion crit( p, c, 1e9/*, 0.1*/ ) ;
    
    LogarithmicCreepWithExternalParameters paste("young_modulus = 12e9, poisson_ratio = 0.2, creep_modulus = 30e9, creep_characteristic_time = 1, creep_poisson = 0.2", &crit, new SpaceTimeFiberBasedIsotropicLinearDamage(0.9999, 1e-9, 0.99)) ;
    box.setBehaviour( &paste );

    LogarithmicCreepWithExternalParameters aggregates("young_modulus = 60e9, poisson_ratio = 0.2") ;

    PSDGenerator::get2DConcrete(&F, &aggregates, 500, 0.01, 0.00001, new PSDBolomeA(), CIRCLE, 1., M_PI, 1000000, 0.9) ;

    F.step() ;

    F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM_AFTER) ) ;
    F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, LEFT_AFTER ) ) ;
    BoundingBoxDefinedBoundaryCondition * load = new BoundingBoxDefinedBoundaryCondition( SET_STRESS_ETA, TOP_AFTER, 1e6 ) ;
    F.addBoundaryCondition(load) ;

    F.setMaxIterationsPerStep(256) ;

    std::string toto = "test_damage";
    std::fstream out ;
    out.open(toto.c_str(), std::ios::out) ;
    int i = 0 ;

    while(F.getCurrentTime() < 100)
    {
	i++ ;
	time_step *= 1.1 ;
	F.setDeltaTime(time_step) ;
//	load->setData(0.0001*F.getCurrentTime()) ;
	while(!F.step()) { 
		i++ ;
		std::cout << std::endl ;
		std::string test = "test_" ;
		test.append(itoa(i)) ;
		TriangleWriter writer(test.c_str(), &F, 1.) ;
		writer.getField(STRAIN_FIELD) ;
		writer.getField(PRINCIPAL_REAL_STRESS_FIELD) ;
		writer.getField(SCALAR_DAMAGE_FIELD) ;
		writer.getField(TWFT_STIFFNESS) ;
		writer.write() ;
	} ;
	std::cout << "\n" << F.getAverageField( STRAIN_FIELD, -1, 1.)[1] << "\t" << F.getAverageField( REAL_STRESS_FIELD, -1, 1.)[1] << std::endl; 
	out << F.getCurrentTime() << "\t"<< F.getAverageField( STRAIN_FIELD, -1, 1.)[1] << "\t" << F.getAverageField( REAL_STRESS_FIELD, -1, 1.)[1] << std::endl; 

	std::string test = "test_" ;
	test.append(itoa(i)) ;
	TriangleWriter writer(test.c_str(), &F, 1.) ;
	writer.getField(STRAIN_FIELD) ;
	writer.getField(PRINCIPAL_REAL_STRESS_FIELD) ;
	writer.getField(SCALAR_DAMAGE_FIELD) ;
	writer.getField(TWFT_STIFFNESS) ;
	writer.write() ;

    }



    return 0 ;

}

