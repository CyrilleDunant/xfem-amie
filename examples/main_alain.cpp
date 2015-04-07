// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../features/features.h"
#include "../physics/physics_base.h"
#include "../physics/stiffness.h"
#include "../physics/dual_behaviour.h"
#include "../physics/finite_difference_viscoelasticity.h"
#include "../physics/logarithmic_creep.h"
#include "../physics/logarithmic_creep_with_external_parameters.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
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
#include "../physics/fracturecriteria/spacetimemultilinearsofteningfracturecriterion.h"

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
	omp_set_num_threads(1) ;

	Sample box(nullptr, 0.01,0.01,0.,0.) ;

	FeatureTree F(&box) ;
	F.setSamplingNumber(1) ;
	double time_step = 1. ;
	F.setDeltaTime(time_step) ;
	F.setMaxIterationsPerStep(1000) ;

	F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, LEFT_AFTER) ) ;
	F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM_AFTER) ) ;
//	F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_STRESS_ETA, TOP_AFTER, 1e6) ) ;
	BoundaryCondition* disp = new BoundingBoxDefinedBoundaryCondition( SET_ALONG_ETA, TOP_AFTER, 0. ) ;
	F.addBoundaryCondition(disp) ;

	std::vector<Point> tension ;
	std::vector<Point> compression ;
	tension.push_back( Point( 0.02, 20e6 ) ) ;
	tension.push_back( Point( 0.12, 0e6 ) ) ;
	FractureCriterion* crit = new AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion( tension, compression, 1e9 ) ;
	DamageModel* dam = new SpaceTimeFiberBasedIsotropicLinearDamage( 0.01, 0.1, 0.999 ) ;

	LogarithmicCreepWithExternalParameters paste("young_modulus = 1e9, poisson_ratio = 0.3", crit, dam) ;
//	box.setBehaviour(new /*FiniteDifference*/Viscoelasticity( GENERALIZED_MAXWELL, paste.C, paste.C*2, paste.C*10, BACKWARD_EULER) ) ;
	box.setBehaviour( &paste ) ;

	F.step() ;
	std::cout << F.getCurrentTime() << "\t" << F.getAverageField( STRAIN_FIELD, -1, 1.)[1] << "\t" << F.getAverageField(REAL_STRESS_FIELD, -1, 1.)[1] << std::endl ;

	while(F.getCurrentTime() < 300)
	{
		time_step += 1. ;
		disp->setData( F.getCurrentTime()*0.00001 ) ;
		F.setDeltaTime(time_step) ;
		F.step() ;
                double strain = F.getAverageField( STRAIN_FIELD, -1, 1.)[1] ;
                double stress = 20e6*(1.-(strain-0.02)/0.1) ;
		std::cout << F.getCurrentTime() << "\t" << F.getAverageField( STRAIN_FIELD, -1, 1.)[1] << "\t" << F.getAverageField(REAL_STRESS_FIELD, -1, 1.)[1] << "\t" << stress << std::endl ;
	}


	return 0 ;
}

