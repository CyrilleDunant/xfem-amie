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
//	omp_set_num_threads(1) ;


	Sample box(nullptr, 0.01,0.01,0.,0.) ;

	FeatureTree F(&box) ;
	F.setSamplingNumber(4) ;
//	F.setOrder(LINEAR_TIME_LINEAR) ;
	double time_step = 1. ;
	F.setDeltaTime(time_step) ;

	F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, LEFT_AFTER) ) ;
	F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM_AFTER) ) ;
	BoundingBoxDefinedBoundaryCondition * disp = new BoundingBoxDefinedBoundaryCondition( SET_ALONG_ETA, TOP_AFTER, 0.00001) ;
	F.addBoundaryCondition(disp) ;

        ViscoElasticOnlyPasteBehaviour falsePaste ;
        falsePaste.freeblocks = 1 ;
	box.setBehaviour(&falsePaste) ;

	LogarithmicCreepWithExternalParameters paste("young_modulus = 15e9, poisson_ratio = 0.3, creep_modulus = 30e9, creep_poisson = 0.3, creep_characteristic_time = 1.") ;
	Inclusion * notch = new Inclusion( 0.001, -0.005, 0.) ;
	notch->setBehaviour(&paste) ;
	F.addFeature(&box, notch) ;

	while(F.getCurrentTime() < 120)
	{
		time_step += 1. ;
		F.setDeltaTime(time_step) ;
		F.step() ;
		std::cout << F.getAssembly()->getNumberOfDegreesOfFreedom() << "\t" << F.getAverageField( STRAIN_FIELD, -1, 1.)[1] << "\t" << F.getAverageField(REAL_STRESS_FIELD, -1, 1)[1] << std::endl ;
	}


	return 0 ;
}

