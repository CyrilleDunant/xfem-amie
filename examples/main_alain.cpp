// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../features/features.h"
#include "../features/microstructuregenerator.h"
#include "../features/polygonSample.h"
#include "../physics/physics_base.h"
#include "../physics/stiffness.h"
#include "../physics/dual_behaviour.h"
#include "../physics/finite_difference_viscoelasticity.h"
#include "../physics/logarithmic_creep.h"
#include "../physics/homogenization/phase.h"
#include "../physics/homogenization/composite.h"
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
#include "../features/polygonSample.h"
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
// 	omp_set_num_threads(1) ;

        Sample rect(nullptr, 0.04,0.04,0,0) ;
//	rect.setBehaviour(new Viscoelasticity(MAXWELL, ElasticOnlyPasteBehaviour(10e9).param, ElasticOnlyPasteBehaviour(10e9).param*10 ) ) ;
	rect.setBehaviour(new LogarithmicCreepWithExternalParameters("young_modulus = 10e9, poisson_ratio = 0.2, creep_characteristic_time = 5, creep_modulus = 20e9, creep_poisson = 0.2, imposed_deformation = 0.1") ) ;

	FeatureTree f(&rect) ;
        f.setOrder( LINEAR_TIME_LINEAR ) ;
	f.setSamplingNumber(32) ;
        f.setDeltaTime(1) ;
	
	f.step() ;

//	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_STRESS_ETA, TOP_AFTER, 1e6 ) ) ;
//	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_PROPORTIONAL_DISPLACEMENT_XI_ETA, TOP, -1 ) ) ; // ux = 0.5 u_y
//        Point n(-0.004,0.008) ;
//	f.addBoundaryCondition( new GeometryAndFaceDefinedSurfaceBoundaryCondition( SET_TANGENT_DISPLACEMENT, dynamic_cast<Polygon*>(&s), n, 0.0001 ) ) ;
	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, BOTTOM_LEFT_AFTER) ) ;
	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM_AFTER ) ) ;
	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_STRESS_ETA, TOP_AFTER, -1e6 ) ) ;

	for(size_t i = 0 ; i < 1 ; i++)
	{
                f.setDeltaTime((double) i+2) ;
		f.step() ;
		std::cout << f.getCurrentTime() << "\t" << f.getAverageField( REAL_STRESS_FIELD, -1,1 )[1] << "\t" << f.getAverageField( STRAIN_FIELD, -1,1 )[1] << std::endl ;
	}

	return 0 ;
}

