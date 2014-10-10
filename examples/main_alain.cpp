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
    Sample box(nullptr, 0.2,0.2,0.,0.) ;

    FeatureTree F(&box) ;
    F.setSamplingNumber(128) ;
    F.setOrder(LINEAR) ;
    double time_step = 0.01 ;
    F.setDeltaTime(time_step) ;
    F.setMinDeltaTime(1e-9) ;
    F.setSamplingRestriction( SAMPLE_RESTRICT_4 ) ;

    ElasticOnlyPasteBehaviour paste ;
    box.setBehaviour( &paste );

    std::vector<Feature *> inc = PSDGenerator::get2DConcrete(&F, new ElasticOnlyAggregateBehaviour(), 1, 0.04, 0.0002, new ConstantSizeDistribution(), ELLIPSE, 0.7, M_PI, 10, 0.8) ;

	std::vector<Geometry *> buffer ;
	for(size_t i = 0 ; i < inc.size() ; i++)
	{
		buffer.push_back(dynamic_cast<Geometry *>(inc[i])) ;
	}

    std::vector<Feature *> inc2 = PSDGenerator::get2DConcrete(&F, new ElasticOnlyAggregateBehaviour(30e9), 1, 0.05, 0.0002, new ConstantSizeDistribution(), CIRCLE, 1., M_PI, 1000, 0.8, nullptr, buffer) ;

	std::vector<Feature *> inc3 = PSDGenerator::get2DEmbeddedInclusions(&F, new VoidForm(), inc2, 3, 0.01, 0.0001, new ConstantSizeDistribution(), CIRCLE, 1, M_PI, 100, 0.8) ;

	std::vector<Geometry *> buffer2 ;
	for(size_t i = 0 ; i < inc3.size() ; i++)
	{
		buffer2.push_back(dynamic_cast<Geometry *>(inc3[i])) ;
	}

	PSDGenerator::get2DEmbeddedInclusions(&F, new VoidForm(), inc2, 20, 0.005, 0.0001, new ConstantSizeDistribution(), ELLIPSE, 0.8, M_PI, 1000, 0.8, nullptr, buffer2) ;

    F.step() ;

	TriangleWriter trg("toto", &F, 1.) ;
	trg.getField(TWFT_STIFFNESS) ;
	trg.write() ;

    return 0 ;

}

