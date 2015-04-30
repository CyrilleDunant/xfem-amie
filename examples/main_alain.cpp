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
	Sample s(nullptr, 0.07, 0.07,0,0) ;
	Inclusion * inc = new Inclusion( 0.003,0.,0. ) ;
	GravelPolygonalInclusionGenerator gravel(1.9, 0.5, 2, 12,0,M_PI, 3) ;
	PolygonalSample * pol = dynamic_cast<PolygonalSample*>(gravel.convert( inc )) ;
	s.setBehaviour(new ElasticOnlyPasteBehaviour() ) ;
	pol->setBehaviour(new ElasticOnlyAggregateBehaviour() ) ;


	FeatureTree f(&s) ;
	f.setSamplingNumber(256) ;

	PSDGenerator::get2DConcrete( &f, new ElasticOnlyAggregateBehaviour(), 50, 0.01, 0.000001, new PSDBolomeA(), &gravel ) ;

	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, LEFT ) ) ;
	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM ) ) ;

	f.step() ;

	TriangleWriter trg( "concrete_256", &f, 1.) ;
	trg.getField( TWFT_STIFFNESS ) ;
	trg.write() ;




	return 0 ;
}

