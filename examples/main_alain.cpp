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
 	omp_set_num_threads(1) ;

        Sample rect(nullptr, 0.04,0.04,0,0) ;
	Inclusion * left = new Inclusion( 0.02,-0.025,0. ) ;
	Inclusion * right = new Inclusion( 0.02,0.025,0. ) ;
        Inclusion * son = new Inclusion( 0.02, 0.00, 0 ) ;
	rect.setBehaviour(new ElasticOnlyPasteBehaviour() ) ;
	left->setBehaviour(new ElasticOnlyAggregateBehaviour() ) ;
	right->setBehaviour(new ElasticOnlyAggregateBehaviour(40e9) ) ;
	son->setBehaviour(new ElasticOnlyAggregateBehaviour(25e9) ) ;
//        son->addToMask( left ) ;
        son->addToMask( right ) ;

	FeatureTree f(&rect) ;
	f.setSamplingNumber(512) ;
	std::vector<Feature *> agg = PSDGenerator::get2DConcrete( &f, new ElasticOnlyAggregateBehaviour(),250, 0.002, 0.0002, nullptr, new GravelPolygonalInclusionGenerator(1.9,0.2,2,10,0,M_PI,3), 10000) ;
        PSDGenerator::get2DMaskedInclusions( &f, new ElasticOnlyAggregateBehaviour(40e9), agg, 250, 0.0008, 0.0001, new ConstantSizeDistribution(), new PolygonalInclusionGenerator(5,0,M_PI,0), 10000,0.5, nullptr, std::vector<Geometry *>(), 20 ) ;

	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_STRESS_ETA, TOP, -1e6 ) ) ;
//	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_PROPORTIONAL_DISPLACEMENT_XI_ETA, TOP, -1 ) ) ; // ux = 0.5 u_y
//        Point n(-0.004,0.008) ;
//	f.addBoundaryCondition( new GeometryAndFaceDefinedSurfaceBoundaryCondition( SET_TANGENT_DISPLACEMENT, dynamic_cast<Polygon*>(&s), n, 0.0001 ) ) ;
	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, BOTTOM_LEFT) ) ;
	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM ) ) ;

	f.step() ;

//        f.getAssembly()->print() ;


	TriangleWriter trg( "larger", &f, 1.) ;
	trg.getField( TWFT_STIFFNESS ) ;
	trg.write() ;




	return 0 ;
}

