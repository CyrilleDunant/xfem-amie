// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../utilities/samplingcriterion.h"
#include "../features/features.h"
#include "../physics/physics_base.h"
#include "../physics/stiffness.h"
#include "../physics/dual_behaviour.h"
#include "../physics/logarithmic_creep.h"
#include "../physics/logarithmic_creep_with_external_parameters.h"
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
	Sample box(nullptr, 0.04,0.04,0.,0.) ;
	FeatureTree F(&box) ;
	box.setBehaviour( new ElasticOnlyPasteBehaviour( 10e9 ) ) ;

	std::map<Form *, double> behaviour ;
	behaviour[ new ElasticOnlyAggregateBehaviour( 45e9 ) ] = 0.1 ;
	behaviour[ new ElasticOnlyAggregateBehaviour( 50e9 ) ] = 0.25 ;
	behaviour[ new ElasticOnlyAggregateBehaviour( 55e9 ) ] = 0.35 ;
	behaviour[ new ElasticOnlyAggregateBehaviour( 70e9 ) ] = 0.45 ;
	behaviour[ new ElasticOnlyAggregateBehaviour( 75e9 ) ] = 0.6 ;
	behaviour[ new ElasticOnlyAggregateBehaviour( 80e9 ) ] = 0.7 ;
	behaviour[ new ElasticOnlyAggregateBehaviour( 90e9 ) ] = 0.75 ;
	behaviour[ new ElasticOnlyAggregateBehaviour( 95e9 ) ] = 0.95 ;
	behaviour[ new ElasticOnlyAggregateBehaviour( 100e9 ) ] = 1. ;

	std::vector<Feature *> incs = PSDGenerator::get2DConcrete( &F, new ElasticOnlyAggregateBehaviour(), 500 ) ;
	std::vector<PolygonalSample *> poly = PSDGenerator::get2DVoronoiPolygons(&F, behaviour, incs, 200, 0.00002, true ) ;
/*	Inclusion * inc = new Inclusion( 0.01,0,0) ;
	inc->setBehaviour( new ElasticOnlyAggregateBehaviour( 42e9 ) ) ;
	F.addFeature(&box, inc) ;
	for(size_t i = 0 ; i < poly.size() ; i++)
	{
		if(inc->in(poly[i]->getCenter()) || inc->intersects(dynamic_cast<Polygon *>(poly[i])))
		{
			ElasticOnlyAggregateBehaviour toto( 10e9*(i+4) ) ;
			poly[i]->setBehaviour( new  ElasticOnlyAggregateBehaviour(10e9*(i*4+4)) ) ;
			F.addFeature(inc, poly[i]) ;
			poly[i]->addToMask( inc ) ;
		}
	}*/


	F.setSamplingNumber(256) ;
	F.step() ;

        TriangleWriter writer("tata", &F, 1.) ;
	writer.getField(TWFT_STIFFNESS) ;
	writer.write() ;

	return 0 ;
}

