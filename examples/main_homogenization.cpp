
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
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
#include "../physics/homogeneised_behaviour.h"
#include "../physics/stiffness.h"
#include "../physics/parallel_behaviour.h"
//#include "../physics/generalized_spacetime_viscoelasticity.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../features/pore.h"
#include "../utilities/writer/triangle_writer.h"
#include "../physics/materials/gel_behaviour.h"
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

#include <fstream>

#include <cmath>
#include <typeinfo>
#include <limits>
#include <sys/time.h>

using namespace Amie ;

int main(int argc, char *argv[])
{
  	Sample box(nullptr, 0.01,0.01,0.,0.) ;
	box.setBehaviour(new ElasticOnlyAggregateBehaviour()) ;
	// no intersection
//	ExpansiveZone * zone = new ExpansiveZone(nullptr, 0.001, 0.001,0.001, new GelBehaviour()) ;
	// intersect without node inside
	ExpansiveZone * zone = new ExpansiveZone(nullptr, 0.0012, -0.0005,0.000, new GelBehaviour(1e6)) ;
	// intersect with node inside
//	ExpansiveZone * zone = new ExpansiveZone(nullptr, 0.001, 0.0015,0.001, new GelBehaviour()) ;
	
	FeatureTree F(&box) ;
	F.setSamplingNumber(4) ;
	F.setOrder(LINEAR) ;
	F.addFeature(&box, zone) ;
	F.addPoint(&(zone->getCenter())) ;

	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT) ) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM) ) ;

	F.step() ;

	Vector strain = F.getAverageField(STRAIN_FIELD) ;
	Vector stress = F.getAverageField(REAL_STRESS_FIELD) ;

	std::cout << strain[0] << "\t" << stress[0] << std::endl ;

	for(auto i = F.get2DMesh()->begin() ; i != F.get2DMesh()->end() ; i++)
	{
		if(dynamic_cast<HomogeneisedBehaviour *>(i->getBehaviour()))
			std::cout << "homogeneised" << std::endl ;
		if(zone->intersects(i->getPrimitive()))
			std::cout << "intersection" << std::endl  ;
		for(size_t j = 0 ; j < i->getBoundingPoints().size() ; j++)	
		{
			if(zone->in(i->getBoundingPoint(j)))
				i->getBoundingPoint(j).print() ;
		}
	}
  
	TriangleWriter writer("test-zone", &F) ;
	writer.getField(STRAIN_FIELD) ;
	writer.getField(REAL_STRESS_FIELD) ;
	writer.getField(TWFT_STIFFNESS) ;
	writer.write() ;

	return 0 ;
}
