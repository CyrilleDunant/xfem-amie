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
#include "../physics/parallel_behaviour.h"
#include "../physics/viscoelasticity.h"
#include "../physics/viscoelasticity_and_fracture.h"
#include "../physics/fracturecriteria/maxstrain.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../features/pore.h"
#include "../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h"
#include "../utilities/writer/triangle_writer.h"
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

#include <fstream>

#include <cmath>
#include <typeinfo>
#include <limits>
#include <sys/time.h>
#define DEBUG
using namespace Mu ;

Sample box(nullptr, 0.02, 0.01,0.,0.) ;


int main(int argc, char *argv[])
{
	double timestep = 1. ;//atof(argv[1]) ;
	FeatureTree F(&box) ;
	F.setSamplingNumber(4) ;

	ViscoElasticOnlyPasteBehaviour dpaste ;
	SpaceTimeFiberBasedIsotropicLinearDamage * dampaste = new SpaceTimeFiberBasedIsotropicLinearDamage( 0.05, 1e-9) ;//atof(argv[1]) ) ;
	ViscoelasticityAndFracture paste( PURE_ELASTICITY, dpaste.param, new SpaceTimeNonLocalLinearSofteningMaximumStrain( 0.0001, 0.0001*dpaste.param[0][0], 0.0001*1.5), dampaste ) ;
	paste.getFractureCriterion()->setMaterialCharacteristicRadius(0.003) ;
	box.setBehaviour(&paste);

	F.setSamplingFactor(&box, 2.) ;

	F.setOrder(LINEAR_TIME_LINEAR) ;
	F.setDeltaTime(timestep) ;
	F.setMaxIterationsPerStep(500) ;
	F.setMinDeltaTime(1e-20) ;

	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_LEFT_AFTER, 0,1)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0,0)) ;
	F.step() ;
	std::vector<DelaunayTriangle *> trg = F.getElements2D() ;
	int count = 0 ;
	for(size_t i = 0 ; i < trg.size() ; i++)
	{
		if(trg[i]->getBehaviour()->getFractureCriterion())
		{
			if(std::abs(trg[i]->getCenter().getX()) < 0.0002)
			{
				SpaceTimeNonLocalLinearSofteningMaximumStrain* dam = dynamic_cast<SpaceTimeNonLocalLinearSofteningMaximumStrain *>(trg[i]->getBehaviour()->getFractureCriterion()) ;
				dam->upVal *= 0.99 ;
				dam->maxstress *= 0.99 ;
				dam->yieldstrain *= 0.99 ;
				count++ ;
			}
		}
	}
	std::cout << "paste weakened in " << count << " elements" << std::endl ;

	std::fstream out ;
	std::string name = "tensile_" ;
//	name.append(argv[1]) ;
	out.open(name.c_str(), std::ios::out) ;

	BoundingBoxDefinedBoundaryCondition * disp = new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, RIGHT_AFTER, 0.0, 0) ;
	F.addBoundaryCondition(disp) ;
	for(size_t i = 0 ; i < 100 ; i++)
	{
		disp->setData(0.0000001*i) ;
		if(!F.step())
			break ;

	}

	std::string nametrg = name ;
	nametrg.append("_trg_0") ;
	TriangleWriter writer(nametrg, &F, 1) ;
	writer.getField(STRAIN_FIELD) ;
	writer.getField(REAL_STRESS_FIELD) ;
	writer.getField(TWFT_STIFFNESS) ;
	writer.getField(TWFT_DAMAGE) ;
	writer.write() ;
		
	return 0 ;
}
