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

Sample box(nullptr, 0.08, 0.08,0.,0.) ;


int main(int argc, char *argv[])
{
	double timestep = 1. ;//atof(argv[1]) ;
	FeatureTree F(&box) ;
	F.setSamplingNumber(50) ;

	ViscoElasticOnlyPasteBehaviour dpaste ;
	SpaceTimeFiberBasedIsotropicLinearDamage * dampaste = new SpaceTimeFiberBasedIsotropicLinearDamage( 0.05, atof(argv[1]) ) ;
	dampaste->setLogitViscousDamageLaw(0.025, 0.3, 2.5) ;
	ViscoelasticityAndFracture paste( PURE_ELASTICITY, dpaste.param, /*dpaste.param*0.3, dpaste.param*10,*/ new SpaceTimeNonLocalMaximumStrain( 0.0001, 0.0001*dpaste.param[0][0]), dampaste ) ;
	paste.getFractureCriterion()->setMaterialCharacteristicRadius(0.003) ;
	box.setBehaviour(&paste);

	F.setSamplingFactor(&box, 2.) ;

	ViscoElasticOnlyAggregateBehaviour dagg ;
	std::vector<Feature *> feats = PSDGenerator::get2DConcrete(&F, new Viscoelasticity(PURE_ELASTICITY, dagg.param, 0), 500, 0.004, 0.0006, new PSDBolomeA(), CIRCLE, 1., M_PI, 50000, 0.5) ;

	
	F.setOrder(LINEAR_TIME_LINEAR) ;
	F.setDeltaTime(timestep) ;
	F.setMaxIterationsPerStep(10) ;
	F.setMinDeltaTime(atof(argv[1])) ;

	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_LEFT_AFTER, 0,0)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0,1)) ;
/*	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0,2)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0,3)) ;*/
	F.step() ;

	BoundingBoxDefinedBoundaryCondition * disp = new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, TOP_AFTER, 0.000008, 1) ;
	F.addBoundaryCondition(disp) ;

	std::vector<DelaunayTriangle *> trg = F.getElements2D() ;
	int count = 0 ;
	for(size_t i = 0 ; i < trg.size() ; i++)
	{
		if(trg[i]->getBehaviour()->getFractureCriterion())
		{
			if(std::abs(trg[i]->getCenter().x) > 0.038 && std::abs(trg[i]->getCenter().y) < 0.002)
			{
				dynamic_cast<Viscoelasticity *>(trg[i]->getBehaviour())->param *= 2. ;
				count++ ;
			}
		}
	}
	std::cout << "paste weakened in " << count << " elements" << std::endl ;

 	
	std::fstream out ;
	std::string name = "tensile_" ;
	name.append(argv[1]) ;
	
	out.open(name.c_str(), std::ios::out) ;

	double time = timestep ;	
	Vector stress = F.getAverageField(REAL_STRESS_FIELD,-1,1) ;
	Vector strain = F.getAverageField(STRAIN_FIELD,-1,1) ;
	Vector rate = F.getAverageField(STRAIN_RATE_FIELD,-1,1) ;
	Vector displ = F.getDisplacements() ;
	out << std::setprecision(16) << 0. << "\t" << displ.max() << "\t" << displ.min() << "\t" <<  stress[0] << "\t" << stress[1] << "\t" << stress[2] << 
		"\t" << strain[0] << "\t" << strain[1] << "\t" << strain[2] << 
		"\t" << rate[0] << "\t" << rate[1] << "\t" << rate[2] << std::endl ;

	std::string nametrg = name ;
	nametrg.append("_trg_0") ;
	TriangleWriter writer(nametrg, &F, 1) ;
	writer.getField(STRAIN_FIELD) ;
	writer.getField(REAL_STRESS_FIELD) ;
	writer.getField(TWFT_STIFFNESS) ;
	writer.getField(TWFT_DAMAGE) ;
	writer.write() ;

	bool goOn = false ;
	int i = 0 ;
	
	while(!goOn)
	{
		i++ ;
		goOn = F.step() ;

		stress = F.getAverageField(REAL_STRESS_FIELD,-1,-1) ;
		strain = F.getAverageField(STRAIN_FIELD,-1,-1) ;
		rate = F.getAverageField(STRAIN_RATE_FIELD,-1,-1) ;
		displ = F.getDisplacements(-1, false) ;
		out << std::setprecision(16) << (F.getElements2D()[0])->getBoundingPoint(0).t << "\t" << displ.max() << "\t" << displ.min() << "\t" << "\t" << stress[0] << "\t" << stress[1] << "\t" << stress[2] << 
			"\t" << strain[0] << "\t" << strain[1] << "\t" << strain[2] << 
			"\t" << rate[0] << "\t" << rate[1] << "\t" << rate[2] << std::endl ;
		
		std::string nametrg = name ;
		nametrg.append("_trg_") ;
		nametrg.append(std::to_string((int) i)) ;
		TriangleWriter writer(nametrg, &F, -1) ;
		writer.getField(STRAIN_FIELD) ;
		writer.getField(REAL_STRESS_FIELD) ;
		writer.getField(TWFT_STIFFNESS) ;
		writer.getField(TWFT_DAMAGE) ;
		writer.write() ;
		std::cout << nametrg << std::endl ;
	}

	stress = F.getAverageField(REAL_STRESS_FIELD,-1,1) ;
	strain = F.getAverageField(STRAIN_FIELD,-1,1) ;
	rate = F.getAverageField(STRAIN_RATE_FIELD,-1,1) ;
	displ = F.getDisplacements(-1, false) ;
	out << std::setprecision(16) << 1. << "\t" << displ.max() << "\t" << displ.min() << "\t" << "\t" << stress[0] << "\t" << stress[1] << "\t" << stress[2] << 
			"\t" << strain[0] << "\t" << strain[1] << "\t" << strain[2] << 
			"\t" << rate[0] << "\t" << rate[1] << "\t" << rate[2] << std::endl ;

	std::string nametrg2 = name ;
	nametrg2.append("_trg_") ;
	nametrg2.append(std::to_string((int) i+1)) ;
	TriangleWriter writer2(nametrg2, &F, 1) ;
	writer2.getField(STRAIN_FIELD) ;
	writer2.getField(REAL_STRESS_FIELD) ;
	writer2.getField(TWFT_STIFFNESS) ;
	writer2.getField(TWFT_DAMAGE) ;
	writer2.write() ;

		
	return 0 ;
}
