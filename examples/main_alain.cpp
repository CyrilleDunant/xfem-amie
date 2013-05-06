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
#include "../physics/fracturecriteria/mohrcoulomb.h"
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
#include <GL/glut.h>
#include <sys/time.h>
#define DEBUG

#define ID_QUIT 1
#define ID_ZOOM 5
#define ID_UNZOOM 6
#define ID_NEXT10 7
#define ID_NEXT100 3
#define ID_NEXT1000 4
#define ID_NEXT 2
#define ID_NEXT_TIME 0
#define ID_REFINE 8
#define ID_AMPLIFY 9
#define ID_DEAMPLIFY 10

#define ID_DISP 11
#define ID_STRAIN_XX 12
#define ID_STRAIN_XY 13
#define ID_STRAIN_YY 14
#define ID_STRESS_XX 15
#define ID_STRESS_XY 16
#define ID_STRESS_YY 17
#define ID_STIFNESS 18
#define ID_ELEM 19
#define ID_VON_MISES 20
#define ID_ANGLE 22
#define ID_ENRICHMENT 21

#define DISPLAY_LIST_DISPLACEMENT 1
#define DISPLAY_LIST_ELEMENTS 2
#define DISPLAY_LIST_STRAIN_XX 3
#define DISPLAY_LIST_STRAIN_YY 4
#define DISPLAY_LIST_STRAIN_XY 5
#define DISPLAY_LIST_STRESS_XX 6
#define DISPLAY_LIST_STRESS_YY 7
#define DISPLAY_LIST_STRESS_XY 8
#define DISPLAY_LIST_CRACK 9
#define DISPLAY_LIST_STIFFNESS 10
#define DISPLAY_LIST_VON_MISES 11
#define DISPLAY_LIST_ANGLE 23
#define DISPLAY_LIST_ENRICHMENT 12
#define DISPLAY_LIST_STIFFNESS_DARK 24

using namespace Mu ;

Sample box(nullptr, 0.08, 0.08,0.,0.) ;


int main(int argc, char *argv[])
{
	double timestep = atof(argv[1]) ;
	int sampling = (int) atof(argv[2]) ;
	int ninc = (int) atof(argv[3]) ;

	double ekv = atof(argv[4]) ;
	double tkv = atof(argv[5]) ;
	double tmx = atof(argv[6]) ;
  
	FeatureTree F(&box) ;
	F.setSamplingNumber(sampling) ;
	
// 	Matrix e = (new ElasticOnlyPasteBehaviour(12e9, 0.3))->param ;
// 	Matrix a = (new ElasticOnlyAggregateBehaviour(70e9, 0.3))->param ;
	box.setBehaviour(new ElasticOnlyPasteBehaviour());
//   	box.setBehaviour(new Viscoelasticity(BURGER, e*ekv, e*ekv*tkv, e, e*tmx)) ;
// 	Viscoelasticity * agg = new Viscoelasticity(PURE_ELASTICITY, a,2) ;
	ElasticOnlyAggregateBehaviour * agg = new ElasticOnlyAggregateBehaviour() ;

	ParticleSizeDistribution::get2DConcrete(&F, agg, ninc, 0.008, 0.001, BOLOME_A, TRIANGLE, 0.25) ;
	
	F.setOrder(LINEAR) ;
	F.setDeltaTime(timestep) ;

	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0,0)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0,1)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0,2)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0,3)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0,4)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0,5)) ;
// 	F.addBoundaryCondition(new TimeContinuityBoundaryCondition()) ;
	F.step() ;

// 	F.getAssembly()->setEpsilon(1e-14) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition( SET_STRESS_ETA, TOP_AFTER, -10e6, 1)) ;
// 	F.step() ;
// 	
// 	std::fstream out ;
// 	std::string name = "visco_" ;
// 	name.append(argv[1]) ;
// 	name.append("_") ;
// 	name.append(argv[2]) ;
// 	name.append("_") ;
// 	name.append(argv[3]) ;
// 	
// 	out.open(name.c_str(), std::ios::out) ;
// 
// 	double time = timestep ;	
// 	Vector stress = F.getAverageField(REAL_STRESS_FIELD,-1,1) ;
// 	Vector strain = F.getAverageField(STRAIN_FIELD,-1,1) ;
// 	Vector rate = F.getAverageField(STRAIN_RATE_FIELD,-1,1) ;
// 	Vector disp = F.getDisplacements() ;
// 	out << std::setprecision(16) << time << "\t" << disp.max() << "\t" << disp.min() << "\t" <<  stress[0] << "\t" << stress[1] << "\t" << stress[2] << 
// 		"\t" << strain[0] << "\t" << strain[1] << "\t" << strain[2] << 
// 		"\t" << rate[0] << "\t" << rate[1] << "\t" << rate[2] << std::endl ;

	if(true)
	{
		std::string nametrg = "test" ;
		nametrg.append("_trg_1") ;
		TriangleWriter writer(nametrg, &F, 1) ;
// 		writer.getField(STRAIN_FIELD) ;
// 		writer.getField(REAL_STRESS_FIELD) ;
		writer.getField(TWFT_STIFFNESS) ;
		writer.write() ;
	}
	
/*	while(time < 401)
	{
		F.step() ;

		if(time > 10)
		{
			timestep++ ;
			F.setDeltaTime(timestep) ;
		}

		time += timestep ;
		stress = F.getAverageField(REAL_STRESS_FIELD,-1,1) ;
		strain = F.getAverageField(STRAIN_FIELD,-1,1) ;
		rate = F.getAverageField(STRAIN_RATE_FIELD,-1,1) ;
		disp = F.getDisplacements() ;
		out << std::setprecision(16) << time << "\t" << disp.max() << "\t" << disp.min() << "\t" << "\t" << stress[0] << "\t" << stress[1] << "\t" << stress[2] << 
			"\t" << strain[0] << "\t" << strain[1] << "\t" << strain[2] << 
			"\t" << rate[0] << "\t" << rate[1] << "\t" << rate[2] << std::endl ;
		
		if(time == 25 || time == 50 || time == 100 || time == 200 || time == 400)
		{
			std::string nametrg = name ;
			nametrg.append("_trg_") ;
			nametrg.append(std::to_string((int) time)) ;
			TriangleWriter writer(nametrg, &F, 1) ;
			writer.getField(STRAIN_FIELD) ;
			writer.getField(REAL_STRESS_FIELD) ;
			writer.getField(TWFT_STIFFNESS) ;
			writer.write() ;
		}
	}*/

		
	return 0 ;
}
