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
#include "../utilities/random.h"
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
	double timestep = 1 ;//atof(argv[1]) ;
	int sampling = 100 ;//(int) atof(argv[2]) ;
	bool fullkv = true ; //(bool) atof(argv[3]) ;
	int axis = 1 ;// (int) atof(argv[1]) ;
	int ninc = 10000 ;//(int) atof(argv[3]) ;
	int micro = 0 ;//(int) atof(argv[2]) ;
	double afraction = atof(argv[1]) ;
	double itz = 0.002 ; //atof(argv[7]) ;
	double aspect = 1. ;
	double orientation = M_PI ;
	GeometryType inclusions = CIRCLE ;
	size_t seed = atof(argv[2]) ;

	if(micro == 1)
	{
		inclusions = ELLIPSE ;
		sampling = 230 ;
		aspect = 0.4 ;
		orientation = 0.001 ;
	}
	if(micro == 2)
	{
		inclusions = TRIANGLE ;
		sampling = 450 ;
	}
	
	FeatureTree F(&box) ;
	F.setSamplingNumber(sampling) ;
	Rectangle placement(0.078,0.078,0.,0.) ;

	box.setBehaviour(new ViscoElasticOnlyPasteBehaviour());
	ElasticOnlyAggregateBehaviour agg ;
	Viscoelasticity toto( PURE_ELASTICITY, (agg.param)*1e-20, 2) ;

	std::vector<Feature *> aggs = ParticleSizeDistribution::get2DConcrete(&F, new ViscoElasticOnlyAggregateBehaviour(), ninc, 0.005, itz, PSD_UNIFORM, inclusions, aspect, orientation, ninc*100,0.2,dynamic_cast<Geometry*>(&placement),seed) ;
//	std::cout << aggs.size() << "\t" << aggs.size()*aggs[0]->area() << std::endl ;	

	std::vector<Circle *> aggregates ;
	for(size_t i = 0 ; i < aggs.size() ; i++)
		aggregates.push_back( dynamic_cast<Circle *>(aggs[i]) ) ;


	double placed = 0. ;
	double radius = 0.001 ;
	RandomNumber rnd ;
	Inclusion * last = nullptr ;
	while(placed < afraction*0.08*0.08)
	{
		double x = rnd.uniform(-0.04,0.04) ;
		double y = rnd.uniform(-0.038,0.038) ; 
		Circle test( radius, x, y) ;
		bool intersects = false ;
		for(size_t i = 0 ; i < aggregates.size() ; i++)
		{
			intersects |= aggregates[i]->intersects(&test) ;
			intersects |= aggregates[i]->in(test.getCenter()) ;
		}
		if(!intersects)
		{
			Inclusion * inc = new Inclusion(radius,x,y) ;
			inc->setBehaviour( &toto ) ;
			if(!last)	
				F.addFeature(&box, inc) ;
			else
				F.addFeature(last, inc) ;
			F.setSamplingFactor(inc, 4.) ;
			last = inc ;
			placed += 0.001*0.001*M_PI ;
		}
	}
	


	F.setOrder(LINEAR_TIME_LINEAR) ;
	F.setDeltaTime(timestep) ;

	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0,0)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0,1)) ;
 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0,2)) ;
 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0,3)) ;
 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0,4)) ;
 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0,5)) ;
	F.step() ;

	F.getAssembly()->setEpsilon(1e-12) ;
//	if(axis == 1)
		F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition( SET_STRESS_ETA, TOP_AFTER, -10e6, 1)) ;
//	else if(axis == 0)
//		F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition( SET_STRESS_XI, RIGHT_AFTER, -10e6, 0)) ;
	F.step() ;
 	
	double realarea = 0. ;
	std::vector<DelaunayTriangle *> trg = F.getElements2D() ;
	for(size_t i = 0 ; i < trg.size() ; i++)
	{
		if(trg[i]->getBehaviour()->param[0][0] < 1)
			realarea += trg[i]->area() ;
	}
	realarea /= (0.08*0.08) ;	

	std::fstream out ;
	std::string name = "visco_porosity_inclusions_" ;
/*	if(fullkv)
		name.append("kv_") ;
	if(axis == 0)
		name.append("xi_") ;
	if(axis == 1)
		name.append("eta_") ;
	if(micro == 0)
		name.append("circle_") ;
	name.append(itoa(aggs.size())) ;
	name.append("_") ;
	if(micro == 1)
	{
		name.append("ellipse_") ;
	}
	if(micro == 2)
	{
		name.append("triangle_") ;
	}*/
	name.append(std::to_string(realarea)) ;
	name.append("_") ;
	name.append(std::string(argv[2])) ;
	
	out.open(name.c_str(), std::ios::out) ;

	double time = timestep ;	
	Vector stress = F.getAverageField(REAL_STRESS_FIELD,-1,1) ;
	Vector strain = F.getAverageField(STRAIN_FIELD,-1,1) ;
	Vector rate = F.getAverageField(STRAIN_RATE_FIELD,-1,1) ;
	Vector disp = F.getDisplacements() ;
	out << std::setprecision(16) << time << "\t" << disp.max() << "\t" << disp.min() << "\t" <<  stress[0] << "\t" << stress[1] << "\t" << stress[2] << 
		"\t" << strain[0] << "\t" << strain[1] << "\t" << strain[2] << 
		"\t" << rate[0] << "\t" << rate[1] << "\t" << rate[2] << std::endl ;

	if(true)
	{
		std::string nametrg = name ;
		nametrg.append("_trg_0") ;
		TriangleWriter writer(nametrg, &F, 1) ;
/*		writer.getField(STRAIN_FIELD) ;
		writer.getField(REAL_STRESS_FIELD) ;*/
		writer.getField(TWFT_STIFFNESS) ;
		writer.write() ;
	}
	
//	exit(0) ;
	while(time < 400)
	{
		F.step() ;

		if(time > 10)
		{
			timestep += (int) (time/10) ;
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
		
		if(time == 25 || time == 55 || time == 101 || time == 200 || time == 388)
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
	}

		
	return 0 ;
}
