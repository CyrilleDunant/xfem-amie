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
#include "../physics/stiffness_and_fracture.h"
#include "../physics/parallel_behaviour.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/maxstrain.h"
#include "../physics/fracturecriteria/creeprupture.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h"
#include "../physics/damagemodels/fiberbasedisotropiclineardamage.h"
#include "../features/pore.h"
#include "../utilities/writer/triangle_writer.h"
#include "../physics/materials/gel_behaviour.h"
#include "../physics/viscoelasticity.h"
#include "../physics/viscoelasticity_and_fracture.h"
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

#include <omp.h>
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

// see Cecot 2001, p98, for sample dimensions

double length = 0.2 ; //5.5 ;
double width = 0.06 ;
double depth = 0.03 ;//2.5 ;
double nnotch = 0.085 ;
double nwidth = 0.0025 ;
Sample box(nullptr, length*0.5, length,length*0.25,length*0.5) ;
Sample top(nullptr, width*0.5, depth, width*0.25, length - depth*0.5) ;
Sample notch(nullptr,nwidth, nnotch, 0.0, length*0.5 - depth - nnotch*0.5) ;

int main(int argc, char *argv[])
{
//	 omp_set_schedule(omp_sched_static, 60) ;
//	omp_set_num_threads(1) ;

	FeatureTree F(&box) ;
	F.setSamplingNumber(64) ;
	F.setOrder(LINEAR_TIME_LINEAR) ;
	F.setMaxIterationsPerStep( 100 ) ;
	double totaltime = atof(argv[1]) ;
	F.setDeltaTime(totaltime/100.) ;
	F.setMinDeltaTime((totaltime/100.)*1e-9) ;
	size_t seed = 0 ;//atof(argv[4]) ;
//	F.setOrder(LINEAR) ;
	
//	Matrix cmx = (new PasteBehaviour(14e9))->param ;
//	Matrix ckv = (new PasteBehaviour(30e9))->param ;
//	Viscoelasticity dummy( BURGER, ckv, ckv*300, cmx, cmx*5000,2) ;

//	Viscoelasticity concrete( dummy.param, dummy.eta, 5) ;

	PseudoBurgerViscoElasticOnlyPasteBehaviour pastenodamage(12e9,0.3,2,10000) ;

	Rectangle * placement= new Rectangle(length*0.7,nnotch, length*(-0.1), nnotch*0.5) ;

//	Sample left(nullptr, (length-width)*.5, length, (length-width)*.25-length*.5, 0) ;
//	left.isVirtualFeature  = true ;
//	left.setBehaviour(&pastenodamage) ;
	
	Sample middle(nullptr, length*0.5, depth+nnotch*0.9, length*0.25, length-(depth+nnotch*0.9)*0.5) ;
	middle.isVirtualFeature = true ;
	middle.setBehaviour(&pastenodamage) ;

	Sample right(nullptr, (length-width)*.5, nnotch*2., length*.5-(length-width)*.25, nnotch) ;
	right.isVirtualFeature  = true ;
	right.setBehaviour(&pastenodamage) ;

	Rectangle refinement( length*0.6,nnotch*1.2, length*0.1, length-(depth+nnotch)*0.5) ;
/*	Rectangle large( width, length, 0.,0.) ;
	Rectangle large2( 0.04, length, 0.,0.) ;
	Rectangle topfine( length, depth*2., 0., length*0.5-depth) ;*/
//	F.addRefinementZone(placement);
//	F.addRefinementZone(&refinement);
//	F.addRefinementZone(&topfine);
// 	F.addRefinementZone(&large);
// 	F.addRefinementZone(&refinement);
	Rectangle support( 0.001, depth, width*0.5, length-depth*0.5) ;
	
	PseudoBurgerViscoDamagePasteBehaviour paste(12e9, 0.3,2, 10000, atof(argv[2]), 0.002) ;
	paste.freeblocks = 0 ;
/*	if(argv[2] == std::string("stress"))
		paste.ctype = STRESS_CRITERION ;
	if(argv[2] == std::string("mixed"))
		paste.ctype = MIXED_CRITERION ;*/

	box.setBehaviour( &paste ) ;

//	box.setBehaviour( new Viscoelasticity(GENERALIZED_KELVIN_VOIGT, c, c*0.3, c*0.3*10)) ;
//	box.setBehaviour( new Viscoelasticity(PURE_ELASTICITY, c) ) ;
//	box.setBehaviour( new StiffnessAndFracture( c, new NonLocalMohrCoulomb( 0.001, -0.008, 15e9) ) ) ;
//	box.setBehaviour( new Stiffness( c ) ) ;
	
//	box.getBehaviour()->getFractureCriterion()->setMaterialCharacteristicRadius(0.00025);	
	top.setBehaviour( new VoidForm() ) ;
	notch.setBehaviour( new VoidForm() ) ;

	/*std::vector<Inclusion *> inclusions = */
//	size_t s = seed ;
	ViscoElasticOnlyAggregateBehaviour agg ;
	agg.freeblocks = 0 ;

	std::vector<Feature *> aggregates = ParticleSizeDistribution::get2DConcrete( &F, &agg, 170, 0.008, 0.000001, BOLOME_A, CIRCLE, 1., M_PI, 100000, 0.9 ) ;
	for(size_t i = 0 ; i < aggregates.size() ; i++)
	{
		if(aggregates[i]->in( Point(0., nnotch) ) )
		{
			F.removeFeature(aggregates[i]) ;
			continue ;
		}

//		if(placement->in(aggregates[i]->getCenter()))
			F.setSamplingFactor(aggregates[i], 4.) ;
/*		else
		{
			if(aggregates[i]->getRadius() < 0.006)
				aggregates[i]->setBehaviour( &pastenodamage ) ;
			else 
				F.setSamplingFactor(aggregates[i], 8.) ;
		}*/
	}
//	for(size_t i = 1000 ; i < 2000 ; i++)
//		F.setSamplingFactor(aggregates[i], 3.) ;

/*	Inclusion * support = new Inclusion(nullptr, 0.008, 0., -length*0.5) ;
	support->setBehaviour( &agg) ;*/

 	F.addFeature(&box, &top) ;
	F.addFeature(&box, &middle) ;
// 	F.addFeature(&middle, &notch) ;
//	F.addFeature(&box, support) ;
//	F.setSamplingFactor(support, 2.) ;
//	F.addFeature(&box, &left) ;
	F.addFeature(&box, &right) ;
	F.setSamplingRestriction( SAMPLE_RESTRICT_4 ) ;

	F.setSamplingFactor( &top, 5. ) ;
//	F.setSamplingFactor( &notch, 2. ) ;
//	F.setSamplingFactor( support, 2. ) ;
	
 	F.addBoundaryCondition(new BoundingBoxAndRestrictionDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, LEFT_AFTER, -1e-5, 1e-5, -1e-5, nnotch+1e-5,0, 0 )) ;
 	F.addBoundaryCondition(new BoundingBoxAndRestrictionDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, LEFT_AFTER, -1e-5, 1e-5, -1e-5, nnotch+1e-5, 0, 2 )) ;
 	F.addBoundaryCondition(new BoundingBoxAndRestrictionDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, LEFT_AFTER, -1e-5, 1e-5, -1e-5, nnotch+1e-5, 0, 4 )) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0, 6 )) ;
 //	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0, 8 )) ;

	F.addBoundaryCondition( new BoundingBoxNearestNodeDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, Point(0., 0.), 0, 1 ) ) ;
 	F.addBoundaryCondition( new BoundingBoxNearestNodeDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, Point(0., 0.), 0, 3 ) ) ;
 	F.addBoundaryCondition( new BoundingBoxNearestNodeDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, Point(0., 0.), 0, 5 ) ) ;
// 	F.addBoundaryCondition( new BoundingBoxNearestNodeDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, Point(0., -length*0.5), 0, 7 ) ) ;
 //	F.addBoundaryCondition( new BoundingBoxNearestNodeDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, Point(0., -length*0.5), 0, 9 ) ) ;

//	BoundingBoxAndRestrictionDefinedBoundaryCondition * disp = new BoundingBoxAndRestrictionDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, TOP_AFTER, -0.1*width, width*0.6, length*0.99, length*1.01, 0., 0 ) ;
	BoundingBoxNearestNodeDefinedBoundaryCondition * dispx = new BoundingBoxNearestNodeDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, TOP_AFTER, Point(width*0.5, length), 0., 0) ;
	BoundingBoxNearestNodeDefinedBoundaryCondition * dispy = new BoundingBoxNearestNodeDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, TOP_AFTER, Point(width*0.5, length), 0., 1) ;
	F.addBoundaryCondition(dispx) ;
	F.addBoundaryCondition(dispy) ;

	
	F.step() ;
	F.getAssembly()->setEpsilon(1e-8) ;
 	Vector x = F.getAverageField(STRAIN_FIELD) ;
 	Vector y = F.getAverageField(REAL_STRESS_FIELD) ;
//	std::cout << 0. << "\t" << x[0] << "\t" << y[0] << std::endl ;

	size_t i = 0 ;
	size_t j = 0 ;

	std::fstream out ;
	std::string tata = "wedge_" ;
	tata.append(argv[1]) ;
	tata.append("_strain_") ;
	tata.append(argv[2]) ;
	tata.append(".txt") ;
	out.open(tata.c_str(), std::ios::out) ;

	double speed = 0.0005/totaltime ;
	double totaldisp = 0.;
	bool goOn = true ;
	std::vector<DelaunayTriangle *> trg = F.getElements2D() ;

	while(speed*F.getCurrentTime() <= 0.0005)
	{
		if(goOn)
		{
			i++ ;
			totaldisp = speed * F.getCurrentTime() ;
			dispx->setData( totaldisp*0.5 ) ;
			dispy->setData( totaldisp*0.25/(-1.8666) ) ;
			F.setDeltaTime(totaltime/100.) ;
			F.setMinDeltaTime((totaltime/100.)*1e-9) ;
			std::cout << totaldisp*0.5 << std::endl ;
		}

		goOn = F.step() ;

		if(goOn)
		{
			std::string tati = tata ; tati.append("_") ;
			tati.append(itoa(i)) ;
			std::cout << tati << std::endl ;
			TriangleWriter writer(tati, &F, 1) ;
	 		writer.getField(STRAIN_FIELD) ;
	 		writer.getField(REAL_STRESS_FIELD) ;
			writer.getField(TWFT_DAMAGE) ;
			writer.getField(TWFT_STIFFNESS) ;
			writer.write() ;
			j = 0 ;
		}
		if(!goOn)
		{
			j++ ;
			std::string tati = tata ; tati.append("_") ;
			tati.append(itoa(i)) ; tati.append("_inter_") ;
			tati.append(itoa(j)) ;
			std::cout << tati << std::endl ;
			TriangleWriter writer(tati, &F, -1) ;
	 		writer.getField(STRAIN_FIELD) ;
	 		writer.getField(REAL_STRESS_FIELD) ;
			writer.getField(TWFT_DAMAGE) ;
			writer.getField(TWFT_STIFFNESS) ;
			writer.write() ;

		}
			
 		x = F.getAverageField(STRAIN_FIELD, -1, -1+2*goOn) ;
 		y = F.getAverageField(REAL_STRESS_FIELD, -1, -1+2*goOn) ;
		out << trg[0]->getBoundingPoint(0+3*goOn).t << "\t" << totaldisp*0.5 << "\t" << x[0] << "\t" << y[0]*length*length << "\t" << F.averageDamage << std::endl ;

	}
	
//	F.getAssembly()->print() ;
	
	return 0 ;
}

