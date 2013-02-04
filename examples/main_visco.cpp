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
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../features/pore.h"
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

Sample box(nullptr, 1., 1.,0.0,0.0) ;

double s2d(double s)
{
	return s/(24*60*60) ;
}

double getMinDisplacement(const Vector & x)
{
	double ret = x[0] ;
	for(size_t i = 0 ; i < x.size()/6 ; i++)
	{
		if(x[i*6] < ret)
			ret = x[i*6] ;
		if(x[i*6+1] < ret)
			ret = x[i*6+1] ;
	}
	return ret ;
}

int main(int argc, char *argv[])
{
	double tau = 10 ;
//	double tauRadius = atof(argv[1]) ;

/*	if(tauRadius == 0)
		tauRadius = 1e-9 ;  */

	FeatureTree F(&box) ;
	F.setSamplingNumber(10) ;

	Matrix e = (new ElasticOnlyPasteBehaviour())->param ;
	box.setBehaviour(new Viscoelasticity(BURGER, e, e*500, e, e*300)) ;
//	std::vector<Inclusion *> inclusions = ParticleSizeDistribution::get2DConcrete( &F, new GeneralizedSpaceTimeViscoelasticity(PURE_ELASTICITY, (new ElasticOnlyAggregateBehaviour())->param, 2), 0.008, 1000) ;

	F.setOrder(LINEAR_TIME_LINEAR) ;
	F.setDeltaTime(tau) ;

	
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT, 0,0)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM, 0,1)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT, 0,2)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM, 0,3)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT, 0,4)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM, 0,5)) ;
	F.addBoundaryCondition(new TimeContinuityBoundaryCondition()) ;
	srand(0) ;
	F.step() ;

	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition( SET_STRESS_ETA, TOP_AFTER, -10e6) ) ;

	std::string toto = "cyrille" ;
	std::fstream out(toto.c_str(), std::ios::out) ;

	F.step() ;
	double time = tau ;

	out << time << "\t" << F.getDisplacements().min() << std::endl ;
	std::cout << time << "\t" << F.getDisplacements().min() << std::endl ;
	
	while(time < 50)
	{
		time += tau ;
		F.step() ;

		out << time <<  "\t" << F.getDisplacements().max() << std::endl ;
		std::cout << time << "\t" << F.getDisplacements().min() << std::endl ;

	}
	
		
	return 0 ;
}

