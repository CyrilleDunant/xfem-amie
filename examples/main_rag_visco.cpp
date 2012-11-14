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
#include "../physics/generalized_spacetime_viscoelasticity.h"
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

Sample box(nullptr, 0.07, 0.07,0.0,0.0) ;
Rectangle rbox(0.07,0.07,0.0,0.0) ;
double tau = 1.;

double s2d(double s)
{
	return s/(24*60*60) ;
}

int main(int argc, char *argv[])
{
	ElasticOnlyPasteBehaviour * paste_e = new ElasticOnlyPasteBehaviour() ;
  
	FeatureTree F(&box) ;
	F.setSamplingNumber(100) ;
	F.setOrder(LINEAR) ;
	F.setDeltaTime(tau) ;
	
	box.setBehaviour( new ViscoElasticOnlyPasteBehaviour() ) ;
	
	std::vector<Inclusion *> inclusions = ParticleSizeDistribution::get2DConcrete(0.008, 0.07, 50) ;//, masseInitiale, BOLOME_A, PSDEndCriteria(-1, 0.001, inclusionNumber)) ;
	std::vector<Feature *> feats ;
	for( size_t i = 0; i < inclusions.size() ; i++ )
		feats.push_back( inclusions[i] ) ;
	inclusions.clear() ;

	int nAgg = 1 ;
	srand(0) ;
	feats = placement( &rbox, feats, &nAgg, 0, 6400 );
	for( size_t i = 0; i < feats.size() ; i++ )
		inclusions.push_back( static_cast<Inclusion *>( feats[i] ) ) ;
	
	ElasticOnlyAggregateBehaviour *stiff = new ElasticOnlyAggregateBehaviour() ;
	GeneralizedSpaceTimeViscoelasticity * agg = new GeneralizedSpaceTimeViscoelasticity( PURE_ELASTICITY, stiff->param, 7) ;
	for( size_t i = 0 ; i < inclusions.size() ; i++ )
	{
// 		inclusions[i]->setBehaviour( agg ) ;
		inclusions[i]->setBehaviour( stiff ) ;
		F.addFeature( &box, inclusions[i] ) ;
	}
	
		
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
	F.addBoundaryCondition(new TimeContinuityBoundaryCondition()) ;
	F.step() ;
 	F.step() ;
 	size_t i = 0 ;
   	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP/*_AFTER*/, -10e6)) ;

	
	double time = 0. ;
	F.step() ;
	Vector x = F.getDisplacements() ;
	double hey = x.max() ;
	std::fstream out ;
	std::string toto = "/home/alain/Code/FitExperiment/fit10_" ;
	toto.append(argv[1]) ;
	out.open(toto.c_str(), std::ios::out) ;
	out << time << "\t" << x.min()/(0.07*10) << std::endl ;
	
	while(time < 365)
	{
		i++ ;
 		F.setDeltaTime( tau*i ) ;
		F.step() ;
		x = F.getDisplacements() ;
		time += tau *i ;
		out << time << "\t" << x.min()/(0.07*10) << std::endl ;
	}
		
	return 0 ;
}

