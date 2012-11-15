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
#include "../physics/materials/gel_behaviour.h"
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

double s2d(double s)
{
	return s/(24*60*60) ;
}

double reaction( std::vector<std::pair<ExpansiveZone *, Inclusion*> > zones)
{
	double r = zones[0].first->getRadius() ;
	return M_PI*r*r*zones.size() ;
}

void updateZones(  std::vector<std::pair<ExpansiveZone *, Inclusion*> > zones, double nexta)
{
	double prevr = zones[0].first->getRadius() ;
	if(nexta > reaction(zones))
	{
		double newr = sqrt(nexta/(M_PI*zones.size())) ;
		std::cout << prevr << "\t" << newr << std::endl ;
		for(size_t i = 0 ; i < zones.size() ; i++)
		{
			zones[i].first->setRadius(newr) ;
		}
	}
}

double getReactiveSurface( std::vector<std::pair<ExpansiveZone *, Inclusion*> > zones )
{
	double area = 0 ; 
	Inclusion * last = nullptr ;
	for(size_t i = 0 ; i < zones.size() ; i++)
	{
		if(zones[i].second != last)
			area += zones[i].second->area() ;
		last = zones[i].second ;
	}
	return area ;
}

int main(int argc, char *argv[])
{
	double timeScale = 100 ;
	double tau = 1.;

	FeatureTree F(&box) ;
	F.setSamplingNumber(50) ;
	F.setOrder(LINEAR) ;
	F.setDeltaTime(tau) ;
	
	box.setBehaviour( new ViscoElasticOnlyPasteBehaviour() ) ;
	std::vector<Inclusion *> inclusions = ParticleSizeDistribution::get2DConcrete( &F, new ElasticOnlyAggregateBehaviour(), 0.008, 10) ;
		
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
	F.step() ;
	Vector x = F.getAverageField(STRAIN_FIELD) ;
	Vector y = F.getAverageField(REAL_STRESS_FIELD) ;
	std::cout << 0. << "\t" << x[0] << "\t" << y[0] << std::endl ;

	size_t i = 0 ;
	double time = 0. ;
	std::vector<std::pair<ExpansiveZone *, Inclusion*> > zones = ParticleSizeDistribution::get2DExpansiveZonesInAggregates( &F, inclusions, new GelBehaviour(), 0.00001, 3000, 1) ;
	F.step() ;
	
	std::vector<DelaunayTriangle *> paste = F.getFeature(0)->getElements2D(&F) ;
	std::vector<DelaunayTriangle *> all = F.getElements2D() ;
	
// 	std::fstream out ;
// 	std::string toto = "rag_visco_" ;
// 	toto.append(argv[1]) ;
// 	out.open(toto.c_str(), std::ios::out ) ;
	
	
	x = F.getAverageField(STRAIN_FIELD) ;
	y = F.getAverageField(EFFECTIVE_STRESS_FIELD, paste) ;
	std::cout << time << "\t" << reaction(zones) << "\t" << x[0] << "\t" << y[0] << std::endl ;

	
	while(time < 100)
	{
		i++ ;
 		F.setDeltaTime( tau*i ) ;
		time += tau *i ;
			
		if(time < timeScale )
			updateZones( zones, getReactiveSurface(zones)*0.03*time/timeScale) ;
		else
			updateZones( zones, getReactiveSurface(zones)*0.03) ;
		F.step() ;
		x = F.getAverageField(STRAIN_FIELD) ;
		y = F.getAverageField(REAL_STRESS_FIELD, paste) ;
		std::cout << time << "\t" << reaction(zones)/getReactiveSurface(zones) << "\t" << x[0] << "\t" << y[0] << std::endl ;
	}
		
	return 0 ;
}

