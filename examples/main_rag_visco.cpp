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
	if(zones.empty())
		return 0 ;
	double r = zones[0].first->getRadius() ;
	return M_PI*r*r*zones.size() ;
}

void updateZones(  std::vector<std::pair<ExpansiveZone *, Inclusion*> > zones, double nexta)
{
	if(zones.empty())
		return ;
		
	double prevr = zones[0].first->getRadius() ;
	if(nexta > reaction(zones))
	{
		double newr = sqrt(nexta/(M_PI*zones.size())) ;
//		std::cout << prevr << "\t" << newr << std::endl ;
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

double getDamagedAggregateArea( std::vector<std::pair<ExpansiveZone *, Inclusion *> > zones, FeatureTree * F)
{
	double area = 0 ;
	Inclusion * last = nullptr ;
	std::vector<DelaunayTriangle *> agg ;
	for(size_t i = 0 ; i < zones.size() ; i++)
	{
		agg.clear() ;
		if(zones[i].second != last)
		{
			agg = zones[i].second->getElements2D(F) ;
			for(size_t j = 0 ; j < agg.size() ; j++)
			{
				if(agg[j]->getBehaviour()->getDamageModel())
				{
					area += (agg[j]->area())*(agg[j]->getBehaviour()->getDamageModel()->getState()[0]) ;
				}
			}
		}
		last = zones[i].second ;
	}
	return area ;
}

double getDamagedArea(std::vector<DelaunayTriangle *> & paste)
{
	double area = 0 ;
	for(size_t i = 0 ; i < paste.size() ; i++)
	{
		if(paste[i]->getBehaviour()->getDamageModel())
		{
			area += (paste[i]->area())*paste[i]->getBehaviour()->getDamageModel()->getState()[0] ;
		}
	}
	return area ;
}

double getCrackedArea(std::vector<DelaunayTriangle *> & paste)
{
	double area = 0 ;
	for(size_t i = 0 ; i < paste.size() ; i++)
	{
		if(paste[i]->getBehaviour()->fractured())
		{
			area += (paste[i]->area()) ;
		}
	}
	return area ;
}

int main(int argc, char *argv[])
{
	double timeScale = 1000 ;
	double tau = 10 ;
	int nzones = atof(argv[1]) ;
	double stress = 0 ;	
	if(argc == 3)
		stress = atof(argv[2])*(-1e6) ;

	FeatureTree F(&box) ;
	F.setSamplingNumber(450) ;
	F.setMaxIterationsPerStep(50000) ;
	F.setOrder(LINEAR) ;
	F.setDeltaTime(tau) ;
	
	box.setBehaviour( new ElasticOnlyPasteBehaviour() ) ;
	std::vector<Feature *> incs = ParticleSizeDistribution::get2DConcrete( &F, new ElasticOnlyAggregateBehaviour(), 6000, 0.008, 0., BOLOME_A) ;
	std::vector<Inclusion *> inclusions ;
	for(size_t i = 0 ; i < incs.size() ; i++)
		inclusions.push_back(dynamic_cast<Inclusion *>(incs[i]));
	
	double aggregate_area = 0 ;
	for(size_t i = 0 ; i < inclusions.size() ; i++)
		aggregate_area += inclusions[i]->area() ;
	
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
//	std::cout << 0. << "\t" << x[0] << "\t" << y[0] << std::endl ;

	if(stress != 0)
		F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition( SET_STRESS_ETA, TOP, stress) ) ;

	size_t i = 0 ;
	double time = 0. ;
	std::vector<std::pair<ExpansiveZone *, Inclusion*> > zones = ParticleSizeDistribution::get2DExpansiveZonesInAggregates( &F, inclusions, new GelBehaviour(), 0.00001, nzones*5, nzones) ;
	F.getAssembly()->setEpsilon(1e-15) ;
	F.step() ;

		TriangleWriter writer("hey", &F) ;
		writer.getField(TWFT_STRAIN) ;
		writer.getField(TWFT_STRESS) ;
		writer.getField(TWFT_DAMAGE) ;
		writer.getField(TWFT_STIFFNESS) ;
		writer.write() ;
	
	
	Vector x = F.getAverageField(STRAIN_FIELD) ;
	Vector y = F.getAverageField(REAL_STRESS_FIELD) ;
	Vector z = F.getAverageField(REAL_STRESS_FIELD) ;
	
	std::vector<DelaunayTriangle *> paste = F.getFeature(0)->getElements2D(&F) ;
	std::vector<DelaunayTriangle *> all = F.getElements2D() ;
	
	std::fstream out ;
	std::string toto = "rag_" ;
	toto.append(argv[1]) ;
	toto.append("_") ;
	if(argc == 3)
		toto.append(argv[2]) ;
	out.open(toto.c_str(), std::ios::out ) ;
	
	
	x = F.getAverageField(STRAIN_FIELD) ;
	y = F.getAverageField(REAL_STRESS_FIELD, paste) ;
	z = F.getAverageField(REAL_STRESS_FIELD) ;
	out << time << "\t"                          // 1
	<< reaction(zones) << "\t"                   // 2
	<< aggregate_area << "\t"                    // 3
	<< getDamagedAggregateArea(zones,&F) << "\t" // 4
	<< x[0] << "\t"                              // 5
	<< x[1] << "\t"                              // 6
	<< y[0] << "\t"                              // 7
	<< y[1] << "\t"                              // 8
	<< z[0] << "\t"                              // 9
	<< z[1] << "\t"                              // 10
	<<  getDamagedArea(paste) << "\t"            // 11
	<<  getDamagedArea(all) << "\t"              // 12 
	<< getCrackedArea(paste) << "\t"             // 13
	<< getCrackedArea(all) << std::endl ;        // 14


	double damage = 0. ;
	int notConverged = 0 ;	
	
	do
	{
		i++ ;
 		F.setDeltaTime( tau ) ;
		time += tau  ;

		notConverged += !F.solverConverged() ;
		if(time > timeScale || notConverged > 4)
			return 0 ;
		std::string tati = toto ;
		tati.append("_") ;
		tati.append(itoa(i)) ;
		if(!F.solverConverged())
			tati.append("_notconverged") ;
		std::cout << tati << std::endl ;
		TriangleWriter writer(tati, &F) ;
		writer.getField(TWFT_STRAIN) ;
		writer.getField(TWFT_STRESS) ;
		writer.getField(TWFT_DAMAGE) ;
		writer.getField(TWFT_STIFFNESS) ;
		writer.write() ;
			
		updateZones( zones, getReactiveSurface(zones)*0.05*time/timeScale) ;
		F.step() ;
		x = F.getAverageField(STRAIN_FIELD) ;
		y = F.getAverageField(REAL_STRESS_FIELD, paste) ;
		z = F.getAverageField(REAL_STRESS_FIELD) ;

		damage = reaction(zones) + getDamagedArea(all)-getDamagedArea(paste) ;

		out << time << "\t" << reaction(zones) << "\t" << aggregate_area << "\t" << damage << "\t" << x[0] << "\t" << x[1] << "\t" << y[0] << "\t" << y[1] << "\t" << z[0] << "\t" << z[1] <<  "\t" << getDamagedArea(paste) << "\t" <<  getDamagedArea(all) << "\t" << getCrackedArea(paste) << "\t" << getCrackedArea(all) << std::endl ;
	} while (time < timeScale) ;
		
	return 0 ;
}

