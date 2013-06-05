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
//#include "../physics/generalized_spacetime_viscoelasticity.h"
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

Sample box(nullptr, 0.08, 0.08,0.0,0.0) ;

int findIndexOfAggregate( DelaunayTriangle * tri, std::vector<Inclusion *> aggregates, FeatureTree * f)
{
	for(size_t i = 0 ; i < aggregates.size() ; i++)
	{
		std::vector<DelaunayTriangle *> aggTriangles = aggregates[i]->getElements2D(f) ;
		if( std::find( aggTriangles.begin(), aggTriangles.end(), tri) != aggTriangles.end() )
			return i ;
	}
  
	return 0 ;
}

int main(int argc, char *argv[])
{
	double timeScale = atof(argv[1]) ;
	double tau = timeScale * 0.01 ;
	double timestep = tau ;
	double appliedLoadEta = atof(argv[2])*(-1e6) ;
	double appliedLoadXi = atof(argv[3])*(-1e6) ;
	
	bool elastic = false ;
	bool noZones = false ;
	
	for(size_t i = 0 ; i < argc ; i++)
	{
		if(std::string(argv[i]) == "--elastic")
			elastic = true ;
		if(std::string(argv[i]) == "--nozones")
			noZones = true ;
	}
	
	if(elastic)
		std::cout << "elastic case" << std::endl ;
	if(noZones)
		std::cout << "disabling asr zones" << std::endl ;
	
	int nzones = 200 ;
	int naggregates = 500 ;

	FeatureTree F(&box) ;
	F.setSamplingNumber(400) ;
	F.setOrder(LINEAR_TIME_LINEAR) ;
	
	F.setDeltaTime(tau) ;
	
// 	box.setBehaviour( new ViscoElasticOnlyPasteBehaviour() ) ;
// 	if(elastic)
//		box.setBehaviour(new ViscoElasticOnlyPasteBehaviour(12e9,0.3,0,0) );
	
	Form * stAggregateBehaviour = new ViscoElasticOnlyAggregateBehaviour() ;
	Form * stPasteBehaviour = new ViscoElasticOnlyPasteBehaviour() ;
	ViscoelasticityAndImposedDeformation * stGelBehaviour = dynamic_cast<ViscoelasticityAndImposedDeformation *>((new ViscoElasticOnlyGelBehaviour())->getCopy()) ;
	
	Form * elasticAggregateBehaviour = new ElasticOnlyAggregateBehaviour() ;
	Form * elasticPasteBehaviour = new ElasticOnlyPasteBehaviour() ;
	GelBehaviour * elasticGelBehaviour = new GelBehaviour() ;
	
	Form * stElasticAggregateBehaviour = new Viscoelasticity( PURE_ELASTICITY, elasticAggregateBehaviour->param) ;
	Form * stElasticPasteBehaviour = new Viscoelasticity( PURE_ELASTICITY, elasticPasteBehaviour->param) ;
	ViscoelasticityAndImposedDeformation * stElasticGelBehaviour = new ViscoelasticityAndImposedDeformation( PURE_ELASTICITY, elasticGelBehaviour->param, elasticGelBehaviour->imposed ) ;
	
	Form * aggBehaviour = stAggregateBehaviour ;
	Form * pasteBehaviour = stPasteBehaviour ;
	ViscoelasticityAndImposedDeformation * gelBehaviour = stGelBehaviour ;
	
	if(elastic)
	{
		aggBehaviour = stElasticAggregateBehaviour ;
		pasteBehaviour = stElasticPasteBehaviour ;
		gelBehaviour = stElasticGelBehaviour ;
	}
	
	box.setBehaviour( pasteBehaviour );
	std::vector<Feature *> feats = ParticleSizeDistribution::get2DConcrete( &F, aggBehaviour, naggregates, 0.008, 0.0001, BOLOME_A) ;
	
	std::vector<Inclusion *> aggregates ;
	for(size_t i = 0 ; i < feats.size() ; i++)
	      aggregates.push_back( dynamic_cast<Inclusion *>(feats[i]) ) ;
	
	double aggregatesArea = 0 ; 
	for(size_t i = 0 ; i < aggregates.size() ; i++)
		aggregatesArea += aggregates[i]->area() ;
	double maxGelArea = aggregatesArea*0.03/nzones ;
	double maxGelRadius = sqrt(maxGelArea/M_PI) ;
	Function radius("t") ;
	radius *= (maxGelRadius / timeScale) ;
	
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0, 0)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0, 1)) ;
	if(!elastic)
	{
		F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0, 2)) ;
		F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0, 3)) ;
		F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0, 4)) ;
		F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0, 5)) ;
	}
	
	F.step() ;
	F.step() ;
	
	if(appliedLoadEta < 0)
		F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition( SET_STRESS_ETA, TOP_AFTER, appliedLoadEta));

	if(appliedLoadXi < 0)
		F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition( SET_STRESS_XI, RIGHT_AFTER, appliedLoadXi));

	
	if(!noZones)
		std::vector<std::pair<GrowingExpansiveZone *, Inclusion*> > zones = ParticleSizeDistribution::get2DGrowingExpansiveZonesInAggregates( &F, aggregates, gelBehaviour, radius, maxGelRadius, nzones*50, nzones) ;
	
	std::string name = "asr_creep/asr_creep_no_damage_" ;
	if(elastic)
		name.append("elastic_") ;
	if(noZones)
		name.append("nozones_") ;
	name.append(argv[1]) ;
	name.append("_") ;
	name.append(argv[2]) ;
	name.append("_") ;
	name.append(argv[3]) ;

	std::string nameshort = name ;
	nameshort.append("_summary") ;
	
	std::ofstream summary ;
	summary.open(nameshort.c_str(), std::ios::out) ;
	
	std::vector<DelaunayTriangle *> triangles ;
	std::map<DelaunayTriangle *, int> map ;
	for(size_t i = 0 ; i < 15 ; i++)
	{
		F.setDeltaTime(tau) ;
		F.step() ;
		
		if(triangles.size() == 0)
		{
			triangles = F.getElements2D() ;
			for(size_t k = 0 ; k < aggregates.size() ; k++)
			{
				std::vector<DelaunayTriangle *> trgAggregates = aggregates[k]->getElements2D( &F ) ;
				for(size_t j = 0 ; j < trgAggregates.size() ; j++)
				{
					map[ trgAggregates[j] ] = k+1 ;
				}
			}
			
			std::string namemap = name ;
			namemap.append("_map") ;
			
			std::ofstream mapstream ;
			mapstream.open(namemap.c_str(), std::ios::out ) ;
					
			for(size_t k = 0 ; k < triangles.size() ; k++)
			{
				mapstream << k << "\t" << map[ triangles[k] ] << "\t" << triangles[k]->area() << std::endl ;
			}
			mapstream.close() ;
			
		}
		
		summary << F.getCurrentTime() << "\t" ;
		double r = VirtualMachine().eval( radius, 0,0,0, F.getCurrentTime() ) ;
		summary << r << "\t" << nzones*r*r*M_PI/aggregatesArea << "\t" ;
		
		Vector strain = F.getAverageField(STRAIN_FIELD) ;
		Vector stress = F.getAverageField(REAL_STRESS_FIELD) ;
		summary << strain[0] << "\t" << strain[1] << "\t" << strain[2] << "\t" ;
		summary << stress[0] << "\t" << stress[1] << "\t" << stress[2] << std::endl ;
		
		std::string namei = name ;
		namei.append("_") ;
		namei.append( itoa(i) ) ;
		
		std::ofstream out ;
		out.open(namei.c_str(), std::ios::out ) ;
		
		for(size_t k = 0 ; k < triangles.size() ; k++)
		{
			triangles[k]->getState().getAverageField( STRAIN_FIELD, strain ) ;
			triangles[k]->getState().getAverageField( REAL_STRESS_FIELD, stress ) ;
			out << k << "\t" << strain[0] << "\t" << strain[1] << "\t" << strain[2] << "\t" ;
			out << stress[0] << "\t" << stress[1] << "\t" << stress[2] << std::endl ;
		}
		out.close() ;
		
		tau += timestep ;
	}
	
	
	return 0 ;
}

