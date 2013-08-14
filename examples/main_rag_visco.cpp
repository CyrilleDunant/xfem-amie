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
#include "../physics/dual_behaviour.h"
#include "../physics/parallel_behaviour.h"
//#include "../physics/generalized_spacetime_viscoelasticity.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../features/pore.h"
#include "../utilities/writer/triangle_writer.h"
#include "../physics/materials/gel_behaviour.h"
#include "../physics/orthotropicstiffness.h"
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
	omp_set_num_threads(4) ;

	double timeScale = atof(argv[1]) ;
	double tau = timeScale * 0.01 ;
	double timestep = tau ;
	double appliedLoadEta = atof(argv[2])*(-1e6) ;
	double appliedLoadXi = atof(argv[3])*(-1e6) ;
		
	int nzones = 100 ;
	int naggregates = 1000 ;

	FeatureTree F(&box) ;
	F.setSamplingNumber(160) ;
	F.setOrder(LINEAR_TIME_LINEAR) ;
	F.setDeltaTime(tau) ;
	F.setMinDeltaTime(tau*1e-6) ;
	F.setMaxIterationsPerStep(5000) ;
		
//	box.setBehaviour( new ViscoElasticOnlyPasteBehaviour() );
	box.setBehaviour( new ViscoDamagePasteBehaviour() );
	std::vector<Feature *> feats = ParticleSizeDistribution::get2DConcrete( &F, new ViscoDamageAggregateBehaviour(), naggregates, 0.008, 0.0001, BOLOME_A) ;
	
	std::vector<Inclusion *> aggregates ;
	for(size_t i = 0 ; i < feats.size() ; i++)
	      aggregates.push_back( dynamic_cast<Inclusion *>(feats[i]) ) ;
	
	double aggregatesArea = 0 ; 
	for(size_t i = 0 ; i < aggregates.size() ; i++)
		aggregatesArea += aggregates[i]->area() ;
	double maxGelArea = aggregatesArea*0.05/nzones ;
	double maxGelRadius = sqrt(maxGelArea/M_PI) ;
	Function radius("t") ;
	radius *= (maxGelRadius / timeScale) ;
	
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0, 0)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0, 1)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0, 2)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0, 3)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0, 4)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0, 5)) ;
	
	if(appliedLoadEta < 0)
		F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition( SET_STRESS_ETA, TOP_AFTER, appliedLoadEta));

	if(appliedLoadXi < 0)
		F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition( SET_STRESS_XI, RIGHT_AFTER, appliedLoadXi));

	ViscoelasticityAndImposedDeformation * gel = dynamic_cast<ViscoelasticityAndImposedDeformation *>( (new ViscoElasticOnlyGelBehaviour())->getCopy() ) ;
	
	std::vector<std::pair<TimeDependentHomogenisingInclusion *, Inclusion*> > zones = ParticleSizeDistribution::get2DGrowingExpansiveZonesInAggregates( &F, aggregates, gel, radius, maxGelRadius, nzones*5, nzones) ;
	
	std::string name = "asr_creep_full_" ;
	name.append(argv[1]) ;
	name.append("_") ;
	name.append(argv[2]) ;
	name.append("_") ;
	name.append(argv[3]) ;

	std::ofstream summary ;
	summary.open(name.c_str(), std::ios::out) ;
	
	for(size_t i = 0 ; i < 15 ; i++)
	{
		F.setDeltaTime(tau) ;
		F.setMinDeltaTime(tau*1e-6) ;
		F.step() ;

		summary << F.getCurrentTime() << "\t" ;
		double r = VirtualMachine().eval( radius, 0,0,0, F.getCurrentTime() ) ;
		summary << r << "\t" << nzones*r*r*M_PI/aggregatesArea << "\t" ;
		
		Vector strain = F.getAverageField(STRAIN_FIELD) ;
		Vector stress = F.getAverageField(REAL_STRESS_FIELD) ;
		summary << strain[0] << "\t" << strain[1] << "\t" << strain[2] << "\t" ;
		summary << stress[0] << "\t" << stress[1] << "\t" << stress[2] << "\t" ;
		summary << F.averageDamage << "\t" << std::endl ;
		
		if(i == 0 || (i+1)%5 == 0)
		{
			std::string nametrg = name ;
			nametrg.append("_trg_") ;
			nametrg.append(itoa(i)) ;
			TriangleWriter w( nametrg, &F, 1) ;
			w.getField( STRAIN_FIELD ) ;
			w.getField( REAL_STRESS_FIELD ) ;
			w.getField( TWFT_STIFFNESS ) ;
			w.getField( TWFT_DAMAGE ) ;
			w.write() ;
		}
		
		tau += timestep ;
	}
	
	
	return 0 ;
}

