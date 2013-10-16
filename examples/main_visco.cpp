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

Sample box(nullptr, 0.17, 0.17,0.0,0.0) ;

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

double avgTensileStressInPaste( std::vector<DelaunayTriangle *> trg)
{
	double ar = 0. ;
	double tension = 0. ;
	for(size_t i = 0 ; i < trg.size() ; i++)
	{
		Viscoelasticity * visc = dynamic_cast<Viscoelasticity *>(trg[i]->getBehaviour()) ;
		if(visc->model == GENERALIZED_KELVIN_VOIGT)
		{
			Vector stress(2) ; stress = 0 ;
			trg[i]->getState().getAverageField( PRINCIPAL_REAL_STRESS_FIELD, stress ) ;
			if(stress.max() > 0)
			{
				double a = trg[i]->area() ;
				tension += stress.max()*a ;
				ar += a ;
			}
		}
	}
	if(ar > 0)
		return tension/ar ;
	return 0 ;
}

double areaTensionInPaste( std::vector<DelaunayTriangle *> trg)
{
	double ar = 0. ;
	for(size_t i = 0 ; i < trg.size() ; i++)
	{
		Viscoelasticity * visc = dynamic_cast<Viscoelasticity *>(trg[i]->getBehaviour()) ;
		if(visc->model == GENERALIZED_KELVIN_VOIGT)
		{
			Vector stress(2) ; stress = 0 ;
			trg[i]->getState().getAverageField( PRINCIPAL_REAL_STRESS_FIELD, stress ) ;
			if(stress.max() > 0)
			{
				ar = trg[i]->area() ;
			}
		}
	}
	return ar ;
}

double maxTensileStressInPaste( std::vector<DelaunayTriangle *> trg)
{
	double tension = 0. ;
	for(size_t i = 0 ; i < trg.size() ; i++)
	{
		Viscoelasticity * visc = dynamic_cast<Viscoelasticity *>(trg[i]->getBehaviour()) ;
		if(visc->model == GENERALIZED_KELVIN_VOIGT)
		{
			Vector stress(2) ; stress = 0 ;
			trg[i]->getState().getAverageField( PRINCIPAL_REAL_STRESS_FIELD, stress ) ;
			if(stress.max() > tension)
				tension = stress.max() ;
		}
	}
	return tension ;
}

Vector dispOnTop( Vector disp, std::vector<Point *> nodes)
{
	Vector top(6) ;
	top = 0. ;
	size_t count = 0 ;
	for(size_t i = 0 ; i < nodes.size() ; i++)
	{
		if(nodes[i]->t > nodes[0]->t && nodes[i]->y >= 0.0849)
		{
			for(size_t j = 0 ; j < 6 ; j++)
			{
				top[j] += disp[ nodes[i]->id*6 + j ] ;
			}
			count++ ;
		}
	}
	top /= count ;
	return top ;
}

int main(int argc, char *argv[])
{
//	omp_set_num_threads(1) ;

	double timeScale = 1000. ;//atof(argv[1]) ;
	double tau = 1. ;
	double timestep = tau ;
	double appliedLoadEta = -1e6 *atof(argv[2]);
	double appliedLoadXi = 0 ;//atof(argv[3])*(-1e6) ;
	double degree = 0.01 ;//atof(argv[1]) ;
	double maxstress = atof(argv[1])*1e6 ;
//	tau /= degree/0.001 ;
		
	double limit = 0.1 ;//atof(argv[1]) ;

	int nzones = 400 ;
	int naggregates = 2000 ;

	FeatureTree F(&box) ;
	F.setSamplingNumber(196) ;
	F.setOrder(LINEAR_TIME_LINEAR) ;
	F.setDeltaTime(tau) ;
	F.setMinDeltaTime(tau*1e-5) ;
	F.setMaxIterationsPerStep(5000) ;
		
	PseudoBurgerViscoDamagePasteBehaviour paste(12e9, 0.3, 2, 10000, maxstress/12e9) ;
	paste.ctype = STRESS_CRITERION ;

	box.setBehaviour( &paste );

	std::vector<Feature *> feats = ParticleSizeDistribution::get2DConcrete( &F, new ViscoElasticOnlyAggregateBehaviour(), naggregates, 0.015, 0.0001, BOLOME_A) ;	
	
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_LEFT_AFTER, 0, 0)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0, 1)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_LEFT_AFTER, 0, 2)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0, 3)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_LEFT_AFTER, 0, 4)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0, 5)) ;
	
	if(appliedLoadEta < 0)
		F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition( SET_STRESS_ETA, TOP_AFTER, appliedLoadEta));

	if(appliedLoadXi < 0)
		F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition( SET_STRESS_XI, RIGHT_AFTER, appliedLoadXi));

	
	std::string name = "creep_nonlinear_" ;
	name.append(argv[1]) ;
	name.append("_") ;
	name.append(argv[2]) ;
	std::ofstream summary ;
	summary.open(name.c_str(), std::ios::out) ;

	std::vector<DelaunayTriangle *> trg ;
	std::vector<Point *> nodes ;
	size_t i = 0 ;
	while(F.getCurrentTime() < timeScale)
	{
		i++ ;
		F.setDeltaTime(tau) ;
		F.setMinDeltaTime(tau*1e-6) ;
		F.step() ;

		F.getAssembly()->setEpsilon(1e-8) ;
		if(trg.size() == 0)
		{
			trg = F.getElements2D() ;
			nodes = F.getNodes() ;
		}
		
		Vector strain = F.getAverageField(STRAIN_FIELD, -1, 1) ;
		Vector stress = F.getAverageField(REAL_STRESS_FIELD, -1, 1) ;
		Vector d =  dispOnTop( F.getDisplacements(), nodes) ;
		summary << F.getCurrentTime() << "\t" << strain[1] << "\t" << stress[1] << "\t" ; 
		for(size_t j = 0 ; j < d.size() ; j++)
			summary << d[j] << "\t" ;
		summary <<  areaTensionInPaste( trg ) << "\t" << avgTensileStressInPaste( trg ) << "\t" << maxTensileStressInPaste( trg ) << std::endl ;
		
//		if(F.getCurrentTime() > 20)
//			timestep *= 10. ;
		tau += timestep ;

		if(i%5 == 0)
		{
			std::string tati = name ;
			tati.append("_pattern") ;
			tati.append(itoa(i)) ;
			TriangleWriter writer(tati, &F, 1) ;
	 		writer.getField(STRAIN_FIELD) ;
	 		writer.getField(REAL_STRESS_FIELD) ;
			writer.getField(TWFT_DAMAGE) ;
			writer.getField(TWFT_STIFFNESS) ;
			writer.write() ;
		}
	}


	


	return 0 ;
}

