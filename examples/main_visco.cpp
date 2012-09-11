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

FeatureTree * featureTree ;
std::vector<DelaunayTriangle *> triangles ;

Sample box(nullptr, 0.08, 0.16,0.0,0.0) ;

double young = 1.25e9 ;
double eta = 100. ;
double stress = -10e6 ;

double tau = 15./(60*24) ;
int sampling = 200 ;
double dmax = 0.008 ;

Matrix C0(3,3) ;
Matrix Cagg(3,3) ;
Matrix C(3,3) ;
Matrix E(3,3) ;


int main(int argc, char *argv[])
{
	C.array() = 0. ;
	C[0][0] = 1. ;
	C[1][1] = 1. ;
	C[2][2] = 1. ;

	Matrix C1(3,3) ;
	Matrix C2(3,3) ;
	Matrix C3(3,3) ;
	
	C1 = C*1.2e9 ; 
	C2 = C*750e9 ;
	C3 = C*1.5e12 ;
	
	Matrix E1(3,3) ;
	Matrix E2(3,3) ;
	Matrix E3(3,3) ;
	
	E1 = C1*1./2 ;
	E2 = C2*1./10 ;
	E3 = C3*1./140 ;
	
	std::vector<std::pair<Matrix, Matrix> > branches ;
	branches.push_back(std::make_pair(C1,E1)) ;
	branches.push_back(std::make_pair(C2,E2)) ;
	branches.push_back(std::make_pair(C3,E3)) ;
		
	ElasticOnlyPasteBehaviour * elasticPaste = new ElasticOnlyPasteBehaviour() ;
	C0 = elasticPaste->param ;
	
	ElasticOnlyAggregateBehaviour * elasticAggregate = new ElasticOnlyAggregateBehaviour() ;
	Cagg = elasticAggregate->param ;
	
	FeatureTree F(&box) ;
	F.setSamplingNumber(sampling) ;
	F.setOrder(LINEAR_TIME_LINEAR) ;
//	F.setOrder(LINEAR) ;
	F.setDeltaTime(tau) ;

	GeneralizedSpaceTimeViscoelasticity * paste = new GeneralizedSpaceTimeViscoelasticity( GENERALIZED_KELVIN_VOIGT, C0, branches ) ;
	GeneralizedSpaceTimeViscoelasticity * paste0= new GeneralizedSpaceTimeViscoelasticity( PURE_ELASTICITY, C0, 0 ) ;
	//GeneralizedSpaceTimeViscoelasticity * paste1= new GeneralizedSpaceTimeViscoelasticity( PURE_ELASTICITY, C, 0 ) ;
	GeneralizedSpaceTimeViscoelasticity * agg   = new GeneralizedSpaceTimeViscoelasticity( PURE_ELASTICITY , Cagg, 3 ) ;
		
	F.addBoundaryCondition(new TimeContinuityBoundaryCondition()) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM, 0, 1)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM, 0, 3)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM, 0, 5)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM, 0, 7)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT, 0, 0)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT, 0, 2)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT, 0, 4)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT, 0, 6)) ;

	Function stressRamp("t") ;
	stressRamp = stressRamp + (tau/2) ;
	stressRamp = stressRamp/tau ;
	stressRamp = stressRamp * stress ;
	
	Function stressConstant("1") ;
	stressConstant = stressConstant * stress ;
	Function stressNull("0") ;
	
	BoundingBoxDefinedBoundaryCondition * stressBCRamp = new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP_AFTER, stressRamp) ;
	BoundingBoxDefinedBoundaryCondition * stressBCConstant = new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP_AFTER, 0.) ;
  	F.addBoundaryCondition(stressBCRamp) ;
//  	F.addBoundaryCondition(stressBCConstant) ;
	
 	box.setBehaviour(paste) ;
//	box.setBehaviour(elasticPaste) ;
	
	Sample * sensor = new Sample(NULL,0.005,0.1,-0.04,0) ;
	sensor->setBehaviour(new VoidForm()) ;

	std::vector<Inclusion *> inclusions = ParticleSizeDistribution::get2DConcrete(0.008, 0.105, 4000 ) ;
	std::vector<Feature *> feats ;
//	feats.push_back(sensor) ;
	for( size_t i = 0; i < inclusions.size() ; i++ )
		feats.push_back( inclusions[i] ) ;
	inclusions.clear() ;
	Rectangle placeGeometry( 0.08, 0.16, 0, 0 ) ;
	int nAgg = 1 ;
	feats = placement( &placeGeometry, feats, &nAgg, 0, 6400 );
	for( size_t i = 0; i < feats.size() ; i++ )
		inclusions.push_back( static_cast<Inclusion *>( feats[i] ) ) ;

//	F.addFeature(&box, sensor) ;
	for( size_t i = 0; i < inclusions.size() ; i++ )
	{
		inclusions[i]->setBehaviour( agg ) ;
//		inclusions[i]->setBehaviour( elasticAggregate ) ;
		F.addFeature(&box, inclusions[i] ) ;
	}
	
	
	srand(0) ;
	int imax = 4*24 ; 
	Vector straino(3) ;
	straino = 0. ;
	
	for(size_t i = 0 ; i < imax ; i++)
	{
		F.step() ;
		Vector x = F.getDisplacements() ;
		
		triangles = F.getElements2D() ;
		Vector strain(3) ;
		Vector bstrain(3) ;
		strain = 0. ; bstrain = 0. ;
		double area = 0. ;
		double xmin = 0. ;
		double xmax = 0. ;
		for(size_t j = 0 ; j < triangles.size() ; j++)
		{
			triangles[j]->getState().getAverageField( STRAIN_FIELD, bstrain ) ;
			strain += bstrain ;
			area += triangles[j]->area() ;
		}
		if(i == 0)
		{
			straino = strain ;
			straino *= 2. ;
		}
		else
		{
			std::cout << (strain[1]-straino[1])/(area*tau) << "\t" << (strain[0]-straino[0])/(area*tau) <<std::endl ;
		}
//		std::cout << x.min()/0.16 << "\t" << x.max()/0.08 << std::endl ;

		
		if(i == 0)
		{
			stressBCRamp->setData( stressConstant ) ;
			stressBCConstant->setData( stress ) ;
		}		
	}
	
// 	TriangleWriter writer( "test", &F, 1) ;
// 	writer.getField(TWFT_STRESS) ;
// 	writer.getField(TWFT_STRAIN) ;
// 	writer.write() ;

	
	
	return 0 ;
}

