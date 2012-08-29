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

Sample box(nullptr, 1.,1.,0.5,0.5) ;

double young = 1.*1e9 ;
double nu = 0.3 ;
double eta = 1. ;
double stress = 1.*1e6 ;
double imposed = 1e-3 ;

double tau = 0.1 ;
int sampling = 100 ;

Matrix C(3,3) ;
Matrix E(3,3) ;
Matrix K(3,3) ;


int main(int argc, char *argv[])
{
	C[0][0] = 1. ; C[1][0] = nu ; C[2][0] = 0. ;
	C[0][1] = nu ; C[1][1] = 1. ; C[2][1] = 0. ;
	C[0][2] = 0. ; C[1][2] = 0. ; C[2][2] = 1.-nu ;
	C *= young/(1.-nu*nu) ;
	
	E = C*eta ;
	
	K = C*3 ;
	
	FeatureTree F(&box) ;
	F.setSamplingNumber(sampling) ;
	F.setOrder(LINEAR_TIME_LINEAR) ;
	F.setDeltaTime(tau) ;

	F.addBoundaryCondition(new TimeContinuityBoundaryCondition()) ;

	Matrix Ci = C*10 ;
	Inclusion * inc = new Inclusion( 0.3,0.5,0.5 ) ;
	inc->setBehaviour(new GeneralizedSpaceTimeViscoelasticity( PURE_ELASTICITY, Ci, 1 ) ) ;
	F.addFeature( &box, inc ) ;
	
// 	box.setBehaviour(new GeneralizedSpaceTimeViscoelasticity(PURE_ELASTICITY, C)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM, 0, 1)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT, 0, 0)) ;
// 
// 	box.setBehaviour(new GeneralizedSpaceTimeViscoelasticity(PURE_VISCOSITY, E)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM, 0, 1)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT, 0, 0)) ;
// 	
// 	box.setBehaviour(new GeneralizedSpaceTimeViscoelasticity(KELVIN_VOIGT, C, E)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM, 0, 1)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT, 0, 0)) ;
// 
// 	box.setBehaviour(new GeneralizedSpaceTimeViscoelasticity(MAXWELL, C, E)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM, 0, 1)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM, 0, 3)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT, 0, 0)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT, 0, 2)) ;
// 		
// 	box.setBehaviour(new GeneralizedSpaceTimeViscoelasticity(BURGER, C, E, C, E)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM, 0, 1)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM, 0, 3)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM, 0, 5)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT, 0, 0)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT, 0, 2)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT, 0, 4)) ;
// 
	box.setBehaviour(new GeneralizedSpaceTimeViscoelasticity(GENERALIZED_KELVIN_VOIGT, C, C, E)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM, 0, 1)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM, 0, 3)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT, 0, 0)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT, 0, 2)) ;
// 	
// 	box.setBehaviour(new GeneralizedSpaceTimeViscoelasticity(GENERALIZED_MAXWELL, K, C, E)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM, 0, 1)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM, 0, 3)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT, 0, 0)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT, 0, 2)) ;

	srand(0) ;
	F.step() ;

	Function stressRamp("t") ;
	stressRamp = stressRamp + tau/2 ;
	stressRamp = stressRamp/tau ;
	stressRamp = stressRamp * stress ;
	Function stressConstant("1") ;
	stressConstant = stressConstant * stress ;
	Function stressNull("1") ;
	stressNull = stressNull * 0. ;
	BoundingBoxDefinedBoundaryCondition * stressBCRamp = new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP_AFTER, stressRamp) ;
	BoundingBoxDefinedBoundaryCondition * stressBCConstant = new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP_AFTER, 0.) ;
//  	F.addBoundaryCondition(stressBCRamp) ;
//  	F.addBoundaryCondition(stressBCConstant) ;
		
	BoundingBoxDefinedBoundaryCondition * imposedBC = new BoundingBoxDefinedBoundaryCondition( SET_ALONG_ETA, TOP_AFTER, imposed) ;
	F.addBoundaryCondition(imposedBC) ;
	
	int imax = 20 ; 
	
	for(size_t i = 0 ; i < imax ; i++)
	{
		F.step() ;
		Vector x = F.getDisplacements() ;
		
		TriangleWriter writer( "test_"+itoa(i), &F, 1) ;
		writer.getField(TWFT_STRESS) ;
		writer.getField(TWFT_STRAIN) ;
		writer.write() ;

		if(i == 0)
		{
			stressBCRamp->setData( stressNull ) ;
			stressBCConstant->setData( stress ) ;
		}		
	}
	
	return 0 ;
}

