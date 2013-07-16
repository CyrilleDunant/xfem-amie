// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../utilities/samplingcriterion.h"
#include "../features/features.h"
#include "../features/growingExpansiveZone.h"
#include "../features/timeDependentEnrichmentInclusion.h"
#include "../physics/physics_base.h"
#include "../physics/kelvinvoight.h"
#include "../physics/maxwell.h"
#include "../physics/stiffness.h"
#include "../physics/parallel_behaviour.h"
#include "../physics/viscoelasticity.h"
#include "../physics/viscoelasticity_and_fracture.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../features/pore.h"
#include "../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h"
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
#include "../geometry/space_time_geometry_2D.h"

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

Sample box(nullptr, 1.,1.,0.,0.) ;


int main(int argc, char *argv[])
{	
  
  	FeatureTree F(&box) ;
	F.setSamplingNumber(32) ;

	Vector alpha(3) ;
	alpha[0] = 0.0 ;
	alpha[1] = 0.0 ;

	
 	Matrix e = (new ElasticOnlyPasteBehaviour(60e9, 0.3))->param ;
	Matrix e2 = (new ElasticOnlyPasteBehaviour(22e9, 0.3))->param ;
  	box.setBehaviour(new Viscoelasticity(PURE_ELASTICITY, e,0)) ;
	
	F.setOrder(LINEAR_TIME_LINEAR) ;
	F.setDeltaTime(0.1) ;
	
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0,0)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0,1)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP_AFTER, -1e6,1)) ;
	
	Function r("0.05 t *") ;
	r += 0.05 ;
	GrowingExpansiveZone *  tarata = new GrowingExpansiveZone( nullptr, r,0,0, new ViscoelasticityAndImposedDeformation( PURE_ELASTICITY, e2, alpha,0 )) ;
//	tarata->setInitialTime(-2.);
	Inclusion * taratatata = new Inclusion( 0.1, 0.,0.) ;
	taratatata->setBehaviour( new ViscoelasticityAndImposedDeformation( PURE_ELASTICITY, e, alpha,0 ));
	F.addFeature(&box, tarata);
	F.step() ;
	{
	  Vector str = F.getAverageField(STRAIN_FIELD, -1, 0.) ;
	  Vector stress = F.getAverageField(REAL_STRESS_FIELD, -1, 0.) ;
	  std::cout << F.getCurrentTime() << "\t" << /*tarata->radiusAtTime(Point(0,0,0,F.getCurrentTime())) << "\t" <<*/ str[0] << "\t" << str[1] << "\t" << str[2] << "\t" << stress[0] << "\t" << stress[1] << "\t" << stress[2] << std::endl ;
	  
	}
//	std::vector<DelaunayTriangle *> mesh = F.getElements2D() ;
	
	std::fstream out ;
	out.open("test_zone_moving_", std::ios::out) ;
	
	while(F.getCurrentTime() < 20)
	{
	  F.step() ;
	  
	  Vector str = F.getAverageField(STRAIN_FIELD, -1, 0.) ;
	  Vector stress = F.getAverageField(REAL_STRESS_FIELD, -1, 0.) ;
	  std::cout << F.getCurrentTime() << "\t" << tarata->radiusAtTime(Point(0,0,0,F.getCurrentTime())) << "\t" << str[0] << "\t" << str[1] << "\t" << str[2] << "\t" << stress[0] << "\t" << stress[1] << "\t" << stress[2] << std::endl ;
	  out << F.getCurrentTime() << "\t" << tarata->radiusAtTime(Point(0,0,0,F.getCurrentTime())) << "\t" << str[0] << "\t" << str[1] << "\t" << str[2] << "\t" << stress[0] << "\t" << stress[1] << "\t" << stress[2] << std::endl ;
	}
	
	F.getNodes()[0]->print() ;
	 
	
	return 0 ;
}
