// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../utilities/samplingcriterion.h"
#include "../features/features.h"
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
	Function x("x") ;
	Function y("y") ;
	Function t("t") ;
	Function one("1") ;
	Function zero("0") ;
	x.setNumberOfDerivatives(4);
	y.setNumberOfDerivatives(4);
	t.setNumberOfDerivatives(4);
	x.setDerivative(XI,one) ;
	y.setDerivative(ETA,one) ;
	t.setDerivative(TIME_VARIABLE,one) ;
	y.setDerivative(XI,zero) ;
	x.setDerivative(ETA,zero) ;
	y.setDerivative(ZETA,zero) ;
	x.setDerivative(ZETA,zero) ;
	y.setDerivative(TIME_VARIABLE,zero) ;
	x.setDerivative(TIME_VARIABLE,zero) ;
	t.setDerivative(XI,zero) ;
	t.setDerivative(ETA,zero) ;
	t.setDerivative(ZETA,zero) ;
  
  
	Vector toto ;
	FeatureTree F(&box) ;
	F.setSamplingNumber(2) ;
	Matrix e = (new ElasticOnlyPasteBehaviour(10e9, 0.3))->param ;
  	box.setBehaviour(new Stiffness(e)) ;
	F.setOrder(LINEAR_TIME_LINEAR) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0,0)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0,1)) ;
	F.step() ;
	
	toto = F.getDisplacements() ;
	
	DelaunayTriangle* tri = F.getElements2D()[0] ;
	Point c(0.27,0.35) ;
	double radius = 0.2 ;
	Circle inc(radius, c) ;
	if(inc.intersects(dynamic_cast<Triangle *>(tri)))
		std::cout << "intersects" << std::endl ;
	Function position = f_sqrt((tri->getXTransform()-c.x)*(tri->getXTransform()-c.x) +
		                           (tri->getYTransform()-c.y)*(tri->getYTransform()-c.y)) ;
	Function hat = 1.-f_abs(position-radius);
	hat *= (t+1) ;
	
	Function otherPosition = f_sqrt( (y - c.y)*(y - c.y) + (x - c.x)*(x - c.x) ) ;
	Function otherHat = 1.-f_abs( otherPosition - radius ) ;
	otherHat *= t ;
// 	otherHat.setNumberOfVariableTransforms(2);
	Function dx = tri->getXTransform() ;
	Function dy = tri->getYTransform() ;
 	otherHat.setVariableTransform(XI, dx);
 	otherHat.setVariableTransform(ETA, dy);
	
	Function dt = t+1 ;
 	otherHat.setVariableTransform(TIME_VARIABLE, dt);
	
	otherHat.makeVariableTransformDerivative() ;
	
	tri->getBoundingPoint(0).print() ;
	tri->getBoundingPoint(1).print() ;
	tri->getBoundingPoint(2).print() ;
	
	Point p(0,1,0,1) ;
	Point dummy ;
	
	std::cout << otherHat.getNumberOfDerivatives() << std::endl ;
	std::cout << hat.getNumberOfDerivatives() << std::endl ;
	
	std::cout << VirtualMachine().ddeval( hat, XI, TIME_VARIABLE, p ) << std::endl ;
	std::cout << VirtualMachine().ddeval( otherHat, XI, TIME_VARIABLE, p ) << std::endl ;	
	
	
  
	return 0 ;
  
  
	
// 	Matrix e = (new ElasticOnlyPasteBehaviour(10e9, 0.3))->param ;
  	box.setBehaviour(new Viscoelasticity(PURE_ELASTICITY, e)) ;

	F.setOrder(LINEAR_TIME_LINEAR) ;
	F.setDeltaTime(1.) ;

	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0,0)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0,1)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0,2)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0,3)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0,4)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0,5)) ;
	F.addBoundaryCondition(new TimeContinuityBoundaryCondition()) ;
	F.step() ;
	F.step() ;
	
	Vector instants(2) ;
	instants[0] = 0.5 ;
	instants[1] = 1.5 ;
	
	Function r("t 0.03 *") ;
	TimeDependentEnrichmentInclusion tarata( r,0,0) ;
	F.addFeature(&box, &tarata);
	Mesh<DelaunayTriangle, DelaunayTreeItem> * mesh = F.get2DMesh() ;
	
	for(size_t i = 0 ; i < 10 ; i++)
	{
	  F.step() ;
// 	  instants += 1 ;
// //	  tarata.setTimeCircles(instants);
// 	  std::vector<DelaunayTriangle *> hop = mesh->getConflictingElements(&tarata) ;
// 	  std::cout << tarata.radiusAtTime(Point(0,0,0,F.getCurrentTime())) << "\t" <<  hop.size() << std::endl ;
// 	  for(size_t j = 0 ; j < hop.size() ; j++)
// 	  {
// 		  hop[j]->getBehaviour()->param += e ;
// 	  }
// 	  std::string name = "hop_" ;
// 	  name.append( itoa(i) ) ;
// 	  TriangleWriter wrt(name, &F, 1) ;
// 	  wrt.getField( TWFT_STIFFNESS) ;
// 	  wrt.write() ;
	  
	}
	 
	
	return 0 ;
}
