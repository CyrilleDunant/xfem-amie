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
// 	Function tata("cst", 2) ;
// 	std::cout << VirtualMachine().eval(tata, 0,1,2,3) << std::endl ;
//   
// 
// 	std::vector<double>  hop ;
// 	hop.push_back(3) ; hop.push_back(1) ;
// 	Function toto("cst cst +", hop) ;
// 	std::cout << VirtualMachine().eval(toto, 0,1,2,3) << std::endl ;
// 
// 	Function tutu("cst cst cst + +", hop) ;
// 	std::cout << VirtualMachine().eval(tutu, 0,1,2,3) << std::endl ;
// 	
// 	exit(0) ;
	
/*	Function x("x") ;
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
	
	
  
	return 0 ;*/
  
  	FeatureTree F(&box) ;
	F.setSamplingNumber(32) ;

	Vector alpha(0.,3) ;
	alpha[0] = 0.1 ;
	alpha[1] = 0.1 ;

	
 	Matrix e = (new ElasticOnlyPasteBehaviour(10e9, 0.3))->param ;
	Matrix e2 = e*0.5 ;
  	box.setBehaviour(new Viscoelasticity(PURE_ELASTICITY, e,0)) ;
//  	box.setBehaviour(new Stiffness(e)) ;
	
	F.setOrder(LINEAR_TIME_LINEAR) ;
	F.setDeltaTime(1.) ;
// 	Inclusion * inc = new Inclusion( 0.35, Point(0,0)) ;
// 	inc->setBehaviour(new Viscoelasticity(PURE_ELASTICITY, e*0.));
// 	F.addFeature(&box, inc);
	
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0,0)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0,1)) ;
// 	F.step() ;
//	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP_AFTER, -1e6,1));
//    	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0,2)) ;
//    	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0,3)) ;
//	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP_AFTER, 1e6));
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0,2)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0,3)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0,4)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0,5)) ;
// 	F.step() ;
// 	F.step() ;
	
	Function r("0.03 t *") ;
	GrowingExpansiveZone *  tarata = new GrowingExpansiveZone( nullptr, r,0,0, new ViscoelasticityAndImposedDeformation( PURE_ELASTICITY, e, alpha,0 )) ;
	tarata->setInitialTime(-2.);
//	ExpansiveZone * tarata = new ExpansiveZone( &box, VirtualMachine().eval(r, 0,0,0,2), 0,0, new StiffnessWithImposedDeformation( e*0.7,alpha)) ;
	Inclusion * taratatata = new Inclusion( 0.1, 0.,0.) ;
	taratatata->setBehaviour( new ViscoelasticityAndImposedDeformation( PURE_ELASTICITY, e, alpha,0 ));
	F.addFeature(&box, tarata);
	std::vector<DelaunayTriangle *> mesh = F.getElements2D() ;
	
	std::fstream out ;
	out.open("test_zone_moving_1", std::ios::out) ;
	
	while(F.getCurrentTime() < 100)
	{
//	  F.setDeltaTime(i+1);
	  F.step() ;
// 	  if(i == 2)
// 		  exit(0) ;
	  
	  Vector str = F.getAverageField(STRAIN_FIELD, -1, 1) ;
	  Vector stress = F.getAverageField(REAL_STRESS_FIELD, -1, 1) ;
	  out << F.getCurrentTime() << "\t" << tarata->radiusAtTime(Point(0,0,0,F.getCurrentTime())) << "\t" << str[0] << "\t" << str[1] << "\t" << str[2] << "\t" << stress[0] << "\t" << stress[1] << "\t" << stress[2] << std::endl ;
// 	  std::cout << "finish\t" << F.getDisplacements(-1, false).size() << std::endl ;
  	  std::cout << F.getCurrentTime() << "\t"  << str[0] << "\t" << str[1] << std::endl ;
// 	  instants += 1 ;
// //	  tarata.setTimeCircles(instants);
// 	  std::vector<DelaunayTriangle *> hop = mesh->getConflictingElements(&tarata) ;
// 	  std::cout << tarata.radiusAtTime(Point(0,0,0,F.getCurrentTime())) << "\t" <<  hop.size() << std::endl ;
// 	  for(size_t j = 0 ; j < hop.size() ; j++)
// 	  {
// 		  hop[j]->getBehaviour()->param += e ;
// 	  }
	  std::string name = "hop_stfem_" ;
	  name.append( itoa(F.getCurrentTime()) ) ;
	  TriangleWriter wrt(name, &F, 1) ;
	  wrt.getField( STRAIN_FIELD ) ;
	  wrt.getField( REAL_STRESS_FIELD ) ;
	  wrt.getField( TWFT_STIFFNESS ) ;
	  wrt.write() ;
	  
/*	  int cIn = 0 ;
	  int cRing = 0 ;
	  int cRing6 = 0 ;
	  int cOut = 0 ;
	  for(size_t j = 0 ; j < mesh.size() ; j++)
	  {
		  if(mesh[j]->getEnrichmentFunctions().size() > 0)
		  {
			  cRing++ ;
			  if(mesh[j]->getEnrichmentFunctions().size() == 6)
			  {
				  cRing6++ ;
			  }
				  
		  }
		  else
		  {
//			  mesh[j]->getBoundingPoint(0).print() ;
			  if(tarata->getPrimitive()->in( mesh[j]->getBoundingPoint(0) ) )
			  {
				  cIn ++ ;
			  }
			  else
			  {
				  cOut++ ;
			  }
		  }
		  
	  }
	  
	  
//	  tarata->setRadius(VirtualMachine().eval(r, 0,0,0,F.getCurrentTime()+1));*/
	  
	}
	
	F.getNodes()[0]->print() ;
	 
	
	return 0 ;
}
