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
#include "../physics/viscoelasticity.h"
#include "../physics/viscoelasticity_and_fracture.h"
#include "../physics/viscoelasticity_and_imposed_deformation.h"
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

Sample box(nullptr, 1, 1,0.,0.) ;


int main(int argc, char *argv[])
{
	Function z2("0") ;

	Function z1("0") ;
	z1.setNumberOfDerivatives(4) ;
	for(int i = 0 ; i < 4 ; i++)
	{
		z1.setDerivative( (const Variable) i, z2) ;
	}
  
	Function zero("0") ;
	zero.setNumberOfDerivatives(4) ;
	for(int i = 0 ; i < 4 ; i++)
	{
		zero.setDerivative( (const Variable) i, z1) ;
	}	
	Function one = zero +1 ;
	Function mone = zero-1 ;
  
	Function s0("y") ;
	Function s1("1 x - y -") ;
	Function s2("x") ;
	
	Function t0("t 2 ^ t - 0.5 *") ;			
	Function t1("1 t 2 ^ -") ;
	Function t2("t 2 ^ t + 0.5 *") ;

	Function tt0("t 0.5 -") ;
	Function tt1("t 2 *") ; tt1 *= -1 ;
	Function tt2("t 0.5 +") ;
	
	Function ttt0 = one ;
	Function ttt1 = mone *2 ;
	Function ttt2 = one ;
	
	Function test = s0/tt1 ;
	
	std::cout << VirtualMachine().eval( test, Point(0,1,0,-1) ) << std::endl ;
	
	s0.setNumberOfDerivatives(4) ;
	s0.setDerivative( XI, zero) ;
	s0.setDerivative( ETA, one) ;
	s0.setDerivative( ZETA, zero) ;
	s0.setDerivative( TIME_VARIABLE, zero) ;
	s1.setNumberOfDerivatives(4) ;
	s1.setDerivative( XI, mone) ;
	s1.setDerivative( ETA, mone) ;
	s1.setDerivative( ZETA, zero) ;
	s1.setDerivative( TIME_VARIABLE, zero) ;
	s2.setNumberOfDerivatives(4) ;
	s2.setDerivative( XI, one) ;
	s2.setDerivative( ETA, zero) ;
	s2.setDerivative( ZETA, zero) ;
	s2.setDerivative( TIME_VARIABLE, zero) ;

	tt0.setNumberOfDerivatives(4) ;
	tt0.setDerivative( XI, zero) ;
	tt0.setDerivative( ETA, zero) ;
	tt0.setDerivative( ZETA, zero) ;
	tt0.setDerivative( TIME_VARIABLE, ttt0 ) ;
	tt1.setNumberOfDerivatives(4) ;
	tt1.setDerivative( XI, zero) ;
	tt1.setDerivative( ETA, zero) ;
	tt1.setDerivative( ZETA, zero) ;
	tt1.setDerivative( TIME_VARIABLE, ttt1 ) ;
	tt2.setNumberOfDerivatives(4) ;
	tt2.setDerivative( XI, zero) ;
	tt2.setDerivative( ETA, zero) ;
	tt2.setDerivative( ZETA, zero) ;
	tt2.setDerivative( TIME_VARIABLE, ttt2 ) ;
	
	t0.setNumberOfDerivatives(4) ;
	t0.setDerivative( XI, zero) ;
	t0.setDerivative( ETA, zero) ;
	t0.setDerivative( ZETA, zero) ;
	t0.setDerivative( TIME_VARIABLE, tt0 ) ;
	t1.setNumberOfDerivatives(4) ;
	t1.setDerivative( XI, zero) ;
	t1.setDerivative( ETA, zero) ;
	t1.setDerivative( ZETA, zero) ;
	t1.setDerivative( TIME_VARIABLE, tt1 ) ;
	t2.setNumberOfDerivatives(4) ;
	t2.setDerivative( XI, zero) ;
	t2.setDerivative( ETA, zero) ;
	t2.setDerivative( ZETA, zero) ;
	t2.setDerivative( TIME_VARIABLE, tt2 ) ;
	
	Function test0 = s0/tt1 ;
	
	std::cout << VirtualMachine().eval( test0, Point(0,1,0,-1) ) << std::endl ;
	
	std::cout << VirtualMachine().eval( test0.d(ETA), Point(0,1,0,-1) ) << std::endl ;
	std::cout << VirtualMachine().eval( test0.d(ETA).d(TIME_VARIABLE), Point(0,1,0,-1) ) << std::endl ;
	std::cout << VirtualMachine().eval( test0.d(ETA).d(TIME_VARIABLE).d(TIME_VARIABLE), Point(0,1,0,-1) ) << std::endl ;
exit(0) ;
/*	test = s0-t0 ;
	
	std::cout << VirtualMachine().eval( test, Point(0,1,0,-1) ) << std::endl ;
	std::cout << VirtualMachine().eval( test.d(ETA), Point(0,1,0,-1) ) << std::endl ;
	std::cout << VirtualMachine().eval( test.d(ETA).d(TIME_VARIABLE), Point(0,1,0,-1) ) << std::endl ;
	std::cout << VirtualMachine().eval( test.d(ETA).d(TIME_VARIABLE).d(TIME_VARIABLE), Point(0,1,0,-1) ) << std::endl ;

	test = s0/t0 ;
	
	std::cout << VirtualMachine().eval( test, Point(0,1,0,-1) ) << std::endl ;
	std::cout << VirtualMachine().eval( test.d(ETA), Point(0,1,0,-1) ) << std::endl ;
	std::cout << VirtualMachine().eval( test.d(ETA).d(TIME_VARIABLE), Point(0,1,0,-1) ) << std::endl ;
	std::cout << VirtualMachine().eval( test.d(ETA).d(TIME_VARIABLE).d(TIME_VARIABLE), Point(0,1,0,-1) ) << std::endl ;*/
	
	
	return 0 ;
	
  
  
  }

