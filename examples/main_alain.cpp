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
#include "../utilities/random.h"
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

Sample box(nullptr, 0.08, 0.08,0.,0.) ;


int main(int argc, char *argv[])
{
	FeatureTree F(&box) ;
	box.setBehaviour(new ViscoDamagePasteBehaviour()) ;
	F.setDeltaTime(0.1) ;
	F.setSamplingNumber(5) ;
	F.setOrder(LINEAR_TIME_LINEAR) ;
	F.setMaxIterationsPerStep(50) ;

	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0,1)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0,3)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0,5)) ;

	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0,0)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0,2)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0,4)) ;
		
	F.step() ;

	BoundingBoxDefinedBoundaryCondition * disp = new BoundingBoxDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, TOP_AFTER, 0.,1) ;
	F.addBoundaryCondition(disp) ;

	size_t j = 0 ;
	for(size_t i = 0 ; i < 10 ; i++)
	{
		disp->setData(1e-5*j) ;
		F.setDeltaTime(0.1) ;
		bool goOn = F.step() ;

		std::string ex("ex_") ;
		ex.append(itoa(j)) ;
		if(!goOn)
		{
			ex.append("_inter_") ;			
			ex.append(itoa(i)) ;
		}
		TriangleWriter w(ex, &F, -1+2*goOn) ;
		w.getField(STRAIN_FIELD) ;
		w.getField(TWFT_DAMAGE) ;
		w.write() ;
		if(goOn)
			j++ ;
	}

	return 0 ;
}
