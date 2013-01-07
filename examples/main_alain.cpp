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
#include "../physics/generalized_spacetime_viscoelasticity.h"
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

double sampleLength = 0.08 ; //5.5 ;
double sampleHeight = 0.04 ;
double supportLever = 0.05 ;//2.5 ;
double platewidth = 0.005 ;
double plateHeight = 0.0025 ;
Sample box(nullptr, sampleLength, sampleHeight,0.0,0.0) ;

int main(int argc, char *argv[])
{
	FeatureTree F(&box) ;
	F.setSamplingNumber(200) ;
//	F.setMaxIterationsPerStep(50000) ;
	F.setOrder(LINEAR) ;
	
	box.setBehaviour( new PasteBehaviour() ) ;

	ParticleSizeDistribution::get2DMortar(&F, new AggregateBehaviour(), 0.0025, 500) ;

	BoundingBoxAndRestrictionDefinedBoundaryCondition * load = new BoundingBoxAndRestrictionDefinedBoundaryCondition( SET_ALONG_ETA, TOP, -sampleLength*0.5, platewidth-sampleLength*0.5, -10, 10, 0. ) ;
	F.addBoundaryCondition(load) ;
	F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, LEFT ) ) ;
	F.addBoundaryCondition( new BoundingBoxNearestNodeDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM, Point(supportLever-sampleLength*0.5,-sampleHeight*0.5) ) ) ;


	F.step() ;
	Vector x = F.getAverageField(STRAIN_FIELD) ;
	Vector y = F.getAverageField(REAL_STRESS_FIELD) ;
	std::cout << 0. << "\t" << x[0] << "\t" << y[0] << std::endl ;

	size_t i = 0 ;

	std::fstream out ;
	out.open("tripoint_curve", std::ios::out) ;
	
	while(i < 100)
	{
		i++ ;
		load->setData(- 0.000001*i ) ;

		F.step() ;
		
		std::string tati = "tripoint" ;
		tati.append("_") ;
		tati.append(itoa(i)) ;
		std::cout << tati << std::endl ;
		TriangleWriter writer(tati, &F) ;
		writer.getField(TWFT_STRAIN) ;
		writer.getField(TWFT_STRESS) ;
		writer.getField(TWFT_DAMAGE) ;
		writer.getField(TWFT_STIFFNESS) ;
		writer.write() ;
			
		x = F.getAverageField(STRAIN_FIELD) ;
		y = F.getAverageField(REAL_STRESS_FIELD) ;
		out << 0.000001*i << "\t" << x[1] << "\t" << y[1]*sampleLength*sampleHeight*sampleHeight << "\t" << F.averageDamage << std::endl ;

	}
		
	return 0 ;
}

