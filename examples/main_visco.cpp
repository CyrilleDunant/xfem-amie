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

Sample box(nullptr, 1., 1.,0.,0.) ;
double decrease = 1. ;

double fr(Point p, double t)
{
	double r = std::sqrt(p.x*p.x+p.y*p.y) ;
	if(r < 0.24 || std::abs(p.x) >= 0.5 || std::abs(p.y) >= 0.5 || decrease <= 0)
		return 0. ;
	double rmax = 0.5 ;
	if(p.y > 0 && p.y > std::abs(p.x))
	{
		double xtmp = std::abs(p.x)*0.5/p.y ;
		rmax = std::sqrt(xtmp*xtmp+0.5*0.5) ;
	}
	else if(p.y < 0 && p.y < -std::abs(p.x))
	{
		double xtmp = std::abs(p.x)*0.5/p.y ;
		rmax = std::sqrt(xtmp*xtmp+0.5*0.5) ;
	}
	else if(p.x > 0 && p.x > std::abs(p.y))
	{
		double ytmp = std::abs(p.y)*0.5/p.x ;
		rmax = std::sqrt(ytmp*ytmp+0.5*0.5) ;
	}
	else if(p.x < 0 && p.x < -std::abs(p.y))
	{
		double ytmp = std::abs(p.y)*0.5/p.x ;
		rmax = std::sqrt(ytmp*ytmp+0.5*0.5) ;
	}
	double x = (rmax-r)/(rmax-0.25) ;
	if(decrease <= 0)
	{
		return 0. ;
	}
	return x*(1.-exp(-t/decrease))*(0.3-0.25)/0.25 ;
	
}

double fx(Point p, double t)
{
	if(std::abs(p.y) == 0.5)
		return p.x ;
	if(decrease == 0)
		return p.x ;
	if(decrease == -1)
		return p.x*0.75 ;
	return p.x*(1.+fr(p,t)) ;
}

double fy(Point p, double t)
{
	if(std::abs(p.y) == 0.5)
		return p.y ;
	if(decrease == 0)
		return p.y ;
	if(decrease == -1)
		return p.y*0.75 ;
	return p.y*(1.+fr(p,t)) ;
}

int main(int argc, char *argv[])
{
	double timestep = atof(argv[1]) ;
	decrease = atof(argv[2]) ;
	int sampling = (int) atof(argv[3]) ;
  
	FeatureTree F(&box) ;
	F.setSamplingNumber(sampling) ;
	
	Matrix e = (new ElasticOnlyPasteBehaviour(1e9, 0.3))->param ;
  	box.setBehaviour(new Viscoelasticity(BURGER, e*30, e*30*300, e*14, e*14*5000)) ;

	Inclusion * hole = new Inclusion(0.25, 0.,0.) ;
	if(decrease < 0)
		hole = new Inclusion(0.3,0.,0.) ;
	hole->setBehaviour( new VoidForm() ) ;
	F.addFeature( &box, hole) ;
	
	F.setOrder(LINEAR_TIME_LINEAR) ;
	F.setDeltaTime(timestep) ;

	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0,0)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0,1)) ;
 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0,2)) ;
 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0,3)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0,4)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0,5)) ;
	F.addBoundaryCondition(new TimeContinuityBoundaryCondition()) ;
	F.step() ;
	std::cout << "hip" << 0. << "\t" << 0.25+(1.-exp(-0./decrease))*(0.3-0.25) << "\t" << F.getAverageField(REAL_STRESS_FIELD,-1,1)[1] << "\t" << F.getAverageField(STRAIN_FIELD, -1,1)[1]  << std::endl ;

	BoundingBoxDefinedBoundaryCondition * disp = new BoundingBoxDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, TOP_AFTER, -0.0005, 1) ;
	F.addBoundaryCondition(disp) ;

	F.step() ;

	double time = timestep ;	
	std::cout << "hap" << time << "\t" << 0.25+(1.-exp(-time/decrease))*(0.3-0.25) << "\t" << F.getAverageField(REAL_STRESS_FIELD,-1,1)[1] << "\t" << F.getAverageField(STRAIN_FIELD, -1,1)[1]  << std::endl ;
	
	disp->setData(-0.001) ;
	
	F.step() ;
	time += timestep ;	
	std::cout << "hop" << time << "\t" << 0.25+(1.-exp(-time/decrease))*(0.3-0.25) << "\t" << F.getAverageField(REAL_STRESS_FIELD,-1,1)[1] << "\t" << F.getAverageField(STRAIN_FIELD, -1,1)[1]  << std::endl ;

	std::vector<Point> nodes = F.getNodes() ;
	
	while(time < 1500)
	{
		F.step() ;

		
		time += timestep ;
		std::cout << time << "\t" << 0.25*(1.+fr(Point(0.25,0.), time))  << "\t" << F.getAverageField(REAL_STRESS_FIELD,-1,0)[1] << "\t" << F.getAverageField(STRAIN_FIELD, -1,0)[1]  << std::endl ;

		
		std::vector<DelaunayTriangle *> tri = F.getElements2D() ;
		std::valarray<bool> done(F.getDisplacements().size()) ;
		done = (decrease <= 0) ;
		
		if(!done[0])
		{
			for(size_t i = 0 ; i < tri.size() ; i++)
			{
				for(size_t j = 0 ; j < 3 ; j++)
				{
					if(!done[tri[i]->getBoundingPoint(j).id])
					{
						done[tri[i]->getBoundingPoint(j).id] = true ;
						tri[i]->getBoundingPoint(j).x = tri[i]->getBoundingPoint(j+3).x ;
						tri[i]->getBoundingPoint(j).y = tri[i]->getBoundingPoint(j+3).y ;
						tri[i]->getBoundingPoint(j+3).x = fx(nodes[tri[i]->getBoundingPoint(j).id], time-timestep) ;
						tri[i]->getBoundingPoint(j+3).y = fy(nodes[tri[i]->getBoundingPoint(j).id], time-timestep);
					}
				}
				tri[i]->clearElementaryMatrix() ;
				tri[i]->moved = true ;
			}
		}
	}

	TriangleWriter writer("yop", &F, 1) ;
	writer.getField(TWFT_STIFFNESS) ;
	writer.write() ;

		
	return 0 ;
}

