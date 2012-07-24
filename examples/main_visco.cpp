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
#include "../physics/stiffness.h"
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

typedef enum
{
	ANALYTICAL = 0,
	FDI = 1,
	FDE = 2,
	STFEM = 3,
	FULL = 4,
	QUAD,
} TimeSteppingScheme ;


FeatureTree * featureTree ;
std::vector<DelaunayTriangle *> triangles ;

Sample box(nullptr, 0.1,0.1,0.,0.) ;
Sample help(nullptr, 0.1,0.1,0.,0.) ;

double young = 3.*1e9 ;
double nu = 0.3 ;
double eta = 60. ;
double stress = 10.*1e6 ;

Matrix C(3,3) ;
Matrix E(3,3) ;

double day = 1 ;

double tau ;
double alpha ;
int sampling ;
TimeSteppingScheme scheme ;

TimeSteppingScheme getTimeSteppingScheme(std::string arg)
{
	if(arg == std::string("fdi"))
			return FDI ;
	if(arg == std::string("fde"))
		return FDE ;
	if(arg == std::string("full"))
		return FULL ;
	if(arg == std::string("stfem"))
		return STFEM ;
	if(arg == std::string("quad"))
		return QUAD ;
	return ANALYTICAL ;
}

double getTimeStep(std::string arg)
{
	return atof(arg.c_str())*day ;
}

int getSampling(std::string arg)
{
	return atoi(arg.c_str()) ;
}

double getAlpha(std::string arg)
{
	double a = atof(arg.c_str()) ;
	if(a < 0)
		a = -a ;
	while(a > 1)
		a *= 0.1 ;
	return a ;
}

Vector getInstants()
{
	double steps = 2.*365./tau*day ;
	Vector instants((int) steps) ;
	instants[0] = tau/2. ;
	for(int i = 1 ; i < instants.size() ; i++)
		instants[i] = instants[i-1] + tau ;
	return instants ;
}

std::string getFileName(TimeSteppingScheme s)
{
	std::string name = "visco_" ;
	switch(s)
	{
		case FDI:
			name.append("fdi_") ;
			break ;
		case FDE:
			name.append("fde_") ;
			break ;
		case FULL:
			name.append("full_") ;
			break ;
		case QUAD:
			name.append("quad_") ;
			break ;
		case STFEM:
			name.append("stfem_") ;
			break ;
	}
	
	name.append(itoa(tau)) ;
	name.append("_") ;
	name.append(itoa(sampling)) ;
	if(s == FDI)
	{
		name.append("_") ;
		name.append(itoa(alpha*100)) ;
	}
	name.append(".txt") ;
	return name ;
}

TimeSteppingScheme getTimeSteppingScheme(int i)
{
	switch(i)
	{
		case 0:
			return ANALYTICAL ;
		case 1:
			return FDI ;
		case 2:
			return FDE ;
		case 3:
			return STFEM ;
	}
	return ANALYTICAL ;
}


int main(int argc, char *argv[])
{
	std::cout << "usage = ./visco <scheme> <time-step> <sampling-number> <alpha>" << std::endl ;
	std::cout << "\t<scheme>\tstring representing the time-stepping scheme among <fdi> <fde> <stfem>" << std::endl ;
	std::cout << "\t<time-step>\tdouble representing the number of days between two time steps" << std::endl ;
	std::cout << "\t<sampling-number>\tinteger representing the number points at the boundary of the sample" << std::endl ;
	std::cout << "\t<alpha>\tdouble between 0 and 1 between (for <fdi> only)" << std::endl ;
	
	
	
	scheme = getTimeSteppingScheme(argv[1]) ;
	tau = getTimeStep(argv[2]) ;
	sampling = getSampling(argv[3]) ;
	alpha = getAlpha(argv[4]) ;
	
	C[0][0] = 1. ; C[1][0] = nu ; C[2][0] = 0. ;
	C[0][1] = nu ; C[1][1] = 1. ; C[2][1] = 0. ;
	C[0][2] = 0. ; C[1][2] = 0. ; C[2][2] = 1.-nu ;
	
	C *= young/(1.-nu*nu) ;
	E = C*eta*day ;
	Vector decay(3) ;
	decay[0] = decay[1] = decay[2] = (eta*day) ;

	Matrix results(getInstants().size(), 5) ;
	
/*	for(int s = 0 ; s < 4 ; s++)
	{
		scheme = getTimeSteppingScheme(s) ;*/
	
		FeatureTree F(&box) ;
		F.setSamplingNumber(sampling) ;
		if(scheme == STFEM)
			F.setOrder(LINEAR_TIME_LINEAR) ;
		else
			F.setOrder(LINEAR) ;
		F.setDeltaTime(tau) ;
	
		switch(scheme)
		{
			case FDI:
				box.setBehaviour( new NewmarkNumeroffKelvinVoigt( C, decay, alpha ) ) ;
				break ;
			case FDE:
				box.setBehaviour( new ExponentiallyPredictedKelvinVoigt( C, decay ) ) ;
				break ;
			case STFEM:
				box.setBehaviour( new KelvinVoight( C, E, eta*day ) ) ;
				break ;
			case ANALYTICAL:
				box.setBehaviour( new Stiffness( C ) ) ;
				break ;
		}
	
		F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
		F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
	
		srand(0) ;
		std::vector<DelaunayTriangle *> tri ;
		std::set<std::pair<std::pair<Point *, Point *>, DelaunayTriangle *> > pointList ;
		std::set<std::pair<DofDefinedBoundaryCondition *, size_t> > pointBC ;
		BoundingBoxDefinedBoundaryCondition * stressBC = nullptr ;
		Function stressRamp("t") ;
		switch(scheme)
		{
			case STFEM:
				F.step() ;
				tri = F.getElements2D() ;
				for(size_t i = 0 ; i < tri.size() ; i++)
				{	  
					pointList.insert(std::make_pair(std::make_pair(&tri[i]->getBoundingPoint(0),&tri[i]->getBoundingPoint(3)), tri[i])) ;
					pointList.insert(std::make_pair(std::make_pair(&tri[i]->getBoundingPoint(1),&tri[i]->getBoundingPoint(4)), tri[i])) ;
					pointList.insert(std::make_pair(std::make_pair(&tri[i]->getBoundingPoint(2),&tri[i]->getBoundingPoint(5)), tri[i])) ;
				}
				for(auto i = pointList.begin() ; i != pointList.end() ; i++)
				{
					pointBC.insert(std::make_pair(new DofDefinedBoundaryCondition(SET_ALONG_XI, i->second, i->first.first->id, 0),i->first.second->id*2)) ;
					pointBC.insert(std::make_pair(new DofDefinedBoundaryCondition(SET_ALONG_ETA,i->second, i->first.first->id, 0),i->first.second->id*2+1)) ;
				}
	
				for(auto i = pointBC.begin() ; i != pointBC.end() ; i++)
					F.addBoundaryCondition(i->first) ;
			
				stressRamp = stressRamp + tau/2 ;
				stressRamp = stressRamp/tau ;
				stressRamp = stressRamp * stress*0.5 ;
				stressBC = new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP_AFTER, stressRamp) ;
		
				F.addBoundaryCondition(stressBC) ;
				break ;
			case FDI:
			case FDE:
			case ANALYTICAL:
				stressBC = new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP, stress) ;
				break ;
			
		}
	
		srand(0) ;
		size_t imax = 5 ;
		imax = getInstants().size() ;
		if(scheme == FDE || scheme == FDI)
			imax += 2 ;
		
		for(size_t i = 0 ; i < imax ; i++)
		{
			if(i == 0 && scheme == ANALYTICAL)
			{
				F.addBoundaryCondition(stressBC) ;			
			}

			if(i == 2 && (scheme == FDI || scheme == FDE ))
			{
				F.addBoundaryCondition(stressBC) ;
			}
	
			F.step() ;
			Vector x = F.getDisplacements() ;
			std::pair<Vector, Vector> all = F.getStressAndStrain() ;
			std::cout << all.first.max() << std::endl ;
			if(scheme == ANALYTICAL)
				x *= 1. - std::exp( - getInstants()[i] / eta ) ;
/*			if(scheme == FDI || scheme == FDE)
			{
				if(i >= 2)
				{
					results[i-2][s+1] = x.max() ;
					results[i-2][0] = getInstants()[i] ;
				}
			}
			else
			{
					results[i][s+1] = x.max() ;
					results[i][0] = getInstants()[i] ;
			}*/

			for(auto j = pointBC.begin() ; j != pointBC.end() ; j++)
				j->first->setData(x[j->second]) ;

			if(i == 1 && scheme == STFEM)
			{
				Function stressConstant("1") ;
				stressConstant = stressConstant * (stress) ;
				stressBC->setData( stressConstant ) ;
			}
		
		
		}
//	}
	
//	results.print() ;
		
	return 0 ;
}

