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
#include <time.h>
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
	ANALYTICAL,
	FDI,
	FDE,
	STFEM,
} TimeSteppingScheme ;


FeatureTree * featureTree ;
std::vector<DelaunayTriangle *> triangles ;

Sample box(NULL, 0.1,0.1,0.,0.) ;
Sample help(NULL, 0.1,0.1,0.,0.) ;

double young = 3.*1e9 ;
double nu = 0.3 ;
double eta = 60 ;
double stress = 10.*1e6 ;

Matrix C(3,3) ;
Matrix E(3,3) ;

double day = 86400. ;

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
	if(arg == std::string("stfem"))
		return STFEM ;
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
	double steps = 365./tau*day ;
	Vector instants((int) steps) ;
	instants[0] = tau/2 ;
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

Vector getFDIResults()
{
	double alphatau = tau * alpha ;
	Matrix CE = C + (E/alphatau) ;
	
	FeatureTree F(&box) ;
	F.setSamplingNumber(sampling) ;
	F.setOrder(LINEAR) ;
	
	FeatureTree H(&help) ;
	H.setSamplingNumber(sampling) ;
	H.setOrder(LINEAR) ;
	
	box.setBehaviour(new Stiffness(CE)) ;
	help.setBehaviour(new Stiffness(C)) ;

	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP, stress)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;

	Vector du ;
	Vector u ;
	Vector v ;
	Vector results(getInstants().size()) ;
	
	srand(0) ;
	F.step() ;
	du = F.getDisplacements() ;

	u.resize(du.size(),0.) ;
	v.resize(du.size(),0.) ;

	u = du ;
	v = du/alphatau ;
	
	results[0] = du.max() ;

	srand(0) ;
	H.step() ;

	for(size_t i = 1 ; i < results.size() ; i++)
	{
		F.resetBoundaryConditions() ;
		F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
		F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;

		Vector k(v.size()) ;
		k = H.getAssembly()->getMatrix()*v ;
		k *= tau ;
		for(size_t j = 0 ; j < k.size()/2 ; j++)
		{
			F.getAssembly()->addForceOn(XI,-k[j*2],j);
			F.getAssembly()->addForceOn(ETA,-k[j*2+1],j);
		}

		F.step() ;
		du = F.getDisplacements() ;

		u = u + v*tau + du ;
		v = v + du/alphatau ;
		std::cout << v.min() << std::endl ;

		results[i] = u.max() ;
	}
	return results ;
}

Vector getFDEResults()
{
	Vector u ;
	Vector ue ;
	Vector uv ;
	Vector v ;
	Vector k ;
	Vector results(getInstants().size()) ;

	FeatureTree F(&box) ;
	F.setSamplingNumber(sampling) ;
	F.setOrder(LINEAR) ;

	FeatureTree H(&help) ;
	H.setSamplingNumber(sampling) ;
	H.setOrder(LINEAR) ;

	box.setBehaviour(new Stiffness(C)) ;
	help.setBehaviour(new Stiffness(E)) ;

	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP, stress)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
	srand(0) ;
	F.step() ;
	ue.resize(F.getDisplacements().size()) ;
	ue = F.getDisplacements() ;

	srand(0) ;
	H.step() ;

	// get viscoelastic forces
	Vector f(ue.size()) ;
	f = H.getAssembly()->getMatrix() * ue ;
	f *= (1./tau) ;
	F.resetBoundaryConditions() ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
	for(size_t j = 0 ; j < f.size()/2 ; j++)
	{
		F.getAssembly()->addForceOn(XI, -f[j*2],  j);
		F.getAssembly()->addForceOn(ETA,-f[j*2+1],j);
	}

	// solve auxiliary problem
	F.step() ;
	uv.resize(F.getDisplacements().size()) ;
	uv = F.getDisplacements() ;
	k.resize(uv.size()) ;
	for(size_t i = 0 ; i < k.size() ; i++)
		k[i] = -ue[i] / uv[i] ;

	u.resize(uv.size()) ;
	v.resize(uv.size()) ;

	for(size_t i = 0 ; i < u.size() ; i++)
	{
		u[i] = (1.-std::exp(-k[i]))*uv[i] + ue[i] ;
		v[i] = ((k[i]*std::exp(-k[i]))*uv[i] + ue[i] ) / tau ;
	}
	
	results[0] = u.max() ;
	
	for(size_t i = 1 ; i < results.size() ; i++)
	{
		for(size_t j = 0 ; j < u.size() ; j++)
		{
			uv[j] = tau*v[j]/k[j] ;
			u[j] = u[j] + (1.-std::exp(-k[j]))*uv[j] ;
			v[j] = ((k[j]*std::exp(-k[j]))*uv[j]) / tau ;
		}
		results[i] = u.max() ;
	}

	return results ;

}


Vector getSTFEMResults()
{
	Vector u ;
	Vector results(getInstants().size()) ;

	FeatureTree F(&box) ;
	F.setSamplingNumber(sampling) ;
	F.setOrder(LINEAR_TIME_LINEAR) ;
	F.setDeltaTime(tau) ;

	box.setBehaviour(new KelvinVoight(C,E)) ;

	srand(0) ;
	F.step() ;
	std::vector<DelaunayTriangle *> tri = F.getElements2D() ;
	std::set<std::pair<std::pair<Point *, Point *>, DelaunayTriangle *> > pointList ;
	for(size_t i = 0 ; i < tri.size() ; i++)
	{	  
		pointList.insert(std::make_pair(std::make_pair(&tri[i]->getBoundingPoint(0),&tri[i]->getBoundingPoint(3)), tri[i])) ;
		pointList.insert(std::make_pair(std::make_pair(&tri[i]->getBoundingPoint(1),&tri[i]->getBoundingPoint(4)), tri[i])) ;
		pointList.insert(std::make_pair(std::make_pair(&tri[i]->getBoundingPoint(2),&tri[i]->getBoundingPoint(5)), tri[i])) ;
	}
	
	std::set<std::pair<DofDefinedBoundaryCondition *, size_t> > pointBC ;
	for(auto i = pointList.begin() ; i != pointList.end() ; i++)
	{
	  pointBC.insert(std::make_pair(new DofDefinedBoundaryCondition(SET_ALONG_XI, i->second, i->first.first->id, 0),i->first.second->id*2)) ;
	  pointBC.insert(std::make_pair(new DofDefinedBoundaryCondition(SET_ALONG_ETA,i->second, i->first.first->id, 0),i->first.second->id*2+1)) ;
	}
	
	for(auto i = pointBC.begin() ; i != pointBC.end() ; i++)
	  F.addBoundaryCondition(i->first) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP_AFTER, stress)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;

	F.step() ;
	u.resize(F.getDisplacements().size()) ;
	u = F.getDisplacements() ;
	
	results[0] = u.max() ;
	
	for(size_t i = 1 ; i < results.size() ; i++)
	{
		for(auto j = pointBC.begin() ; j != pointBC.end() ; j++)
		    j->first->setData(u[j->second]) ;
		
		F.step() ;
		u = F.getDisplacements() ;

		results[i] = u.max() ;
	}

	return results ;
}

Vector getAnalyticalResults()
{
	FeatureTree F(&box) ;
	F.setSamplingNumber(sampling) ;
	F.setOrder(LINEAR) ;
	
	Vector u ;
	
	box.setBehaviour(new Stiffness(C)) ;
	
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP, stress)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
	
	srand(0) ;
	F.step() ;
	u.resize(F.getDisplacements().size()) ;
	u = F.getDisplacements() ;
	
	double umax = u.max() ;
	
	Vector instants = getInstants() ;
	Vector results(instants.size()) ;
	for(int i = 0 ; i < results.size() ; i++)
		results[i] = umax*(1.-exp(-instants[i]/eta)) ;
	
	return results ;
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
		if(scheme == FDI)
			alpha = getAlpha(argv[4]) ;
		else
			alpha = 0 ;
		
		C[0][0] = 1. ; C[1][0] = nu ; C[2][0] = 0. ;
		C[0][1] = nu ; C[1][1] = 1. ; C[2][1] = 0. ;
		C[0][2] = 0. ; C[1][2] = 0. ; C[2][2] = 1.-nu ;
		C *= young/(1.-nu*nu) ;
		
		E = C*eta*day ;
		
		Vector instants = getInstants() ;
		Vector analytical = getAnalyticalResults() ;
		Vector fem ;
		switch(scheme)
		{
			case FDI:
				fem = getFDIResults() ;
				break ;
			case FDE:
				fem = getFDEResults() ;
				break ;
			case STFEM:
				fem = getSTFEMResults() ;
				break ;
		}
		

		double error = 0. ;
		for(size_t i = 0 ; i < analytical.size() ; i++)
			error += (analytical[i]-fem[i])*(analytical[i]-fem[i]) ;

		std::ofstream out ;
		std::string filename = getFileName(scheme) ;
		out.open(filename.c_str(), std::ios::out) ;
		out << error << std::endl ;
		out << std::endl ;
		for(int i = 0 ; i < fem.size() ; i++)
			out << instants[i]/day << "\t" << fem[i] << std::endl ;
		out.close() ;
		
		return 0 ;
}

