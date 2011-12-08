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
	ANALYTICAL,
	FDI,
	FDE,
	FULL,
	QUAD,
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

Vector getFDIResults()
{
	timeval t1, t2 ;
	
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
	gettimeofday(&t1, NULL) ;
	du = F.getDisplacements() ;

	u.resize(du.size(),0.) ;
	v.resize(du.size(),0.) ;

	u = du ;
	v = du/alphatau ;
	
	results[0] = du.max() ;
	gettimeofday(&t2, NULL) ;
	double delta = t2.tv_sec*1000000 - t1.tv_sec*1000000 + t2.tv_usec - t1.tv_usec ;
	std::cerr << "Post-processing " << delta/1e6 << std::endl ;
	
	srand(0) ;
	H.step() ;

	for(size_t i = 1 ; i < results.size() ; i++)
	{
		F.resetBoundaryConditions() ;
		F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
		F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;

		Vector k(v.size()) ;
		gettimeofday(&t1, NULL) ;
		k = H.getAssembly()->getMatrix()*v ;
		k *= tau ;
		for(size_t j = 0 ; j < k.size()/2 ; j++)
		{
			F.getAssembly()->addForceOn(XI,-k[j*2],j);
			F.getAssembly()->addForceOn(ETA,-k[j*2+1],j);
		}
		gettimeofday(&t2, NULL) ;
		delta = t2.tv_sec*1000000 - t1.tv_sec*1000000 + t2.tv_usec - t1.tv_usec ;
		std::cerr << "Boundary conditions " << delta/1e6 << std::endl ;
		
		F.step() ;
		du = F.getDisplacements() ;

		gettimeofday(&t1, NULL) ;
		u = u + v*tau + du ;
		v = v + du/alphatau ;
		gettimeofday(&t2, NULL) ;
		delta = t2.tv_sec*1000000 - t1.tv_sec*1000000 + t2.tv_usec - t1.tv_usec ;
		std::cerr << "Post-processing " << delta/1e6 << std::endl ;
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
	
	timeval t1, t2 ;

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
	std::cout << "MAIN" << std::endl ;
	srand(0) ;
	F.step() ;
	ue.resize(F.getDisplacements().size()) ;
	ue = F.getDisplacements() ;

	std::cout << "HELP" << std::endl ;
	H.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
	H.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
	srand(0) ;
	H.step() ;

	// get viscoelastic forces
	Vector f(ue.size()) ;
	gettimeofday(&t1, NULL) ;
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
	gettimeofday(&t2, NULL) ;
	
	double delta = t2.tv_sec*1000000 - t1.tv_sec*1000000 + t2.tv_usec - t1.tv_usec ;
	std::cerr << "Additional BC time " << delta/1e6 << std::endl ;
	
	
	// solve auxiliary problem
	F.step() ;
	gettimeofday(&t1, NULL) ;
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
	gettimeofday(&t2, NULL) ;
	delta = t2.tv_sec*1000000 - t1.tv_sec*1000000 + t2.tv_usec - t1.tv_usec ;
	std::cerr << "Post-processing " << delta/1e6 << std::endl ;
	
	return results ;

}


Vector getSTFEMResults(bool quad = true)
{
	Vector u ;
	Vector results(getInstants().size()) ;

	FeatureTree F(&box) ;
	F.setSamplingNumber(sampling) ;
	F.setOrder(LINEAR_TIME_LINEAR) ;
	if(quad)
		F.setOrder(QUADRATIC_TIME_QUADRATIC) ;
	F.setDeltaTime(tau) ;

	box.setBehaviour(new KelvinVoight(C,E)) ;

	srand(0) ;
	F.step() ;
	std::vector<DelaunayTriangle *> tri = F.getElements2D() ;
	std::set<std::pair<std::pair<Point *, Point *>, DelaunayTriangle *> > pointList ;
	if(quad)
	{
		std::cout << "quadratiques" << std::endl ;
		for(size_t i = 0 ; i < tri.size() ; i++)
		{	  
			pointList.insert(std::make_pair(std::make_pair(&tri[i]->getBoundingPoint(0),&tri[i]->getBoundingPoint(12)), tri[i])) ;
			pointList.insert(std::make_pair(std::make_pair(&tri[i]->getBoundingPoint(1),&tri[i]->getBoundingPoint(13)), tri[i])) ;
			pointList.insert(std::make_pair(std::make_pair(&tri[i]->getBoundingPoint(2),&tri[i]->getBoundingPoint(14)), tri[i])) ;
			pointList.insert(std::make_pair(std::make_pair(&tri[i]->getBoundingPoint(3),&tri[i]->getBoundingPoint(15)), tri[i])) ;
			pointList.insert(std::make_pair(std::make_pair(&tri[i]->getBoundingPoint(4),&tri[i]->getBoundingPoint(16)), tri[i])) ;
			pointList.insert(std::make_pair(std::make_pair(&tri[i]->getBoundingPoint(5),&tri[i]->getBoundingPoint(17)), tri[i])) ;
		}
	}
	else
	{
		for(size_t i = 0 ; i < tri.size() ; i++)
		{	  
			pointList.insert(std::make_pair(std::make_pair(&tri[i]->getBoundingPoint(0),&tri[i]->getBoundingPoint(3)), tri[i])) ;
			pointList.insert(std::make_pair(std::make_pair(&tri[i]->getBoundingPoint(1),&tri[i]->getBoundingPoint(4)), tri[i])) ;
			pointList.insert(std::make_pair(std::make_pair(&tri[i]->getBoundingPoint(2),&tri[i]->getBoundingPoint(5)), tri[i])) ;
		}
	}
	
	std::set<std::pair<DofDefinedBoundaryCondition *, size_t> > pointBC ;
	for(auto i = pointList.begin() ; i != pointList.end() ; i++)
	{
	  pointBC.insert(std::make_pair(new DofDefinedBoundaryCondition(SET_ALONG_XI, i->second, i->first.first->id, 0),i->first.second->id*2)) ;
	  pointBC.insert(std::make_pair(new DofDefinedBoundaryCondition(SET_ALONG_ETA,i->second, i->first.first->id, 0),i->first.second->id*2+1)) ;
	}
	
	for(auto i = pointBC.begin() ; i != pointBC.end() ; i++)
	  F.addBoundaryCondition(i->first) ;
	
	Function stressRamp("t") ;
	stressRamp = stressRamp + tau/2 ;
	stressRamp = stressRamp/tau ;
	stressRamp = stressRamp * stress*0.5 ;
	std::cout << VirtualMachine().eval(stressRamp, Point(0,0,0,-tau/2)) << std::endl ;
	std::cout << VirtualMachine().eval(stressRamp, Point(0,0,0, 0)) << std::endl ;
	std::cout << VirtualMachine().eval(stressRamp, Point(0,0,0, tau/2)) << std::endl ;
	
	Function stressConstant("1") ;
	stressConstant = stressConstant * (stress) ;
	
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP_AFTER, stressRamp)) ;
	if(quad)
		F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP_NOW, stressRamp)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;

	F.step() ;
	u.resize(F.getDisplacements().size()) ;
	u = F.getDisplacements() ;
	
	results[0] = u.max() ;
	F.resetBoundaryConditions() ;
	F.getAssembly()->clear() ;
	for(auto i = pointBC.begin() ; i != pointBC.end() ; i++)
		F.addBoundaryCondition(i->first) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP_AFTER, stressConstant)) ;
	if(quad)
		F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP_NOW, stressConstant)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
	
	
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

Vector getFullSTFEMResults()
{
	Vector u ;
	Vector results(getInstants().size()) ;

	FeatureTree F(&box) ;
	F.setSamplingNumber(sampling) ;
	F.setOrder(LINEAR_TIME_LINEAR) ;
	F.setDeltaTime(tau) ;
	F.instants.resize(results.size()+1) ;
	F.instants[0] = -tau/2 ;
	for(size_t i = 1 ; i < F.instants.size() ; i++)
		F.instants[i] = getInstants()[i-1] ;
	
	box.setBehaviour(new KelvinVoight(C,C*0.)) ;

	Function stressRamp("t") ;
	stressRamp = stressRamp + tau/2 ;
	stressRamp = stressRamp/tau ;
	stressRamp = stressRamp * stress ;
	
	Function stressConstant("1") ;
	stressConstant = stressConstant * (stress) ;
	
	Function time("t") ;
	time = time - tau*0.51 ;
	
	Function plus = f_positivity(time) ;
	Function minus = f_negativity(time) ;
	
	Function stressFunction = plus * stressConstant + minus * stressRamp ;
	
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP, stressFunction)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BEFORE)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, BEFORE)) ;

	srand(0) ;
	F.step() ;
/*	F.getAssembly()->print() ;*/
	u.resize(F.getDisplacements().size()) ;
	u = F.getDisplacements() ;

	std::vector<DelaunayTriangle *> tri = F.getElements2D() ;
	size_t pointsPerTriangle = tri[0]->getBoundingPoints().size() ;
	for(size_t i = 0 ; i < tri.size() ; i++)
	{
		if(i%1000 == 0)
			std::cout << "processing triangle " << i << "/" << tri.size() << std::endl ;
	
		for(size_t j = 0 ; j < pointsPerTriangle ; j++)
		{
			Point p = tri[i]->getBoundingPoint(j) ;
			size_t t = 1 ;
			bool found = false ;
			while( t < F.instants.size() && !found)
			{
				if(p.t == F.instants[t])
					found = true ;
				else t++ ;
			}
			if(found)
			{
				t-- ;
				double umax = u[p.id*2+1] ;
				if(umax > results[t])
				{
					results[t] = umax ;
				}
			}
		}
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
	{
		results[i] = (1.-std::exp(-instants[i]/eta))*umax ;
	}
	
	return results ;
}

int main(int argc, char *argv[])
{
/*		Function t("t t *") ;
		TriElement * element = new TriElement(LINEAR_TIME_LINEAR) ;
		std::valarray<Point *> points(6) ;
		points[0] = new Point(0,1,0,-1) ;
		points[1] = new Point(0,0,0,-1) ;
		points[2] = new Point(1,0,0,-1) ;
		points[3] = new Point(0,1,0,1) ;
		points[4] = new Point(0,0,0,1) ;
		points[5] = new Point(1,0,0,1) ;
		element->setBoundingPoints(points) ;
		std::cout << element->getBoundingPoints().size() << std::endl ;
		
		std::cout << element->jacobianAtPoint(Point(0,0,0,-0.5)) << std::endl ;
		std::cout << element->jacobianAtPoint(Point(0,0,0,+0.5)) << std::endl ;
		
		std::cout << VirtualMachine().ieval(t, element) << std::endl ;
		
		return 0 ;*/
	
	
		std::cout << "usage = ./visco <scheme> <time-step> <sampling-number> <alpha>" << std::endl ;
		std::cout << "\t<scheme>\tstring representing the time-stepping scheme among <fdi> <fde> <stfem>" << std::endl ;
		std::cout << "\t<time-step>\tdouble representing the number of days between two time steps" << std::endl ;
		std::cout << "\t<sampling-number>\tinteger representing the number points at the boundary of the sample" << std::endl ;
		std::cout << "\t<alpha>\tdouble between 0 and 1 between (for <fdi> only)" << std::endl ;
		
		scheme = getTimeSteppingScheme(argv[1]) ;
		tau = getTimeStep(argv[2]) ;
		sampling = getSampling(argv[3]) ;
//		if(scheme == FDI)
//			alpha = getAlpha(argv[4]) ;
//		else
			alpha = 0.5 ;
		
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
			case FULL:
				fem = getFullSTFEMResults() ;
				break ;
			case QUAD:
				fem = getSTFEMResults(true) ;
				break ;
			case STFEM:
				fem = getSTFEMResults(false) ;
				break ;
		}
		

		double error = 0. ;
		for(size_t i = 0 ; i < analytical.size() ; i++)
			error += (analytical[i]-fem[i])*(analytical[i]-fem[i]) ;

		std::ofstream out ;
		std::string filename = getFileName(scheme) ;
		out.open(filename.c_str(), std::ios::out) ;
		for(int i = 0 ; i < analytical.size() ; i++)
			out << std::setprecision(16) << instants[i]/day << "\t" << std::setprecision(16) << fem[i] << std::endl ;
		out.close() ;
		
		return 0 ;
}

