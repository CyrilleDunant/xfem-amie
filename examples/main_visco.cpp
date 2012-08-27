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
#include "../physics/generalized_spacetime_viscoelasticity.h"
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
	KELVINVOIGHT,
	MX,
	STANDARDLINEARSOLID,
} ViscoModel ;

FeatureTree * featureTree ;
std::vector<DelaunayTriangle *> triangles ;

Sample box(nullptr, 1.,1.,0.5,0.5) ;

double young = 1.*1e9 ;
double nu = 0.3 ;
double eta = 1. ;
double stress = 1.*1e6 ;

double tau = 0.1 ;
int sampling = 100 ;

Matrix C(3,3) ;
Matrix E(3,3) ;
Matrix K(3,3) ;

ViscoModel model ;

ViscoModel getViscoModel(std::string arg)
{
	if(arg == std::string("kv"))
		return KELVINVOIGHT ;
	if(arg == std::string("mx"))
		return MX ;
	if(arg == std::string("sls"))
		return STANDARDLINEARSOLID ;
}

double getTimeStep(std::string arg)
{
	return atof(arg.c_str()) ;
}

int getSampling(std::string arg)
{
	return atoi(arg.c_str()) ;
}

Vector getInstants()
{
	double steps = 20./tau ;
	Vector instants((int) steps) ;
	instants[0] = tau ;
	for(int i = 1 ; i < instants.size() ; i++)
		instants[i] = instants[i-1] + tau ;
	return instants ;
}

Form * getBehaviour( ViscoModel v , Matrix & C, Matrix & E, Matrix & K)
{
	switch(v)
	{
		case KELVINVOIGHT:
			return new KelvinVoight(C,E) ;
		case MX:
			return new Maxwell(C,E) ;
		case STANDARDLINEARSOLID:
			return new StandardLinearSolid(K, C, E) ;
	}
	return new KelvinVoight(C,E) ;
}

std::string getFileName(ViscoModel v)
{
	std::string name = "visco_" ;
	switch(v)
	{
		case KELVINVOIGHT:
			name.append("kv_") ;
			break ;
		case MX:
			name.append("mx_") ;
			break ;
		case STANDARDLINEARSOLID:
			name.append("sls_") ;
			break ;
	}
	
	name.append(itoa(getInstants().size())) ;
	name.append("_") ;
	name.append(itoa(sampling)) ;
	name.append(".txt") ;
	return name ;
}

Vector getAnalytical(ViscoModel v)
{
	Vector instants = getInstants() ;
	Vector result( 0., instants.size()) ;
	for(size_t i = 0 ; i < result.size() ; i++)
	{
		double t = instants[i] ;
		switch(v)
		{
			case KELVINVOIGHT:
				result[i] = 0.001*(1-exp(-t/eta)) ;
				break ;
			case MX:
				result[i] = 0.001*(1.+t/eta) ;
				break ;
			case STANDARDLINEARSOLID:
				result[i] = 0.001*(1./3.-(1./3.-1./4.)*exp(-3*t/(4*eta))) ;
				break ;
		}
	  
	}
	return result ;
}

std::string getTWFileName(ViscoModel v, int i)
{
	std::string name = "trg/" ;
	switch(v)
	{
		case KELVINVOIGHT:
			name.append("kv_") ;
			break ;
		case MX:
			name.append("mx_") ;
			break ;
		case STANDARDLINEARSOLID:
			name.append("sls_") ;
			break ;
	}
	name.append(itoa(i)) ;
	return name ;
}

int main(int argc, char *argv[])
{
	std::cout << "usage = ./visco <model> <time-step> <sampling-number>" << std::endl ;
	std::cout << "\t<model>\tstring representing the viscoelastic model among <kv> <mx> <sls>" << std::endl ;
	std::cout << "\t<time-step>\tdouble representing the number of days between two time steps" << std::endl ;
	std::cout << "\t<sampling-number>\tinteger representing the number points at the boundary of the sample" << std::endl ;
	
	model = getViscoModel( argv[1]) ;
	tau = getTimeStep(argv[2]) ;
	sampling = getSampling(argv[3]) ;
	
	C[0][0] = 1. ; C[1][0] = nu ; C[2][0] = 0. ;
	C[0][1] = nu ; C[1][1] = 1. ; C[2][1] = 0. ;
	C[0][2] = 0. ; C[1][2] = 0. ; C[2][2] = 1.-nu ;
	
	C *= young/(1.-nu*nu) ;
	E = C*eta ;
	K = C ;
	
	FeatureTree F(&box) ;
	F.setSamplingNumber(sampling) ;
	F.setOrder(LINEAR_TIME_LINEAR) ;
//	F.setOrder(LINEAR) ;
	F.setDeltaTime(tau) ;

	box.setBehaviour( getBehaviour( model, C, E, K ) ) ;
	box.setBehaviour( new GeneralizedSpaceTimeViscoelasticity( PURE_ELASTICITY, C, 1))  ;
//	box.setBehaviour( new Stiffness( K ) ) ;
	
	Inclusion * hole = new Inclusion( 0.2, 0., 0.) ;
	hole->setBehaviour(new VoidForm()) ;
//	F.addFeature(&box, hole) ;
	
//	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM, 0, 0)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM, 0, 1)) ;
//	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM, 0, 2)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM, 0, 3)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT, 0, 0)) ;
//	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT, 0, 1)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT, 0, 2)) ;
//	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT, 0, 3)) ;
//	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP, stress)) ;
	
	srand(0) ;
	std::vector<DelaunayTriangle *> tri ;
	std::set<std::pair<std::pair<Point *, Point *>, DelaunayTriangle *> > pointList ;
	std::set<std::pair<DofDefinedBoundaryCondition *, size_t> > pointBC ;
	BoundingBoxDefinedBoundaryCondition * stressBCRamp = nullptr ;
	BoundingBoxDefinedBoundaryCondition * stressBCConstant = nullptr ;
	Function stressRamp("t") ;
	stressRamp = stressRamp + tau/2 ;
	stressRamp = stressRamp/tau ;
	stressRamp = stressRamp * stress ;
	Function stressConstant("1") ;
	stressConstant = stressConstant * stress ;
	Function stressNull("1") ;
	stressNull = stressNull * 0. ;

	F.step() ;

	stressBCRamp = new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP_AFTER, stressRamp) ;
	stressBCConstant = new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP_AFTER, 0.) ;
	
	
	tri = F.getElements2D() ;
	for(size_t i = 0 ; i < tri.size() ; i++)
	{	  
		pointList.insert(std::make_pair(std::make_pair(&tri[i]->getBoundingPoint(0),&tri[i]->getBoundingPoint(3)), tri[i])) ;
		pointList.insert(std::make_pair(std::make_pair(&tri[i]->getBoundingPoint(1),&tri[i]->getBoundingPoint(4)), tri[i])) ;
		pointList.insert(std::make_pair(std::make_pair(&tri[i]->getBoundingPoint(2),&tri[i]->getBoundingPoint(5)), tri[i])) ;
	}
	for(auto i = pointList.begin() ; i != pointList.end() ; i++)
	{
		GaussPointArray gp = i->second->getGaussPoints() ;
		std::valarray<Matrix> Jinv( Matrix(), i->second->getGaussPoints().gaussPoints.size() ) ;
						
		for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
		{
			i->second->getInverseJacobianMatrix( gp.gaussPoints[j].first, Jinv[j] ) ;
		}
		
		pointBC.insert(std::make_pair(new DofDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS,i->second,gp,Jinv, i->first.first->id, 0, 0),i->first.second->id*4)) ;
		pointBC.insert(std::make_pair(new DofDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS,i->second,gp,Jinv, i->first.first->id, 0, 1),i->first.second->id*4+1)) ;
		pointBC.insert(std::make_pair(new DofDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS,i->second,gp,Jinv, i->first.first->id, 0, 2),i->first.second->id*4+2)) ;
		pointBC.insert(std::make_pair(new DofDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS,i->second,gp,Jinv, i->first.first->id, 0, 3),i->first.second->id*4+3)) ;
	}

	for(auto i = pointBC.begin() ; i != pointBC.end() ; i++)
		F.addBoundaryCondition(i->first) ;
		
	F.addBoundaryCondition(stressBCRamp) ;
	F.addBoundaryCondition(stressBCConstant) ;


	size_t imax = getInstants().size() ;
	Vector analytical = getAnalytical(model) ;
	Vector fem( 0., analytical.size() ) ;
//	imax = 1 ;
	
	for(size_t i = 0 ; i < imax ; i++)
	{
		double t = getInstants()[i] ;
		F.step() ;
		Vector x = F.getDisplacements() ;
		double femi = x[0] ;
		for(size_t j = 0 ; j < x.size()/4 ; j++)
		{
			if(x[j*4+0] > femi)
				femi = x[j*4+0] ;
			if(x[j*4+1] > femi)
				femi = x[j*4+1] ;
		}
		fem[i] = femi ;
		
		
		TriangleWriter w( getTWFileName( model, i ), &F, 1) ;
		w.getField(TWFT_STRESS) ;
/*		w.getField(TWFT_STRESS) ;
		w.write() ;
		std::cout << "yes ?" << std::endl ;*/

		for(auto j = pointBC.begin() ; j != pointBC.end() ; j++)
			j->first->setData(x[j->second]) ;

		if(i == 0)
		{
			stressBCRamp->setData( stressNull ) ;
			stressBCConstant->setData( stress ) ;
		}		
	}
	
	
	std::string file = getFileName( model ) ;
	std::ofstream out ;
	out.open(file.c_str(), std::ios::out) ;
	
	for(size_t i = 0 ; i < analytical.size() ; i++)
		out << std::setprecision(16) << analytical[i] << "\t" << std::setprecision(16) << fem[i] << std::endl ;
	
	return 0 ;
}

