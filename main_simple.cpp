// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "samplingcriterion.h"
#include "features/features.h"
#include "features/pore.h"
#include "features/sample.h"
#include "features/sample3d.h"
#include "features/inclusion.h"
#include "features/inclusion3d.h"
#include "physics/mohrcoulomb.h"

#include <fstream>

#include <cmath>
#include <typeinfo>
#include <limits>
#include <GL/glut.h>
#include <time.h> 
#define DEBUG 


using namespace Mu ;

FeatureTree * featureTree ;
std::vector<DelaunayTetrahedron *> triangles ;
DelaunayTree *dt ; //(pts) ;

double timepos = 0.1 ;
bool firstRun = true ;

Vector b(0) ;
Vector x(0) ;
Vector sigma(0) ; 
Vector sigma11(0) ; 
Vector sigma22(0) ; 
Vector sigma12(0) ; 
Vector epsilon(0) ; 
Vector epsilon11(0) ; 
Vector epsilon22(0) ; 
Vector epsilon12(0) ; 
Vector vonMises(0) ; 

double nu = 0.3 ;
double E = 20 ;

bool nothingToAdd = false ;
bool dlist = false ;
int count = 0 ;


void setBC()
{
	triangles = featureTree->getTetrahedrons() ;
	
	for(size_t k = 0 ; k < triangles.size() ;k++)
	{
		for(size_t c = 0 ;  c < triangles[k]->getBoundingPoints().size() ; c++ )
		{
			
			if (triangles[k]->getBoundingPoint(c).x < -3.999)
			{
				featureTree->getAssembly()->setPoint(-timepos,0 ,0,triangles[k]->getBoundingPoint(c).id) ;
			}
			if(triangles[k]->getBoundingPoint(c).x > 3.999)
			{
				featureTree->getAssembly()->setPoint( timepos,0, 0,triangles[k]->getBoundingPoint(c).id) ;
			}
			if (triangles[k]->getBoundingPoint(c).y < -2.999)
			{
				featureTree->getAssembly()->setPointAlong( ETA, 0, triangles[k]->getBoundingPoint(c).id) ;
			}
			if(triangles[k]->getBoundingPoint(c).y > 2.999)
			{
				featureTree->getAssembly()->setPointAlong( ETA, 0, triangles[k]->getBoundingPoint(c).id) ;
			}
			if (triangles[k]->getBoundingPoint(c).z < -2.999)
			{
				featureTree->getAssembly()->setPointAlong( ZETA, 0, triangles[k]->getBoundingPoint(c).id) ;
			}
			if(triangles[k]->getBoundingPoint(c).z > 2.999)
			{
				featureTree->getAssembly()->setPointAlong( ZETA, 0, triangles[k]->getBoundingPoint(c).id) ;
			}
		}
	}
	
}


void step()
{
	
	for(size_t i = 0 ; i < 10 ; i++)
	{
		std::cout << "\r iteration " << i << "/800" << std::flush ;
		setBC() ;
		featureTree->step(timepos) ;
		
// 		timepos+= 0.01 ;
	}
	std::cout << "\r iteration " << "800/800 ... done" << std::endl ;
	x.resize(featureTree->getDelaunayTree()->numPoints()*2) ;
	x = featureTree->getDisplacements() ;
	dt = featureTree->getDelaunayTree() ;
	sigma.resize(triangles.size()*3*3) ;

	epsilon.resize(triangles.size()*3*3) ;

// 	sigma = F.strainFromDisplacements() ;
// 	epsilon = F.stressFromDisplacements() ;
	std::pair<Vector, Vector > sigma_epsilon = featureTree->getStressAndStrain() ;
	sigma.resize(sigma_epsilon.first.size()) ;
	sigma = sigma_epsilon.first ;
	epsilon.resize(sigma_epsilon.second.size()) ;
	epsilon = sigma_epsilon.second ;

	sigma11.resize(sigma.size()/4) ;
	sigma22.resize(sigma.size()/4) ;
	sigma12.resize(sigma.size()/4) ;
	epsilon11.resize(sigma.size()/4) ;
	epsilon22.resize(sigma.size()/4) ;
	epsilon12.resize(sigma.size()/4) ;
	vonMises.resize(sigma.size()/4) ;

	std::cout << "unknowns :" << x.size() << std::endl ;
	std::cout << "max value :" << x.max() << std::endl ;
	std::cout << "min value :" << x.min() << std::endl ;


	for(size_t k = 0 ; k < triangles.size() ; k++)
	{

		if(!triangles[k]->getBehaviour()->fractured())
		{
			sigma11[k*3] = sigma[k*3*3];
			sigma22[k*3] = sigma[k*3*3+1];
			sigma12[k*3] = sigma[k*3*3+2];
			sigma11[k*3+1] = sigma[k*3*3+3];
			sigma22[k*3+1] = sigma[k*3*3+4];
			sigma12[k*3+1] = sigma[k*3*3+5];
			sigma11[k*3+2] = sigma[k*3*3+6];
			sigma22[k*3+2] = sigma[k*3*3+7];
			sigma12[k*3+2] = sigma[k*3*3+8];

			epsilon11[k*3] = epsilon[k*3*3];
			epsilon22[k*3] = epsilon[k*3*3+1];
			epsilon12[k*3] = epsilon[k*3*3+2];
			epsilon11[k*3+1] = epsilon[k*3*3+3];
			epsilon22[k*3+1] = epsilon[k*3*3+4];
			epsilon12[k*3+1] = epsilon[k*3*3+5];
			epsilon11[k*3+2] = epsilon[k*3*3+6];
			epsilon22[k*3+2] = epsilon[k*3*3+7];
			epsilon12[k*3+2] = epsilon[k*3*3+8];

			Vector vm0 = triangles[k]->getState().getPrincipalStresses(*triangles[k]->first) ;
			Vector vm1 = triangles[k]->getState().getPrincipalStresses(*triangles[k]->second) ;
			Vector vm2 = triangles[k]->getState().getPrincipalStresses(*triangles[k]->third) ;
			vonMises[k*3]   = sqrt(((vm0[0]-vm0[1])*(vm0[0]-vm0[1]))/2.) ;
			vonMises[k*3+1] = sqrt(((vm1[0]-vm1[1])*(vm1[0]-vm1[1]))/2.) ;
			vonMises[k*3+2] = sqrt(((vm2[0]-vm2[1])*(vm2[0]-vm2[1]))/2.) ;

// 			vonMises[k*3]   = sigma11[k*3]*epsilon11[k*3]     + sigma22[k*3]*epsilon22[k*3]     + sigma12[k*3]*epsilon12[k*3];
// 			vonMises[k*3+1] = sigma11[k*3+1]*epsilon11[k*3+1] + sigma22[k*3+1]*epsilon22[k*3+1] + sigma12[k*3+1]*epsilon12[k*3+2];
// 			vonMises[k*3+2] = sigma11[k*3+2]*epsilon11[k*3+2] + sigma22[k*3+2]*epsilon22[k*3+2] + sigma12[k*3+2]*epsilon12[k*3+2];
		}
		else
		{
			sigma11[k*3] = 0;
			sigma22[k*3] = 0;
			sigma12[k*3] = 0;
			sigma11[k*3+1] = 0;
			sigma22[k*3+1] = 0;
			sigma12[k*3+1] = 0;
			sigma11[k*3+2] = 0;
			sigma22[k*3+2] = 0;
			sigma12[k*3+2] = 0;

			epsilon11[k*3] = 0;
			epsilon22[k*3] = 0;
			epsilon12[k*3] = 0;
			epsilon11[k*3+1] = 0;
			epsilon22[k*3+1] = 0;
			epsilon12[k*3+1] = 0;
			epsilon11[k*3+2] = 0;
			epsilon22[k*3+2] = 0;
			epsilon12[k*3+2] = 0;

			vonMises[k*3]   = 0 ;
			vonMises[k*3+1] = 0 ;
			vonMises[k*3+2] = 0 ;
		}
	}
	std::cout << "max sigma11 :" << sigma11.max() << std::endl ;
	std::cout << "min sigma11 :" << sigma11.min() << std::endl ;
	std::cout << "max sigma12 :" << sigma12.max() << std::endl ;
	std::cout << "min sigma12 :" << sigma12.min() << std::endl ;
	std::cout << "max sigma22 :" << sigma22.max() << std::endl ;
	std::cout << "min sigma22 :" << sigma22.min() << std::endl ;

	std::cout << "max epsilon11 :" << epsilon11.max() << std::endl ;
	std::cout << "min epsilon11 :" << epsilon11.min() << std::endl ;
	std::cout << "max epsilon12 :" << epsilon12.max() << std::endl ;
	std::cout << "min epsilon12 :" << epsilon12.min() << std::endl ;
	std::cout << "max epsilon22 :" << epsilon22.max() << std::endl ;
	std::cout << "min epsilon22 :" << epsilon22.min() << std::endl ;

	std::cout << "max von Mises :" << vonMises.max() << std::endl ;
	std::cout << "min von Mises :" << vonMises.min() << std::endl ;
	
}


int main(int argc, char *argv[])
{

	Sample3D sample(NULL, 8,6,6,0,0,0) ;
	FeatureTree F(&sample) ;
	featureTree = &F ;

	Matrix m0(3,3) ;
	m0[0][0] = E/(1-nu*nu) ; m0[0][1] =E/(1-nu*nu)*nu ; m0[0][2] = 0 ;
	m0[1][0] = E/(1-nu*nu)*nu ; m0[1][1] = E/(1-nu*nu) ; m0[1][2] = 0 ; 
	m0[2][0] = 0 ; m0[2][1] = 0 ; m0[2][2] = E/(1-nu*nu)*(1.-nu)/2. ; 


	Inclusion3D * inc = new Inclusion3D(1, -1.5,0,0) ;
	F.addFeature(&sample,inc) ;
	

	
	inc->setBehaviour(new Stiffness(m0)) ;
	sample.setBehaviour(new StiffnessAndFracture(m0*0.5, new MohrCoulomb(30, -60))) ;

	F.sample(128) ;
	F.setOrder(QUADRATIC) ;
	F.generateElements() ;

	step() ;
	
	return 0 ;
}
