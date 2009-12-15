// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../utilities/samplingcriterion.h"
#include "../features/features.h"
#include "../features/pore.h"
#include "../features/sample.h"
#include "../features/sample3d.h"
#include "../features/inclusion.h"
#include "../features/inclusion3d.h"
#include "../physics/mohrcoulomb.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"

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
DelaunayTree3D *dt ; //(pts) ;

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
			
			if (triangles[k]->getBoundingPoint(c).x < .00001)
			{
				featureTree->getAssembly()->setPoint(-timepos,0 ,0,triangles[k]->getBoundingPoint(c).id) ;
			}
			if(triangles[k]->getBoundingPoint(c).x > .03999)
			{
				featureTree->getAssembly()->setPoint( timepos,0, 0,triangles[k]->getBoundingPoint(c).id) ;
			}
			if (triangles[k]->getBoundingPoint(c).y < .00001)
			{
				featureTree->getAssembly()->setPointAlong( ETA, 0, triangles[k]->getBoundingPoint(c).id) ;
			}
			if(triangles[k]->getBoundingPoint(c).y > .03999)
			{
				featureTree->getAssembly()->setPointAlong( ETA, 0, triangles[k]->getBoundingPoint(c).id) ;
			}
			if (triangles[k]->getBoundingPoint(c).z < .00001)
			{
				featureTree->getAssembly()->setPointAlong( ZETA, 0, triangles[k]->getBoundingPoint(c).id) ;
			}
			if(triangles[k]->getBoundingPoint(c).z > .03999)
			{
				featureTree->getAssembly()->setPointAlong( ZETA, 0, triangles[k]->getBoundingPoint(c).id) ;
			}
		}
	}
	
}


void step()
{
	
	for(size_t i = 0 ; i < 1 ; i++)
	{
		std::cout << "\r iteration " << i << "/2" << std::flush ;
		setBC() ;
		featureTree->step(timepos) ;
		
// 		timepos+= 0.01 ;
	}
	std::cout << "\r iteration " << "2/2 ... done" << std::endl ;
	x.resize(featureTree->getDelaunayTree3D()->numPoints()*3) ;
	x = featureTree->getDisplacements() ;
	dt = featureTree->getDelaunayTree3D() ;
	sigma.resize(triangles.size()*6*4) ;

	epsilon.resize(triangles.size()*6*4) ;

// 	sigma = F.strainFromDisplacements() ;
// 	epsilon = F.stressFromDisplacements() ;
	std::pair<Vector, Vector > sigma_epsilon = featureTree->getStressAndStrain() ;
	sigma.resize(sigma_epsilon.first.size()) ;
	sigma = sigma_epsilon.first ;
	epsilon.resize(sigma_epsilon.second.size()) ;
	epsilon = sigma_epsilon.second ;

	sigma11.resize(sigma.size()/6) ;
	sigma22.resize(sigma.size()/6) ;
	sigma12.resize(sigma.size()/6) ;
	epsilon11.resize(sigma.size()/6) ;
	epsilon22.resize(sigma.size()/6) ;
	epsilon12.resize(sigma.size()/6) ;
	vonMises.resize(sigma.size()/6) ;

	std::cout << "unknowns :" << x.size() << std::endl ;
	std::cout << "max value :" << x.max() << std::endl ;
	std::cout << "min value :" << x.min() << std::endl ;


	for(size_t k = 0 ; k < triangles.size() ; k++)
	{

		if(!triangles[k]->getBehaviour()->fractured())
		{
			sigma11[k*4] = sigma[k*6*4];
			sigma22[k*4] = sigma[k*6*4+1];
			sigma12[k*4] = sigma[k*6*4+3];
			sigma11[k*4+1] = sigma[k*6*4+6];
			sigma22[k*4+1] = sigma[k*6*4+7];
			sigma12[k*4+1] = sigma[k*6*4+9];
			sigma11[k*4+2] = sigma[k*6*4+12];
			sigma22[k*4+2] = sigma[k*6*4+13];
			sigma12[k*4+2] = sigma[k*6*4+15];

			epsilon11[k*4] = epsilon[k*6*4];
			epsilon22[k*4] = epsilon[k*6*4+1];
			epsilon12[k*4] = epsilon[k*6*4+3];
			epsilon11[k*3+1] = epsilon[k*6*4+6];
			epsilon22[k*4+1] = epsilon[k*6*4+7];
			epsilon12[k*4+1] = epsilon[k*6*4+8];
			epsilon11[k*4+2] = epsilon[k*6*4+12];
			epsilon22[k*4+2] = epsilon[k*6*4+13];
			epsilon12[k*4+2] = epsilon[k*6*4+15];

			Vector vm0 = triangles[k]->getState().getPrincipalStresses(*triangles[k]->first) ;
			Vector vm1 = triangles[k]->getState().getPrincipalStresses(*triangles[k]->second) ;
			Vector vm2 = triangles[k]->getState().getPrincipalStresses(*triangles[k]->third) ;
			vonMises[k*4]   = sqrt(((vm0[0]-vm0[1])*(vm0[0]-vm0[1]))/2.) ;
			vonMises[k*4+1] = sqrt(((vm1[0]-vm1[1])*(vm1[0]-vm1[1]))/2.) ;
			vonMises[k*4+2] = sqrt(((vm2[0]-vm2[1])*(vm2[0]-vm2[1]))/2.) ;
			vonMises[k*4+3] = sqrt(((vm2[0]-vm2[1])*(vm2[0]-vm2[1]))/2.) ;

// 			vonMises[k*3]   = sigma11[k*3]*epsilon11[k*3]     + sigma22[k*3]*epsilon22[k*3]     + sigma12[k*3]*epsilon12[k*3];
// 			vonMises[k*3+1] = sigma11[k*3+1]*epsilon11[k*3+1] + sigma22[k*3+1]*epsilon22[k*3+1] + sigma12[k*3+1]*epsilon12[k*3+2];
// 			vonMises[k*3+2] = sigma11[k*3+2]*epsilon11[k*3+2] + sigma22[k*3+2]*epsilon22[k*3+2] + sigma12[k*3+2]*epsilon12[k*3+2];
		}
		else
		{
			sigma11[k*4] = 0;
			sigma22[k*4] = 0;
			sigma12[k*4] = 0;
			sigma11[k*4+1] = 0;
			sigma22[k*4+1] = 0;
			sigma12[k*4+1] = 0;
			sigma11[k*4+2] = 0;
			sigma22[k*4+2] = 0;
			sigma12[k*4+2] = 0;

			epsilon11[k*4] = 0;
			epsilon22[k*4] = 0;
			epsilon12[k*4] = 0;
			epsilon11[k*3+1] = 0;
			epsilon22[k*4+1] = 0;
			epsilon12[k*4+1] = 0;
			epsilon11[k*4+2] = 0;
			epsilon22[k*4+2] = 0;
			epsilon12[k*4+2] = 0;

// 			Vector vm0 = triangles[k]->getState().getPrincipalStresses(*triangles[k]->first) ;
// 			Vector vm1 = triangles[k]->getState().getPrincipalStresses(*triangles[k]->second) ;
// 			Vector vm2 = triangles[k]->getState().getPrincipalStresses(*triangles[k]->third) ;
			vonMises[k*4]   = 0;
			vonMises[k*4+1] = 0;
			vonMises[k*4+2] = 0;
			vonMises[k*4+3] = 0;
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

	Sample3D sample(NULL, 0.04,0.04,0.04,0.02,0.02,0.02) ;
	FeatureTree F(&sample) ;
	featureTree = &F ;

	Matrix m0(6,6) ;
	m0[0][0] = 1. - nu ; m0[0][1] = nu ; m0[0][2] = nu ;
	m0[1][0] = nu ; m0[1][1] = 1. - nu ; m0[1][2] = nu ;
	m0[2][0] = nu ; m0[2][1] = nu ; m0[2][2] = 1. - nu ;
	m0[3][3] = 0.5 - nu ;
	m0[4][4] = 0.5 - nu ;
	m0[5][5] = 0.5 - nu ;
	m0 *= E/((1.+nu)*(1.-2.*nu)) ;
	

	double itzSize = 0.00001;
			int inclusionNumber = 3200 ;

			std::vector<Inclusion3D *> inclusions = GranuloBolome(0.0000384, 1, BOLOME_D)(true, .002, .0001, inclusionNumber, itzSize);

			if(inclusionNumber)
							itzSize = inclusions[inclusions.size()-1]->getRadius() ;
			for(size_t i = 0; i < inclusions.size() ; i++)
							delete inclusions[i] ;

			inclusions = GranuloBolome(0.0000384, 1, BOLOME_D)(true, .002, .0001, inclusionNumber, itzSize);

			std::vector<Feature *> feats ;
			for(size_t i = 0; i < inclusions.size() ; i++)
							feats.push_back(inclusions[i]) ;

			int nAgg = 3200 ;
			feats=placement(sample.getPrimitive(), feats, &nAgg, 64000);
	
	for(size_t i = 0 ; i < feats.size() ; i++)
	{
		F.addFeature(&sample, feats[i]) ;
		feats[i]->setBehaviour(new Stiffness(m0*4)) ;
//		std::cout << feats[i]->getRadius() << "   " << feats[i]->getCenter().x << "   "  << feats[i]->getCenter().y << "   " << feats[i]->getCenter().z << std::endl ; 
	}
// 	Inclusion3D inc(1, -1.5,0,0) ;
// 	F.addFeature(&sample, &inc) ;
// 	
// 
// 	Stiffness * s = new Stiffness(m0) ;
// 	inc.setBehaviour(s) ;
	MohrCoulomb * mc = new MohrCoulomb(30, -60) ;
	StiffnessAndFracture * sf = new StiffnessAndFracture(m0*0.5, mc) ;
	sample.setBehaviour(new Stiffness(m0)) ;

	F.sample(2048) ;
	F.setOrder(LINEAR) ;
	F.generateElements() ;
	step() ;
	
	return 0 ;
}
