// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//


#include "main.h"
#include "samplingcriterion.h"
#include "features/features.h"
#include "physics/physics.h"
#include "physics/mohrcoulomb.h"
#include "features/pore.h"
#include "features/sample.h"
#include "features/inclusion.h"
#include "features/expansiveZone.h"
#include "features/crack.h"
#include "features/enrichmentInclusion.h"
#include "delaunay_3d.h"
#include "solvers/assembly.h"

#include <fstream>

#include <cmath>
#include <typeinfo>
#include <limits>
#include <GL/glut.h>
#include <time.h> 
#define DEBUG 


using namespace Mu ;

FeatureTree * featureTree ;
std::vector<DelaunayTriangle *> triangles ;
DelaunayTree *dt ; //(pts) ;
std::vector<bool> cracked ;
std::vector<Crack *> crack ;

double E_min = 10;
double E_max = 0;

double x_max = 0 ;
double y_max = 0 ;

double x_min = 0 ;
double y_min = 0 ;

std::vector<DelaunayTriangle *> tris__ ;

std::pair<std::vector<Inclusion * >, std::vector<Pore * > > i_et_p ;

std::vector<ExpansiveZone *> zones ;



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
Vector angle(0) ; 



double nu = 0.3 ;
double E = 20 ;

bool nothingToAdd = false ;
bool dlist = false ;
int count = 0 ;


void setBC()
{
	triangles = featureTree->getTriangles() ;
	
	Matrix m0(3,3) ;
	m0[0][0] = E/(1-nu*nu) ; m0[0][1] =E/(1-nu*nu)*nu ; m0[0][2] = 0 ;
	m0[1][0] = E/(1-nu*nu)*nu ; m0[1][1] = E/(1-nu*nu) ; m0[1][2] = 0 ; 
	m0[2][0] = 0 ; m0[2][1] = 0 ; m0[2][2] = E/(1-nu*nu)*(1.-nu)/2. ; 
	
	Matrix m1 = m0*0.1 ;
	
	if(firstRun)
	{
		for(size_t k = 0 ; k < triangles.size() ;k++)
		{
// 			if(triangles[k]->getBehaviour()->type !=VOID_BEHAVIOUR)
// 				triangles[k]->setBehaviour(new StiffnessAndFracture(m0, 0.1)) ;
// 			if(triangles[k]->getBehaviour()->type !=VOID_BEHAVIOUR && 
// 			   !featureTree->getFeature(0)->getChild(0)->in(triangles[k]->getCenter()))
// 			{
// 				triangles[k]->setBehaviour(new StiffnessAndFracture(m0, 1)) ;
// // 				triangles[k]->setBehaviour(new VoidForm()) ;
// 			}
// 			if(triangles[k]->getBehaviour()->type !=VOID_BEHAVIOUR && (
// 				featureTree->getFeature(0)->getChild(1)->in(triangles[k]->getCenter()) ||
// 				featureTree->getFeature(0)->getChild(0)->in(triangles[k]->getCenter())
// 				)
// 			  )
// 				triangles[k]->setBehaviour(new VoidForm()) ;
// 			if(triangles[k]->getBehaviour()->type !=VOID_BEHAVIOUR && 
// 			   featureTree->getFeature(0)->getChild(2)->in(triangles[k]->getCenter()))
// 				triangles[k]->setBehaviour(new VoidForm()) ;
		}
		
		firstRun = false ;
	}
	
	for(size_t k = 0 ; k < triangles.size() ;k++)
	{
		for(size_t c = 0 ;  c < triangles[k]->getBoundingPoints().size() ; c++ )
		{
			
// 			if ((std::abs(triangles[k]->getBoundingPoint(c).x +3.5) < 0.1 
// 			     || std::abs(triangles[k]->getBoundingPoint(c).x -3.5) < 0.1) 
// 			    && triangles[k]->getBoundingPoint(c).y < -.9999)
// 			{
// 				featureTree->getAssembly()->setPointAlong( ETA,0 ,triangles[k]->getBoundingPoint(c).id) ;
// 			}
// 			
// 			if (std::abs(triangles[k]->getBoundingPoint(c).x ) < 0.1 && triangles[k]->getBoundingPoint(c).y > .9999)
// 			{
// 				featureTree->getAssembly()->setPointAlong( ETA,-timepos ,triangles[k]->getBoundingPoint(c).id) ;
// 			}

// 			if (std::abs(triangles[k]->getBoundingPoint(c).x) < .1 && triangles[k]->getBoundingPoint(c).y > 0.9999 )
// 			{
// 				featureTree->getAssembly()->setForceOn(ETA, -timepos- 0.00001, triangles[k]->getBoundingPoint(c).id) ;
// // 				featureTree->getAssembly()->setPointAlong( ETA,-timepos ,triangles[k]->getBoundingPoint(c).id) ;
// 			}

			if(triangles[k]->getBoundingPoint(c).x < -3.999)
			{
				featureTree->getAssembly()->setPointAlong( XI,0, triangles[k]->getBoundingPoint(c).id) ;
			}
			if (triangles[k]->getBoundingPoint(c).y < -0.999 )
			{
				featureTree->getAssembly()->setPointAlong( ETA,0 ,triangles[k]->getBoundingPoint(c).id) ;
			}
// 			if(triangles[k]->getBoundingPoint(c).y > 2.999 && triangles[k]->getBoundingPoint(c).x > 2.999)
// 			{
// 				featureTree->getAssembly()->setPoint( 0,0, triangles[k]->getBoundingPoint(c).id) ;
// 			}
// 			if (triangles[k]->getBoundingPoint(c).y < -2.999)
// 			{
// 				featureTree->getAssembly()->setPointAlong( ETA, 0, triangles[k]->getBoundingPoint(c).id) ;
// 			}
// 			if(triangles[k]->getBoundingPoint(c).y > 2.999)
// 			{
// 				featureTree->getAssembly()->setPointAlong( ETA, 0, triangles[k]->getBoundingPoint(c).id) ;
// 			}
			
// 			if (triangles[k]->getBoundingPoint(c).y <-2.9999)
// 			{
// 				F.getAssembly()->setPoint( 0, -.500,triangles[k]->getBoundingPoint(c).id) ;
// 			}
// 			if(triangles[k]->getBoundingPoint(c).y > 2.99999 )
// 			{
// 				force++ ;
// 			}
		}
// 		if(force == 3)
// 			for(size_t c = 0 ;  c < triangles[k]->getBoundingPoints().size() ; c++ )
// 			{
// 				if(triangles[k]->getBoundingPoint(c).y > 2.99999 )
// 				{
// 					F.getAssembly()->setForceOn(ETA, -0.01, triangles[k]->getBoundingPoint(c).id) ;
// 				}
// 			}
		
	}
	
// 	for(size_t i = 0 ; i < i_et_p.first.size() ;i++)
// 	{
// 		applyForcesInInclusion(4, i_et_p.first[i]) ;
// 	}
}

void step()
{
	
	for(size_t i = 0 ; i < 1 ; i++)
	{
		std::cout << "\r iteration " << i << "/10" << std::flush ;
		setBC() ;
		while(!featureTree->step(timepos))
		{
			break ;
// 			timepos-= 0.0001 ;
			setBC() ;
			
		}
// 		
// 		
		timepos+= 0.0001 ;
	
	
	std::cout << "\r iteration " << "100/100 ... done" << std::endl ;
	x.resize(featureTree->getDisplacements().size()) ;
	x = featureTree->getDisplacements() ;
	dt = featureTree->getDelaunayTree() ;
	sigma.resize(triangles.size()*triangles[0]->getBoundingPoints().size()*3) ;
	epsilon.resize(triangles.size()*triangles[0]->getBoundingPoints().size()*3) ;
	
// 	sigma = F.strainFromDisplacements() ;
// 	epsilon = F.stressFromDisplacements() ;
	std::pair<Vector, Vector > sigma_epsilon = featureTree->getStressAndStrain() ;
	sigma.resize(sigma_epsilon.first.size()) ;
	sigma = sigma_epsilon.first ;
	epsilon.resize(sigma_epsilon.second.size()) ;
	epsilon = sigma_epsilon.second ;
	
	sigma11.resize(sigma.size()/3) ;
	sigma22.resize(sigma.size()/3) ;
	sigma12.resize(sigma.size()/3) ;
	epsilon11.resize(sigma.size()/3) ;
	epsilon22.resize(sigma.size()/3) ;
	epsilon12.resize(sigma.size()/3) ;
	vonMises.resize(sigma.size()/3) ;
	angle.resize(sigma.size()/3) ;
	
	std::cout << "unknowns :" << x.size() << std::endl ;
	
	if(crack.size() > 0)
		tris__ = crack[0]->getIntersectingTriangles(dt) ;
	
	for(size_t k = 1 ; k < crack.size() ; k++)
	{
		std::vector<DelaunayTriangle *> temp = crack[k]->getIntersectingTriangles(dt) ;
		if(tris__.empty())
			tris__ = temp ;
		else if(!temp.empty())
			tris__.insert(tris__.end(), temp.begin(), temp.end() ) ;
	}
	cracked.clear() ;
	
	int npoints = triangles[0]->getBoundingPoints().size() ;
	
	double area = 0 ;
	double avg_e_xx = 0;
	double avg_e_yy = 0;
	double avg_e_xy = 0;
	double avg_s_xx = 0;
	double avg_s_yy = 0;
	double avg_s_xy = 0;
	
	for(size_t k = 0 ; k < triangles.size() ; k++)
	{
/*		bool in = !triangles[k]->getEnrichmentFunctions().empty() ;*/
		bool in = false ;
		for(size_t m = 0 ; m < tris__.size() ; m++)
		{
			if(triangles[k] == tris__[m])
			{
				in = true ;
				break ;
			}
		}
		cracked.push_back(in) ;
		
		if(!in && !triangles[k]->getBehaviour()->fractured())
		{
			for(size_t p = 0 ;p < triangles[k]->getBoundingPoints().size() ; p++)
			{
				if(x[triangles[k]->getBoundingPoint(p).id*2] > x_max)
					x_max = x[triangles[k]->getBoundingPoint(p).id*2];
				if(x[triangles[k]->getBoundingPoint(p).id*2] < x_min)
					x_min = x[triangles[k]->getBoundingPoint(p).id*2];
				if(x[triangles[k]->getBoundingPoint(p).id*2+1] > y_max)
					y_max = x[triangles[k]->getBoundingPoint(p).id*2+1];
				if(x[triangles[k]->getBoundingPoint(p).id*2+1] < y_min)
					y_min = x[triangles[k]->getBoundingPoint(p).id*2+1];
			}
			area += triangles[k]->area() ;
			if(triangles[k]->getBehaviour()->type != VOID_BEHAVIOUR)
			{
				if(triangles[k]->getBehaviour()->param[0][0] > E_max)
					E_max = triangles[k]->getBehaviour()->param[0][0] ;
				if(triangles[k]->getBehaviour()->param[0][0] < E_min)
					E_min = triangles[k]->getBehaviour()->param[0][0] ;
			}
// 			std::cout << sigma[k*npoints*3] << "   " 
// 				<< sigma[k*npoints*3+1] << "   " 
// 				<< sigma[k*npoints*3+2] << "   " 
// 				<< sigma[k*npoints*3+3] << "   " 
// 				<< sigma[k*npoints*3+4] << "   " 
// 				<< sigma[k*npoints*3+5] << "   " 
// 				<< sigma[k*npoints*3+6] << "   " 
// 				<< sigma[k*npoints*3+7] << "   " 
// 				<< sigma[k*npoints*3+8] << "   " 
// 				<< sigma[k*npoints*3+9] << "   " 
// 				<< sigma[k*npoints*3+10] << "   " 
// 				<< sigma[k*npoints*3+11] << "   " 
// 				<< sigma[k*npoints*3+12] << "   " 
// 				<< sigma[k*npoints*3+13] << "   " 
// 				<< sigma[k*npoints*3+14] << "   " 
// 				<< sigma[k*npoints*3+15] << "   " 
// 				<< sigma[k*npoints*3+16] << "   " 
// 				<< sigma[k*npoints*3+17] << std::endl ;
				
			sigma11[k*npoints] = sigma[k*npoints*3];
			sigma22[k*npoints] = sigma[k*npoints*3+1];
			sigma12[k*npoints] = sigma[k*npoints*3+2];
			sigma11[k*npoints+1] = sigma[k*npoints*3+3];
			sigma22[k*npoints+1] = sigma[k*npoints*3+4];
			sigma12[k*npoints+1] = sigma[k*npoints*3+5];
			sigma11[k*npoints+2] = sigma[k*npoints*3+6];
			sigma22[k*npoints+2] = sigma[k*npoints*3+7];
			sigma12[k*npoints+2] = sigma[k*npoints*3+8];
			
			if(npoints >3)
			{
				sigma11[k*npoints+3] = sigma[k*npoints*3+9];
				sigma22[k*npoints+3] = sigma[k*npoints*3+10];
				sigma12[k*npoints+3] = sigma[k*npoints*3+11];
				sigma11[k*npoints+4] = sigma[k*npoints*3+12];
				sigma22[k*npoints+4] = sigma[k*npoints*3+13];
				sigma12[k*npoints+4] = sigma[k*npoints*3+14];
				sigma11[k*npoints+5] = sigma[k*npoints*3+15];
				sigma22[k*npoints+5] = sigma[k*npoints*3+16];
				sigma12[k*npoints+5] = sigma[k*npoints*3+17];
			}
			
			epsilon11[k*npoints] = epsilon[k*npoints*3];
			epsilon22[k*npoints] = epsilon[k*npoints*3+1];
			epsilon12[k*npoints] = epsilon[k*npoints*3+2];
			epsilon11[k*npoints+1] = epsilon[k*npoints*3+3];
			epsilon22[k*npoints+1] = epsilon[k*npoints*3+4];
			epsilon12[k*npoints+1] = epsilon[k*npoints*3+5];
			epsilon11[k*npoints+2] = epsilon[k*npoints*3+6];
			epsilon22[k*npoints+2] = epsilon[k*npoints*3+7];
			epsilon12[k*npoints+2] = epsilon[k*npoints*3+8];
			
			if(npoints > 3)
			{
				epsilon11[k*npoints+3] = epsilon[k*npoints*3+9];
				epsilon22[k*npoints+3] = epsilon[k*npoints*3+10];
				epsilon12[k*npoints+3] = epsilon[k*npoints*3+11];
				epsilon11[k*npoints+4] = epsilon[k*npoints*3+12];
				epsilon22[k*npoints+4] = epsilon[k*npoints*3+13];
				epsilon12[k*npoints+4] = epsilon[k*npoints*3+14];
				epsilon11[k*npoints+5] = epsilon[k*npoints*3+15];
				epsilon22[k*npoints+5] = epsilon[k*npoints*3+16];
				epsilon12[k*npoints+5] = epsilon[k*npoints*3+17];
			}  
			
			for(size_t l = 0 ; l < triangles[k]->getBoundingPoints().size() ; l++)
			{
				Vector vm0 = triangles[k]->getState()->getPrincipalStresses(triangles[k]->getBoundingPoint(l)) ;
				vonMises[k*triangles[k]->getBoundingPoints().size()+l]  = sqrt(((vm0[0]-vm0[1])*(vm0[0]-vm0[1]))/2.) ;

				double agl = triangles[k]->getState()->getPrincipalAngle(triangles[k]->getBoundingPoint(l)) ;
				angle[k*triangles[k]->getBoundingPoints().size()+l]  = agl ;
			}
			
			double ar = triangles[k]->area() ;
			for(size_t l = 0 ; l < npoints ;l++)
			{
				avg_e_xx += (epsilon11[k*npoints+l]/npoints)*ar;
				avg_e_yy += (epsilon22[k*npoints+l]/npoints)*ar;
				avg_e_xy += (epsilon12[k*npoints+l]/npoints)*ar;
				avg_s_xx += (sigma11[k*npoints+l]/npoints)*ar;
				avg_s_yy += (sigma22[k*npoints+l]/npoints)*ar;
				avg_s_xy += (sigma12[k*npoints+l]/npoints)*ar;
			}
				

			
// 			vonMises[k*3]   = sigma11[k*3]*epsilon11[k*3]     + sigma22[k*3]*epsilon22[k*3]     + sigma12[k*3]*epsilon12[k*3];
// 			vonMises[k*3+1] = sigma11[k*3+1]*epsilon11[k*3+1] + sigma22[k*3+1]*epsilon22[k*3+1] + sigma12[k*3+1]*epsilon12[k*3+2];
// 			vonMises[k*3+2] = sigma11[k*3+2]*epsilon11[k*3+2] + sigma22[k*3+2]*epsilon22[k*3+2] + sigma12[k*3+2]*epsilon12[k*3+2];
		}
		else
		{
			sigma11[k*npoints] = 0 ;
			sigma22[k*npoints] = 0 ;
			sigma12[k*npoints] = 0 ;
			sigma11[k*npoints+1] = 0 ;
			sigma22[k*npoints+1] = 0 ;
			sigma12[k*npoints+1] = 0 ;
			sigma11[k*npoints+2] = 0 ;
			sigma22[k*npoints+2] = 0 ;
			sigma12[k*npoints+2] = 0 ;
			
			if(npoints >3)
			{
				sigma11[k*npoints+3] = 0 ;
				sigma22[k*npoints+3] = 0 ;
				sigma12[k*npoints+3] = 0 ;
				sigma11[k*npoints+4] = 0 ;
				sigma22[k*npoints+4] = 0 ;
				sigma12[k*npoints+4] = 0 ;
				sigma11[k*npoints+5] = 0 ;
				sigma22[k*npoints+5] = 0 ;
				sigma12[k*npoints+5] =0 ;
			}
			
			epsilon11[k*npoints] = 0 ;
			epsilon22[k*npoints] = 0 ;
			epsilon12[k*npoints] = 0 ;
			epsilon11[k*npoints+1] = 0 ;
			epsilon22[k*npoints+1] = 0 ;
			epsilon12[k*npoints+1] = 0 ;
			epsilon11[k*npoints+2] = 0 ;
			epsilon22[k*npoints+2] = 0 ;
			epsilon12[k*npoints+2] = 0 ;
			
			if(npoints > 3)
			{
				epsilon11[k*npoints+3] = 0 ;
				epsilon22[k*npoints+3] = 0 ;
				epsilon12[k*npoints+3] =0 ;
				epsilon11[k*npoints+4] = 0 ;
				epsilon22[k*npoints+4] = 0 ;
				epsilon12[k*npoints+4] =0 ;
				epsilon11[k*npoints+5] = 0 ;
				epsilon22[k*npoints+5] =0 ;
				epsilon12[k*npoints+5] = 0 ;
			}  
			
			for(size_t l = 0 ; l < triangles[k]->getBoundingPoints().size() ; l++)
			{
				vonMises[k*triangles[k]->getBoundingPoints().size()+l]  = 0 ;
				angle[k*triangles[k]->getBoundingPoints().size()+l]  = 0 ;
			}
		}
	}
		
	std::cout << "max value :" << x_max << std::endl ;
	std::cout << "min value :" << x_min << std::endl ;
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
	
	std::cout << "average sigma11 : " << avg_s_xx/area << std::endl ;
	std::cout << "average sigma22 : " << avg_s_yy/area << std::endl ;
	std::cout << "average sigma12 : " << avg_s_xy/area << std::endl ;
	std::cout << "average epsilon11 : " << avg_e_xx/area << std::endl ;
	std::cout << "average epsilon22 : " << avg_e_yy/area << std::endl ;
	std::cout << "average epsilon12 : " << avg_e_xy/area << std::endl ;
		
		double reactedArea = 0 ;
		
		for(size_t z = 0 ; z < zones.size() ; z++)
		{
			zones[z]->setRadius(zones[z]->getGeometry()->getRadius()+0.0001) ;
			reactedArea += zones[z]->area() ;
			zones[z]->reset() ;
		}
		
		std::cout << "reacted Area : " << reactedArea << std::endl ;
	}
}


int main(int argc, char *argv[])
{

  // create a sample of 8 by 2
	Sample sample(NULL, 12,2,0,0) ;


	FeatureTree F(&sample) ;
	featureTree = &F ;
	
	// material behaviour
	Matrix m0(3,3) ;
	m0[0][0] = E/(1-nu*nu) ; m0[0][1] =E/(1-nu*nu)*nu ; m0[0][2] = 0 ;
	m0[1][0] = E/(1-nu*nu)*nu ; m0[1][1] = E/(1-nu*nu) ; m0[1][2] = 0 ; 
	m0[2][0] = 0 ; m0[2][1] = 0 ; m0[2][2] = E/(1-nu*nu)*(1.-nu)/2. ; 

	// create inclusion centered 0,0 radius 0.5
	Inclusion * inc = new Inclusion(.5, 0,0) ;
	std::vector<Inclusion *> inclusions ;
	inclusions.push_back(inc) ;

	F.addFeature(&sample,inc) ;

	// set behaviour for whole sample
// 	sample.setBehaviour(new WeibullDistributedStiffness(m0*0.125, 0.02)) ;	
	sample.setBehaviour(new StiffnessAndFracture(m0*0.5, new MohrCoulomb(1,-2))) ;
	// set behaviour for inclusion
	inc->setBehaviour(new Stiffness(m0)) ;
	// specify point density
	F.sample(64) ;
	F.setOrder(LINEAR) ; // element type
	F.generateElements() ;// mesh

	step() ; // t -> t+dt
	
	return 0 ;
}
