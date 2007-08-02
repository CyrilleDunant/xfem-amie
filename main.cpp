// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
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

using namespace Mu ;

FeatureTree * featureTree ;
std::vector<DelaunayTriangle *> triangles ;
std::vector<bool> cracked ;
DelaunayTree *dt ; //(pts) ;
std::vector<Crack *> crack ;

double E_min = 10;
double E_max = 0;

double timepos = 0.0 ;

bool firstRun = true ;

std::vector<DelaunayTriangle *> tris__ ;

std::pair<std::vector<Inclusion * >, std::vector<Pore * > > i_et_p ;

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
double E = 2 ;

size_t current_list = DISPLAY_LIST_STRAIN_XX ;
double factor = 1 ;
MinimumAngle cri(M_PI/6.) ;
bool nothingToAdd = false ;
bool dlist = false ;
int count = 0 ;

// void makeDisplacementsGlobal(std::vector<DelaunayTriangle *> * tri, size_t start, Vector * eps)
// {
// 	(*tri)[start]->displace(eps) ;
// }

/** Les fonctions de Mohsen...
 * 
 */

double generateE( double surface_fissure, double rayon_fissure, double nombre_fissure, double nu ,double E_0, double sigma_max, double sigma, double sigma_moy)
{
	double Q = 225.*(1.-2.*nu)*(nu-2.)*(nu-2.)*(1.+nu)/(16.*(1.-nu)) ;
	double E_max = 1./(Q*nombre_fissure*rayon_fissure*rayon_fissure*rayon_fissure/surface_fissure - 1./E_0) ;
	return E_max - (E_max - E_0)*(sigma/sigma_max -1)*(sigma/sigma_max -1)*sigma_moy ;
}

struct fissureSimple
{
	fissureSimple(double r,Point c,double p,double e)
	{
		rayon = r ;
		centre  = c ; 
		phi = p ;
		epaisseur = e ;
	}
	
	fissureSimple()
	{
		rayon = 0 ;
		centre  = Point() ; 
		phi = 0 ;
		epaisseur = 0 ;
	}
	
	
	double rayon ;
	Point centre ;
	double phi ;
	double epaisseur ;
	
	double surface() const
	{
		return M_PI*rayon*epaisseur ;
	}
	
	void applyForces(Vector * stresses, Vector *b)
	{
		std::cout << "demander la fonction magique de Mohsen" << std::endl ;
		exit(0) ;
	}
} ;

void applyForcesInInclusion(size_t n, Inclusion * i)
{
	size_t tot = 0 ;
	for(size_t j =0 ; j < n ; j++)
	{
		if(i->getRadius() > 0.045 + 0.045)
		{
			Point center = i->Circle::getCenter() + 
				Point(
				       (2.*random()/RAND_MAX-1.),
				       (2.*random()/RAND_MAX-1.)
				     )*0.55*i->getRadius() ;
			double radius = (0.045 + 0.035*random()/RAND_MAX) ;
			
			Circle c(radius, center) ;
			
			bool allin = true ;
			
			
			std::vector<DelaunayTriangle *> tris =featureTree->getDelaunayTree()->conflicts(&c) ;
			
			
			for(size_t k = 0 ; k < tris.size() ; k++)
			{
				if(c.intersects(dynamic_cast<Triangle *>(tris[k])))
				{
					allin = false ;
					break ;
				}
			}
			if(allin) ;
			{
				tot++ ;
				for(size_t k = 0 ; k < tris.size() ; k++)
				{
					for(size_t l = 0 ; l < tris[k]->getBoundingPoints().size() ; l++)
					{
						if(c.in(tris[k]->getBoundingPoint(l)))
						{
							Point vector = (c.getCenter() - tris[k]->getBoundingPoint(l))*timepos ;
							featureTree->getAssembly()->setForceOn(XI,  vector.x, tris[k]->getBoundingPoint(l).id) ;
							featureTree->getAssembly()->setForceOn(ETA, vector.y, tris[k]->getBoundingPoint(l).id) ;
						}
					}
				}
			}
		}
	}
	
// 	std::cout << "applied " << tot << " force zones" << std::endl ;
} ;

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
			
			if (std::abs(triangles[k]->getBoundingPoint(c).x +4) < 0.001)
			{
				featureTree->getAssembly()->setPointAlong( XI,0 ,triangles[k]->getBoundingPoint(c).id) ;
			}
			
			if (std::abs(triangles[k]->getBoundingPoint(c).y +2) < 0.001)
			{
				featureTree->getAssembly()->setPointAlong( ETA,0 ,triangles[k]->getBoundingPoint(c).id) ;
			}

// 			if (std::abs(triangles[k]->getBoundingPoint(c).x) < .1 && triangles[k]->getBoundingPoint(c).y > 0.9999 )
// 			{
// 				featureTree->getAssembly()->setForceOn(ETA, -timepos- 0.00001, triangles[k]->getBoundingPoint(c).id) ;
// // 				featureTree->getAssembly()->setPointAlong( ETA,-timepos ,triangles[k]->getBoundingPoint(c).id) ;
// 			}

// 			if(triangles[k]->getBoundingPoint(c).x < -2.999)
// 			{
// 				featureTree->getAssembly()->setPointAlong( XI,0, triangles[k]->getBoundingPoint(c).id) ;
// 			}
// 			if (triangles[k]->getBoundingPoint(c).y < -2.999 && triangles[k]->getBoundingPoint(c).x > 2.999)
// 			{
// 				featureTree->getAssembly()->setPoint( 0,0 ,triangles[k]->getBoundingPoint(c).id) ;
// 			}
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
	
	for(size_t i = 0 ; i < 100 ; i++)
	{
		std::cout << "\r iteration " << i << "/100" << std::flush ;
		setBC() ;
		featureTree->step(timepos) ;
		
// 		timepos+= 0.01 ;
	}
	
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
	std::cout << "max value :" << x.max() << std::endl ;
	std::cout << "min value :" << x.min() << std::endl ;
	
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

std::vector<fissureSimple> generateInInclusion(size_t n, Inclusion * i)
{
	std::vector<fissureSimple> ret ;
	for(size_t j =0 ; j < n ; j++)
	{
		Point center = i->Circle::getCenter() + Point(i->getRadius()*0.9*random()/RAND_MAX, i->getRadius()*0.9*random()/RAND_MAX) ;
		double radius = 0.09*i->getRadius()*random()/RAND_MAX ;
		double phi = 2*M_PI*random()/RAND_MAX - M_PI ;
		double e = 0.001*radius ;
		
		bool alone  = true ;
		
		for(size_t k = 0 ; k < ret.size() ; k++ )
		{
			if (squareDist(center, ret[k].centre) < (radius+ret[k].rayon)*(radius+ret[k].rayon))
			{
				alone = false ;
				break ;
			}
		}
		if (alone)
			ret.push_back(fissureSimple(radius,center, phi, e)) ;
		else
			j-- ;
	}
	
	return ret ;
} ;

std::vector<Pore * > generatePores(size_t n)
{
	std::vector<Pore *> ret ;
	for(size_t j =0 ; j < n ; j++)
	{

		double radius = 0.1 + 0.6*random()/RAND_MAX ;
		
		Point center = Point(
		                      (2.*random()/RAND_MAX-1.),
		                      (2.*random()/RAND_MAX-1.)
		                    )*(3./*-2.*radius-0.05*/) ; 
		
		bool alone  = true ;
		
		for(size_t k = 0 ; k < ret.size() ; k++ )
		{
			if (squareDist(center, ret[k]->Circle::getCenter()) < (radius+ret[k]->Circle::getRadius()+0.02)*(radius+ret[k]->Circle::getRadius()+0.02))
			{
				alone = false ;
				break ;
			}
		}
		if (alone)
		{
			ret.push_back(new Pore(radius, center)) ;
			//(*ret.rbegin())->setStrainTensor(tensor) ;
		}
		else
			j-- ;
	}
	
	return ret ;
} ;

std::vector<Inclusion * > generateInclusions(size_t n, Matrix * tensor)
{
	std::vector<Inclusion *> ret ;
	for(size_t j =0 ; j < n ; j++)
	{
		
		double radius = 0.06 + 0.6*random()/RAND_MAX ;
		
		Point center = Point(
		                      (2.*random()/RAND_MAX-1.),
		                      (2.*random()/RAND_MAX-1.)
		                    )*(3./*-2.*radius-0.05*/) ; 
		
		bool alone  = true ;
		
		for(size_t k = 0 ; k < ret.size() ; k++ )
		{
			if (squareDist(center, ret[k]->Circle::getCenter()) < (radius+ret[k]->Circle::getRadius()+0.06)*(radius+ret[k]->Circle::getRadius()+0.06))
			{
				alone = false ;
				break ;
			}
		}
		if (alone)
		{
			ret.push_back(new Inclusion(radius, center)) ;
			(*ret.rbegin())->setBehaviour(new Stiffness(*tensor)) ;
		}
		else
			j-- ;
	}
	
	return ret ;
} ;

void generateExpansiveZones(int n, std::vector<Inclusion * > & incs , FeatureTree & F)
{
	double E = .1 ;
	double nu = .45 ;
	Matrix m0(3,3) ;
	m0[0][0] = E/(1-nu*nu) ; m0[0][1] =E/(1-nu*nu)*nu ; m0[0][2] = 0 ;
	m0[1][0] = E/(1-nu*nu)*nu ; m0[1][1] = E/(1-nu*nu) ; m0[1][2] = 0 ; 
	m0[2][0] = 0 ; m0[2][1] = 0 ; m0[2][2] = E/(1-nu*nu)*(1.-nu)/2. ; 
	
	
	for(size_t i = 0 ; i < incs.size() ; i++)
	{
		std::vector<Inclusion *> ret ;
		for(int j = 0 ; j < n ; j++)
		{
			double radius = 0.0005 ;
			
			Point center = incs[i]->getCenter()+Point(
			                      (2.*random()/RAND_MAX-1.),
			                      (2.*random()/RAND_MAX-1.)
			                                         )*(incs[i]->getRadius()*.6) ; 
			
			bool alone  = true ;
			
			for(size_t k = 0 ; k < ret.size() ; k++ )
			{
				if (squareDist(center, ret[k]->Circle::getCenter()) < (radius+ret[k]->Circle::getRadius()+0.06)*(radius+ret[k]->Circle::getRadius()+0.06))
				{
					alone = false ;
					break ;
				}
			}
			if (alone)
			{
				Vector a(double(0), 3) ;
				a[0] = 0.5 ;
				a[1] = 0.5 ;
				a[2] = 0.00 ;
				
				ret.push_back(new Inclusion(radius, center)) ;
				(*ret.rbegin())->setBehaviour(new StiffnessWithImposedDeformation(m0,a)) ;
				F.addFeature(incs[i],*ret.rbegin()) ; 
			}
		}
	}
}

std::vector<Crack *> generateCracks(size_t n)
{
	std::vector<Crack *> ret ;
	std::vector<Circle *> pos ;
	size_t nit = 0 ;
	for(size_t j =0 ; j < n && nit < 2048; j++)
	{
		nit++ ;
		double radius = 0.1 + 0.5*random()/RAND_MAX;
		
		Point center = Point(
		                      (2.*random()/RAND_MAX-1.),
		                      (2.*random()/RAND_MAX-1.)
		                    )*(3. - .2*radius ) ; 

		
		bool alone  = true ;
		
		for(size_t k = 0 ; k < pos.size() ; k++ )
		{
			if (squareDist(center, pos[k]->getCenter()) <
			    (radius+pos[k]->getRadius()+0.05)*(radius+pos[k]->getRadius()+0.05))
			{
				alone = false ;
				break ;
			}
		}
		if (alone)
		{
			pos.push_back(new Circle(radius, center)) ;
		}
		else
			j-- ;
		
// 		pos.push_back(new Circle(radius, center)) ;

	}
	
	for(size_t j = 0 ; j < pos.size() ; j++)
	{
		std::valarray<Point *> ptset1(2) ;
		double angle = (2.*random()/RAND_MAX-1.)*M_PI ;
		double x_0 = pos[j]->getCenter().x + pos[j]->getRadius()*cos(angle);
		double y_0 = pos[j]->getCenter().y + pos[j]->getRadius()*sin(angle);
		double x_1 = pos[j]->getCenter().x + pos[j]->getRadius()*cos(angle+M_PI) ;
		double y_1 = pos[j]->getCenter().y + pos[j]->getRadius()*sin(angle+M_PI);
		
		ptset1[0] = new Point(x_0, y_0) ;
		ptset1[1] = new Point(x_1, y_1) ;
		ret.push_back(new Crack(ptset1, 0.02)) ;
	}
	std::cout << "placed " << ret.size() << " cracks" << std::endl ;
	return ret ;
} ;

std::pair<std::vector<Inclusion * >, std::vector<Pore * > > generateInclusionsAndPores(size_t n, double fraction, Matrix * tensor, Feature * father, FeatureTree * F)
{
// 	srandom(time(NULL)) ;
	size_t nombre_de_pores = static_cast<size_t>(round(n*fraction)) ;
	size_t nombre_d_inclusions = static_cast<size_t>(round(n*(1. - fraction))) ;
	
	std::pair<std::vector<Inclusion * >, std::vector<Pore * > > ret ;
	ret.first = std::vector<Inclusion * >() ;
	ret.second = std::vector<Pore * >() ;
	
	std::vector<Circle *> cercles ;
	for(size_t j =0 ; j < n ; j++)
	{
		
		double radius = 0.01 + 0.5*random()/RAND_MAX ;
		
		Point center = Point(
		                      (2.*random()/RAND_MAX-1.)*(4.-2.*radius-0.001),
		                      (2.*random()/RAND_MAX-1.)*(2.-2.*radius-0.001)
		                    ); 
		bool alone  = true ;
		
		for(size_t k = 0 ; k < cercles.size() ; k++ )
		{
			if (squareDist(center, cercles[k]->getCenter()) < (radius+cercles[k]->getRadius()+0.001)*(radius+cercles[k]->getRadius()+0.001))
			{
				alone = false ;
				break ;
			}
		}
		if (alone)
		{
			cercles.push_back(new Circle(radius, center)) ;
		}
		else
			j-- ;
		
	}
	for(size_t j =0 ; j < nombre_d_inclusions ; j++)
	{
		Vector imp(double(0),3) ;
		imp[0] = 0.01 ;
		imp[1] = 0.01 ;
		Inclusion * temp = new Inclusion(cercles[j]->getRadius(), cercles[j]->getCenter()) ;
		ret.first.push_back(temp) ;
		
		(*ret.first.rbegin())->setBehaviour(new WeibullDistributedStiffness(*tensor, 0.04)) ;
		F->addFeature(father, temp) ;
	}
	
	for(size_t j =0 ; j < nombre_de_pores ; j++)
	{
		Pore * temp = new Pore(cercles[j+nombre_d_inclusions]->getRadius(), cercles[j+nombre_d_inclusions]->getCenter()) ;
		ret.second.push_back(temp) ;
		F->addFeature(father, temp) ;
	}
	
	for(size_t k = 0 ; k < cercles.size() ; k++ )
	{
		delete cercles[k] ;
	}
	
	return ret ;
}

std::vector<std::vector<Inclusion *> > generateFissuresIn(std::vector<Inclusion * > inclusions, Matrix * cg)
{
	std::vector<std::vector<Inclusion *> > ret ;
	for(size_t j =0 ; j < inclusions.size() ; j++)
	{
		size_t num_fissures = (size_t)std::ceil(inclusions[j]->Circle::getRadius()/0.2) ;
		std::vector<Inclusion *> incs ;
		for(size_t k = 0 ; k < num_fissures ; k++)
		{

			Point center = inclusions[j]->Circle::getCenter() + 
				Point( (2.*random()/RAND_MAX-1.), (2.*random()/RAND_MAX-1.) )*0.55*inclusions[j]->getRadius() ;
			double radius = (0.055 + 0.045*random()/RAND_MAX) ;
			
			Inclusion * new_inc = new Inclusion(inclusions[j], radius, center) ;
			incs.push_back(new_inc) ;
		}
		ret.push_back(incs) ;
	}
	
	return ret ;
} ;

void HSVtoRGB( double *r, double *g, double *b, double h, double s, double v )
{
	int i;
	double f, p, q, t;
	if( s == 0 ) {
                // achromatic (grey)
		*r = *g = *b = v;
		return;
	}
	h /= 60.;                        // sector 0 to 5
	i = (int)floor( h );
	f = h - i;                      // factorial part of h
	p = v * ( 1. - s );
	q = v * ( 1. - s * f );
	t = v * ( 1. - s * ( 1. - f ) );
	switch( i ) {
	case 0:
		*r = v;
		*g = t;
		*b = p;
		break;
	case 1:
		*r = q;
		*g = v;
		*b = p;
		break;
	case 2:
		*r = p;
		*g = v;
		*b = t;
		break;
	case 3:
		*r = p;
		*g = q;
		*b = v;
		break;
	case 4:
		*r = t;
		*g = p;
		*b = v;
		break;
	default:                // case 5:
		*r = v;
		*g = p;
		*b = q;
		break;
	}
}

void init(void) 
{
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glShadeModel(GL_SMOOTH);   // Enables Smooth Shading
	glEnable(GL_LINE_SMOOTH) ;
	glEnable(GL_POLYGON_SMOOTH);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST) ;
	
// 	glPointSize(std::max(0.4*(double)width()/(double)columns, 1.));
	glClearColor(0.0f,0.0f,0.0f,0.0f);                                      // Black Background

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE); 
}

void Menu(int selection)
{
	switch (selection)
	{
	case ID_NEXT:
		{
			step() ;
			dlist = false ;
			break ;
		}
	case ID_NEXT_TIME:
		{
			timepos +=0.0001 ;
		}
	case ID_DISP : 
		{
			current_list = DISPLAY_LIST_DISPLACEMENT ;
			break ;
		}
	case ID_STIFNESS : 
	{
		current_list = DISPLAY_LIST_STIFFNESS ;
		break ;
	}
	case ID_STRAIN_XX : 
		{
			current_list = DISPLAY_LIST_STRAIN_XX ;
			break ;
		}
	case ID_STRAIN_YY : 
		{
			current_list = DISPLAY_LIST_STRAIN_YY ;
			break ;
		}
	case ID_STRAIN_XY : 
		{
			current_list = DISPLAY_LIST_STRAIN_XY ;
			break ;
		}
	case ID_STRESS_XX : 
		{
			current_list = DISPLAY_LIST_STRESS_XX ;
			break ;
		}
	case ID_STRESS_YY : 
		{
			current_list = DISPLAY_LIST_STRESS_YY ;
			break ;
		}
	case ID_STRESS_XY : 
		{
			current_list = DISPLAY_LIST_STRESS_XY ;
			break ;
		}
	case ID_ELEM : 
		{
			current_list = DISPLAY_LIST_ELEMENTS ;
			break ;
		}
	case ID_VON_MISES: 
		{
			current_list = DISPLAY_LIST_VON_MISES ;
			break ;
		}
	case ID_ANGLE: 
		{
			current_list = DISPLAY_LIST_ANGLE ;
			break ;
		}	
	case ID_ENRICHMENT: 
		{
			current_list = DISPLAY_LIST_ENRICHMENT ;
			break ;
		}

	case ID_QUIT : exit(0) ;
		
	case ID_ZOOM :
		{
			factor *= 1.5 ;
			break ;
		}
	case ID_UNZOOM :
		{
			factor /= 1.5 ;
			break ;
		}
		
	case ID_AMPLIFY :
		{
			x *= 1.5 ;
			sigma11 *= 1.5 ;
			sigma22 *= 1.5 ;
			sigma12 *= 1.5 ;
			dlist = false ;
			break ;
		}
	case ID_DEAMPLIFY :
		{
			x /= 1.5 ;
			sigma11 /= 1.5 ;
			sigma22 /= 1.5 ;
			sigma12 /= 1.5 ;
			dlist = false ;
			break ;
		}
	}
}

void reshape(int w, int h)
{
	if (h==0)                                          // Prevent A Divide By Zero By
		h=1;                                           // Making Height Equal One
	
	glViewport(0, 0, (int)(w*factor), (int)(h*factor));
	gluPerspective((double)h/(double)w,1.,1.f,45.0f);
}

void Display(void)
{
	
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT) ;
	glFlush();
	glutSwapBuffers();
	glMatrixMode(GL_PROJECTION) ;
	glLoadIdentity() ;
	glOrtho(-4.5/factor, 4.5/factor, -4.5/factor, 4.5/factor, -4.5, 4.5);
// 	glEnable( GL_POLYGON_OFFSET_FILL );
// 	glPolygonOffset( 0.5, 0.5 );
	
	//std::cout << x.max() << std::endl ;
	//std::cout << x.min() << std::endl ;
	
	double x_max = std::abs(x).min() ;
	double y_max = std::abs(x).min() ;
	
	for(size_t k = 0 ; k < x.size()/2 ; k++)
	{
		if(x[k*2]*x[k*2]+x[k*2+1]*x[k*2+1] > x_max*x_max+y_max*y_max )
		{
			x_max = x[k*2] ;
			y_max = x[k*2+1] ;
		}
	}
	
	if(!dlist)
	{
		
		
		glNewList( DISPLAY_LIST_DISPLACEMENT,  GL_COMPILE ) ;
			for (unsigned int j=0 ; j< triangles.size() ; j++ )
			{

				if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !cracked[j] && !triangles[j]->getBehaviour()->fractured())
				{
					double c1 ;
					double c2 ;
					double c3 ;
		
					double vx = x[triangles[j]->getBoundingPoint(0).id*2]; 
					double vy = x[triangles[j]->getBoundingPoint(0).id*2+1]; 
					
					glBegin(GL_TRIANGLE_FAN);
					HSVtoRGB( &c1, &c2, &c3, 300. - sqrt((vx*vx + vy*vy)/(x_max*x_max + y_max*y_max))*300., 1., 1. ) ;
						glColor3f(c1, c2, c3) ;
							
						glVertex2f(double(triangles[j]->getBoundingPoint(0).x + vx) , double(triangles[j]->getBoundingPoint(0).y + vy) );
			
						for(size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++)
						{
							vx = x[triangles[j]->getBoundingPoint(k).id*2];
							vy = x[triangles[j]->getBoundingPoint(k).id*2+1]; 
						
							HSVtoRGB( &c1, &c2, &c3, 300. - sqrt((vx*vx + vy*vy)/(x_max*x_max + y_max*y_max))*300., 1., 1. ) ;
							glColor3f(c1, c2, c3) ;
							
							glVertex2f( double(triangles[j]->getBoundingPoint(k).x + vx) ,  double(triangles[j]->getBoundingPoint(k).y + vy) );
							
						}
					glEnd() ;
				}
			}
		glEndList() ;
		
		double sigma11_min = sigma11.min() ;
		double sigma11_max = sigma11.max() ;
		glNewList(  DISPLAY_LIST_STRAIN_XX,  GL_COMPILE ) ;
		
			for (unsigned int j=0 ; j< triangles.size() ; j++ )
			{
				
				if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !cracked[j]&& !triangles[j]->getBehaviour()->fractured())
				{
					 double c1 ;
					 double c2 ;
					 double c3 ;
					
					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(sigma11[j*triangles[j]->getBoundingPoints().size()]-sigma11_min)/(sigma11_max-sigma11_min), 1., 1.) ;
					glColor3f(c1, c2, c3) ;
					
					double vx = x[triangles[j]->first->id*2]; 
					double vy = x[triangles[j]->first->id*2+1]; 
					
					glBegin(GL_TRIANGLE_FAN);
					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(sigma11[j*triangles[j]->getBoundingPoints().size()]-sigma11_min)/(sigma11_max-sigma11_min), 1., 1.) ;
					
					glVertex2f(double(triangles[j]->getBoundingPoint(0).x + vx) , double(triangles[j]->getBoundingPoint(0).y + vy) );
					
					for(size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++)
					{
						vx = x[triangles[j]->getBoundingPoint(k).id*2];
						vy = x[triangles[j]->getBoundingPoint(k).id*2+1]; 
						
						HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(sigma11[j*triangles[j]->getBoundingPoints().size()+k]-sigma11_min)/(sigma11_max-sigma11_min), 1., 1.) ;
						glColor3f(c1, c2, c3) ;
						glVertex2f( double(triangles[j]->getBoundingPoint(k).x + vx) ,  double(triangles[j]->getBoundingPoint(k).y + vy) );
						
					}
					glEnd() ;
				}
			}
		glEndList() ;
		
		double vonMises_max = vonMises.max() ;
		double vonMises_min = vonMises.min() ;
		
		glNewList(  DISPLAY_LIST_VON_MISES,  GL_COMPILE ) ;
		for (unsigned int j=0 ; j< triangles.size() ; j++ )
		{
			
			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !cracked[j]&& !triangles[j]->getBehaviour()->fractured())
			{
				double c1 ;
				double c2 ;
				double c3 ;

				double vx = x[triangles[j]->first->id*2]; 
				double vy = x[triangles[j]->first->id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(vonMises[j*triangles[j]->getBoundingPoints().size()]-vonMises_min)/(vonMises_max-vonMises_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				
				glVertex2f(double(triangles[j]->getBoundingPoint(0).x + vx) , double(triangles[j]->getBoundingPoint(0).y + vy) );
				
				for(size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++)
				{
					vx = x[triangles[j]->getBoundingPoint(k).id*2];
					vy = x[triangles[j]->getBoundingPoint(k).id*2+1]; 
					
					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(vonMises[j*triangles[j]->getBoundingPoints().size()+k]-vonMises_min)/(vonMises_max-vonMises_min), 1., 1.) ;
					glColor3f(c1, c2, c3) ;
					glVertex2f( double(triangles[j]->getBoundingPoint(k).x + vx) ,  double(triangles[j]->getBoundingPoint(k).y + vy) );
					
				}
				glEnd() ;
			}
		}
		glEndList() ;
		
		
		double angle_max = angle.max() ;
		double angle_min = angle.min() ;
		glNewList(  DISPLAY_LIST_ANGLE,  GL_COMPILE ) ;
		for (unsigned int j=0 ; j< triangles.size() ; j++ )
		{
			
			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !cracked[j]&& !triangles[j]->getBehaviour()->fractured())
			{
				double c1 ;
				double c2 ;
				double c3 ;

				double vx = x[triangles[j]->first->id*2]; 
				double vy = x[triangles[j]->first->id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(angle[j*triangles[j]->getBoundingPoints().size()]-angle_min)/(angle_max-angle_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				
				glVertex2f(double(triangles[j]->getBoundingPoint(0).x + vx) , double(triangles[j]->getBoundingPoint(0).y + vy) );
				
				for(size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++)
				{
					vx = x[triangles[j]->getBoundingPoint(k).id*2];
					vy = x[triangles[j]->getBoundingPoint(k).id*2+1]; 
					
					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(angle[j*triangles[j]->getBoundingPoints().size()+k]-angle_min)/(angle_max-angle_min), 1., 1.) ;
					glColor3f(c1, c2, c3) ;
					glVertex2f( double(triangles[j]->getBoundingPoint(k).x + vx) ,  double(triangles[j]->getBoundingPoint(k).y + vy) );
					
				}
				glEnd() ;
			}
		}
		glEndList() ;
		
		double sigma22_min = sigma22.min() ;
		double sigma22_max = sigma22.max() ;
		
		glNewList(  DISPLAY_LIST_STRAIN_YY,  GL_COMPILE ) ;
		for (unsigned int j=0 ; j< triangles.size() ; j++ )
		{
			
			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !cracked[j]&& !triangles[j]->getBehaviour()->fractured())
			{
				double c1 ;
				double c2 ;
				double c3 ;
				
				double vx = x[triangles[j]->first->id*2]; 
				double vy = x[triangles[j]->first->id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(sigma22[j*triangles[j]->getBoundingPoints().size()]-sigma22_min)/(sigma22_max-sigma22_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				glVertex2f(double(triangles[j]->getBoundingPoint(0).x + vx) , double(triangles[j]->getBoundingPoint(0).y + vy) );
				
				for(size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++)
				{
					vx = x[triangles[j]->getBoundingPoint(k).id*2];
					vy = x[triangles[j]->getBoundingPoint(k).id*2+1]; 
					
					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(sigma22[j*triangles[j]->getBoundingPoints().size()+k]-sigma22_min)/(sigma22_max-sigma22_min), 1., 1.) ;
					glColor3f(c1, c2, c3) ;
					glVertex2f( double(triangles[j]->getBoundingPoint(k).x + vx) ,  double(triangles[j]->getBoundingPoint(k).y + vy) );
					
				}
				glEnd() ;
			}
		}
		glEndList() ;
		
		double sigma12_min = sigma12.min() ;
		double sigma12_max = sigma12.max() ;
		glNewList(  DISPLAY_LIST_STRAIN_XY,  GL_COMPILE ) ;
		
		for (unsigned int j=0 ; j< triangles.size() ; j++ )
		{
			
			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !cracked[j]&& !triangles[j]->getBehaviour()->fractured())
			{
				double c1 ;
				double c2 ;
				double c3 ;
				
				double vx = x[triangles[j]->first->id*2]; 
				double vy = x[triangles[j]->first->id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(sigma12[j*triangles[j]->getBoundingPoints().size()]-sigma12_min)/(sigma12_max-sigma12_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				glVertex2f(double(triangles[j]->getBoundingPoint(0).x + vx) , double(triangles[j]->getBoundingPoint(0).y + vy) );
				
				for(size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++)
				{
					vx = x[triangles[j]->getBoundingPoint(k).id*2];
					vy = x[triangles[j]->getBoundingPoint(k).id*2+1]; 
					
					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(sigma12[j*triangles[j]->getBoundingPoints().size()+k]-sigma12_min)/(sigma12_max-sigma12_min), 1., 1.) ;
					glColor3f(c1, c2, c3) ;
					glVertex2f( double(triangles[j]->getBoundingPoint(k).x + vx) ,  double(triangles[j]->getBoundingPoint(k).y + vy) );
					
				}
				glEnd() ;
			}
		}
		glEndList() ;
		
		glNewList(  DISPLAY_LIST_STIFFNESS,  GL_COMPILE ) ;
		for (unsigned int j=0 ; j< triangles.size() ; j++ )
		{
			
			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !cracked[j]&& !triangles[j]->getBehaviour()->fractured())
			{
				double c1 ;
				double c2 ;
				double c3 ;
				
				double vx = x[triangles[j]->first->id*2]; 
				double vy = x[triangles[j]->first->id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				
				Point a = triangles[j]->inLocalCoordinates(triangles[j]->getBoundingPoint(0)) ;
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(triangles[j]->getBehaviour()->getTensor(a)[0][0]-E_min)/(E_max-E_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				glVertex2f(double(triangles[j]->getBoundingPoint(0).x + vx) , double(triangles[j]->getBoundingPoint(0).y + vy) );
				
				for(size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++)
				{
					vx = x[triangles[j]->getBoundingPoint(k).id*2];
					vy = x[triangles[j]->getBoundingPoint(k).id*2+1]; 
					a = triangles[j]->inLocalCoordinates(triangles[j]->getBoundingPoint(k)) ;
					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(triangles[j]->getBehaviour()->getTensor(a)[0][0]-E_min)/(E_max-E_min), 1., 1.) ;
					glColor3f(c1, c2, c3) ;
					glVertex2f( double(triangles[j]->getBoundingPoint(k).x + vx) ,  double(triangles[j]->getBoundingPoint(k).y + vy) );
					
				}
				glEnd() ;
			}
		}
		glEndList() ;
		
		double epsilon11_min = epsilon11.min() ;
		double epsilon11_max = epsilon11.max() ;
		glNewList(  DISPLAY_LIST_STRESS_XX,  GL_COMPILE ) ;
			
		for (unsigned int j=0 ; j< triangles.size() ; j++ )
		{
			
			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !cracked[j]&& !triangles[j]->getBehaviour()->fractured())
			{
				double c1 ;
				double c2 ;
				double c3 ;
				
				double vx = x[triangles[j]->first->id*2]; 
				double vy = x[triangles[j]->first->id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(epsilon11[j*triangles[j]->getBoundingPoints().size()]-epsilon11_min)/(epsilon11_max-epsilon11_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				glVertex2f(double(triangles[j]->getBoundingPoint(0).x + vx) , double(triangles[j]->getBoundingPoint(0).y + vy) );
				
				for(size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++)
				{
					vx = x[triangles[j]->getBoundingPoint(k).id*2];
					vy = x[triangles[j]->getBoundingPoint(k).id*2+1]; 
					
					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(epsilon11[j*triangles[j]->getBoundingPoints().size()+k]-epsilon11_min)/(epsilon11_max-epsilon11_min), 1., 1.) ;
					glColor3f(c1, c2, c3) ;
					glVertex2f( double(triangles[j]->getBoundingPoint(k).x + vx) ,  double(triangles[j]->getBoundingPoint(k).y + vy) );
					
				}
				glEnd() ;
			}
		}
		glEndList() ;
		
				
		double epsilon22_min = epsilon22.min() ;
		double epsilon22_max = epsilon22.max() ;
		
		glNewList(  DISPLAY_LIST_STRESS_YY,  GL_COMPILE ) ;
		for (unsigned int j=0 ; j< triangles.size() ; j++ )
		{
			
			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !cracked[j]&& !triangles[j]->getBehaviour()->fractured())
			{
				double c1 ;
				double c2 ;
				double c3 ;
				
				double vx = x[triangles[j]->first->id*2]; 
				double vy = x[triangles[j]->first->id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(epsilon22[j*triangles[j]->getBoundingPoints().size()]-epsilon22_min)/(epsilon22_max-epsilon22_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				glVertex2f(double(triangles[j]->getBoundingPoint(0).x + vx) , double(triangles[j]->getBoundingPoint(0).y + vy) );
				
				for(size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++)
				{
					vx = x[triangles[j]->getBoundingPoint(k).id*2];
					vy = x[triangles[j]->getBoundingPoint(k).id*2+1]; 
					
					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(epsilon22[j*triangles[j]->getBoundingPoints().size()+k]-epsilon22_min)/(epsilon22_max-epsilon22_min), 1., 1.) ;
					glColor3f(c1, c2, c3) ;
					glVertex2f( double(triangles[j]->getBoundingPoint(k).x + vx) ,  double(triangles[j]->getBoundingPoint(k).y + vy) );
					
				}
				glEnd() ;
			}
		}
		
		glEndList() ;
		
		double epsilon12_min = epsilon12.min() ;
		double epsilon12_max = epsilon12.max() ;
		
		glNewList(  DISPLAY_LIST_STRESS_XY,  GL_COMPILE ) ;
		for (unsigned int j=0 ; j< triangles.size() ; j++ )
		{
			
			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !cracked[j]&& !triangles[j]->getBehaviour()->fractured())
			{
				double c1 ;
				double c2 ;
				double c3 ;
				
				double vx = x[triangles[j]->first->id*2]; 
				double vy = x[triangles[j]->first->id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(epsilon12[j*triangles[j]->getBoundingPoints().size()]-epsilon12_min)/(epsilon12_max-epsilon12_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				glVertex2f(double(triangles[j]->getBoundingPoint(0).x + vx) , double(triangles[j]->getBoundingPoint(0).y + vy) );
				
				for(size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++)
				{
					vx = x[triangles[j]->getBoundingPoint(k).id*2];
					vy = x[triangles[j]->getBoundingPoint(k).id*2+1]; 
					
					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(epsilon12[j*triangles[j]->getBoundingPoints().size()+k]-epsilon12_min)/(epsilon12_max-epsilon12_min), 1., 1.) ;
					glColor3f(c1, c2, c3) ;
					glVertex2f( double(triangles[j]->getBoundingPoint(k).x + vx) ,  double(triangles[j]->getBoundingPoint(k).y + vy) );
					
				}
				glEnd() ;
			}
		}
		glEndList() ;
		
		glNewList(  DISPLAY_LIST_ENRICHMENT,  GL_COMPILE ) ;
		glBegin(GL_TRIANGLES);
		for (unsigned int j=0 ; j< triangles.size() ; j++ )
		{

			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR )
			{
				double c1 ;
				double c2 ;
				double c3 ;
				
				double enrichment = triangles[j]->getEnrichmentFunctions().size() ;
				//HSVtoRGB( &c1, &c2, &c3, 180. + 180.*(sigma12[j]-sigma12.min())/(sigma12.max()-sigma12.min()), 1., 1. ) 
				
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*enrichment/20., 1., 1.) ;
				if(enrichment)
					glColor3f(c1, c2, c3) ;
				else
					glColor3f(.25, .25, .25) ;
				
				double vx = x[triangles[j]->first->id*2]; 
				double vy = x[triangles[j]->first->id*2+1]; 
				
				glVertex2f( double(triangles[j]->first->x + vx) ,
				            double(triangles[j]->first->y + vy) );
				
				vx = x[triangles[j]->second->id*2];
				vy = x[triangles[j]->second->id*2+1]; 

				glVertex2f( double(triangles[j]->second->x + vx) ,
				            double(triangles[j]->second->y + vy) );
				
				
				vx = x[triangles[j]->third->id*2]; 
				vy = x[triangles[j]->third->id*2+1]; 
				

				glVertex2f( double(triangles[j]->third->x + vx) ,
				            double(triangles[j]->third->y + vy) );
			}
		}
		glEnd();
		glEndList() ;
		
		
		glNewList(  DISPLAY_LIST_ELEMENTS,  GL_COMPILE ) ;
		glColor3f(1, 1, 1) ;
		for(unsigned int j=0 ; j< triangles.size() ; j++ )
		{
			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR)
			{
				if(triangles[j]->getBehaviour()->fractured())
					glColor3f(1, 0, 0) ;
				else
					glColor3f(1, 1, 1) ;
				
				glBegin(GL_LINE_LOOP);
				for(size_t k = 0 ; k < triangles[j]->getBoundingPoints().size() ; k++)
				{
					double vx = x[triangles[j]->getBoundingPoint(k).id*2]; 
					double vy = x[triangles[j]->getBoundingPoint(k).id*2+1]; 
				
					glVertex2f( double(triangles[j]->getBoundingPoint(k).x+vx) ,  double(triangles[j]->getBoundingPoint(k).y+vy) );
					
				}
				glEnd();
			}
			
			glColor3f(1, 1, 1) ;
		}
		glEndList() ;
		
		
		glNewList(  DISPLAY_LIST_CRACK,  GL_COMPILE ) ;
		glLineWidth(4) ;
		for(size_t k  = 0 ; k < crack.size() ; k++)
		{
			glColor3f(1, 0, 0) ;
// 			for(unsigned int j=0 ; j< tris__.size() ; j++ )
// 			{
// 				glBegin(GL_LINE_LOOP);
// 				double vx = x[tris__[j]->first->id*2]; 
// 				double vy = x[tris__[j]->first->id*2+1]; 
// 				
// 				glVertex2f( double(tris__[j]->first->x/*+ vx*/) ,
// 				            double(tris__[j]->first->y/*+ vy*/) );
// 				
// 				vx = x[tris__[j]->second->id*2]; 
// 				vy = x[tris__[j]->second->id*2+1]; 
// 				
// 				glVertex2f( double(tris__[j]->second->x/*+ vx*/) ,
// 				            double(tris__[j]->second->y/*+ vy*/) );
// 				
// 				vx = x[tris__[j]->third->id*2]; 
// 				vy = x[tris__[j]->third->id*2+1]; 
// 				
// 				glVertex2f( double(tris__[j]->third->x/*+ vx*/) ,
// 				            double(tris__[j]->third->y/*+ vy*/) );
// 				glEnd();
// 			}
// 			
// 			glColor3f(0, 1, 1) ;
			glBegin(GL_LINES) ;
			for(size_t j=0 ; j< crack[k]->getBoundingPoints().size()-1 ; j++ )
			{
				glVertex2f( double(crack[k]->getBoundingPoint(j).x) ,
				            double(crack[k]->getBoundingPoint(j).y) );
				glVertex2f( double(crack[k]->getBoundingPoint(j+1).x) ,
				            double(crack[k]->getBoundingPoint(j+1).y) );
			}
			glEnd();
		}
		
// 		for(unsigned int j=0 ; j< triangles.size() ; j++ )
// 		{
// 			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR)
// 			{
// 				
// 				
// 				Vector t = triangles[j]->getState()->getPrincipalStresses(triangles[j]->getCenter()) ;
// 				glBegin(GL_LINE_LOOP);
// 				
// 				glColor3f(1, 1, 1) ;
// 				glVertex2f( triangles[j]->getCenter().x ,  triangles[j]->getCenter().y  );
// 				glColor3f(1, 1, 1) ;
// 				glVertex2f( triangles[j]->getCenter().x +5.*t[0],  triangles[j]->getCenter().y +5.*t[1] );
// 				
// 				glEnd();
// 			}
// 			
// 			glColor3f(1, 1, 1) ;
// 		}
		glLineWidth(1) ;
		glEndList() ;
		
		dlist = true ;
		glCallList(current_list) ;
	}
	else
	{
		//glCallList(DISPLAY_LIST_DISPLACEMENT) ;
		//glCallList(DISPLAY_LIST_STRAIN) ;
		double c1, c2, c3 = 0;
		HSVtoRGB( &c1, &c2, &c3, 180. + 0, 1., 1.) ;
// 		glBegin(GL_LINE) ;
// 		glVertex2f(3.5 ,
// 		           3. );
// 		glVertex2f(3.5 ,
// 		           -3. );
// 		glEnd() ;
		
		glCallList(current_list) ;
		glCallList(DISPLAY_LIST_CRACK) ;
// 		if(current_list == DISPLAY_LIST_ELEMENTS)
// 			glCallList(DISPLAY_LIST_CRACK) ;
		if(current_list == DISPLAY_LIST_ENRICHMENT)
		{
			glCallList(DISPLAY_LIST_ELEMENTS) ;
			
		}
		
		glColor3f(1, 1, 1) ;

		
	}
	glColor3f(1, 0, 0) ;
	glFlush();
	glutSwapBuffers();
}

int main(int argc, char *argv[])
{
	
	
	Matrix m0(3,3) ;
	m0[0][0] = E/(1-nu*nu) ; m0[0][1] =E/(1-nu*nu)*nu ; m0[0][2] = 0 ;
	m0[1][0] = E/(1-nu*nu)*nu ; m0[1][1] = E/(1-nu*nu) ; m0[1][2] = 0 ; 
	m0[2][0] = 0 ; m0[2][1] = 0 ; m0[2][2] = E/(1-nu*nu)*(1.-nu)/2. ; 
	
	Sample sample(NULL, 8,4,0,0) ;
	
	Sample reinforcement0(NULL, 8,.15,0,.5) ;
	reinforcement0.setBehaviour(new Stiffness(m0*5)) ;
	
	Sample reinforcement1(NULL, 8,.15,0,-.5) ;
	reinforcement1.setBehaviour(new Stiffness(m0*5)) ;
	
	FeatureTree F(&sample) ;
	featureTree = &F ;
// 	F.addFeature(&sample,&reinforcement0) ;
// 	F.addFeature(&sample,&reinforcement1) ;
	Point A(0,-.7); Point b(-0.1,-1.5); Point c(0.1,-1.5) ;
// 	TriangularPore * pore = new TriangularPore(A, b, c) ;
// 	F.addFeature(&sample,pore) ;
	
	
	PointSet ptset(5) ;
	
	std::valarray<Point *> centerpoint(2) ;
	
	centerpoint[0] = new Point(0, -1.5 ) ;
	centerpoint[1] = new Point(0, -.7) ;
	
	std::valarray<Point *> centerpoint2(2) ;
	
	centerpoint2[0] = new Point(5,1.5 ) ;
	centerpoint2[1] = new Point(3.5, 1.5) ;
	
	
	std::valarray<Point *> side0(2) ;
	
	side0[1] = new Point(-0.7,-4 ) ;
	side0[0] = new Point(-0.7, 4 ) ;
	
	std::valarray<Point *> side1(2) ;
	
	side1[1] = new Point(-4,0.75 ) ;
	side1[0] = new Point(4,0.75 ) ;
	
	std::valarray<Point *> side2(2) ;
	
	side2[1] = new Point(0.7,4 ) ;
	side2[0] = new Point(0.7,-4 ) ;
	
// 	PointSet side3(2) ;
// 	
// 	side3.set(1, new Point(1.5,-0.75 )) ;
// 	side3.set(0, new Point(-1.5,-0.75 )) ;
	
// 	crack = generateCracks(10) ;
// 	for(size_t j = 0 ; j < crack.size() ; j++)
// 		F.addFeature(&sample, crack[j]) ;
	
// 	for(size_t j = 0 ; j < crack.size() ; j++)
// 	{
// 		crack[j]->getHead()->print() ; std::cout << " <- " << j << std::endl ;
// 	}
// 	Crack cr(&sample, centerpoint, 0.03) ;
// 	crack.push_back(&cr) ;
// 	F.addFeature(&sample, crack[0]) ;
// 	
// 	Crack cr2(&sample, centerpoint2, 0.03) ;
// 	crack.push_back(&cr2) ;
// 	F.addFeature(&sample, crack[1]) ;

// 	crack.push_back(new Crack(&sample, &side0, 0.1)) ;
// 	F.addFeature(sample, crack[0]) ;
// 	crack.push_back(new Crack(&sample, &side1, 0.1)) ;
// 	F.addFeature(sample, crack[1]) ;
// 	crack.push_back(new Crack(&sample, &side2, 0.1)) ;
// 	F.addFeature(sample, crack[2]) ;
// 	crack.push_back(new Crack(&sample, &side3, 0.1)) ;
// 	F.addFeature(sample, crack[3]) ;
	
	i_et_p = generateInclusionsAndPores(1024, .05, &m0, &sample, &F) ;
// 	Inclusion * inc = new Inclusion(1, 0,0) ;
// 	F.addFeature(&sample,inc) ;
// 	inc->setBehaviour(new Stiffness(m0)) ;
// 	F.addFeature(&sample,new TriangularPore(Point(-1.25, -4), Point(-1,-2.5), Point(-0.75,-4))) ;
// 	F.addFeature(&sample,new TriangularPore(Point(2.25, -4), Point(2,-1.5), Point(1.75,-4))) ;
// 	F.addFeature(&sample,new Pore(1, 1.5,1)) ;
// 	Inclusion * inc = new Inclusion(1, 0,0) ;
// 	F.addFeature(&sample,inc) ;
	
	
	Circle cercle(1, 0,0) ;
	
	
// 	sample.setBehaviour(new BimaterialInterface(&cercle, m0,  m0*4)) ;
	Vector a(double(0), 3) ;
	a[0] = 0.001 ;
	a[1] = 0.00 ;
	a[1] = 0.00 ;
// 	inc->setBehaviour(new Stiffness(m0*4)) ;
	sample.setBehaviour(new WeibullDistributedStiffness(m0*0.125, 0.005)) ;
	generateExpansiveZones(10, i_et_p.first, F) ;
// 	sample.setBehaviour(new Stiffness(m0*0.35)) ;
// 	sample.setBehaviour(new StiffnessAndFracture(m0, 0.03)) ;
// 	F.addFeature(&sample,new EnrichmentInclusion(1, 0,0)) ;
// 	F.addFeature(&sample,new Pore(1, 0,0)) ;
// 	F.addFeature(&sample,new Pore(0.75, 1,-1)) ;
// 	F.addFeature(&sample,new Pore(0.75, -1,-1)) ;
// 	F.addFeature(&sample,new Pore(0.75, -1,1)) ;
	
	F.sample(1200) ;
	F.setOrder(LINEAR) ;

	F.generateElements() ;
	
	for(size_t j = 0 ; j < crack.size() ; j++)
		crack[j]->setInfluenceRadius(0.03) ;
// 	
	step() ;
	
	glutInit(&argc, argv) ;	
	glutInitDisplayMode(GLUT_RGBA) ;
	glutInitWindowSize(600, 600) ;
	glutReshapeFunc(reshape) ;
	glutCreateWindow("coucou !") ;
	
	int submenu = glutCreateMenu(Menu) ;
	
	glutAddMenuEntry(" Displacements ", ID_DISP);
	glutAddMenuEntry(" Strain (s) xx ", ID_STRAIN_XX);
	glutAddMenuEntry(" Strain (s) yy ", ID_STRAIN_YY);
	glutAddMenuEntry(" Strain (s) xy ", ID_STRAIN_XY);
	glutAddMenuEntry(" Stress (e) xx ", ID_STRESS_XX);
	glutAddMenuEntry(" Stress (e) yy ", ID_STRESS_YY);
	glutAddMenuEntry(" Stress (e) xy ", ID_STRESS_XY);
	glutAddMenuEntry(" Elements      ", ID_ELEM);
	glutAddMenuEntry(" Stiffness     ", ID_STIFNESS);
	glutAddMenuEntry(" Von Mises     ", ID_VON_MISES);
	glutAddMenuEntry(" Princ. angle  ", ID_ANGLE);
	glutAddMenuEntry(" Enrichment    ", ID_ENRICHMENT);
	
	glutCreateMenu(Menu) ;

 	glutAddMenuEntry(" Step          ", ID_NEXT);
	glutAddMenuEntry(" Step time     ", ID_NEXT_TIME);
	glutAddMenuEntry(" Zoom in       ", ID_ZOOM);
	glutAddMenuEntry(" Zoom out      ", ID_UNZOOM);
	glutAddMenuEntry(" Amplify       ", ID_AMPLIFY);
	glutAddMenuEntry(" Deamplify     ", ID_DEAMPLIFY);
	glutAddSubMenu(  " Display       ", submenu);
	glutAddMenuEntry(" Quit          ", ID_QUIT) ;
	
	
	glutAttachMenu(GLUT_RIGHT_BUTTON) ;
	
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glShadeModel(GL_SMOOTH);
	
	glutDisplayFunc(Display) ;
	glutMainLoop() ;
	
// 	delete dt ;
	
	return 0 ;
}
