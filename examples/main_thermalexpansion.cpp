// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../utilities/samplingcriterion.h"
#include "../features/features.h"
#include "../physics/physics.h"
#include "../physics/mohrcoulomb.h"
#include "../physics/ruptureenergy.h"
#include "../physics/kelvinvoight.h"
#include "../physics/vonmises.h"
#include "../physics/spatially_distributed_stiffness.h"
#include "../features/pore.h"
#include "../features/sample.h"
#include "../features/inclusion.h"
#include "../features/expansiveZone.h"
#include "../features/crack.h"
#include "../features/enrichmentInclusion.h"
#include "../features/expansiveZone.h"
#include "../mesher/delaunay_3d.h"
#include "../solvers/assembly.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"

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
#define ID_FRAC_CRIT 23

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
#define DISPLAY_LIST_FRAC_CRIT 25
#include "../features/sample3d.h"

using namespace Mu ;
using namespace std;

FeatureTree * featureTree ;
std::vector<DelaunayTriangle *> triangles ;
DelaunayTree *dt ; //(pts) ;
std::vector<bool> cracked ;
std::vector<BranchedCrack *> crack ;

double E_min = 10;
double E_max = 0;

double x_max = 0 ;
double y_max = 0 ;

double x_min = 0 ;
double y_min = 0 ;

double timepos = 0.001 ;

bool firstRun = true ;

std::vector<DelaunayTriangle *> tris__ ;

std::pair<std::vector<Inclusion * >, std::vector<Pore * > > i_et_p ;

std::vector<std::pair<ExpansiveZone *, Inclusion *> > zones ;

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
Vector fracCrit(0) ;

Vector g_count(0) ;

double width = 1. ;
double height = 1. ;
Sample sample(NULL, width, height, 0.5, 0.5) ;	

//double temperature = 50 ;

size_t current_list = DISPLAY_LIST_STRAIN_XX ;
double factor = 200 ;
MinimumAngle cri(M_PI/6.) ;
bool nothingToAdd = false ;
bool dlist = false ;
int count = 0 ;
double aggregateArea = 0;

std::string outfilename("results") ;
std::fstream outresultfile ;


int totit = 5 ;

std::vector<double> energy ;

std::string itoa(int value, int base) {

	enum { kMaxDigits = 35 };
	std::string buf;
	buf.reserve( kMaxDigits ); // Pre-allocate enough space.

	// check that the base if valid
	if (base < 2 || base > 16) return buf;

	int quotient = value;
	
	
	// Translating number to string with base:
	do {
		buf += "0123456789abcdef"[ std::abs( quotient % base ) ];
		quotient /= base;
	} while ( quotient );
	
	// Append the negative sign for base 10
	if ( value < 0 && base == 10) buf += '-';
	
	std::reverse( buf.begin(), buf.end() );
	
	return buf;
	
}




void setBC()
{

}


Vector step(int nInc)
{
	
  
	double avg_e_xx = 0;
	double avg_e_yy = 0;
	double avg_e_xy = 0;
	double avg_s_xx = 0;
	double avg_s_yy = 0;
	double avg_s_xy = 0;
	
	bool cracks_did_not_touch = true;
	size_t max_growth_steps = 1;
	size_t countit = 0;	
	
	while ( (cracks_did_not_touch) && (countit < max_growth_steps) )
	{
	  
		countit++;
		std::cout << "\r iteration " << countit << "/" << max_growth_steps << std::flush ;
		setBC() ;
      
		int limit = 0 ;
		while(!featureTree->step(timepos) && limit < 50)//as long as we can update the features
		{
			std::cout << "." << std::flush ;
// 			timepos-= 0.0001 ;
			setBC() ;
			limit++ ;
	  // check if the two cracks did not touch
			cracks_did_not_touch = true;
			for(size_t j = 0 ; j < crack.size() ; j++)
			{
				for(size_t k = j+1 ; k < crack.size() ; k++)
				{
		  
					if (static_cast<SegmentedLine *>(crack[j])->intersects(static_cast<SegmentedLine *>(crack[k])))
					{
						cracks_did_not_touch = false;
						break;
					}
		  //Circle headj(radius, crack[j]->getHead());
		  
		  //	      head0.intersects(crack[1]);
				}
			}
		}
	
  // Prints the crack geo to a file for each crack
		if (cracks_did_not_touch == false) // if cracks touched
		{
			std::cout << "** Cracks touched exporting file **" << endl;	
				// Print the state of the cracks to a file  
				std::string filename = "crackGeo.txt";
				fstream filestr;
				filestr.open (filename.c_str(), fstream::in | fstream::out | fstream::app);
				filestr << "Crack vertices" << std::endl;
				filestr << "x" << " " << "y" << std::endl ; 
		}
  

		timepos+= 0.0001 ;

		triangles = featureTree->getTriangles() ;
		x.resize(featureTree->getDisplacements().size()) ;
		x = featureTree->getDisplacements() ;

		dt = featureTree->getDelaunayTree() ;
		sigma.resize(triangles.size()*triangles[0]->getBoundingPoints().size()*3) ;
		epsilon.resize(triangles.size()*triangles[0]->getBoundingPoints().size()*3) ;
	
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
		fracCrit.resize(sigma.size()/3) ;
		g_count.resize(sigma.size()/3) ;
		std::cout << "unknowns :" << x.size() << std::endl ;
		
		if(crack.size() > 0)
			tris__ = crack[0]->getTriangles(dt) ;
		
		for(size_t k = 1 ; k < crack.size() ; k++)
		{
			std::vector<DelaunayTriangle *> temp = crack[k]->getTriangles(dt) ;
			if(tris__.empty())
				tris__ = temp ;
			else if(!temp.empty())
				tris__.insert(tris__.end(), temp.begin(), temp.end() ) ;
		}
		cracked.clear() ;
		
		int npoints = triangles[0]->getBoundingPoints().size() ;
	
		double area = 0 ;
		avg_e_xx = 0;
		avg_e_yy = 0;
		avg_e_xy = 0;
		avg_s_xx = 0;
		avg_s_yy = 0;
		avg_s_xy = 0;
		double e_xx_max = 0 ;
		double e_xx_min = 0 ;
		double ex_count = 0 ;
		double enr = 0 ;
		for(size_t k = 0 ; k < triangles.size() ; k++)
		{
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
		
		
		
			if(!in && !triangles[k]->getBehaviour()->fractured() && triangles[k]->getBehaviour()->type != VOID_BEHAVIOUR)
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
					if(triangles[k]->getBoundingPoint(p).x > sample.width()*.9999)
					{
						if(e_xx_max < x[triangles[k]->getBoundingPoint(p).id*2])
							e_xx_max=x[triangles[k]->getBoundingPoint(p).id*2] ;
// 						ex_count++ ;
					}
					if(triangles[k]->getBoundingPoint(p).x < sample.width()*.0001)
					{
						if(e_xx_min > x[triangles[k]->getBoundingPoint(p).id*2])
							e_xx_min=x[triangles[k]->getBoundingPoint(p).id*2] ;
// 						ex_count++ ;
					}
				}
				area += triangles[k]->area() ;
				if(triangles[k]->getBehaviour()->type != VOID_BEHAVIOUR)
				{
					if(triangles[k]->getBehaviour()->param[0][0] > E_max)
						E_max = triangles[k]->getBehaviour()->param[0][0] ;
					if(triangles[k]->getBehaviour()->param[0][0] < E_min)
						E_min = triangles[k]->getBehaviour()->param[0][0] ;
				}
				
				g_count[k*npoints]++ ;
				g_count[k*npoints+1]++ ;
				g_count[k*npoints+2]++ ;
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
					g_count[k*npoints+3]++ ;
					g_count[k*npoints+4]++ ;
					g_count[k*npoints+5]++ ;
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
				
				for(size_t l = 0 ; l < npoints ; l++)
				{
					Vector vm0 = triangles[k]->getState().getPrincipalStresses(triangles[k]->getBoundingPoint(l)) ;
					vonMises[k*npoints+l]  += sqrt(((vm0[0]-vm0[1])*(vm0[0]-vm0[1]))/2.) ;
	
					double agl = triangles[k]->getState().getPrincipalAngle(triangles[k]->getBoundingPoint(l)) ;
					angle[k*npoints+l]  += agl ;
					if(triangles[k]->getBehaviour()->getFractureCriterion())
					{
						fracCrit[k*npoints+l] = triangles[k]->getBehaviour()->getFractureCriterion()->grade(triangles[k]->getState()) ;
						enr += triangles[k]->getState().elasticEnergy() ;
					}
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
	
			}
			else
			{
				sigma11[k*npoints] += 0 ;
				sigma22[k*npoints] += 0 ;
				sigma12[k*npoints] += 0 ;
				sigma11[k*npoints+1] += 0 ;
				sigma22[k*npoints+1] += 0 ;
				sigma12[k*npoints+1] += 0 ;
				sigma11[k*npoints+2] += 0 ;
				sigma22[k*npoints+2] += 0 ;
				sigma12[k*npoints+2] += 0 ;
				
				if(npoints >3)
				{
					sigma11[k*npoints+3] += 0 ;
					sigma22[k*npoints+3] += 0 ;
					sigma12[k*npoints+3] += 0 ;
					sigma11[k*npoints+4] += 0 ;
					sigma22[k*npoints+4] += 0 ;
					sigma12[k*npoints+4] += 0 ;
					sigma11[k*npoints+5] += 0 ;
					sigma22[k*npoints+5] += 0 ;
					sigma12[k*npoints+5] += 0 ;
				}
				
				epsilon11[k*npoints] += 0 ;
				epsilon22[k*npoints] += 0 ;
				epsilon12[k*npoints] += 0 ;
				epsilon11[k*npoints+1] += 0 ;
				epsilon22[k*npoints+1] += 0 ;
				epsilon12[k*npoints+1] += 0 ;
				epsilon11[k*npoints+2] += 0 ;
				epsilon22[k*npoints+2] += 0 ;
				epsilon12[k*npoints+2] += 0 ;
				
				if(npoints > 3)
				{
					epsilon11[k*npoints+3] += 0 ;
					epsilon22[k*npoints+3] += 0 ;
					epsilon12[k*npoints+3] += 0 ;
					epsilon11[k*npoints+4] += 0 ;
					epsilon22[k*npoints+4] += 0 ;
					epsilon12[k*npoints+4] += 0 ;
					epsilon11[k*npoints+5] += 0 ;
					epsilon22[k*npoints+5] += 0 ;
					epsilon12[k*npoints+5] += 0 ;
				}  
				
				for(size_t l = 0 ; l < triangles[k]->getBoundingPoints().size() ; l++)
				{
					vonMises[k*triangles[k]->getBoundingPoints().size()+l]  += 0 ;
					angle[k*triangles[k]->getBoundingPoints().size()+l]  += 0 ;
				}
			}
		}
		
		std::string filename("triangles") ;
		filename.append(itoa(totit++, 10)) ;
		std::cout << filename << std::endl ;
		std::fstream outfile  ;
		outfile.open(filename.c_str(), std::ios::out) ;
		
		outfile << "TRIANGLES" << std::endl ;
		outfile << triangles.size() << std::endl ;
		outfile << 3 << std::endl ;
		outfile << 8 << std::endl ;
		
		for(size_t j = 0 ; j < triangles.size() ;j++)
		{
			for(size_t l = 0 ; l < triangles[j]->getBoundingPoints().size() ; l++)
			{
				outfile << triangles[j]->getBoundingPoint(l).x << " " << triangles[j]->getBoundingPoint(l).y << " ";
			}

			for(size_t l = 0 ; l < triangles[j]->getBoundingPoints().size() ; l++)
			{
				outfile <<  epsilon11[j*3+l] << " ";
			}
			for(size_t l = 0 ; l < triangles[j]->getBoundingPoints().size() ; l++)
			{
				outfile <<  epsilon22[j*3+l] << " " ;
			}
			for(size_t l = 0 ; l < triangles[j]->getBoundingPoints().size() ; l++)
			{
				outfile <<   epsilon12[j*3+l]<< " " ;
			}
			for(size_t l = 0 ; l < triangles[j]->getBoundingPoints().size() ; l++)
			{
				outfile <<  sigma11[j*3+l]<< " " ;
			}
			for(size_t l = 0 ; l < triangles[j]->getBoundingPoints().size() ; l++)
			{
				outfile <<  sigma22[j*3+l]<< " ";
			}
			for(size_t l = 0 ; l < triangles[j]->getBoundingPoints().size() ; l++)
			{
				outfile <<  sigma12[j*3+l] << " ";
			}
			for(size_t l = 0 ; l < triangles[j]->getBoundingPoints().size() ; l++)
			{
				outfile << vonMises[j*3+l]<< " " ;
			}
			for(size_t l = 0 ; l < triangles[j]->getBoundingPoints().size() ; l++)
			{
				outfile <<  triangles[j]->getBehaviour()->getTensor(Point(.3, .3))[0][0] << " ";
			}
			outfile << "\n" ;
		}
		
		std::cout << std::endl ;
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
			
		std::cout << "energy index :" << enr << std::endl ;
		energy.push_back(enr) ;
		
		if(limit < 2)
			break ;
		
	}
	for(size_t i = 0 ; i < energy.size() ; i++)
		std::cout << energy[i] << std::endl ;
	
/*	std::vector<Point *> c ;
	Inclusion * inc_ ;
	double radius = featureTree->getFeature(1)->getRadius() ;
	for(int i=1 ; i < (nInc * nInc) + 1 ; i++)
	{
		inc_ = featureTree->getFeature(i) ;
		c = inc_->getBoundingPoints() ;
	}*/
	
	
	
	Vector results(6) ;
	results[0] = avg_e_xx ;
	results[1] = avg_e_xy ;
	results[2] = avg_e_yy ;
	results[3] = avg_s_xx ;
	results[4] = avg_s_xy ;
	results[5] = avg_s_yy ;
	return results ;
}

Matrix buildStiffnessMatrix(Vector p)
{
	Matrix prop(3,3) ;
	double E = p[1] * 1e9 ;
	double nu = p[2] ;
	prop[0][0] = E/(1.-nu*nu) ; prop[0][1] =E/(1.-nu*nu)*nu ; prop[0][2] = 0 ;
	prop[1][0] = E/(1.-nu*nu)*nu ; prop[1][1] = E/(1.-nu*nu) ; prop[1][2] = 0 ; 
	prop[2][0] = 0 ; prop[2][1] = 0 ; prop[2][2] = E/(1.-nu*nu)*(1.-nu)/2. ; 
	return prop ;  
}

std::vector<Inclusion *> buildInclusion(int nInc, double pInc)
{
	std::vector<Inclusion *> inc ;
	double radius = std::sqrt(pInc / (nInc * nInc * 3.1415926)) ;
	std::cout << radius << std::endl ;
	double dist = (1 - nInc * radius * 2) / (nInc + 1) ;
	std::cout << dist << std::endl ;
	for(int i = 0 ; i < nInc ; i++)
	{
		for(int j = 0 ; j < nInc ; j++)
		{
			inc.push_back(new Inclusion(radius,(i +1) * dist + radius * (2 * i + 1),(j +1) * dist  + radius * (2 * j + 1))) ;
//			std::cout << (i +1) * dist + radius * (2 * i + 1) << ";" << (j +1) * dist  + radius * (2 * j + 1) << std::endl ;
		}
	}
	return inc ;
}


Vector thermalDeformation(Vector matrixProp, Vector inclusionProp, int nInc, double pInc, double temperature)
{
	Vector result(6) ;
	FeatureTree F(&sample) ;
	featureTree = &F ;
	
	Matrix mK = buildStiffnessMatrix(matrixProp) ;
	Matrix iK = buildStiffnessMatrix(inclusionProp) ;
	double mA = matrixProp[0] * 1e-6 ;
	double iA = inclusionProp[0] * 1e-6 ;
	
	Vector mDef(3) ;
	mDef[0] = mA * temperature ;
	mDef[1] = mA * temperature ;
	mDef[2] = 0 ;
	Vector iDef(3) ;
	iDef[0] = iA * temperature ;
	iDef[1] = iA * temperature ;
	iDef[2] = 0 ;
	
	StiffnessWithImposedDeformation * mTDef = new StiffnessWithImposedDeformation(mK,mDef) ;
	StiffnessWithImposedDeformation * iTDef = new StiffnessWithImposedDeformation(iK,iDef) ;
	
	std::vector<Inclusion *> inc = buildInclusion(nInc,pInc) ;
	
	sample.setBehaviour(mTDef) ;
	for(int i = 0 ; i < inc.size() ; i++)
	{
		  inc[i]->setBehaviour(iTDef) ;
		  F.addFeature(&sample, inc[i]) ;
	}
	F.sample(1024) ;
	F.setOrder(LINEAR) ;
	F.generateElements() ;
	result = step(nInc) ;
	return result ;
		
}


int main(int argc, char *argv[])
{
  
	outresultfile.open(outfilename.c_str(), std::ios::out) ;
	Vector exp(6) ;
	
	Vector paste_prop(3) ;
	Vector similipaste_prop(3) ;
	Vector quartz_prop(3) ;
	Vector limestone_prop(3) ;

	paste_prop[0] = 9.69 ;
	paste_prop[1] = 25 ;
	paste_prop[2] = 0.22 ;

	similipaste_prop[0] = 4.34 ;
	similipaste_prop[1] = 25 ;
	similipaste_prop[2] = 0.22 ;

	quartz_prop[0] = 11.38 ;
	quartz_prop[1] = 70 ;
	quartz_prop[2] = 0.2 ;

	limestone_prop[0] = 4.34 ;
	limestone_prop[1] = 50 ;
	limestone_prop[2] = 0.22 ;
    
	exp = thermalDeformation(paste_prop,similipaste_prop,1,0.75,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,similipaste_prop,3,0.75,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,similipaste_prop,10,0.75,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,similipaste_prop,33,0.75,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,quartz_prop,1,0.75,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,quartz_prop,3,0.75,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,quartz_prop,10,0.75,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;	
	exp = thermalDeformation(paste_prop,quartz_prop,33,0.75,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
 	exp = thermalDeformation(paste_prop,limestone_prop,1,0.75,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,limestone_prop,3,0.75,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,limestone_prop,10,0.75,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;	
	exp = thermalDeformation(paste_prop,limestone_prop,33,0.75,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,similipaste_prop,1,0.5,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,similipaste_prop,3,0.5,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,similipaste_prop,10,0.5,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;	
	exp = thermalDeformation(paste_prop,similipaste_prop,33,0.5,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,quartz_prop,1,0.5,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,quartz_prop,3,0.5,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,quartz_prop,10,0.5,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;	
	exp = thermalDeformation(paste_prop,quartz_prop,33,0.5,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,limestone_prop,1,0.5,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,limestone_prop,3,0.5,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,limestone_prop,10,0.5,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;	
	exp = thermalDeformation(paste_prop,limestone_prop,33,0.5,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,similipaste_prop,1,0.25,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,similipaste_prop,3,0.25,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,similipaste_prop,10,0.25,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;	
	exp = thermalDeformation(paste_prop,similipaste_prop,33,0.25,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,quartz_prop,1,0.25,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,quartz_prop,3,0.25,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,quartz_prop,10,0.25,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;	
	exp = thermalDeformation(paste_prop,quartz_prop,33,0.25,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,limestone_prop,1,0.25,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,limestone_prop,3,0.25,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,limestone_prop,10,0.25,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;	
	exp = thermalDeformation(paste_prop,limestone_prop,33,0.25,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	return 0 ;
}
