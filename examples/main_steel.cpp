// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../features/features.h"
#include "../physics/stiffness.h"
#include "../physics/stiffness_and_fracture.h"
#include "../physics/materials/gel_behaviour.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../features/sample.h"
#include "../features/inclusion.h"
#include "../features/expansiveZone.h"
#include "../utilities/writer/triangle_writer.h"



#include <fstream>
#ifdef HAVE_SSE3
#include <pmmintrin.h>
#endif
#ifdef HAVE_SSE4
#include <smmintrin.h>
#endif
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h> 
#define DEBUG 


using namespace Amie ;




FeatureTree * featureTree ;
std::vector<bool> cracked ;

double E_min = 10;
double E_max = 0;

double x_max = 0 ;
double y_max = 0 ;

double x_min = 0 ;
double y_min = 0 ;

double timepos = 0.05e-07 ;

double delta_displacement =  1e-5 ;
double displacement_tolerance = 0.05*delta_displacement ; 
double softeningFactor = 1. ;

double percent = 0.01 ;
double displacement  = 0 ;
double prescribedDisplacement = 0;
double derror = 0 ;
double ierror = 0 ;
double preverror = 0 ;
bool firstRun = true ;

std::vector<DelaunayTriangle *> tris__ ;
double apriori_command = 0 ;

std::vector<std::pair<ExpansiveZone *, Inclusion *> > zones ;

std::vector<std::pair<double, double> > expansion_reaction ;
std::vector<std::pair<double, double> > expansion_stress ;
std::vector<std::pair<double, double> > load_displacement ;
std::vector< double > loads ;
std::vector< double > displacements ;

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

double nu = 0.2 ;
double E_agg = 58.9e9 ;
double E_paste = 31.5e9 ;
// BoundingBoxAndRestrictionDefinedBoundaryCondition * load = new BoundingBoxAndRestrictionDefinedBoundaryCondition(SET_STRESS_ETA, TOP, -.15, .15, -10, 10, -10.) ;
// BoundingBoxDefinedBoundaryCondition * load = new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP,0) ;
BoundingBoxDefinedBoundaryCondition * load = new BoundingBoxDefinedBoundaryCondition(SET_ALONG_XI, RIGHT, 0) ;
double factor = 25 ;
MinimumAngle cri(M_PI/6.) ;
bool nothingToAdd = false ;
bool dlist = false ;
int count = 0 ;
double aggregateArea = 0;

GelBehaviour * gel = new GelBehaviour() ;

void computeDisplacement()
{
	x.resize(featureTree->getDisplacements().size()) ;
	x = featureTree->getDisplacements() ;
	Circle C(.01, 0, 0.05) ;
	std::vector<DelaunayTriangle *> t = featureTree->get2DMesh()->getConflictingElements(&C) ;
	std::vector<int> indices ;
	for(size_t i = 0 ; i < t.size() ; i++)
	{
		for(size_t c = 0 ;  c < t[i]->getBoundingPoints().size() ; c++ )
		{
		if(C.in(t[i]->getBoundingPoint(c)))
			indices.push_back(t[i]->getBoundingPoint(c).getId()) ;
		}
	}
	
	std::sort(indices.begin(), indices.end()) ;
	auto e = std::unique(indices.begin(), indices.end()) ;
	displacement = 0 ;
	for(auto i = indices.begin() ; i != e ; i++)
	{
		displacement+=x[(*i)*2+1]/(e-indices.begin()) ;
	}
	
	
// 	Circle C(.000001, -0.0025, -0.05) ;
// 	std::vector<DelaunayTriangle *> t = featureTree->get2DMesh()->getConflictingElements(&C) ;
// 	std::vector<int> indices ;
// 	for(size_t i = 0 ; i < t.size() ; i++)
// 	{
// 		for(size_t c = 0 ;  c < t[i]->getBoundingPoints().size() ; c++ )
// 		{
// 		if(C.in(t[i]->getBoundingPoint(c)))
// 			indices.push_back(t[i]->getBoundingPoint(c).getId()) ;
// 		}
// 	}
// 	
// 	std::sort(indices.begin(), indices.end()) ;
// 	std::vector<int>::iterator e = std::unique(indices.begin(), indices.end()) ;
// 	displacement = 0 ;
// 	for(std::vector<int>::iterator i = indices.begin() ; i != e ; i++)
// 	{
// 		displacement+=x[(*i)*2.]/(e-indices.begin()) ;
// 	}
// 
// 	Circle C0(.000001, 0.0025, -0.05) ;
// 	std::vector<DelaunayTriangle *> t0 = featureTree->get2DMesh()->getConflictingElements(&C) ;
// 	std::vector<int> indices0 ;
// 	for(size_t i = 0 ; i < t0.size() ; i++)
// 	{
// 		for(size_t c = 0 ;  c < t0[i]->getBoundingPoints().size() ; c++ )
// 		{
// 		if(C0.in(t0[i]->getBoundingPoint(c)))
// 			indices0.push_back(t0[i]->getBoundingPoint(c).getId()) ;
// 		}
// 	}
// 	
// 	std::sort(indices0.begin(), indices0.end()) ;
// 	std::vector<int>::iterator e0 = std::unique(indices0.begin(), indices0.end()) ;
// 	double displacement0 = 0 ;
// 	for(std::vector<int>::iterator i = indices0.begin() ; i != e0 ; i++)
// 	{
// 		displacement0+=x[(*i)*2.]/(e0-indices0.begin()) ;
// 	}
// 
// displacement = displacement-displacement0 ;
}

void step()
{
	
	size_t nsteps = 4000 ; //16*10;
	size_t ntries = 5;
	size_t dsteps = 60 ;
	size_t tries = 0 ;
	size_t dit = 0 ;
	featureTree->setMaxIterationsPerStep(dsteps) ;
	for(size_t v = 0 ; v < nsteps ; v++)
	{
		tries = 0 ;
		while(tries < ntries)
		{
			tries++ ;

			featureTree->step() ;

			if(dit < dsteps)
			{
				load->setData(load->getData()+5e-7) ;
				break ;
			}
		}
		
		
		x.resize(featureTree->getDisplacements().size()) ;
		x = featureTree->getDisplacements() ;
	
		
// 		displacement_tolerance = 0.01*(std::abs(delta_displacement)+std::abs(displacement)) ;

		sigma.resize(featureTree->get2DMesh()->begin().size()*featureTree->get2DMesh()->begin()->getBoundingPoints().size()*3) ;
		epsilon.resize(featureTree->get2DMesh()->begin().size()*featureTree->get2DMesh()->begin()->getBoundingPoints().size()*3) ;
		
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
		Vector forces(featureTree->getAssembly()->getForces()) ;
		double appliedForce = std::accumulate(&forces[0], &forces[forces.size()], 0.) ;
		std::cerr << "unknowns :" << x.size() << std::endl ;
		
		
		int npoints = featureTree->get2DMesh()->begin()->getBoundingPoints().size() ;
		
		double area = 0 ;
		double avg_e_xx = 0;
		double avg_e_yy = 0;
		double avg_e_xy = 0;
		double avg_s_xx = 0;
		double avg_s_yy = 0;
		double avg_s_xy = 0;
		double e_xx = 0 ;
		double ex_count = 1 ;
		double avg_e_xx_nogel = 0;
		double avg_e_yy_nogel = 0;
		double avg_e_xy_nogel = 0;
		double avg_s_xx_nogel = 0;
		double avg_s_yy_nogel = 0;
		double avg_s_xy_nogel = 0;
		double nogel_area = 0 ;
		int tsize = 0 ;
		for(auto k = featureTree->get2DMesh()->begin() ; k != featureTree->get2DMesh()->end() ; k++)
		{
	/*		bool in = !k->getEnrichmentFunctions().empty() ;*/
			bool in = false ;
			for(size_t m = 0 ; m < tris__.size() ; m++)
			{
				if(k == tris__[m])
				{
					in = true ;
					break ;
				}
			}
			cracked.push_back(in) ;
			
			
			
			if(!in )
			{
				
				for(size_t p = 0 ;p < k->getBoundingPoints().size() ; p++)
				{
					if(k->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						if(x[k->getBoundingPoint(p).getId()*2] > x_max)
							x_max = x[k->getBoundingPoint(p).getId()*2];
						if(x[k->getBoundingPoint(p).getId()*2] < x_min)
							x_min = x[k->getBoundingPoint(p).getId()*2];
						if(x[k->getBoundingPoint(p).getId()*2+1] > y_max)
							y_max = x[k->getBoundingPoint(p).getId()*2+1];
						if(x[k->getBoundingPoint(p).getId()*2+1] < y_min)
							y_min = x[k->getBoundingPoint(p).getId()*2+1];
						if(k->getBoundingPoint(p).getX() > .9999)
						{
							e_xx=x[k->getBoundingPoint(p).getId()*2] ;
							ex_count = 1 ;
						}
					}
				}
				area += k->area() ;
				if(k->getBehaviour()->type != VOID_BEHAVIOUR)
				{
					tsize++ ;
					if(k->getBehaviour()->param[0][0] > E_max)
						E_max = k->getBehaviour()->param[0][0] ;
					if(k->getBehaviour()->param[0][0] < E_min)
						E_min = k->getBehaviour()->param[0][0] ;
				}
					
				sigma11[k.getPosition()*npoints] = sigma[k.getPosition()*npoints*3];
				sigma22[k.getPosition()*npoints] = sigma[k.getPosition()*npoints*3+1];
				sigma12[k.getPosition()*npoints] = sigma[k.getPosition()*npoints*3+2];
				sigma11[k.getPosition()*npoints+1] = sigma[k.getPosition()*npoints*3+3];
				sigma22[k.getPosition()*npoints+1] = sigma[k.getPosition()*npoints*3+4];
				sigma12[k.getPosition()*npoints+1] = sigma[k.getPosition()*npoints*3+5];
				sigma11[k.getPosition()*npoints+2] = sigma[k.getPosition()*npoints*3+6];
				sigma22[k.getPosition()*npoints+2] = sigma[k.getPosition()*npoints*3+7];
				sigma12[k.getPosition()*npoints+2] = sigma[k.getPosition()*npoints*3+8];
				
				if(npoints >3)
				{
					sigma11[k.getPosition()*npoints+3] = sigma[k.getPosition()*npoints*3+9];
					sigma22[k.getPosition()*npoints+3] = sigma[k.getPosition()*npoints*3+10];
					sigma12[k.getPosition()*npoints+3] = sigma[k.getPosition()*npoints*3+11];
					sigma11[k.getPosition()*npoints+4] = sigma[k.getPosition()*npoints*3+12];
					sigma22[k.getPosition()*npoints+4] = sigma[k.getPosition()*npoints*3+13];
					sigma12[k.getPosition()*npoints+4] = sigma[k.getPosition()*npoints*3+14];
					sigma11[k.getPosition()*npoints+5] = sigma[k.getPosition()*npoints*3+15];
					sigma22[k.getPosition()*npoints+5] = sigma[k.getPosition()*npoints*3+16];
					sigma12[k.getPosition()*npoints+5] = sigma[k.getPosition()*npoints*3+17];
				}
				
				epsilon11[k.getPosition()*npoints] = epsilon[k.getPosition()*npoints*3];
				epsilon22[k.getPosition()*npoints] = epsilon[k.getPosition()*npoints*3+1];
				epsilon12[k.getPosition()*npoints] = epsilon[k.getPosition()*npoints*3+2];
				epsilon11[k.getPosition()*npoints+1] = epsilon[k.getPosition()*npoints*3+3];
				epsilon22[k.getPosition()*npoints+1] = epsilon[k.getPosition()*npoints*3+4];
				epsilon12[k.getPosition()*npoints+1] = epsilon[k.getPosition()*npoints*3+5];
				epsilon11[k.getPosition()*npoints+2] = epsilon[k.getPosition()*npoints*3+6];
				epsilon22[k.getPosition()*npoints+2] = epsilon[k.getPosition()*npoints*3+7];
				epsilon12[k.getPosition()*npoints+2] = epsilon[k.getPosition()*npoints*3+8];
				
				if(npoints > 3)
				{
					epsilon11[k.getPosition()*npoints+3] = epsilon[k.getPosition()*npoints*3+9];
					epsilon22[k.getPosition()*npoints+3] = epsilon[k.getPosition()*npoints*3+10];
					epsilon12[k.getPosition()*npoints+3] = epsilon[k.getPosition()*npoints*3+11];
					epsilon11[k.getPosition()*npoints+4] = epsilon[k.getPosition()*npoints*3+12];
					epsilon22[k.getPosition()*npoints+4] = epsilon[k.getPosition()*npoints*3+13];
					epsilon12[k.getPosition()*npoints+4] = epsilon[k.getPosition()*npoints*3+14];
					epsilon11[k.getPosition()*npoints+5] = epsilon[k.getPosition()*npoints*3+15];
					epsilon22[k.getPosition()*npoints+5] = epsilon[k.getPosition()*npoints*3+16];
					epsilon12[k.getPosition()*npoints+5] = epsilon[k.getPosition()*npoints*3+17];
				}  
				
				for(size_t l = 0 ; l < k->getBoundingPoints().size() ; l++)
				{
					Vector vm0(0., 3) ;
					k->getState().getField( PRINCIPAL_REAL_STRESS_FIELD, k->getBoundingPoint(l), vm0, false) ;
					vonMises[k.getPosition()*k->getBoundingPoints().size()+l]  = sqrt(((vm0[0]-vm0[1])*(vm0[0]-vm0[1]))/2.) ;
	
					Vector agl(0., 1) ;
					k->getState().getField( PRINCIPAL_STRESS_ANGLE_FIELD, k->getBoundingPoint(l), agl, false) ;
					angle[k.getPosition()*k->getBoundingPoints().size()+l]  = agl[0] ;
				}
				
				double ar = k->area() ;
				for(int l = 0 ; l < npoints ;l++)
				{
					avg_e_xx += (epsilon11[k.getPosition()*npoints+l]/npoints)*ar;
					avg_e_yy += (epsilon22[k.getPosition()*npoints+l]/npoints)*ar;
					avg_e_xy += (epsilon12[k.getPosition()*npoints+l]/npoints)*ar;
					avg_s_xx += (sigma11[k.getPosition()*npoints+l]/npoints)*ar;
					avg_s_yy += (sigma22[k.getPosition()*npoints+l]/npoints)*ar;
					avg_s_xy += (sigma12[k.getPosition()*npoints+l]/npoints)*ar;
				}
				
				if(k->getEnrichmentFunctions().size() == 0)
				{
					for(int l = 0 ; l < npoints ;l++)
					{
						avg_e_xx_nogel += (epsilon11[k.getPosition()*npoints+l]/npoints)*ar;
						avg_e_yy_nogel += (epsilon22[k.getPosition()*npoints+l]/npoints)*ar;
						avg_e_xy_nogel += (epsilon12[k.getPosition()*npoints+l]/npoints)*ar;
						avg_s_xx_nogel += (sigma11[k.getPosition()*npoints+l]/npoints)*ar;
						avg_s_yy_nogel += (sigma22[k.getPosition()*npoints+l]/npoints)*ar;
						avg_s_xy_nogel += (sigma12[k.getPosition()*npoints+l]/npoints)*ar;
						
					}
					nogel_area+= ar ;
				}

			}
			else
			{
				sigma11[k.getPosition()*npoints] = 0 ;
				sigma22[k.getPosition()*npoints] = 0 ;
				sigma12[k.getPosition()*npoints] = 0 ;
				sigma11[k.getPosition()*npoints+1] = 0 ;
				sigma22[k.getPosition()*npoints+1] = 0 ;
				sigma12[k.getPosition()*npoints+1] = 0 ;
				sigma11[k.getPosition()*npoints+2] = 0 ;
				sigma22[k.getPosition()*npoints+2] = 0 ;
				sigma12[k.getPosition()*npoints+2] = 0 ;
				
				if(npoints >3)
				{
					sigma11[k.getPosition()*npoints+3] = 0 ;
					sigma22[k.getPosition()*npoints+3] = 0 ;
					sigma12[k.getPosition()*npoints+3] = 0 ;
					sigma11[k.getPosition()*npoints+4] = 0 ;
					sigma22[k.getPosition()*npoints+4] = 0 ;
					sigma12[k.getPosition()*npoints+4] = 0 ;
					sigma11[k.getPosition()*npoints+5] = 0 ;
					sigma22[k.getPosition()*npoints+5] = 0 ;
					sigma12[k.getPosition()*npoints+5] =0 ;
				}
				
				epsilon11[k.getPosition()*npoints] = 0 ;
				epsilon22[k.getPosition()*npoints] = 0 ;
				epsilon12[k.getPosition()*npoints] = 0 ;
				epsilon11[k.getPosition()*npoints+1] = 0 ;
				epsilon22[k.getPosition()*npoints+1] = 0 ;
				epsilon12[k.getPosition()*npoints+1] = 0 ;
				epsilon11[k.getPosition()*npoints+2] = 0 ;
				epsilon22[k.getPosition()*npoints+2] = 0 ;
				epsilon12[k.getPosition()*npoints+2] = 0 ;
				
				if(npoints > 3)
				{
					epsilon11[k.getPosition()*npoints+3] = 0 ;
					epsilon22[k.getPosition()*npoints+3] = 0 ;
					epsilon12[k.getPosition()*npoints+3] =0 ;
					epsilon11[k.getPosition()*npoints+4] = 0 ;
					epsilon22[k.getPosition()*npoints+4] = 0 ;
					epsilon12[k.getPosition()*npoints+4] =0 ;
					epsilon11[k.getPosition()*npoints+5] = 0 ;
					epsilon22[k.getPosition()*npoints+5] =0 ;
					epsilon12[k.getPosition()*npoints+5] = 0 ;
				}  
				
				for(size_t l = 0 ; l < k->getBoundingPoints().size() ; l++)
				{
					vonMises[k.getPosition()*k->getBoundingPoints().size()+l]  = 0 ;
					angle[k.getPosition()*k->getBoundingPoints().size()+l]  = 0 ;
				}
			}
		}
		
		
		if(dit < dsteps)
			loads.push_back(avg_s_xx/area);
		if(dit < dsteps)
			displacements.push_back(avg_e_xx/area);
		if(v%25 == 0)
		{
			std::cout << std::endl ;
			std::cout << "load :" << appliedForce/1000. << std::endl ;
			std::cout << "displacement :" << 1000.*e_xx/(double)ex_count << std::endl ;
			std::cout << "max value :" << x_max << std::endl ;
			std::cout << "min value :" << x_min << std::endl ;
			std::cout << "max sigma11 :" << sigma11.max()/1000000. << std::endl ;
			std::cout << "min sigma11 :" << sigma11.min()/1000000. << std::endl ;
			std::cout << "max sigma12 :" << sigma12.max()/1000000. << std::endl ;
			std::cout << "min sigma12 :" << sigma12.min()/1000000. << std::endl ;
			std::cout << "max sigma22 :" << sigma22.max()/1000000. << std::endl ;
			std::cout << "min sigma22 :" << sigma22.min()/1000000. << std::endl ;
			
			std::cout << "max epsilon11 :" << epsilon11.max() << std::endl ;
			std::cout << "min epsilon11 :" << epsilon11.min() << std::endl ;
			std::cout << "max epsilon12 :" << epsilon12.max() << std::endl ;
			std::cout << "min epsilon12 :" << epsilon12.min() << std::endl ;
			std::cout << "max epsilon22 :" << epsilon22.max() << std::endl ;
			std::cout << "min epsilon22 :" << epsilon22.min() << std::endl ;
			
			std::cout << "max von Mises :" << vonMises.max()/1000000. << std::endl ;
			std::cout << "min von Mises :" << vonMises.min()/1000000. << std::endl ;
			
			std::cout << "average sigma11 : " << avg_s_xx/area/1000000. << std::endl ;
			std::cout << "average sigma22 : " << avg_s_yy/area/1000000. << std::endl ;
			std::cout << "average sigma12 : " << avg_s_xy/area/1000000. << std::endl ;
			std::cout << "average epsilon11 : " << avg_e_xx/area<< std::endl ;
			std::cout << "average epsilon22 : " << avg_e_yy/area << std::endl ;
			std::cout << "average epsilon12 : " << avg_e_xy/area << std::endl ;
		}
// 		for (int i = 0 ; i < displacements.size() ; i++)
// 		{
// 			std::cout << loads[i] << "  " << displacements[i] << std::endl ;
// 		}
		
		if(dit < dsteps)	
			std::cout << displacements.back() << "  "<< appliedForce/1000.  <<  std::endl ;
		
		std::fstream ldfile  ;
		ldfile.open("ldn", std::ios::out) ;
		for(size_t j = 0 ; j < loads.size() ; j++)
		{
			ldfile << displacements[j] << "   " << loads[j] << "\n" ;
		}
		ldfile.close();
		
		if(v%1 == 0)
		{
			std::stringstream filename ;
			if(dit >= dsteps)
				filename << "intermediate-" ;
			
			filename << "triangles-" ;
			filename << round(appliedForce/1000.) ;
			filename << "-" ;
			filename << 1000.*e_xx/(double)ex_count ;
			
		TriangleWriter writer(filename.str(), featureTree) ;
		writer.getField(STRAIN_FIELD) ;
		writer.getField(REAL_STRESS_FIELD) ;
		writer.getField(TWFT_VON_MISES) ;
		writer.getField(TWFT_STIFFNESS) ;
		writer.write() ;

		}
		//(1./epsilon11.getX())*( stressMoyenne.getX()-stressMoyenne.getY()*modulePoisson);
		
		double delta_r = sqrt(aggregateArea*0.03/((double)zones.size()*M_PI))/nsteps ;
		double reactedArea = 0 ;
			
// 		if (tries < ntries)
			for(size_t z = 0 ; z < zones.size() ; z++)
			{
				zones[z].first->setRadius(zones[z].first->getRadius()+delta_r) ;	
		// 		zones[z].first->reset() ;
				reactedArea += zones[z].first->area() ;
			}
		
		std::cerr << "reacted Area : " << reactedArea << std::endl ;
		
// 		if (tries < ntries)
// 		{
			expansion_reaction.push_back(std::make_pair(reactedArea, avg_e_xx/area)) ;
			expansion_stress.push_back(std::make_pair(avg_e_xx_nogel/nogel_area, avg_s_xx_nogel/nogel_area)) ;
// 		}
		
// 		if (tries >= ntries)
// 			break ;
	}
}

std::vector<std::pair<ExpansiveZone *, Inclusion *> > generateExpansiveZones(int n, std::vector<Inclusion * > & incs , FeatureTree & F)
{	
	std::vector<std::pair<ExpansiveZone *, Inclusion *> > ret ;
	aggregateArea = 0 ;
	for(size_t i = 0 ; i < incs.size() ; i++)
	{
		aggregateArea += incs[i]->area() ;
		for(int j = 0 ; j < n ; j++)
		{
			double radius = 0.000001 ;
			double rangle = (2.*rand()/RAND_MAX-1.)*M_PI ;
			double rradius = (double)rand()/RAND_MAX*(incs[i]->getRadius()-2.*radius*64) ;
			Point center = incs[i]->getCenter()+Point(rradius*cos(rangle), rradius*sin(rangle)) ; 
			
			bool alone  = true ;
			
			for(size_t k = 0 ; k < ret.size() ; k++ )
			{
				if (squareDist(center, ret[k].first->Circle::getCenter()) < 2.*(2.*radius)*(2.*radius)*64*64)
				{
					alone = false ;
					break ;
				}
			}
			if (alone)
			{
				ExpansiveZone * z = new ExpansiveZone(incs[i], radius, center.getX(), center.getY(), gel) ;
				ret.push_back(std::make_pair(z, incs[i])) ;
				F.addFeature(incs[i],z) ; 
			}
		}
	}
	std::cout << "initial Reacted Area = " << M_PI*0.000001*0.000001*ret.size() << " in "<< ret.size() << " zones"<< std::endl ;
	return ret ;	
}

int main(int argc, char *argv[])
{
	
	Matrix m0_steel(3,3) ;
	double E_steel = 200e9 ;
	double nu_steel = 0.3 ; 
	m0_steel[0][0] = E_steel/(1.-nu_steel*nu_steel) ;          m0_steel[0][1] = E_steel/(1.-nu_steel*nu_steel)*nu_steel ; m0_steel[0][2] = 0 ;
	m0_steel[1][0] = E_steel/(1.-nu_steel*nu_steel)*nu_steel ; m0_steel[1][1] = E_steel/(1.-nu_steel*nu_steel) ;          m0_steel[1][2] = 0 ; 
	m0_steel[2][0] = 0 ;                                       m0_steel[2][1] = 0 ;                                       m0_steel[2][2] = E_steel/(1.-nu_steel*nu_steel)*(1.-nu_steel)/2. ; 
	
	Matrix m0_barSteel(3,3) ;
	m0_barSteel[0][0] = E_steel*(1-nu_steel)/((1+nu_steel)*(1.-2.*nu_steel)) ; m0_barSteel[0][1] =E_steel*nu_steel/((1.+nu_steel)*(1.-2.*nu_steel)) ;      m0_barSteel[0][2] = 0 ;
	m0_barSteel[1][0] = E_steel*nu_steel/((1+nu_steel)*(1.-2.*nu_steel)) ;     m0_barSteel[1][1] = E_steel*(1-nu_steel)/((1.+nu_steel)*(1.-2.*nu_steel)) ; m0_barSteel[1][2] = 0 ; 
	m0_barSteel[2][0] = 0 ;                                                    m0_barSteel[2][1] = 0 ;                                                     m0_barSteel[2][2] = E_steel/((1.+nu_steel)) ; 
	
	Matrix m0_paste(3,3) ;
	m0_paste[0][0] = E_paste/(1-nu*nu) ; m0_paste[0][1] =E_paste/(1-nu*nu)*nu ; m0_paste[0][2] = 0 ;
	m0_paste[1][0] = E_paste/(1-nu*nu)*nu ; m0_paste[1][1] = E_paste/(1-nu*nu) ; m0_paste[1][2] = 0 ; 
	m0_paste[2][0] = 0 ; m0_paste[2][1] = 0 ; m0_paste[2][2] = E_paste/(1-nu*nu)*(1.-nu)/2. ; 

	Sample sample(nullptr, 3.9*.5, 1.2, 0,0) ;
	
	Sample rebar0(1, 0.025, (3.9*.5)*.5,  -1.2*.5+0.064) ; 
	rebar0.setBehaviour(new Stiffness(m0_steel*0.15));
	Sample rebar1(1, 0.025, (3.9*.5)*.5,  -1.2*.5+0.064+0.085) ; 
	rebar1.setBehaviour(new Stiffness(m0_steel*0.15));
	
	FeatureTree F(&sample) ;
	featureTree = &F ;

	sample.setBehaviour(new WeibullDistributedStiffness(E_paste, nu, SPACE_TWO_DIMENSIONAL, -37.0e6, 4.3e6)) ;
	dynamic_cast<WeibullDistributedStiffness *>(sample.getBehaviour())->variability = 0. ;
	dynamic_cast<WeibullDistributedStiffness *>(sample.getBehaviour())->materialRadius = .3 ;

	F.addBoundaryCondition(load) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT) );
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM) );
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, LEFT) );

	

// 	F.addFeature(&sample,&rebar0) ;
// 	F.addFeature(&sample,&rebar1) ;
	
	
	F.setSamplingNumber(400) ;
	F.setOrder(LINEAR) ;
	
	step() ;
	
// 	delete dt ;
	
	return 0 ;
}
