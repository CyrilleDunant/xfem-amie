// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../utilities/samplingcriterion.h"
#include "../features/features.h"
#include "../physics/radialstiffnessgradient.h"
#include "../physics/physics_base.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/stiffness.h"
#include "../physics/materials/gel_behaviour.h"
#include "../physics/stiffness_and_fracture.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../features/pore.h"
#include "../features/sample.h"
#include "../features/inclusion.h"
#include "../features/layeredinclusion.h"
#include "../features/expansiveZone.h"
#include "../features/crack.h"
#include "../features/enrichmentInclusion.h"
#include "../mesher/delaunay_3d.h"
#include "../solvers/assembly.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"

#include <fstream>

#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h> 
#define DEBUG 

using namespace Amie ;

FeatureTree * featureTree ;
std::vector<DelaunayTriangle *> triangles ;
std::vector<bool> cracked ;

GelBehaviour * gel = new GelBehaviour() ;

double E_min = 10;
double E_max = 0;

double x_max = 0 ;
double y_max = 0 ;

double x_min = 0 ;
double y_min = 0 ;

double timepos = 0.05e-07 ;
double delta_displacement =  0.5e-6 ;
double displacement_tolerance = 0.01*delta_displacement ; 
double softeningFactor = 1. ;

double percent = 0.01 ;
double load = 0 ;
double displacement  = 0 ;
double prescribedDisplacement = 0 ;// -1.15e-10 ; //-1.25e-7;
double derror = 0 ;
double ierror = 0 ;
double preverror = 0 ;
bool firstRun = true ;

std::vector<DelaunayTriangle *> tris__ ;
double apriori_command = 0 ;
std::pair<std::vector<Inclusion * >, std::vector<Pore * > > i_et_p ;

std::vector<std::pair<ExpansiveZone *, Inclusion *> > zones ;

std::vector<std::pair<double, double> > expansion_reaction ;
std::vector<std::pair<double, double> > expansion_stress ;
std::vector<std::pair<double, double> > load_displacement ;

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
double E_agg = 58900000000 ;
double E_paste = 12000000000 ;

double factor = 200 ;
MinimumAngle cri(M_PI/6.) ;
bool nothingToAdd = false ;
bool dlist = false ;
int count = 0 ;
double aggregateArea = 0;

void computeDisplacement()
{
	x.resize(featureTree->getDisplacements().size()) ;
	x = featureTree->getDisplacements() ;
	Circle c(.0001, 0., 0.02) ;
	std::vector<DelaunayTriangle *> t = featureTree->getElements2D(&c) ;
	std::vector<int> indices ;
	for(size_t i = 0 ; i < t.size() ; i++)
	{
		for(size_t c = 0 ;  c < t[i]->getBoundingPoints().size() ; c++ )
		{
			indices.push_back(t[i]->getBoundingPoint(c).getId()) ;
		}
	}
	
	std::sort(indices.begin(), indices.end()) ;
	auto e = std::unique(indices.begin(), indices.end()) ;
	displacement = 0 ;
	for(auto i = indices.begin() ; i != e ; i++)
	{
		displacement+=x[(*i)*2.+1]/(e-indices.begin()) ;
	}
	
}

double pidUpdate()
{
	double currentLoad = load ;
	if(load_displacement.size() > 2 && std::abs(load_displacement.back().second) > 1e-15)
		apriori_command = load_displacement.back().first
		/ load_displacement.back().second
				* prescribedDisplacement ;
	double error = prescribedDisplacement-displacement ;
	derror =  error-preverror ;
	ierror += (error+preverror)*.5 ;
	double K_p = 100000000. ;
	
	if(load_displacement.size() > 2 && std::abs(load_displacement.back().second) > 1e-12 && std::abs(load_displacement[2].second) > 1e-12)
	{
		softeningFactor = (load_displacement.back().first 
			/ load_displacement.back().second )
			/ (load_displacement[2].first 
			/ load_displacement[2].second ) ;
	}
	softeningFactor = std::min(1., softeningFactor) ;
	K_p *= softeningFactor ;
	load = /*apriori_command +*/ K_p*error + K_p*ierror+ K_p* .5 *derror;
// 	if(load > 0)
// 		load = -load ;
// 	if(std::abs(error) > std::abs(preverror) && std::abs(load) > std::abs(currentLoad))
// 	{
// 		std::cout << "windup ! Resetting integral" << std::endl ;
// 		ierror = 0 ;
// 		load = /*apriori_command +*/ K_p*error + K_p*ierror+ K_p* .5 *derror;
// 		if(load > 0)
// 			load = -load ;
// 	}
	
	preverror = error ;

	return error ;
}

void step()
{
	
	size_t nsteps = 64;
	size_t nit = 20 ;
	size_t ntries = 25;

	for(size_t i = 0 ; i < nit ; i++)
	{
		size_t tries = 0 ;
		bool go_on = true ;
		std::vector<std::pair<double,double> > saved_load_displacement = load_displacement ;
		double originalPrescribedDisplacement = prescribedDisplacement ;

		featureTree->step() ;

		std::cout << " " << tries << " tries." << std::endl ;
		
// 		saved_load_displacement.push_back(load_displacement.back()) ;
// 		load_displacement = saved_load_displacement ;
		prescribedDisplacement -= delta_displacement ;
		x.resize(featureTree->getDisplacements().size()) ;
		x = featureTree->getDisplacements() ;
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

		cracked.clear() ;
		
		int npoints = triangles[0]->getBoundingPoints().size() ;
		
		double area = 0 ;
		double avg_e_xx = 0;
		double avg_e_yy = 0;
		double avg_e_xy = 0;
		double avg_s_xx = 0;
		double avg_s_yy = 0;
		double avg_s_xy = 0;
		double e_xx = 0 ;
		double ex_count = 0 ;
		double avg_e_xx_nogel = 0;
		double avg_e_yy_nogel = 0;
		double avg_e_xy_nogel = 0;
		double avg_s_xx_nogel = 0;
		double avg_s_yy_nogel = 0;
		double avg_s_xy_nogel = 0;
		double nogel_area = 0 ;
		
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
			
			
			
			if(!in /*&& !triangles[k]->getBehaviour()->fractured()*/)
			{
				
				for(size_t p = 0 ;p < triangles[k]->getBoundingPoints().size() ; p++)
				{
					if(x[triangles[k]->getBoundingPoint(p).getId()*2] > x_max)
						x_max = x[triangles[k]->getBoundingPoint(p).getId()*2];
					if(x[triangles[k]->getBoundingPoint(p).getId()*2] < x_min)
						x_min = x[triangles[k]->getBoundingPoint(p).getId()*2];
					if(x[triangles[k]->getBoundingPoint(p).getId()*2+1] > y_max)
						y_max = x[triangles[k]->getBoundingPoint(p).getId()*2+1];
					if(x[triangles[k]->getBoundingPoint(p).getId()*2+1] < y_min)
						y_min = x[triangles[k]->getBoundingPoint(p).getId()*2+1];
					if(triangles[k]->getBoundingPoint(p).getX() > 0.0799)
					{
						e_xx+=x[triangles[k]->getBoundingPoint(p).getId()*2] ;
						ex_count++ ;
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
					Vector vm0(0., 3) ;
					triangles[k]->getState().getField( PRINCIPAL_REAL_STRESS_FIELD, triangles[k]->getBoundingPoint(l), vm0, false) ;
					vonMises[k*triangles[k]->getBoundingPoints().size()+l]  = sqrt(((vm0[0]-vm0[1])*(vm0[0]-vm0[1]))/2.) ;
	
					Vector agl(0., 1) ;
					triangles[k]->getState().getField( PRINCIPAL_ANGLE_FIELD, triangles[k]->getBoundingPoint(l), agl, false) ;
					angle[k*triangles[k]->getBoundingPoints().size()+l]  = agl[0] ;
				}
				
				double ar = triangles[k]->area() ;
				for(int l = 0 ; l < npoints ;l++)
				{
					avg_e_xx += (epsilon11[k*npoints+l]/npoints)*ar;
					avg_e_yy += (epsilon22[k*npoints+l]/npoints)*ar;
					avg_e_xy += (epsilon12[k*npoints+l]/npoints)*ar;
					avg_s_xx += (sigma11[k*npoints+l]/npoints)*ar;
					avg_s_yy += (sigma22[k*npoints+l]/npoints)*ar;
					avg_s_xy += (sigma12[k*npoints+l]/npoints)*ar;
				}
				
				if(triangles[k]->getEnrichmentFunctions().size() == 0)
				{
					for(int l = 0 ; l < npoints ;l++)
					{
						avg_e_xx_nogel += (epsilon11[k*npoints+l]/npoints)*ar;
						avg_e_yy_nogel += (epsilon22[k*npoints+l]/npoints)*ar;
						avg_e_xy_nogel += (epsilon12[k*npoints+l]/npoints)*ar;
						avg_s_xx_nogel += (sigma11[k*npoints+l]/npoints)*ar;
						avg_s_yy_nogel += (sigma22[k*npoints+l]/npoints)*ar;
						avg_s_xy_nogel += (sigma12[k*npoints+l]/npoints)*ar;
						
					}
					nogel_area+= ar ;
				}

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
		
		std::cout << "average sigma11 (no gel): " << avg_s_xx_nogel/nogel_area << std::endl ;
		std::cout << "average sigma22 (no gel): " << avg_s_yy_nogel/nogel_area << std::endl ;
		std::cout << "average sigma12 (no gel): " << avg_s_xy_nogel/nogel_area << std::endl ;
		std::cout << "average epsilon11 (no gel): " << avg_e_xx_nogel/nogel_area << std::endl ;
		std::cout << "average epsilon22 (no gel): " << avg_e_yy_nogel/nogel_area << std::endl ;
		std::cout << "average epsilon12 (no gel): " << avg_e_xy_nogel/nogel_area << std::endl ;
		
		std::cout << "apparent extension " << e_xx/ex_count << std::endl ;
		//(1./epsilon11.getX())*( stressMoyenne.getX()-stressMoyenne.getY()*modulePoisson);
		
		double delta_r = sqrt(aggregateArea*0.03/((double)zones.size()*M_PI))/nsteps ;
		double reactedArea = 0 ;
			
		if (tries < ntries)
			for(size_t z = 0 ; z < zones.size() ; z++)
			{
				zones[z].first->setRadius(zones[z].first->getRadius()+delta_r) ;	
		// 		zones[z].first->reset() ;
				reactedArea += zones[z].first->area() ;
			}
		
		std::cout << "reacted Area : " << reactedArea << std::endl ;
		
		if (tries < ntries)
		{
			expansion_reaction.push_back(std::make_pair(reactedArea, avg_e_xx/area)) ;
			expansion_stress.push_back(std::make_pair(avg_e_xx_nogel/nogel_area, avg_s_xx_nogel/nogel_area)) ;
		}
		
		if (tries >= 100)
			break ;

	for(size_t i = 0 ; i < expansion_reaction.size() ; i++)
		std::cout << expansion_reaction[i].first << "   " 
		<< expansion_reaction[i].second << "   " 
		<< expansion_stress[i].first << "   " 
		<< expansion_stress[i].second << "   " << std::endl ;
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

std::pair<std::vector<Inclusion * >, std::vector<Pore * > > generateInclusionsAndPores(size_t n, double fraction, double E, double nu, Feature * father, FeatureTree * F)
{
	size_t nombre_de_pores = static_cast<size_t>(round(n*fraction)) ;
	size_t nombre_d_inclusions = static_cast<size_t>(round(n*(1. - fraction))) ;
	
	std::pair<std::vector<Inclusion * >, std::vector<Pore * > > ret ;
	ret.first = std::vector<Inclusion * >() ;
	ret.second = std::vector<Pore * >() ;
	double v = 0 ;
	std::vector<Circle *> cercles ;
	for(size_t j =0 ; j < n ; j++)
	{
		
		double radius = .0005 + .0025*rand()/RAND_MAX ;
		
		Point center = Point(
		                      (2.*rand()/RAND_MAX-1.)*(.08-2.*radius-0.00001),
		                      (2.*rand()/RAND_MAX-1.)*(.02-2.*radius-0.00001)
		                    ); 
		bool alone  = true ;
		
		for(size_t k = 0 ; k < cercles.size() ; k++ )
		{
			if (squareDist(center, cercles[k]->getCenter()) < (radius+cercles[k]->getRadius()+0.00001)*(radius+cercles[k]->getRadius()+0.00001))
			{
				alone = false ;
				break ;
			}
		}
		if (alone)
		{
			cercles.push_back(new Circle(radius, center)) ;
			v+= M_PI*radius*radius ;
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
// 		(*ret.first.rbegin())->setBehaviour(new StiffnessAndFracture(*tensor, new MohrCoulomb(1000000, -10000000))) ;
		(*ret.first.rbegin())->setBehaviour(new WeibullDistributedStiffness(E, nu, SPACE_TWO_DIMENSIONAL, -8000000, 1000000)) ;
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
	
	std::cout << "initial aggregate volume was : " << v << std::endl ;
	aggregateArea = v ;
	return ret ;
}

int main(int argc, char *argv[])
{
// 	srandom(time(nullptr)) ;

	Matrix m0_agg(3,3) ;
	m0_agg[0][0] = E_agg/(1-nu*nu) ; m0_agg[0][1] =E_agg/(1-nu*nu)*nu ; m0_agg[0][2] = 0 ;
	m0_agg[1][0] = E_agg/(1-nu*nu)*nu ; m0_agg[1][1] = E_agg/(1-nu*nu) ; m0_agg[1][2] = 0 ; 
	m0_agg[2][0] = 0 ; m0_agg[2][1] = 0 ; m0_agg[2][2] = E_agg/(1-nu*nu)*(1.-nu)/2. ; 
	
	Matrix m0_dest(3,3) ;
	m0_agg[0][0] = E_agg/(1-nu*nu) ; m0_agg[0][1] =0 ; m0_agg[0][2] = 0 ;
	m0_agg[1][0] = 0 ; m0_agg[1][1] = E_agg/(1-nu*nu) ; m0_agg[1][2] = 0 ; 
	m0_agg[2][0] = 0 ; m0_agg[2][1] = 0 ; m0_agg[2][2] = 1 ; 
	
	Matrix m0_paste(3,3) ;
	m0_paste[0][0] = E_paste/(1-nu*nu) ; m0_paste[0][1] =E_paste/(1-nu*nu)*nu ; m0_paste[0][2] = 0 ;
	m0_paste[1][0] = E_paste/(1-nu*nu)*nu ; m0_paste[1][1] = E_paste/(1-nu*nu) ; m0_paste[1][2] = 0 ; 
	m0_paste[2][0] = 0 ; m0_paste[2][1] = 0 ; m0_paste[2][2] = E_paste/(1-nu*nu)*(1.-nu)/2. ; 

	Sample sample(nullptr, 0.16, 0.04,0,0) ;
	Amie::Rectangle box( 0.28, 0.07,0,0) ;
	
	Inclusion inclusion(.00001, 0.02, -.02) ;
// return 0 ; 
	FeatureTree F(&sample) ;
	featureTree = &F ;


	double itzSize = 0.00004;
	int inclusionNumber = 4092 ;
	double masseInitiale = .00000743*4;
	double densite = 1.;
//	std::vector<Inclusion *> inclusions = GranuloBolome(4.79263e-07*4, 1, BOLOME_D)(.0025, .0001, inclusionNumber, itzSize);
	std::vector<Inclusion *> inclusions = PSDGenerator::get2DInclusions(.0025, 4.79263e-07, new PSDBolomeD(), PSDEndCriteria(-1, 0.001, inclusionNumber)) ;


	std::vector<Feature *> feats ;
	for(size_t i = 0; i < inclusions.size() ; i++)
		feats.push_back(inclusions[i]) ;

	int nAgg = 0 ;
	feats=placement(sample.getPrimitive(), feats, &nAgg, 0, 64000);

	double volume = 0 ;
	for(size_t i = 0 ; i < feats.size() ; i++)
		volume += feats[i]->area() ;
	if(!feats.empty())
		std::cout << "largest r = " << feats.back()->getRadius() 
		<< ", smallest r =" << feats.front()->getRadius() 
		<< ", filling = " << volume/sample.area()*100.<< "%"<< std::endl ; 
// 	feats=placement(new Circle(.01, 0,0), feats, &nAgg, 32000);
	inclusions.clear() ;
	for(size_t i = 0; i < feats.size() ; i++)
		inclusions.push_back(dynamic_cast<Inclusion *>(feats[i])) ;
	
	std::cout << "incs : " << inclusions.size() << std::endl ;
	double placed_area = 0 ;
	sample.setBehaviour(new WeibullDistributedStiffness(E_paste, nu, SPACE_TWO_DIMENSIONAL, -40000*8., 40000)) ;
// 	sample.setBehaviour(new StiffnessAndFracture(m0_paste, new MohrCoulomb(40000, -40000*8))) ;
	Inclusion * pore = new Inclusion(0.01, 0.458, -0.153) ;
	pore->setBehaviour(new Stiffness(m0_paste)) ;
	F.addFeature(&sample,pore) ;
	Inclusion * pore0 = new Inclusion(0.01, -0.061, -0.153) ;
	pore0->setBehaviour(new Stiffness(m0_paste)) ;
	F.addFeature(pore,pore0) ;
	
	Inclusion * pore1 = new Inclusion(0.01, -0.458, 0.153) ;
	pore1->setBehaviour(new Stiffness(m0_paste)) ;
	F.addFeature(pore0,pore1) ;
	Inclusion * pore2 = new Inclusion(0.01, 0.061, 0.153) ;
	pore2->setBehaviour(new Stiffness(m0_paste)) ;
	Inclusion * pore3 = new Inclusion(0.01, -0.458, -0.153) ;
	pore3->setBehaviour(new Stiffness(m0_paste)) ;
	F.addFeature(pore1,pore2) ;
	F.addFeature(pore2,pore3) ;
	TriangularPore * pore4 = new TriangularPore(Point(0, -0.071), Point(-0.005, -0.154), Point(0.005, -0.154)) ;
	F.addFeature(pore2,pore4) ;
	for(size_t i = 0 ; i < inclusions.size(); i++)
	{
		std::vector<double> radii ;
		std::vector<Form *> behavs ;
		radii.push_back(inclusions[i]->getRadius()) ;
		if(inclusions[i]->getRadius()-itzSize > 0)
			radii.push_back(inclusions[i]->getRadius()-itzSize) ;
// 		behavs.push_back(new WeibullDistributedStiffness(m0_agg,80000)) ;
		behavs.push_back(new StiffnessAndFracture(m0_agg, new MohrCoulomb(80000, -80000*8))) ;
// 		behavs.push_back(new WeibullDistributedStiffness(m0_paste*.5, 40000*.5)) ;
// 		behavs.push_back(new WeibullDistributedStiffness(m0_paste, 40000)) ;
		RadialStiffnessGradient * b = new RadialStiffnessGradient(.75*E_paste, nu, inclusions[i]->getRadius()-itzSize, 
			E_paste, nu, inclusions[i]->getRadius(),
			inclusions[i]->getCenter()
			) ;
		b->setFractureCriterion(new MohrCoulomb(40000*.75, -40000*4)) ;
		behavs.push_back(b) ;

		LayeredInclusion * newinc = new LayeredInclusion(radii, inclusions[i]->getCenter()) ;
		newinc->setBehaviours(behavs) ;
		F.addFeature(pore1,newinc) ;
		
// 		inclusions[i]->setRadius(inclusions[i]->getRadius()-itzSize) ;
// 		inclusions[i]->setBehaviour(new WeibullDistributedStiffness(m0_agg,80000)) ;
// 		F.addFeature(pore1,inclusions[i]) ;
// 		F.addFeature(pore1,new Pore(inclusions[i]->getRadius()-itzSize*.75, inclusions[i]->getCenter())) ;
		placed_area += inclusions[i]->area() ;
	}


	if(!inclusions.empty())
	{
		std::cout << "largest inclusion with r = " << inclusions.front()->getRadius() << std::endl ;
		std::cout << "smallest inclusion with r = " << inclusions.back()->getRadius() << std::endl ;
		std::cout << "placed area = " <<  placed_area << std::endl ;
	}
// 	inclusions.erase(inclusions.begin()+1, inclusions.end()) ;
// 	zones = generateExpansiveZones(1, inclusions, F) ;

	F.setSamplingNumber(400) ;

	F.setOrder(LINEAR) ;

// 	F.refine(2, new MinimumAngle(M_PI/8.)) ;
	
// 	
	step() ;
	
// 	delete dt ;
	
	return 0 ;
}
