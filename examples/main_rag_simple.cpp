// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../utilities/samplingcriterion.h"
#include "../features/features.h"
#include "../physics/physics_base.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/stiffness_with_variable_imposed_deformation.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../physics/stiffness.h"
#include "../physics/materials/gel_behaviour.h"
#include "../features/pore.h"
#include "../features/sample.h"
#include "../features/inclusion.h"
#include "../features/expansiveZone.h"
#include "../features/expansiveRing.h"
#include "../features/crack.h"
#include "../features/enrichmentInclusion.h"
#include "../mesher/delaunay_3d.h"
#include "../solvers/assembly.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"
#include "../utilities/writer/triangle_writer.h"

#include <fstream>

#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h> 
#define DEBUG 

using namespace Mu ;

FeatureTree * featureTree ;
std::vector<DelaunayTriangle *> triangles ;
std::vector<bool> cracked ;

double E_min = 10;
double E_max = 0;

double x_max = 0 ;
double y_max = 0 ;

double x_min = 0 ;
double y_min = 0 ;

double timepos = 0.00 ;
double percent = 0.10 ;
double placed_area = 0 ;

bool firstRun = true ;

std::vector<DelaunayTriangle *> tris__ ;

std::pair<std::vector<Inclusion * >, std::vector<Pore * > > i_et_p ;

std::vector<std::pair<ExpansiveZone *, Inclusion *> > zones ;
std::vector<Inclusion *> inclusions ;
std::vector<ExpansiveRing *> reactionRims ;

std::vector<std::pair<double, double> > expansion_reaction ;
std::vector<std::pair<double, double> > expansion_stress ;
std::vector<double> apparent_extension ;

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
double E_agg = 58.9e9 ;
double E_paste = 12e9 ;
Matrix m0_agg(3,3) ;

double factor = 200 ;
MinimumAngle cri(M_PI/6.) ;
bool nothingToAdd = false ;
bool dlist = false ;
int count = 0 ;
double aggregateArea = 0;

GelBehaviour * gel = new GelBehaviour() ;

void step()
{
	featureTree->setMaxIterationsPerStep(2000) ;
	int nsteps = 1;
	for(size_t i = 0 ; i < nsteps ; i++)
	{

		bool go_on = featureTree->step() ;

		triangles = featureTree->getElements2D() ;

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
		
		int npoints = triangles[0]->getBoundingPoints().size() ;
		
		double area = 0 ;
		double avg_e_xx = 0;
		double avg_e_yy = 0;
		double avg_e_xy = 0;
		double avg_s_xx = 0;
		double avg_s_yy = 0;
		double avg_s_xy = 0;
		double e_xx = 0 ;
		double e_yy = 0 ;
		double ex_count = 0 ;
		double ey_count = 0 ;
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
					if(triangles[k]->getBoundingPoint(p).x > 0.0199)
					{
						e_xx+=x[triangles[k]->getBoundingPoint(p).id*2] ;
						ex_count++ ;
					}
					if(triangles[k]->getBoundingPoint(p).y > 0.0199)
					{
						e_yy+=x[triangles[k]->getBoundingPoint(p).id*2+1] ;
						ey_count++ ;
					}
				}
				area += triangles[k]->area() ;
				if(triangles[k]->getBehaviour()->type != VOID_BEHAVIOUR)
				{
					if(!triangles[k]->getBehaviour()->param.isNull() && triangles[k]->getBehaviour()->param[0][0] > E_max)
						E_max = triangles[k]->getBehaviour()->param[0][0] ;
					if(!triangles[k]->getBehaviour()->param.isNull() && triangles[k]->getBehaviour()->param[0][0] < E_min)
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
				for(size_t l = 0 ; l < npoints ;l++)
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
					for(size_t l = 0 ; l < npoints ;l++)
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
		
		double reactedArea = 0 ;
		double final_radius = sqrt(.90)*0.001 ;
		double delta_radius = 2e-5 ;/*(inclusions[inclusions.size()/2]->getRadius()-final_radius)/nsteps ;*/
		for(size_t m = 0 ; m < reactionRims.size() ; m++)
		{
			double new_in_radius = std::max(reactionRims[m]->getInRadius()-delta_radius, 0.) ;
			double new_area = M_PI*(reactionRims[m]->getRadius()*reactionRims[m]->getRadius()-new_in_radius*new_in_radius) ;
			double inc_area = M_PI*reactionRims[m]->getRadius()*reactionRims[m]->getRadius() ;
			std::cout << "..." << reactionRims[m]->getInRadius()<< " -> "<< new_in_radius << "..." << std::endl ;
			if(new_in_radius >= 1e-6)
			{
				reactedArea += (reactionRims[m]->getRadius()*reactionRims[m]->getRadius() - new_in_radius*new_in_radius)*M_PI ;
				reactionRims[m]->setInRadius(new_in_radius) ;
			}
			else
			{
				reactedArea += (reactionRims[m]->getRadius()*reactionRims[m]->getRadius() - reactionRims[m]->getInRadius()*reactionRims[m]->getInRadius())*M_PI ;
			}
			
			
		}
		
		std::string filename("rag_simple_out") ;
// 		filename.append(itoa(totit++, 10)) ;
		std::cout << filename << std::endl ;

		TriangleWriter writer(filename, featureTree) ;
		writer.getField(STRAIN_FIELD) ;
		writer.getField(REAL_STRESS_FIELD) ;
		writer.getField(TWFT_VON_MISES) ;
		writer.getField(TWFT_STIFFNESS) ;
		writer.write() ;

		
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
		
		std::cout << "apparent extension x" << e_xx/ex_count << std::endl ;
		std::cout << "apparent extension y" << e_yy/ey_count << std::endl ;
		//(1./epsilon11.x)*( stressMoyenne.x-stressMoyenne.y*modulePoisson);
		
		if (go_on)
		{
			expansion_reaction.push_back(std::make_pair(reactedArea/placed_area, avg_e_xx/area)) ;
			expansion_reaction.push_back(std::make_pair(reactedArea/placed_area, avg_e_yy/area)) ;
			expansion_stress.push_back(std::make_pair((avg_e_xx_nogel+avg_e_yy_nogel)/(2.*nogel_area), (avg_s_xx_nogel+avg_s_yy_nogel)/(2.*nogel_area))) ;
			apparent_extension.push_back(e_xx/ex_count) ;
			apparent_extension.push_back(e_yy/ey_count) ;
		}
		
		if (!go_on >= 10000)
			break ;

		for(size_t l = 0 ; l < expansion_reaction.size()/2 ; l++)
			std::cout << expansion_reaction[l*2].first << "   " 
			<< expansion_reaction[l*2].second << "   " 
			<< expansion_reaction[l*2+1].second << "   " 
			<< expansion_stress[l].first << "   " 
			<< expansion_stress[l].second << "   " 
			<< apparent_extension[l*2]  << "   " 
			<< apparent_extension[l*2+1]  << "   " 
			<< std::endl ;
	}
	
	for(size_t l = 0 ;l < expansion_reaction.size()/2 ; l++)
			std::cout << expansion_reaction[l*2].first << "   " 
			<< expansion_reaction[l*2].second << "   " 
			<< expansion_reaction[l*2+1].second << "   " 
			<< expansion_stress[l].first << "   " 
			<< expansion_stress[l].second << "   " 
			<< apparent_extension[l*2]  << "   " 
			<< apparent_extension[l*2+1]  << "   " 
			<< std::endl ;
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
			double radius = 0.00005 ;
			
			Point pos((2.*rand()/RAND_MAX-1.),(2.*rand()/RAND_MAX-1.)) ;
			pos /= pos.norm() ;
			pos *= (2.*rand()/RAND_MAX-1.)*(incs[i]->getRadius() - 0.00003) ;
			Point center = incs[i]->getCenter()+pos ; 
			
			bool alone  = true ;
			
			for(size_t k = 0 ; k < ret.size() ; k++ )
			{
				if (squareDist(center, ret[k].first->Circle::getCenter()) < 32.*(radius+radius)*(radius+radius))
				{
					alone = false ;
					break ;
				}
			}
			if (alone)
			{
				ExpansiveZone * z = new ExpansiveZone(incs[i], radius, center.x, center.y, gel) ;
				ret.push_back(std::make_pair(z, incs[i])) ;
			}
		}
	}

	for(size_t i = 0 ; i < ret.size() ; i++)
	{
		ret[i].first->setRadius(0.000001) ;
		F.addFeature(ret[i].second, ret[i].first) ;
	}
	std::cout << "initial Reacted Area = " << M_PI*0.000001*0.000001*ret.size() << " in "<< ret.size() << " zones"<< std::endl ;
	return ret ;	
}

std::pair<std::vector<Inclusion * >, std::vector<Pore * > > generateInclusionsAndPores(size_t n, double fraction, double E, double nu, Feature * father, FeatureTree * F)
{
// 	srandom(time(nullptr)) ;
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
		(*ret.first.rbegin())->setBehaviour(new WeibullDistributedStiffness(E, nu, SPACE_TWO_DIMENSIONAL,-8.*1000000, 1000000)) ;
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
	
// 	TriElement father(LINEAR) ;
// 	Function x = father.getXTransform() ;
// 	Function y = father.getYTransform() ;
// 	double angle = 0 ;
// 	double rotatedSingularityX = .3*cos ( angle ) + .3*sin ( angle ) ;
// 	double rotatedSingularityY = -.3*sin ( angle ) + .3*cos ( angle ) ;
// 	Function rotatedX = x*cos ( angle ) + y*sin ( angle ) ;
// 	Function rotatedY = x*sin ( -angle ) + y*cos ( angle ) ;
// 	Function x_ = x - .3 ;
// 	Function y_ = y - .3 ;
// 	Function theta = f_atan2 ( rotatedY-rotatedSingularityY, rotatedX-rotatedSingularityX );
// 	Function r = f_sqrt ( ( x_^2 ) + ( y_^2 ) );
// 	
// 	Function f0 = f_sqrt ( r ) *f_sin ( theta/2 );
// 	Function f1 = f_sqrt ( r ) *f_cos ( theta/2 );
// 	Function f2 = f_sqrt ( r ) *f_sin ( theta/2 ) *f_cos ( theta );
// 	Function f3 = f_sqrt ( r ) *f_cos ( theta/2 ) *f_cos ( theta );
// 	f3.compile() ;
// 	VirtualMachine vm ;
// 	for(size_t i = 0 ; i < 100000 ; i++)
// 		/*std::cout <<*/ vm.eval(f3, .333, .333, .333) /*<< std::endl*/ ;
// 	return 0 ;
	
	percent = atof(argv[1]) ;

	double E_csh = 31e9 ;
	double nu_csh = .28 ;
	double nu_incompressible = .499997 ;
	
	double E = percent*E_csh ;
	double nu = nu_incompressible ;
	
	Matrix m0(3,3) ;
	m0[0][0] = E/(1.-nu*nu) ; m0[0][1] =E/(1.-nu*nu)*nu ; m0[0][2] = 0 ;
	m0[1][0] = E/(1.-nu*nu)*nu ; m0[1][1] = E/(1.-nu*nu) ; m0[1][2] = 0 ; 
	m0[2][0] = 0 ; m0[2][1] = 0 ; m0[2][2] = E/(1.-nu*nu)*(1.-nu)/2. ; 
	
	m0_agg[0][0] = E_agg/(1-nu*nu) ; m0_agg[0][1] =E_agg/(1-nu*nu)*nu ; m0_agg[0][2] = 0 ;
	m0_agg[1][0] = E_agg/(1-nu*nu)*nu ; m0_agg[1][1] = E_agg/(1-nu*nu) ; m0_agg[1][2] = 0 ; 
	m0_agg[2][0] = 0 ; m0_agg[2][1] = 0 ; m0_agg[2][2] = E_agg/(1-nu*nu)*(1.-nu)/2. ; 
	
	Matrix m0_paste(3,3) ;
	m0_paste[0][0] = E_paste/(1-nu*nu) ; m0_paste[0][1] =E_paste/(1-nu*nu)*nu ; m0_paste[0][2] = 0 ;
	m0_paste[1][0] = E_paste/(1-nu*nu)*nu ; m0_paste[1][1] = E_paste/(1-nu*nu) ; m0_paste[1][2] = 0 ; 
	m0_paste[2][0] = 0 ; m0_paste[2][1] = 0 ; m0_paste[2][2] = E_paste/(1-nu*nu)*(1.-nu)/2. ; 

	Sample sample(nullptr, 0.04, 0.04, 0, 0) ;
	
	FeatureTree F(&sample) ;
	featureTree = &F ;


	double itzSize = 4e-5;
	int inclusionNumber = 5000 ; // 10 100 500 1000 2000 4000

	double masseInitiale = .00000743;
	double densite = 1.;
// 	inclusions = GranuloBolome(masseInitiale, densite, BOLOME_A)(.008, .0001, inclusionNumber, itzSize);
//	inclusions = GranuloBolome(4.79263e-07/.47, 1, BOLOME_D)(.0025, .0001, inclusionNumber, itzSize);
	inclusions = PSDGenerator::get2DInclusions(0.0025, 4.79263e-07/.47, new PSDBolomeD(), PSDEndCriteria(-1, 0.001, inclusionNumber)) ;

	std::vector<Feature *> feats ;
	for(size_t i = 0; i < inclusions.size() ; i++)
		feats.push_back(inclusions[i]) ;


	inclusions.clear() ;
	for(size_t i = 0; i < feats.size() ; i++)
		inclusions.push_back(static_cast<Inclusion *>(feats[i])) ;

	
	int nAgg = 0 ;
	feats=placement(sample.getPrimitive(), feats, &nAgg, 0, 6400);
	double volume = 0 ;
	for(size_t i = 0 ; i < feats.size() ; i++)
		volume += feats[i]->area() ;
	if(!feats.empty())
		std::cout << "n = " << feats.size() << ", largest r = " << feats.front()->getRadius() 
		<< ", smallest r =" << feats.back()->getRadius() 
		<< ", filling = " << volume/sample.area()*100.<< "%"<< std::endl ; 

	sample.setBehaviour(new WeibullDistributedStiffness(E_paste, nu, SPACE_TWO_DIMENSIONAL, -8.*12500000, 12500000)) ;
// 	sample.setBehaviour(new Stiffness(m0_paste)) ;
	Vector a(3) ;
	a[0] = .5 ;
	a[1] = .5 ;
// 	ExpansiveRing er(&sample, .01, .0099, 0,0,  m0_paste, a) ;
// 	F.addFeature(&sample,&er) ;
	double reactedArea = 0 ;
	std::random_shuffle(inclusions.begin(), inclusions.end());
// 	for(size_t i = 0 ; i < inclusions.size()/10 ; i++)
// 	{
// 		inclusions[i]->setRadius(inclusions[i]->getRadius()-itzSize) ;
// 		inclusions[i]->setBehaviour(new WeibullDistributedStiffness(m0_agg,57000000)) ;
// // 		inclusions[i]->setBehaviour(new Stiffness(m0_agg)) ;
// 		double rPlus = inclusions[i]->getRadius()+itzSize/120000. ;
// 		double rMinus = inclusions[i]->getRadius()-itzSize/120000. ;
// 		double cx = inclusions[i]->getCenter().x ;
// 		double cy = inclusions[i]->getCenter().y ;
// 		ExpansiveRing * er = new ExpansiveRing(&sample, rPlus, rMinus, cx, cy,  m0, a) ;
// 		reactionRims.push_back(er) ;
// // 		inclusions[i]->setBehaviour(new StiffnessWithVariableImposedDeformation(m0_agg, a)) ;
// 		F.addFeature(&sample,inclusions[i]) ;
// 		F.addFeature(inclusions[i],reactionRims.back()) ;
// 		placed_area += inclusions[i]->area() ;
// 		reactedArea += (reactionRims.back()->getRadius()*reactionRims[i]->getRadius() - reactionRims.back()->getInRadius()*reactionRims.back()->getInRadius())*M_PI ;
// 	}
// 	for(size_t i = inclusions.size()/10 ; i < inclusions.size() ; i++)
// 	{
// 		inclusions[i]->setRadius(inclusions[i]->getRadius()-itzSize) ;
// 		inclusions[i]->setBehaviour(new WeibullDistributedStiffness(m0_agg,57000000)) ;
// // 		inclusions[i]->setBehaviour(new Stiffness(m0_agg)) ;
// // 		inclusions[i]->setBehaviour(new StiffnessWithVariableImposedDeformation(m0_agg, a)) ;
// 		F.addFeature(&sample,inclusions[i]) ;
// 		placed_area += inclusions[i]->area() ;
// 	}
	Inclusion * inc = new Inclusion(0.010, 0, 0) ;
// 	inc->setBehaviour(new Stiffness(m0_agg)) ;
	inc->setBehaviour(new WeibullDistributedStiffness(E_agg, nu, SPACE_TWO_DIMENSIONAL,-8.*57000000,57000000)) ;
	double rPlus = 0.010+itzSize ;
	double rMinus = 0.009 ;
	ExpansiveRing * er = new ExpansiveRing(&sample, rPlus, rMinus, 0, 0,  m0, a) ;
	reactionRims.push_back(er) ;
	F.addFeature(&sample,inc) ;
	F.addFeature(inc,er) ;
	
	if(!inclusions.empty())
	{
		std::cout << "largest inclusion with r = " << (*inclusions.begin())->getRadius() << std::endl ;
		std::cout << "smallest inclusion with r = " << (*inclusions.rbegin())->getRadius() << std::endl ;
		std::cout << "placed area = " <<  placed_area << std::endl ;
		std::cout << "initial reacted area = " <<  reactedArea << std::endl ;
	}
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, BOTTOM_LEFT)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM_LEFT)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, RIGHT)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_XI, RIGHT, -5e6)) ;
	F.setSamplingNumber(256) ;

	F.setOrder(LINEAR) ;

// 	
	step() ;

// 	delete dt ;
	
	return 0 ;
}
