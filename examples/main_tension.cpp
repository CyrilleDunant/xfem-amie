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
#include "../physics/fracturecriteria/mcft.h"
#include "../physics/fracturecriteria/fractionmcft.h"
#include "../physics/fracturecriteria/vonmises.h"
#include "../physics/fracturecriteria/nonlocalvonmises.h"
#include "../physics/fracturecriteria/druckerprager.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/fracturecriteria/boundedvonmises.h"
#include "../physics/damagemodels/plasticstrain.h"
#include "../physics/stiffness.h"
#include "../physics/materials/aggregate_behaviour.cpp"
#include "../physics/stiffness_and_fracture.h"
#include "../physics/fraction_stiffness_and_fracture.h"
#include "../physics/void_form.h"
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
#include "../utilities/optimizer.h"
#include "../utilities/itoa.h"
#include "../utilities/random.h"
#include "../utilities/writer/triangle_writer.h"
#include "../physics/materials/concrete_behaviour.h"
#include "../physics/materials/paste_behaviour.h"
#include "../physics/orthotropicstiffness.h"


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

double timepos = 0.05e-07 ;

double effectiveRadius = 0.13506098 ; //.5*.165*sqrt(M_PI*.25) ;//
double rebarDiametre = 0.0254*sqrt(M_PI*.25) ;

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

double sampleLength = 5.5 ;
double sampleHeight = 1.2 ;
double supportLever = 2.5 ; 
double supportMidPointToEndClearance = 0.25 ;
double platewidth = 0.15 ;
double plateHeight = 0.051 ;
double rebarEndCover = 0.047 ;

std::vector<DelaunayTriangle *> tris__ ;
double apriori_command = 0 ;
std::pair<std::vector<Inclusion * >, std::vector<Pore * > > i_et_p ;

std::vector<std::pair<ExpansiveZone *, Inclusion *> > zones ;

std::vector<std::pair<double, double> > expansion_reaction ;
std::vector<std::pair<double, double> > expansion_stress ;
std::vector<std::pair<double, double> > load_displacement ;
std::vector< double > loads ;
std::vector< double > displacements ;
std::vector< double > loadsx ;
std::vector< double > displacementsx ;
std::vector< double > damages ;
std::vector< double > cmod ;
std::vector< double > idx ;

Vector fracCrit(0) ;

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

// BoundingBoxNearestNodeDefinedBoundaryCondition * loadr = new BoundingBoxNearestNodeDefinedBoundaryCondition(SET_ALONG_ETA, TOP, Point(+0.04445, .0672*.5), -0.000005) ;
// BoundingBoxNearestNodeDefinedBoundaryCondition * loadt = new BoundingBoxNearestNodeDefinedBoundaryCondition(SET_ALONG_ETA, TOP, Point(-0.04445, .0672*.5), -0.000005) ;


// BoundingBoxAndRestrictionDefinedBoundaryCondition * load = new BoundingBoxAndRestrictionDefinedBoundaryCondition(SET_STRESS_ETA, TOP, -.15, .15, -10, 10, -10.) ;
BoundingBoxDefinedBoundaryCondition * loadt = new BoundingBoxDefinedBoundaryCondition(SET_ALONG_XI, RIGHT,0) ;
// BoundingBoxDefinedBoundaryCondition * loadr = new BoundingBoxDefinedBoundaryCondition(SET_STRESS_XI, RIGHT,1e4) ;
BoundingBoxDefinedBoundaryCondition * loadr = new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA, TOP,0) ;
// BoundingBoxNearestNodeDefinedBoundaryCondition * loadr = new BoundingBoxNearestNodeDefinedBoundaryCondition(SET_FORCE_XI, RIGHT, Point(1.3*.5+.225, 0)) ;
// BoundingBoxDefinedBoundaryCondition * loadl = new BoundingBoxDefinedBoundaryCondition(SET_STRESS_XI, LEFT,0) ;
// BoundingBoxNearestNodeDefinedBoundaryCondition * load = new BoundingBoxNearestNodeDefinedBoundaryCondition(SET_FORCE_ETA, TOP, Point(0., 1.2), 0) ;
double factor = 25 ;
MinimumAngle cri(M_PI/6.) ;
bool nothingToAdd = false ;
bool dlist = false ;
int count = 0 ;
double aggregateArea = 0;

MultiTriangleWriter writer( "triangles_head", "triangles_layers", nullptr ) ;
MultiTriangleWriter writerc( "triangles_converged_head", "triangles_converged_layers", nullptr ) ;

void step(size_t nsteps)
{
	
	size_t nit = 2 ;
	size_t ntries = 5;
	size_t dsteps = 60 ;
	size_t tries = 0 ;
	size_t dit = 0 ;
	int totit = 0 ;
	for(size_t v = 0 ; v < nsteps ; v++)
	{
 		tries = 0 ;

		tries++ ;

		bool go_on = featureTree->step() ;
		double appliedForce = loadr->getData()*effectiveRadius*2.*rebarDiametre;
		if(go_on)
		{
// 			loadr->setData(sin(double(count)/40.)*5e-5) ;
// 			if(count < 80)
// 				loadr->setData(loadr->getData()-1e-5) ;
// 			else
// 				loadr->setData(loadr->getData()+1e-7) ;
			count++ ;
// 			loadt->setData(loadt->getData()+6e-6) ;
			loadr->setData(loadr->getData()-6e-5) ;
			
// 			loadt->setData(0) ;
		}
		
		triangles = featureTree->getActiveElements2D() ;
		x.resize( featureTree->getDisplacements().size() ) ;
		x = featureTree->getDisplacements() ;



		sigma.resize( triangles.size()*triangles[0]->getBoundingPoints().size()*3 ) ;
		epsilon.resize( triangles.size()*triangles[0]->getBoundingPoints().size()*3 ) ;
		std::pair<Vector, Vector > sigma_epsilon = featureTree->getStressAndStrainInAllLayers() ;
		sigma.resize( sigma_epsilon.first.size() ) ;
		sigma = sigma_epsilon.first ;
		epsilon.resize( sigma_epsilon.second.size() ) ;
		epsilon = sigma_epsilon.second ;
		sigma11.resize( sigma.size() / 3 ) ;
		sigma22.resize( sigma.size() / 3 ) ;
		sigma12.resize( sigma.size() / 3 ) ;
		epsilon11.resize( sigma.size() / 3 ) ;
		fracCrit.resize( sigma.size() / 3 ) ;
		epsilon22.resize( sigma.size() / 3 ) ;
		epsilon12.resize( sigma.size() / 3 ) ;
		vonMises.resize( sigma.size() / 3 ) ;
		angle.resize( sigma.size() / 3 ) ;
		
		std::cerr << "unknowns :" << x.size() << std::endl ;
		
		
		int npoints = triangles[0]->getBoundingPoints().size() ;
		
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
		double leftDisp = 0;
		double rightDisp = 0 ;
		int tsize = 0 ;
		std::vector< std::pair<double, double> > pos_strain ;
		for(size_t k = 0 ; k < triangles.size() ; k++)
		{
			if(!triangles[k]->getBehaviour())
				continue ;
// 			if(triangles[k]->getBehaviour()->fractured())
// 				continue ;
			if(triangles[k]->getBehaviour()->type == VOID_BEHAVIOUR)
				continue ;
			bool in = false ;

			if(std::abs(triangles[k]->getCenter().y) < 0.0125)
			{
				Point where(triangles[k]->getCenter().x, 0) ;
				Vector stra(0., 3) ;
				if(triangles[k]->in(where))
				{
					triangles[k]->getState().getField(STRAIN_FIELD, where, stra, false) ;
					pos_strain.push_back(std::make_pair( triangles[k]->getCenter().x , stra[0] )) ;
				}
			}
			if(dist(triangles[k]->getCenter(), Point(0.0025, -.1*.5+.03*.5)) < 1.5*triangles[k]->getRadius())
			{
				Vector tmpd(2) ;
				triangles[k]->getState().getField(DISPLACEMENT_FIELD, triangles[k]->getCenter(),tmpd, false) ;
				leftDisp = tmpd[0] ;
			}
			if(dist(triangles[k]->getCenter(), Point(-0.0025, -.1*.5+.03*.5)) < 1.5*triangles[k]->getRadius())
			{
				Vector tmpd(2) ;
				triangles[k]->getState().getField(DISPLACEMENT_FIELD,triangles[k]->getCenter(),tmpd, false) ;
				rightDisp = tmpd[0] ;
			}
			
			for(size_t p = 0 ;p < triangles[k]->getBoundingPoints().size() ; p++)
			{
				if(triangles[k]->getBehaviour()->type != VOID_BEHAVIOUR)
				{
					if(x[triangles[k]->getBoundingPoint(p).id*2] > x_max)
						x_max = x[triangles[k]->getBoundingPoint(p).id*2];
					if(x[triangles[k]->getBoundingPoint(p).id*2] < x_min)
						x_min = x[triangles[k]->getBoundingPoint(p).id*2];
					if(x[triangles[k]->getBoundingPoint(p).id*2+1] > y_max)
						y_max = x[triangles[k]->getBoundingPoint(p).id*2+1];
					if(x[triangles[k]->getBoundingPoint(p).id*2+1] < y_min)
						y_min = x[triangles[k]->getBoundingPoint(p).id*2+1];
					if(triangles[k]->getBoundingPoint(p).y > (sampleHeight+plateHeight*2)*.5*.9999)
					{
						e_xx+=x[triangles[k]->getBoundingPoint(p).id*2+1] ;
						ex_count++ ;
					}
						
				}
			}
			area += triangles[k]->area() ;
			if(triangles[k]->getBehaviour()->type != VOID_BEHAVIOUR)
			{
				tsize++ ;
				if(triangles[k]->getBehaviour()->getTensor(triangles[k]->getCenter())[0][0] > E_max)
					E_max = triangles[k]->getBehaviour()->getTensor(triangles[k]->getCenter())[0][0] ;
				if(triangles[k]->getBehaviour()->getTensor(triangles[k]->getCenter())[0][0] < E_min)
					E_min = triangles[k]->getBehaviour()->getTensor(triangles[k]->getCenter())[0][0] ;
			}
			
//			Vector se = triangles[k]->getState().getPrincipalStresses(triangles[k]->getBoundingPoint(0)) ;
			sigma11[k*npoints] = sigma[k*npoints*3];//se[0] ;
			sigma22[k*npoints] = sigma[k*npoints*3+1];//se[1] ;
			sigma12[k*npoints] = sigma[k*npoints*3+2];
//			se = triangles[k]->getState().getPrincipalStresses(triangles[k]->getBoundingPoint(1)) ;
			sigma11[k*npoints+1] = sigma[k*npoints*3+3];//se[0] ;
			sigma22[k*npoints+1] = sigma[k*npoints*3+4];//se[1] ;
			sigma12[k*npoints+1] = sigma[k*npoints*3+5];
//			se = triangles[k]->getState().getPrincipalStresses(triangles[k]->getBoundingPoint(2)) ;
			sigma11[k*npoints+2] = sigma[k*npoints*3+6];//se[0] ;
			sigma22[k*npoints+2] = sigma[k*npoints*3+7];//se[1] ;
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
//			Vector pe =triangles[k]->getState().getPrincipalStrains(triangles[k]->getBoundingPoint(0)) ;
			epsilon11[k*npoints] = epsilon[k*npoints*3];//pe[0] ;
			epsilon22[k*npoints] = epsilon[k*npoints*3+1];//pe[1] ; 
			epsilon12[k*npoints] = epsilon[k*npoints*3+2];
//			pe =triangles[k]->getState().getPrincipalStrains(triangles[k]->getBoundingPoint(1)) ;
			epsilon11[k*npoints+1] = epsilon[k*npoints*3+3];//pe[0] ;
			epsilon22[k*npoints+1] = epsilon[k*npoints*3+4];//pe[1] ;
			epsilon12[k*npoints+1] = epsilon[k*npoints*3+5];
//			pe =triangles[k]->getState().getPrincipalStrains(triangles[k]->getBoundingPoint(2)) ;
			epsilon11[k*npoints+2] = epsilon[k*npoints*3+6];//pe[0] ;
			epsilon22[k*npoints+2] = epsilon[k*npoints*3+7];//pe[1] ;
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
				
				if(triangles[k]->getBehaviour()->getFractureCriterion())
				{
					fracCrit[k*triangles[k]->getBoundingPoints().size()+l] = triangles[k]->getBehaviour()->getFractureCriterion()->getNonLocalScoreAtState() ;
				}
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
		
		if(go_on)
		{
			displacements.push_back(1e6*avg_e_yy/area);
			loads.push_back((avg_s_yy/area)/1e6);
			displacementsx.push_back(1e6*avg_e_xx/area);
			loadsx.push_back((avg_s_xx/area)/1e6);
			damages.push_back(featureTree->averageDamage);
			cmod.push_back(rightDisp-leftDisp);
			idx.push_back(loadr->getData());
		}
		if(v%5 == 0 || true)
		{

			std::cout << std::endl ;
			if(!loads.empty())
				std::cout << "load : " <<  loads.back() << std::endl ;
			std::cout << "max value : " << x_max*1000 << std::endl ;
			std::cout << "min value : " << x_min*1000 << std::endl ;
			std::cout << "max sigma11 : " << sigma11.max()/1e6 << std::endl ;
			std::cout << "min sigma11 : " << sigma11.min()/1e6 << std::endl ;
			std::cout << "max sigma12 : " << sigma12.max()/1e6 << std::endl ;
			std::cout << "min sigma12 : " << sigma12.min()/1e6 << std::endl ;
			std::cout << "max sigma22 : " << sigma22.max()/1e6 << std::endl ;
			std::cout << "min sigma22 : " << sigma22.min()/1e6 << std::endl ;
			
			std::cout << "max epsilon11 : " << epsilon11.max()*1e6 << std::endl ;
			std::cout << "min epsilon11 : " << epsilon11.min()*1e6 << std::endl ;
			std::cout << "max epsilon12 : " << epsilon12.max()*1e6 << std::endl ;
			std::cout << "min epsilon12 : " << epsilon12.min()*1e6 << std::endl ;
			std::cout << "max epsilon22 : " << epsilon22.max()*1e6 << std::endl ;
			std::cout << "min epsilon22 : " << epsilon22.min()*1e6 << std::endl ;
			
			std::cout << "max von Mises : " << vonMises.max()/1e6 << std::endl ;
			std::cout << "min von Mises : " << vonMises.min()/1e6 << std::endl ;
			
			std::cout << "average sigma11 : " << (avg_s_xx/area)/1e6 << std::endl ;
			std::cout << "average sigma22 : " << (avg_s_yy/area)/1e6 << std::endl ;
			std::cout << "average sigma12 : " << (avg_s_xy/area)/1e6 << std::endl ;
			std::cout << "average epsilon11 : " << 1e6*avg_e_xx/area<< std::endl ;
			std::cout << "average epsilon22 : " << 1e6*avg_e_yy/area << std::endl ;
			std::cout << "average epsilon12 : " << 1e6*avg_e_xy/area << std::endl ;

		}

		
		if(go_on)
			std::cout << appliedForce/1000. << std::endl ;
		
		std::fstream ldfile  ;
		ldfile.open("ldn", std::ios::out) ;
		for(int j = 0 ; j < loads.size() ; j++)
		{
			ldfile << displacements[j] << "   " << loads[j] << "   " <<  displacementsx[j] << "   " << loadsx[j] << "  " << cmod[j] <<  "  "<<  idx[j]<< "\n" ;
		}
		ldfile.close();
		
		std::fstream strfile  ;
		std::sort(pos_strain.begin(), pos_strain.end()) ;
		ldfile.open("str", std::ios::out) ;
		for(int j = 0 ; j < pos_strain.size() ; j++)
		{
			ldfile << pos_strain[j].first << "   " << pos_strain[j].second << "\n" ;
		}
		strfile.close();
		
		if ( true )
		{
			writer.reset( featureTree ) ;
			writer.getField( TWFT_PRINCIPAL_STRESS ) ;
			writer.getField( TWFT_PRINCIPAL_STRAIN ) ;
			writer.getField( TWFT_CRITERION ) ;
			writer.getField( TWFT_PRINCIPAL_ANGLE ) ;
			writer.getField( TWFT_STIFFNESS_X ) ;
			writer.getField( TWFT_STIFFNESS_Y ) ;
			writer.append() ;
		}
		if(go_on)
		{
			writerc.reset( featureTree ) ;
			writerc.getField( TWFT_PRINCIPAL_STRESS ) ;
			writerc.getField( TWFT_PRINCIPAL_STRAIN ) ;
			writerc.getField( TWFT_CRITERION ) ;
			writerc.getField( TWFT_PRINCIPAL_ANGLE ) ;
			writerc.getField( TWFT_STIFFNESS_X ) ;
			writerc.getField( TWFT_STIFFNESS_Y ) ;
			writerc.append() ;
		}
		//(1./epsilon11.x)*( stressMoyenne.x-stressMoyenne.y*modulePoisson);
	}

}




int main(int argc, char *argv[])
{
	OrthotropicStiffness orth(1, 2, 1,  .2, 0) ;
// 	std::vector<std::pair<std::string, double> >  vars ;
// 	vars.push_back(std::make_pair("L", 300)); Function L("L") ;
// 	vars.push_back(std::make_pair("f", 1.)); Function f("f") ;
// 	vars.push_back(std::make_pair("m", 2)); Function m("m") ;
// 	vars.push_back(std::make_pair("E", 200)); Function E("E") ;
// 	
// 	double I=70.*pow(1.02,3) ;
// 	Function EI = E*I*1e6 ;
// 	Function stif = (EI*2)/(L^3) ;
// 	Function mass = m+1.080 ; 
// 	std::vector<std::pair<double, double> >   bounds ;
// 	bounds.push_back(std::make_pair(450, 450));
// 	bounds.push_back(std::make_pair(1.2, 1.2));
// 	bounds.push_back(std::make_pair(0., .5));
// 	bounds.push_back(std::make_pair(200, 200));
// 	Function objectiveFunction = f_abs(f/0.64-(1./(2.*M_PI))*f_sqrt(stif/mass)) ;
// 	GeneticAlgorithmOptimizer ga(vars, bounds, objectiveFunction) ;
// 	ga.optimize(1e-8, 100, 5000, .1, .65) ;
// 	std::vector<std::pair<std::string, double> > vals = ga.getValues() ;
// 	for(size_t i = 0 ; i < vals.size() ; i++)
// 		std::cout << vals[i].first << " = " << vals[i].second << std::endl ;
// 	std::cout << "effective f = " << VirtualMachine().eval((0.64/(2.*M_PI))*f_sqrt(stif/mass),ga.getValues() ) << std::endl;
// 	exit(0) ;
// 	
// 	
// 	
	
	
	double compressionCrit = -32.6e6 ; 
	double steelfraction = 0.5*rebarDiametre/effectiveRadius ;
	std::cout << "steel fraction = " << steelfraction << std::endl ;
	double mradius = .05 ; // .010 ;//
// 	double mradius = .25 ;
	double length =  0.3048 ; //1.300*.5;//
	double E_steel = 193e9/1.7 ;
	double nu_steel = 0. ; 
	
	double nu = 0.2 ;
	double E_paste = 30e9 ;
	
	Matrix m0_paste= Material::cauchyGreen(std::make_pair(E_paste,nu), true,SPACE_TWO_DIMENSIONAL) ;
	Matrix m0_steel = Material::cauchyGreen(std::make_pair(E_steel,nu_steel), true,SPACE_TWO_DIMENSIONAL) ;

	Sample box(1.300*.5+.225, effectiveRadius, (1.300*.5+.225)*.5, effectiveRadius*0.5) ;
	box.setBehaviour(new VoidForm()) ;  
	Sample sample(1.300*.5, effectiveRadius-rebarDiametre*.5, 1.300*.25, rebarDiametre*.5+(effectiveRadius-rebarDiametre*.5)*0.5) ;
// 	Sample samplef(length, effectiveRadius, 1.300*.25, (effectiveRadius)*0.5) ;
// 	Sample samplef(.1200, .1200, 0., 0.) ;
// 	Sample samplef(0.3048, .0672, 0., 0.) ;
	Sample samplef(0.6, 0.6, 0., 0.) ;
	Sample notch(.005, .03, 0., -.1*.5+.03*.5) ;
	
	Sample toprightvoid(0.225, effectiveRadius-rebarDiametre*.5, 1.300*.5+0.225*0.5, rebarDiametre*.5+(effectiveRadius-rebarDiametre*.5)*0.5) ;     
	toprightvoid.setBehaviour(new VoidForm()) ;  

	Sample rebarright(0.225, rebarDiametre*.5, 1.300*.5+.225*.5, rebarDiametre*0.25) ; 
	rebarright.setBehaviour(new VoidForm());
	
// 	Sample rebarright((1.300)*.5+0.225, rebarDiametre*.5, (1.300)*.25+0.225*.5, rebarDiametre*0.25) ; 
// 	rebarright.setBehaviour(new StiffnessAndFracture(m0_steel, new VonMises(500e6, MIRROR_XY)));
// 	rebarright.getBehaviour()->getFractureCriterion()->setMaterialCharacteristicRadius(mradius);
	
	Sample rebarinternal((1.300)*.5+0.225, rebarDiametre*.5, (1.300)*.25+0.225*.5, rebarDiametre*0.25) ; 
	rebarinternal.setBehaviour(new StiffnessAndFracture(m0_steel*.5, new VonMises(500e6, MIRROR_XY)));
	rebarinternal.getBehaviour()->getFractureCriterion()->setMaterialCharacteristicRadius(mradius);
// 	rebarinternal.getBehaviour()->getFractureCriterion()->setNeighbourhoodRadius(nradius);
// 	



// 	sample.setBehaviour(new VoidForm()) ;  
		
// 	samplef.setBehaviour(new ConcreteBehaviour()) ;
// 	dynamic_cast<ConcreteBehaviour *>(samplef.getBehaviour())->materialRadius = mradius ;
// 	dynamic_cast<ConcreteBehaviour *>(samplef.getBehaviour())->neighbourhoodRadius = nradius ;
// // 	dynamic_cast<ConcreteBehaviour *>(sample.getBehaviour() )->variability = 0.03 ;
	
// 	samplef.setBehaviour(new ConcreteBehaviour(E_paste, nu, compressionCrit,PLANE_STRESS)) ;
// 	dynamic_cast<ConcreteBehaviour *>(samplef.getBehaviour())->materialRadius = mradius ;
	
// 		samplef.setBehaviour(new ConcreteBehaviour()) ;
		notch.setBehaviour(new VoidForm());
	
// 	samplef.setBehaviour(new Stiffness(Material::cauchyGreen(std::make_pair(E_paste,nu), true,SPACE_TWO_DIMENSIONAL, PLANE_STRESS))) ;
// 		samplef.setBehaviour(new StiffnessAndFracture(Material::cauchyGreen(std::make_pair(E_paste,nu), true,SPACE_TWO_DIMENSIONAL, PLANE_STRESS) ,new NonLocalVonMises(20e6, E_paste, mradius), new NullDamage())) ;
// 	samplef.setBehaviour(new PasteBehaviour()) ;
// 	dynamic_cast<PasteBehaviour *>(samplef.getBehaviour())->materialRadius = mradius ;
	
// 	dynamic_cast<ConcreteBehaviour *>(samplef.getBehaviour())->neighbourhoodRadius = nradius ;
// // 	dynamic_cast<ConcreteBehaviour *>(samplef.getBehaviour() )->variability = 0.03 ;
// // 	
	
		
	FeatureTree F(&samplef) ;
// 	F.addFeature(&samplef,&notch);
	featureTree = &F ;
	
// 	Inclusion inc(.06, -samplef.width()*.5, samplef.width()*.5) ;
// 	F.addFeature(&samplef, &inc);
// 	PlasticStrain * incdamagemodel = new PlasticStrain() ;
// 	incdamagemodel->factor = .90 ;
// 	inc.setBehaviour( new StiffnessAndFracture(Material::cauchyGreen(std::make_pair(E_paste,nu), true,SPACE_TWO_DIMENSIONAL, PLANE_STRAIN) , new DruckerPrager(-12.315e6, 12.315e6,E_paste,0.1 , mradius),incdamagemodel));
// 	inc.isVirtualFeature = true ;
// // 	
// 	Inclusion inc0(.06, samplef.width()*.5,  -samplef.width()*.5) ;
// 	F.addFeature(&samplef, &inc0);
// 	PlasticStrain * inc0damagemodel = new PlasticStrain() ;
// 	inc0damagemodel->factor = .90 ;
// 	inc0.setBehaviour( new StiffnessAndFracture(Material::cauchyGreen(std::make_pair(E_paste,nu), true,SPACE_TWO_DIMENSIONAL, PLANE_STRAIN) , new DruckerPrager(-12.315e6, 12.315e6,E_paste,0.1 , mradius),inc0damagemodel));
// 	inc0.isVirtualFeature = true ;
// 	
// 	Inclusion inc2(.06, -samplef.width()*.25, samplef.width()*.25) ;
// 	F.addFeature(&samplef, &inc2);
// 	PlasticStrain * inc2damagemodel = new PlasticStrain() ;
// 	inc2damagemodel->factor = .90 ;
// 	inc2.setBehaviour( new StiffnessAndFracture(Material::cauchyGreen(std::make_pair(E_paste,nu), true,SPACE_TWO_DIMENSIONAL, PLANE_STRAIN) , new DruckerPrager(-12.315e6, 12.315e6,E_paste,0.1 , mradius),inc2damagemodel));
// 	inc2.isVirtualFeature = true ;
// // 	
// 	Inclusion inc3(.06, samplef.width()*.25,  -samplef.width()*.25) ;
// 	F.addFeature(&samplef, &inc3);
// 	PlasticStrain * inc3damagemodel = new PlasticStrain() ;
// 	inc3damagemodel->factor = .90 ;
// 	inc3.setBehaviour( new StiffnessAndFracture(Material::cauchyGreen(std::make_pair(E_paste,nu), true,SPACE_TWO_DIMENSIONAL, PLANE_STRAIN) , new DruckerPrager(-12.315e6, 12.315e6,E_paste,0.1 , mradius),inc3damagemodel));
// 	inc3.isVirtualFeature = true ;
// 	
// 	Inclusion inc1(.06, 0., 0.) ;
// 	F.addFeature(&samplef, &inc1);
// 	DruckerPrager * inc1criterion = new DruckerPrager(-20e6*.9, -20e6*.9,E_paste,0.1 , mradius) ; 
// 	PlasticStrain * inc1damagemodel = new PlasticStrain() ;
// // 	inc1criterion->cap = .99 ;
// // 	inc1damagemodel->factor = .90 ;
// 	StiffnessAndFracture * inc1Behaviour = new StiffnessAndFracture(Material::cauchyGreen(std::make_pair(E_paste,nu), true,SPACE_TWO_DIMENSIONAL, PLANE_STRAIN) , inc1criterion , inc1damagemodel) ; 
// 	inc1.setBehaviour( inc1Behaviour );
// 	inc1.setBehaviourSource(&samplef);
// 	inc1.isVirtualFeature = true ;
	
	double verticaloffset = 0 ; //1./14.*samplef.height()
	TriangularInclusion t0(Point(samplef.width()*.5,-samplef.width()*.5-mradius+verticaloffset),
				 Point(samplef.width()*.5,-samplef.width()*.5+mradius+verticaloffset), 
				 Point(-samplef.width()*.5,samplef.width()*.5+mradius+verticaloffset));
	TriangularInclusion t1(Point(samplef.width()*.5,-samplef.width()*.5-mradius+verticaloffset), 
				 Point(-samplef.width()*.5,samplef.width()*.5+mradius+verticaloffset), 
				 Point(-samplef.width()*.5,samplef.width()*.5-mradius+verticaloffset));
	t0.isVirtualFeature = true ;
	t1.isVirtualFeature = true ;
	PlasticStrain * t0damagemodel = new PlasticStrain() ;
	PlasticStrain * t1damagemodel = new PlasticStrain() ;
// 	t0.setBehaviour( new StiffnessAndFracture(Material::cauchyGreen(std::make_pair(E_paste,nu), true,SPACE_TWO_DIMENSIONAL, PLANE_STRAIN) , new DruckerPrager(-20e6*.95, -20e6*.95,E_paste,0.1 , mradius),t0damagemodel));
// 	t1.setBehaviour( new StiffnessAndFracture(Material::cauchyGreen(std::make_pair(E_paste,nu), true,SPACE_TWO_DIMENSIONAL, PLANE_STRAIN) , new DruckerPrager(-20e6*.95, -20e6*.95,E_paste,0.1 , mradius),t1damagemodel));
	
	t0.setBehaviour(  new ConcreteBehaviour( E_paste, nu, compressionCrit*.95,PLANE_STRAIN, UPPER_BOUND, SPACE_TWO_DIMENSIONAL ));
	t1.setBehaviour(  new ConcreteBehaviour( E_paste, nu, compressionCrit*.95,PLANE_STRAIN, UPPER_BOUND, SPACE_TWO_DIMENSIONAL ));
	
// 	F.addFeature(&samplef, &t0);
// 	F.addFeature(&samplef, &t1);
	t0.setBehaviourSource(&samplef);
	t1.setBehaviourSource(&samplef);
	
	
// 	
// 	
// 	inc.setBehaviourSource(&samplef);
// 	samplef.setBehaviour( new StiffnessAndFracture(Material::cauchyGreen(std::make_pair(E_paste,nu), true,SPACE_TWO_DIMENSIONAL, PLANE_STRAIN) , new DruckerPrager(-20e6, -20e6,E_paste,0.1 , mradius),new PlasticStrain())) ;
// 	samplef.setBehaviour( new ConcreteBehaviour( E_paste, nu, compressionCrit,PLANE_STRAIN, UPPER_BOUND, SPACE_TWO_DIMENSIONAL ) ) ;
	samplef.setBehaviour(new Stiffness(Material::cauchyGreen(std::make_pair(E_paste,nu), true,SPACE_TWO_DIMENSIONAL, PLANE_STRAIN))) ;
	
// 	inc.setBehaviour(new VoidForm()) ;
// 	F.addFeature(&samplef, new Pore(samplef.height()*.15, samplef.getCenter().x, samplef.getCenter().y+samplef.height()*.5));
// 	F.addFeature(&samplef, new Pore(samplef.height()*.1, samplef.getCenter().x,samplef.getCenter().y-samplef.height()*.5));
	
// 	FeatureTree F(&box) ;
// 	featureTree = &F ;
// 	F.addFeature(NULL,&sample) ;        F.setSamplingFactor(&sample, 2.) ;
// 	F.addFeature(NULL,&rebarinternal) ; F.setSamplingFactor(&rebarinternal, 2.) ;
// 	F.addFeature(NULL,&samplef, 0,1.-steelfraction) ;
// 	F.addFeature(NULL,&rebarright,0,1.-steelfraction) ;
// 	F.addFeature(NULL,&toprightvoid) ;

// 	BranchedCrack * bc = new BranchedCrack(new Point(0, samplef.height()*.15), new Point(0, -samplef.height()*.15));
// 	bc->setEnrichementRadius(samplef.height()*0.) ;
// 	F.addFeature(&sample, bc) ; //add the crack to the feature tree
// 	
// 	BranchedCrack * bd = new BranchedCrack(new Point(samplef.height()*.15, 0), new Point(-samplef.height()*.15, 0));
// 	bd->setEnrichementRadius(samplef.height()*0.) ;
// 	F.addFeature(&sample, bd) ; //add the crack to the feature tree
	
// 	Vector def(3) ; 
// 	def[0] = 0.1 ;
// 	def[1] = 0.1 ;
// 	F.addFeature(&sample, new ExpansiveZone(&sample, samplef.height()*0.05, 0, 0, m0_paste*.5, def)) ;
	
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM_RIGHT)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM_LEFT)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, BOTTOM_LEFT)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, TOP_LEFT)) ;
	
	F.addBoundaryCondition(loadr);
	F.addBoundaryCondition(loadt);

// 	F.addBoundaryCondition(new BoundingBoxNearestNodeDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM, Point(-0.13335, -.0672*.5))) ;
// 	F.addBoundaryCondition(new BoundingBoxNearestNodeDefinedBoundaryCondition(FIX_ALONG_XI, BOTTOM, Point(+0.13335, -.0672*.5))) ;
// 	F.addBoundaryCondition(new BoundingBoxNearestNodeDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM, Point(+0.13335, -.0672*.5))) ;
// 	F.addBoundaryCondition(new BoundingBoxNearestNodeDefinedBoundaryCondition(FIX_ALONG_XI, BOTTOM, Point(+0.13335, -.0672*.5))) ;
	
// 	F.addBoundaryCondition( new BoundingBoxNearestNodeDefinedBoundaryCondition(FIX_ALONG_XI, TOP, Point(+0.04445, .0672*.5))) ;
// 	F.addBoundaryCondition( new BoundingBoxNearestNodeDefinedBoundaryCondition(FIX_ALONG_XI, TOP, Point(-0.04445, .0672*.5))) ;
	
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA,TOP, 0)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA,BOTTOM)) ;

	F.setSamplingNumber(atoi(argv[1])) ;
// 	F.setSamplingFactor(&rebarinternal, .5) ;
	
// 	F.addRefinementZone(new Circle(.02, 0., 0.)) ;
// 	F.addRefinementZone(new Circle(.01, 0., 0.)) ;
// 	F.addRefinementZone(new Rectangle(.06, .03, 0., 0.)) ;
	
	
// 	F.addRefinementZone(new Rectangle(0.04445*3, .0672, 0., 0.)) ;
// 	F.addRefinementZone(new Rectangle(0.04445*2, .0672, 0., 0.)) ;

	F.setOrder(LINEAR) ;
// F.addPoint(new Point(0, 0)) ;
	triangles = F.getElements2D() ;
	F.setMaxIterationsPerStep(3400);
// 	F.addPoint( new Point( 0, -sampleHeight*.5 - plateHeight ) ) ;
// 	F.addPoint(new Point(-0.04445, .0672*.5)) ;
// 	F.addPoint(new Point(+0.04445, .0672*.5)) ;
// 	F.addPoint(new Point(-0.13335, -.0672*.5)) ;
// 	F.addPoint(new Point(+0.13335, -.0672*.5)) ;
	
	
// 	F.addPoint(new Point(1.300*.5+.225, effectiveRadius*.5)) ;
// 	F.addPoint(new Point(-1.300*.5-.225, effectiveRadius*.5)) ;
	
	step(2) ;

// 	delete dt ;
	
	return 0 ;
}
