// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../utilities/samplingcriterion.h"
#include "../features/features.h"
#include "../physics/physics_base.h"
#include "../utilities/optimizer.h"
#include "../mesher/structuredmesh.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/nonlocalvonmises.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/kelvinvoight.h"
#include "../physics/fracturecriteria/vonmises.h"
#include "../physics/fracturecriteria/maxstrain.h"
#include "../physics/stiffness_and_indexed_fracture.h"
#include "../physics/damagemodels/isotropiclineardamage.h"
#include "../physics/spatially_distributed_stiffness.h"
#include "../physics/stiffness.h"
#include "../physics/void_form.h"
#include "../physics/homogeneised_behaviour.h"
#include "../physics/damagemodels/damageindexeddamage.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../physics/stiffness_and_fracture.h"
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
#include "../physics/stiffness_with_imposed_deformation.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

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
#define ID_BACK -2
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

using namespace Mu ;
using namespace std;

FeatureTree * featureTree ;
std::vector<DelaunayTriangle *> triangles ;
std::vector<bool> cracked ;
std::vector<BranchedCrack *> crack ;

double E_min = 10;
double E_max = 0;

double x_max = 0 ;
double y_max = 0 ;
double disp = 0 ;
BoundingBoxDefinedBoundaryCondition * imposeddisp = new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA, TOP, disp) ;
double width = 30;
double height = 30;
Sample sample(NULL, width , height, 0, 0) ;
double x_min = 0 ;
double y_min = 0 ;

double timepos = 0.00 ;
int grid = -1 ;
bool firstRun = true ;


int samplingnumber ;
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

double nu = 0.3 ;
double E_agg = 30000.;//softest
double E_paste = 147.25 ; //E_agg/4. ;//stiff
double E_stiff = E_agg*10. ;//stiffer
double E_soft = E_agg/10.; //stiffest

size_t current_list = DISPLAY_LIST_STRAIN_XX ;
double factor = .5 ;
MinimumAngle cri(M_PI/6.) ;
bool nothingToAdd = false ;
bool dlist = false ;
int count = 0 ;
double aggregateArea = 0;

std::vector<double> energy ;
SingleElementMesh<DelaunayTriangle, DelaunayTreeItem> * mesh = new SingleElementMesh<DelaunayTriangle, DelaunayTreeItem>(new DelaunayTriangle(NULL, NULL, 
																	  new Point(sample.getCenter().x - sample.width()/2, sample.getCenter().y - sample.height()/2), 
																	  new Point(sample.getCenter().x - sample.width()/2, sample.getCenter().y + sample.height()/2), 
																	  new Point(sample.getCenter().x + sample.width()/2, sample.getCenter().y + sample.height()/2), NULL)) ;
Point *pa = new Point(sample.getCenter().x - 2, sample.getCenter().y) ;
Point *pb = new Point(sample.getCenter().x + 2, sample.getCenter().y) ;

BranchedCrack * crack0 = new BranchedCrack(pa, pb) ;

std::vector<DelaunayTriangle *> supertris = mesh->getElements() ;

std::set<Point *> points0 ;

LeastSquaresApproximation * ls0 ;
LeastSquaresApproximation * ls1 ;
std::map<size_t, size_t> trans0 ;
std::map<size_t, size_t> trans1 ;
std::map<Point *, Point *> coincidentPoints ;
std::set<int> idsall ;

Vector origXdisps0 ;
Vector origYdisps0 ;
Vector origXdisps1 ;
Vector origYdisps1 ;

// void enrichedEquivalentElements()
// {
// 	std::vector<Point *> newTips ;
// 	newTips.push_back(pb);
// 	newTips.push_back(pc);
// 	newTips.push_back(pd);
// 	Segment ab(*pa, *pb) ;
// 	Segment cd(*pc, *pd) ;
// 	Point inter = ab.intersection(cd) ;
// 	pi->x = inter.x ;
// 	pi->y = inter.y ;
// 	crack0->branch(pi, newTips);
// 	std::cout << "setup equivalent elements...  " << std::flush ;
// 	crack0->setEnrichementRadius(sample.height()*0.0001) ;
// 	TriElement * father = new TriElement(QUADRATIC) ;
// 	father->compileAndPrecalculate() ;
// 	mesh->setElementOrder(QUADRATIC);
// 	supertris[0]->setBehaviour(new HomogeneisedBehaviour(featureTree, supertris[0])) ;
// 	supertris[0]->getState().initialize() ;
// 	supertris[0]->refresh(father) ;
// 	for(int k = 0 ; k < supertris[0]->getBoundingPoints().size() ; k++)
// 		supertris[0]->getBoundingPoint(k).id = k ;
// 	mesh->getLastNodeId()  =  supertris[0]->getBoundingPoints().size()+1;
// 	
// 	crack0->enrich(mesh->getLastNodeId(), mesh);
// 	
// 	for(size_t i = 0 ; i < triangles.size() ; i++)
// 	{
// 		if(triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR)
// 		{
// 
// 			for(size_t j = 0 ; j < triangles[i]->getBoundingPoints().size() ; j++)
// 			{
// 				Point test(triangles[i]->getBoundingPoint(j)) ;
// 				supertris[0]->project(&test) ;
// 				if(supertris[0]->in(triangles[i]->getBoundingPoint(j)) || dist(test, triangles[i]->getBoundingPoint(j)) < 100.*POINT_TOLERANCE)
// 				{
// 					points0.insert(&triangles[i]->getBoundingPoint(j)) ;
// 				}
// 				for(int k = 0 ; k < supertris[0]->getBoundingPoints().size() ; k++)
// 				{
// 					if( supertris[0]->getBoundingPoint(k) == triangles[i]->getBoundingPoint(j))
// 					{
// 						coincidentPoints[&supertris[0]->getBoundingPoint(k)] = &triangles[i]->getBoundingPoint(j) ;
// 						points0.insert(&triangles[i]->getBoundingPoint(j)) ;
// 					}
// 				}
// 				
// 				
// 			}
// 		}
// 	}
// 	std::cout << "...done." << std::endl ;
// 	std::cout << "found " << coincidentPoints.size() << " coincident points" << std::endl ;
// }

void setLSBoundaryConditions(LeastSquaresApproximation * ls, int xy)
{
	ls->clearParameterValues();
	for(auto i = coincidentPoints.begin() ; i != coincidentPoints.end() ; ++i)
	{
		ls->setParameterValue(i->first->id, x[i->second->id*2+xy]);
	}
}

void setupLeastSquares()
{
	trans0.clear();
	idsall.clear();
	std::vector< size_t > ids0 (supertris[0]->getState().getInterpolatingFactors(supertris[0]->getCenter(), false).size()) ; //= supertris[0]->getDofIds() ;
	for(size_t j = 0 ; j < ids0.size() ; j++)
	{
		trans0[j] = j ;
	}

	Matrix X0(ids0.size(), points0.size()+coincidentPoints.size()) ;
		
	origXdisps0.resize(points0.size()+coincidentPoints.size()) ;
	origYdisps0.resize(points0.size()+coincidentPoints.size()) ;
	
	int indexj = 0 ;
	for(std::set<Point *>::const_iterator i = points0.begin() ; i != points0.end() ; i++)
	{
		origXdisps0[indexj] = x[(*i)->id*2] ;
		origYdisps0[indexj] = x[(*i)->id*2+1] ;
		std::vector<double> interp = supertris[0]->getState().getInterpolatingFactors(*(*i), false) ;
		
		for(size_t j = 0 ; j < interp.size() ; j++)
		{
			X0[trans0[j]][indexj] = interp[j] ;
		}
		indexj++ ;
	}
	
	delete ls0 ;
	ls0 = new LeastSquaresApproximation(origXdisps0, X0) ;
}


// double distanceBetweenMeshes()
// {
// 
// 	crack0->enrich(mesh->getLastNodeId(), mesh);
// 	setupLeastSquares() ;
// 	setLSBoundaryConditions(ls0, 0);
// 	ls0->setMeasures(origXdisps0) ;
// 	ls0->optimize() ;
// 	Vector x0disp = ls0->getParameters() ;
// 	
// 	setLSBoundaryConditions(ls0, 1);
// 	ls0->setMeasures(origYdisps0) ;
// 	ls0->optimize() ;
// 	Vector y0disp = ls0->getParameters() ;
// 	
// 	double distancex = 0;
// 	double distancey = 0;
// 	int indexj = 0 ;
// 
// 	for(std::set<Point *>::const_iterator i = points0.begin() ; i != points0.end() ; i++)
// 	{
// 		double dispx = 0 ;
// 		double dispy = 0 ;
// 		for(size_t j = 0 ; j < ls0->getLinearModel().numRows() ; j++)
// 		{
// 			dispx += ls0->getLinearModel()[j][indexj]*x0disp[j] ;
// 			dispy += ls0->getLinearModel()[j][indexj]*y0disp[j] ;
// 		}
// 
// 		indexj++ ;
// 		if(!isnan(dispx))
// 			distancex += std::abs(dispx-x[(*i)->id*2]) ;
// 		else
// 			distancex += 1e6 ;
// 		if(!isnan(dispy))
// 			distancey += std::abs(dispy-x[(*i)->id*2+1])  ;
// 		else
// 			distancey += 1e6 ;
// 	}
// 
// 	return (distancex+distancey)/points0.size();
// }

// void optimize()
// {
// 	enrichedEquivalentElements() ;
// 
// 	std::vector<double *> vars ;
// 	vars.push_back( &pa->y) ;
// 	vars.push_back( &pb->y) ;
// 	vars.push_back( &pc->x) ;
// 	vars.push_back( &pd->x) ;
// 	vars.push_back( &pi->x) ;
// 	vars.push_back( &pi->y) ;
// 	
// 	std::vector<std::pair<double, double> > limits ;
// 	limits.push_back(std::make_pair(sample.getCenter().y-sample.height()/2.1, sample.getCenter().y+sample.height()/2.1)) ;
// 	limits.push_back(std::make_pair(sample.getCenter().y-sample.height()/2.1, sample.getCenter().y+sample.height()/2.1)) ;
// 	limits.push_back(std::make_pair(sample.getCenter().x-sample.width()/2.1, sample.getCenter().x+sample.width()/2.1)) ;
// 	limits.push_back(std::make_pair(sample.getCenter().x-sample.width()/2.1, sample.getCenter().x+sample.width()/2.1)) ;
// 	limits.push_back(std::make_pair(sample.getCenter().x-sample.width()/2, sample.getCenter().x+sample.width()/2)) ;
// 	limits.push_back(std::make_pair(sample.getCenter().x-sample.height()/2, sample.getCenter().x+sample.height()/2)) ;
// 	
// 	GeneticAlgorithmOptimizer ga( vars, limits, &distanceBetweenMeshes) ;
// 	ga.optimize(1e-12, 60, 150,  .1, .1) ;
// 
// 	pa->y = ga.getValues()[0].second ;
// 	pb->y = ga.getValues()[1].second ;
// 	pc->x = ga.getValues()[2].second ;
// 	pd->x = ga.getValues()[3].second ;
// 	pi->x = ga.getValues()[4].second ;
// 	pi->y = ga.getValues()[5].second ;
// 	crack0->enrich(mesh->getLastNodeId(), mesh);
// 	
// 	crack0->print();
// 
// 	ls0->setMeasures(origXdisps0) ;
// 	setLSBoundaryConditions(ls0, 0);
// 	ls0->optimize() ;
// 	Vector dispx = ls0->getApproximation() ;
// 	ls0->setMeasures(origYdisps0) ;
// 	setLSBoundaryConditions(ls0, 1);
// 	ls0->optimize() ;
// 	Vector dispy = ls0->getApproximation() ;
// 	
// 	ls0->printParameters() ;
// 	int indexj = 0 ;
// 	for(std::set<Point *>::const_iterator i = points0.begin() ; i != points0.end() ; i++)
// 	{
// 		std::cout << (*i)->x << "  " << (*i)->y << "  " << x[(*i)->id*2] << "  "<< x[(*i)->id*2+1] << "  " << dispx[indexj] << "  " << dispy[indexj] << std::endl;
// 		indexj++ ;
// 	}
// 	
// }

bool go = true ;
void step()
{
	
	size_t max_growth_steps = 1;
	size_t max_limit = 2000 ;
	int limit = 0 ;
	

		go = featureTree->step() ;

// 		if(go)
// 		{
// 			imposeddisp->setData(imposeddisp->getData()+.1);
// 			go = true ;
// 		}
		double da = 0 ;
		
		triangles = featureTree->getElements2D(grid) ;
		x.resize(featureTree->getDisplacements(grid).size()) ;
		x = featureTree->getDisplacements(grid) ;
		sigma.resize(triangles.size()*triangles[0]->getBoundingPoints().size()*3) ;
		epsilon.resize(triangles.size()*triangles[0]->getBoundingPoints().size()*3) ;
	
		std::pair<Vector, Vector > sigma_epsilon = featureTree->getStressAndStrain(grid) ;
		sigma.resize(sigma_epsilon.first.size()) ;
		sigma = sigma_epsilon.first ;
		epsilon.resize(sigma_epsilon.second.size()) ;
		epsilon = sigma_epsilon.second ;
		
		Vector avgdisplacement(2) ;
		double avgdisplacementarea(0) ;
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
// 		if(countit%100 == 0)
// 			std::cout << "unknowns :" << x.size() << std::endl ;
		
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
		
		
			if(triangles[k]->getBehaviour() && !in && triangles[k]->getBehaviour()->type != VOID_BEHAVIOUR)
			{
				avgdisplacement += triangles[k]->getState().getAverageDisplacement()*triangles[k]->area() ;
				avgdisplacementarea += triangles[k]->area() ;
			}
			
			
			if(triangles[k]->getBehaviour() && !in && !triangles[k]->getBehaviour()->fractured() && triangles[k]->getBehaviour()->type != VOID_BEHAVIOUR)
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
					if(triangles[k]->getBoundingPoint(p).x > 0.0799)
					{
						e_xx+=x[triangles[k]->getBoundingPoint(p).id*2] ;
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
					}
					enr += triangles[k]->getState().elasticEnergy() ;
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
		avgdisplacement /= avgdisplacementarea ;
		
		if(go)
		{
			
			for(double k = 0 ; k < 4 ; k += 4./800)
			{
				for(double l = -2 ; l < 2 ; l += 4./800)
				{
					Point p(k, l) ;
					auto tri = featureTree->getElements2D(&p) ;
					
					for(size_t j = 0 ; j  < tri.size() ; j++ )
					{
						if(tri[j]->in(p))
						{
							Vector s = tri[j]->getState().getStress(p, false) ;
							Vector e = tri[j]->getState().getStrain(p, false) ;
							std::cout << 0.5*(s[0]*e[0] + s[1]*e[1]+ s[2]*e[2])<< "  " << std::flush ;
// 							std::cout << tri[j]->getState().getStress(p, false)[0] << "  " << std::flush ;
// 							std::cout << tri[j]->getState().getDisplacements(p, false)[0] << "  " << std::flush ;
							break ;
						}
					}
				}
				std::cout << std::endl ;
			}
			std::cout << std::endl ;
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
			
			std::cout << "average sigma11 : " << (avg_s_xx/area)/1000000. << std::endl ;
			std::cout << "average sigma22 : " << (avg_s_yy/area)/1000000. << std::endl ;
			std::cout << "average sigma12 : " << (avg_s_xy/area)/1000000. << std::endl ;
			std::cout << "average epsilon11 : " << avg_e_xx/area<< std::endl ;
			std::cout << "average epsilon22 : " << avg_e_yy/area << std::endl ;
			std::cout << "average epsilon12 : " << avg_e_xy/area << std::endl ;
		}
		energy.push_back(enr) ;

// 	for(size_t i = 0 ; i < energy.size() ; i++)
// 		std::cout << energy[i] << std::endl ;
}

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
// 			imposeddisp->setData(imposeddisp->getData()+0.025);
		step() ;
			
			dlist = false ;
			break ;
		}
	case ID_BACK:
	{
// 		samplingnumber *= 1.5 ;
// 		featureTree->setSamplingNumber(samplingnumber);
		dynamic_cast<Inclusion *>(featureTree->getFeature(1))->setRadius(featureTree->getFeature(1)->getRadius()+0.1 ) ;
		dynamic_cast<Inclusion *>(featureTree->getFeature(2))->setRadius(featureTree->getFeature(2)->getRadius()+0.1 ) ;
		step() ;
// 		imposeddisp->setData(imposeddisp->getData()-.1);
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
	case ID_FRAC_CRIT: 
		{
			current_list = DISPLAY_LIST_FRAC_CRIT ;
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
// 			sigma11 *= 1.5 ;
// 			sigma22 *= 1.5 ;
// 			sigma12 *= 1.5 ;
			dlist = false ;
			break ;
		}
	case ID_DEAMPLIFY :
		{
			x /= 1.5 ;
// 			sigma11 /= 1.5 ;
// 			sigma22 /= 1.5 ;
// 			sigma12 /= 1.5 ;
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
	
// 	double x_max = std::abs(x).min() ;
// 	double y_max = std::abs(x).min() ;
// 	
// 	for(size_t k = 0 ; k < x.size()/2 ; k++)
// 	{
// 		if(x[k*2]*x[k*2]+x[k*2+1]*x[k*2+1] > x_max*x_max+y_max*y_max )
// 		{
// 			x_max = x[k*2] ;
// 			y_max = x[k*2+1] ;
// 		}
// 	}
	int start = 0 ;
	for (unsigned int j=0 ; j< triangles.size() ; j++ )
	{
		if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR)
		{
			int timePlanes = triangles[j]->timePlanes() ;
			start = ((triangles[j]->getBoundingPoints().size()-1)*(timePlanes-1))/timePlanes ;
			break ;
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
				
				double vx = x[triangles[j]->getBoundingPoint(start).id*2]; 
				double vy = x[triangles[j]->getBoundingPoint(start).id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - sqrt(((vx-x_min)*(vx-x_min) + (vy-y_min)*(vy-y_min))/((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min)))*300., 1., 1. ) ;
				glColor3f(c1, c2, c3) ;
				
				glVertex2f(double(triangles[j]->getBoundingPoint(start).x + vx) , double(triangles[j]->getBoundingPoint(start).y + vy) );
				
				for(size_t k = start+1 ; k < triangles[j]->getBoundingPoints().size() ; k++)
				{
					vx = x[triangles[j]->getBoundingPoint(k).id*2];
					vy = x[triangles[j]->getBoundingPoint(k).id*2+1]; 
					
					HSVtoRGB( &c1, &c2, &c3, 300. - sqrt(((vx-x_min)*(vx-x_min) + (vy-y_min)*(vy-y_min))/((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min)))*300., 1., 1. ) ;
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
				
				double vx = x[triangles[j]->getBoundingPoint(start).id*2]; 
				double vy = x[triangles[j]->getBoundingPoint(start).id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(sigma11[j*triangles[j]->getBoundingPoints().size()]-sigma11_min)/(sigma11_max-sigma11_min), 1., 1.) ;
				
				glVertex2f(double(triangles[j]->getBoundingPoint(start).x + vx) , double(triangles[j]->getBoundingPoint(start).y + vy) );
				
				for(size_t k = start+1 ; k < triangles[j]->getBoundingPoints().size() ; k++)
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
				
				double vx = x[triangles[j]->getBoundingPoint(start).id*2]; 
				double vy = x[triangles[j]->getBoundingPoint(start).id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(vonMises[j*triangles[j]->getBoundingPoints().size()]-vonMises_min)/(vonMises_max-vonMises_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				
				glVertex2f(double(triangles[j]->getBoundingPoint(start).x + vx) , double(triangles[j]->getBoundingPoint(start).y + vy) );
				
				for(size_t k = 1+start ; k < triangles[j]->getBoundingPoints().size() ; k++)
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
				
				double vx = x[triangles[j]->getBoundingPoint(start).id*2]; 
				double vy = x[triangles[j]->getBoundingPoint(start).id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(angle[j*triangles[j]->getBoundingPoints().size()]-angle_min)/(angle_max-angle_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				
				glVertex2f(double(triangles[j]->getBoundingPoint(start).x + vx) , double(triangles[j]->getBoundingPoint(start).y + vy) );
				
				for(size_t k = 1+start ; k < triangles[j]->getBoundingPoints().size() ; k++)
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

		double crit_max = fracCrit.max() ;
		double crit_min = fracCrit.min() ;
		glNewList(  DISPLAY_LIST_FRAC_CRIT,  GL_COMPILE ) ;
		for (unsigned int j=0 ; j< triangles.size() ; j++ )
		{
			
			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !cracked[j]&& !triangles[j]->getBehaviour()->fractured())
			{
				double c1 ;
				double c2 ;
				double c3 ;
				
				double vx = x[triangles[j]->getBoundingPoint(start).id*2]; 
				double vy = x[triangles[j]->getBoundingPoint(start).id*2+1]; 
				double s = 1. ;
				if(fracCrit[j*triangles[j]->getBoundingPoints().size()] <= 0)
					s = .2 ;
					
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(fracCrit[j*triangles[j]->getBoundingPoints().size()]-crit_min)/(crit_max-crit_min), s, 1.) ;
				glColor3f(c1, c2, c3) ;
				
				glVertex2f(double(triangles[j]->getBoundingPoint(start).x + vx) , double(triangles[j]->getBoundingPoint(start).y + vy) );
				
				for(size_t k = 1+start ; k < triangles[j]->getBoundingPoints().size() ; k++)
				{
					vx = x[triangles[j]->getBoundingPoint(k).id*2];
					vy = x[triangles[j]->getBoundingPoint(k).id*2+1]; 
					
					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(fracCrit[j*triangles[j]->getBoundingPoints().size()+k]-crit_min)/(crit_max-crit_min), s, 1.) ;
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
				
				double vx = x[triangles[j]->getBoundingPoint(start).id*2]; 
				double vy = x[triangles[j]->getBoundingPoint(start).id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(sigma22[j*triangles[j]->getBoundingPoints().size()]-sigma22_min)/(sigma22_max-sigma22_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				glVertex2f(double(triangles[j]->getBoundingPoint(start).x + vx) , double(triangles[j]->getBoundingPoint(start).y + vy) );
				
				for(size_t k = 1+start ; k < triangles[j]->getBoundingPoints().size() ; k++)
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
				
				double vx = x[triangles[j]->getBoundingPoint(start).id*2]; 
				double vy = x[triangles[j]->getBoundingPoint(start).id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(sigma12[j*triangles[j]->getBoundingPoints().size()]-sigma12_min)/(sigma12_max-sigma12_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				glVertex2f(double(triangles[j]->getBoundingPoint(start).x + vx) , double(triangles[j]->getBoundingPoint(start).y + vy) );
				
				for(size_t k = 1+start ; k < triangles[j]->getBoundingPoints().size() ; k++)
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
				
				double vx = x[triangles[j]->getBoundingPoint(start).id*2]; 
				double vy = x[triangles[j]->getBoundingPoint(start).id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				
				Point a = triangles[j]->inLocalCoordinates(triangles[j]->getBoundingPoint(start)) ;
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(triangles[j]->getBehaviour()->getTensor(a)[2][2]-0)/(E_max-0), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				glVertex2f(double(triangles[j]->getBoundingPoint(start).x + vx) , double(triangles[j]->getBoundingPoint(start).y + vy) );
				
				for(size_t k = 1+start ; k < triangles[j]->getBoundingPoints().size() ; k++)
				{
					vx = x[triangles[j]->getBoundingPoint(k).id*2];
					vy = x[triangles[j]->getBoundingPoint(k).id*2+1]; 
					a = triangles[j]->inLocalCoordinates(triangles[j]->getBoundingPoint(k)) ;
					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(triangles[j]->getBehaviour()->getTensor(a)[2][2]-0)/(E_max-0), 1., 1.) ;
					glColor3f(c1, c2, c3) ;
					glVertex2f( double(triangles[j]->getBoundingPoint(k).x + vx) ,  double(triangles[j]->getBoundingPoint(k).y + vy) );
					
				}
				glEnd() ;
			}
		}
		glEndList() ;
		
		glNewList(  DISPLAY_LIST_STIFFNESS_DARK,  GL_COMPILE ) ;
		for (unsigned int j=0 ; j< triangles.size() ; j++ )
		{
			
			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !cracked[j]&& !triangles[j]->getBehaviour()->fractured())
			{
				double c1 ;
				double c2 ;
				double c3 ;
				
				double vx = x[triangles[j]->getBoundingPoint(start).id*2]; 
				double vy = x[triangles[j]->getBoundingPoint(start).id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				
				Point a = triangles[j]->inLocalCoordinates(triangles[j]->getBoundingPoint(start)) ;
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(triangles[j]->getBehaviour()->getTensor(a)[0][0]-E_min)/(E_max-E_min), 1., .2) ;
				glColor3f(c1, c2, c3) ;
				glVertex2f(double(triangles[j]->getBoundingPoint(start).x + vx) , double(triangles[j]->getBoundingPoint(start).y + vy) );
				
				for(size_t k = 1+start ; k < triangles[j]->getBoundingPoints().size() ; k++)
				{
					vx = x[triangles[j]->getBoundingPoint(k).id*2];
					vy = x[triangles[j]->getBoundingPoint(k).id*2+1]; 
					a = triangles[j]->inLocalCoordinates(triangles[j]->getBoundingPoint(k)) ;
					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(triangles[j]->getBehaviour()->getTensor(a)[0][0]-E_min)/(E_max-E_min), 1., .2) ;
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
				
				double vx = x[triangles[j]->getBoundingPoint(start).id*2]; 
				double vy = x[triangles[j]->getBoundingPoint(start).id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(epsilon11[j*triangles[j]->getBoundingPoints().size()]-epsilon11_min)/(epsilon11_max-epsilon11_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				glVertex2f(double(triangles[j]->getBoundingPoint(start).x + vx) , double(triangles[j]->getBoundingPoint(start).y + vy) );
				
				for(size_t k = 1+start ; k < triangles[j]->getBoundingPoints().size() ; k++)
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
		double epsilon22_max =  epsilon22.max() ;
		
		glNewList(  DISPLAY_LIST_STRESS_YY,  GL_COMPILE ) ;
		for (unsigned int j=0 ; j< triangles.size() ; j++ )
		{
			
			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !cracked[j]&& !triangles[j]->getBehaviour()->fractured())
			{
				double c1 ;
				double c2 ;
				double c3 ;
				
				double vx = x[triangles[j]->getBoundingPoint(start).id*2]; 
				double vy = x[triangles[j]->getBoundingPoint(start).id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(epsilon22[j*triangles[j]->getBoundingPoints().size()]-epsilon22_min)/(epsilon22_max-epsilon22_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				glVertex2f(double(triangles[j]->getBoundingPoint(start).x + vx) , double(triangles[j]->getBoundingPoint(start).y + vy) );
				
				for(size_t k = 1+start ; k < triangles[j]->getBoundingPoints().size() ; k++)
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
		double epsilon12_max =  epsilon12.max() ;
		
		glNewList(  DISPLAY_LIST_STRESS_XY,  GL_COMPILE ) ;
		for (unsigned int j=0 ; j< triangles.size() ; j++ )
		{
			
			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !cracked[j]&& !triangles[j]->getBehaviour()->fractured())
			{
				double c1 ;
				double c2 ;
				double c3 ;
				
				double vx = x[triangles[j]->getBoundingPoint(start).id*2]; 
				double vy = x[triangles[j]->getBoundingPoint(start).id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(epsilon12[j*triangles[j]->getBoundingPoints().size()]-epsilon12_min)/(epsilon12_max-epsilon12_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				glVertex2f(double(triangles[j]->getBoundingPoint(start).x + vx) , double(triangles[j]->getBoundingPoint(start).y + vy) );
				
				for(size_t k = 1+start ; k < triangles[j]->getBoundingPoints().size() ; k++)
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
				
				int enrichment = triangles[j]->getEnrichmentFunctions().size() ;
				//HSVtoRGB( &c1, &c2, &c3, 180. + 180.*(sigma12[j]-sigma12.min())/(sigma12.max()-sigma12.min()), 1., 1. ) 
				
				
				if(enrichment)
				{
					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*enrichment/20., 1., 1.) ;
					glColor3f(c1, c2, c3) ;
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
// 				else
// 				{
// 					double vx = x[triangles[j]->first->id*2]; 
// 					double vy = x[triangles[j]->first->id*2+1]; 
// 					Point a = triangles[j]->inLocalCoordinates(triangles[j]->getBoundingPoint(0)) ;
// 					
// 					HSVtoRGB( &c1, &c2, &c3, 0,0, 300. - 300.*(triangles[j]->getBehaviour()->getTensor(a)[0][0]-E_min)/(E_max-E_min)) ;
// 					glColor3f(c1, c2, c3) ;
// 					
// 					glVertex2f( double(triangles[j]->first->x + vx) ,
// 					            double(triangles[j]->first->y + vy) );
// 					
// 					vx = x[triangles[j]->second->id*2];
// 					vy = x[triangles[j]->second->id*2+1]; 
// 					
// 					glVertex2f( double(triangles[j]->second->x + vx) ,
// 					            double(triangles[j]->second->y + vy) );
// 					
// 					
// 					vx = x[triangles[j]->third->id*2]; 
// 					vy = x[triangles[j]->third->id*2+1]; 
// 					
// 					
// 					glVertex2f( double(triangles[j]->third->x + vx) ,
// 					            double(triangles[j]->third->y + vy) );
// 					
// 
// 				}
			}
		}
		glEnd();
		glEndList() ;
		
		
		glNewList( DISPLAY_LIST_ELEMENTS,  GL_COMPILE ) ;
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
				for(size_t k = start ; k < triangles[j]->getBoundingPoints().size() ; k++)
				{
					double vx = 0 ; //x[triangles[j]->getBoundingPoint(k).id*2]; 
					double vy =  0 ; //x[triangles[j]->getBoundingPoint(k).id*2+1]; 
					
					glVertex2f( double(triangles[j]->getBoundingPoint(k).x+vx) ,  double(triangles[j]->getBoundingPoint(k).y+vy) );
					
				}
				glEnd();
			}
			else
			{
				glColor3f(0, 0, 1) ;

				
				glBegin(GL_LINE_LOOP);
				for(size_t k = start ; k < triangles[j]->getBoundingPoints().size() ; k++)
				{
					double vx =  0 ; //x[triangles[j]->getBoundingPoint(k).id*2]; 
					double vy =  0 ; //x[triangles[j]->getBoundingPoint(k).id*2+1]; 
					
					glVertex2f( double(triangles[j]->getBoundingPoint(k).x+vx) ,  double(triangles[j]->getBoundingPoint(k).y+vy) );
					
				}
				glEnd();
			}
			
			glColor3f(1, 1, 1) ;
		}
		glEndList() ;
		
		
		glNewList(  DISPLAY_LIST_CRACK,  GL_COMPILE ) ;

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
		
		if(current_list != DISPLAY_LIST_ENRICHMENT)
			glCallList(current_list) ;
		if(current_list == DISPLAY_LIST_ENRICHMENT)
		{
			glCallList(DISPLAY_LIST_STIFFNESS_DARK) ;
			glCallList(current_list) ;
			
		}
		
		glCallList(DISPLAY_LIST_CRACK) ;
// 		if(current_list == DISPLAY_LIST_ELEMENTS)
// 			glCallList(DISPLAY_LIST_CRACK) ;
		
		glColor3f(1, 1, 1) ;
	}
	glColor3f(1, 0, 0) ;
	glFlush();
	glutSwapBuffers();
}


double llx = 0 ;
double lly = 0 ;
// double simpleObjective()
// {
// 	return (1.-llx)*(1.-llx) + 100.*(lly-llx*llx)*(lly-llx*llx) ;
// }


struct Block
{
	Block()
	{
		series = rand()%2 ;
		stiff = 1. ;
	}
	
	Block(int nblocks)
	{
		int selfblocks = 1 ;
		if(nblocks>1)
			selfblocks = rand()%(nblocks)+1 ;
		if(selfblocks > 1)
		{
			std::vector<int> subblocks(selfblocks) ;
			for(int i = 0 ; i < selfblocks ; i++)
			{
				subblocks[i]++ ;
				nblocks-- ;
			}
			
			while(nblocks)
			{
				subblocks[rand()%selfblocks]++ ;
				nblocks-- ;
			}
			for(int i = 0 ; i < selfblocks ; i++)
			{
				blocks.push_back(new Block(subblocks[i])) ;
			}
		}
		else if(nblocks > 1)
		{
			blocks.push_back(new Block(nblocks)) ;
		}
		
		series = rand()%2 ;
		stiff = 1. ;
	}
	
	~Block()
	{
		for(int i = 0 ; i < blocks.size() ; i++)
			delete blocks[i] ;
	}
	std::vector<Block *> blocks ;
	double stiff ;
	bool series ;
	double equivStifness() const
	{
		if(blocks.size() == 0)
			return stiff ;
		if(series)
		{
			double k = 0 ;
			for(int i = 0 ; i < blocks.size() ; i++)
			{
				k += 1./blocks[i]->equivStifness() ;
			}
			
			return 1./k ;
		}
		else
		{
			double k = 0 ;
			for(int i = 0 ; i < blocks.size() ; i++)
			{
				k += blocks[i]->equivStifness() ;
			}
			
			return k ;
		}
	}
} ;

int main(int argc, char *argv[])
{

  // Material behaviour of the matrix
	Matrix m0_paste(3,3) ;
	m0_paste[0][0] = E_paste/(1.-nu*nu) ; m0_paste[0][1] =E_paste/(1.-nu*nu)*nu ; m0_paste[0][2] = 0 ;
	m0_paste[1][0] = E_paste/(1.-nu*nu)*nu ; m0_paste[1][1] = E_paste/(1.-nu*nu) ; m0_paste[1][2] = 0 ; 
	m0_paste[2][0] = 0 ; m0_paste[2][1] = 0 ; m0_paste[2][2] = .99*E_paste/(1.-nu*nu)*(1.-nu)/2. ; 

	// Material behaviour of the fibres
	Matrix m0_agg(3,3) ;
	m0_agg[0][0] = E_agg/(1-nu*nu) ; m0_agg[0][1] =E_agg/(1-nu*nu)*nu ; m0_agg[0][2] = 0 ;
	m0_agg[1][0] = E_agg/(1-nu*nu)*nu ; m0_agg[1][1] = E_agg/(1-nu*nu) ; m0_agg[1][2] = 0 ; 
	m0_agg[2][0] = 0 ; m0_agg[2][1] = 0 ; m0_agg[2][2] = E_agg/(1-nu*nu)*(1.-nu)/2. ; 

	// Material behaviour for the "very" stiff inclusion
	Matrix m0_stiff(3,3) ;
	m0_stiff[0][0] = E_stiff/(1.-nu*nu) ; m0_stiff[0][1] =E_stiff/(1.-nu*nu)*nu ; m0_stiff[0][2] = 0 ;
	m0_stiff[1][0] = E_stiff/(1.-nu*nu)*nu ; m0_stiff[1][1] = E_stiff/(1.-nu*nu) ; m0_stiff[1][2] = 0 ; 
	m0_stiff[2][0] = 0 ; m0_stiff[2][1] = 0 ; m0_stiff[2][2] = E_stiff/(1.-nu*nu)*(1.-nu)/2. ; 

	// Material behaviour for the "very" soft inclusion
	Matrix m0_soft(3,3) ;
	m0_soft[0][0] = E_soft/(1.-nu*nu) ; m0_soft[0][1] =E_soft/(1.-nu*nu)*nu ; m0_soft[0][2] = 0 ;
	m0_soft[1][0] = E_soft/(1.-nu*nu)*nu ; m0_soft[1][1] = E_soft/(1.-nu*nu) ; m0_soft[1][2] = 0 ; 
	m0_soft[2][0] = 0 ; m0_soft[2][1] = 0 ; m0_soft[2][2] = E_soft/(1.-nu*nu)*(1.-nu)/2. ; 

	Matrix d(3,3) ;
	d[0][0] = .1*E_paste ;
	d[1][1] = .1*E_paste ;
	d[2][2] = .1*E_paste ;
	FeatureTree F(&sample, 2) ;
	featureTree = &F ;

	Sample sm(0.2, 0.1, 0, 0) ;
// 	std::vector<Feature *> incs = AggregateDistribution2DGenerator(sm.area(), 0.016, 0.00002, .72, 65).getFeatures(sm.getPrimitive()) ;
// 	exit (0);
// 	for(int i = 0 ; i < incs.size() ; i++)
// 	{
// 		incs[i]->setBehaviour(new Stiffness(m0_paste*4)) ;
// 		dynamic_cast<Inclusion * >(incs[i])->setRadius(incs[i]->getRadius()*500./0.2);
// 		incs[i]->setCenter(incs[i]->getCenter()* 500./0.2)  ;
// 		F.addFeature(&sample, incs[i]);
// 	}
//  	sample.setBehaviour(new WeibullDistributedStiffness(m0_paste, 50./8)) ;

	double cradius = 10 ;
	double mradius = 5 ;
	double tdamage = .999 ;
	double dincrement = .01 ;
	IsotropicLinearDamage * dfunc = new IsotropicLinearDamage(2, .01) ;
	dfunc->setMaterialCharacteristicRadius(mradius) ;
	dfunc->setThresholdDamageDensity(tdamage);
	dfunc->setDamageDensityIncrement(dincrement);
	
	PseudoPlastic * psp = new PseudoPlastic(m0_paste, 0.156, mradius*.175) ;
// 	psp->crit->setNeighbourhoodRadius(cradius);
// 	psp->crit->setMaterialCharacteristicRadius(mradius);
// 	StiffnessAndFracture * saf = new StiffnessAndFracture(m0_paste, new VonMises(35), cradius) ; //1.5640 ; 5625 too low ; 5650 too high
	StiffnessAndFracture * saf = new StiffnessAndFracture(m0_paste, new MohrCoulomb(10, -10*8) , cradius) ; 
	saf->dfunc->setMaterialCharacteristicRadius(mradius) ;
	saf->criterion->setMaterialCharacteristicRadius(mradius);
	saf->criterion->setNeighbourhoodRadius(cradius);
	saf->dfunc->setThresholdDamageDensity(tdamage);
	saf->dfunc->setDamageDensityIncrement(dincrement);
	StiffnessAndIndexedFracture * saif = new StiffnessAndIndexedFracture(m0_paste, new /*MohrCoulomb(0.01, -0.01)*/ VonMises(0.01), cradius) ; //1.5640; 5625 too low ; 5650 too high
	saif->dfunc->setMaterialCharacteristicRadius(mradius) ;
	saif->criterion->setMaterialCharacteristicRadius(mradius);
	saif->criterion->setNeighbourhoodRadius(cradius);
	saif->dfunc->setThresholdDamageDensity(tdamage);
	saif->dfunc->setDamageDensityIncrement(dincrement);
	Stiffness * sf = new Stiffness(m0_paste) ;

// 	sample.setBehaviour(saf) ;
// 	sample.setBehaviour(saif) ;
// 		sample.setBehaviour(new WeibullDistributedStiffness(m0_paste, -37.0e6, 2)) ;
// 	dynamic_cast<WeibullDistributedStiffness *>(sample.getBehaviour())->variability = 0. ;
// 	dynamic_cast<WeibullDistributedStiffness *>(sample.getBehaviour())->materialRadius = mradius ;
// 	dynamic_cast<WeibullDistributedStiffness *>(sample.getBehaviour())->neighbourhoodRadius =  3.;
// 	sample.setBehaviour(psp) ;
	sample.setBehaviour(sf) ;
//	sample.setBehaviour(new StiffnessAndFracture(m0_paste, new VonMises(25))) ;
// 	sample.setBehaviour(new KelvinVoight(m0_paste, m0_paste*100.)) ;
// 	F.addFeature(&sample, new Pore(2, -7,2) );
// 	F.addFeature(&sample, new Pore(2, 7,-2) );
// 	F.addFeature(&sample, new Pore(20, 0, 0) );
// 	F.addFeature(&sample, new Pore(20, 85, -85) );
// 	F.addFeature(&sample, new Pore(20, 155, -155) );
// 	
// 	F.addFeature(&sample, new Pore(20, -155, -155) );
// 	F.addFeature(&sample, new Pore(20, -85, -85) );
// 	F.addFeature(&sample, new Pore(20, 85, 85) );
// 	F.addFeature(&sample, new Pore(20, 155, 155) );
	
// 	F.addFeature(&sample, new Pore(20, -250, 0) );
// 	F.addFeature(&sample, new Pore(20, -200, 0) );
// 	F.addFeature(&sample, new Pore(20, -150, 0) );
// 	F.addFeature(&sample, new Pore(20, -100, 0) );
// 	F.addFeature(&sample, new Pore(20, -50, 0) );
// 	F.addFeature(&sample, new Pore(20, 0, 0) );
// 	F.addFeature(&sample, new Pore(20, 50, -0) );
// 	F.addFeature(&sample, new Pore(20, 100, -0) );
// 	F.addFeature(&sample, new Pore(20, 150, -0) );
// 	F.addFeature(&sample, new Pore(20, 200, -0) );
// 	F.addFeature(&sample, new Pore(20, 250, -0) );

	Vector a(0., 3) ; a[0] = 1 ; a[1] = 1 ; a[2] = 0 ;
// 	ExpansiveZone * inc0 = new ExpansiveZone(&sample, 1, 0, 0.,m0_paste*2., a) ;

// 	ExpansiveZone * inc0 = new ExpansiveZone(&sample, 0.1, 0.5, 0.,m0_paste*2., a) ;
// 	ExpansiveZone * inc1 = new ExpansiveZone(&sample, 0.1, -0.5, 0.,m0_paste*2., a) ;
	Inclusion * inc0 = new Inclusion(1, 0, 0.) ;
// 	Inclusion * inc1 = new Inclusion(0.1, -0.5, 0.) ;
// 	
// 	inc0->setBehaviour(new PseudoPlastic(m0_paste*2., new MohrCoulomb(20./8, -20), new IsotropicLinearDamage(2, .01))) ;
// 	inc0->setBehaviour(new VoidForm()) ;
	inc0->setBehaviour(new StiffnessWithImposedDeformation(m0_paste*2., a)) ;
// 	inc1->setBehaviour(new StiffnessWithImposedDeformation(m0_paste*2., a)) ;
	F.addFeature(&sample, inc0) ;
// 	F.addFeature(&inc0, inc1) ;
	

// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA , TOP, -1)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI , TOP/*_LEFT*/)) ;
// 	F.addBoundaryCondition(imposeddisp) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA , BOTTOM_RIGHT)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA , BOTTOM_LEFT)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI , BOTTOM_LEFT)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI , TOP_LEFT)) ;
	
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA , TOP)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA, BOTTOM, -10)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI , LEFT)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI , RIGHT)) ;

// 	std::vector<Point *> newTips ;
// 	newTips.push_back(pb);
// 	newTips.push_back(pc);
// 	newTips.push_back(pd);
// 
// 	crack0->branch(pi, newTips);
	crack0->setEnrichementRadius(sample.height()*0.1) ;
// 	F.addFeature(&sample, crack0);
	
	samplingnumber = atoi(argv[1]);
	F.setSamplingNumber(samplingnumber) ;
	F.setOrder(QUADRATIC) ;
	F.setMaxIterationsPerStep(200) ;
	F.setDeltaTime(0.1);

	std::cout << "# max value x ; " << "mean value x ; " <<  "min value x ; " << "max value y ; " << "mean value y ;" << "min value y ; " << "max sigma11 ; " << "min sigma11 ; " << "max sigma12 ; " << "min sigma12 ; " << "max sigma22 ; " << "min sigma22 ; " << "max epsilon11 ; " << "min epsilon11 ; " << "max epsilon12 ; " << "min epsilon12 ; " << "max epsilon22 ; " << "min epsilon22 ; " << "max von Mises : " << "min von Mises : " << "average sigma11 ; " << "average sigma22 ; " << "average sigma12 ; " << "average epsilon11 ; " << "average epsilon22 ; " << "average epsilon12 ; " << "energy index ;" <<  std::endl ;
	step() ;
	
	
	glutInit(&argc, argv) ;	
	glutInitDisplayMode(GLUT_RGBA) ;
	glutInitWindowSize(600, 600) ;
	glutReshapeFunc(reshape) ;
	glutCreateWindow("S&F !") ;
	
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
	glutAddMenuEntry(" Frac Crit     ", ID_FRAC_CRIT);
	
	glutCreateMenu(Menu) ;

 	glutAddMenuEntry(" Step (fix damage) ", ID_NEXT);
	glutAddMenuEntry(" Grow pore         ", ID_BACK);
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
// 	delete saf ;
	
	return 0 ;
}
