
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
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
#include "../physics/fracturecriteria/nonlocalvonmises.h"
#include "../physics/fracturecriteria/maxstrain.h"
#include "../physics/stiffness_and_indexed_fracture.h"
#include "../physics/damagemodels/isotropiclineardamage.h"
#include "../physics/damagemodels/plasticstrain.h"
#include "../physics/spatially_distributed_stiffness.h"
#include "../physics/stiffness.h"
#include "../physics/orthotropicstiffness.h"
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
#include "../utilities/writer/triangle_writer.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include <fstream>

#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h> 
#define DEBUG 
using namespace Amie ;
using namespace std;

FeatureTree * featureTree ;
std::vector<DelaunayTriangle *> triangles ;
std::vector<bool> cracked ;
std::vector<BranchedCrack *> crack ;

MultiTriangleWriter writer("triangles_head","triangles_layers",nullptr) ;

double E_min = 10;
double E_max = 0;

double x_max = 0 ;
double y_max = 0 ;
double disp = 1 ;//.350 .355
BoundingBoxDefinedBoundaryCondition * imposeddispright = new BoundingBoxDefinedBoundaryCondition(SET_STRESS_XI, RIGHT, 0) ;
BoundingBoxDefinedBoundaryCondition * imposeddisptop = new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP, 0) ;

double width = 20;
double height = 10;
Sample sample(nullptr, width , height, 0, 0) ;
double x_min = 0 ;
double y_min = 0 ;

double timepos = 0.00 ;
int grid = -1 ;
bool firstRun = true ;


int samplingnumber ;
std::vector<DelaunayTriangle *> tris__ ;

std::pair<std::vector<Inclusion * >, std::vector<Pore * > > i_et_p ;

std::vector<std::pair<ExpansiveZone *, Inclusion *> > zones ;

std::vector<double> scales ;
std::vector<double> disps ;

Vector g_count(0) ;

double nu = 0.3 ;
double E_agg = 30000.;//softest
double E_paste = 147.25 ; //E_agg/4. ;//stiff
double E_stiff = E_agg*10. ;//stiffer
double E_soft = E_agg/10.; //stiffest

double factor = .5 ;
MinimumAngle cri(M_PI/6.) ;
bool nothingToAdd = false ;
bool dlist = false ;
int count = 0 ;
double aggregateArea = 0;

std::vector<double> energy ;
SingleElementMesh<DelaunayTriangle, DelaunayTreeItem> * mesh = new SingleElementMesh<DelaunayTriangle, DelaunayTreeItem>(new DelaunayTriangle(nullptr, nullptr, 
																	  new Point(sample.getCenter().getX() - sample.width()/2, sample.getCenter().getY() - sample.height()/2), 
																	  new Point(sample.getCenter().getX() - sample.width()/2, sample.getCenter().getY() + sample.height()/2), 
																	  new Point(sample.getCenter().getX() + sample.width()/2, sample.getCenter().getY() + sample.height()/2), nullptr), SPACE_TWO_DIMENSIONAL) ;
Point *pa = new Point(sample.getCenter().getX() - 2, sample.getCenter().getY()) ;
Point *pb = new Point(sample.getCenter().getX() + 2, sample.getCenter().getY()) ;

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
// 	pi->getX() = inter.getX() ;
// 	pi->getY() = inter.getY() ;
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
// 		supertris[0]->getBoundingPoint(k).getId() = k ;
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

void setLSBoundaryConditions(LeastSquaresApproximation * ls, int xy, const std::valarray<double> & x)
{
	ls->clearParameterValues();
	for(auto i = coincidentPoints.begin() ; i != coincidentPoints.end() ; ++i)
	{
		ls->setParameterValue(i->first->getId(), x[i->second->getId()*2+xy]);
	}
}

void setupLeastSquares(const std::valarray<double> & x)
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
		origXdisps0[indexj] = x[(*i)->getId()*2] ;
		origYdisps0[indexj] = x[(*i)->getId()*2+1] ;
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
// 			distancex += std::abs(dispx-x[(*i)->getId()*2]) ;
// 		else
// 			distancex += 1e6 ;
// 		if(!isnan(dispy))
// 			distancey += std::abs(dispy-x[(*i)->getId()*2+1])  ;
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
// 	vars.push_back( &pa->getY()) ;
// 	vars.push_back( &pb->getY()) ;
// 	vars.push_back( &pc->getX()) ;
// 	vars.push_back( &pd->getX()) ;
// 	vars.push_back( &pi->getX()) ;
// 	vars.push_back( &pi->getY()) ;
// 	
// 	std::vector<std::pair<double, double> > limits ;
// 	limits.push_back(std::make_pair(sample.getCenter().getY()-sample.height()/2.1, sample.getCenter().getY()+sample.height()/2.1)) ;
// 	limits.push_back(std::make_pair(sample.getCenter().getY()-sample.height()/2.1, sample.getCenter().getY()+sample.height()/2.1)) ;
// 	limits.push_back(std::make_pair(sample.getCenter().getX()-sample.width()/2.1, sample.getCenter().getX()+sample.width()/2.1)) ;
// 	limits.push_back(std::make_pair(sample.getCenter().getX()-sample.width()/2.1, sample.getCenter().getX()+sample.width()/2.1)) ;
// 	limits.push_back(std::make_pair(sample.getCenter().getX()-sample.width()/2, sample.getCenter().getX()+sample.width()/2)) ;
// 	limits.push_back(std::make_pair(sample.getCenter().getX()-sample.height()/2, sample.getCenter().getX()+sample.height()/2)) ;
// 	
// 	GeneticAlgorithmOptimizer ga( vars, limits, &distanceBetweenMeshes) ;
// 	ga.optimize(1e-12, 60, 150,  .1, .1) ;
// 
// 	pa->getY() = ga.getValues()[0].second ;
// 	pb->getY() = ga.getValues()[1].second ;
// 	pc->getX() = ga.getValues()[2].second ;
// 	pd->getX() = ga.getValues()[3].second ;
// 	pi->getX() = ga.getValues()[4].second ;
// 	pi->getY() = ga.getValues()[5].second ;
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
// 		std::cout << (*i)->getX() << "  " << (*i)->getY() << "  " << x[(*i)->getId()*2] << "  "<< x[(*i)->getId()*2+1] << "  " << dispx[indexj] << "  " << dispy[indexj] << std::endl;
// 		indexj++ ;
// 	}
// 	
// }
bool go = true ;
int countit = 0 ;
std::vector<double> ang ;
std::vector<double> sig ;
std::vector<double> eps ;
void step()
{
	
	size_t max_growth_steps = 20;
	size_t max_limit = 2000 ;
	int limit = 0 ;
	
	for(int s = 0 ; s < max_growth_steps ; s++)
	{
		go = featureTree->step() ;
		featureTree->setDeltaTime(.1);
		scales.push_back(imposeddisptop->getScale());

		if(go)
		{
// 			imposeddispright->setData(1e1);
			imposeddisptop->setData(1e1);
		}

		double da = 0 ;
		
		std::vector<DelaunayTriangle *> triangles = featureTree->getElements2D() ;

        Vector x = featureTree->getDisplacements() ;

        std::cout << "unknowns :" << x.size() << std::endl ;


        int npoints = triangles[0]->getBoundingPoints().size() ;

        std::string filename( "triangles" ) ;

        if( !featureTree->solverConverged() )
            filename = std::string( "failed-triangles" ) ;

        filename.append( itoa( s, 10 ) ) ;
        std::cout << filename << std::endl ;

        TriangleWriter writer(filename.c_str(), featureTree) ;
//      writer.reset( featureTree ) ;
//      writer.getField(TWFT_STRAIN_AND_STRESS ) ;
        writer.getField(TWFT_PRINCIPAL_STRESS ) ;
        writer.getField(TWFT_PRINCIPAL_STRAIN ) ;
        writer.getField( TWFT_VON_MISES ) ;
        writer.getField( TWFT_STIFFNESS ) ;
        writer.getField( TWFT_DAMAGE ) ;

        writer.write() ;

        double volume = 0 ;
        double xavg = 0 ;
        double yavg = 0 ;
        for(size_t k = 0 ; k < triangles.size() ; k++)
        {
            if(triangles[k]->getBehaviour()->type != VOID_BEHAVIOUR )
            {       
                double ar = triangles[k]->area() ;
                volume += ar ;
                for(size_t l = 0 ; l < npoints ;l++)
                {
                    xavg += x[triangles[k]->getBoundingPoint(l).getId()*2]*ar/npoints ;
                    yavg += x[triangles[k]->getBoundingPoint(l).getId()*2+1]*ar/npoints ;
                }

            }
        }
            
        xavg /= volume ;
        yavg /= volume ;
        std::pair<Vector, Vector> stempm = featureTree->getFieldMinMax(REAL_STRESS_FIELD) ;
        std::pair<Vector, Vector> etempm = featureTree->getFieldMinMax(STRAIN_FIELD) ;
        std::pair<Vector, Vector> vmm = featureTree->getFieldMinMax(VON_MISES_REAL_STRESS_FIELD) ;
        Vector stemp = featureTree->getAverageField(REAL_STRESS_FIELD) ;
        Vector etemp = featureTree->getAverageField(STRAIN_FIELD) ;
        
        std::cout << std::endl ;
        std::cout << "max value :" << x.max() << std::endl ;
        std::cout << "min value :" << x.min() << std::endl ;
        std::cout << "avg x value :" << xavg << std::endl ;
        std::cout << "avg y value :" << xavg << std::endl ;

        std::cout << "max sigma11 :" << stempm.second[0]  << std::endl ;
        std::cout << "min sigma11 :" << stempm.first[0]   << std::endl ;
        std::cout << "max sigma12 :" << stempm.second[2]  << std::endl ;
        std::cout << "min sigma12 :" << stempm.first[2]   << std::endl ;
        std::cout << "max sigma22 :" << stempm.second[1]  << std::endl ;
        std::cout << "min sigma22 :" << stempm.first[1]   << std::endl ;
        
        std::cout << "max epsilon11 :" << etempm.second[0] << std::endl ;
        std::cout << "min epsilon11 :" << etempm.first[0]  << std::endl ;
        std::cout << "max epsilon12 :" << etempm.second[2] << std::endl ;
        std::cout << "min epsilon12 :" << etempm.first[2]  << std::endl ;
        std::cout << "max epsilon22 :" << etempm.second[1] << std::endl ;
        std::cout << "min epsilon22 :" << etempm.first[1]  << std::endl ;
        
        std::cout << "max von Mises :" << vmm.second[0] << std::endl ;
        std::cout << "min von Mises :" << vmm.first[0] << std::endl ;
        
        std::cout << "average sigma11 : " << stemp[0] << std::endl ;
        std::cout << "average sigma22 : " << stemp[1] << std::endl ;
        std::cout << "average sigma12 : " << stemp[2] << std::endl ;
        std::cout << "average epsilon11 : " << etemp[0] << std::endl ;
        std::cout << "average epsilon22 : " << etemp[1] << std::endl ;
        std::cout << "average epsilon12 : " << etemp[2] << std::endl ;
        
        std::cout << std::endl ;;
		
		writer.reset(featureTree) ;
		writer.getField(TWFT_PRINCIPAL_STRESS ) ;
		writer.getField(TWFT_PRINCIPAL_STRAIN ) ;
		writer.getField(TWFT_CRITERION) ;
		writer.getField(TWFT_STIFFNESS) ;
		writer.getField(TWFT_DAMAGE) ;
		writer.append() ;
	}
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
	Matrix m0_paste=Material::cauchyGreen(std::make_pair(E_paste,nu), true,SPACE_TWO_DIMENSIONAL) ;

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
	
	Matrix m0_steely(3,3) ;
	m0_steely[0][0] =E_paste/(1.-2.*nu*nu) ; m0_steely[0][1] = nu*sqrt(E_stiff*E_paste)/(1.-2.*nu*nu) ; m0_steely[0][2] = 0 ; 
	m0_steely[1][0] = nu*sqrt(E_stiff*E_paste)/(1.-2.*nu*nu) ; m0_steely[1][1] =  E_stiff/(1.-2.*nu*nu) ; m0_steely[1][2] = 0 ; 
	m0_steely[2][0] = 0 ; m0_steely[2][1] = 0 ; m0_steely[2][2] = 0.25*(E_paste+E_stiff-2.*nu*sqrt(E_stiff*E_paste))/(1.-2.*nu*nu) ; 

	Matrix m0_steelx(3,3) ;
	m0_steelx[0][0] = E_stiff/(1.-2.*nu*nu) ; m0_steelx[0][1] = nu*sqrt(E_stiff*E_paste)/(1.-2.*nu*nu) ; m0_steelx[0][2] = 0 ; 
	m0_steelx[1][0] = nu*sqrt(E_stiff*E_paste)/(1.-2.*nu*nu) ; m0_steelx[1][1] = E_paste/(1.-2.*nu*nu) ; m0_steelx[1][2] = 0 ; 
	m0_steelx[2][0] = 0 ; m0_steelx[2][1] = 0 ; m0_steelx[2][2] = 0.25*(E_paste+E_stiff-2.*nu*sqrt(E_stiff*E_paste))/(1.-2.*nu*nu) ; 
	
	double nu_concreteSteel = nu ;
	double nu_steelConcrete = nu* E_paste/E_stiff ;
	
	Matrix complianceSteelx(3,3) ;
	complianceSteelx[0][0] = 1./E_stiff ;                complianceSteelx[0][1] = -nu_steelConcrete/E_paste ;
	complianceSteelx[1][0] = -nu_concreteSteel/E_stiff ; complianceSteelx[1][1] =  1./E_paste ;
	                                                                                                         complianceSteelx[2][2] =  E_paste/(1.-nu*nu)*(1.-nu)*.5 ;
	
	m0_steelx = inverse3x3Matrix(complianceSteelx) ;
	
	// Material behaviour for the "very" soft inclusion
	Matrix m0_soft(3,3) ;
	m0_soft[0][0] = E_soft/(1.-nu*nu) ; m0_soft[0][1] =E_soft/(1.-nu*nu)*nu ; m0_soft[0][2] = 0 ;
	m0_soft[1][0] = E_soft/(1.-nu*nu)*nu ; m0_soft[1][1] = E_soft/(1.-nu*nu) ; m0_soft[1][2] = 0 ; 
	m0_soft[2][0] = 0 ; m0_soft[2][1] = 0 ; m0_soft[2][2] = E_soft/(1.-nu*nu)*(1.-nu)/2. ; 

	Matrix d(3,3) ;
	d[0][0] = .1*E_paste ;
	d[1][1] = .1*E_paste ;
	d[2][2] = .1*E_paste ;
	FeatureTree F(&sample) ;
	featureTree = &F ;

	double cradius = 1*4 ;
	double mradius = 0.5 ;
	IsotropicLinearDamage * dfunc = new IsotropicLinearDamage() ;
	
	PseudoPlastic * psp = new PseudoPlastic(m0_paste, E_paste, 20, mradius) ;
	StiffnessAndFracture * saf = new StiffnessAndFracture(m0_paste, new NonLocalVonMises(20,E_paste, mradius), /*new NonLocal*//*IsotropicLinearDamage()*/new PlasticStrain()) ; 
	saf->criterion->setMaterialCharacteristicRadius(mradius);
	Stiffness * sf = new Stiffness(m0_steelx) ;

// 	sample.setBehaviour(sf) ;	
	sample.setBehaviour(saf) ;
// 	sample.setBehaviour(psp) ;
// 	sample.setBehaviour(new OrthothropicStiffness(E_paste, E_paste*.5, E_paste*.75, 0., 0)) ;
// 	F.addFeature(&sample, new Pore(2, -7,2) );
// 	F.addFeature(&sample, new Pore(2, 7,-2) );


	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA , BOTTOM)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI , BOTTOM_LEFT)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI , RIGHT)) ;
	F.addBoundaryCondition(imposeddispright) ;
	F.addBoundaryCondition(imposeddisptop) ;

	samplingnumber = atof(argv[1]);
	F.setSamplingNumber(samplingnumber) ;
	F.setOrder(LINEAR) ;
	F.setMaxIterationsPerStep(800) ;
	F.setDeltaTime(0.1);

	
// 	delete dt ;
// 	delete saf ;
	
	return 0 ;
}
