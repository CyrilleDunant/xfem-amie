// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../utilities/samplingcriterion.h"
#include "../features/features.h"
#include "../physics/physics_base.h"
#include "../physics/kelvinvoight.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/laplacian.h"
#include "../physics/materials/aggregate_behaviour.h"
#include "../physics/stiffness_and_fracture.h"
#include "../physics/materials/paste_behaviour.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../features/pore.h"
#include "../features/sample.h"
#include "../features/sample3d.h"
#include "../features/inclusion.h"
#include "../geometry/geometry_3D.h"
#include "../features/inclusion3d.h"
#include "../features/expansiveZone.h"
#include "../features/crack.h"
#include "../features/enrichmentInclusion.h"
#include "../mesher/delaunay_3d.h"
#include "../features/expansiveZone3d.h"
#include "../solvers/assembly.h"
#include "../utilities/granulo.h"
#include "../utilities/itoa.h"
#include "../utilities/writer/triangle_writer.h"
#include "../utilities/random.h"
#include "../utilities/placement.h"
#include "../physics/stiffness.h"
#include "../physics/void_form.h"
#include "../utilities/writer/voxel_writer.h"
#include <sys/time.h>

#include <fstream>

#include <cmath>
#include <typeinfo>
#include <limits>
#include <GL/glut.h>
#include <time.h> 
#define DEBUG 

#define ID_QUIT 1
#define ID_ZOOM 2
#define ID_UNZOOM 3
#define ID_NEXT10 4
#define ID_NEXT100 5
#define ID_NEXT1000 6
#define ID_NEXT 7
#define ID_NEXT_TIME 8
#define ID_REFINE 9
#define ID_AMPLIFY 10
#define ID_DEAMPLIFY 11

#define ID_DISP 12
#define ID_STRAIN_XX 13
#define ID_STRAIN_XY 14
#define ID_STRAIN_XZ 15
#define ID_STRAIN_YZ 16
#define ID_STRAIN_YY 17
#define ID_STRAIN_ZZ 18
#define ID_STRESS_XX 19
#define ID_STRESS_XY 20
#define ID_STRESS_YY 21
#define ID_STRESS_ZZ 22
#define ID_STRESS_XZ 23
#define ID_STRESS_YZ 24

#define ID_STIFNESS 25
#define ID_ELEM 26
#define ID_VON_MISES 27
#define ID_ANGLE 28
#define ID_ENRICHMENT 29

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

double stress = 15e6 ;

bool firstRun = true ;

std::vector<DelaunayTriangle *> tris__ ;

std::pair<std::vector<Inclusion * >, std::vector<Pore * > > i_et_p ;

std::vector<std::pair<ExpansiveZone *, Inclusion *> > zones ;

std::vector<Inclusion *> inclusions ;

std::vector<std::pair<double, double> > expansion_reaction ;
std::vector<std::pair<double, double> > expansion_stress_xx ;
std::vector<std::pair<double, double> > expansion_stress_yy ;
std::vector<std::pair<double, double> > apparent_extension ;
std::vector<double> cracked_volume ;
std::vector<double> damaged_volume ;

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

VoidForm* black = new VoidForm() ;

int totit = 5 ;

double factor = 200 ;
MinimumAngle cri(M_PI/6.) ;
bool nothingToAdd = false ;
bool dlist = false ;
int count = 0 ;
double aggregateArea = 0;

double shape ;
double orientation ;
double spread ;

void step()
{
	bool go_on = featureTree->step() ;
	if(featureTree->solverConverged())
	{
		cracked_volume.push_back(featureTree->crackedVolume) ;
		damaged_volume.push_back(featureTree->damagedVolume) ;
	}

	triangles = featureTree->getElements2D() ;
	x.resize(featureTree->getDisplacements().size()) ;
	x = featureTree->getDisplacements() ;
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
	double e_xx_max = 0 ;
	double e_xx_min = 0 ;
        double e_yy_max = 0 ;
        double e_yy_min = 0 ;
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
			
/*			for(size_t l = 0 ; l < triangles[k]->getBoundingPoints().size() ; l++)
			{
				Vector vm0 = triangles[k]->getState().getPrincipalStresses(triangles[k]->getBoundingPoint(l)) ;
				vonMises[k*triangles[k]->getBoundingPoints().size()+l]  = sqrt(((vm0[0]-vm0[1])*(vm0[0]-vm0[1]))/2.) ;

				double agl = triangles[k]->getState().getPrincipalAngle(triangles[k]->getBoundingPoint(l))[0] ;
				agl = (vm0[0] <= 0 && vm0[1] <= 0)*180. ;
				angle[k*triangles[k]->getBoundingPoints().size()+l]  = agl ;
			}*/
			
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
	
        std::cout << "apparent extension X " << e_xx_max-e_xx_min << std::endl ;
        std::cout << "apparent extension Y " << e_yy_max-e_yy_min << std::endl ;

}

int main(int argc, char *argv[])
{

	double E = 35*1e9 ;
	double nu = 0.3 ;
	Matrix m0(3,3) ;
	m0[0][0] = E/(1.-nu*nu) ; m0[0][1] =E/(1.-nu*nu)*nu ; m0[0][2] = 0 ;
	m0[1][0] = E/(1.-nu*nu)*nu ; m0[1][1] = E/(1.-nu*nu) ; m0[1][2] = 0 ;
	m0[2][0] = 0 ; m0[2][1] = 0 ; m0[2][2] = E/(1.-nu*nu)*(1.-nu)/2. ;

	Sample box(NULL, 0.04,0.04,0.,0.) ;
	box.setBehaviour(new PasteBehaviour()) ;

	Inclusion* inc = new Inclusion(0.01, 0.003,0.003) ;
	inc->setBehaviour(new AggregateBehaviour()) ;
	

	double alpha = 0.1 ;
	Vector a(3) ;
	a[0] = alpha ; a[1] = alpha ; a[2] = 0. ;

	ExpansiveZone* exp1 = new ExpansiveZone(NULL, 0.0001, 0.003, 0.003, m0, a) ;
	ExpansiveZone* exp2 = new ExpansiveZone(NULL, 0.0001, 0.009, 0.0025, m0, a) ;
	ExpansiveZone* exp3 = new ExpansiveZone(NULL, 0.0001, -0.005, 0.003, m0, a) ;

	FeatureTree F(&box) ;


	featureTree = &F ;
	F.addFeature(&box, inc) ;
	F.addFeature(inc, exp1) ;
//	F.addFeature(inc, exp2) ;
//	F.addFeature(inc, exp3) ;
	F.setSamplingNumber(64) ;
	F.setDeltaTime(0.0001) ;
	F.setMaxIterationsPerStep(200) ;
	F.setOrder(LINEAR) ;


	TriangleWriter writer("test",NULL,0) ;

	for(size_t i = 0 ; i < 100 ; i++)
	{
	    F.resetBoundaryConditions() ;
	    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
	    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_XI, RIGHT, -0.001*i)) ;
	    step() ;

	    writer.reset(featureTree, 0);
	    writer.getField(TWFT_STIFFNESS) ;
	    writer.getField(TWFT_STRAIN) ;
	    writer.getField(TWFT_STRESS) ;
	    writer.append() ;

	}


	
	return 0 ;
}


