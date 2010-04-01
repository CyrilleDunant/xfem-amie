// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
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
#include "../physics/kelvinvoight.h"
#include "../physics/fracturecriteria/vonmises.h"
#include "../physics/stiffness.h"
#include "../physics/spatially_distributed_stiffness.h"
#include "../features/pore.h"
#include "../features/sample3d.h"
#include "../features/inclusion.h"
#include "../features/expansiveZone.h"
#include "../features/crack.h"
#include "../features/enrichmentInclusion.h"
#include "../features/expansiveZone3d.h"
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
std::vector<DelaunayTetrahedron *> tets ;

double E_min = 10;
double E_max = 0;

double x_max = 0 ;
double y_max = 0 ;

double x_min = 0 ;
double y_min = 0 ;

Vector b(0) ;
Vector x(0) ;
Vector sigma(0) ; 
Vector sigma11(0) ; 
Vector sigma22(0) ; 
Vector sigma33(0) ; 
Vector sigma12(0) ; 
Vector sigma13(0) ; 
Vector sigma23(0) ; 
Vector epsilon(0) ; 
Vector epsilon11(0) ; 
Vector epsilon22(0) ; 
Vector epsilon33(0) ; 
Vector epsilon12(0) ; 
Vector epsilon13(0) ; 
Vector epsilon23(0) ; 
Vector vonMises(0) ; 
Vector stiffness(0) ; 
Vector angle(0) ; 
Vector damage(0) ; 

double aggregateVolume = 0;

std::vector<double> energy ;
double percent = 0.7 ;

Matrix makeStiffnessMatrix(Vector E, double nu)
{
	Matrix m(6,6) ;
	m[0][0] = 1. - nu ; m[0][1] = nu ; m[0][2] = nu ;
	m[1][0] = nu ; m[1][1] = 1. - nu ; m[1][2] = nu ;
	m[2][0] = nu ; m[2][1] = nu ; m[2][2] = 1. - nu ;
	m[3][3] = 0.5 - nu ;
	m[4][4] = 0.5 - nu ;
	m[5][5] = 0.5 - nu ;
	m *= 1/((1.+nu)*(1.-2.*nu)) ;

	for(int i = 0 ; i < 3 ; i++)
	{
		for(int j = 0 ; j < 3 ; j++)
			m[i][j] *= E[i] ;
		m[i+3][i+3] *= E[i] ;
	}

	return m ;
}

std::vector<std::pair<ExpansiveZone3D *, Inclusion3D *> > generateExpansiveZones(int n, std::vector<Inclusion3D * > & incs , FeatureTree & F)
{
	double E_csh = 31e9 ;
	double nu_csh = .28 ;
	
	Vector csh(3) ;
	csh[0] = percent*E_csh ;
	csh[1] = percent*E_csh ;
	csh[2] = percent*E_csh ;
	
	Matrix m_csh = makeStiffnessMatrix(csh,nu_csh) ;
	
	std::vector<std::pair<ExpansiveZone3D *, Inclusion3D *> > ret ;
	aggregateVolume = 0 ;
	double radius = 0.5 ;
	for(size_t i = 0 ; i < incs.size() ; i++)
	{
		aggregateVolume += incs[i]->volume() ;
		for(int j = 0 ; j < n ; j++)
		{
				
			Point pos((2.*rand()/RAND_MAX-1.),(2.*rand()/RAND_MAX-1.),(2.*rand()/RAND_MAX-1.)) ;
			pos /= pos.norm() ;
			pos *= (2.*rand()/RAND_MAX-1.)*(incs[i]->getRadius() - 3.) ;
			Point center = incs[i]->getCenter()+pos ; 
			
			bool alone  = true ;
			
			for(size_t k = 0 ; k < ret.size() ; k++ )
			{
				if (squareDist(center, ret[k].first->Sphere::getCenter()) < (radius*60.+radius*60.)*(radius*60.+radius*60.))
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
				a[2] = 0.5 ;
				
				ExpansiveZone3D * z = new ExpansiveZone3D(incs[i], radius, center.x, center.y, center.z, m_csh, a) ;
				ret.push_back(std::make_pair(z, incs[i])) ;
			}
		}
	}

	for(size_t i = 0 ; i < ret.size() ; i++)
	{
		ret[i].first->setRadius(radius) ;
		F.addFeature(ret[i].second, ret[i].first) ;
	}
	std::cout << "initial Reacted Area = " << 4/3*M_PI*radius*radius*radius*ret.size() << " in "<< ret.size() << " zones"<< std::endl ;
	std::cout << "Reactive aggregate Area = " << aggregateVolume << std::endl ;
	return ret ;	
}



double uni_directional_step(int dir)
{
	bool go_on = true ;
	int tries = 0 ;
	int maxtries = 200 ;

	while(go_on && tries < maxtries)
	{
		featureTree->step(0.00) ;
		go_on = featureTree->solverConverged() &&  (featureTree->meshChanged() || featureTree->enrichmentChanged());
		tries++ ;
	}

	std::cout << std::endl ;
	std::cout << std::endl ;
	std::cout << std::endl ;
	std::cout << std::endl ;
	std::cout << tries << " iterations before convergence" << std::endl ;
	std::cout << std::endl ;
	std::cout << std::endl ;
	std::cout << std::endl ;
	std::cout << std::endl ;

	tets= featureTree->getTetrahedrons() ;
	x.resize(featureTree->getDisplacements().size()) ;
	x = featureTree->getDisplacements() ;

	std::pair<Vector, Vector > sigma_epsilon ;
	sigma_epsilon.first.resize(24*tets.size()) ;
	sigma_epsilon.second.resize(24*tets.size()) ;
	sigma_epsilon = featureTree->getStressAndStrain(tets) ;
	sigma.resize(sigma_epsilon.first.size()) ;
	sigma = sigma_epsilon.first ;
	epsilon.resize(sigma_epsilon.second.size()) ;
	epsilon = sigma_epsilon.second ;
	sigma11.resize(sigma.size()/6, 0.) ;
	sigma22.resize(sigma.size()/6, 0.) ;
	sigma33.resize(sigma.size()/6, 0.) ;
	sigma12.resize(sigma.size()/6, 0.) ;
	sigma13.resize(sigma.size()/6, 0.) ;
	sigma23.resize(sigma.size()/6, 0.) ;
	
	epsilon11.resize(sigma.size()/6, 0.) ;
	epsilon22.resize(sigma.size()/6, 0.) ;
	epsilon33.resize(sigma.size()/6, 0.) ;
	epsilon12.resize(sigma.size()/6, 0.) ;
	epsilon13.resize(sigma.size()/6, 0.) ;
	epsilon23.resize(sigma.size()/6, 0.) ;
	stiffness.resize(sigma.size()/6, 0.) ;
	vonMises.resize(sigma.size()/6, 0.) ;
	angle.resize(sigma.size()/6, 0.) ;
	damage.resize(sigma.size()/6, 0.) ;

	std::cout << "unknowns :" << x.size() << std::endl ;

	
	int npoints = 4 ;
	
	double volume = 0 ;
	double avg_e_xx = 0;
	double avg_e_yy = 0;
	double avg_e_zz = 0;
	double avg_e_xy = 0;
	double avg_e_xz = 0;
	double avg_e_yz = 0;
	double avg_s_xx = 0;
	double avg_s_zz = 0;
	double avg_s_yy = 0;
	double avg_s_xy = 0;
	double avg_s_xz = 0;
	double avg_s_yz = 0;
	double e_xx = 0 ;
	double ex_count = 0 ;
	double xavg = 0 ;
	
	double n_void = 0 ;
	for(size_t k = 0 ; k < tets.size() ; k++)
	{
			
		if(tets[k]->getBehaviour()->type == VOID_BEHAVIOUR )
			n_void++ ;
		if(tets[k]->getBehaviour()->type != VOID_BEHAVIOUR )
		{
			volume += tets[k]->volume() ;
			if(tets[k]->getBehaviour()->type != VOID_BEHAVIOUR)
			{
				if(tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0] > E_max)
					E_max = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0] ;
				if(tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0] < E_min)
					E_min = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0] ;
				
				stiffness[k*npoints] = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0] ;
				stiffness[k*npoints+1] = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0] ;
				stiffness[k*npoints+2] = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0] ;
				stiffness[k*npoints+3] = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0] ;
				damage[k*npoints] = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0]/tets[k]->getBehaviour()->param[0][0] ;
				damage[k*npoints+1] = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0]/tets[k]->getBehaviour()->param[0][0] ;
				damage[k*npoints+2] = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0]/tets[k]->getBehaviour()->param[0][0] ;
				damage[k*npoints+3] = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0]/tets[k]->getBehaviour()->param[0][0] ;
			}
				
			sigma11[k*npoints] = sigma[k*npoints*6];
			sigma22[k*npoints] = sigma[k*npoints*6+1];
			sigma33[k*npoints] = sigma[k*npoints*6+2];
			sigma12[k*npoints] = sigma[k*npoints*6+3];
			sigma13[k*npoints] = sigma[k*npoints*6+4];
			sigma23[k*npoints] = sigma[k*npoints*6+5];
			
			sigma11[k*npoints+1] = sigma[k*npoints*6+6];
			sigma22[k*npoints+1] = sigma[k*npoints*6+7];
			sigma33[k*npoints+1] = sigma[k*npoints*6+8];
			sigma12[k*npoints+1] = sigma[k*npoints*6+9];
			sigma13[k*npoints+1] = sigma[k*npoints*6+10];
			sigma23[k*npoints+1] = sigma[k*npoints*6+11];
			
			sigma11[k*npoints+2] = sigma[k*npoints*6+12];
			sigma22[k*npoints+2] = sigma[k*npoints*6+13];
			sigma33[k*npoints+2] = sigma[k*npoints*6+14];
			sigma12[k*npoints+2] = sigma[k*npoints*6+15];
			sigma13[k*npoints+2] = sigma[k*npoints*6+16];
			sigma23[k*npoints+2] = sigma[k*npoints*6+17];
			
			sigma11[k*npoints+3] = sigma[k*npoints*6+18];
			sigma22[k*npoints+3] = sigma[k*npoints*6+19];
			sigma33[k*npoints+3] = sigma[k*npoints*6+20];
			sigma12[k*npoints+3] = sigma[k*npoints*6+21];
			sigma13[k*npoints+3] = sigma[k*npoints*6+22];
			sigma23[k*npoints+3] = sigma[k*npoints*6+23];
			
			epsilon11[k*npoints] = epsilon[k*npoints*6];
			epsilon22[k*npoints] = epsilon[k*npoints*6+1];
			epsilon33[k*npoints] = epsilon[k*npoints*6+2];
			epsilon12[k*npoints] = epsilon[k*npoints*6+3];
			epsilon13[k*npoints] = epsilon[k*npoints*6+4];
			epsilon23[k*npoints] = epsilon[k*npoints*6+5];
			
			epsilon11[k*npoints+1] = epsilon[k*npoints*6+6];
			epsilon22[k*npoints+1] = epsilon[k*npoints*6+7];
			epsilon33[k*npoints+1] = epsilon[k*npoints*6+8];
			epsilon12[k*npoints+1] = epsilon[k*npoints*6+9];
			epsilon13[k*npoints+1] = epsilon[k*npoints*6+10];
			epsilon23[k*npoints+1] = epsilon[k*npoints*6+11];
			
			epsilon11[k*npoints+2] = epsilon[k*npoints*6+12];
			epsilon22[k*npoints+2] = epsilon[k*npoints*6+13];
			epsilon33[k*npoints+2] = epsilon[k*npoints*6+14];
			epsilon12[k*npoints+2] = epsilon[k*npoints*6+15];
			epsilon13[k*npoints+2] = epsilon[k*npoints*6+16];
			epsilon23[k*npoints+2] = epsilon[k*npoints*6+17];
			
			epsilon11[k*npoints+3] = epsilon[k*npoints*6+18];
			epsilon22[k*npoints+3] = epsilon[k*npoints*6+19];
			epsilon33[k*npoints+3] = epsilon[k*npoints*6+20];
			epsilon12[k*npoints+3] = epsilon[k*npoints*6+21];
			epsilon13[k*npoints+3] = epsilon[k*npoints*6+22];
			epsilon23[k*npoints+3] = epsilon[k*npoints*6+23];
			
			double vm0 = 0 ;
			double agl = 0 ;
			if(tets[k]->getBehaviour()->type != VOID_BEHAVIOUR)
				vm0 = tets[k]->getState().getMaximumVonMisesStress() ;
			if(tets[k]->getBehaviour()->type != VOID_BEHAVIOUR)
				agl = tets[k]->getState().getPrincipalAngle(tets[k]->getCenter()) ;
			for(size_t l = 0 ; l < 4 ; l++)
			{
				vonMises[k*npoints+l]  = vm0 ;
				angle[k*npoints+l]  = agl ;
			}
			
			double ar = tets[k]->volume() ;
			for(size_t l = 0 ; l < npoints ;l++)
			{
				avg_e_xx += (epsilon11[k*npoints+l]/npoints)*ar;
				avg_e_yy += (epsilon22[k*npoints+l]/npoints)*ar;
				avg_e_zz += (epsilon33[k*npoints+l]/npoints)*ar;
				avg_e_xy += (epsilon12[k*npoints+l]/npoints)*ar;
				avg_e_xz += (epsilon13[k*npoints+l]/npoints)*ar;
				avg_e_yz += (epsilon23[k*npoints+l]/npoints)*ar;
				avg_s_xx += (sigma11[k*npoints+l]/npoints)*ar;
				avg_s_yy += (sigma22[k*npoints+l]/npoints)*ar;
				avg_s_zz += (sigma33[k*npoints+l]/npoints)*ar;
				avg_s_xy += (sigma12[k*npoints+l]/npoints)*ar;
				avg_s_xz += (sigma13[k*npoints+l]/npoints)*ar;
				avg_s_yz += (sigma23[k*npoints+l]/npoints)*ar;
				xavg += x[tets[k]->getBoundingPoint(l).id]*ar/npoints ;
			}

		}
	}
	
	switch(dir)
	{
		case 0:
			return avg_s_xx/avg_e_xx ;
		case 1:
			return avg_s_yy/avg_e_yy ;
		case 2:
			return avg_s_zz/avg_e_zz ;
	}
	std::cout << "unknwon direction" << std::endl ;
	return -1 ;
}

Vector tri_directional_step(double box_dim, Matrix m_mat, Matrix m_inc, std::vector<double> radius)
{
	Vector Young(3) ;

	Sample3D box(NULL, box_dim , box_dim, box_dim, 0, 0, 0) ;
	box.setBehaviour(new Stiffness(m_mat)) ;

	std::vector<Inclusion3D *> inc ;
	for(size_t i = 0 ; i < radius.size() ; i++)
	{
		inc.push_back(new Inclusion3D(radius[i],0,0,0)) ;
		inc[i]->setBehaviour(new Stiffness(m_inc)) ;
	}

	std::vector<Feature *> feats ;
	for(size_t i = 0; i < inc.size() ; i++)
		feats.push_back(inc[i]) ;
	inc.clear() ;

	int n_agg = (int) radius.size() ;

	feats = placement(box.getPrimitive(), feats, &n_agg, n_agg * 1000);

	for(size_t i = 0; i < feats.size() ; i++)
		inc.push_back(static_cast<Inclusion3D *>(feats[i])) ;
	feats.clear() ;

	int n_zones = 10 * (int) box_dim ;

	for(int dir = 0 ; dir < 3 ; dir++)
	{
		FeatureTree F(&box) ;
		featureTree = &F ;
		for(size_t i = 0 ; i < inc.size() ; i++)
			F.addFeature(&box,inc[i]) ;

		std::vector<std::pair<ExpansiveZone3D *, Inclusion3D *> > zones = generateExpansiveZones(n_zones, inc , F) ;

		switch(dir)
		{
			case 0:
			{
				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_XI, LEFT, 0)) ;
				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, LEFT)) ;
				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, LEFT)) ;
				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_XI, RIGHT, 100)) ;
				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, RIGHT)) ;
				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, RIGHT)) ;
				break ;
			}
			case 1:
			{
				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, BOTTOM)) ;
				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA, BOTTOM,0)) ;
				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, BOTTOM)) ;
				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, TOP)) ;
				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA, TOP,100)) ;
				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, TOP)) ;
				break ;
			}
			case 2:
			{
				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, FRONT)) ;
				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, FRONT)) ;
				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ZETA, FRONT,0)) ;
				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, BACK)) ;
				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BACK)) ;
				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ZETA, BACK,100)) ;
				break ;
			}
		}
	
		F.sample(256) ;
	
		F.setOrder(LINEAR) ;
		F.generateElements() ;
	
		Young[dir] = uni_directional_step(dir) ;
		std::cout << Young[dir] << std::endl ;

	}
	
	return Young ;	
}

int main(int argc, char *argv[])
{
	double agg_fraction = 0.2 ;

	std::vector<std::vector<double> > radius ;
	std::vector<double> r0 ;
	std::vector<double> r1 ;
	for(int i = 1 ; i < 5 ; i++)
	{
		for(int j = 0 ; j < i ; j++)
			r0.push_back(100 / (i*i)) ;
	}
	r1.push_back(250) ;
	radius.push_back(r0) ;
	radius.push_back(r1) ;

	double agg_volume = 0 ;
	std::vector<double> agg_volume_by_class ;
	for(size_t i = 0 ; i < radius.size() ; i++)
	{
		agg_volume = 0 ;
		for(size_t j = 0 ; j < radius[i].size() ; j++)
			agg_volume += radius[i][j]*radius[i][j]*radius[i][j] ;
		agg_volume *= 4/3 * M_PI ;
		agg_volume_by_class.push_back(agg_volume) ;
	}
	agg_volume = 0 ;
	for(size_t i = 0 ; i < agg_volume_by_class.size() ; i++)
		agg_volume += agg_volume_by_class[i] ;

	double total_volume = agg_volume / agg_fraction ;
	double cement_volume = total_volume - agg_volume ;

	Vector cement(3) ;
	cement[0] = 25*1e9 ;
	cement[1] = 25*1e9 ;
	cement[2] = 25*1e9 ;

	Vector aggregates(3) ;
	aggregates[0] = 70*1e9 ;
	aggregates[1] = 70*1e9 ;
	aggregates[2] = 70*1e9 ;

	double iter_volume = cement_volume ;

	std::vector<Vector> homogenized ;
	homogenized.push_back(cement) ;

	for(size_t i = 0 ; i < radius.size() ; i++)
	{
		iter_volume += agg_volume_by_class[i] ;
		double box_dim = std::pow(iter_volume, 0.33333333) ;
		std::cout << box_dim << std::endl ;

		Vector hom = tri_directional_step(box_dim,makeStiffnessMatrix(homogenized.back(),0.2),makeStiffnessMatrix(aggregates,0.2),radius[i]) ;

		homogenized.push_back(hom) ;
	}
	
	return 0 ;
}
