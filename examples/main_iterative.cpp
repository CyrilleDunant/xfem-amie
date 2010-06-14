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
#include "../physics/homogenization/elastic_homogenization.h"
#include "../physics/homogenization/elastic_bounds.h"
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

	double E_ = E[0]+E[1]+E[2] ;
	E_ /= 3 ;

	for(int i = 0 ; i < 3 ; i++)
	{
		for(int j = 0 ; j < 3 ; j++)
			m[i][j] *= E_ ;
		m[i+3][i+3] *= E_ ;
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



std::pair<double, double> uni_directional_step(int dir)
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

	std::cout << tries << " iterations before convergence" << std::endl ;

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
	
	double Exx = avg_s_xx/avg_e_xx ;
	double Eyy = avg_s_yy/avg_e_yy ;
	double Ezz = avg_s_zz/avg_e_zz ;

	double nuxx = std::sqrt(avg_e_yy*avg_e_yy+avg_e_zz*avg_e_zz)/(std::sqrt(2)*(avg_e_xx)) ;
	double nuyy = std::sqrt(avg_e_xx*avg_e_xx+avg_e_zz*avg_e_zz)/(std::sqrt(2)*(avg_e_yy)) ;
	double nuzz = std::sqrt(avg_e_xx*avg_e_xx+avg_e_yy*avg_e_yy)/(std::sqrt(2)*(avg_e_zz)) ;

	if(nuxx > 1 || nuyy > 1 || nuzz > 1)
	{
		std::cout << std::endl ;
		std::cout << "sigma" << std::endl ;
		std::cout << avg_s_xx << std::endl ;
		std::cout << avg_s_yy << std::endl ;
		std::cout << avg_s_zz << std::endl ;
		std::cout << "epsilon" << std::endl ;
		std::cout << avg_e_xx << std::endl ;
		std::cout << avg_e_yy << std::endl ;
		std::cout << avg_e_zz << std::endl ;
		std::cout << std::endl ;
	}

	std::pair<double, double> Enu ;
	Enu.first = -1 ;
	Enu.second = -1 ;

	switch(dir)
	{
		case 0:
		{
			Enu = std::make_pair(Exx,nuxx) ;
			break ;
		}
		case 1:
		{
			Enu = std::make_pair(Eyy,nuyy) ;
			break ;
		}
		case 2:
		{
			Enu = std::make_pair(Ezz,nuzz) ;
			break ;
		}
	}
	return Enu ;
}

Matrix tri_directional_step(double box_dim, Matrix m_mat, Matrix m_inc, double radius, size_t n_agg, double disp)
{
	std::pair<double,double> Enu_xx ;
	std::pair<double,double> Enu_yy ;
	std::pair<double,double> Enu_zz ;

	Sample3D box(NULL, box_dim , box_dim, box_dim, 0, 0, 0) ;
	box.setBehaviour(new Stiffness(m_mat)) ;

	std::vector<Inclusion3D *> inc ;
	for(size_t i = 0 ; i < n_agg ; i++)
	{
		inc.push_back(new Inclusion3D(radius,0,0,0)) ;
		inc[i]->setBehaviour(new Stiffness(m_inc)) ;
	}

	std::vector<Feature *> feats ;
	for(size_t i = 0; i < inc.size() ; i++)
		feats.push_back(inc[i]) ;
	inc.clear() ;

	int nAgg = n_agg ;

	feats = placement(box.getPrimitive(), feats, &nAgg, nAgg * 1000);

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

//		std::vector<std::pair<ExpansiveZone3D *, Inclusion3D *> > zones = generateExpansiveZones(n_zones, inc , F) ;
		F.resetBoundaryConditions() ;
		switch(dir)
		{
			case 0:
			{
				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_XI, LEFT, 0)) ;
//				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, LEFT)) ;
//				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, LEFT)) ;
				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_XI, RIGHT, disp)) ;
//				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, RIGHT)) ;
//				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, RIGHT)) ;
				break ;
			}
			case 1:
			{
//				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, BOTTOM)) ;
				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA, BOTTOM, 0)) ;
//				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, BOTTOM)) ;
//				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, TOP)) ;
				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA, TOP, disp)) ;
//				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, TOP)) ;
				break ;
			}
			case 2:
			{
//				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, BACK)) ;
//				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BACK)) ;
				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ZETA, BACK, 0)) ;
//				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, FRONT)) ;
//				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, FRONT)) ;
				F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ZETA, FRONT, disp)) ;
				break ;
			}
		}
	
		F.sample(16000) ;
	
		F.setOrder(LINEAR) ;
		F.generateElements() ;
	
		std::pair<double,double> Enu = uni_directional_step(dir) ;
		std::cout << Enu.first << ";" << Enu.second << std::endl ;
		switch(dir)
		{
			case 0 :
			{
				Enu_xx = Enu ;
				break ;
			}
			case 1 :
			{
				Enu_yy = Enu ;
				break ;
			}
			case 2 :
			{
				Enu_zz = Enu ;
				break ;
			}
		}

	}
	double nu = Enu_xx.second + Enu_yy.second + Enu_zz.second ;
	nu /= 3 ;

	std::cout << std::endl ;
	std::cout << std::endl ;
	std::cout << std::endl ;
	std::cout << (Enu_xx.first + Enu_yy.first + Enu_zz.first)/3 << std::endl ;
	std::cout << nu << std::endl ;
	std::cout << std::endl ;
	std::cout << std::endl ;
	std::cout << std::endl ;

	Vector Young(3) ;
	Young[0] = Enu_xx.first ;
	Young[1] = Enu_yy.first ;
	Young[2] = Enu_zz.first ;
	
	return makeStiffnessMatrix(Young, nu) ;	
}

int main(int argc, char *argv[])
{
	double agg_fraction = 0.7 ;

	std::vector<std::pair<double, double > > radius ;
	radius.push_back(std::make_pair(53, 3.5)) ;
	radius.push_back(std::make_pair(63, 4.7-3.5)) ;
	radius.push_back(std::make_pair(100, 7.1-4.7)) ;
	radius.push_back(std::make_pair(160, 10.0-7.1)) ;
	radius.push_back(std::make_pair(250, 14.1-10.0)) ;
	radius.push_back(std::make_pair(400, 22.1-14.1)) ;
	radius.push_back(std::make_pair(630, 34.1-22.1)) ;
	radius.push_back(std::make_pair(1000, 49.0-34.1)) ;
	radius.push_back(std::make_pair(1600, 64.2-49.0)) ;
	radius.push_back(std::make_pair(2500, 80.5-64.2)) ;
	radius.push_back(std::make_pair(4000, 97.8-80.5)) ;
	radius.push_back(std::make_pair(5000, 100.0-97.8)) ;

	double agg_mass = 0 ;
	for(size_t i = 0 ; i < radius.size() ; i++)
		agg_mass += radius[i].second ;
	double total_mass = agg_mass / agg_fraction ;
	double cem_mass = total_mass - agg_mass ;

	double agg_volume = agg_mass/2.2 ;
	double cement_volume = cem_mass/3.1 ;

	double total_volume = agg_volume + cement_volume ;
	double vfraction = agg_volume / total_volume ;
	std::cout << vfraction << std::endl ;
//	double cement_volume = total_volume - agg_volume ;

	Vector cement(3) ;
	cement[0] = 25 ;
	cement[1] = 25 ;
	cement[2] = 25 ;

	Matrix stiff_cem = makeStiffnessMatrix(cement, 0.2) ;
//	SimpleMaterial sm_cem(25, 0.2) ;

	Vector aggregates(3) ;
	aggregates[0] = 70 ;
	aggregates[1] = 70 ;
	aggregates[2] = 70 ;

	Matrix stiff_agg = makeStiffnessMatrix(aggregates, 0.2) ;
//	SimpleMaterial sm_agg(70, 0.2) ;

	double iter_volume = cement_volume ;

	std::vector<Matrix> stiff_hom ;
	stiff_hom.push_back(stiff_cem) ;

/*	for(size_t i = 0 ; i < radius.size() ; i++)
	{

		std::cout << std::endl ;
		std::cout << std::endl ;
		std::cout << std::endl ;
		std::cout << i << std::endl ;
		std::cout << std::endl ;
		std::cout << std::endl ;
		std::cout << std::endl ;

		iter_volume += radius[i].second/2.2 ;
		double box_dim = std::pow(iter_volume, 0.33333333) ;
		std::cout << box_dim << std::endl ;

		size_t n_agg = 150  ;
		double this_agg_volume = radius[i].first*radius[i].first*radius[i].first ;
		this_agg_volume *= 4/3 * M_PI * n_agg ;
		double this_iter_volume = iter_volume * this_agg_volume * 2.2 / radius[i].second ;

//		std::cout << std::endl ;
//		std::cout << radius[i].second/iter_volume << std::endl ;
//		std::cout << this_agg_volume/this_iter_volume << std::endl ;

		box_dim = std::pow(this_iter_volume, 0.33333333) ;
		double effective_radius = radius[i].first ;

		while(box_dim > 1000)
		{
			if(box_dim / effective_radius > 10 || n_agg==1)
				effective_radius *= 0.75 ;
			else
			{
				n_agg *= 9 ;
				n_agg /= 10 ;
				if(n_agg == 0)
				{
					n_agg = 1 ;
				}
			}
			this_agg_volume = effective_radius*effective_radius*effective_radius ;
			this_agg_volume *= 4/3 * M_PI * n_agg ;
			this_iter_volume = iter_volume * this_agg_volume * 2.2 / radius[i].second ;
			box_dim = std::pow(this_iter_volume, 0.33333333) ;
		}

		double disp = std::min(effective_radius, box_dim/(10)) ;

		std::cout << box_dim << ";" << effective_radius << ";" << disp << std::endl ;
		std::cout << std::endl ;
		

		Matrix hom = tri_directional_step(box_dim,stiff_hom.back(),stiff_agg,effective_radius,n_agg, disp) ;

		stiff_hom.push_back(hom) ;
	}

	Matrix hom_final = stiff_hom.back() ;

	double nu_final = hom_final[0][1] / (hom_final[0][0] + hom_final[0][1]) ;
	double Exx_final = hom_final[0][0] * (1+nu_final)*(1-2*nu_final) /(1 - nu_final) ;
	double Eyy_final = hom_final[0][0] * (1+nu_final)*(1-2*nu_final) /(1 - nu_final) ;
	double Ezz_final = hom_final[0][0] * (1+nu_final)*(1-2*nu_final) /(1 - nu_final) ;

	std::cout << std::endl ;
	std::cout << std::endl ;
	std::cout << std::endl ;
	std::cout << std::endl ;
	std::cout << Exx_final << ";" << Eyy_final << ";" << Ezz_final << ";" << nu_final << std::endl ;
	std::cout << std::endl ;
	std::cout << std::endl ;
	std::cout << std::endl ;
	std::cout << std::endl ;*/


//	std::cout << std::endl ;
//	std::cout << "Incremental" << std::endl ;
//	SimpleMaterial hom_inc(INCREMENTAL, std::make_pair(vfraction, new SimpleMaterial(70,0.2)), new SimpleMaterial(25,0.2)) ;
//	std:: cout << hom_inc.getE() << ";" << hom_inc.getnu() << std::endl ;
//	std::cout << std::endl ;


/*	for(size_t i = 1 ; i < 99 ; i++) {
	Properties cem_Enu(HOOKE,std::make_pair(25,0.2)) ;
	Properties agg_Enu(HOOKE,std::make_pair(70,0.2)) ;

	Properties cem_frac(FRACTION,(double) i/99.) ;
	Properties agg_frac(FRACTION,1. - cem_frac.getValue(0)) ;

	Material m_cement(cem_frac) ;
	Material m_aggregates(agg_frac) ;

	std::pair<bool,Properties> cem_test = cem_Enu.convert(BULK_SHEAR) ;
	std::pair<bool,Properties> agg_test = agg_Enu.convert(BULK_SHEAR) ;
	
	if(cem_test.first)
	{
		m_cement.push_back(cem_test.second) ;
		if(agg_test.first)
		{
			m_aggregates.push_back(agg_test.second) ;
			std::vector<Material> mat ;
			mat.push_back(m_cement) ;
			mat.push_back(m_aggregates) ;
			std::pair<bool, Material> hom_mt = MoriTanaka().apply(mat) ;
			std::pair<bool, Properties> conv_u = hom_mt.second[0].convert(HOOKE) ;
//			std::pair<bool, Properties> conv_l = hom_mt.second[1].convert(HOOKE) ;
//			std::cout << std::endl ;
//			std::cout << std::endl ;
//			std::cout << "UPPER BOUNDS" << std::endl ;
			std::cout << conv_u.second.getValue(0) << /*";" << conv_u.second[0].getValue(1) <<*/ //std::endl ;
//			std::cout << std::endl ;
//			std::cout << std::endl ;
//			std::cout << std::endl ;
//			std::cout << std::endl ;
//			std::cout << "LOWER BOUNDS" << std::endl ;
//			std::cout << hom_mt.second[1].getValue(0) << ";" << hom_mt.second[1].getValue(1) << std::endl ;
//			std::cout << std::endl ;
//			std::cout << std::endl ;
//			Properties voigt("BULK_SHEAR",


/*		} else {
			std::cout << "agg shit" <<std::endl ;
		}
	} else {
		std::cout << "shit..." << std::endl ;
	}
	}*/

	return 0 ;
}
