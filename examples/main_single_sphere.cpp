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
#include "../physics/laplacian.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/stiffness.h"
#include "../physics/stiffness_and_fracture.h"
#include "../features/pore.h"
#include "../features/sample.h"
#include "../features/sample3d.h"
#include "../features/inclusion.h"
#include "../features/inclusion3d.h"
#include "../features/expansiveZone.h"
#include "../features/crack.h"
#include "../features/enrichmentInclusion.h"
#include "../mesher/delaunay_3d.h"
#include "../utilities/writer/voxel_writer.h"
#include "../features/expansiveZone3d.h"
#include "../solvers/assembly.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"
#include <sys/time.h>

#include <fstream>

#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h> 
#define DEBUG 


using namespace Mu ;


FeatureTree * featureTree ;

void step()
{
	
  int nsteps = 1;// number of steps between two clicks on the opengl thing
	featureTree->setMaxIterationsPerStep(50) ;

	for(size_t i = 0 ; i < nsteps ; i++)
	{
		featureTree->step() ;
		
		std::vector<DelaunayTetrahedron *> tets= featureTree->getElements3D() ;
		Vector x = featureTree->getDisplacements() ;
	
		std::cout << "unknowns :" << x.size() << std::endl ;
	
		int npoints = 4 ;
	
		double volume = 0 ;

		double xavg = 0 ;
		
		for(size_t k = 0 ; k < tets.size() ; k++)
		{
			
			if(tets[k]->getBehaviour()->type != VOID_BEHAVIOUR )
			{
				double ar = tets[k]->volume() ;
				volume += ar ;
				for(size_t l = 0 ; l < npoints ;l++)
				{
					xavg += x[tets[k]->getBoundingPoint(l).id]*ar/npoints ;
				}
			}
		}
			
		xavg /= volume ;
		std::pair<Vector, Vector> stempm = featureTree->getFieldMinMax(REAL_STRESS_FIELD) ;
		std::pair<Vector, Vector> etempm = featureTree->getFieldMinMax(STRAIN_FIELD) ;
		std::pair<Vector, Vector> vmm = featureTree->getFieldMinMax(VON_MISES_REAL_STRESS_FIELD) ;
		Vector stemp = featureTree->getAverageField(REAL_STRESS_FIELD) ;
		Vector etemp = featureTree->getAverageField(STRAIN_FIELD) ;
		
		std::cout << std::endl ;
		std::cout << "max value :" << x.max() << std::endl ;
		std::cout << "min value :" << x.min() << std::endl ;
		std::cout << "avg value :" << xavg << std::endl ;

		std::cout << "max sigma11 :" << stempm.second[0]  << std::endl ;
		std::cout << "min sigma11 :" << stempm.first[0]   << std::endl ;
		std::cout << "max sigma12 :" << stempm.second[3]  << std::endl ;
		std::cout << "min sigma12 :" << stempm.first[3]   << std::endl ;
		std::cout << "max sigma13 :" << stempm.second[4]  << std::endl ;
		std::cout << "min sigma13 :" << stempm.first[4]   << std::endl ;
		std::cout << "max sigma22 :" << stempm.second[1]  << std::endl ;
		std::cout << "min sigma22 :" << stempm.first[1]   << std::endl ;
		std::cout << "max sigma23 :" << stempm.second[5]  << std::endl ;
		std::cout << "min sigma23 :" << stempm.first[5]   << std::endl ;
		std::cout << "max sigma33 :" << stempm.second[2]  << std::endl ;
		std::cout << "min sigma33 :" << stempm.first[2]   << std::endl ;
		
		std::cout << "max epsilon11 :" << etempm.second[0] << std::endl ;
		std::cout << "min epsilon11 :" << etempm.first[0]  << std::endl ;
		std::cout << "max epsilon12 :" << etempm.second[3] << std::endl ;
		std::cout << "min epsilon12 :" << etempm.first[3]  << std::endl ;
		std::cout << "max epsilon13 :" << etempm.second[4] << std::endl ;
		std::cout << "min epsilon13 :" << etempm.first[4]  << std::endl ;
		std::cout << "max epsilon22 :" << etempm.second[1] << std::endl ;
		std::cout << "min epsilon22 :" << etempm.first[1]  << std::endl ;
		std::cout << "max epsilon23 :" << etempm.second[5] << std::endl ;
		std::cout << "min epsilon23 :" << etempm.first[5]  << std::endl ;
		std::cout << "max epsilon33 :" << etempm.second[2] << std::endl ;
		std::cout << "min epsilon33 :" << etempm.first[2]  << std::endl ;
		
		std::cout << "max von Mises :" << vmm.second[0] << std::endl ;
		std::cout << "min von Mises :" << vmm.first[0] << std::endl ;
		
		std::cout << "average sigma11 : " << stemp[0] << std::endl ;
		std::cout << "average sigma22 : " << stemp[1] << std::endl ;
		std::cout << "average sigma33 : " << stemp[2] << std::endl ;
		std::cout << "average sigma12 : " << stemp[3] << std::endl ;
		std::cout << "average sigma13 : " << stemp[4] << std::endl ;
		std::cout << "average sigma23 : " << stemp[5] << std::endl ;
		std::cout << "average epsilon11 : " << etemp[0] << std::endl ;
		std::cout << "average epsilon22 : " << etemp[1] << std::endl ;
		std::cout << "average epsilon33 : " << etemp[2] << std::endl ;
		std::cout << "average epsilon12 : " << etemp[3] << std::endl ;
		std::cout << "average epsilon13 : " << etemp[4] << std::endl ;
		std::cout << "average epsilon23 : " << etemp[5] << std::endl ;
		
	}
// 	VoxelWriter vw1("sphere_stiffness", 100) ;
// 	vw1.getField(featureTree, VWFT_STIFFNESS) ;
// 	vw1.write();
	
	VoxelWriter vw("sphere_stress", 150) ;
	vw.getField(featureTree, VWFT_STRESS) ;
	vw.write();
// 	VoxelWriter vw0("sphere_strain", 50) ;
// 	vw0.getField(featureTree, VWFT_STRAIN) ;
// 	vw0.write();
	exit(0) ;
}


std::pair<double, double> centile(const Vector & v)
{
	Vector vs(v) ;
	std::sort(&vs[0], &vs[vs.size()]) ;
	return std::make_pair(vs[(vs.size()-1)*0.01], vs[(vs.size()-1)*0.99]) ;
}


int main(int argc, char *argv[])
{

	double nu = 0.2 ;
	double E = 1 ;
   Sample3D samplers(nullptr, 400,400,400,200,200,200) ;

	FeatureTree F(&samplers) ;
	featureTree = &F ;

	Matrix m0(6,6) ;
	m0[0][0] = 1. - nu ; m0[0][1] = nu ; m0[0][2] = nu ;
	m0[1][0] = nu ; m0[1][1] = 1. - nu ; m0[1][2] = nu ;
	m0[2][0] = nu ; m0[2][1] = nu ; m0[2][2] = 1. - nu ;
	m0[3][3] = 0.5 - nu ;
	m0[4][4] = 0.5 - nu ;
	m0[5][5] = 0.5 - nu ;
	m0 *= E/((1.+nu)*(1.-2.*nu)) ;

	nu = 0.2 ;
	E = 10 ;
	Matrix m1(6,6) ;
	m1[0][0] = 1. - nu ; m1[0][1] = nu ;      m1[0][2] = nu ;
	m1[1][0] = nu ;      m1[1][1] = 1. - nu ; m1[1][2] = nu ;
	m1[2][0] = nu ;      m1[2][1] = nu ;      m1[2][2] = 1. - nu ;
	m1[3][3] = 0.5 - nu ;
	m1[4][4] = 0.5 - nu ;
	m1[5][5] = 0.5 - nu ;
	m1 *= E/((1.+nu)*(1.-2.*nu)) ;



	
	MohrCoulomb * mc = new MohrCoulomb(30, -60) ;
	StiffnessAndFracture * sf = new StiffnessAndFracture(m0*0.5, mc) ;
	Stiffness * s = new Stiffness(m0) ;
	Stiffness * ss = new Stiffness(m1) ;

	samplers.setBehaviour(new Stiffness(m0)) ;
	Vector a(0.,6) ;// a[0] = 1 ; a[1] = 1 ; a[2] = 1 ; 
// 	ExpansiveZone3D inc(&samplers,100, 200, 200, 200, m1*4, a) ;
	Inclusion3D inc(100, 200, 200, 200) ;
	
// 	inc->setBehaviour(new StiffnessWithImposedDeformation(m1*4.,a)) ;
// 	inc.setBehaviour(new Stiffness(m1*4)) ;
	inc.setBehaviour(new VoidForm()) ;
	
	
    F.addFeature(&samplers, &inc) ;
// 	F.addFeature(&samplers, inc0) ;
    F.setSamplingNumber(atof(argv[1])) ;
// 		F.setProjectionOnBoundaries(false) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM_RIGHT_BACK)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, BOTTOM_RIGHT_BACK)) ;
// 	
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, TOP_LEFT_BACK)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, TOP_LEFT_BACK)) ;
// 	
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, BOTTOM_LEFT_FRONT)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM_LEFT_FRONT)) ;
	
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ZETA, FRONT, 1.)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_XI, RIGHT, 1.)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP, 1.)) ;
	F.addBoundaryCondition(new GeometryDefinedBoundaryCondition(SET_NORMAL_STRESS, inc.getPrimitive(), -1)) ;

// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, BOTTOM_LEFT_BACK)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM_LEFT_BACK)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, BOTTOM_LEFT_BACK)) ;
	
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, BOTTOM_RIGHT_BACK)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, TOP_LEFT_BACK)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, BOTTOM_LEFT_FRONT)) ;
	
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, BOTTOM_LEFT_BACK)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM_LEFT_BACK)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, BOTTOM_LEFT_BACK)) ;
	
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, BACK)) ;
	F.setOrder(QUADRATIC) ;

	step() ;
// 	delete dt ;
	
	return 0 ;
}
