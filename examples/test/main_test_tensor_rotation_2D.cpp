// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../main.h"
#include "../../features/features.h"
#include "../../physics/stiffness.h"
#include "../../physics/viscoelasticity_and_fracture.h"
#include "../../physics/logarithmic_creep_with_external_parameters.h"
#include "../../utilities/writer/triangle_writer.h"
//#include "../../physics/damagemodels/spacetimefiberbasedtensileshearlineardamage.h"
#include "../../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h"
#include "../../physics/damagemodels/spacetimebadisotropiclineardamage.h"
#include "../../physics/damagemodels/spacetimefiberbasedfixedcrack.h"
//#include "../../physics/fracturecriteria/spacetimemultisurfacefracturecriterion.h"
#include "../../physics/fracturecriteria/spacetimelimitsurfacefracturecriterion.h"
#include "../../physics/fracturecriteria/maxstrain.h"
#include "../../physics/homogenization/composite.h"
#include "../../features/sample.h"
#include "../../features/polygonSample.h"
#include "../../features/inclusion.h"
#include "../../utilities/parser/command_line_parser.h"
#include "../../utilities/postprocessor.h"
#include "../../utilities/itoa.h"

#include <fstream>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include <cmath>
#include <typeinfo>
#include <limits>
#include <sys/time.h>
#define DEBUG


using namespace Amie ;

void print(Vector & v)
{
	for(size_t i = 0 ; i < v.size() ; i++)
		std::cout << v[i] << "\t" ;
	std::cout << std::endl ;
}

void printVector(std::vector<int> v)
{
	for(size_t i = 0 ; i < v.size() ; i++)
		std::cout << v[i] << "\t" ;
	std::cout << std::endl ;
}

double invariant_1( Vector & v ) 
{
	return v[0]+v[1] ; 
}

double invariant_2( Vector & v ) 
{
	return v[0]*v[1] - v[2]*v[2] ; 
}

Vector invariant( Vector & v ) 
{
	Vector ret(2) ;
	ret[0] = invariant_1(v) ;
	ret[1] = invariant_2(v) ;
	return ret ;
}

int main(int argc, char *argv[])
{
	CommandLineParser parser("Test the 2D rotation of 4th-order tensors") ;
	parser.addFlag("--renew-base", "renew the base of results") ;
	parser.addString("--output-directory","../examples/test/","directory where the results are stored", "-D") ;
	parser.parseCommandLine(argc, argv) ;
	bool renew = parser.getFlag("--renew-base") ;
	std::string outdir = parser.getString("--output-directory") ;

	std::ofstream out ;
	if(renew)
		out.open(outdir+"/test_tensor_rotation_2D_base", std::ios::out) ;
	else
		out.open(outdir+"/test_tensor_rotation_2D_current", std::ios::out) ;

	Vector eigen(3) ; eigen[0] = 2e6 ; eigen[1] = -1e6 ; eigen[2] = 1e6 ;
	Vector reference = invariant(eigen) ;
	Vector eigenRotated(3) ;
	Vector invar(2) ;
	for(double theta = 0 ; theta < 2.01 ; theta += 0.1)
	{
		eigenRotated = Tensor::rotate2ndOrderTensor2D( eigen, theta*M_PI ) ;

		invar = invariant( eigenRotated) - reference ;

		out << theta ;
		for(size_t i = 0 ; i < 2 ; i++)
			out << "\t" << invar[i]/reference[i] ;
		out << std::endl ;


	}

	Matrix stiffness2D = Tensor::cauchyGreen( 10e9, 0.25, SPACE_TWO_DIMENSIONAL, PLANE_STRESS ) ;
	Vector shear(6) ; shear[0] = 2e-4 ; shear[1] = -1e-4 ; shear[2] = 1e-4 ;
	double lastTheta = 0. ;
	double dtheta = 0. ;

	Matrix incrementalStiffnessRotated(3,3) ;
	Matrix stiffnessRotated(3,3) ;

	Vector shearRotated(3) ;
	Vector stressRotated(3) ;
	Vector stressFromRotated(3) ;
	Vector stressUnRotated(3) ;

	Vector difference(3) ;

	for(double theta = 0 ; theta < 2.01 ; theta += 0.1 )
	{
		dtheta = theta-lastTheta ;		

		stiffnessRotated = Tensor::rotate4thOrderTensor2D( stiffness2D, theta*M_PI, 1 ) ;
		if( theta == 0 )
			incrementalStiffnessRotated = Tensor::rotate4thOrderTensor2D( stiffness2D, dtheta*M_PI, 1 ) ;
		else
			incrementalStiffnessRotated = Tensor::rotate4thOrderTensor2D( incrementalStiffnessRotated, dtheta*M_PI, 1 ) ;
		lastTheta = theta ;

//		stiffnessRotated.print() ;

		shearRotated = Tensor::rotate2ndOrderTensor2D( shear, -theta*M_PI ) ;

		stressRotated = stiffnessRotated * shear ;
		stressUnRotated = Tensor::rotate2ndOrderTensor2D( stressRotated, -theta*M_PI ) ;

		stressFromRotated = stiffness2D * shearRotated ;

		difference = stressUnRotated -  stressFromRotated ;

		out << theta << "\t" << abs((stiffnessRotated-incrementalStiffnessRotated).array()).max()/1e9 ;
		for(size_t i = 0 ; i < 3 ; i++)
			out << "\t" << difference[i]/1e6 ;
		out << std::endl ;
	}




	Vector components(21) ;
	components[0] = 12e9 ;
	components[1] = 10e9 ;
	components[2] = 13e9 ;
	components[3] = 6e9 ;
	components[4] = 4e9 ;
	components[5] = 5e9 ;
	components[6] = 0.5e9 ;
	components[7] = 0.3e9 ;
	components[8] = 1e9 ;
	components[9] = 2e9 ;
	components[10] = -1e9 ;
	components[11] = 3e9 ;
	components[12] = 0.2e9 ;
	components[13] = 0.5e9 ;
	components[14] = 1e9 ;
	components[15] = -3e9 ;
	components[16] = 0.4e9 ;
	components[17] = 4e9 ;
	components[18] = 3e9 ;
	components[19] = -2e9 ;
	components[20] = 0.5e9 ;

	Matrix stiffness3D = Tensor::orthotropicCauchyGreen( components, SYMMETRY_TRICLINIC, true ) ;
	stiffness2D = Tensor::to2D( stiffness3D, PLANE_STRAIN, XI ) ;
	lastTheta = 0. ;

	for(double theta = 0 ; theta < 2.01 ; theta += 0.1 )
	{
		dtheta = theta-lastTheta ;		

		stiffnessRotated = Tensor::rotate4thOrderTensor2D( stiffness2D, theta*M_PI, 1 ) ;
		if( theta == 0 )
			incrementalStiffnessRotated = Tensor::rotate4thOrderTensor2D( stiffness2D, dtheta*M_PI, 1 ) ;
		else
			incrementalStiffnessRotated = Tensor::rotate4thOrderTensor2D( incrementalStiffnessRotated, dtheta*M_PI, 1 ) ;
		lastTheta = theta ;

		shearRotated = Tensor::rotate2ndOrderTensor2D( shear, -theta*M_PI ) ;

		Matrix realStiffnessRotated = stiffnessRotated ;
		Matrix realStiffness = stiffness2D ;
		realStiffnessRotated[2][0] *= 0.5 ;
		realStiffnessRotated[2][1] *= 0.5 ;
		realStiffness[2][0] *= 0.5 ;
		realStiffness[2][1] *= 0.5 ;

		stressRotated = realStiffnessRotated * shear ;
		stressUnRotated = Tensor::rotate2ndOrderTensor2D( stressRotated, -theta*M_PI ) ;

		stressFromRotated = realStiffness * shearRotated ;

		difference = stressUnRotated -  stressFromRotated ;

		out << theta << "\t" << abs((stiffnessRotated-incrementalStiffnessRotated).array()).max()/1e9 ;
		for(size_t i = 0 ; i < 3 ; i++)
			out << "\t" << difference[i]/1e6 ;
		out << std::endl ;
	}


	return 0 ;
}

