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
	return v[0]+v[1]+v[2] ; 
}

double invariant_2( Vector & v ) 
{
	return v[0]*v[1] + v[1]*v[2] + v[2]*v[0] - v[3]*v[3] - v[4]*v[4] - v[5]*v[5] ; 
}

double invariant_3( Vector & v ) 
{
	return v[0]*v[1]*v[2]+2.*v[3]*v[4]*v[5]-v[5]*v[5]*v[2]-v[3]*v[3]*v[0]-v[4]*v[4]*v[1] ;
}

Vector invariant( Vector & v ) 
{
	Vector ret(3) ;
	ret[0] = invariant_1(v) ;
	ret[1] = invariant_2(v) ;
	ret[2] = invariant_3(v) ;
	return ret ;
}

int main(int argc, char *argv[])
{
	CommandLineParser parser("Test the 3D rotation of 4th-order tensors") ;
	parser.addFlag("--renew-base", "renew the base of results") ;
	parser.addString("--output-directory","../examples/test/","directory where the results are stored", "-D") ;
	parser.parseCommandLine(argc, argv) ;
	bool renew = parser.getFlag("--renew-base") ;
	std::string outdir = parser.getString("--output-directory") ;

	std::ofstream out ;
	if(renew)
		out.open(outdir+"/test_tensor_rotation_3D_base", std::ios::out) ;
	else
		out.open(outdir+"/test_tensor_rotation_3D_current", std::ios::out) ;

	Vector eigen(6) ; eigen[0] = 1e6 ; eigen[1] = -1e6 ; eigen[2] = 1e6 ; eigen[3] = 0.5e6 ; eigen[4] = 0 ; eigen[5] = 0.5e6 ;
	Vector reference = invariant(eigen) ;
	Vector eigenRotated(6) ;
	Vector invar(3) ;
	for(double theta = 0 ; theta < 2.01 ; theta += 0.1)
	{
		Point p( theta*M_PI, 0, 0 ) ;
		eigenRotated = Tensor::rotate2ndOrderTensor3D( eigen, p, 1. ) ;

		invar = invariant( eigenRotated) - reference ;

		out << theta ;
		for(size_t i = 0 ; i < 3 ; i++)
			out << "\t" << invar[i]/reference[i] ;
		out << std::endl ;


	}



	Matrix stiffness3D = Tensor::cauchyGreen( 10e9, 0.25, SPACE_THREE_DIMENSIONAL ) ;
	Vector shear(6) ; shear[0] = 0 ; shear[1] = 0 ; shear[2] = 0 ; shear[3] = 1e-4 ; shear[4] = 1e-4 ; shear[5] = 1e-4 ;
	double lastTheta = 0. ;

	Matrix incrementalStiffnessRotated(6,6) ;
	Matrix stiffnessRotated(6,6) ;

	Vector shearRotated(6) ;
	Vector stressRotated(6) ;
	Vector stressFromRotated(6) ;
	Vector stressUnRotated(6) ;

	Vector difference(6) ;

	for(double theta = 0 ; theta < 2.01 ; theta += 0.1 )
	{
		Point p( theta*M_PI, 0, 0 ) ;
		Point q( (theta-lastTheta)*M_PI, 0, 0 ) ;
		lastTheta = theta ;

		stiffnessRotated = Tensor::rotate4thOrderTensor3D( stiffness3D, theta*M_PI, 0, 0, 1 ) ;
		if( theta == 0 )
			incrementalStiffnessRotated = Tensor::rotate4thOrderTensor3D( stiffness3D, (theta-lastTheta)*M_PI, 0, 0, 1 ) ;
		else
			incrementalStiffnessRotated = Tensor::rotate4thOrderTensor3D( incrementalStiffnessRotated, (theta-lastTheta)*M_PI, 0, 0, 1 ) ;

		shearRotated = Tensor::rotate2ndOrderTensor3D( shear, -p, -1. ) ;

		stressRotated = stiffnessRotated * shear ;
		stressUnRotated = Tensor::rotate2ndOrderTensor3D( stressRotated, -p, 1.) ;

		stressFromRotated = stiffness3D * shearRotated ;

		difference = stressUnRotated -  stressFromRotated ;

		out << theta << "\t" << abs((stiffnessRotated-incrementalStiffnessRotated).array()).max()/1e9 ;
		for(size_t i = 0 ; i < 6 ; i++)
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

	stiffness3D = Tensor::orthotropicCauchyGreen( components, SYMMETRY_TRICLINIC, true ) ;
	lastTheta = 0. ;

	for(double theta = 0 ; theta < 2.01 ; theta += 0.1 )
	{
		Point p( theta*M_PI, 0, 0 ) ;
		Point q( (theta-lastTheta)*M_PI, 0, 0 ) ;
		lastTheta = theta ;

		stiffnessRotated = Tensor::rotate4thOrderTensor3D( stiffness3D, p, -1 ) ;
		if( theta == 0 )
			incrementalStiffnessRotated = Tensor::rotate4thOrderTensor3D( stiffness3D, q, -1 ) ;
		else
			incrementalStiffnessRotated = Tensor::rotate4thOrderTensor3D( incrementalStiffnessRotated, q, -1 ) ;

		Matrix realStiffnessRotated = stiffnessRotated ;
		Matrix realStiffness = stiffness3D ;
		for(size_t i = 0 ; i < 3 ; i++)
		{
			for(size_t j = 0 ; j < 3 ; j++)
			{
				realStiffnessRotated[i+3][j] *= 0.5 ;
				realStiffness[i+3][j] *= 0.5 ;
			}

		}

		shearRotated = Tensor::rotate2ndOrderTensor3D( shear, -p, -1. ) ;

		stressRotated = realStiffnessRotated * shear ;
		stressUnRotated = Tensor::rotate2ndOrderTensor3D( stressRotated, -p , -1.) ;

		stressFromRotated = realStiffness * shearRotated ;

		difference = stressUnRotated -  stressFromRotated ;
		out << theta << "\t" << abs((stiffnessRotated-incrementalStiffnessRotated).array()).max()/1e9 ;
		for(size_t i = 0 ; i < 6 ; i++)
			out << "\t" << difference[i]/1e6 ;
		out << std::endl ;
	}


	return 0 ;
}

