// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../main.h"
#include "../../features/features.h"
#include "../../physics/stiffness.h"
#include "../../physics/logarithmic_creep_with_external_parameters.h"
#include "../../utilities/writer/triangle_writer.h"
#include "../../features/sample.h"
#include "../../features/polygonSample.h"
#include "../../features/inclusion.h"
#include "../../utilities/parser.h"
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

int main(int argc, char *argv[])
{
	double E = 10e9 ;
	double nu = 0.2 ;
	double K = E/(3*(1.-2.*nu)) ;
	double G = E/(2*(1.+nu)) ;

	std::cout << E << "\t" << nu <<  "\t" << K <<  "\t" << G << std::endl ;

	Matrix C_strain = Tensor::cauchyGreen( E, nu, SPACE_TWO_DIMENSIONAL, PLANE_STRAIN, YOUNG_POISSON) ;
	Matrix C_stress = Tensor::cauchyGreen( E, nu, SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON) ;

	C_strain.print() ;
	C_stress.print() ;

	std::pair<double,double> e_nu_strain = Tensor::getIsotropicMaterialParameters( C_strain, YOUNG_POISSON, PLANE_STRAIN ) ;
	std::pair<double,double> e_nu_stress = Tensor::getIsotropicMaterialParameters( C_stress, YOUNG_POISSON, PLANE_STRESS ) ;
	std::pair<double,double> k_g_strain = Tensor::getIsotropicMaterialParameters( C_strain, BULK_SHEAR, PLANE_STRAIN ) ;
	std::pair<double,double> k_g_stress = Tensor::getIsotropicMaterialParameters( C_stress, BULK_SHEAR, PLANE_STRESS ) ;
	std::pair<double,double> e_g_strain = Tensor::getIsotropicMaterialParameters( C_strain, YOUNG_SHEAR, PLANE_STRAIN ) ;
	std::pair<double,double> e_g_stress = Tensor::getIsotropicMaterialParameters( C_stress, YOUNG_SHEAR, PLANE_STRESS ) ;

	std::cout << E << "\t" << e_nu_strain.first << std::endl ;
	std::cout << nu << "\t" << e_nu_strain.second << std::endl ;
	std::cout << E << "\t" << e_nu_stress.first << std::endl ;
	std::cout << nu << "\t" << e_nu_stress.second << std::endl ;
	std::cout << K << "\t" << k_g_strain.first << std::endl ;
	std::cout << G << "\t" << k_g_strain.second << std::endl ;
	std::cout << K << "\t" << k_g_stress.first << std::endl ;
	std::cout << G << "\t" << k_g_stress.second << std::endl ;
	std::cout << E << "\t" << e_g_strain.first << std::endl ;
	std::cout << G << "\t" << e_g_strain.second << std::endl ;
	std::cout << E << "\t" << e_g_stress.first << std::endl ;
	std::cout << G << "\t" << e_g_stress.second << std::endl ;




	return 0 ;
}

