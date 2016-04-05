// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../main.h"
#include "../../features/features.h"
#include "../../physics/viscoelasticity.h"
#include "../../physics/viscoelasticity_and_fracture.h"
#include "../../physics/fracturecriteria/spacetimemultilinearsofteningfracturecriterion.h"
#include "../../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h"
#include "../../physics/materials/paste_behaviour.h"
#include "../../utilities/writer/triangle_writer.h"
#include "../../features/sample.h"
#include "../../utilities/parser/command_line_parser.h"

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
	CommandLineParser parser("Test the construction of Cauchy-Green stiffness matrices") ;
	parser.addFlag("--renew-base", "renew the base of results") ;
	parser.addString("--output-directory","../examples/test/","directory where the results are stored", "-D") ;
	parser.parseCommandLine(argc, argv) ;
	bool renew = parser.getFlag("--renew-base") ;
	std::string outdir = parser.getString("--output-directory") ;

	double E = 10e9 ;
	double nu = 0.2 ;
	double K = E/(3*(1.-2.*nu)) ;
	double G = E/(2*(1.+nu)) ;

	std::vector< IsotropicMaterialParameters > types ;
	types.push_back( YOUNG_POISSON ) ;
	types.push_back( BULK_SHEAR ) ;
	types.push_back( YOUNG_SHEAR ) ;

	Vector firsts(3) ;
	Vector seconds(3) ;
	firsts[0] = E ;
	firsts[1] = K ;
	firsts[2] = E ;
	seconds[0] = nu ;
	seconds[1] = G ;
	seconds[2] = G ;

	std::vector< planeType > planes ;
	planes.push_back( PLANE_STRAIN) ;
	planes.push_back( PLANE_STRESS) ;

	std::ofstream out ;
	if(renew)
		out.open(outdir+"/test_cauchy_green_base", std::ios::out) ;
	else
		out.open(outdir+"/test_cauchy_green_current", std::ios::out) ;

	out << E << "\t" << nu << "\t" << K << "\t" << G << std::endl ;

	for(size_t i = 0 ; i < types.size() ; i++)
	{
		for(size_t j = 0 ; j < planes.size() ; j++)
		{
			Matrix C = Tensor::cauchyGreen( firsts[i], seconds[i], SPACE_TWO_DIMENSIONAL, planes[j], types[i]) ;
			out << C[0][0] << "\t" << C[0][1] << "\t" << C[2][2] ;
			for(size_t k = 0 ; k < types.size() ; k++)
			{
				std::pair<double,double> props = Tensor::getIsotropicMaterialParameters( C, types[k], planes[j] ) ;
				out << "\t" << props.first << "\t" << props.second ;
			}
			out << std::endl ;
		}

	}

	return 0 ;
}

