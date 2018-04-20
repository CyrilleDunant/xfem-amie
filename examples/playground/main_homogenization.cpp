// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../main.h"
#include "../../features/features.h"
#include "../../physics/stiffness.h"
#include "../../physics/homogenization/phase.h"
#include "../../physics/homogenization/composite.h"
#include "../../physics/logarithmic_creep_with_external_parameters.h"
#include "../../utilities/writer/triangle_writer.h"
#include "../../features/sample.h"
#include "../../features/polygonSample.h"
#include "../../features/inclusion.h"
#include "../../utilities/parser/command_line_parser.h"
#include "../../utilities/parser/config_parser.h"
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

Vector readVector( ConfigTreeItem * problem, std::string arg, std::string def ) 
{
	std::string ret = problem->getStringData( arg, def ) ;
	if(ret.length() > 0)
		return  ConfigTreeItem::readLineAsVector( problem->getStringData( arg, def ) ) ;
	Vector v(1) ;
	v[0] = problem->getData( arg, atof(def.c_str()) ) ; 
	return v ;
}


Composite * getComposite( std::string scheme, Phase matrix, Phase inclusion )
{
	if( scheme == "Dilute" ) { return new DiluteMatrixInclusionComposite( matrix, inclusion ) ; }
	if( scheme == "Voigt" ) { return new VoigtMatrixInclusionComposite( matrix, inclusion ) ; }
	if( scheme == "Reuss" ) { return new ReussMatrixInclusionComposite( matrix, inclusion ) ; }
	if( scheme == "MoriTanaka" ) { return new MoriTanakaMatrixInclusionComposite( matrix, inclusion ) ; }
	if( scheme == "InverseMoriTanaka" ) { return new InverseMoriTanakaMatrixInclusionComposite( matrix, inclusion ) ; }
	if( scheme == "BiphasicSelfConsistent" ) { return new BiphasicSelfConsistentComposite( matrix, inclusion ) ; }
	return new MoriTanakaMatrixInclusionComposite(matrix, inclusion) ;
}

int main(int argc, char *argv[])
{
	CommandLineParser parser("Makes a series of homogenization analyses on a biphasic elastic material") ;
	parser.addArgument("input_file", "../examples/playground/homogenization.ini", "*.ini file containing the description of the analysis") ;
	parser.addString("--scheme","MoriTanaka","homogenization scheme to be used", "-s") ;
	parser.addString("--output","homogenization.txt","files where the results are stored", "-o") ;
	parser.disableFeatureTreeArguments() ;
	parser.parseCommandLine(argc, argv) ;
	std::string file = parser.getStringArgument("input_file") ;

	std::string input = parser.getStringArgument("input_file") ;
	std::map<std::string, std::string> direct = parser.getDirectConfiguration() ;
	ConfigTreeItem * problem = ConfigParser::readFile(input, nullptr, false ) ;
	problem->configure( direct ) ;

	Vector E_mat = readVector( problem, "matrix.young_modulus","1e9") ;
	Vector nu_mat = readVector( problem, "matrix.poisson_ratio","0.2") ;
	Vector a_mat = readVector( problem, "matrix.imposed_deformation","0") ;
//	Vector f_mat = ConfigTreeItem::readLineAsVector( problem->getStringData("matrix.volume_fraction","0.5") ) ;
		
	Vector E_inc = readVector( problem, "inclusions.young_modulus","1e9") ;
	Vector nu_inc = readVector( problem, "inclusions.poisson_ratio","0.2") ;
	Vector a_inc = readVector( problem, "inclusions.imposed_deformation","0") ;
	Vector f_inc = readVector( problem, "inclusions.volume_fraction","0.5") ;

	std::string scheme = problem->getStringData("scheme", parser.getString("--scheme")) ;
	std::string output = problem->getStringData("output",parser.getString("--output")) ;
	std::fstream out ;
	out.open( output.c_str(), std::ios::out ) ;
	SpaceDimensionality d3 = SPACE_THREE_DIMENSIONAL ;
	
	for(size_t imE = 0 ; imE < E_mat.size() ; imE++)
	{
		for(size_t imn = 0 ; imn < nu_mat.size() ; imn++)
		{
			for(size_t ima = 0 ; ima < a_mat.size() ; ima++)
			{
				for(size_t iiE = 0 ; iiE < E_inc.size() ; iiE++)
				{
					for(size_t iin = 0 ; iin < nu_inc.size() ; iin++)
					{
						for(size_t iia = 0 ; iia < a_inc.size() ; iia++)
						{
							for(size_t iif = 0 ; iif < f_inc.size() ; iif++)
							{
								                        StiffnessWithImposedStrain matrixBehaviour( E_mat[imE], nu_mat[imn], a_mat[ima], d3 ) ;
								                        StiffnessWithImposedStrain inclusionBehaviour( E_inc[iiE], nu_inc[iin], a_inc[iia], d3 ) ;
								Phase matrix( &matrixBehaviour, 1.-f_inc[iif] ) ;
								Phase inclusion( &inclusionBehaviour, f_inc[iif] ) ;

								Composite * composite = getComposite( scheme, matrix, inclusion ) ;
								composite->apply() ;
								Matrix compliance = composite->C ;
								Composite::invertTensor( compliance ) ;
								Vector alpha = compliance*(composite->beta) ;

								double E_hom = 1./compliance[0][0] ;
								double nu_hom = -compliance[0][1]*E_hom ;
								double a_hom = alpha[0] ;

								std::cout << f_inc[iif] << "\t" << E_mat[imE] << "\t" << nu_mat[imn] << "\t" << a_mat[ima] << "\t" << E_inc[iiE] << "\t" << nu_inc[iin] << "\t" << a_inc[iia] << "\t" << E_hom << "\t" << nu_hom << "\t" << a_hom << std::endl ;
								out << f_inc[iif] << "\t" << E_mat[imE] << "\t" << nu_mat[imn] << "\t" << a_mat[ima] << "\t" << E_inc[iiE] << "\t" << nu_inc[iin] << "\t" << a_inc[iia] << "\t" << E_hom << "\t" << nu_hom << "\t" << a_hom << std::endl ;

//								delete composite ;
							}
						}
					}
				}
			}
		}

	}




	return 0 ;
}

