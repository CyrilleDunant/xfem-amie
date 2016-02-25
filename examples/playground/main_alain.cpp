// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../main.h"
#include "../../features/features.h"
#include "../../features/sample.h"
#include "../../features/polygonSample.h"
#include "../../physics/stiffness.h"
#include "../../physics/logarithmic_creep_with_external_parameters.h"
#include "../../physics/material_laws/mechanical_material_laws.h"
#include "../../physics/material_laws/material_laws.h"
#include "../../utilities/parser.h"
#include "../../utilities/itoa.h"
#include "../../utilities/writer/triangle_writer.h"
#include "../../geometry/sampler/gradient_sampler.h" 
#include "../../geometry/sampler/regular_sampler.h" 
#include "../../utilities/mineral.h" 


#include <dirent.h>
#include <fstream>
#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h>
#include <sys/time.h>


using namespace Amie ;

double bulk( double E, double nu ) { return (1./3.)*E/(1.-2.*nu) ; }
double shear(double E, double nu ) { return (1./2.)*E/(1.+nu) ; }
double young(double k, double mu ) { return (9.*k*mu)/(3.*k+mu) ; }

double bulk_mt(double fa, double ka, double km, double mum)
{
    return km + (fa*(ka-km))/(1.+(1.-fa)*((ka-km)/km+(4./3.)*mum)) ;
}

double shear_mt(double fa, double mua, double mum, double km)
{
    return mum+(fa*(mua-mum))/(1.+(1.-fa)*( (mua-mum)/mum+(mum*(9.*km+8.*mum))/(6.*(km+2.*mum)))) ;
}

double beta_mt(double fa, double ka, double km, double mum)
{
    return (fa*ka*(3.*km+4.*mum))/(km*(3.*ka+4.*mum)-4.*fa*mum*(km-ka)) ;
}

double beta_simple(double fa, double ka, double km)
{
    return (2.*fa*ka)/(km+ka+fa*(ka-km)) ;
}

int main( int argc, char *argv[] )
{
    CommandLineParser parser("Run an analytical simulation with parameters found in a *.ini file", true, true) ;
    parser.addString( "--directory", std::string() , "directory to read additional input and write output", "-D") ;
    parser.addArgument( "file_name", "*.ini", "relative path to *.ini file to run") ;
    parser.parseCommandLine( argc, argv ) ;
    std::string file = parser.getStringArgument(0) ;
    std::string path = parser.getString("--directory") ;

    ConfigTreeItem * define = parser.getLocalConfiguration() ;
    std::map<std::string, std::string> direct = parser.getDirectConfiguration() ;
    std::vector<std::string> flags = parser.getActiveFlags() ;

    ConfigTreeItem * problem = ConfigParser::readFile(file, define, true, true, flags, path ) ;

    double E_agg = problem->getData("aggregate.young_modulus", 70e9) ;
    double nu_agg = 0.2 ;
    double epsilon_max = problem->getData("aggregate.maximum_radiation_expansion", 0.06) ;
    double phi_L = problem->getData("aggregate.fluence_latency", 15 ) ;
    double phi_C = problem->getData("aggregate.fluence_characteristic", 3 ) ;
    double alpha_agg = 8e-6 ;
    double f_agg = 0.5 ;

    double E_paste = 20e9 ;
    double nu_paste = 0.2 ;
    double strength_paste = problem->getData("paste.tensile_strength", 5e6 ) ;
    double fracture_energy = 0 ;
    double eta_paste = problem->getData("paste.creep_viscosity", 40e9 ) ;
    double tau_paste = 2 ;
    double humidity_coefficient = 0.2 ;
    double alpha_paste = 10e-6 ;
    double drying_paste = 1000e-6 ;

    double T = problem->getData("conditions.temperature", 318 ) ;
    double h = problem->getData("conditions.relative_humidity", 0.95 ) ;
    double phi_80 = problem->getData("conditions.final_fluence", 10 ) ;

    double k_agg = bulk( E_agg, nu_agg ) ;
    double mu_agg = shear( E_agg, nu_agg ) ;
    double k_paste = bulk( E_paste, nu_paste ) ;
    double mu_paste = shear( E_paste, nu_paste ) ;

    double k_c = bulk_mt( f_agg, k_agg, k_paste, mu_paste ) ;
    double mu_c = shear_mt( f_agg, mu_agg, mu_paste, k_paste) ;

    double xhi_k = (16./9.)*(1.-nu_paste*nu_paste)/(1.-2.*nu_paste) ;
    double xhi_m = (32./45.)*(1.-nu_paste)*(5.-nu_paste)/(2.-nu_paste) ;

    double E_c = young(k_c, mu_c ) ;

    double damage = 0 ;

    for(double i = 0 ; i < 81 ; i+=1)
    {
        double phi = phi_80*((double) i / 80.) ;

        double epsilon_paste = (T-293)*alpha_paste + (h-1)*drying_paste ;
        double epsilon_agg = (T-293)*alpha_agg + epsilon_max*(1.-exp(-phi/phi_C))/(1.+exp(-(phi-phi_L)/phi_C)) ;

        if( E_agg * (epsilon_agg - epsilon_paste) > strength_paste )
        {
            damage = std::max(damage, std::pow(((E_agg * epsilon_agg / strength_paste-1.))/((eta_paste/E_paste)*i/tau_paste), 5. )) ;
            damage = std::min( damage, 9./16. ) ;
        }

        double k_paste_damaged = std::max(0., k_paste / (1.+xhi_k*damage)) ;
        double mu_paste_damaged = std::max(0., mu_paste / (1.+xhi_m*damage)) ;

        double beta = beta_simple( f_agg, k_agg, k_paste_damaged ) ;

        double epsilon_concrete = beta*epsilon_agg + (1.-beta)*epsilon_paste ;
        double k_concrete = bulk_mt( f_agg, k_agg, k_paste_damaged, mu_paste_damaged ) ;
        double mu_concrete = shear_mt( f_agg, mu_agg, mu_paste_damaged, k_paste_damaged) ;

        double E_c_damaged = young( k_concrete, mu_concrete ) ;

        std::cout << phi << "\t" << epsilon_concrete << "\t" << 1.-E_c_damaged/E_c << std::endl ;

    }


    return 0 ;
}
