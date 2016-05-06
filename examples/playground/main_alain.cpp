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
#include "../../physics/fracturecriteria/spacetimemultisurfacefracturecriterion.h"
#include "../../physics/fracturecriteria/spacetimelimitsurfacefracturecriterion.h"
#include "../../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h"
#include "../../physics/damagemodels/spacetimefiberbasedfixedcrack.h"
#include "../../physics/damagemodels/spacetimebifurcation.h"
#include "../../physics/material_laws/mechanical_material_laws.h"
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
        CommandLineParser parser("Test a fixed crack damage behaviour on two elements") ;
        parser.addFlag("--renew-base", "renew the base of results") ;
        parser.addString("--output-directory","../examples/test/","directory where the results are stored", "-D") ;
        parser.parseCommandLine(argc, argv) ;
        bool renew = parser.getFlag("--renew-base") ;
        std::string outdir = parser.getString("--output-directory") ;

        std::ofstream out ;
        if(renew)
                out.open(outdir+"/test_spacetime_fixed_crack_base", std::ios::out) ;
        else
                out.open(outdir+"/test_spacetime_fixed_crack_current", std::ios::out) ;

        Sample rect(nullptr, 0.01,0.01,0,0) ;

        SpaceTimeMultiSurfaceFractureCriterion * crit = new SpaceTimeMultiSurfaceFractureCriterion() ;
        crit->add(new SpaceTimeLimitFirstStrainInvariant(-1.8e-3)) ;
        crit->add(new SpaceTimeNonLocalMaximumStress(15e6)) ;

        LogarithmicCreepWithExternalParameters * toto = new LogarithmicCreepWithExternalParameters("young_modulus = 10e9, poisson_ratio = 0.2, imposed_deformation = 0, microcracking = -1e-3", crit, new SpaceTimeBifurcationAndDamage(0.5, new SpaceTimeFiberBasedIsotropicLinearDamage(0.001,0.00001,0.99))) ;

        rect.setBehaviour( toto ) ;

        FeatureTree f(&rect) ;
        f.setSamplingNumber(0) ;
        f.setDeltaTime(0.001) ;
        f.setMinDeltaTime(1e-9) ;

        f.step() ;
        f.step() ;

        f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, LEFT_AFTER) ) ;
        f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM_AFTER ) ) ;
        BoundingBoxDefinedBoundaryCondition * top = new BoundingBoxDefinedBoundaryCondition( SET_ALONG_ETA, TOP_AFTER, 0) ;
        f.addBoundaryCondition( top ) ;

        Vector stress(3) ;
        Vector strain(3) ;
        Vector alpha(3) ;

        for(double i = 0. ; i < 10 ; i++)
        {

            top->setData( - 0.00001*i )   ;

                f.step() ;

                stress = f.getAverageField( REAL_STRESS_FIELD, 1. ) ;
                strain = f.getAverageField( MECHANICAL_STRAIN_FIELD, 1. ) ;
                alpha = f.getAverageField( IMPOSED_STRAIN_FIELD, 1. ) ;
                for(size_t i = 0 ; i < stress.size() ; i++)
                {
                    if(std::abs(stress[i]) < 1e3) { stress[i] = 0 ; }
                    if(std::abs(strain[i]) < 1e-6) { strain[i] = 0 ; }
                }
                std::cout << top->getData()*1e6 << "\t" << stress[0]/1e6  << "\t" << stress[1]/1e6 << "\t" << stress[2]/1e6 << "\t" << strain[0]*1e3  << "\t" << strain[1]*1e3 << "\t" << strain[2]*1e3 << "\t" << alpha[0]*1e3 << "\t" << alpha[1]*1e3 << "\t" << alpha[2]*1e3 << std::endl ;
        }

        for(double i = 10. ; i > -20.5 ; i--)
        {

            top->setData( - 0.00001*i )   ;

                f.step() ;

                stress = f.getAverageField( REAL_STRESS_FIELD, 1. ) ;
                strain = f.getAverageField( MECHANICAL_STRAIN_FIELD, 1. ) ;
                alpha = f.getAverageField( IMPOSED_STRAIN_FIELD, 1. ) ;
                for(size_t i = 0 ; i < stress.size() ; i++)
                {
                    if(std::abs(stress[i]) < 1e3) { stress[i] = 0 ; }
                    if(std::abs(strain[i]) < 1e-6) { strain[i] = 0 ; }
                }
                std::cout << top->getData()*1e6 << "\t" << stress[0]/1e6  << "\t" << stress[1]/1e6 << "\t" << stress[2]/1e6 << "\t" << strain[0]*1e3  << "\t" << strain[1]*1e3 << "\t" << strain[2]*1e3 << "\t" << alpha[0]*1e3 << "\t" << alpha[1]*1e3 << "\t" << alpha[2]*1e3 << std::endl ;
        }

        for(double i = -19 ; i < 0.5 ; i++)
        {

            top->setData( - 0.00001*i )   ;

                f.step() ;

                stress = f.getAverageField( REAL_STRESS_FIELD, 1. ) ;
                strain = f.getAverageField( MECHANICAL_STRAIN_FIELD, 1. ) ;
                alpha = f.getAverageField( IMPOSED_STRAIN_FIELD, 1. ) ;
                for(size_t i = 0 ; i < stress.size() ; i++)
                {
                    if(std::abs(stress[i]) < 1e3) { stress[i] = 0 ; }
                    if(std::abs(strain[i]) < 1e-6) { strain[i] = 0 ; }
                }
                std::cout << top->getData()*1e6 << "\t" << stress[0]/1e6  << "\t" << stress[1]/1e6 << "\t" << stress[2]/1e6 << "\t" << strain[0]*1e3  << "\t" << strain[1]*1e3 << "\t" << strain[2]*1e3 << "\t" << alpha[0]*1e3 << "\t" << alpha[1]*1e3 << "\t" << alpha[2]*1e3 << std::endl ;
        }

        return 0 ;
}

