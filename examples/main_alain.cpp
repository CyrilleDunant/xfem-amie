// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "./main.h"
#include "./../features/features.h"
#include "./../features/sample.h"
#include "./../physics/logarithmic_creep_with_external_parameters.h"
#include "./../physics/materials/aggregate_behaviour.h"
#include "./../utilities/granulo.h"
#include "./../utilities/parser.h"
#include "./../utilities/writer/triangle_writer.h"


#include <fstream>
#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h>
#include <sys/time.h>


using namespace Amie ;

int main( int argc, char *argv[] )
{
    CommandLineParser parser("Make a tensile test on a concrete microstructure with a small notch") ;
    parser.addFlag("--mesh-only", false, "stop simulation after meshing" ) ;
    parser.addValue("--inclusions", 1000, "number of aggregates (default 1,000)" ) ;
    parser.addValue("--notch",0.01,"size of the notch (default 0.01)") ;
    parser.addValue("--seed",1,"start of the random number sequence (default 1)") ;
    parser.addValue("--stress",0,"applied stress") ;
    parser.addValue("--strain-rate",0.0001,"applied strain rate") ;
    parser.parseCommandLine(argc, argv) ;

    bool bc = !parser.getFlag("--mesh-only") ;
    int inc = parser.getValue("--inclusions") ;
    int seed = parser.getValue("--seed") ;
    double length = parser.getValue("--notch") ;
    double stress = parser.getValue("--stress") ;
    double rate = parser.getValue("--strain-rate") ;
    std::string file = "tension_notch_trg" ;

    Form * paste = new LogarithmicCreepWithExternalParameters("young_modulus = 12e9, poisson_ratio = 0.2, creep_modulus = 40e9, creep_characteristic_time = 1") ;
    Form * agg = new LogarithmicCreepWithExternalParameters("young_modulus = 60e9, poisson_ratio = 0.2") ;
    
    Sample box(0.1,0.1,0,0) ;
    box.setBehaviour( paste ) ;

    Sample notch(length, length*0.2, -0.05+length*0.5, 0) ;
    notch.setBehaviour(new VoidForm()) ;
    std::vector<Geometry *> geom ; geom.push_back( dynamic_cast<Rectangle *>(&notch) ) ;

    FeatureTree F(&box) ;
    F.setSamplingNumber(256) ;
    F.setOrder(LINEAR_TIME_LINEAR) ;
    F.setSamplingRestriction( 8 ) ;
    F.setDeltaTime( 0.01 ) ;
    F.setSolverPrecision( 1e-6) ;
    parser.setFeatureTree( &F ) ;

    F.addFeature( &box, &notch ) ;
    F.setSamplingFactor( &notch, 1.5 ) ;
    PSDGenerator::get2DConcrete(&F, agg, inc, 0.012, 0.00005, new PSDBolomeA(), nullptr , 1e8, 0.8, nullptr, geom, seed ) ;

    F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, BOTTOM_LEFT_AFTER ) ) ;
    F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM_AFTER ) ) ;
    BoundingBoxDefinedBoundaryCondition * strainBC = new BoundingBoxDefinedBoundaryCondition( SET_ALONG_ETA, TOP_AFTER, 0 ) ;
    if(bc)
    {
        if(stress > POINT_TOLERANCE)
             F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_STRESS_ETA, TOP_AFTER, stress ) ) ;
        else
             F.addBoundaryCondition(strainBC) ;
    }

    F.step() ;
    Vector epsilon = F.getAverageField( STRAIN_FIELD, -1, 1. ) ;
    Vector sigma = F.getAverageField( REAL_STRESS_FIELD, -1, 1. ) ;
    Vector dmg = F.getAverageField( SCALAR_DAMAGE_FIELD , -1, 1. ) ;

    std::cout << F.getCurrentTime() << "\t" << epsilon[0] << "\t" << epsilon[1] << "\t" << epsilon[2] << "\t" << sigma[0] << "\t" << sigma[1] << "\t" << sigma[2] << "\t" << dmg[0] << std::endl ;

    TriangleWriter trg( file.c_str() , &F, 1.) ;
    trg.getField( STRAIN_FIELD ) ;
    trg.getField( REAL_STRESS_FIELD ) ;
    trg.getField( SCALAR_DAMAGE_FIELD ) ;
    trg.getField( TWFT_STIFFNESS ) ;
    trg.getField( TWFT_CRITERION ) ;
    trg.write() ;

    if(bc)
    {
        for(size_t i = 0 ; i < 100 ; i++)
        {
            strainBC->setData( rate*0.01*i ) ;

            F.step() ;
            epsilon = F.getAverageField( STRAIN_FIELD, -1, 1. ) ;
            sigma = F.getAverageField( REAL_STRESS_FIELD, -1, 1. ) ;
            dmg = F.getAverageField( SCALAR_DAMAGE_FIELD , -1, 1. ) ;

            std::cout << F.getCurrentTime() << "\t" << epsilon[0] << "\t" << epsilon[1] << "\t" << epsilon[2] << "\t" << sigma[0] << "\t" << sigma[1] << "\t" << sigma[2] << "\t" << dmg[0] << std::endl ;
        }
    }


    return 0 ;
}
