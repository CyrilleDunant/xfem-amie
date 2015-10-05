// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "./main.h"
#include "./../features/features.h"
#include "./../features/sample.h"
#include "./../physics/materials/paste_behaviour.h"
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
    CommandLineParser parser("Make an elastic tensile test on a finely discretized concrete microstructure") ;
    parser.addFlag("--damage", false, "the test accounts for damage in the cement paste and aggregates") ;
    parser.addFlag("--viscous", false, "the test accounts for visco-elasticity of the cement paste") ;
    parser.addFlag("--pores", false, "use empty pores instead of aggregates") ;
    parser.addFlag("--no-export", false, "disable export of the triangles" ) ;
    parser.addFlag("--mesh-only", false, "stop simulation after meshing" ) ;
    parser.addValue("--inclusions", 10000, "number of aggregates (default 10,000)" ) ;
    parser.addString("--export-file", std::string(), "name of the file to export the triangles (auto-generated if not specified)" ) ;
    parser.parseCommandLine(argc, argv) ;

    bool damage = parser.getFlag("--damage") ;
    bool viscous = parser.getFlag("--viscous") ;
    bool pores = parser.getFlag("--pores") ;
    bool exp = !parser.getFlag("--no-export") ;
    bool bc = !parser.getFlag("--mesh-only") ;
    int inc = parser.getValue("--inclusions") ;
    std::string file = parser.getString("--export-file") ;
    if(file.size() == 0 && exp)
    {
        file = "performance" ;
        if(viscous)
            file += "_visco" ;
        if(damage)
            file += "_damage" ;
        if(pores)
            file += "_pores" ;
        file += "_trg" ;
    }
    

    Form * paste = new ElasticOnlyPasteBehaviour() ;
    Form * agg = new ElasticOnlyAggregateBehaviour() ;

    if(viscous)
    {
        if(damage)
        {
            paste = new ViscoDamagePasteBehaviour() ;
            agg = new ViscoDamageAggregateBehaviour() ;
        }
        else
        {
            paste = new ViscoElasticOnlyPasteBehaviour() ;
            agg = new ViscoElasticOnlyAggregateBehaviour() ;
        }
    }
    else if(damage)
    {
        paste = new PasteBehaviour() ;
        agg = new AggregateBehaviour() ;
    }

    if(pores)
        agg = new VoidForm() ;


    timeval time0, time1, time2 ;
    gettimeofday ( &time0, nullptr );

    Sample box(0.1,0.1,0,0) ;
    box.setBehaviour( paste ) ;

    FeatureTree F(&box) ;
    F.setSamplingRestriction( 8 ) ;
    F.setDeltaTime( 0.01 ) ;
    parser.setFeatureTree( &F ) ;

    if(inc > 0)   
        PSDGenerator::get2DConcrete(&F, agg, inc, 0.012, 0.00005, new PSDBolomeA(), nullptr , 1e8  ) ;

    if(!viscous)
    {
        F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, BOTTOM_LEFT ) ) ;
        F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM ) ) ;
        if(bc)
            F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_ALONG_ETA, TOP, 5e-5 ) ) ;
    }
    else
    {
        F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, BOTTOM_LEFT_AFTER ) ) ;
        F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM_AFTER ) ) ;
        if(bc)
            F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_ALONG_ETA, TOP_AFTER, 5e-5 ) ) ;
    }

    F.step() ;
    Vector strain = F.getAverageField( STRAIN_FIELD, -1, 1. ) ;
    Vector stress = F.getAverageField( REAL_STRESS_FIELD, -1, 1. ) ;
    Vector dmg = F.getAverageField( SCALAR_DAMAGE_FIELD , -1, 1. ) ;

    std::cout << std::endl ;
    std::cout << std::endl ;
    F.printReport() ;

//    std::cout << F.getCurrentTime() << "\t" << strain[0] << "\t" << strain[1] << "\t" << strain[2] << "\t" << stress[0] << "\t" << stress[1] << "\t" << stress[2] << "\t" << dmg[0] << std::endl ;

    gettimeofday ( &time2, nullptr );

    if(exp)
    {
        TriangleWriter trg( file.c_str() , &F, 1.) ;
        trg.getField( STRAIN_FIELD ) ;
        trg.getField( REAL_STRESS_FIELD ) ;
        trg.getField( SCALAR_DAMAGE_FIELD ) ;
        trg.getField( TWFT_STIFFNESS ) ;
        trg.getField( TWFT_CRITERION ) ;
        trg.write() ;
    }

    gettimeofday ( &time1, nullptr );
    double dt1 = time1.tv_sec * 1000000 - time0.tv_sec * 1000000 + time1.tv_usec - time0.tv_usec ;
    double dt2 = time1.tv_sec * 1000000 - time2.tv_sec * 1000000 + time1.tv_usec - time2.tv_usec ;
    std::cout << std::endl ;
    std::cout << "run time: " << dt1/1000000 << " seconds, including " << dt2/1000000 << " seconds for export" <<  std::endl ;


    return 0 ;
}
