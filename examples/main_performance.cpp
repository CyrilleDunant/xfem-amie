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
#include "./../utilities/parser/command_line_parser.h"
#include "./../utilities/writer/triangle_writer.h"
#include "./../utilities/writer/exodus_writer.h"


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
    parser.addFlag("--damage", "the test accounts for damage in the cement paste and aggregates", "-d") ;
    parser.addFlag("--viscous", "the test accounts for visco-elasticity of the cement paste", "-V") ;
    parser.addFlag("--pores", "use empty pores instead of aggregates", "-p") ;
    parser.addFlag("--no-export", "disable export of the triangles" ) ;
    parser.addFlag("--mesh-only", "stop simulation after meshing", "-m" ) ;
    parser.addValue("--inclusions", 10000, "number of aggregates (default 10,000)", "-i" ) ;
    parser.addValue("--seed", 1, "seed for the random placement of aggregates" ) ;
    parser.addString("--matrix-behaviour", "", "path to a *.ini file containing the matrix mechanical behaviour") ;
    parser.addString("--inclusion-behaviour", "", "path to a *.ini file containing the inclusion mechanical behaviour") ;
    parser.addString("--export-file", std::string(), "name of the file to export the triangles (auto-generated if not specified)", "-e" ) ;
    parser.parseCommandLine(argc, argv) ;

    bool damage = parser.getFlag("--damage") ;
    bool viscous = parser.getFlag("--viscous") ;
    bool pores = parser.getFlag("--pores") ;
    bool exp = !parser.getFlag("--no-export") ;
    bool bc = !parser.getFlag("--mesh-only") ;
    int inc = parser.getValue("--inclusions") ;
    size_t seed = parser.getValue("--seed") ;
    std::string matrixBehaviour = parser.getString("--matrix-behaviour") ;
    std::string inclusionBehaviour = parser.getString("--inclusion-behaviour") ;
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
    

    Form * paste = new PasteBehaviour(true) ;
    Form * agg = new AggregateBehaviour(true) ;

    if(viscous)
    {
        if(damage)
        {
            paste = new PasteBehaviour(false, true) ;
            agg = new AggregateBehaviour(false, true) ;
        }
        else
        {
            paste = new PasteBehaviour(true, true) ;
            agg = new AggregateBehaviour(true, true) ;
        }
    }
    else if(damage)
    {
        paste = new PasteBehaviour() ;
        agg = new AggregateBehaviour() ;
    }

    if(pores)
        agg = new VoidForm() ;

    paste = parser.getBehaviour("--matrix-behaviour", paste, SPACE_TWO_DIMENSIONAL) ;
    agg = parser.getBehaviour("--inclusion-behaviour", agg, SPACE_TWO_DIMENSIONAL) ;

    timeval time0, time1, time2 ;
    gettimeofday ( &time0, nullptr );

    RectangularFeature box(0.1,0.1,0,0) ;
    box.setBehaviour( paste ) ;

    FeatureTree F(&box) ;
    F.setSamplingRestriction( 0.002 ) ;
    F.setSamplingNumber(64) ;
    F.setDeltaTime( 0.01 ) ;
    parser.setFeatureTree( &F ) ;

    std::vector<Feature *> incs ;
    if(inc > 0)   
        incs = (PSDGenerator::get2DConcrete(&F, agg, inc, 0.02, 0.00005, new PSDBolomeA(), nullptr , 1e8 , 0.8, nullptr, std::vector<Geometry *>(), seed )) ;

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

//    MeshWriter::exportExodus2D("test.e.txt", &F, all) ;
    Vector strain = F.getAverageField( TOTAL_STRAIN_FIELD ) ;
    Vector stress = F.getAverageField( REAL_STRESS_FIELD ) ;
    Vector dmg = F.getAverageField( SCALAR_DAMAGE_FIELD ) ;

    std::cout << std::endl ;
    std::cout << std::endl ;
    F.printReport() ;

//    std::cout << F.getCurrentTime() << "\t" << strain[0] << "\t" << strain[1] << "\t" << strain[2] << "\t" << stress[0] << "\t" << stress[1] << "\t" << stress[2] << "\t" << dmg[0] << std::endl ;

    gettimeofday ( &time2, nullptr );

    if(exp)
    {
        TriangleWriter trg( file.c_str() , &F, 1.) ;
        trg.getField( TOTAL_STRAIN_FIELD ) ;
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
