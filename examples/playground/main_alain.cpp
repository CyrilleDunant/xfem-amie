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
#include "../../utilities/parser.h"
#include "../../utilities/writer/triangle_writer.h"
#include "../../geometry/level_set.h" 


#include <fstream>
#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h>
#include <sys/time.h>


using namespace Amie ;


int main( int argc, char *argv[] )
{
    Sample box(0.064,0.04,0,0) ;
    box.setBehaviour( new Stiffness( 10e9, 0.2 ) ) ;

    FeatureTree F(&box) ;

    PolygonGranuloFromFile reader("/home/ag3/Code/cv/test.poly") ;
    std::vector<Feature *> polys =  reader.getFeatures( SPACE_TWO_DIMENSIONAL, &box, 2 ) ;
    for(size_t i = 0 ; i < polys.size() ; i++)
    {
        polys[i]->setBehaviour( new Stiffness( 20e9 + i*10e9, 0.2 ) ) ;
        if(reader.getData( 1, i ) >= 0 && reader.getData( 1, i ) < polys.size() )
            polys[i]->setFather( polys[ reader.getData( 1, i ) ] ) ;
//	if(i == 11)
            F.addFeature( polys[i]->getFather(), polys[i] ) ;
    }

    F.setSamplingNumber( 16 ) ;
    F.setMinimumMeshDensity( 0.4 ) ;
    CommandLineParser::setFeatureTree( &F, argc, argv ) ;
    F.step() ;

/*    int cache = F.get2DMesh()->generateCache( dynamic_cast<Geometry *>( polys[7] ) ) ;
    std::vector<int> elements = F.get2DMesh()->getCache(cache) ;
    Vector y(elements.size()) ;
    for(size_t i = 0 ; i < elements.size() ; i++)
    {
	y[i] = F.get2DMesh()->getElement( cache, i )->getCenter().getY() ;
    }
    size_t j = 0 ;
    while( y[j] > y.min() ) { j++ ; } 
    std::cout << j << " " ; F.get2DMesh()->getElement( cache, j )->getCenter().print() ;

    std::valarray<Point> pts = dynamic_cast<PolygonalSample *>(polys[7])->getOriginalPoints() ;
    for(size_t i = 0 ; i < pts.size() ; i++)
        pts[i].print() ;*/


    TriangleWriter trg("test_poly", &F, 1) ;
    trg.getField( TWFT_STIFFNESS) ; 
    trg.write() ;


    return 0 ;
}
