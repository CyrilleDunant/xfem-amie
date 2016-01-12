// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Ruzena Chamrova <ruzena.chamrova@epfl.ch>, (C) 2007
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution

#include "regular_sampler.h"

namespace Amie
{

std::vector<Point> RegularSampler::sampleInnerSurface( const Geometry * geom, double linearDensity ) 
{
    if(force && geom->getGeometryType() != RECTANGLE && geom->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
    {
        std::vector<Point> bounds = geom->getBoundingBox() ;
	std::cout << "  " << bounds.size() << std::endl ;
        Rectangle * box = new Rectangle( bounds ) ;
        std::vector<Point> tries = this->sampleRectangleInnerSurface( box, linearDensity ) ;
        std::vector<Point> pts ;
        for(size_t i = 0 ; i < tries.size() ; i++)
        {
            if(geom->in(tries[i]))
               pts.push_back( tries[i] ) ;
        }
        delete box ;
        return pts ;
    }


    switch( geom->getGeometryType() )
    {
    case CIRCLE:
        return this->sampleCircleInnerSurface( dynamic_cast<const Circle *>(geom), linearDensity ) ;
    case TRIANGLE:
        return this->sampleTriangleInnerSurface( dynamic_cast<const Triangle *>(geom), linearDensity ) ;
    case RECTANGLE:
        return this->sampleRectangleInnerSurface( dynamic_cast<const Rectangle *>(geom), linearDensity ) ;
    case POLYGON:
        return this->samplePolygonInnerSurface( dynamic_cast<const Polygon *>(geom), linearDensity ) ;
    case HEXAHEDRON:
        return this->sampleHexahedronInnerSurface( dynamic_cast<const Hexahedron *>(geom), linearDensity ) ;
    case SPHERE:
        return this->sampleSphereInnerSurface( dynamic_cast<const Sphere *>(geom), linearDensity ) ;
    case POLYGON_PRISM:
        return this->samplePolygonPrismInnerSurface( dynamic_cast<const PolygonPrism *>(geom), linearDensity ) ;
    case LOFTED_POLYGON:
        return this->sampleLoftedPolygonInnerSurface( dynamic_cast<const LoftedPolygonPrism *>(geom), linearDensity ) ;
    case ELLIPSE:
        return this->sampleEllipseInnerSurface( dynamic_cast<const Ellipse *>(geom), linearDensity ) ;
    default:
        break ;
    }
    return std::vector<Point>() ;
}

std::vector<Point> RegularSampler::sampleRectangleBoundingSurface( const Rectangle * geom, double linearDensity ) 
{
    std::vector<Point> pts ;

    double size_x = geom->width() ;
    double size_y = geom->height() ;
    Point c = geom->getCenter() ;

    double perimeter = 2.*(size_x+size_y) ;
    double realDensity = std::max(linearDensity, 2./std::min(size_x,size_y)) ;
    size_t num_points = round(perimeter*realDensity) ;

    double distanceBetweenPointsx = std::min(perimeter/num_points, size_x) ;
    double distanceBetweenPointsy = std::min(perimeter/num_points, size_y) ;
    double dy = distanceBetweenPointsy ;
    double dx = distanceBetweenPointsx ;
    distanceBetweenPointsy = std::min(dx*1.5, dy) ;
    distanceBetweenPointsx = std::min(dy*1.5, dx) ;

    size_t numberOfPointsAlongX = round(1.5*static_cast<size_t>(std::ceil(size_x/distanceBetweenPointsx) + 1));
    double distanceBetweenPointsAlongX = size_x/(numberOfPointsAlongX-1) ;

    size_t numberOfPointsAlongY = round(1.5*static_cast<size_t>(std::ceil(size_y/distanceBetweenPointsy) + 1));
    double distanceBetweenPointsAlongY = size_y/(numberOfPointsAlongY-1) ;

    num_points = ((numberOfPointsAlongX)*2 + (numberOfPointsAlongY)*2 - 4) ;
    
    for (size_t i = 0 ; i < numberOfPointsAlongY; i++)
        pts.push_back(Point(c.getX()-0.5*size_x, c.getY() + 0.5*size_y - i*distanceBetweenPointsAlongY) );

    for (size_t i = 1 ; i < numberOfPointsAlongX ; i++)
        pts.push_back(Point(c.getX()-0.5*size_x + i*distanceBetweenPointsAlongX, c.getY() - 0.5*size_y ) ) ;

    for (size_t i = 1 ; i < numberOfPointsAlongY ; i++)
        pts.push_back(Point(c.getX()+0.5*size_x, c.getY() - 0.5*size_y + i*distanceBetweenPointsAlongY) );

    for (size_t i = 1 ; i < numberOfPointsAlongX-1 ; i++)
        pts.push_back(Point(c.getX()+0.5*size_x- i*distanceBetweenPointsAlongX, c.getY() + 0.5*size_y ) );

    return pts ;
}

std::vector<Point> RegularSampler::sampleRectangleInnerSurface( const Rectangle * geom, double linearDensity ) 
{
    std::vector<Point> pts ;

    double size_x = geom->width() ;
    double size_y = geom->height() ;
    Point c = geom->getCenter() ;

    double perimeter = 2.*(size_x+size_y) ;
    double realDensity = std::max(linearDensity, 2./std::min(size_x,size_y)) ;
    size_t num_points = round(perimeter*realDensity) ;

    double distanceBetweenPointsx = std::min(perimeter/num_points, size_x) ;
    double distanceBetweenPointsy = std::min(perimeter/num_points, size_y) ;
    double dy = distanceBetweenPointsy ;
    double dx = distanceBetweenPointsx ;
    distanceBetweenPointsy = std::min(dx*1.5, dy) ;
    distanceBetweenPointsx = std::min(dy*1.5, dx) ;

    size_t numberOfPointsAlongX = round(1.5*static_cast<size_t>(std::ceil(size_x/distanceBetweenPointsx) + 1));
    double distanceBetweenPointsAlongX = size_x/(numberOfPointsAlongX-1) ;

    size_t numberOfPointsAlongY = round(1.5*static_cast<size_t>(std::ceil(size_y/distanceBetweenPointsy) + 1));
    double distanceBetweenPointsAlongY = size_y/(numberOfPointsAlongY-1) ;

    num_points = ((numberOfPointsAlongX)*2 + (numberOfPointsAlongY)*2 - 4) ;
    
    for (size_t i = 1 ; i < numberOfPointsAlongY-1; i++)
    {
        for (size_t j = 1 ; j < numberOfPointsAlongX-1; j++)
        {
            pts.push_back( Point( c.getX()-0.5*size_x + j*distanceBetweenPointsAlongX, c.getY()-0.5*size_y + i*distanceBetweenPointsAlongY ) ) ;
        }
    }

    return pts ;

}


}

