// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Ruzena Chamrova <ruzena.chamrova@epfl.ch>, (C) 2007
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution

#include "sampler.h"
#include "../../features/feature_base.h"

namespace Amie
{

double Sampler::getEquivalentDensity( Feature * f, double linearDensity, double surfaceDensityFactor ) const 
{
    if(f->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
    {
        Rectangle * box ;
        if(dynamic_cast<Rectangle *>(f))
            box = new Rectangle( dynamic_cast<Rectangle *>(f)->width(), dynamic_cast<Rectangle *>(f)->height(), f->getCenter() ) ;
        else
            box = new Rectangle( f->getBoundingBox() ) ;
        if(box == nullptr || box->area() < 0)
            return -1 ;

        box->sampleSurface( linearDensity, surfaceDensityFactor ) ;
        double factor = sqrt( box->area() / ( box->getBoundingPoints().size() + box->getInPoints().size() ) ) ;
        delete box ;
        return factor ;
    }
    return -1 ;
}


std::vector<Point> Sampler::sampleBoundingSurface( const Geometry * geom, double linearDensity ) 
{
    std::vector<Point> pts ;
    switch( geom->getGeometryType() )
    {
    case CIRCLE:
        pts = this->sampleCircleBoundingSurface( dynamic_cast<const Circle *>(geom), linearDensity ) ;
        break ;
    case TRIANGLE:
        pts = this->sampleTriangleBoundingSurface( dynamic_cast<const Triangle *>(geom), linearDensity ) ;
        break ;
    case RECTANGLE:
        pts = this->sampleRectangleBoundingSurface( dynamic_cast<const Rectangle *>(geom), linearDensity ) ;
        break ;
    case POLYGON:
        pts = this->samplePolygonBoundingSurface( dynamic_cast<const Polygon *>(geom), linearDensity ) ;
        break ;
    case HEXAHEDRON:
        pts = this->sampleHexahedronBoundingSurface( dynamic_cast<const Hexahedron *>(geom), linearDensity ) ;
        break ;
    case SPHERE:
        pts = this->sampleSphereBoundingSurface( dynamic_cast<const Sphere *>(geom), linearDensity ) ;
        break ;
    case POLYGON_PRISM:
        pts = this->samplePolygonPrismBoundingSurface( dynamic_cast<const PolygonPrism *>(geom), linearDensity ) ;
        break ;
    case LOFTED_POLYGON:
        pts = this->sampleLoftedPolygonBoundingSurface( dynamic_cast<const LoftedPolygonPrism *>(geom), linearDensity ) ;
        break ;
    case ELLIPSE:
        pts = this->sampleEllipseBoundingSurface( dynamic_cast<const Ellipse *>(geom), linearDensity ) ;
        break ;
    default:
        break ;
    }

    if(pts.size() == 0)
        pts = geom->getSamplingBoundingPoints( linearDensity ) ;

    return pts ;
}

std::vector<Point> Sampler::sampleInnerSurface( const Geometry * geom, double linearDensity ) 
{
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


}

