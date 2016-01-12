// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Ruzena Chamrova <ruzena.chamrova@epfl.ch>, (C) 2007
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution

#include "gradient_sampler.h"

namespace Amie
{

double GradientSampler::getLinearDensity( Point p ) 
{
    Point p_ = direction.project( p ) ;
    return start + (end-start)*( p_-direction.first() ).norm() / direction.norm() ;
}

std::vector<Point> GradientSampler::sampleRectangleBoundingSurface( const Rectangle * geom, double linearDensity ) 
{
    std::vector<Point> pts ;
    if(!valid)
        return pts ;

    double size_x = geom->width() ;
    double size_y = geom->height() ;
    Point c = geom->getCenter() ;

    pts.push_back( Point( c.getX()-0.5*size_x, c.getY()-0.5*size_y ) ) ;
    double next = getLinearDensity( pts[0] )  / linearDensity ;
    std::cout << next << std::endl ;
    while( pts[ pts.size()-1 ].getX() < c.getX()+0.5*size_x - next*1.5 )
    {
        double nextx = pts[ pts.size()-1 ].getX() + next ;
        double nexty = pts[ pts.size()-1 ].getY() ;
        pts.push_back( Point( nextx, nexty ) ) ;
        next = getLinearDensity( pts[ pts.size()-1 ] )  / linearDensity ;
    }

    pts.push_back( Point( c.getX()+0.5*size_x, c.getY()-0.5*size_y ) ) ;
    next = getLinearDensity( pts[ pts.size()-1 ] )  / linearDensity ;
    std::cout << next << std::endl ;
    while( pts[ pts.size()-1 ].getY() < c.getY()+0.5*size_y - next*1.5 )
    {
        double nextx = pts[ pts.size()-1 ].getX() ;
        double nexty = pts[ pts.size()-1 ].getY() + next ;
        pts.push_back( Point( nextx, nexty ) ) ;
        next = getLinearDensity( pts[ pts.size()-1 ] )  / linearDensity ;
    }

    pts.push_back( Point( c.getX()+0.5*size_x, c.getY()+0.5*size_y ) ) ;
    next = getLinearDensity( pts[ pts.size()-1 ] )  / linearDensity ;
    std::cout << next << std::endl ;
    while( pts[ pts.size()-1 ].getX() > c.getX()-0.5*size_x + next*1.5 )
    {
        double nextx = pts[ pts.size()-1 ].getX() - next ;
        double nexty = pts[ pts.size()-1 ].getY() ;
        pts.push_back( Point( nextx, nexty ) ) ;
        next = getLinearDensity( pts[ pts.size()-1 ] )  / linearDensity ;
    }

    pts.push_back( Point( c.getX()-0.5*size_x, c.getY()+0.5*size_y ) ) ;
    next = getLinearDensity( pts[ pts.size()-1 ] )  / linearDensity ;
    std::cout << next << std::endl ;
    while( pts[ pts.size()-1 ].getY() > c.getY()-0.5*size_y + next*1.5 )
    {
        double nextx = pts[ pts.size()-1 ].getX() ;
        double nexty = pts[ pts.size()-1 ].getY() - next ;
        pts.push_back( Point( nextx, nexty ) ) ;
        next = getLinearDensity( pts[ pts.size()-1 ] )  / linearDensity ;
    }

    return pts ;
}

std::vector<Point> GradientSampler::sampleRectangleInnerSurface( const Rectangle * geom, double linearDensity ) 
{
    std::vector<Point> pts ;
    if(!valid)
        return pts ;

    Point dir = direction.vector()*(-1./direction.norm()) ;
    Point norm = direction.normal() ;

    double size_x = geom->width() ;
    double size_y = geom->height() ;
    Point c = geom->getCenter() ;

    Point start = Point(c.getX()-size_x*0.5, c.getY()-size_y*0.5) ;
    if( dir.getX() < 0 && dir.getY() >= 0)
        start = Point(c.getX()+size_x*0.5, c.getY()-size_y*0.5) ;
    if( dir.getX() < 0 && dir.getY() < 0)
        start = Point(c.getX()+size_x*0.5, c.getY()+size_y*0.5) ;
    if( dir.getX() >= 0 && dir.getY() < 0)
        start = Point(c.getX()-size_x*0.5, c.getY()+size_y*0.5) ;

    double next = getLinearDensity( start ) / linearDensity ;
    bool in = true ;
    int index = 0 ;
    start.print() ;
    start += dir*next ;
    while(in)
    {
        Line current( start, norm ) ;
        if(current.intersects(geom))
        {
            std::vector<Point> local = current.intersection( geom ) ;
            if( local.size() > 1 )
            {
                 Point localStart = local[index] ;
                 double localNext = getLinearDensity( localStart ) / linearDensity ;
                 double up = 1. ;
                 if(!geom->in(localStart + norm*localNext))
                     up = -1. ;
                 localStart += norm*(localNext*up) ;
                 while( geom->in(localStart) && dist(localStart, local[1-index]) > next*0.75 )
                 {
                     pts.push_back( localStart ) ;
                     localNext = getLinearDensity( localStart ) / linearDensity ;
                     localStart += norm*(localNext*up) ;
                 }
            }
            else
                in = false ;
        }
        else
            in = false ;
        index = 1-index ;
        next = getLinearDensity( start ) / linearDensity ;
        start += dir*next ;
    }

    return pts ;

}


}

