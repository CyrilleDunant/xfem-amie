
// Author: Jérôme Krebs <jerome.krebs@epfl.ch>, (C) 2007
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//

#include "placement.h"
#include "random.h"
#include "../geometry/geometry_base.h"


using namespace Amie ;


bool intersections( Feature * feat, const std::vector<Geometry *> & exclusionZones)
{
    for(size_t i = 0 ; i < exclusionZones.size() ; i++)
    {
        if(feat->intersects(exclusionZones[i]))
            return true ;
    }
    return false ;
}

bool isInside( Feature * feat, const std::vector<Geometry *> & base)
{
    for(size_t j = 0 ; j < base.size() ; j++)
    {
        if(feat->intersects(base[j]))
        {
            Point p = feat->getCenter() ;
            base[j]->project(&p) ;
            return !feat->in(p) ;
        }
        if(base[j]->in(feat->getCenter()))
            return true ;
    }
    return false ;
}

void transform2D( Feature * inc, RandomDistribution & xDistribution, RandomDistribution & yDistribution, RandomDistribution & rDistribution)
{
    Point c( xDistribution.draw(), yDistribution.draw() ) ;
    Point theta( 0,0, rDistribution.draw() ) ;
    inc->setCenter( c ) ;
    transform(inc, ROTATE, theta) ;
}

void transform3D( Feature * inc, RandomDistribution & xDistribution, RandomDistribution & yDistribution,  RandomDistribution & zDistribution)
{
    Point c( xDistribution.draw(), yDistribution.draw() , zDistribution.draw() ) ;
    inc->setCenter( c ) ;
}

std::vector<Feature *> Amie::placement2D(const Geometry* box, std::vector<Feature *> inclusions, double minDist, int placedAggregates, int triesMax, double orientation,  std::vector<Geometry *> exclusionZones)
{
    std::vector<Feature *> ret ;
    int tries = 0 ;

    std::vector<Point> boundingBox = box->getBoundingBox() ;
    UniformDistribution xDistribution( boundingBox[0].getX(), boundingBox[2].getX() ) ;
    UniformDistribution yDistribution( boundingBox[0].getY(), boundingBox[2].getY() ) ;
    UniformDistribution rDistribution( -orientation, orientation ) ;
    Grid grid(boundingBox[2].getX()-boundingBox[0].getX(), boundingBox[0].getY()-boundingBox[2].getY(), 10, box->getCenter()) ;

    for(int i = 0 ; i < placedAggregates ; i++)
    {
        ret.push_back(inclusions[i]);
        grid.add(inclusions[i]) ;
    }

    for(size_t i = placedAggregates ; i < inclusions.size() && tries < triesMax ; i++)
    {
        tries++ ;

        double scale = 1. ;
        if(minDist > POINT_TOLERANCE)
        {
            scale = (inclusions[i]->getRadius()+minDist)/inclusions[i]->getRadius() ;
            Point s( scale, scale ) ;
            transform(inclusions[i], SCALE, s) ;
        }

        transform2D( inclusions[i], xDistribution, yDistribution, rDistribution);
        std::vector<Point> bbox = inclusions[i]->getBoundingBox() ;
        while((!box->in(inclusions[i]->getCenter()) || !(box->in(bbox[0]) && box->in(bbox[1]) && box->in(bbox[2]) && box->in(bbox[3])) || (intersections(inclusions[i], exclusionZones))) && tries < triesMax)
        {
            tries++ ;
            transform2D( inclusions[i], xDistribution, yDistribution, rDistribution);
            bbox = inclusions[i]->getBoundingBox() ;
        }

        while(!grid.add(inclusions[i]) && tries < triesMax)
        {
            tries++ ;

            transform2D( inclusions[i], xDistribution, yDistribution, rDistribution);
            std::vector<Point> bbox = inclusions[i]->getBoundingBox() ;
            while((!box->in(inclusions[i]->getCenter()) || !(box->in(bbox[0]) && box->in(bbox[1]) && box->in(bbox[2]) && box->in(bbox[3]))|| (intersections(inclusions[i], exclusionZones))) && tries < triesMax )
            {
                tries++ ;
                transform2D( inclusions[i], xDistribution, yDistribution, rDistribution);
                bbox = inclusions[i]->getBoundingBox() ;
            }

        }

        if(tries < triesMax)
        {
            if(i%100 == 0)
                std::cout << "\rplaced " << i << " particles (tries " << tries << "/" << triesMax << ")" << std::flush ;
            if(scale > 1.)
            {
                Point s(1./scale, 1./scale) ;
                transform(inclusions[i], SCALE , s) ;
            }
            ret.push_back(inclusions[i]);
        }

    }

    std::cout << "\n" << ret.size() << " inclusions placed after " << tries << " tries" << std::endl ;
    double area = 0. ;
    for(size_t i = 0. ; i < ret.size() ; i++)
        area += ret[i]->area() ;
    std::cout << "surface covered by the inclusions: " << area << std::endl ;

    return ret ;
}

std::vector<Feature *> Amie::placement3D(const Geometry* box, std::vector<Feature *> inclusions, double minDist, int placedAggregates, int triesMax, double orientation,  std::vector<Geometry *> exclusionZones)
{
    std::vector<Feature *> ret ;
    int tries = 0 ;

    std::vector<Point> boundingBox = box->getBoundingBox() ;
    double min_x = box->getCenter().getX(), min_y = box->getCenter().getY(), max_x = box->getCenter().getX(), max_y = box->getCenter().getY(), max_z = box->getCenter().getZ(), min_z = box->getCenter().getZ();

    for ( size_t j  =  0 ; j <  boundingBox.size() ; j++ )
    {
        if ( boundingBox[j].getY() < min_y )
        {
            min_y = boundingBox[j].getY() ;
        }

        if ( boundingBox[j].getY() > max_y )
        {
            max_y = boundingBox[j].getY() ;
        }

        if ( boundingBox[j].getX() < min_x )
        {
            min_x = boundingBox[j].getX() ;
        }

        if ( boundingBox[j].getX() > max_x )
        {
            max_x = boundingBox[j].getX() ;
        }

        if ( boundingBox[j].getZ() < min_z )
        {
            min_z = boundingBox[j].getZ() ;
        }

        if ( boundingBox[j].getZ() > max_z )
        {
            max_z = boundingBox[j].getZ() ;
        }
    }
    UniformDistribution xDistribution( min_x, max_x ) ;
    UniformDistribution yDistribution( min_y, max_y ) ;
    UniformDistribution zDistribution( min_z, max_z ) ;
    Grid3D grid(max_x- min_x, max_y- min_y, max_z- min_z, 50, box->getCenter()) ;

    for(int i = 0 ; i < placedAggregates ; i++)
    {
        ret.push_back(inclusions[i]);
        grid.add(inclusions[i]) ;
    }

    for(size_t i = placedAggregates ; i < inclusions.size() && tries < triesMax ; i++)
    {
        tries++ ;

        double scale = 1. ;
        if(minDist > POINT_TOLERANCE)
        {
            scale = (inclusions[i]->getRadius()+minDist)/inclusions[i]->getRadius() ;
            Point s( scale, scale, scale ) ;
            transform(inclusions[i],SCALE, s) ;
        }

        transform3D( inclusions[i], xDistribution, yDistribution, zDistribution);

        while((box->intersects(inclusions[i]) || intersections(inclusions[i], exclusionZones)) && tries < triesMax)
        {
            tries++ ;
            transform3D( inclusions[i], xDistribution, yDistribution, zDistribution);
        }

        while(!grid.add(inclusions[i]) && tries < triesMax)
        {
            tries++ ;

            transform3D( inclusions[i], xDistribution, yDistribution, zDistribution);
            while((box->intersects(inclusions[i])|| intersections(inclusions[i], exclusionZones)) && tries < triesMax )
            {
                tries++ ;
                transform3D( inclusions[i], xDistribution, yDistribution, zDistribution);
            }

        }

        if(tries < triesMax)
        {
            if(i%100 == 0)
                std::cerr << "\rplaced " << i << " particles (tries " << tries << "/" << triesMax << ")" << std::flush ;
            if(scale > 1.)
            {
                Point s(1./scale, 1./scale, 1./scale) ;
                transform(inclusions[i], SCALE , s) ;
            }
            ret.push_back(inclusions[i]);
            tries = 0 ;
        }

    }

    std::cerr << "\n" << ret.size() << " inclusions placed after " << tries << " tries" << std::endl ;
    double area = 0. ;
    for(size_t i = 0. ; i < ret.size() ; i++)
        area += ret[i]->volume() ;
    std::cerr << "volume covered by the inclusions: " << area << std::endl ;

    return ret ;
}


std::vector<Feature *> Amie::placement2DInInclusions(const Geometry* box, std::vector<Geometry *> base, std::vector<Feature *> inclusions, double minDist, int placedAggregates, int triesMax, double orientation,  std::vector<Geometry *> exclusionZones)
{
    std::vector<Feature *> ret ;
    int tries = 0 ;

    std::vector<Point> boundingBox = box->getBoundingBox() ;
    UniformDistribution xDistribution( boundingBox[0].getX(), boundingBox[2].getX() ) ;
    UniformDistribution yDistribution( boundingBox[0].getY(), boundingBox[2].getY() ) ;
    UniformDistribution rDistribution( -orientation, orientation ) ;
    Grid grid(boundingBox[2].getX()-boundingBox[0].getX(), boundingBox[0].getY()-boundingBox[2].getY(), 10, box->getCenter()) ;

    for(int i = 0 ; i < placedAggregates ; i++)
    {
        ret.push_back(inclusions[i]);
        grid.add(inclusions[i]) ;
    }


    for(int i = placedAggregates ; i < (int)inclusions.size() && tries < triesMax ; i++)
    {
        tries++ ;

        double scale = 1. ;
        if(minDist > POINT_TOLERANCE)
        {
            scale = (inclusions[i]->getRadius()+minDist)/inclusions[i]->getRadius() ;
            Point s( scale, scale ) ;
            transform(inclusions[i], SCALE, s) ;
        }

        transform2D( inclusions[i], xDistribution, yDistribution, rDistribution);
        std::vector<Point> bbox = inclusions[i]->getBoundingBox() ;
        while(!box->in(inclusions[i]->getCenter()) || !(box->in(bbox[0]) && box->in(bbox[1]) && box->in(bbox[2]) && box->in(bbox[3])) || !isInside(inclusions[i], base) || (intersections(inclusions[i], exclusionZones) && tries < triesMax))
        {
            tries++ ;
            transform2D( inclusions[i], xDistribution, yDistribution, rDistribution);
            bbox = inclusions[i]->getBoundingBox() ;
        }

        while(!grid.add(inclusions[i]) && tries < triesMax)
        {
            tries++ ;

            transform2D( inclusions[i], xDistribution, yDistribution, rDistribution);
            std::vector<Point> bbox = inclusions[i]->getBoundingBox() ;
            while(!box->in(inclusions[i]->getCenter()) || !(box->in(bbox[0]) && box->in(bbox[1]) && box->in(bbox[2]) && box->in(bbox[3]))  || !isInside(inclusions[i], base) || (intersections(inclusions[i], exclusionZones) && tries < triesMax))
            {
                tries++ ;
                transform2D( inclusions[i], xDistribution, yDistribution, rDistribution);
                bbox = inclusions[i]->getBoundingBox() ;
            }

        }

        if(tries < triesMax)
        {
            if(i%100 == 0)
                std::cout << "\rplaced " << i << " particles (tries " << tries << "/" << triesMax << ")" << std::flush ;
            if(scale > 1.)
            {
                Point s(1./scale, 1./scale) ;
                transform(inclusions[i], SCALE , s) ;
            }
            ret.push_back(inclusions[i]);
        }


    }

    std::cout << "\n" << ret.size() << " inclusions placed after " << tries << " tries" << std::endl ;

    return ret ;
}



