
// Author: Jérôme Krebs <jerome.krebs@epfl.ch>, (C) 2007
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//

#include "placement.h"
#include "../geometry/geometry_base.h"


using namespace Amie ;


bool intersections( Feature * feat, const std::vector<Geometry *> & exclusionZones)
{
    for(size_t i = 0 ; i < exclusionZones.size() ; i++)
    {
        if(feat->intersects(exclusionZones[i]) || exclusionZones[i]->in(feat->getCenter()))
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

void projectOnEdge( Feature * feat, const std::vector<Geometry *> & base, bool vertex)
{
    for(size_t j = 0 ; j < base.size() ; j++)
    {
        if(feat->intersects(base[j]) || base[j]->in( feat->getCenter()) )
        {
            Point p = feat->getCenter() ;
            if(vertex && base[j]->getGeometryType() == POLYGON)
            {
                std::valarray<Point> pts = dynamic_cast<Polygon *>(base[j])->getOriginalPoints() ;
                size_t k = 0 ;
                double d = dist(p, pts[0] ) ;
                for(size_t i = 1 ; i < pts.size() ; i++)
                {
                    if( dist(p, pts[i]) < d)
                    {
                        d = dist( p, pts[i] ) ;
                        k = i ;
                    }
                }
                p = pts[k] ;
            }
            else
                base[j]->project(&p) ;

            feat->setCenter(p) ;
            return ;
        }
    }
}

void transform2D( Feature * inc, double x, double y, double r)
{
    Point c( x,y ) ;
    Point theta( 0,0, r ) ;
    inc->setCenter( c ) ;
    transform(inc, ROTATE, theta) ;
}

void transform3D( Feature * inc, double x, double y, double z)
{
    Point c( x,y,z ) ;
    inc->setCenter( c ) ;
}

std::vector<Feature *> Amie::placement2D(const Geometry* box, std::vector<Feature *> inclusions, double minDist, int placedAggregates, int triesMax, double orientation,  std::vector<Geometry *> exclusionZones)
{
    std::vector<Feature *> ret ;
    int tries = 0 ;

    std::vector<Point> boundingBox = box->getBoundingBox() ;
    std::default_random_engine generator ;
    std::uniform_real_distribution< double > xDistribution( boundingBox[0].getX(), boundingBox[2].getX() ) ;
    std::uniform_real_distribution< double > yDistribution( boundingBox[0].getY(), boundingBox[2].getY() ) ;
    std::uniform_real_distribution< double > rDistribution( -orientation, orientation ) ;
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

        transform2D( inclusions[i], xDistribution(generator), yDistribution(generator), rDistribution(generator));
        std::vector<Point> bbox = inclusions[i]->getBoundingBox() ;
        while((!box->in(inclusions[i]->getCenter()) || !(box->in(bbox[0]) && box->in(bbox[1]) && box->in(bbox[2]) && box->in(bbox[3])) || (intersections(inclusions[i], exclusionZones))) && tries < triesMax)
        {
            tries++ ;
            transform2D( inclusions[i], xDistribution(generator), yDistribution(generator), rDistribution(generator));
            bbox = inclusions[i]->getBoundingBox() ;
        }

        while(!grid.add(inclusions[i]) && tries < triesMax)
        {
            tries++ ;

            transform2D( inclusions[i], xDistribution(generator), yDistribution(generator), rDistribution(generator));
            std::vector<Point> bbox = inclusions[i]->getBoundingBox() ;
            while((!box->in(inclusions[i]->getCenter()) || !(box->in(bbox[0]) && box->in(bbox[1]) && box->in(bbox[2]) && box->in(bbox[3]))|| (intersections(inclusions[i], exclusionZones))) && tries < triesMax )
            {
                tries++ ;
                transform2D( inclusions[i], xDistribution(generator), yDistribution(generator), rDistribution(generator));
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
    std::default_random_engine generator ;
    std::uniform_real_distribution< double > xDistribution( min_x, max_x ) ;
    std::uniform_real_distribution< double > yDistribution( min_y, max_y ) ;
    std::uniform_real_distribution< double > zDistribution( min_z, max_z ) ;
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

        transform3D( inclusions[i], xDistribution(generator), yDistribution(generator), zDistribution(generator));

        while((box->intersects(inclusions[i]) || intersections(inclusions[i], exclusionZones)) && tries < triesMax)
        {
            tries++ ;
            transform3D( inclusions[i], xDistribution(generator), yDistribution(generator), zDistribution(generator));
        }

        while(!grid.add(inclusions[i]) && tries < triesMax)
        {
            tries++ ;

            transform3D( inclusions[i], xDistribution(generator), yDistribution(generator), zDistribution(generator));
            while((box->intersects(inclusions[i])|| intersections(inclusions[i], exclusionZones)) && tries < triesMax )
            {
                tries++ ;
                transform3D( inclusions[i], xDistribution(generator), yDistribution(generator), zDistribution(generator));
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
    std::default_random_engine generator ;
    std::uniform_real_distribution< double > xDistribution( boundingBox[0].getX(), boundingBox[2].getX() ) ;
    std::uniform_real_distribution< double > yDistribution( boundingBox[0].getY(), boundingBox[2].getY() ) ;
    std::uniform_real_distribution< double > rDistribution( -orientation, orientation ) ;
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

        transform2D( inclusions[i], xDistribution(generator), yDistribution(generator), rDistribution(generator));
        std::vector<Point> bbox = inclusions[i]->getBoundingBox() ;
        while(!box->in(inclusions[i]->getCenter()) || !(box->in(bbox[0]) && box->in(bbox[1]) && box->in(bbox[2]) && box->in(bbox[3])) || !isInside(inclusions[i], base) || (intersections(inclusions[i], exclusionZones) && tries < triesMax))
        {
            tries++ ;
            transform2D( inclusions[i], xDistribution(generator), yDistribution(generator), rDistribution(generator));
            bbox = inclusions[i]->getBoundingBox() ;
        }

        while(!grid.add(inclusions[i]) && tries < triesMax)
        {
            tries++ ;

            transform2D( inclusions[i], xDistribution(generator), yDistribution(generator), rDistribution(generator));
            std::vector<Point> bbox = inclusions[i]->getBoundingBox() ;
            while(!box->in(inclusions[i]->getCenter()) || !(box->in(bbox[0]) && box->in(bbox[1]) && box->in(bbox[2]) && box->in(bbox[3]))  || !isInside(inclusions[i], base) || (intersections(inclusions[i], exclusionZones) && tries < triesMax))
            {
                tries++ ;
                transform2D( inclusions[i], xDistribution(generator), yDistribution(generator), rDistribution(generator));
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

std::vector<Feature *> Amie::placement2DOnEdge(const Geometry* box, std::vector<Geometry *> base, std::vector<Feature *> inclusions, bool vertex, double minDist, int placedAggregates, int triesMax, double orientation,  std::vector<Geometry *> exclusionZones)
{
    std::vector<Feature *> ret ;
    int tries = 0 ;

    std::vector<Point> boundingBox = box->getBoundingBox() ;
    std::default_random_engine generator ;
    std::uniform_real_distribution< double > xDistribution( boundingBox[0].getX(), boundingBox[2].getX() ) ;
    std::uniform_real_distribution< double > yDistribution( boundingBox[0].getY(), boundingBox[2].getY() ) ;
    std::uniform_real_distribution< double > rDistribution( -orientation, orientation ) ;
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

        transform2D( inclusions[i], xDistribution(generator), yDistribution(generator), rDistribution(generator));
        projectOnEdge( inclusions[i], base, vertex) ;
        std::vector<Point> bbox = inclusions[i]->getBoundingBox() ;
        while((!box->in(inclusions[i]->getCenter()) || !(box->in(bbox[0]) && box->in(bbox[1]) && box->in(bbox[2]) && box->in(bbox[3]))) && tries < triesMax)
        {
            tries++ ;
            transform2D( inclusions[i], xDistribution(generator), yDistribution(generator), rDistribution(generator));
            projectOnEdge( inclusions[i], base, vertex) ;
            bbox = inclusions[i]->getBoundingBox() ;
        }

        while(!grid.add(inclusions[i]) && tries < triesMax)
        {
            tries++ ;

            transform2D( inclusions[i], xDistribution(generator), yDistribution(generator), rDistribution(generator));
            projectOnEdge( inclusions[i], base, vertex) ;
            std::vector<Point> bbox = inclusions[i]->getBoundingBox() ;
            while((!box->in(inclusions[i]->getCenter()) || !(box->in(bbox[0]) && box->in(bbox[1]) && box->in(bbox[2]) && box->in(bbox[3]))) && tries < triesMax)
            {
                tries++ ;
                transform2D( inclusions[i], xDistribution(generator), yDistribution(generator), rDistribution(generator));
                projectOnEdge( inclusions[i], base, vertex) ;
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


