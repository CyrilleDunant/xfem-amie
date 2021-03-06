//
// C++ Implementation: granulo
//
// Description:
//
//
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2007-2011
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include <cmath>
#include <vector>
#include <cstring>
#include "granulo.h"
#include <iostream>  // I/O 
#include <fstream>   // file I/O
#include "placement.h"
#include "configuration.h"
#include "parser/parser.h"
#include "../geometry/geometry_base.h"
#include "../features/sample.h"
#include "../features/microstructuregenerator.h"

using namespace Amie ;

PSDGenerator::PSDGenerator()
{

}

// ParticleSizeDistribution * ParticleSizeDistribution::getPSD(PSDType type)
// {
// 	switch(type)
// 	{
// 	  case BOLOME_A:
// 	    return new PSDBolomeA() ;
// 	  case BOLOME_B:
// 	    return new PSDBolomeB() ;
// 	  case BOLOME_C:
// 	    return new PSDBolomeC() ;
// 	  case BOLOME_D:
// 	    return new PSDBolomeD() ;
// 	  case PSD_UNIFORM:
// 			return new ParticleSizeDistribution() ;
// 	}
// 	return new ParticleSizeDistribution() ;
// }
//
std::vector<Inclusion *> PSDGenerator::get2DInclusions(double rmax, double mass, ParticleSizeDistribution * type, PSDEndCriteria crit)
{
    std::vector<double> radii ;
    double diameter = rmax*2. ;
    double remainingMass = mass ;
    double remainingFraction = 1. ;
    int tries = 0 ;
    while(!crit.meets(diameter*0.5, remainingFraction, tries++))
    {
        radii.push_back(diameter*0.5) ;
        remainingMass -= diameter*diameter*M_PI*.25 ;
        if(remainingMass < 0)
        {
            radii.pop_back() ;
            break ;
        }
        remainingFraction = remainingMass / mass ;
        diameter = type->getNext2DDiameter(diameter, remainingFraction, rmax*2.) ;
    }
// 	crit.print(diameter*0.5, remainingFraction, radii.size()) ;

    std::sort(radii.begin(), radii.end()) ;
    std::reverse(radii.begin(), radii.end());

    std::vector<Inclusion *> incs ;
    if(radii.size() == 0)
    {
        return incs ;
    }
    for(size_t i = 0 ; i < radii.size()-1 ; i++)
    {
        incs.push_back(new Inclusion(radii[i],0.,0.)) ;
    }
    return incs ;
}

std::vector<Inclusion3D *> PSDGenerator::get3DInclusions(double rmax, double mass, ParticleSizeDistribution * type, PSDEndCriteria crit)
{
    ParticleSizeDistribution * psd = type ;
    std::vector<double> radii ;
    double diameter = rmax*2. ;
    double remainingMass = mass ;
    double remainingFraction = 1. ;
    int tries = 0 ;

    while(!crit.meets(diameter*0.5, remainingFraction, tries++))
    {
        radii.push_back(diameter*0.5) ;
        remainingMass -= diameter*diameter*diameter*M_PI*1./6. ;
        if(remainingMass < 0)
        {
            radii.pop_back() ;
            break ;
        }
        remainingFraction = remainingMass / mass ;
        diameter = psd->getNext3DDiameter(diameter, remainingFraction, rmax*2.) ;
    }

    std::sort(radii.begin(), radii.end()) ;
    std::reverse(radii.begin(), radii.end());

    std::vector<Inclusion3D *> incs ;
    for(size_t i = 0 ; i < radii.size() ; i++)
        incs.push_back(new Inclusion3D(radii[i],0.,0.,0.)) ;
    delete psd ;
    return incs ;
}

std::vector<Inclusion *> PSDGenerator::get2DMortar(double rmax, double width, size_t n, ParticleSizeDistribution * type)
{
    if(!type)
        type = new PSDBolomeD() ;
    return get2DInclusions(rmax, width*width*0.65, type, PSDEndCriteria(0.5*0.15e-3, 0.01, n)) ;
}

std::vector<Inclusion *> PSDGenerator::get2DConcrete(double rmax, double width, size_t n, ParticleSizeDistribution * type, double percent)
{
    if(!type)
        type = new PSDBolomeA() ;
    return get2DInclusions(rmax, width*width*percent, type, PSDEndCriteria(0.5*0.15e-3, 0.01, n)) ;
}

std::vector<Feature *> PSDGenerator::get2DConcrete(FeatureTree * F, Form * behaviour, size_t n, double rmax, double itz, ParticleSizeDistribution * type, InclusionGenerator * geometry, size_t tries,double fraction, Geometry * placement, std::vector<Geometry *> exclusionZones, size_t seed)
{
    InclusionGenerator * converter = (geometry ? geometry : new InclusionGenerator()  ) ;


    if(!type)
        type = new PSDBolomeA() ;
    Feature * box = F->getFeature(0) ;
    double ar = sqrt(box->area()) ;
    if(placement != nullptr)
        ar = sqrt(placement->area()) ;

    std::vector<Inclusion *> inc = get2DConcrete(rmax, ar, n, type,fraction) ;
    std::vector<Feature *> feats = converter->convert( inc ) ;
    std::vector<Feature *> real ;
    inc.clear() ;

    if(placement)
        feats = placement2D( placement, feats, itz, 0, tries, converter->authorizeRotationsDuringPlacement, exclusionZones, seed ) ;
    else
        feats = placement2D( dynamic_cast<Rectangle *>(box), feats, itz, 0, tries, converter->authorizeRotationsDuringPlacement, exclusionZones, seed ) ;
    double area = 0 ;
    for(size_t i = 0 ; i < feats.size() ; i++)
    {
        area += feats[i]->area() ;
        if(behaviour)
            feats[i]->setBehaviour(behaviour) ;
        real.push_back(feats[i]) ;
// 		feats[i]->isVirtualFeature = true ;
//		feats[i]->isUpdated = false ;
        F->addFeature(box, feats[i]) ;
    }
    return real ;
}

std::vector<Feature *> PSDGenerator::get2DInclusionsOnEdge(FeatureTree * F, Form * behaviour, std::vector<Feature *> base, bool checkMask, bool vertex, size_t n, double rmax, double itz, ParticleSizeDistribution * type, InclusionGenerator * geometry, size_t tries,double fraction, Geometry * placement, std::vector<Geometry *> exclusionZones, size_t seed)
{
    InclusionGenerator * converter = (geometry ? geometry : new InclusionGenerator()  ) ;

    if(!type)
        type = new PSDBolomeA() ;
    Feature * box = F->getFeature(0) ;
    double ar = sqrt(box->area()) ;
    if(placement != nullptr)
        ar = sqrt(placement->area()) ;

    std::vector<Inclusion *> inc = get2DConcrete(rmax, ar, n, type,fraction) ;
    std::vector<Feature *> feats = converter->convert( inc ) ;
    std::vector<Feature *> ret ;
    inc.clear() ;

    std::vector<Geometry *> geom ;
    for(size_t i = 0 ; i < base.size() ; i++)
        geom.push_back( dynamic_cast<Geometry *>(base[i])) ;
    if(placement)
        feats = placement2DOnEdge( placement, geom, feats, vertex, itz, 0, tries, converter->authorizeRotationsDuringPlacement, exclusionZones, seed ) ;
    else
        feats = placement2DOnEdge( dynamic_cast<Rectangle *>(box), geom, feats, vertex, itz, 0, tries, converter->authorizeRotationsDuringPlacement, exclusionZones, seed ) ;
    double area = 0 ;
    for(size_t i = 0 ; i < feats.size() ; i++)
    {
        bool inMask = true ;
        if(checkMask)
        {
            for(size_t j = 0 ; j < base.size() ; j++)
            {
                if( base[j]->in( feats[i]->getCenter() ) || base[j]->intersects( feats[i] ) )
                {
                    inMask = base[j]->inMask( feats[i]->getCenter() ) ;
                    break ;
                }
            }
        }

        if(inMask)
        {
            area += feats[i]->area() ;
            if(behaviour)
                feats[i]->setBehaviour(behaviour) ;
            for(size_t j = 0 ; j < base.size() ; j++)
            {
                if(base[j]->in(feats[i]->getCenter()))
                {
                    F->addFeature(base[j], feats[i]) ;
                    ret.push_back( feats[i] ) ;
                    break ;
                }
            }
        }
    }
    return ret ;
}

std::vector<Feature *> PSDGenerator::get2DEmbeddedInclusions(FeatureTree * F, Form * behaviour, std::vector<Feature *> base, size_t n, double rmax, double itz, ParticleSizeDistribution * type, InclusionGenerator * geometry, size_t tries,double fraction, Geometry * placement, std::vector<Geometry *> exclusionZones, size_t seed)
{
    InclusionGenerator * converter = (geometry ? geometry : new InclusionGenerator()  ) ;

    if(!type)
        type = new PSDBolomeA() ;
    Feature * box = F->getFeature(0) ;
    double ar = sqrt(box->area()) ;
    if(placement != nullptr)
        ar = sqrt(placement->area()) ;

    std::vector<Inclusion *> inc = get2DConcrete(rmax, ar, n, type,fraction) ;
    std::vector<Feature *> feats = converter->convert( inc ) ;
    inc.clear() ;

    std::vector<Geometry *> geom ;
    for(size_t i = 0 ; i < base.size() ; i++)
        geom.push_back( dynamic_cast<Geometry *>(base[i])) ;
    if(placement)
        feats = placement2DInInclusions( placement, geom, feats, itz, 0, tries, converter->authorizeRotationsDuringPlacement, exclusionZones, seed ) ;
    else
        feats = placement2DInInclusions( dynamic_cast<Rectangle *>(box), geom, feats, itz, 0, tries, converter->authorizeRotationsDuringPlacement, exclusionZones, seed ) ;
    double area = 0 ;
    for(size_t i = 0 ; i < feats.size() ; i++)
    {
        area += feats[i]->area() ;
        if(behaviour)
            feats[i]->setBehaviour(behaviour) ;
        for(size_t j = 0 ; j < base.size() ; j++)
        {
            if(base[j]->in(feats[i]->getCenter()))
            {
                F->addFeature(base[j], feats[i]) ;
                break ;
            }
        }
    }
    return feats ;
}

std::vector<Feature *> PSDGenerator::get2DMaskedInclusions(FeatureTree * F, Form * behaviour, std::vector<Feature *> base, size_t n, double rmax, double itz, ParticleSizeDistribution * type, InclusionGenerator * geometry, size_t tries,double fraction, Geometry * placement, std::vector<Geometry *> exclusionZones, size_t seed)
{
    InclusionGenerator * converter = (geometry ? geometry : new InclusionGenerator()  ) ;

    if(!type)
        type = new PSDBolomeA() ;
    Feature * box = F->getFeature(0) ;
    double ar = sqrt(box->area()) ;
    if(placement != nullptr)
        ar = sqrt(placement->area()) ;

    std::vector<Inclusion *> inc = get2DConcrete(rmax, ar, n, type,fraction) ;
    std::vector<Feature *> feats = converter->convert( inc ) ;
    std::vector<Feature *> ret ;
    inc.clear() ;

    if(placement)
        feats = placement2D( placement, feats, itz, 0, tries, converter->authorizeRotationsDuringPlacement, exclusionZones, seed ) ;
    else
        feats = placement2D( dynamic_cast<Rectangle *>(box), feats, itz, 0, tries, converter->authorizeRotationsDuringPlacement, exclusionZones, seed ) ;
    for(size_t i = 0 ; i < feats.size() ; i++)
    {
        std::vector< Feature *> potentials ;
        if(behaviour)
            feats[i]->setBehaviour(behaviour) ;
        for(size_t j = 0 ; j < base.size() ; j++)
        {
            if(base[j]->intersects(feats[i]) || base[j]->in( feats[i]->getCenter() ) )
                potentials.push_back( base[j] ) ;
        }
        if(potentials.size() > 0)
        {
            F->addFeature( potentials[0], feats[i] ) ;
            feats[i]->setMask( potentials ) ;
            ret.push_back( feats[i] ) ;
        }
    }
    return ret ;
}

double getMaximumRadius( Feature * f, std::vector<Geometry *> geom)
{
    if(geom.size() == 0)
        return f->getRadius() ;

    double d = dist( f->getCenter(), geom[0]->getCenter() )-geom[0]->getRadius() ;
    for(size_t i = 1 ; i < geom.size() ; i++)
    {
        double di = dist( f->getCenter(), geom[i]->getCenter() )-geom[i]->getRadius() ;
        if(di < d)
            d = di ;
    }
    return d ;
}

std::vector<std::vector<PolygonalSample *> > PSDGenerator::get2DSourceVoronoiPolygons(Rectangle * box, std::vector<VoronoiGrain> & morphology, size_t truen, double minDist, double delta, size_t seed)
{
    std::vector<Geometry *> exclusion ;
    std::vector<std::pair<Point, size_t > > grains ;

    double r0 = (minDist < 0 ? morphology[morphology.size()-1].radius : minDist) ;

    std::vector<Inclusion *> inclusions = PSDGenerator::get2DInclusions( r0*0.5, box->area(), new ConstantSizeDistribution(), PSDEndCriteria( r0*0.25, box->area(), truen) ) ;
    std::vector<Feature *> circles ;
    for(size_t i = 0 ; i < inclusions.size() ; i++)
        circles.push_back(inclusions[i]) ;
    double tries = 100 ;
/*    if(morphology.size() == 1)
    {
       tries *= 10 ;
    }*/
    std::vector<Feature *> placed = placement2D(box, circles, 0,0, tries*truen, M_PI, std::vector<Geometry* >(), seed) ;
    std::random_shuffle( placed.begin(), placed.end() ) ;

    std::vector<int> index(placed.size(), -1) ;
    size_t j = 0 ;
    size_t i = 0 ;
    double current = 0. ;
    double radius = morphology[i].radius ;
    double correctionFactor = morphology[i].correctionFactor ;
    double area = box->area()*morphology[i].fraction ;
    size_t count = 0 ;
    double target_count = morphology[i].n ; 

    while(j < index.size())
    {
        if(index[j] > -1)
        {
            j++ ;
            continue ;
        }

        Point p = placed[j]->getCenter() ;
        Circle c( radius, p ) ;
        int count_in = 0 ;
        int count_out = 0 ;
        if(radius > r0)
        {
            for(size_t k = 0 ; k < index.size() ; k++)
            {
                if( dist(p,placed[k]->getCenter()) < radius )
                {
                    if(index[k] > -1)
                        count_out++ ;
                    else
                    {
                        index[k] = i ;
                        count_in++ ;
                    }
                }
            }
            count++ ;
        }
        else
        {
            index[j] = i ;
            count_in = 1 ;
            count++ ;
        }

        if(radius > r0)
        {
            double totalcount = /*correctionFactor*/(count_in+count_out) ;
            current += c.area()*correctionFactor*((double) count_in)/((double) totalcount) ;
        }
//        else
//            current += c.area() ;
        grains.push_back( std::make_pair( p, i ) ) ;
        j++ ;
       

        if( (morphology.size() > 1) && ((radius > r0 && current > area ) || (count > target_count) ))
        {
            i++ ;
            if( i < morphology.size())
            {
                correctionFactor = morphology[i].correctionFactor ;
                radius = morphology[i].radius ;
                area = box->area()*morphology[i].fraction ;
                current = 0. ;
                count = 0 ;
                target_count = morphology[i].fraction*placed.size()*r0/(radius) ; 
            }
            else
            {
                break ;
            }
        }
        else if(morphology.size() == 1 && grains.size() > target_count)
            break ;


    }

    std::vector<std::vector<PolygonalSample *> > poly( morphology.size() ) ;
    if(grains.size() < 2)
        return poly ;

    double n = grains.size() ;
    Mesh<DelaunayTriangle, DelaunayTreeItem> * test = new DelaunayTree ( &box->getBoundingPoint(0), &box->getBoundingPoint(1), &box->getBoundingPoint(2) ) ;
    test->insert( &box->getBoundingPoint(3), 0 ) ;
    double divx = box->width()/sqrt(n) ;
    Point c = box->getCenter() ;
    for(double x = c.getX()-box->width()*0.5+divx ; x < c.getX()+box->width()*0.5-divx/2 ; x += divx )
    {
        test->insert( new Point( x, c.getY()-box->height()*0.5 ), 0 ) ;
        test->insert( new Point( x, c.getY()+box->height()*0.5 ), 0 ) ;
    }
    double divy = box->height()/sqrt(n) ;
    for(double y = c.getY()-box->height()*0.5+divx ; y < c.getY()+box->height()*0.5-divy/2 ; y += divy )
    {
        test->insert( new Point( c.getX()-box->width()*0.5, y ), 0 ) ;
        test->insert( new Point( c.getX()+box->width()*0.5, y ), 0 ) ;
    }
    std::vector<Point *> nodes ;
    for(size_t i = 0 ; i < grains.size() ; i++)
    {
        nodes.push_back( new Point(grains[i].first) ) ;
        test->insert( nodes[i], 0 ) ;
    }

    std::vector<DelaunayTriangle *> connectivity = test->getConflictingElements( box ) ;
    std::cout << connectivity.size() << " elements in voronoi sub-mesh" << std::endl ;
    size_t discarded = 0 ;

    std::vector<Point> nextvertex ;
    for(size_t i = 0 ; i < connectivity.size() ; i++)
    {
        double w = 0 ;
        Point r ;
        Vector wk( connectivity[i]->getBoundingPoints().size() ) ;
        for(size_t j = 0 ; j < connectivity[i]->getBoundingPoints().size() ; j++)
        {
            int k = -1 ;
            wk[j] = r0 ;
            for(size_t l = 0 ; l < nodes.size() ; l++)
            {
                if( dist(connectivity[i]->getBoundingPoint(j), *(nodes[l])) < POINT_TOLERANCE )
                {
                    k = l ;
                    break ;
                }
            }
            if(k > -1)
                 wk[j] = morphology[grains[k].second].weight(r0, connectivity[i]->getBoundingPoint(j), connectivity[i]->getCenter() ) ;
            w += wk[j] ;
        }
        double rw = 0 ;
        for(size_t j = 0 ; j < connectivity[i]->getBoundingPoints().size() ; j++)
        {
           r += connectivity[i]->getBoundingPoint(j)*(1.-wk[j]/w) ;
           rw += (1.-wk[j]/w) ;
        }
        r /= rw ;
        nextvertex.push_back(r) ;
    }

    for(size_t i = 0 ; i < nodes.size() ; i++)
    {
        std::vector< std::pair<Point, std::pair< Point,  Point> > > next ;
        for(size_t j = 0 ; j < connectivity.size() ; j++)
        {
            int vertex = -1 ;
            for(size_t k = 0 ; k < connectivity[j]->getBoundingPoints().size() ; k++)
            {
                if( dist(connectivity[j]->getBoundingPoint(k), *(nodes[i])) < POINT_TOLERANCE )
                {
                    vertex = k ;
                    break ;
                }
            }
            if(vertex > -1)
            {
                Point a ;
                Point b ;
                bool first = true ;
                for(size_t k = 0 ; k < connectivity[j]->getBoundingPoints().size() ; k++)
                {
                    if( k != (size_t) vertex )
                    {
                        if(first)
                        {
                            a = connectivity[j]->getBoundingPoint(k) ;
                            first = false ;
                        }
                        else
                            b = connectivity[j]->getBoundingPoint(k) ;
                    }
                }
                next.push_back( std::make_pair( nextvertex[ j ], std::make_pair( a-connectivity[j]->getBoundingPoint(vertex), b-connectivity[j]->getBoundingPoint(vertex) ) ) ) ;
            }
        }
        if(next.size() < 3)
        {
            discarded++ ;
            continue ;
        }
        std::vector<size_t> vertex ;
        std::valarray<bool> done(next.size()) ;
        done = false ;
        size_t start = 0 ;
        while(vertex.size() < next.size())
        {
            bool found = false ;
            if(vertex.size() == 0)
            {
                vertex.push_back( start ) ;
                done[start] = true ;
                found = true ;
            }
            for(size_t j = 0 ; j < next.size() ; j++)
            {
                if(done[j])
                    continue ;
                if(dist( next[j].second.first, next[ vertex[vertex.size()-1] ].second.first ) < POINT_TOLERANCE ||
                        dist( next[j].second.first, next[ vertex[vertex.size()-1] ].second.second ) < POINT_TOLERANCE ||
                        dist( next[j].second.second, next[ vertex[vertex.size()-1] ].second.first ) < POINT_TOLERANCE ||
                        dist( next[j].second.second, next[ vertex[vertex.size()-1] ].second.second ) < POINT_TOLERANCE)
                {
                    vertex.push_back(j) ;
                    found = true ;
                    done[j] = true ;
                }
            }
            if(!found)
            {
                vertex.clear() ;
                done = false ;
                start++ ;
            }
            if(start == next.size())
                break ;
        }
        bool open = true ;
        if(next.size() > 16)
            std::cout << next.size() << std::endl ;
        for(size_t j = 0 ; j < next.size() ; j++)
        {
            if(j == vertex.size()-1 || j == vertex.size()-2)
                continue ;
            if(dist( next[j].second.first, next[ vertex[vertex.size()-1] ].second.first ) < POINT_TOLERANCE ||
                    dist( next[j].second.first, next[ vertex[vertex.size()-1] ].second.second ) < POINT_TOLERANCE ||
                    dist( next[j].second.second, next[ vertex[vertex.size()-1] ].second.first ) < POINT_TOLERANCE ||
                    dist( next[j].second.second, next[ vertex[vertex.size()-1] ].second.second ) < POINT_TOLERANCE)
            {
                open=false ;
            }
        }
        if(open && next.size() > 3)
        {
            discarded++ ;
            continue ;
        }
        std::valarray<Point *> corners(vertex.size()) ;
        for(size_t j = 0 ; j < vertex.size() ; j++)
        {
            corners[j] = new Point( next[vertex[j]].first ) ;
        }
        std::vector<Segment> edges ;
        for(size_t j = 0 ; j < vertex.size() ; j++)
        {
            Segment s( next[vertex[j]].first, next[vertex[(j+1)%vertex.size()]].first ) ;
            edges.push_back( s ) ;
        }

        if(delta > POINT_TOLERANCE)
        {
            Polygon p(corners) ;
            std::vector<Polygon> inscr = p.getInscribedPolygons( delta ) ;
            for(size_t j = 0 ; j < inscr.size() ; j++)
            {
                std::valarray<Point> opts = inscr[j].getOriginalPoints() ;
                std::valarray<Point *> pts( opts.size()  ) ;
                for(size_t k = 0 ; k < opts.size() ; k++)
                    pts[k] = new Point( opts[k] ) ;
                poly[grains[i].second].push_back( new PolygonalSample( nullptr, pts ) ) ;
            }
        }
        else
        {
            poly[grains[i].second].push_back( new PolygonalSample( nullptr, corners) ) ;
        }
    }
    delete test ;
    return poly ;
}

std::vector<std::vector<Feature *> > PSDGenerator::get2DVoronoiPolygons(FeatureTree * F, std::vector<VoronoiGrain> & grains, size_t n, double minDist, double border, size_t nmax, bool copy, double delta, size_t seed)
{
    RectangularFeature * sample = dynamic_cast<RectangularFeature *>(F->getFeature(0)) ;
    Rectangle * placement = new Rectangle( sample->width()+minDist*2.+border*2., sample->height()+minDist*2.+border*2., sample->getCenter().getX(), sample->getCenter().getY() ) ;
    double realn = placement->area()/(minDist*minDist*M_PI) ;
    std::vector<std::vector<PolygonalSample *> > poly = PSDGenerator::get2DSourceVoronoiPolygons( placement, grains, (n==0?realn : n), minDist, delta, seed) ;
    std::vector<std::vector<Feature *> > ret( std::max(1, (int) grains.size()) ) ;

    for(size_t i = 0 ; i < poly.size() ; i++)
    {
        for(size_t j = 0 ; j < poly[i].size() ; j++)
        {
            poly[i][j]->setBehaviour( copy ? grains[i].behaviour->getCopy() : grains[i].behaviour ) ;
            F->addFeature( sample, poly[i][j] ) ;
            ret[i].push_back(poly[i][j]) ;
        }
    }
    delete placement ;
    return ret ;
}

std::vector<std::vector<Feature *> > PSDGenerator::get2DVoronoiPolygons(Rectangle * placement, std::vector<VoronoiGrain> & grains, size_t n, double minDist, double border, size_t nmax, bool copy, double delta, size_t seed)
{
    double r0 = minDist ;
    if(r0 < 0)
        r0 = grains[grains.size()-1].radius ;
    Rectangle * realplacement = new Rectangle( placement->width()+r0*2.+border*2., placement->height()+r0*2.+border*2., placement->getCenter().getX(), placement->getCenter().getY() ) ;
    double realn = realplacement->area()/(r0*r0*M_PI) ;
    std::vector<std::vector<PolygonalSample *> > poly = PSDGenerator::get2DSourceVoronoiPolygons( realplacement, grains, (n==0?realn : n), minDist, delta, seed) ;
    std::vector<std::vector<Feature *> > ret( std::max(1, (int) grains.size()) ) ;

    for(size_t i = 0 ; i < poly.size() ; i++)
    {
        for(size_t j = 0 ; j < poly[i].size() ; j++)
        {
            poly[i][j]->setBehaviour( copy ? grains[i].behaviour->getCopy() : grains[i].behaviour ) ;
//            F->addFeature( sample, poly[i][j] ) ;
            ret[i].push_back(poly[i][j]) ;
        }
    }
    delete realplacement ;
    return ret ;
}

std::vector<std::vector<Feature *> > PSDGenerator::get2DVoronoiPolygons(Feature * feat, std::vector<VoronoiGrain> & grains, size_t n, double minDist, double border, size_t nmax, bool copy, double delta, size_t seed)
{
    std::vector<Point> box = feat->getBoundingBox() ;
    Rectangle * placement = new Rectangle( box ) ;
    Rectangle * realbox = new Rectangle( placement->width()+minDist+border*2., placement->height()+minDist+border*2., placement->getCenter().getX(), placement->getCenter().getY() ) ;
    std::vector<std::vector<PolygonalSample *> > poly = PSDGenerator::get2DSourceVoronoiPolygons( realbox, grains, n, minDist, delta, seed) ;
    std::vector<std::vector<Feature *> > ret( std::max(1, (int) grains.size()) ) ;
    for(size_t i = 0 ; i < poly.size() ; i++)
    {
        for(size_t j = 0 ; j < poly[i].size() ; j++)
        {
            if( poly[i][j]->getOriginalPoints().size() < nmax+1 && (feat->in(poly[i][j]->getCenter()) || feat->intersects( dynamic_cast<Polygon *>(poly[i][j]) ) ))
            {
                poly[i][j]->setBehaviour(copy ? grains[i].behaviour->getCopy() : grains[i].behaviour) ;
                poly[i][j]->addToMask( feat ) ;
                ret[i].push_back(poly[i][j]) ;
            }
        }
    }
    return ret ;
}

std::vector<std::vector<Feature *> > PSDGenerator::get2DVoronoiPolygons(FeatureTree * F, std::vector<VoronoiGrain> & grains, std::vector<Feature *> feats, size_t n, double minDist, double border, size_t nmax, bool copy, double delta, size_t seed)
{
    std::vector<std::vector<Feature *> > ret( std::max(1, (int) grains.size()) ) ;
    double realn = (n == 0 ? F->getFeature(0)->area()/(minDist*minDist*M_PI) : n) ;
    for(size_t i = 0 ; i < feats.size() ; i++)
    {
        double num = ((double) realn)*sqrt(feats[i]->area()/F->getFeature(0)->area()) ;
        if(num > 1 && feats[i]->getRadius() > minDist*1.5)
        {
            std::vector<std::vector<Feature *> > poly = PSDGenerator::get2DVoronoiPolygons( feats[i], grains, std::max(4., num), minDist, border, nmax, copy, delta, seed) ;
            for(size_t j = 0 ; j < poly.size() ; j++)
            {
                for(size_t k = 0 ; k < poly[j].size() ; k++)
                {
                    F->addFeature(feats[i], poly[j][k]) ;
                    ret[j].push_back(poly[j][k]) ;
                }
            }
        }
    }

    return ret ;
}

std::vector<std::vector<Feature *> > PSDGenerator::get2DVoronoiPolygons(Rectangle * box, std::vector<VoronoiGrain> & grains, std::vector<Feature *> feats, size_t n, double minDist, double border, size_t nmax, bool copy, double delta, size_t seed)
{
    std::vector<std::vector<Feature *> > ret( std::max(1, (int) grains.size()) ) ;
    double r0 = minDist ;
    if(r0 < 0) { r0 = grains[grains.size()-1].radius ; }
    double realn = (n == 0 ? box->area()/(r0*r0*M_PI) : n) ;
    for(size_t i = 0 ; i < feats.size() ; i++)
    {
        double num = ((double) realn)*sqrt(feats[i]->area()/box->area()) ;
        if(num > 1 && feats[i]->getRadius() > r0*1.5)
        {
            std::vector<std::vector<Feature *> > poly = PSDGenerator::get2DVoronoiPolygons( feats[i], grains, std::max(4., num), r0, border, nmax, copy, delta, seed) ;
            for(size_t j = 0 ; j < poly.size() ; j++)
            {
                for(size_t k = 0 ; k < poly[j].size() ; k++)
                {
//                    F->addFeature(feats[i], poly[j][k]) ;
                    poly[j][k]->setFather( feats[i] ) ;
                    ret[j].push_back(poly[j][k]) ;
                }
            }
        }
    }

    return ret ;
}

std::vector<Inclusion *> PSDGenerator::get2DMortar(FeatureTree * F, Form * behaviour, double rmax, size_t n, ParticleSizeDistribution * type, size_t tries, size_t seed)
{
    if(!type)
        type = new PSDBolomeD() ;
    Feature * box = F->getFeature(0) ;
    std::vector<Inclusion *> inc = get2DMortar(rmax, sqrt(box->area()), n, type) ;
    std::vector<Feature *> feats ;
    for(size_t i = 0 ; i < inc.size() ; i++)
    {
        feats.push_back(inc[i]) ;
    }
    inc.clear() ;
    srand(seed) ;
    feats = placement2D( dynamic_cast<Rectangle *>(box), feats, 0.00001, 0, tries, M_PI, std::vector<Geometry *>(), seed ) ;
    for(size_t i = 0 ; i < feats.size() ; i++)
    {
        inc.push_back(dynamic_cast<Inclusion *>(feats[i])) ;
    }
    for(size_t i = 0 ; i < inc.size() ; i++)
    {
        inc[i]->setBehaviour(behaviour) ;
        F->addFeature(box, inc[i]) ;
    }
    return inc ;
}

std::vector<std::pair<ExpansiveZone *, Inclusion *> > PSDGenerator::get2DExpansiveZonesInAggregates(FeatureTree * F, std::vector<Inclusion *> incs, StiffnessWithImposedStrain * behaviour, double radius, size_t n, size_t max, int maxPerAgg)
{
    Feature * box = F->getFeature(0) ;
    RectangularFeature * sample = dynamic_cast<RectangularFeature *>(box) ;
    double w = sample->width()*0.5-radius*60 ;
    double h = sample->height()*0.5-radius*60 ;
    std::default_random_engine gen ;
    std::uniform_real_distribution< double > distributionx(-w,w) ;
    std::uniform_real_distribution< double > distributiony(-h,h) ;
    std::vector<std::pair<ExpansiveZone *, Inclusion *> > ret ;
    double aggregateArea = 0 ;

    std::vector<ExpansiveZone *> zonesToPlace ;

    for(size_t i = 0 ; i < n ; i++)
    {
        Point pos( distributionx(gen), distributiony(gen) ) ;
        pos += sample->getCenter() ;
        bool alone  = true ;
        for(size_t j = 0 ; j< zonesToPlace.size() ; j++)
        {
            if (squareDist(pos, zonesToPlace[j]->Circle::getCenter()) < (radius*60.+radius*60.)*(radius*60.+radius*60.))
            {
                alone = false ;
                break ;
            }
        }
        if (alone)
            zonesToPlace.push_back(new ExpansiveZone(nullptr, radius, pos.getX(), pos.getY(), behaviour)) ;
    }

    std::map<Inclusion *, int> zonesPerIncs ;
    for(size_t i = 0 ; i < zonesToPlace.size() ; i++)
    {
        bool placed = false ;
        if(ret.size() < max)
        {
            for(size_t j = 0 ; j < incs.size() ; j++)
            {
                Circle circle(incs[j]->getRadius()*0.95 , incs[j]->getCenter()) ;
                if(circle.in(zonesToPlace[i]->getCenter()))
                {
                    if((maxPerAgg < 0 || zonesPerIncs[incs[j]] < maxPerAgg) )
                    {
                        zonesPerIncs[incs[j]]++ ;
                        F->addFeature(incs[j],zonesToPlace[i]) ;
//					  F->addPoint(&(zonesToPlace[i]->getCenter())) ;
                        ret.push_back(std::make_pair(zonesToPlace[i],incs[j])) ;
                        placed = true ;
                        break ;
                    }
                }
            }
        }
        if(!placed)
            delete zonesToPlace[i] ;

    }

// 	exit(0) ;
    int count = 0 ;
    for(auto i = zonesPerIncs.begin() ; i != zonesPerIncs.end() ; ++i)
    {
        aggregateArea+= i->first->area() ;
        count+= i->second ;
    }

    std::cerr << ret.size() << " zones placed on reactive aggregate area of " << aggregateArea << std::endl ;

    return ret ;

}

std::vector<std::pair<TimeDependentHomogenisingInclusion *, Inclusion *> > PSDGenerator::get2DGrowingExpansiveZonesInAggregates(FeatureTree * F, std::vector<Inclusion *> incs, ViscoelasticityAndImposedDeformation * behaviour, Function radius, double rmax, size_t n, size_t max, int maxPerAgg)
{
    Feature * box = F->getFeature(0) ;
    RectangularFeature * sample = dynamic_cast<RectangularFeature *>(box) ;
    double w = sample->width()*0.5 - 2*rmax ;
    double h = sample->height()*0.5 - 2*rmax ;
    std::default_random_engine gen ;
    std::uniform_real_distribution< double > distributionx(-w,w) ;
    std::uniform_real_distribution< double > distributiony(-h,h) ;
    std::vector<std::pair<TimeDependentHomogenisingInclusion *, Inclusion *> > ret ;
    double aggregateArea = 0 ;

    std::vector<TimeDependentHomogenisingInclusion *> zonesToPlace ;

    for(size_t i = 0 ; i < n ; i++)
    {
        Point pos( distributionx(gen), distributiony(gen) ) ;
        pos += sample->getCenter() ;
        bool alone  = true ;
        for(size_t j = 0 ; j< zonesToPlace.size() ; j++)
        {
            if (squareDist(pos, zonesToPlace[j]->getCenter()) < (rmax*rmax*9))
            {
                alone = false ;
                break ;
            }
        }
        if (alone)
            zonesToPlace.push_back(new TimeDependentHomogenisingInclusion(nullptr, radius, pos.getX(), pos.getY(), behaviour)) ;
    }

    std::map<Inclusion *, int> zonesPerIncs ;
    for(size_t i = 0 ; i < zonesToPlace.size() ; i++)
    {
        bool placed = false ;
        if(ret.size() < max)
        {
            for(size_t j = 0 ; j < incs.size() ; j++)
            {
                if(incs[j]->getRadius() > rmax*4.)
                {
                    Circle circle(incs[j]->getRadius() - rmax*3. , incs[j]->getCenter()) ;
                    if(circle.in(zonesToPlace[i]->getCenter()))
                    {
                        if((maxPerAgg < 0 || zonesPerIncs[incs[j]] < maxPerAgg) )
                        {
                            zonesPerIncs[incs[j]]++ ;
                            F->addFeature(incs[j],zonesToPlace[i]) ;
                            ret.push_back(std::make_pair(zonesToPlace[i],incs[j])) ;
                            placed = true ;
                            break ;
                        }
                    }
                }
            }
        }
        if(!placed)
            delete zonesToPlace[i] ;

    }

// 	exit(0) ;
    int count = 0 ;
    for(auto i = zonesPerIncs.begin() ; i != zonesPerIncs.end() ; ++i)
    {
        aggregateArea+= i->first->area() ;
        count+= i->second ;
    }

    std::cerr << ret.size() << " zones placed on reactive aggregate area of " << aggregateArea << std::endl ;

    return ret ;

}


double PSDBolomeA::getNext2DDiameter(double diameter, double fraction, double dmax)
{
    double b = -(4.*fraction+1.) ;
    double delta = 8.*fraction + 1. ;
    return std::max(15.e-5,(- b - std::sqrt(delta))/2.*dmax/**2./M_PI*/) ;

}

double PSDBolomeA::getNext3DDiameter(double diameter, double fraction, double dmax)
{
    double b = 1.+fraction/.25 ;
    return dmax*(b - std::sqrt(b*b-fraction*fraction/(0.25*0.25))) ;
}

double PSDBolomeB::getNext2DDiameter(double diameter, double fraction, double dmax)
{
    return fraction*fraction*dmax ;
}

double PSDBolomeB::getNext3DDiameter(double diameter, double fraction, double dmax)
{
    return fraction*fraction*dmax ;
}

double PSDBolomeC::getNext2DDiameter(double diameter, double fraction, double dmax)
{
    double next = fraction*fraction*dmax ;
    if(next > 0.004)
        next *= 1.05 ;
    return next ;
}

double PSDBolomeC::getNext3DDiameter(double diameter, double fraction, double dmax)
{
    double next = fraction*fraction*dmax ;
    if(next > 0.004)
        next *= 1.05 ;
    return next ;
}

double PSDBolomeD::getNext2DDiameter(double diameter, double fraction, double dmax)
{
    double m;//pente
    double b;//ordonn??e en x=0
    if (fraction> 0. && fraction<0.1)
    {
        Point A(.000150,0.);
        Point B(.000315,0.1);
        m= (B.getY()-A.getY())/(B.getX()-A.getX());
        b=A.getY()-(m*A.getX());
        return std::max(15.e-05, (fraction-b)/m/**2./M_PI*/) ;
    }
    if (fraction<0.2 && fraction>= 0.1)
    {
        Point A(.000315,0.1);
        Point B(.00063,0.2);
        m= (B.getY()-A.getY())/(B.getX()-A.getX());
        b=A.getY()-(m*A.getX());
        return (fraction-b)/m/**2./M_PI*/ ;
    }
    if (fraction<0.45 && fraction>= 0.2)
    {
        Point A(.00063,0.2);
        Point B(.00125,0.45);
        m= (B.getY()-A.getY())/(B.getX()-A.getX());
        b=A.getY()-(m*A.getX());
        return (fraction-b)/m/**2./M_PI*/ ;
    }
    if (fraction<0.7 && fraction>= 0.45)
    {
        Point A(.00125,0.45);
        Point B(.0025,0.70);
        m= (B.getY()-A.getY())/(B.getX()-A.getX());
        b=A.getY()-(m*A.getX());
        return (fraction-b)/m/**2./M_PI*/;
    }
    if (fraction<1. && fraction>= 0.7)
    {
        Point A(.0025,.7);
        Point B(.005,1.);
        m= (B.getY()-A.getY())/(B.getX()-A.getX());
        b=A.getY()-(m*A.getX());
        return (fraction-b)/m/**2./M_PI*/ ;
    }
    return 15e-5 ;
}

double PSDBolomeD::getNext3DDiameter(double diameter, double fraction, double dmax)
{
    double m;//pente
    double b;//ordonn??e en x=0
    if (fraction> 0. && fraction<0.1)
    {
        Point A(.000150,0.);
        Point B(.000315,0.1);
        m= (B.getY()-A.getY())/(B.getX()-A.getX());
        b=A.getY()-(m*A.getX());
        return std::max(7.5e-05, (fraction-b)/m);
    }
    if (fraction<0.2 && fraction>= 0.1)
    {
        Point A(.000315,0.1);
        Point B(.00063,0.2);
        m= (B.getY()-A.getY())/(B.getX()-A.getX());
        b=A.getY()-(m*A.getX());
        return (fraction-b)/m*.5;
    }
    if (fraction<0.45 && fraction> 0.2)
    {
        Point A(.00063,0.2);
        Point B(.00125,0.45);
        m= (B.getY()-A.getY())/(B.getX()-A.getX());
        b=A.getY()-(m*A.getX());
        return (fraction-b)/m;
    }
    if (fraction<0.7 && fraction> 0.45)
    {
        Point A(.00125,0.45);
        Point B(.0025,0.70);
        m= (B.getY()-A.getY())/(B.getX()-A.getX());
        b=A.getY()-(m*A.getX());
        return (fraction-b)/m;
    }
    if (fraction<1. && fraction> 0.7)
    {
        Point A(.0025,.7);
        Point B(.005,1.);
        m= (B.getY()-A.getY())/(B.getX()-A.getX());
        b=A.getY()-(m*A.getX());
        return (fraction-b)/m;
    }
    return 7.5e-5 ;

}


double PSDFuller::getNext2DDiameter(double diameter, double fraction, double dmax)
{
    if(diameter < dmin + POINT_TOLERANCE)
        return dmin ;
    return dmin + (dmax-dmin)*std::pow(fraction, 1./exponent) ;
}

double PSDFuller::getNext3DDiameter(double diameter, double fraction, double dmax)
{
    if(diameter < dmin + POINT_TOLERANCE)
        return dmin ;
    return dmin + (dmax-dmin)*std::pow(fraction, 1./exponent) ;
}

GranuloFromFile::GranuloFromFile(const std::string & fname, std::vector<std::string> columns)
{
    std::cerr << "importing file: " << fname << std::endl ;
    filename = fname ;
    for(size_t i = 0 ; i < columns.size() ; i++)
        fields.push_back(columns[i]) ;

    std::fstream filereader ;
    filereader.open(this->filename.c_str(),std::ios::in);

    while(!filereader.eof())
    {
        double buff ;
        filereader >> buff ;
        values.push_back(buff) ;
    }
    filereader.close() ;
    std::cerr << "done..." << std::endl ;
}

bool GranuloFromFile::verifyField(std::vector<std::string> columns)
{
    for(size_t i = 0 ; i < columns.size() ; i++)
    {
        if(getFieldNumber(columns[i]) == -1)
            return false ;
    }
    return true ;
}

int GranuloFromFile::getFieldNumber(std::string column)
{
    for(size_t j = 0 ; j < fields.size() ; j++)
    {
        if(strcmp(column.c_str(),fields[j].c_str()) == 0)
        {
            return j ;
        }
    }
    std::cerr << "error: field not found" << std::endl ;
    return -1 ;
}

std::vector<double> GranuloFromFile::getFieldValues(std::string column)
{
    std::vector<double> val ;
    int f = getFieldNumber(column) ;
    if(f > -1)
        val = getFieldValues(f) ;
    return val ;
}

std::vector<double> GranuloFromFile::getFieldValues(int f)
{
    std::vector<double> val ;
    int nv = values.size() ;
    int nf = fields.size() ;
    int i = 0 ;
    while(i * nf < nv)
    {
        val.push_back(values[i * nf + f]) ;
        i++ ;
    }
    return val ;
}

std::vector<Feature *> GranuloFromFile::getFeatures(GeometryType type, int ninc)
{
    std::vector<Feature *> inc ;
    std::vector<std::string> columns ;
    switch(type)
    {
    case CIRCLE:
        // inclusions
        columns.push_back("radius") ;
        columns.push_back("center_x") ;
        columns.push_back("center_y") ;
        break ;
    case SPHERE:
        // inclusions 3D
        columns.push_back("radius") ;
        columns.push_back("center_x") ;
        columns.push_back("center_y") ;
        columns.push_back("center_z") ;
        break ;
    case ELLIPSE:
        // ellipses
        columns.push_back("radius_a") ;
        columns.push_back("radius_b") ;
        columns.push_back("center_x") ;
        columns.push_back("center_y") ;
        columns.push_back("axis_x") ;
        columns.push_back("axis_y") ;
        break ;
    default:
        break ;
    }
    std::vector< std::vector<double> > fieldvalues ;
    std::cout << "extracting columns" ;
    for(size_t i = 0 ; i < columns.size() ; i++)
    {
        std::cout << " ... " << columns[i] ;
        std::vector<double> val = this->getFieldValues(columns[i]) ;
        fieldvalues.push_back(val) ;
    }
    std::cout << std::endl ;
    switch(type)
    {
    case CIRCLE:
        // inclusions
        std::cout << "creating inclusions..." << std::endl ;
        for(size_t i = 0 ; i < fieldvalues[0].size() && (int)i < ninc+1 ; i++)
            inc.push_back(new Inclusion(fieldvalues[0][i], fieldvalues[1][i], fieldvalues[2][i])) ;
        break ;
    case SPHERE:
        // inclusions 3D
        std::cout << "creating 3D inclusions..." << std::endl ;
        for(size_t i = 0 ; i < fieldvalues[0].size() && (int)i < ninc+1 ; i++)
            inc.push_back(new Inclusion3D(fieldvalues[0][i], fieldvalues[1][i], fieldvalues[2][i], fieldvalues[3][i])) ;
        break ;
    case ELLIPSE:
        // ellipses
        std::cout << "creating ellipses..." << std::endl ;
        for(size_t i = 0 ; i < fieldvalues[0].size() && (int)i < ninc+1 ; i++)
        {
            Point center(fieldvalues[2][i], fieldvalues[3][i]) ;
            Point a(fieldvalues[4][i], fieldvalues[5][i]) ;
            a /= a.norm() ;
            a *= fieldvalues[0][i] ;
            double b = fieldvalues[1][i]/fieldvalues[0][i] ;
            inc.push_back(new EllipsoidalInclusion(center,a,b)) ;
        }
        break ;
    default:
        break;
    }
    inc.pop_back() ;
    return inc ;
}


std::vector<Inclusion3D *> GranuloFromFile::getInclusion3D(int ninc, double scale)
{
    std::vector<Inclusion3D *> inc ;
    std::vector<std::string> columns ;
    columns.push_back("radius") ;
    columns.push_back("center_x") ;
    columns.push_back("center_y") ;
    columns.push_back("center_z") ;
    std::vector< std::vector<double> > fieldvalues ;
    std::cout << "extracting columns" ;
    for(size_t i = 0 ; i < columns.size() ; i++)
    {
        std::cout << " ... " << columns[i] ;
        std::vector<double> val = this->getFieldValues(columns[i]) ;
        fieldvalues.push_back(val) ;
    }
    std::cout << "creating 3D inclusions..." << std::endl ;
    for(size_t i = 0 ; i < fieldvalues[0].size() && (int)i < ninc ; i++)
        inc.push_back(new Inclusion3D(fieldvalues[0][i]*scale, fieldvalues[1][i]*scale, fieldvalues[2][i]*scale, fieldvalues[3][i]*scale)) ;
    std::cout << "done" << std::endl ;
    inc.pop_back() ;
    return inc ;
}

std::vector<Feature *> PolygonGranuloFromFile::getFeatures(SpaceDimensionality dim, Feature * father, int fields)
{
    std::vector<Feature *> ret ;
    for(int j = 0 ; j < fields ; j++)
        data.push_back( std::vector<double>() ) ;
    size_t i = 0 ;
    while(i < values.size())
    {
        int sides = values[i] ; if(sides < 0) { sides = 0 ; }
        std::valarray<Point *> pts( sides ) ;
        size_t ndim = (dim == SPACE_TWO_DIMENSIONAL ? 2 : 3) ;
        bool valid = (sides>2) && i+fields+sides*ndim < values.size() ;
        if(!valid)
        {
            i += fields+sides*ndim+1 ;
            continue ;
        }
        for(int j = 0 ; j < fields ; j++)
        {
            data[j].push_back( values[i+j+1] ) ;
        }
        for(int j = 0 ; j < sides  ; j++)
        {
           if(dim == SPACE_TWO_DIMENSIONAL)
               pts[j] = new Point( values[i+fields+j*ndim+1], values[i+fields+j*ndim+2] ) ;
           else
               pts[j] = new Point( values[i+fields+j*ndim+1], values[i+fields+j*ndim+2], values[i+fields+j*ndim+3] ) ;
        }
        if(dim == SPACE_TWO_DIMENSIONAL)
            ret.push_back(new PolygonalSample( father, pts) ) ;
        i += fields+sides*ndim+1 ;
    }

    return ret ;
}


GranuloFromCumulativePSD::GranuloFromCumulativePSD(const std::string & filename, PSDSpecificationType t, double factor, double cutOffUp, double cutOffDown)
{
    std::fstream file(filename) ;

    if(!file.is_open())
    {
        std::cout << "file " << filename << " doesn't exist!" << std::endl ;
        exit(0) ;
    }
    double maxfrac = 0 ;
    double minfrac = 100 ;
    do {
        double frac ;
        double rad ;
        file >> frac >> rad ;

        if(((cutOffUp > 0 && rad*factor <= cutOffUp) || cutOffUp < 0)  &&
                ((cutOffDown > 0 && rad*factor >= cutOffDown) || cutOffDown < 0)
          )
        {
// 			if(frac < 1e-13 )
// 			{
// 				foundZeroFrac = true ;
// 				zeroFracRad = std::max(rad*factor, zeroFracRad) ;
// 				minfrac = std::min(minfrac, frac) ;
// 			}
// 			else
// 			{
            fraction.push_back(frac);
            radius.push_back(rad*factor);
            maxfrac = std::max(maxfrac, frac) ;
            minfrac = std::min(minfrac, frac) ;
// 			}
        }
    } while(!file.eof()) ;

    fraction.pop_back();
    radius.pop_back();

    //collate possible list of zero-radius or fraction
// 	if(radius.front() > radius.back())
// 	{
// 		fraction.push_back(0);
// 		radius.push_back(zeroFracRad);
// 	}
// 	else
// 	{
// 		fraction.insert(fraction.begin(), 0.);
// 		radius.insert(radius.begin(), zeroFracRad) ;
// 	}

    if(cutOffUp > 0  || cutOffDown > 0)
    {
        if(t == CUMULATIVE_PERCENT)
            t = CUMULATIVE_ABSOLUTE ;
        if(t == CUMULATIVE_PERCENT_REVERSE)
            t = CUMULATIVE_ABSOLUTE_REVERSE ;
    }

    for(size_t i = 0 ; i < fraction.size() ; i++)
    {
        fraction[i] -= minfrac ;
    }

    switch(t)
    {
    case CUMULATIVE_PERCENT:
    {
        for(size_t i = 0 ; i < fraction.size() ; i++)
            fraction[i] /= 100. ;
        break ;
    }
    case CUMULATIVE_FRACTION:
    {
        break ;
    }
    case CUMULATIVE_ABSOLUTE:
    {
        double maxv = std::max(fraction.front(), fraction.back()) ;
        for(size_t i = 0 ; i < fraction.size() ; i++)
            fraction[i] /= maxv ;
        break ;
    }
    case CUMULATIVE_PERCENT_REVERSE:
    {
        for(size_t i = 0 ; i < fraction.size() ; i++)
            fraction[i] = (100. - fraction[i])/100. ;
        break ;
    }
    case CUMULATIVE_FRACTION_REVERSE:
    {
        for(size_t i = 0 ; i < fraction.size() ; i++)
            fraction[i] = 1. - fraction[i] ;
        break ;
    }
    case CUMULATIVE_ABSOLUTE_REVERSE:
    {
        double maxv = std::max(fraction.front(), fraction.back()) ;
        for(size_t i = 0 ; i < fraction.size() ; i++)
        {
            fraction[i] /= maxv ;
            fraction[i] = 1. - fraction[i] ;
        }
        break ;
    }
    }
    if(fraction.back() < fraction.front())
    {
        std::reverse(fraction.begin(), fraction.end());
        std::reverse(radius.begin(), radius.end());
    }
//
// 	exit(0) ;

}

double GranuloFromCumulativePSD::getNext2DDiameter(double diameter, double frac, double dmax)
{

    if(diameter < 2.*radius.front()+POINT_TOLERANCE)
    {
        return 2.*radius.front() ;
    }

    for( size_t i = 0  ; i < radius.size()-1 ; i++)
    {
        if(fraction[i] > frac )
        {
            double df = (exp(frac)-exp(fraction[i-1]))/(exp(fraction[i])-exp(fraction[i-1])) ;
            double v = std::max(2.*(radius[i-1]*(1.-df) + radius[i]*df), 2.*radius.front())  ;
            return v ;
        }
    }
    double df = (exp(frac)-exp(fraction[radius.size()-2]))/(exp(fraction[radius.size()-1])-exp(fraction[radius.size()-2])) ;
    double v = std::max(2.*(radius[radius.size()-2]*(1.-df) + radius[radius.size()-1]*df), 2.*radius.front())  ;
    return v ;


}
double GranuloFromCumulativePSD::getNext3DDiameter(double diameter, double frac, double dmax)
{
    if(diameter < 2.*radius.front()+POINT_TOLERANCE)
    {
        return 2.*radius.front() ;
    }

    for( size_t i = 0  ; i < radius.size()-1 ; i++)
    {
        if(fraction[i] > frac )
        {
            double df = (frac-fraction[i-1])/(fraction[i]-fraction[i-1]) ;
            double v = std::max(2.*(radius[i-1]*(1.-df) + radius[i]*df), 2.*radius.front())  ;
            return v ;
        }
    }
    double df = (frac-fraction[radius.size()-2])/(fraction[radius.size()-1]-fraction[radius.size()-2]) ;
    double v = std::max(2.*(radius[radius.size()-2]*(1.-df) + radius[radius.size()-1]*df), 2.*radius.front())  ;
    return v ;
}



