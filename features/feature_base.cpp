// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//


#include "feature_base.h"
#include "../physics/void_form.h"
#include "../mesher/delaunay.h"
#include "../mesher/delaunay_3d.h"
#include "features.h"

using namespace Amie ;


Feature::Feature(Feature * father) :  m_f(father), behaviour(nullptr)
{
    isEnrichmentFeature = false ;
    isCompositeFeature = false ;
    isVirtualFeature = false ;
    isUpdated = false ;
    behaviourSource = nullptr ;

    layer = -1 ;
    if(father != nullptr)
        father->addChild(this) ;

}

std::vector<Point *> Feature::doubleSurfaceSampling()
{
    std::vector<Point *> ret ;
    Amie::PointArray newboundingPoints(getBoundingPoints().size()*2) ;
    for(size_t i = 0 ; i < getBoundingPoints().size()-1 ; i++)
    {
        newboundingPoints[i*2] = &getBoundingPoint(i) ;
        Point * p = new Point(getBoundingPoint(i)*0.5 + getBoundingPoint(i+1)*0.5) ;
        project(p) ;
        p->setId(-1) ;
        newboundingPoints[i*2+1] = p ;
        ret.push_back(newboundingPoints[i*2+1]) ;
    }

    newboundingPoints[(getBoundingPoints().size()-1)*2] = &getBoundingPoint(getBoundingPoints().size()-1) ;
    newboundingPoints[getBoundingPoints().size()*2-1] = new Point(getBoundingPoint(getBoundingPoints().size()-1)*0.5 + getBoundingPoint(0)*0.5) ;
    newboundingPoints[getBoundingPoints().size()*2-1]->setId(-1) ;
    project(newboundingPoints[getBoundingPoints().size()*2-1]) ;
    ret.push_back(newboundingPoints[getBoundingPoints().size()*2-1]) ;

    dynamic_cast<Geometry *>(this)->setBoundingPoints(newboundingPoints) ;

    return ret ;
}


const Feature * Feature::getBehaviourSource() const
{
    return behaviourSource ;
}

Feature * Feature::getBehaviourSource()
{
    return behaviourSource ;
}

void Feature::setBehaviourSource(Feature * const f)
{
    behaviourSource = f ;
}

void  Feature::addChild( Feature*const f) 
{
    if(getChildren().empty() || std::find(getChildren().begin(), getChildren().end(), f) == getChildren().end())
    {
        getChildren().push_back(f) ;
    }
}

void  Feature::removeChild(Feature * const feat)
{
    auto c = std::find(m_c.begin(), m_c.end(), feat) ;
    if(c != m_c.end())
        getChildren().erase(c) ;
}

void Feature::setBehaviour(Form * const f)
{
    /*	if(behaviour)
    		delete this->behaviour ;*/
    behaviour = f ;
}

Form * Feature::getBehaviour( const Point & p) const
{
    return this->behaviour ;
}

Form * Feature::getBehaviour() const
{
    return this->behaviour ;
}

Feature * Feature::getChild(size_t i) const
{
    return m_c[i] ;
}

const std::vector<Feature *> & Feature::getChildren() const
{
    return m_c ;
}

std::vector<Feature *> & Feature::getChildren()
{
    return m_c ;
}

bool Feature::inMask(const Point &p, double d) const
{
    if(mask.size() == 0)
       return true ;

    bool isin = false ;
    for(size_t i = 0 ; !isin && i < mask.size() ; i++)
        isin = mask[i]->inBoundary( p, d ) ;

    return isin ;
}

bool Feature::onMaskBoundary(const Point &p, double d) const
{
    if(mask.size() == 0)
       return false ;

    for(size_t i = 0 ; i < mask.size() ; i++)
    {
        if(mask[i]->onBoundary( p, d ))
            return true ;
    }

    return false ;
}

bool Feature::inBoundary(const Point &p, double d) const
{

    bool inMaskBoundary = false ;
    for(size_t i = 0 ; i < mask.size() ; i++)
       inMaskBoundary |= mask[i]->inBoundary( p, d ) ;

    if(inMaskBoundary && in(p))
        return true ;

    Point proj(p) ;
    project(&proj) ;
    return ((squareDist3D(proj, p) < d*d) || in(p));
}

void Feature::setMask(std::vector<Feature *> & m) 
{
    mask.clear() ;
    mask = m ;
    if(m_f && !isMaskedBy(m_f))
        mask.push_back(m_f) ;
}

bool Feature::onBoundary(const Point &p, double d) const
{
    if(p == getCenter())
        return false ;

//    if(!inMask(p,d))
//        return false ;

    bool onMaskBoundary = false ;
    for(size_t i = 0 ; i < mask.size() ; i++)
       onMaskBoundary |= mask[i]->onBoundary( p, d ) ;

    if(onMaskBoundary && in(p))
        return true ;

    Point proj(p) ;
    project(&proj) ;
    return (dist(proj, p) < d) ;
}

void Feature::makeMaskBoundary() 
{
    if(mask.size() == 0)
        return ;

    if(m_f && !isMaskedBy(m_f))
        mask.push_back( m_f ) ;

    std::vector<Point> newPoints ;
    double density = -1 ;
    if(getBoundingPoints().size() > 1)
        density = dist( getBoundingPoint(0), getBoundingPoint(1) ) ;

    if(density < 0)
        return ;

    for(size_t i = 0 ; i < mask.size() ; i++)
    {
        std::vector<Point> inter = this->intersection( mask[i] ) ;
        for(size_t j = 0 ; j < inter.size() ; j++)
            newPoints.push_back( inter[j] ) ;

        for(size_t j = 0 ; j < getBoundingPoints().size() ; j++)
        {
            if( mask[i]->in( getBoundingPoint(j) ) )
            {
                Point p = getBoundingPoint(j) ;
                mask[i]->project( &p ) ;
                if( dist( p, getBoundingPoint(j) ) < density*0.25 )
                {
                    Point next = getBoundingPoint( j==getBoundingPoints().size()-1 ? 0 : j+1 ) ;
                    Point prev = getBoundingPoint( j==0 ? getBoundingPoints().size()-1 : j-1 ) ;
                    Point target = p ;
                    if(mask[i]->in( next ))
                        target = next ;
                    else if(mask[i]->in( prev ))
                        target = prev ;
                    p += target ;
                    p *= 0.5 ;
                    this->project(&p) ;
                    mask[i]->project( &target ) ;
                    newPoints.push_back( p ) ;
                    newPoints.push_back( target ) ;
                }
                else
                    newPoints.push_back( getBoundingPoint( j ) ) ;
            }
        }       

        for(size_t j = 0 ; j < mask[i]->getBoundingPoints().size() ; j++)
        {
            if( in( mask[i]->getBoundingPoint(j) ) )
                newPoints.push_back( mask[i]->getBoundingPoint( j ) ) ;
        }       
 
    }

    PointArray next( newPoints.size() ) ;
    for(size_t i = 0 ; i < newPoints.size() ; i++)
        next[i] = new Point(newPoints[i].getX(), newPoints[i].getY(), newPoints[i].getZ(), newPoints[i].getT()) ;
    setBoundingPoints( next ) ;

    newPoints.clear() ;
    for(size_t j = 0 ; j < getInPoints().size() ; j++)
    {
        for(size_t i = 0 ; i < mask.size() ; i++)
        {
            if( mask[i]->in( getInPoint(j) ) )
            {
                Point p = getInPoint(j) ;
                mask[i]->project( &p ) ;
                if( dist( p, getInPoint(j) ) > density*0.1 )
                {
                     newPoints.push_back( getInPoint(j) ) ;
                     break ;
                }       
            }
        }
    }

    PointArray nextIn( newPoints.size() ) ;
    for(size_t i = 0 ; i < newPoints.size() ; i++)
        nextIn[i] = new Point(newPoints[i].getX(), newPoints[i].getY(), newPoints[i].getZ(), newPoints[i].getT()) ;
    setInPoints( nextIn ) ;

}


std::vector<Feature *> Feature::getDescendants() const
{
    std::vector<Feature *> childrenToCheck = getChildren() ;
    std::vector<Feature *> ret ;

    while(!childrenToCheck.empty())
    {
        std::vector<Feature *> newChildren ;

        for(size_t i = 0 ; i< childrenToCheck.size() ; i++)
        {
            ret.push_back(childrenToCheck[i]) ;
            newChildren.insert(newChildren.end(), childrenToCheck[i]->getChildren().begin(), childrenToCheck[i]->getChildren().end()) ;
        }

        childrenToCheck = newChildren ;
    }

    return ret ;
}

void  Feature::setFather(Feature * const f)
{
    m_f = f ;
}


CompositeFeature::~CompositeFeature()
{
    for(size_t i = 0 ; i < components.size() ; i++)
        delete components[i] ;
}

std::vector<VirtualFeature *> & CompositeFeature::getComponents()
{
    return components ;
}

const std::vector<VirtualFeature *> & CompositeFeature::getComponents() const
{
    return components ;
}

const Feature * Feature::getFather() const
{
    return m_f ;
}

Feature * Feature::getFather()
{
    return m_f ;
}



Feature::~Feature()
{
// 	for(size_t i = 0 ; i < this->boundary->size() ; i++)
// 		delete this->boundary->getPoint(i) ;

    delete this->behaviour ;
}

std::vector<DelaunayTriangle *> Feature::getBoundingElements2D( FeatureTree * dt) const
{
    std::vector<DelaunayTriangle *> tri = dt->get2DMesh()->getConflictingElements(dynamic_cast<const Geometry *>(this)) ;

    std::vector<DelaunayTriangle *> ret  ;

    for(size_t i = 0 ; i < tri.size() ; i++)
    {
        if(tri[i]->intersects(dynamic_cast<const Geometry *>(this)))
            ret.push_back(tri[i]) ;
    }

    return ret ;

}

std::vector<DelaunayTetrahedron *> Feature::getBoundingElements3D( FeatureTree * dt) const
{
    std::vector<DelaunayTetrahedron *> tri = dt->get3DMesh()->getConflictingElements(dynamic_cast<const Geometry *>(this)) ;

    std::vector<DelaunayTetrahedron *> ret  ;

    for(size_t i = 0 ; i < tri.size() ; i++)
    {
        if(tri[i]->Tetrahedron::intersects(dynamic_cast<const Geometry *>(this)))
            ret.push_back(tri[i]) ;
    }

    return ret ;

}

