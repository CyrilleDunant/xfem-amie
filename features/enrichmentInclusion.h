
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef __ENR_INCLUSION_H__
#define __ENR_INCLUSION_H__

#include "features.h"

namespace Amie
{

/** \brief Enrich the mesh to simulate a circular inclusion.
 *
 * This enrichment feature will add a soft discontinuity 
 * enrichment to the elements it crosses. If the inclusion
 * is smaller than, or fully contained by an element,
 * The enrichment will still be performed, yielding a 
 * composite field within the element.
 * 
 * This Inclusion will set no behaviour in the elements
 * thus derived classes should be used for actual simulation
 * see for example expansiveZone.h
 */
class EnrichmentInclusion :  public EnrichmentFeature,  public Circle
{
protected:
	bool updated ;
	std::vector<DelaunayTriangle *> cache ;
	std::map<Mesh<DelaunayTriangle, DelaunayTreeItem> *,std::set<size_t> > freeIds ;
    Function getBlendingFunction(const std::map<const Point *, int> & dofIds, const DelaunayTriangle * t)
    {
    //  return Function("1") ;

    // if(t->getOrder() == QUADRATIC)
    // {
    //  TriElement father(QUADRATIC) ;
    //  if(dofIds.find(t->first) != dofIds.end() && dofIds.find(t->second) == dofIds.end() && dofIds.find(t->third) == dofIds.end())
    //  {
    //      return father.getShapeFunction(0) + 0.25*father.getShapeFunction(1)+ 0.25*father.getShapeFunction(5);
    //  }
    //  
    //  if(dofIds.find(t->first) == dofIds.end() && dofIds.find(t->second) != dofIds.end() && dofIds.find(t->third) == dofIds.end())
    //  {
    //      return father.getShapeFunction(2) + 0.25*father.getShapeFunction(1)+ 0.25*father.getShapeFunction(3);
    //  }
    //  
    //  if(dofIds.find(t->first) == dofIds.end() && dofIds.find(t->second) == dofIds.end() && dofIds.find(t->third) != dofIds.end())
    //  {
    //      return father.getShapeFunction(4) + 0.25*father.getShapeFunction(3)+ 0.25*father.getShapeFunction(5);
    //  }
    //  
    //  if(dofIds.find(t->first) == dofIds.end() && dofIds.find(t->second) != dofIds.end() && dofIds.find(t->third) != dofIds.end())
    //  {
    //      return father.getShapeFunction(2)+father.getShapeFunction(3)+father.getShapeFunction(4) + 0.25*father.getShapeFunction(1)+ 0.25*father.getShapeFunction(5);
    //  }
    //  
    //  if(dofIds.find(t->first) != dofIds.end() && dofIds.find(t->second) == dofIds.end() && dofIds.find(t->third) != dofIds.end())
    //  {
    //      return father.getShapeFunction(0) + father.getShapeFunction(5) + father.getShapeFunction(4) + 0.25*father.getShapeFunction(1) +0.25*father.getShapeFunction(3);
    //  }
    //  
    //  if(dofIds.find(t->first) != dofIds.end() && dofIds.find(t->second) != dofIds.end() && dofIds.find(t->third) == dofIds.end())
    //  {
    //      return father.getShapeFunction(1)+father.getShapeFunction(0)+father.getShapeFunction(2) + 0.25*father.getShapeFunction(3) + 0.25*father.getShapeFunction(5);
    //  }
    // }


        Amie::TriElement father(Amie::LINEAR) ;
    //  Function f ;
    //  for(size_t i = 0 ; i < t->getBoundingPoints().size() ; i++)
    //  {
    //      if(dofIds.find(&(t->getBoundingPoint(i))) != dofIds.end())
    //          f += father.getShapeFunction(i) ;
    //  }
    //  return f ;
        
        if(dofIds.find(t->first) != dofIds.end() && dofIds.find(t->second) == dofIds.end() && dofIds.find(t->third) == dofIds.end())
        {
            return father.getShapeFunction(0) ;
        }
        
        if(dofIds.find(t->first) == dofIds.end() && dofIds.find(t->second) != dofIds.end() && dofIds.find(t->third) == dofIds.end())
        {
            return father.getShapeFunction(1) ;
        }
        
        if(dofIds.find(t->first) == dofIds.end() && dofIds.find(t->second) == dofIds.end() && dofIds.find(t->third) != dofIds.end())
        {
            return father.getShapeFunction(2) ;
        }
        
        if(dofIds.find(t->first) == dofIds.end() && dofIds.find(t->second) != dofIds.end() && dofIds.find(t->third) != dofIds.end())
        {
            return 1-father.getShapeFunction(0) ;
        }
        
        if(dofIds.find(t->first) != dofIds.end() && dofIds.find(t->second) == dofIds.end() && dofIds.find(t->third) != dofIds.end())
        {
            return 1-father.getShapeFunction(1) ;
        }
        
        if(dofIds.find(t->first) != dofIds.end() && dofIds.find(t->second) != dofIds.end() && dofIds.find(t->third) == dofIds.end())
        {
            return 1-father.getShapeFunction(2) ;
        }
        
        return Amie::Function("1") ;
    } ;

public:

/** \brief Constructor. Construct the inclusion from a supporting feature, a radius and center coordinates */
	EnrichmentInclusion(Feature *father, double radius, double x, double y) ;

/** \brief Constructor. Construct the inclusion from a radius and center coordinates */
	EnrichmentInclusion(double radius, double x, double y) ;
	virtual ~EnrichmentInclusion() ;
	
/** \brief return true if the argument is cut by the circle*/
	virtual bool enrichmentTarget(DelaunayTriangle * t) ;

/** \brief Enrich elements cut by the feature*/
	virtual void enrich(size_t & , Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree) ;
	
/** \brief return false*/
	virtual bool interacts(Feature * f, double d) const ;

/** \brief return false*/
	virtual bool inBoundary(const Point v) const ;

/** \brief return false*/
	virtual bool inBoundary(const Point *v) const ;
	
/** \brief return the list of elements cut by the feature*/
	virtual std::vector<DelaunayTriangle *> getElements2D( FeatureTree * dt)  ;

/** \brief return empty vector*/
	virtual std::vector<DelaunayTetrahedron *> getElements3D(FeatureTree * dt) {return std::vector<DelaunayTetrahedron *>(0) ;} 
	
/** \brief return list of elements cut by the feature*/
	virtual std::vector<DelaunayTriangle *> getIntersectingTriangles( FeatureTree * dt) ;
	
/** \brief do nothing*/
	virtual void setInfluenceRadius(double r) ;
	
/** \brief return empty vector*/
	virtual std::vector<Point *> getSamplingPoints() const { return std::vector<Point *>(0) ; }
	
/** \brief do nothing, return null*/
	virtual Point * pointAfter(size_t i) { return nullptr ; }
	
/** \brief return empty vector*/
	virtual std::vector<Geometry *> getRefinementZones( size_t level) const ;
	
	virtual void print() const
	{
		std::cout << "I am an enriched inclusion" << std::endl ;
	}
	
/** \brief return false*/
	virtual bool isVoid( const Point &) const {return false ;}

/** \brief set the radius*/
	virtual void setRadius(double newR) ;
	
// 	virtual void setSingularityHints(const Point & i, const Point & s, std::vector<Point> * hints) const ;
	
public:
	GEO_DERIVED_OBJECT(Circle) ;
	
	
	/** \brief Sample the surface. Does nothing.
	 * 
	 * This is necessary as we need to implement the interface. Of course, for a segmented line, it makes no sense.
	 * 
	 * @param n 
	 */
	virtual void sample(size_t n)
	{
		
	}
	
/** \brief do nothing*/
	virtual void step(double dt, std::valarray<double> *, const  Mesh <DelaunayTriangle, DelaunayTreeItem > * dtree);
	
/** \brief return true if the radius has changed*/
	virtual bool moved() const ;

/** \brief compute and cache the elements to enrich*/
	void update(Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree) ;

protected:
	
	bool changed ;
} ;


} ;



#endif
