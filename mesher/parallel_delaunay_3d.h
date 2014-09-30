
//
// C++ Interface: delaunay
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//


#ifndef _PARA_DELAUNAY_3D_H_
#define  _PARA_DELAUNAY_3D_H_
#include "delaunay_3d.h"

namespace Amie
{
    
class ParallelDelaunayTree3D :public Mesh<DelaunayTetrahedron, DelaunayTreeItem3D>
{
    int getDomain(const DelaunayTetrahedron * tet) const;
    int getMesh(const DelaunayTreeItem3D * self) const;
    bool inDomain(int domain_index, const DelaunayTetrahedron * tet) const;
    bool isSame(const DelaunayTreeItem3D * i0, const DelaunayTreeItem3D * i1) const;
    int getDomain(const Point & center) const;
protected:
    std::vector<const Geometry *> domains ;
    std::vector<DelaunayTree3D *> meshes ;
    std::vector<Point *> additionalPoints ;
    std::vector<DelaunayTreeItem3D *> tree ;
    void addSharedNodes(size_t nodes_per_side, size_t time_planes = 2, double timestep = 2) ;
   
//     std::vector<double> maxRadius ;
    int global_counter ;
    std::vector< std::vector<DelaunayTetrahedron *> > caches ;
    std::vector<Vector> coefs ;
public:
    virtual std::vector<DelaunayTreeItem3D *> & getTree() ;
    virtual const std::vector<DelaunayTreeItem3D *> & getTree() const;
    virtual std::vector<Point * > & getAdditionalPoints()  ;
    virtual const std::vector<Point * > & getAdditionalPoints() const  ;
    virtual void extrude(double dt);
    virtual void extrude(const Vector & dt) ;
    virtual double getInternalScale() const {
        return meshes[0]->getInternalScale() ;
    } ;
public:
    ParallelDelaunayTree3D(Point * p0,  Point *p1,  Point *p2,  Point *p3, const std::vector<const Geometry *> & domains) ;
    virtual ~ParallelDelaunayTree3D() {} ;
    virtual std::vector<DelaunayTetrahedron *> getElements() ;
    virtual std::vector<DelaunayTetrahedron *> getConflictingElements(const Point  * p) ;
    virtual std::vector<DelaunayTetrahedron *> getConflictingElements(const Geometry * g) ;

    virtual void setElementOrder(Order o, double dt = 0.) ;
    virtual void insert(Point *) ;

    virtual size_t getLastNodeId() const {
        return global_counter ;
    };

    virtual int addToTree(DelaunayTreeItem3D * toAdd)
    {
        return 0 ;
    }
    
    virtual std::vector<DelaunayTetrahedron *> getNeighbourhood(DelaunayTetrahedron * element) const ;

    virtual DelaunayTreeItem3D * getInTree(int index) const
    {
        return nullptr ;
    }
    
    virtual size_t size() const
    {
        size_t s = 0 ;
        for( const auto & m : meshes)
            s += m->size() ;
        return s ;
    }
    
    virtual unsigned int generateCache(const std::vector<DelaunayTetrahedron *> original) 
    { 
        caches.push_back(original);
        coefs.push_back(Vector());
        return caches.size()-1 ;
    } ;
    virtual std::vector<DelaunayTetrahedron *> getCache(unsigned int cacheID) 
    {
        return caches[cacheID] ; 
    } ;
    virtual unsigned int generateCache(const Geometry * locus, const Geometry * source = nullptr, Function smoothing = Function("1"))
    {
        VirtualMachine vm ;
        std::vector<double> co ;
        std::vector<DelaunayTetrahedron *> elems= getConflictingElements(locus) ;
        for(auto & element : elems)
        {
            if(source && element->getBehaviour()->getSource() != source)
                continue ;
            
            if(locus->in(element->getCenter()))
            {
                       
                Function x = element->getXTransform() ;
                Function y = element->getYTransform() ;
                Function z = element->getZTransform() ;
                Function t = element->getTTransform() ;
                for(size_t i = 0 ; i < co.size() ; i++)
                {
                    double xx= vm.eval(x, element->getGaussPoints().gaussPoints[i].first) ;
                    double xy = vm.eval(y, element->getGaussPoints().gaussPoints[i].first) ;
                    double xz = vm.eval(z, element->getGaussPoints().gaussPoints[i].first) ;
                    double xt = vm.eval(t, element->getGaussPoints().gaussPoints[i].first) ;

                    co.push_back(vm.eval(smoothing, xx, xy, xz, xt));
                }
            }
        }
        
        Vector cf(co.size()) ;
        std::copy(co.begin(), co.end(), &cf[0]) ;
        coefs.push_back(cf);
        caches.push_back(elems);
        return coefs.size()-1 ;
    } ;
} ;
} ;



#endif  //_PARA_DELAUNAY_3D_H_
