
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
    std::vector<std::vector<int>> elementMap ;
    void addSharedNodes(size_t nodes_per_side, size_t time_planes = 2, double timestep = 2) ;
   
//     std::vector<double> maxRadius ;
    int global_counter ;
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
    virtual std::vector<DelaunayTetrahedron *> getElements() const ;
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
    
    virtual unsigned int generateCache(const Geometry * locus, const Geometry * source = nullptr, Function smoothing = Function("1")) ;
    
    virtual unsigned int generateCache () ;
    
    Vector getField( FieldType f, unsigned int cacheID, int dummy = 0, double t = 0) const ;

    Vector getField( FieldType f, int dummy = 0, double t = 0) const ;
    
    Vector getSmoothedField (  FieldType f0, unsigned int cacheID, IntegrableEntity * e,int dummy = 0, double t = 0 ) const ;
    
    std::pair<Vector, Vector> getSmoothedFields ( FieldType f0, FieldType f1, unsigned int cacheID, IntegrableEntity * e ,int dummy = 0, double t = 0 ) const ;  
      
    class iterator
    {
    private:
        ParallelDelaunayTree3D * msh ;
        size_t cacheID ;
        size_t position ;

    public:
        iterator( ParallelDelaunayTree3D * msh, size_t cacheID, size_t position ) : msh(msh), cacheID(cacheID), position(position) { } ;
        iterator( ParallelDelaunayTree3D * msh, size_t position) : msh(msh), position(position) 
        { 
            if(msh->allElementsCacheID == -1)
                msh->allElementsCacheID = msh->generateCache() ;
            cacheID = msh->allElementsCacheID ;
        } ;
                bool operator ==(const iterator & i) const
        {
            return i.position == position ;
        }
        
        bool operator <=(const iterator & i)const
        {
            return position <= i.position ;
        }
        bool operator <(const iterator & i) const
        {
            return position < i.position ;
        }
        bool operator >=(const iterator & i)const
        {
            return position >= i.position ;
        }
        bool operator >(const iterator & i) const
        {
            return position > i.position ;
        }
        
        
        iterator& operator++() {
            position++ ;
            // actual increment takes place here
            return *this;
        }
        
        iterator operator++(int) {
            iterator tmp(*this); // copy
            operator++(); // pre-increment
            return tmp;   // return old value
        }
        
        iterator& operator+=(int i) {
            position +=i ;
            return *this;
        }
        
        friend iterator operator+(iterator lhs,  int i) 
        {
            return lhs += i; 
        }
        
        iterator& operator--() {
            position-- ;
            return *this;
        }
        
        iterator operator--(int) {
            iterator tmp(*this); // copy
            operator--(); // pre-increment
            return tmp;   // return old value
        }
        
        iterator& operator-=(int i) {
            position -=i ;
            return *this;
        }
        
        friend iterator operator-(iterator lhs,  int i) 
        {
            return lhs -= i; 
        }
        
        DelaunayTetrahedron * operator-> ( ) {
            int meshId = msh->elementMap[cacheID][position] ;
            int elemID = msh->caches[cacheID][position] ;
            return static_cast<DelaunayTetrahedron *>(msh->meshes[meshId]->getInTree(elemID)) ;
        }
        
    } ;
    
    class const_iterator
    {
    private:
        ParallelDelaunayTree3D * msh ;
        size_t cacheID ;
        size_t position ;
    public:
        const_iterator( ParallelDelaunayTree3D * msh, size_t cacheID, size_t position) : msh(msh), cacheID(cacheID), position(position) { } ;
        
        const_iterator( ParallelDelaunayTree3D * msh, size_t position) : msh(msh), position(position) 
        { 
            if(msh->allElementsCacheID == -1)
                msh->allElementsCacheID = msh->generateCache() ;
            cacheID = msh->allElementsCacheID ;
        } ;
        
        bool operator ==(const const_iterator & i) const
        {
            return i.position == position ;
        }
        bool operator <=(const const_iterator & i)const
        {
            return position <= i.position ;
        }
        bool operator <(const const_iterator & i) const
        {
            return position < i.position ;
        }
        bool operator >=(const const_iterator & i)const
        {
            return position >= i.position ;
        }
        bool operator >(const const_iterator & i) const
        {
            return position > i.position ;
        }
        
        const_iterator& operator++() {
            position++ ;
            // actual increment takes place here
            return *this;
        }
        
        const_iterator operator++(int) {
            const_iterator tmp(*this); // copy
            operator++(); // pre-increment
            return tmp;   // return old value
        }
        
        const_iterator& operator+=(int i) {
            position +=i ;
            return *this;
        }
        
        friend const_iterator operator+(const_iterator lhs,  int i) 
        {
            return lhs += i; 
        }
        
        const_iterator& operator--() {
            position-- ;
            return *this;
        }
        
        const_iterator operator--(int) {
            const_iterator tmp(*this); // copy
            operator--(); // pre-increment
            return tmp;   // return old value
        }
        
        const_iterator& operator-=(int i) {
            position -=i ;
            return *this;
        }
        
        friend const_iterator operator-(const_iterator lhs,  int i) 
        {
            return lhs -= i; 
        }
        
        const DelaunayTetrahedron * operator-> ( ) const { 
            return static_cast<const DelaunayTetrahedron *> (msh->meshes[msh->elementMap[cacheID][position]]->getInTree(msh->caches[cacheID][position])) ;
        }
    } ;
    
    
    iterator begin()
    {
        return iterator(this, 0) ;
    }
    const_iterator cbegin()
    {
        return const_iterator(this, 0) ;
    }
    iterator end()
    {
        if(allElementsCacheID == -1)
            allElementsCacheID = generateCache() ;
        return iterator(this,allElementsCacheID, caches[allElementsCacheID].size()) ;
    }
    const_iterator cend()
    {
        if(allElementsCacheID == -1)
            allElementsCacheID = generateCache() ;
        return const_iterator(this,allElementsCacheID, caches[allElementsCacheID].size()) ;
    }
    
    iterator begin( size_t cacheID)
    {
        return iterator(this, cacheID, 0) ;
    }
    const_iterator cbegin(size_t cacheID)
    {
        return const_iterator(this, cacheID, 0) ;
    }
    iterator end(size_t cacheID)
    {
        return iterator(this,cacheID, caches[cacheID].size()) ;
    }
    const_iterator cend(size_t cacheID)
    {
        return const_iterator(this,cacheID, caches[cacheID].size()) ;
    }

    
} ;
} ;



#endif  //_PARA_DELAUNAY_3D_H_
