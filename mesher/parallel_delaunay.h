
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


#ifndef _PARA_DELAUNAY_2D_H_
#define  _PARA_DELAUNAY_2D_H_
#include "delaunay.h"

namespace Amie
{

class ParallelDelaunayTree :public Mesh<DelaunayTriangle, DelaunayTreeItem>
{
    int getDomain(const DelaunayTriangle * tet) const ;
    int getMesh(const DelaunayTreeItem * self) const;
    bool inDomain(int domain_index, const DelaunayTriangle * tet) const;
    bool isSame(const DelaunayTreeItem * i0, const DelaunayTreeItem * i1) const;
    int getDomain(const Point & center) const;
protected:
    std::vector< std::vector<int> > elementMap ; //the negative ids indicate elements not valid for the mesh
    std::vector<const Geometry *> domains ;
    std::vector<DelaunayTree *> meshes ;
    std::vector<Point *> additionalPoints ;
//     std::vector<DelaunayTreeItem *> tree ;
    int global_counter ;
public:
    virtual std::vector<DelaunayTreeItem *> & getTree() ;
    virtual const std::vector<DelaunayTreeItem *> & getTree() const;
    virtual std::vector<Point * > & getAdditionalPoints()  ;
    virtual const std::vector<Point * > & getAdditionalPoints() const  ;
    void addSharedNodes( size_t nodes_per_side, size_t time_planes, double timestep) ;
    virtual void extrude(double dt);
    virtual void extrude(const Vector & dt) ;
    virtual double getInternalScale() const {
        return meshes[0]->getInternalScale() ;
    } ;
     virtual std::vector<DelaunayTriangle *> getElements() const ;
public:
    ParallelDelaunayTree(Point * p0,  Point *p1,  Point *p2, const std::vector<const Geometry *> & domains) ;
    virtual ~ParallelDelaunayTree() {} ;
    virtual std::vector<DelaunayTriangle *> getConflictingElements(const Point  * p) ;
    virtual std::vector<DelaunayTriangle *> getConflictingElements(const Geometry * g) ;

    virtual void setElementOrder(Order o, double dt = 0.) ;
    virtual void insert(Point *) ;

    virtual size_t getLastNodeId() const {
        return global_counter ;
    };

    virtual std::vector<DelaunayTriangle *> getNeighbourhood(DelaunayTriangle * element) const ;
    
    virtual int addToTree(DelaunayTreeItem * toAdd)
    {
        return 0 ;
    }

    virtual DelaunayTreeItem * getInTree(int index) const 
    {
    }
    
    virtual size_t size() const
    {
        size_t s = 0 ;
        for( const auto & m : meshes)
            s += m->size() ;
        return s ;
    }
    virtual void clearCaches()
    {
        caches.clear();
        coefs.clear();
        elementMap.clear();
        allElementsCacheID = -1 ;
    }
    
    virtual unsigned int generateCache(const Geometry * locus, const Geometry * source = nullptr, Function smoothing = Function("1")) ;
    
    virtual unsigned int generateCache () ;
    
    Vector getField( FieldType f, unsigned int cacheID, int dummy = 0, double t = 0) const;

    Vector getField( FieldType f, int dummy = 0, double t = 0) const;
    
    Vector getSmoothedField (  FieldType f0, unsigned int cacheID, IntegrableEntity * e,int dummy = 0, double t = 0 ) const ;
    
    std::pair<Vector, Vector> getSmoothedFields ( FieldType f0, FieldType f1, unsigned int cacheID, IntegrableEntity * e ,int dummy = 0, double t = 0 ) const ;  
      

    class iterator
    {
    private:
        ParallelDelaunayTree * msh ;
        size_t cacheID ;
        size_t position ;

    public:
        iterator( ParallelDelaunayTree * msh, size_t cacheID, size_t position) : msh(msh), cacheID(cacheID), position(position) { } ;
        iterator( ParallelDelaunayTree * msh, size_t position) : msh(msh), position(position) 
        { 
            if(msh->allElementsCacheID == -1)
                msh->allElementsCacheID = msh->generateCache() ;
            cacheID = msh->allElementsCacheID ;
        } ;
        
        bool operator ==(const iterator & i) const
        {
            return i.position == position ;
        }
        bool operator !=(const iterator & i) const
        {
            return i.position != position ;
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
        
        DelaunayTriangle * operator-> ( ) { 
            return dynamic_cast<DelaunayTriangle *>(msh->meshes[msh->elementMap[cacheID][position]]->getInTree(msh->caches[cacheID][position])) ;
        }
        operator DelaunayTriangle * (){ 
            return static_cast< DelaunayTriangle *>(msh->getInTree(msh->caches[cacheID][position])) ;
        } 
        
        size_t size() const
        {
            return msh->caches[cacheID].size() ;
        }
        size_t getPosition() const
        {
            return position ;
        }
        size_t getId() const {return cacheID ; }
        
    } ;
    
    class const_iterator
    {
    private:
        const ParallelDelaunayTree * msh ;
        size_t cacheID ;
        size_t position ;
    public:
        const_iterator( const ParallelDelaunayTree * msh, size_t cacheID, size_t position) : msh(msh), cacheID(cacheID), position(position) { } ;
        
        const_iterator( const ParallelDelaunayTree * msh, size_t position) : msh(msh), position(position) 
        { 
            cacheID = msh->allElementsCacheID ;
        } ;
        
        bool operator ==(const const_iterator & i) const
        {
            return i.position == position ;
        }
        bool operator !=(const const_iterator & i) const
        {
            return i.position != position ;
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
        
        const DelaunayTriangle * operator-> ( ) const { 
              return dynamic_cast<const DelaunayTriangle *>(msh->meshes[msh->elementMap[cacheID][position]]->getInTree(msh->caches[cacheID][position])) ;
        }
        
        operator const DelaunayTriangle * () const { 
            return static_cast<const DelaunayTriangle *>(msh->getInTree(msh->caches[cacheID][position])) ;
        } 
        size_t size() const
        {
            return msh->caches[cacheID].size() ;
        }
        size_t getPosition() const
        {
            return position ;
        }
        size_t getId() const {return cacheID ; }
        
    } ;
    
    
    iterator begin()
    {
        return iterator(this, 0) ;
    }
    const_iterator cbegin() const
    {
        return const_iterator(this, 0) ;
    }
    iterator end()
    {
        if(allElementsCacheID == -1)
            allElementsCacheID = generateCache() ;
        return iterator(this,allElementsCacheID, caches[allElementsCacheID].size()) ;
    }
    const_iterator cend() const 
    {
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
