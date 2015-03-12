// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011

#ifndef STRUCTURED_MESH_H
#define STRUCTURED_MESH_H

#include "mesh.h"
#include "delaunay.h"
#include "delaunay_3d.h"
#include "../filters/voxelfilter.h"
#include "../utilities/grid.h"

namespace Amie
{
class voxelfilter ;

class StructuredMesh : public Mesh<DelaunayTriangle, DelaunayTreeItem>
{
protected:
    std::vector< Point *> points ;
    Grid grid ;
    void addSharedNodes ( size_t nodes_per_side, size_t time_planes, double timestep ) ;
    size_t global_counter ;
    std::vector <DelaunayTriangle *> tree ;
    virtual std::vector<DelaunayTriangle *> getElements() ;

public:
    StructuredMesh ( double sizeX, double sizeY, int div, const Point & center ) ;
    virtual size_t size() const {
        return tree.size() ;
    } ;
    virtual ~StructuredMesh() ;

    virtual std::vector<DelaunayTriangle *> getConflictingElements ( const Point  * p ) ;
    virtual std::vector<DelaunayTriangle *> getConflictingElements ( const Geometry * g ) ;
    virtual std::vector<Point * > & getAdditionalPoints() {
        return points ;
    } ;
    virtual const std::vector<Point * > & getAdditionalPoints() const {
        return points ;
    } ;
    virtual void setElementOrder ( Order o, double dt = 0 );
    virtual void insert ( Point * ) ;
    virtual size_t getLastNodeId() const;


} ;

class MicDerivedMesh : public Mesh<DelaunayTetrahedron, DelaunayTreeItem3D>
{
protected:
    size_t global_counter ;
    std::vector<DelaunayTetrahedron *> tree ;
    std::vector<double> times ;
    std::vector< Point *> additionalPoints ;
    std::map<unsigned char,Form *> behaviourMap ;
    virtual std::vector<DelaunayTetrahedron *> getElements()  {return tree ; };
    bool single ;
    int currentMesh ;
    double currentTime ;
    std::string voxelSource ;
public:

    virtual int addToTree ( DelaunayTreeItem3D * toAdd ) ;
    
    virtual DelaunayTreeItem3D * getInTree ( int index ) const { return tree[index] ;} ;
    virtual std::vector<Point * > & getAdditionalPoints() { return additionalPoints ;} ;
    virtual const std::vector<Point * > & getAdditionalPoints() const { return additionalPoints ;} ;
    virtual void extrude ( double dt ) ;
    virtual void extrude ( const Vector & dt ) ;

    void addSharedNodes( size_t nodes_per_side, size_t time_planes, double timestep, const TetrahedralElement *father = nullptr) ;
    
public:
    MicDerivedMesh(const char * voxelSource, std::map<unsigned char,Form *> behaviourMap) ;
    MicDerivedMesh(const char * voxelSource, std::map<unsigned char,Form *> behaviourMap, std::vector<double> times) ;
    
    virtual ~MicDerivedMesh() {} ;

    virtual std::vector<DelaunayTetrahedron *> getConflictingElements ( const Point  * p );
    
    virtual std::vector<DelaunayTetrahedron *> getConflictingElements ( const Geometry * g ) ;
    
    virtual std::vector<DelaunayTetrahedron *> getNeighbourhood ( DelaunayTetrahedron * element ) const  ;

    virtual void setElementOrder ( Order elemOrder, double dt = 0. ) ;

    virtual void insert ( Point * ) { } ;

    virtual size_t getLastNodeId() const {return global_counter ;};
    virtual size_t size() const {return tree.size() ;} ;
    
    virtual bool step( double dt ) ;
} ;

}

#endif // STRUCTURED_MESH_H
