// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011

#ifndef STRUCTURED_MESH_H
#define STRUCTURED_MESH_H

#include "mesh.h"
#include "delaunay.h"
#include "../utilities/grid.h"

namespace Amie
{
	class StructuredMesh : public Mesh<DelaunayTriangle, DelaunayTreeItem>
	{
	protected:
		std::vector< Point *> points ;
		Grid grid ;
		void addSharedNodes(size_t nodes_per_side, size_t time_planes, double timestep) ;
		size_t global_counter ;
		std::vector <DelaunayTriangle *> tree ;
        std::vector< std::vector<DelaunayTriangle *> > caches ;
	public:
		StructuredMesh(double sizeX, double sizeY, int div, const Point & center ) ;
        virtual size_t size() const { return tree.size() ; } ;
		virtual ~StructuredMesh() ;
		virtual std::vector<DelaunayTriangle *> getElements();
		virtual std::vector<DelaunayTriangle *> getConflictingElements(const Point  * p) ;
		virtual std::vector<DelaunayTriangle *> getConflictingElements(const Geometry * g) ;
		virtual std::vector<Point * > & getAdditionalPoints(){ return points ;} ;
		virtual const std::vector<Point * > & getAdditionalPoints() const {return points ;} ;
		virtual void setElementOrder(Order o, double dt = 0);
		virtual void insert(Point *) ;
		virtual size_t getLastNodeId() const;
		
        virtual unsigned int generateCache(const std::vector<DelaunayTriangle *> original) 
        { 
            caches.push_back(original);
            return caches.size()-1 ;
        } ;
        virtual std::vector<DelaunayTriangle *> getCache(unsigned int cacheID) 
        {
            return caches[cacheID] ; 
        } ;
	} ;
}

#endif // STRUCTURED_MESH_H
