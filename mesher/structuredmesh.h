#ifndef STRUCTURED_MESH_H
#define STRUCTURED_MESH_H

#include "mesh.h"
#include "delaunay.h"
#include "../utilities/grid.h"

namespace Mu
{
	class StructuredMesh : public Mesh<DelaunayTriangle, DelaunayTreeItem>
	{
	protected:
		std::vector< Point *> points ;
		Grid grid ;
		void addSharedNodes(size_t nodes_per_side, size_t time_planes, double timestep) ;
		size_t global_counter ;
		std::vector <DelaunayTreeItem *> tree ;
	public:
		StructuredMesh(double sizeX, double sizeY, int div, const Point & center ) ;
		virtual ~StructuredMesh() ;
		virtual std::vector<DelaunayTriangle *> getElements();
		virtual std::vector<DelaunayTriangle *> getConflictingElements(const Point  * p) const;
		virtual std::vector<DelaunayTriangle *> getConflictingElements(const Geometry * g) const;
		virtual std::vector<Point * > & getAdditionalPoints(){ return points ;} ;
		virtual const std::vector<Point * > & getAdditionalPoints() const {return points ;} ;
		virtual void setElementOrder(Order o);
		virtual void insert(Point *) ;
		virtual size_t & getLastNodeId();
		virtual const size_t & getLastNodeId() const;
		virtual std::vector<DelaunayTreeItem *> & getTree() { return tree ;}
		virtual const std::vector<DelaunayTreeItem *> & getTree() const { return tree ;}
		
	} ;
}

#endif // STRUCTURED_MESH_H
