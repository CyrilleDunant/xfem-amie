#ifndef STRUCTURED_MESH_H
#define STRUCTURED_MESH_H

#include "mesh.h"
#include "delaunay.h"
#include "../utilities/grid.h"

namespace Mu
{
	class StructuredMesh : public Mesh<DelaunayTriangle>
	{
	protected:
		std::vector<DelaunayTriangle *> triangles ;
		std::vector<Point *> points ;
		Grid grid ;
		void addSharedNodes(size_t nodes_per_side, size_t time_planes, double timestep) ;
		int global_counter ;
	public:
		StructuredMesh(double sizeX, double sizeY, int div, const Point & center ) ;
		virtual ~StructuredMesh() ;
		virtual std::vector<DelaunayTriangle *> getElements();
		virtual std::vector<DelaunayTriangle *> getConflictingElements(const Point  * p);
		virtual std::vector<DelaunayTriangle *> getConflictingElements(const Geometry * g);
		virtual void setElementOrder(Order o);
		virtual void insert(Point *) ;
		
	} ;
}

#endif // STRUCTURED_MESH_H
