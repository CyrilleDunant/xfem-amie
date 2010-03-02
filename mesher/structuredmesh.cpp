#include "structuredmesh.h"

using namespace Mu ;

StructuredMesh::StructuredMesh(double sizeX, double sizeY, int div, const Point & center )
{
}
StructuredMesh::~StructuredMesh() 
{
}
std::vector<DelaunayTriangle *> StructuredMesh::getElements()
{
	return triangles ;
}
std::vector<DelaunayTriangle *> StructuredMesh::getConflictingElements(const Point  * p) 
{
	std::vector<DelaunayTriangle *> ret ;
	return ret ;
}
std::vector<DelaunayTriangle *> StructuredMesh::getConflictingElements(const Geometry * g) 
{
	std::vector<DelaunayTriangle *> ret ;
	return ret ;
}
void StructuredMesh::setElementOrder(Order o)
{
}
void StructuredMesh::insert(Point *) 
{
}