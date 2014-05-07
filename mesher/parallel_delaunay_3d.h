
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

namespace Mu
{

	class ParallelDelaunayTree3D :public Mesh<DelaunayTetrahedron, DelaunayTreeItem3D>
	{
	protected:
		std::vector<Geometry *> domains ;
		std::vector<DelaunayTree3D *> meshes ;
	public:
		virtual std::vector<DelaunayTreeItem3D *> & getTree() = 0;
		virtual const std::vector<DelaunayTreeItem3D *> & getTree() const = 0 ;
		virtual std::vector<Point * > & getAdditionalPoints() = 0 ;
		virtual const std::vector<Point * > & getAdditionalPoints() const = 0 ;
		virtual void extrude(double dt) = 0 ;
		virtual void extrude(const Vector & dt) = 0 ;
		virtual double getInternalScale() const { return 1. ;} ;
	public:
		ParallelDelaunayTree3D(Point * p0,  Point *p1,  Point *p2, Point *p3, const std::vector<Geometry *> & domains) ;
		virtual ~ParallelDelaunayTree3D() {} ;
		virtual std::vector<DelaunayTetrahedron *> getElements() = 0;
		virtual std::vector<DelaunayTetrahedron *> getConflictingElements(const Point  * p)  = 0;
		virtual std::vector<DelaunayTetrahedron *> getConflictingElements(const Geometry * g) = 0;

		virtual void setElementOrder(Order o, double dt = 0.) = 0;
		virtual void insert(Point *) ;

		virtual const size_t & getLastNodeId() const = 0;
	} ;
} ;



#endif  //_PARA_DELAUNAY_3D_H_
