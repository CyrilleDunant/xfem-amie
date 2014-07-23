
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

namespace Mu
{

	class ParallelDelaunayTree :public Mesh<DelaunayTriangle, DelaunayTreeItem>
	{
	protected:
		std::vector< std::vector<int> > elementMap ; //the negative ids indicate elements not valid for the mesh
		std::vector<Geometry *> domains ;
		std::vector<DelaunayTree *> meshes ;
		std::vector<Point *> additionalPoints ;
		int global_counter ;
	public:
		virtual std::vector<DelaunayTreeItem *> & getTree() = 0;
		virtual const std::vector<DelaunayTreeItem *> & getTree() const = 0 ;
		virtual std::vector<Point * > & getAdditionalPoints() = 0 ;
		virtual const std::vector<Point * > & getAdditionalPoints() const = 0 ;
		virtual void extrude(double dt);
		virtual void extrude(const Vector & dt) ;
		virtual double getInternalScale() const { return 1. ;} ;
	public:
		ParallelDelaunayTree(Point * p0,  Point *p1,  Point *p2, const std::vector<Geometry *> & domains) ;
		virtual ~ParallelDelaunayTree() {} ;
		virtual std::vector<DelaunayTriangle *> getElements() ;
		virtual std::vector<DelaunayTriangle *> getConflictingElements(const Point  * p) ;
		virtual std::vector<DelaunayTriangle *> getConflictingElements(const Geometry * g) ;

		virtual void setElementOrder(Order o, double dt = 0.) = 0;
		virtual void insert(Point *) ;

		virtual const size_t & getLastNodeId() const = 0;
	} ;
} ;



#endif  //_PARA_DELAUNAY_3D_H_
