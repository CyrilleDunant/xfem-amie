
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
		std::vector< std::vector<int> > elementMap ; //the negative ids indicate elements not valid for the mesh
		std::vector<Geometry *> domains ;
		std::vector<DelaunayTree3D *> meshes ;
		std::vector<Point *> additionalPoints ;
		std::vector<DelaunayTreeItem3D *> tree ;
		int global_counter ;
	public:
		virtual std::vector<DelaunayTreeItem3D *> & getTree() ;
		virtual const std::vector<DelaunayTreeItem3D *> & getTree() const;
		virtual std::vector<Point * > & getAdditionalPoints()  ;
		virtual const std::vector<Point * > & getAdditionalPoints() const  ;
		virtual void extrude(double dt);
		virtual void extrude(const Vector & dt) ;
		virtual double getInternalScale() const { return meshes[0]->getInternalScale() ;} ;
	public:
		ParallelDelaunayTree3D(Point * p0,  Point *p1,  Point *p2,  Point *p3, const std::vector<Geometry *> & domains) ;
		virtual ~ParallelDelaunayTree3D() {} ;
		virtual std::vector<DelaunayTetrahedron *> getElements() ;
		virtual std::vector<DelaunayTetrahedron *> getConflictingElements(const Point  * p) ;
		virtual std::vector<DelaunayTetrahedron *> getConflictingElements(const Geometry * g) ;

		virtual void setElementOrder(Order o, double dt = 0.) ;
		virtual void insert(Point *) ;

		virtual size_t getLastNodeId() const {return global_counter ;};
		
		virtual size_t addToTree(DelaunayTreeItem3D * toAdd)
		{
			return 0 ;
		}
		
		virtual DelaunayTreeItem3D * getInTree(int index) 
		{
		}
	} ;
} ;



#endif  //_PARA_DELAUNAY_3D_H_
