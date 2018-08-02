//
// C++ Header: contact
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2018-
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef CONTACT_CONDITION
#define CONTACT_CONDITION

#include "../../geometry/geometry_base.h"
#include "../../mesher/mesh.h"
#include "../../mesher/delaunay.h"
#include "../../features/boundarycondition.h"


namespace Amie
{
    class  ContactBoundaryCondition
    {
        
    protected:
        Geometry * baseGeometry ;
        std::map<const Point *, double> contactPointsAndTributary ;
        std::map<const Point *, Point> referencePoints ;
        std::map<const Point *, Point> normalVectors ;
        std::map<const Point *, double> distances ;
        std::map<const Point *, double> stiffnesses ;
        std::map<const Point *, DelaunayTriangle * > affectedElements ;
        
        std::vector<DelaunayTriangle *> edgeElements ;
        
        bool conv ;
        double baselength = -1. ;
        
    public: 
        ContactBoundaryCondition(Mesh<DelaunayTriangle,DelaunayTreeItem> * mesh, Geometry * geo) ;
        
        void update(Mesh<DelaunayTriangle,DelaunayTreeItem> * mesh) ;
        std::vector<BoundaryCondition *> getBoundaryconditions(Mesh<DelaunayTriangle,DelaunayTreeItem> * mesh)  ;
        
        bool converged() const ;
    } ;
}

#endif
