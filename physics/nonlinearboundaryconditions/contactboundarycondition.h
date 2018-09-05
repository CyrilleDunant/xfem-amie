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

#ifndef CONTACT_CONDITION_H
#define CONTACT_CONDITION_H

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
        std::map<const Point *, double> previousStiffnesses ;
        std::map<std::pair<const Point *, const Point *>, DelaunayTriangle * > affectedElements ;
        
        std::map<const Point *, double> errors ;
        std::map<const Point *, double> previousErrors ;
        std::map<const Point *, bool> positions ;
        
        std::vector<DelaunayTriangle *> edgeElements ;
        
        bool conv = true ;
        double scale = 1. ;
        double currentError = 0 ;
        bool active = true ;
        double threshold = 1e-6 ;
        
    public: 
        ContactBoundaryCondition(Geometry * geo) ;
        
        void initialise(Mesh<DelaunayTriangle,DelaunayTreeItem> * mesh) ;
        void reInitialise() ;
        void update() ;
        void applyBoundaryConditions(Assembly * a, Mesh<DelaunayTriangle,DelaunayTreeItem> * mesh)  ;
        
        void setScale(double s) ;
        
        void setActive(bool act) ;
        
        
        bool verifyConvergence() const ;
        bool converged() const ;
        double error() const ;
        
        void print() ;
    } ;
}

#endif
