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
        std::map<const Point *, double> referenceX ;
        std::map<const Point *, double> referenceY ;
        std::map<const Point *, Point> normalVectors ;
        std::map<const Point *, double> stiffnessesX ;
        std::map<const Point *, double> stiffnessesY ;
        std::map<const Point *, double> previousStiffnessesX ;
        std::map<const Point *, double> previousStiffnessesY ;
        std::map<std::pair<const Point *, const Point *>, DelaunayTriangle * > affectedElements ;
        
        std::map<const Point *, double> errorsX ;
        std::map<const Point *, double> previousErrorsX ;
        std::map<const Point *, double> errorsY ;
        std::map<const Point *, double> previousErrorsY ;
        std::map<const Point *, bool> positions ;
        
        std::vector<DelaunayTriangle *> edgeElements ;
        
        bool conv = true ;
        double scale = 1. ;
        double currentError = 0 ;
        bool active = true ;
        double threshold = 1e-5 ;
        int counter = 0 ;
        
    public: 
        ContactBoundaryCondition(Geometry * geo) ;
        
        void initialise(Mesh<DelaunayTriangle,DelaunayTreeItem> * mesh) ;
        void reInitialise() ;
        void update() ;
        void applyBoundaryConditions(Assembly * a, Mesh<DelaunayTriangle,DelaunayTreeItem> * mesh)  ;
        
        void setScale(double s) ;
        
        void setActive(bool act) ;
        virtual void postProcess() ;
        
        bool verifyConvergence() const ;
        bool converged() const ;
        double error() const ;
        
        void print() const;
    } ;
}

#endif
