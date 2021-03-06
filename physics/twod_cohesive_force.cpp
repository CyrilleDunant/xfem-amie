
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//

#include "twod_cohesive_force.h"
#include "../mesher/delaunay.h"
#include "../polynomial/vm_base.h"
#include "../features/boundarycondition.h"

namespace Amie {



TwoDCohesiveForces::TwoDCohesiveForces(IntegrableEntity *s, IntegrableEntity *t, const SegmentedLine * sl)
{
    this->time_d = false ;
    this->type = NON_LINEAR ;
    this->target = t ;
    this->source = s ;
    active = false ;
    startArea = s->area() ;
    for(size_t i = 0 ; i < sl->getBoundingPoints().size()-1 ; i++)
    {
        Segment test(sl->getBoundingPoint(i), sl->getBoundingPoint(i+1)) ;
        normals.push_back(test.normal()) ;
    }

}

Form * TwoDCohesiveForces::getCopy() const
{

    TwoDCohesiveForces * copy = new TwoDCohesiveForces(*this) ;

    return copy ;
}

TwoDCohesiveForces::~TwoDCohesiveForces() { }

void TwoDCohesiveForces::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix &ret, VirtualMachine * vm) const
{
}

bool TwoDCohesiveForces::hasInducedForces() const
{
    return true ;
}

bool TwoDCohesiveForces::hasInducedMatrix() const
{
    return false ;
}


void TwoDCohesiveForces::getForces( ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector & f)
{


    bool enrichedDof = !(p_i.getDofID() == -1) ;

    if(!enrichedDof)
        return ;


    Vector apparentStress(0., gp.gaussPoints.size()*(3+3*(source->spaceDimensions() == 3))) ;
    source->getState().getField( NON_ENRICHED_REAL_STRESS_FIELD, gp.gaussPoints, apparentStress, true) ;
    for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
    {
        if(enrichedDof)
        {
            Vector stress(3) ;
            for(size_t j = 0 ; j < 3 ; j++)
                stress[j] = apparentStress[i*3+j] ;
            std::vector<Variable> v ;
            v.push_back(XI);
            v.push_back(ETA);

            Matrix grad (VirtualMachine().geval(p_i, Jinv[i],v, gp.gaussPoints[i].first, true)) ;
            Vector force = (Vector)(grad*stress) ;

            double normalAmplitude = force[0]*normals[0].getX() + force[1]*normals[0].getY();
            double tangeantAmplitude = -force[0]*normals[0].getY() + force[1]*normals[0].getX();

            Vector normalForce(2) ;
            normalForce[0] = normals[0].getX()*normalAmplitude ;
            normalForce[1] = normals[0].getY()*normalAmplitude ;

            Vector tangeantForce(2) ;
            tangeantForce[0] = -normals[0].getY()*tangeantAmplitude ;
            tangeantForce[1] = normals[0].getX()*tangeantAmplitude ;

            f += normalForce*gp.gaussPoints[i].second ;
            f += tangeantForce*gp.gaussPoints[i].second ;

        }

    }


}

std::vector<BoundaryCondition * > TwoDCohesiveForces::getBoundaryConditions( const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{


    Vector 	f(2, 0.) ;
    bool enrichedDof = !(p_i.getDofID() == -1) ;

    if(!enrichedDof)
        return std::vector<BoundaryCondition * >();


    Vector apparentStress(0., gp.gaussPoints.size()*(3+3*(source->spaceDimensions() == 3))) ;
    source->getState().getField( NON_ENRICHED_REAL_STRESS_FIELD, gp.gaussPoints, apparentStress, true) ;
    for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
    {
        if(enrichedDof)
        {
            Vector stress(3) ;
            for(size_t j = 0 ; j < 3 ; j++)
                stress[j] = apparentStress[i*3+j] ;
            std::vector<Variable> v ;
            v.push_back(XI);
            v.push_back(ETA);

            Matrix grad (VirtualMachine().geval(p_i, Jinv[i],v, gp.gaussPoints[i].first, true)) ;
            Vector force = (Vector)(grad*stress) ;

            double normalAmplitude = force[0]*normals[0].getX() + force[1]*normals[0].getY();
            double tangeantAmplitude = -force[0]*normals[0].getY() + force[1]*normals[0].getX();

            Vector normalForce(2) ;
            normalForce[0] = normals[0].getX()*normalAmplitude ;
            normalForce[1] = normals[0].getY()*normalAmplitude ;

            Vector tangeantForce(2) ;
            tangeantForce[0] = -normals[0].getY()*tangeantAmplitude ;
            tangeantForce[1] = normals[0].getX()*tangeantAmplitude ;

            f += normalForce*gp.gaussPoints[i].second ;
            f += tangeantForce*gp.gaussPoints[i].second ;

        }

    }

    std::vector<BoundaryCondition * > ret ;
    if(f.size() == 2)
    {
        ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, f[0]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, f[1]));
    }
    if(f.size() == 3)
    {
        ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, f[0]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, f[1]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, f[2]));
    }
    return ret ;
}

void TwoDCohesiveForces::step(double timestep, ElementState & s, double maxscore)
{

}

bool TwoDCohesiveForces::isActive() const
{
    Vector u = source->getState().getDisplacements() ;
    Point a = source->getBoundingPoint(0) + Point(u[0], u[1]) ;
    Point b = source->getBoundingPoint(source->getBoundingPoints().size()/3) + Point(u[2], u[3]) ;
    Point c = source->getBoundingPoint(2*source->getBoundingPoints().size()/3) + Point(u[4], u[5]) ;

    double newArea = Triangle(a,b,c).area() ;

    if(newArea< startArea)
        return true ;

    return false ;

}


}


