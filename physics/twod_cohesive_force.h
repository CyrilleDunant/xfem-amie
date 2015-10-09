// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __TWOD_COHESIVE_FORCE_H_
#define __TWOD_COHESIVE_FORCE_H_

#include "physics_base.h"
#include "void_form.h"
#include "stiffness.h"
#include "diffusion.h"
#include "weibull_distributed_stiffness.h"
#include "stiffness_and_fracture.h"
#include "stiffness_with_imposed_deformation.h"

namespace Amie
{



// class Diffusion2D :public LinearForm
// {
// public:
//
// 	double c ;
//
// 	double alpha ;
// 	double tau   ;
//
// 	Diffusion2D(double capacity, double conductivityX, double conductivityY) ;
//
// 	virtual ~Diffusion2D() ;
//
//
// 	virtual Matrix apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const;
//
// 	virtual Matrix apply(const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point,double> > &gp, const std::valarray<Matrix> &Jinv, Matrix & ret) const;
//
// 	virtual void step(double timestep, ElementState * currentState) ;
//
// 	/** Check for fracture state
// 	 *
// 	 * @return true if the element is fractured
// 	 */
// 	virtual bool fractured() const;
//
//
// 	/** get Copy of the behaviour
// 	 *
// 	 * @return pointer to the copy. Caller is responsible for cleaning memory
// 	 */
// 	virtual Form * getCopy() const ;
//
// 	virtual bool hasInducedForces() const ;
//
// 	virtual void getForces(const ElementState & s, const Function & p_i, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv, Vector &v) const ;
//
// } ;

class TwoDCohesiveForces : public NonLinearForm
{
public:

    std::vector<Point> normals ;
    IntegrableEntity * source ;
    IntegrableEntity * target ;
    double startArea ;

    bool active ;

    TwoDCohesiveForces(IntegrableEntity *s, IntegrableEntity *t, const SegmentedLine * sl) ;

    virtual Form * getCopy() const ;

    virtual ~TwoDCohesiveForces() ;

    virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;

    virtual bool hasInducedForces() const ;

    virtual bool hasInducedMatrix() const ;

    virtual void getForces( ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector &v)  ;
    std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;

    virtual void step(double timestep, ElementState & currentState, double maxscore) ;

    virtual bool isActive() const ;
} ;


}

#endif // __ PHYSICS_H_



