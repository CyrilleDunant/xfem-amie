// C++ Interface: logarithmic_creep
//
// Description: Logarithmic visco-elastic behaviour for the Space-Time Finite Element Method
//
//
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2007-2014
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __LOGARITHMIC_CREEP_H_
#define __LOGARITHMIC_CREEP_H_

#include "physics_base.h"
#include "viscoelasticity.h"

namespace Mu
{

struct LogarithmicCreep : public Viscoelasticity
{
    // material parameters
    Matrix C ; //stiffness
    Matrix E ; //viscosity
    double tau ; //characteristic time

	double accumulatedStress ;
	double currentStress ;


    bool isPurelyElastic ;
    bool updated ;
	
    // constructor for pure elasticity
    LogarithmicCreep( const Matrix & rig) ;
    // standard constructor
    LogarithmicCreep( const Matrix & rig, const Matrix & v, double e) ;

    virtual ~LogarithmicCreep() { } ;

	virtual void apply( const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;
	virtual void applyViscous( const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm ) const ;

	virtual Form * getCopy() const ;

	virtual void print() const ;

	virtual void preProcess( double timeStep, ElementState & currentState ) ;

} ;


} ;


#endif
