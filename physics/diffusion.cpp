//
// C++ Implementation: stiffness
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "diffusion.h"

using namespace Mu ;

Diffusion::Diffusion(const Matrix & rig) : LinearForm(rig, false, false, 1) 
{
} ;

Diffusion::~Diffusion() { } ;

Matrix Diffusion::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{
	Matrix ret(1,1) ;
	
	ret[0][0] = VirtualMachine().ieval(VectorGradient(p_i) * param * VectorGradient(p_j, true), e) +  VirtualMachine().ieval(Differential(p_i, TIME_VARIABLE)*p_j, e) ;
	return ret ;
}

Matrix Diffusion::apply(const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point,double> > &gp, const std::valarray<Matrix> &Jinv) const
{
	Matrix ret(1,1) ;
	ret[0][0] = VirtualMachine().ieval(VectorGradient(p_i) * param * VectorGradient(p_j, true),  gp, Jinv) +  VirtualMachine().ieval(Differential(p_i, TIME_VARIABLE)*p_j, gp, Jinv) ;
	return ret ;
}

bool Diffusion::fractured() const
{
	return false ;
}

Form * Diffusion::getCopy() const 
{
	return new Diffusion(*this) ;
}

Vector Diffusion::getForces(const ElementState * s, const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv) const 
{
	return Vector(0) ;
}

