//
// C++ Implementation: kelvinvoight
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "kelvinvoight.h"

using namespace Mu ;

KelvinVoight::KelvinVoight(const Matrix & rig, const Matrix & e) : LinearForm(rig, false, false, rig.numRows()/3+1), eta(e)
{
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() > 9)
		v.push_back(ZETA);
	v.push_back(TIME_VARIABLE);
} ;

KelvinVoight::~KelvinVoight() { } ;

Matrix KelvinVoight::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{
 /** \TODO implement GDtMtGD ieval in VirtualMachine */
	VirtualMachine vm ;
	return vm.ieval(Gradient(p_i) * param * Gradient(p_j, true), e,v) 
   /* +    vm.ieval(GradientDot(p_i) * eta * GradientDot(p_j, true), e,v)*/;
}

void KelvinVoight::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{

	Matrix temp(ret) ;
	Matrix temp0(ret) ;
	Matrix temp1(ret) ;
	Matrix temp2(ret) ;
	vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v,ret) ;
// 	vm->ieval(GradientDot(p_i) * eta * GradientDot(p_j, true), gp, Jinv,v,temp);
	vm->ieval(GradientDot(p_i) * eta * Gradient(p_j, true), gp, Jinv,v,temp);
	vm->ieval(Gradient(p_i) * eta * GradientDot(p_j, true), gp, Jinv,v,temp1);
/*	vm->print(p_i) ;
	vm->print(p_j) */;
	ret += temp+temp1;
// 	ret.print() ;
}

bool KelvinVoight::fractured() const
{
	return false ;
}

bool KelvinVoight::changed() const
{
	return false ;
} 

Form * KelvinVoight::getCopy() const 
{
	return new KelvinVoight(*this) ;
}

void KelvinVoight::getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector & f) const 
{
}

