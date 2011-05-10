//
// C++ Implementation: kelvinvoight
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2008-2011
//
// Copyright: See COPYING file that comes with this distribution
//
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

void KelvinVoight::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{

	Matrix temp(ret) ;
	Matrix temp0(ret) ;
	vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v,ret) ;
	vm->ieval(GradientDot(p_i) * eta * Gradient(p_j, true), gp, Jinv,v,temp);
	vm->ieval(Gradient(p_i) * eta * GradientDot(p_j, true), gp, Jinv,v,temp0);

/*	std::cerr << "---" << std::endl ;
	Jinv[0].print();
	std::cerr << " " << std::endl ;
	temp0.print();*/

	ret += temp+temp0;
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


