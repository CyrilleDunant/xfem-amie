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

#include "laplacian.h"

using namespace Mu ;

Laplacian::Laplacian(const Matrix & rig) : LinearForm(rig, false, false, 1) 
{
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 9)
		v.push_back(ZETA);
} ;

Laplacian::~Laplacian() { } ;

void Laplacian::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine *vm) const
{
	ret[0][0] = vm->ieval(VectorGradient(p_i) * param * VectorGradient(p_j, true),  gp, Jinv, v)  ;

}

bool Laplacian::fractured() const
{
	return false ;
}

Form * Laplacian::getCopy() const 
{
	return new Laplacian(*this) ;
}
