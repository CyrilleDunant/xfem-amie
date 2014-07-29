//
// C++ Implementation: mass
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "mass.h"
#include <valarray>


using namespace Amie ;

Mass::Mass(double rho, SpaceDimensionality dim) : LinearForm(Matrix(3+3*(dim == SPACE_THREE_DIMENSIONAL), 3+3*(dim == SPACE_THREE_DIMENSIONAL)), false, false, 2+(dim == SPACE_THREE_DIMENSIONAL)), density(rho) 
{
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() > 9)
		v.push_back(ZETA);
} ;

Mass::~Mass() { } ;

void Mass::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	double m = density*vm->ieval(p_i*p_j, gp) ;
	for(size_t i = 0 ; i < v.size() ; i++)
	{
		ret[i][i] = m ;
	}
}

Form * Mass::getCopy() const 
{
	Mass * copy =  new Mass(*this) ;
	
	return copy ; 
}


