//
// C++ Implementation: stiffness
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "wave.h"

using namespace Amie ;

Wave::Wave(const Matrix & rig) : LinearForm(rig, false, false, 1) 
{	
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 9)
		v.push_back(ZETA);
	v.push_back(TIME_VARIABLE);

} ;

Wave::~Wave() { } ;

void Wave::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix &ret, VirtualMachine * vm) const
{

	ret[0][0] =  vm->ieval(Differential(p_i, TIME_VARIABLE)*Differential(p_j, TIME_VARIABLE), gp, Jinv, v) 
		+ vm->ieval(VectorGradient(p_i) * param * VectorGradient(p_j, true),  gp, Jinv, v) ;

}

bool Wave::fractured() const
{
	return false ;
}

Form * Wave::getCopy() const 
{
	
	Wave * copy = new Wave(*this) ;
	
	if(getExtra2dMeshes())
	{
		for(size_t i = 0 ; i < getExtra2dMeshes()->size() ; i++)
			copy->addMesh((*getExtra2dMeshes())[i]);
	}
	if(getExtra3dMeshes())
	{
		for(size_t i = 0 ; i < getExtra3dMeshes()->size() ; i++)
			copy->addMesh((*getExtra3dMeshes())[i]);
	}
	return copy ; 
}

