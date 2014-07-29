//
// C++ Interface: stiffness_and_fracture
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "stiffness_and_indexed_fracture.h"
#include "../mesher/delaunay.h"
#include "fracturecriteria/ruptureenergy.h"
#include "fracturecriteria/mcft.h"
#include "fracturecriteria/mohrcoulomb.h"

using namespace Amie ;


StiffnessAndIndexedFracture::StiffnessAndIndexedFracture(const Matrix & rig, FractureCriterion * crit, double eps) : LinearForm(rig, false, true, rig.numRows()/3+1), eps(eps)
{
	dfunc = new IndexedLinearDamage(rig.numRows()-1,1., crit) ; 
	criterion = crit ;
	crit->setMaterialCharacteristicRadius(eps) ;

	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36 )
	{
		v.push_back(ZETA);
	}
// 	v.push_back(TIME_VARIABLE);
} ;

void StiffnessAndIndexedFracture::setNeighbourhoodRadius(double d)
{
	criterion->setMaterialCharacteristicRadius(d);
	eps = d ;
}

StiffnessAndIndexedFracture::~StiffnessAndIndexedFracture() 
{ 
	delete criterion ;
	delete dfunc ;
} ;

FractureCriterion * StiffnessAndIndexedFracture::getFractureCriterion() const
{
	return criterion ;
}

DamageModel * StiffnessAndIndexedFracture::getDamageModel() const
{
	return dfunc ;
}

void StiffnessAndIndexedFracture::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix &ret, VirtualMachine * vm) const
{
	vm->ieval(Gradient(p_i) * dfunc->apply(param) * Gradient(p_j, true), gp, Jinv,v, ret) ;
}

void StiffnessAndIndexedFracture::step(double timestep, ElementState & currentState, double maxscore) 
{
	dfunc->step(currentState, maxscore) ;

}

bool StiffnessAndIndexedFracture::changed() const
{
	return dfunc->changed() ;
} 

bool StiffnessAndIndexedFracture::fractured() const
{
	return dfunc->fractured() ;
}

Form * StiffnessAndIndexedFracture::getCopy() const 
{
	StiffnessAndIndexedFracture * copy = new StiffnessAndIndexedFracture(param, criterion->getCopy(), criterion->getMaterialCharacteristicRadius()) ;
	copy->criterion->setMaterialCharacteristicRadius(criterion->getMaterialCharacteristicRadius()) ;
	copy->dfunc->setThresholdDamageDensity(dfunc->getThresholdDamageDensity());
	
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

Matrix StiffnessAndIndexedFracture::getTensor(const Point & p, IntegrableEntity * e, int g) const
{
	return dfunc->apply(param) ;
}
