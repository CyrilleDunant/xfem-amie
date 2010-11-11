//
// C++ Implementation: stiffness_with_imposed_deformation_and_fracture
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007 ; Alain Giorla <alain.giorla@epfl.ch>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "stiffness_with_diffusion_deformation_and_fracture.h"
#include "../features/boundarycondition.h"

using namespace Mu ;

StiffnessWithDiffusionDeformationAndFracture::StiffnessWithDiffusionDeformationAndFracture(const Matrix & rig, Vector imposedDef, FractureCriterion * crit) : LinearForm(rig, true, false, rig.numRows()/3+1), imposed(imposedDef), criterion(crit), eps(0.005)
{
	dfunc = new IsotropicLinearDamage(rig.numRows()-1, .005) ;
	criterion->setNeighbourhoodRadius(eps*6) ;
	criterion->setMaterialCharacteristicRadius(.005) ;
	frac = false ;
	init = param[0][0] ;
	change  = false ;
	previouschange = false ;
	previousDamage.resize(dfunc->damageState().size()) ; previousDamage =0 ;
	intermediateDamage.resize(dfunc->damageState().size()) ;intermediateDamage = 0 ;
	count = 0 ;
	previousPreviousDamage.resize(dfunc->damageState().size()) ;previousPreviousDamage = 0 ;
	damage = 0 ;	
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36)
		v.push_back(ZETA);
	phi = 0 ;
	previousphi = 0 ;
	diffusivity = 1e-2 ;
	reactionRate = 1e-2;
	accumulatedPhi = 0;
	ageing = false ;
} ;

StiffnessWithDiffusionDeformationAndFracture::~StiffnessWithDiffusionDeformationAndFracture() 
{ 
	delete dfunc ;
	delete criterion ;
}

FractureCriterion * StiffnessWithDiffusionDeformationAndFracture::getFractureCriterion() const
{
	return criterion ;
}

DamageModel * StiffnessWithDiffusionDeformationAndFracture::getDamageModel() const 
{
	return dfunc ;
}

void StiffnessWithDiffusionDeformationAndFracture::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	vm->ieval(Gradient(p_i) * dfunc->apply(param* (1-ageing*(1- accumulatedPhi)*.5)) * Gradient(p_j, true), gp, Jinv,v, ret) ;
}

void StiffnessWithDiffusionDeformationAndFracture::stepBack()
{
// 	if(change)
// 	{
// 		dynamic_cast<MohrCoulomb *>(criterion)->upVal /= .95 ;
// 		dynamic_cast<MohrCoulomb *>(criterion)->downVal /= .95 ;
// 	}
	change = previouschange ;
	damage.resize(previousDamage.size()) ;
	damage = previousDamage ;
	dfunc->damageState() = damage ;
	frac = dfunc->fractured() ;

	previousDamage.resize(previousPreviousDamage.size()) ;
	previousDamage = previousPreviousDamage ;
}

void StiffnessWithDiffusionDeformationAndFracture::step(double timestep, ElementState & currentState) 
{
	previouschange = change ;
	change = false ;
	currentState.getParent()->behaviourUpdated = false ;
	previousphi = phi ;
	//first, we diffuse
	DelaunayTriangle * tri = dynamic_cast<DelaunayTriangle *>(currentState.getParent()) ;
	DelaunayTetrahedron * tet = dynamic_cast<DelaunayTetrahedron *>(currentState.getParent()) ;
	double increment ;
	if(tri)
	{
		for(size_t i = 0 ; i < tri->neighbour.size() ; i++)
		{
			
			if(tri->getNeighbour(i)->isTriangle)
			{
				StiffnessWithDiffusionDeformationAndFracture * b = dynamic_cast<StiffnessWithDiffusionDeformationAndFracture *>(dynamic_cast<DelaunayTriangle *>(tri->getNeighbour(i))->getBehaviour()) ;
				if(b)
				{
				//the flux is the gradient (difference) of the 
				//phis along the surface
					std::pair<Point *, Point *> edge = tri->commonEdge(tri->getNeighbour(i)) ;
					double delta = b->previousphi - previousphi ;
					increment += delta*dist(edge.first, edge.second) ;
				}
				else if(dynamic_cast<DelaunayTriangle *>(tri->getNeighbour(i))->getBehaviour()->type == VOID_BEHAVIOUR)
				{
					std::pair<Point *, Point *> edge = tri->commonEdge(tri->getNeighbour(i)) ;
					double delta = 1. - previousphi ;
					increment += delta*dist(edge.first, edge.second) ;
				}
			}
			else // we assume we are in a source of value 1
			{
				std::pair<Point *, Point *> edge = tri->commonEdge(tri->getNeighbour(i)) ;
				double delta = 1. - previousphi ;
				increment += delta*dist(edge.first, edge.second) ;
			}
		}
	}
	if (tet)
	{
		for(size_t i = 0 ; i < tet->neighbour.size() ; i++)
		{
			if(tet->getNeighbour(i)->isTetrahedron())
			{
				StiffnessWithDiffusionDeformationAndFracture * b = dynamic_cast<StiffnessWithDiffusionDeformationAndFracture *>(dynamic_cast<DelaunayTriangle *>(tet->getNeighbour(i))->getBehaviour()) ;
				if(b)
				{
				//the flux is the gradient (difference) of the 
				//phis along the surface
					std::vector<Point *> edge = tet->commonSurface(tet->getNeighbour(i)) ;
					TriPoint t(edge[0], edge[1], edge[2]) ;
					double delta = b->previousphi - previousphi ;
					increment += delta*t.area() ;
				}
			}
			else // we assume we are in a source of value 1
			{
				std::vector<Point *> edge = tet->commonSurface(tet->getNeighbour(i)) ;
				TriPoint t(edge[0], edge[1], edge[2]) ;
				double delta = 1 - previousphi ;
				increment += delta*t.area() ;
			}
		}
	}
	
	double volume = 0 ;
	if(tri)
		volume = currentState.getParent()->area() ;
	else
		volume = currentState.getParent()->volume() ;
	
	if(std::abs( increment*timestep) > std::numeric_limits<double>::epsilon() )
		change = true ;
	phi += diffusivity*increment*timestep/volume ;
	//now phi gets consumed at the prescribed rate.
	
	if(!ageing)
	{
		accumulatedPhi += diffusivity*increment*timestep/volume ;
	}
	else
	{
		accumulatedPhi -= reactionRate*timestep*volume ;
		if(accumulatedPhi < 0)
			accumulatedPhi = 0 ;
	}
	
	if(accumulatedPhi >= 0.2)
	{
		ageing = true ;
	}
	
	if(!frac && criterion->met(currentState) )
	{
		
		dfunc->step(currentState) ;
// 		dynamic_cast<MohrCoulomb *>(criterion)->upVal *= .95 ;
// 		dynamic_cast<MohrCoulomb *>(criterion)->downVal *= .95 ;
		change = true ;
		currentState.getParent()->behaviourUpdated = true ;
		frac = dfunc->fractured() ;
	}
	previousPreviousDamage.resize(previousDamage.size()) ;
	previousPreviousDamage = previousDamage ;
	previousDamage.resize(damage.size()) ;
	previousDamage = damage ;

	Vector d = dfunc->damageState() ;
	damage.resize(d.size()) ;
	damage = d ;
}

void StiffnessWithDiffusionDeformationAndFracture::artificialDamageStep(double d)
{
	previouschange = change ;
	change = false ;

	dfunc->artificialDamageStep(d) ;
	change = true ;
	frac = dfunc->fractured() ;
	previousPreviousDamage.resize(previousDamage.size()) ;
	previousPreviousDamage = previousDamage ;
	previousDamage.resize(damage.size()) ;
	previousDamage = damage ;

	Vector d_ = dfunc->damageState() ;
	damage.resize(d_.size()) ;
	damage = d ;
}


void StiffnessWithDiffusionDeformationAndFracture::artificialPreviousDamage(Vector previous, Vector previousprevious)
{
	previousDamage.resize(damage.size()) ;
	if(previous.size() < previousDamage.size())
	{
		for(size_t i = 0 ; i < previous.size() ; i++)
			previousDamage[i] = std::min(damage[i],previous[i]) ;
		for(size_t j = previous.size() ; j < previousDamage.size() ; j++)
			previousDamage[j] = std::min(damage[j],previous[previous.size() - 1]) ;
	} else {
		for(size_t i = 0 ; i < previousDamage.size() ; i++)
			previousDamage[i] = std::min(damage[i],previous[i]) ;
	}
	previousPreviousDamage.resize(damage.size()) ;
	if(previousprevious.size() < previousPreviousDamage.size())
	{
		for(size_t i = 0 ; i < previousprevious.size() ; i++)
			previousPreviousDamage[i] = std::min(previousDamage[i],previousprevious[i]) ;
		for(size_t j = previous.size() ; j < previousPreviousDamage.size() ; j++)
			previousPreviousDamage[j] = std::min(previousDamage[j],previousprevious[previousprevious.size() - 1]) ;
	} else {
		for(size_t i = 0 ; i < previousPreviousDamage.size() ; i++)
			previousPreviousDamage[i] = std::min(previousDamage[i],previousprevious[i]) ;
	}
}


bool StiffnessWithDiffusionDeformationAndFracture::changed() const
{
	return change ;
} 

bool StiffnessWithDiffusionDeformationAndFracture::fractured() const
{
	return dfunc->fractured() ;
}


Form * StiffnessWithDiffusionDeformationAndFracture::getCopy() const 
{
	StiffnessWithDiffusionDeformationAndFracture * copy = new StiffnessWithDiffusionDeformationAndFracture(param, imposed, criterion->getCopy()) ;
	copy->damage.resize(damage.size());
	copy->damage = damage ;
	copy->criterion->setNeighbourhoodRadius(criterion->getNeighbourhoodRadius()) ;
	copy->dfunc->setMaterialCharacteristicRadius(dfunc->getCharacteristicRadius());
	return copy ;
}


Matrix StiffnessWithDiffusionDeformationAndFracture::getTensor(const Point & p) const
{
	return dfunc->apply(param) ;
}


Vector StiffnessWithDiffusionDeformationAndFracture::getImposedStress(const Point & p) const
{
	return ((param  * (1-ageing*(1- accumulatedPhi)*.5)* accumulatedPhi) * imposed) ;
}

std::vector<BoundaryCondition * > StiffnessWithDiffusionDeformationAndFracture::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
	Vector f = VirtualMachine().ieval(Gradient(p_i) * ((param * accumulatedPhi) * imposed ), gp, Jinv,v) ;
	
	std::vector<BoundaryCondition * > ret ;
	if(f.size() == 2)
	{
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementarySurface *>(s.getParent()), id, f[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementarySurface *>(s.getParent()), id, f[1]));
	}
	if(f.size() == 3)
	{
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementaryVolume *>(s.getParent()), id, f[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()), id, f[1]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()), id, f[2]));
	}
	return ret ;
}
