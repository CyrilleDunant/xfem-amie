#include "parallel_behaviour.h"
#include "stiffness_with_imposed_deformation.h"
#include "../features/boundarycondition.h"
#include "homogenization/composite.h"

using namespace Mu ;

ParallelBehaviour::ParallelBehaviour( std::vector<Form *> b ) : LinearForm( b[0]->param, false, false, b[0]->getNumberOfDegreesOfFreedom())
{
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() > 9)
		v.push_back(ZETA);
	for(size_t i = 0 ; i < b.size() ; i++)
		branches.push_back(b[i]) ;
}

ParallelBehaviour::ParallelBehaviour(  Form * b1, Form * b2 )  : LinearForm( b1->param, false, false, b1->getNumberOfDegreesOfFreedom())
{
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() > 9)
		v.push_back(ZETA);
	branches.push_back(b1) ;
	branches.push_back(b2) ;
}

ParallelBehaviour::~ParallelBehaviour() 
{
	for(size_t i = 0 ; i < branches.size() ; i++)
		delete branches[i] ;
  
}

void ParallelBehaviour::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const 
{
	std::vector<Matrix> tensorAtGaussPoints ;
	for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
		tensorAtGaussPoints.push_back(this->getTensor( gp.gaussPoints[g].first, NULL, g )) ;
		
	vm->ieval(Gradient(p_i) * tensorAtGaussPoints * Gradient(p_j, true), gp, Jinv,v,ret) ;
}

bool ParallelBehaviour::fractured() const 
{
	for(size_t i = 0 ; i < branches.size() ; i++)
	{
		if(branches[i]->fractured())
			return true ;
	}
	return false ;
}

Form * ParallelBehaviour::getCopy() const 
{
	std::vector<Form *> copy ;
	for(size_t i = 0 ; i < branches.size() ; i++)
		copy.push_back( branches[i]->getCopy() ) ;
	return new ParallelBehaviour( copy ) ;
}

bool ParallelBehaviour::hasInducedForces() const 
{
	for(size_t i = 0 ; i < branches.size() ; i++)
	{
		if(branches[i]->hasInducedForces())
			return true ;
	}
	return false ;
}

Vector ParallelBehaviour::getImposedStress(const Point & p, IntegrableEntity * e, int g) const 
{
	Vector v(0., param.numRows()) ;
	for(size_t i = 0 ; i < branches.size() ; i++)
		v += branches[i]->getImposedStress(p, e, g) ;
	return v ; 
}

Vector ParallelBehaviour::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const 
{
	Matrix m = this->getTensor(p, e, g) ;
	Composite::invertTensor(m) ;
	return (Vector) (m * this->getImposedStress(p,e,g)) ; 
}

std::vector<BoundaryCondition * > ParallelBehaviour::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const 
{
	std::vector<Vector> imposedAtGaussPoints ;
	for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
		imposedAtGaussPoints.push_back( this->getImposedStress( gp.gaussPoints[g].first, NULL, g ) ) ;
  
	Vector f = VirtualMachine().ieval(Gradient(p_i) * imposedAtGaussPoints , gp, Jinv,v) ;
	
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

void ParallelBehaviour::step(double timestep, ElementState & currentState) 
{
	for(size_t i = 0 ; i < branches.size() ; i++)
		branches[i]->step( timestep, dynamic_cast<ParallelElementState &>(currentState).getState(i) ) ;  
}

ElementState * ParallelBehaviour::createElementState( IntegrableEntity * e) 
{
	std::vector<ElementState *> s ;
	for(size_t i = 0 ; i < branches.size() ; i++)
		s.push_back(branches[i]->createElementState( e )) ;
	return new ParallelElementState(e, s) ;
}

void ParallelBehaviour::updateElementState(double timestep, ElementState & currentState) const 
{
	LinearForm::updateElementState(timestep, currentState) ;
	for(size_t i = 0 ; i < branches.size() ; i++)
		branches[i]->updateElementState( timestep, dynamic_cast<ParallelElementState &>(currentState).getState(i) ) ;  
}

Matrix ParallelBehaviour::getTensor(const Point & p, IntegrableEntity * e, int g) const 
{
	Matrix m = branches[0]->getTensor(p, e, g) ;
	for(size_t i = 1 ; i < branches.size() ; i++)
	{
		m += branches[i]->getTensor(p, e, g) ;
	}
	return m ;
}

void ParallelBehaviour::preProcess( double timeStep, ElementState & currentState ) 
{
	for(size_t i = 0 ; i < branches.size() ; i++)
		branches[i]->preProcess( timeStep, dynamic_cast<ParallelElementState &>(currentState).getState(i) ) ;
	param = this->getTensor( currentState.getParent()->inLocalCoordinates( currentState.getParent()->getCenter() ), NULL, -1) ;
}
