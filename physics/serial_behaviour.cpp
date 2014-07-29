#include "serial_behaviour.h"
#include "stiffness_with_imposed_deformation.h"
#include "../features/boundarycondition.h"
#include "homogenization/composite.h"

using namespace Amie ;

SerialBehaviour::SerialBehaviour( std::vector<Form *> b ) : LinearForm( b[0]->param, false, false, b[0]->getNumberOfDegreesOfFreedom())
{
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() > 9)
		v.push_back(ZETA);
	for(size_t i = 0 ; i < b.size() ; i++)
		branches.push_back(b[i]) ;
}

SerialBehaviour::SerialBehaviour(  Form * b1, Form * b2 )  : LinearForm( b1->param, false, false, b1->getNumberOfDegreesOfFreedom())
{
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() > 9)
		v.push_back(ZETA);
	branches.push_back(b1) ;
	branches.push_back(b2) ;
}

SerialBehaviour::~SerialBehaviour() 
{
	for(size_t i = 0 ; i < branches.size() ; i++)
		delete branches[i] ;
  
}

void SerialBehaviour::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const 
{
	std::vector<Matrix> tensorAtGaussPoints ;
	for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
		tensorAtGaussPoints.push_back(this->getTensor( gp.gaussPoints[g].first, nullptr, g )) ;

	vm->ieval(Gradient(p_i) * tensorAtGaussPoints * Gradient(p_j, true), gp, Jinv,v,ret) ;
}

bool SerialBehaviour::fractured() const 
{
	for(size_t i = 0 ; i < branches.size() ; i++)
	{
		if(branches[i]->fractured())
			return true ;
	}
	return false ;
}

Form * SerialBehaviour::getCopy() const 
{
	std::vector<Form *> cop ;
	for(size_t i = 0 ; i < branches.size() ; i++)
		cop.push_back( branches[i]->getCopy() ) ;
	
	SerialBehaviour* copy = new SerialBehaviour( cop ) ;
	
	return copy ; 
}

bool SerialBehaviour::hasInducedForces() const 
{
	for(size_t i = 0 ; i < branches.size() ; i++)
	{
		if(branches[i]->hasInducedForces())
			return true ;
	}
	return false ;
}

Vector SerialBehaviour::getImposedStress(const Point & p, IntegrableEntity * e, int g) const 
{
	return (Vector) (this->getTensor(p, e, g) * this->getImposedStrain( p, e, g )) ;
}

Vector SerialBehaviour::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const 
{
	Vector strain = branches[0]->getImposedStrain(p,e,g) ;
	Matrix m = branches[0]->getTensor(p,e,g) ;
	Composite::invertTensor(m) ;
	Vector ret = (Vector) (m*strain) ;
	for(size_t i = 1 ; i < branches.size() ; i++)
	{
		strain = branches[i]->getImposedStrain(p,e,g) ;
		m = branches[i]->getTensor(p,e,g) ;
		Composite::invertTensor(m) ;
		ret += (Vector) (m*strain) ;
	}
	return ret ;
}

std::vector<BoundaryCondition * > SerialBehaviour::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const 
{
	std::vector<Vector> imposedAtGaussPoints ;
	for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
		imposedAtGaussPoints.push_back( this->getImposedStress( gp.gaussPoints[g].first, nullptr, g )) ;
  
	Vector f = VirtualMachine().ieval(Gradient(p_i) * imposedAtGaussPoints , gp, Jinv,v) ;
	
	std::vector<BoundaryCondition * > ret ;
	if(f.size() == 2)
	{
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, f[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, f[1]));
	}
	if(f.size() == 3)
	{
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, f[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, f[1]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, f[2]));
	}
	return ret ;
}

void SerialBehaviour::step(double timestep, ElementState & currentState, double maxscore) 
{
	for(size_t i = 0 ; i < branches.size() ; i++)
		branches[i]->step( timestep, dynamic_cast<SerialElementState &>(currentState).getState(i), maxscore ) ;  
}

ElementState * SerialBehaviour::createElementState( IntegrableEntity * e) 
{
	std::vector<ElementState *> s ;
	for(size_t i = 0 ; i < branches.size() ; i++)
		s.push_back(branches[i]->createElementState( e )) ;
	return new SerialElementState(e, s) ;
}

void SerialBehaviour::updateElementState(double timestep, ElementState & currentState) const 
{
	LinearForm::updateElementState(timestep, currentState) ;
	for(size_t i = 0 ; i < branches.size() ; i++)
		branches[i]->updateElementState( timestep, dynamic_cast<SerialElementState &>(currentState).getState(i) ) ;  
}

Matrix SerialBehaviour::getTensor(const Point & p, IntegrableEntity * e, int g) const 
{
	Matrix m = branches[0]->getTensor( p, e , g) ;
	Composite::invertTensor(m) ;
	Matrix ret = m ;
	for(size_t i = 1 ; i < branches.size() ; i++)
	{
		m = branches[i]->getTensor( p, nullptr , g) ;
		Composite::invertTensor(m) ;
		ret += m ;
	}
	Composite::invertTensor(ret) ;
	return ret ;
}

void SerialBehaviour::preProcess( double timeStep, ElementState & currentState ) 
{
	for(size_t i = 0 ; i < branches.size() ; i++)
		branches[i]->preProcess( timeStep, dynamic_cast<SerialElementState &>(currentState).getState(i) ) ;
	param = this->getTensor( currentState.getParent()->inLocalCoordinates( currentState.getParent()->getCenter() ), nullptr, -1) ;
}
