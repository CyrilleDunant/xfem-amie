//
// C++ Implementation: physics_base
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "physics_base.h"

using namespace Mu ;


LinearForm::LinearForm(const Matrix & p, bool t, bool s, size_t numdof ) : Form(p, t, s, numdof)
{
	this->type = PURE_LINEAR ;
}

LinearForm::~LinearForm() { } ;

bool LinearForm::fractured() const { return false  ;} ; 

void LinearForm::step(double timestep, ElementState & s)
{

}

void LinearForm::updateElementState(double timestep, ElementState & s) const
{
	if( type == VOID_BEHAVIOUR)
		return ;

	s.getPreviousPreviousDisplacements() = s.getPreviousDisplacements() ;
	s.getPreviousDisplacements() = s.getDisplacements() ;
	s.getPreviousPreviousEnrichedDisplacements() = s.getPreviousEnrichedDisplacements() ;
	s.getPreviousEnrichedDisplacements() = s.getEnrichedDisplacements() ;
	
	size_t ndofs = s.getParent()->getBehaviour()->getNumberOfDegreesOfFreedom() ;
	if(s.getParent()->getBoundingPoints().size()*ndofs != s.getDisplacements().size())
	{
		s.initialize() ;
		std::cout << "uninitialized element" << std::endl ;
	}
	
	if(s.getEnrichedDisplacements().size() != s.getParent()->getEnrichmentFunctions().size()*ndofs)
		s.getEnrichedDisplacements().resize(s.getParent()->getEnrichmentFunctions().size()*ndofs) ;
	
	for(size_t i = 0 ; i < s.getParent()->getShapeFunctions().size() ; i++)
	{
		for(size_t j = 0 ; j < ndofs ; j++)
		{
			s.getDisplacements()[i*ndofs + j] = s.getBuffer()[i*ndofs + j] ;
		}
	}
	
	for(size_t i = 0 ; i < s.getParent()->getEnrichmentFunctions().size() ; i++)
	{
		for(size_t j = 0 ; j < ndofs ; j++)
			s.getEnrichedDisplacements()[i*ndofs+j] = s.getBuffer()[(i+s.getParent()->getBoundingPoints().size())*ndofs+j] ;
	}
}

void NonLinearForm::updateElementState(double timestep, ElementState & s) const
{
	
	s.getPreviousPreviousDisplacements() = s.getPreviousDisplacements() ;
	s.getPreviousDisplacements() = s.getDisplacements() ;
	s.getPreviousPreviousEnrichedDisplacements() = s.getPreviousEnrichedDisplacements() ;
	s.getPreviousEnrichedDisplacements() = s.getEnrichedDisplacements() ;
	
	size_t ndofs = s.getParent()->getBehaviour()->getNumberOfDegreesOfFreedom() ;
	int offset = ndofs-1 ;
	
	if(s.getEnrichedDisplacements().size() != s.getParent()->getEnrichmentFunctions().size()*ndofs)
		s.getEnrichedDisplacements().resize(s.getParent()->getEnrichmentFunctions().size()*ndofs) ;
	
	for(size_t i = 0 ; i < s.getParent()->getBoundingPoints().size() ; i++)
	{
		s.getDisplacements()[i*ndofs] = s.getBuffer()[i*ndofs] ;
		s.getDisplacements()[i*ndofs+offset] = s.getBuffer()[i*ndofs+offset] ;
	}
	
	for(size_t i = 0 ; i < s.getParent()->getEnrichmentFunctions().size() ; i++)
	{
		s.getEnrichedDisplacements()[i*ndofs] = s.getBuffer()[(i+s.getParent()->getBoundingPoints().size())*ndofs] ;
		s.getEnrichedDisplacements()[i*ndofs+offset] = s.getBuffer()[(i+s.getParent()->getBoundingPoints().size())*ndofs+offset] ;
	}
}

NonLinearForm::NonLinearForm()
{
	this->type = NON_LINEAR ;
}

NonLinearForm::~NonLinearForm() { } ;

std::vector<Point> NonLinearForm::getIntegrationHints()
{
	return hints ;
}

Point NonLinearForm::getIntegrationHint(size_t i)
{
	return hints[i] ;
}

void NonLinearForm::setIntegrationHints(std::vector<Point> h)
{
	hints = h ;
}

void NonLinearForm::step(double timestep, ElementState & s)
{
}

bool NonLinearForm::fractured() const 
{
	return false ;
}


LinearFormAndConstant::LinearFormAndConstant(const Matrix & p, Matrix c) : Form(p, false, false, p.numRows()/3+1)
{
	this->type = LINEAR_AND_CONSTANT ;
}

LinearFormAndConstant::~LinearFormAndConstant() { } ;

