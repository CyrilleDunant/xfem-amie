
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//

#include "twod_cohesive_force.h"
#include "../mesher/delaunay.h" 
#include "../polynomial/vm_base.h" 
#include "../features/boundarycondition.h"

using namespace Mu ;



TwoDCohesiveForces::TwoDCohesiveForces(const IntegrableEntity *s, const IntegrableEntity *t, const SegmentedLine * sl) 
{
	this->time_d = false ;
	this->type = NON_LINEAR ;
	this->target = t ;
	this->source = s ;
	active = false ;
	startArea = s->area() ;
	for(size_t i = 0 ; i < sl->getBoundingPoints().size()-1 ; i++)
	{
		Segment test(sl->getBoundingPoint(i), sl->getBoundingPoint(i+1)) ;
		normals.push_back(test.normal()) ;
	}
	
}
	
Form * TwoDCohesiveForces::getCopy() const
{
	return new TwoDCohesiveForces(*this) ;
}

TwoDCohesiveForces::~TwoDCohesiveForces() { } ;

void TwoDCohesiveForces::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix &ret, VirtualMachine * vm) const
{
}
	
bool TwoDCohesiveForces::hasInducedForces() const
{
	return true ;
}

bool TwoDCohesiveForces::hasInducedMatrix() const 
{
	return false ;
}
	
	
void TwoDCohesiveForces::getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector & f) const 
{

	
	bool enrichedDof = !(p_i.getDofID() == -1) ;
	
	if(!enrichedDof)
		return ;
	

	Vector apparentStress = source->getState().getNonEnrichedStress(gp.gaussPoints,Jinv) ; 
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		if(enrichedDof)
		{
			Vector stress(3) ;
			for(size_t j = 0 ; j < 3 ; j++)
				stress[j] = apparentStress[i*3+j] ;
			std::vector<Variable> v ;
			v.push_back(XI);
			v.push_back(ETA);
			
			Matrix grad (VirtualMachine().geval(p_i, Jinv[i],v, gp.gaussPoints[i].first, true)) ;
			Vector force = (Vector)(grad*stress) ;
				
			double normalAmplitude = force[0]*normals[0].x + force[1]*normals[0].y;
			double tangeantAmplitude = -force[0]*normals[0].y + force[1]*normals[0].x;
				
			Vector normalForce(2) ;
			normalForce[0] = normals[0].x*normalAmplitude ;
			normalForce[1] = normals[0].y*normalAmplitude ;
				
			Vector tangeantForce(2) ;
			tangeantForce[0] = -normals[0].y*tangeantAmplitude ;
			tangeantForce[1] = normals[0].x*tangeantAmplitude ;

			f += normalForce*gp.gaussPoints[i].second ;
			f += tangeantForce*gp.gaussPoints[i].second ;
			
		}

	}
	

}

std::vector<BoundaryCondition * > TwoDCohesiveForces::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
	

	Vector 	f(2, 0.) ;
	bool enrichedDof = !(p_i.getDofID() == -1) ;
	
	if(!enrichedDof)
		return std::vector<BoundaryCondition * >();
	

	Vector apparentStress = source->getState().getNonEnrichedStress(gp.gaussPoints,Jinv) ; 
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		if(enrichedDof)
		{
			Vector stress(3) ;
			for(size_t j = 0 ; j < 3 ; j++)
				stress[j] = apparentStress[i*3+j] ;
			std::vector<Variable> v ;
			v.push_back(XI);
			v.push_back(ETA);
			
			Matrix grad (VirtualMachine().geval(p_i, Jinv[i],v, gp.gaussPoints[i].first, true)) ;
			Vector force = (Vector)(grad*stress) ;
				
			double normalAmplitude = force[0]*normals[0].x + force[1]*normals[0].y;
			double tangeantAmplitude = -force[0]*normals[0].y + force[1]*normals[0].x;
				
			Vector normalForce(2) ;
			normalForce[0] = normals[0].x*normalAmplitude ;
			normalForce[1] = normals[0].y*normalAmplitude ;
				
			Vector tangeantForce(2) ;
			tangeantForce[0] = -normals[0].y*tangeantAmplitude ;
			tangeantForce[1] = normals[0].x*tangeantAmplitude ;

			f += normalForce*gp.gaussPoints[i].second ;
			f += tangeantForce*gp.gaussPoints[i].second ;
			
		}

	}
	
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
	
void TwoDCohesiveForces::step(double timestep, ElementState & s) 
{
// 		this->state = currentState ;
	
	
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
		std::cout << i*ndofs << " vs " << s.getEnrichedDisplacements().size() << std::endl ;
		s.getEnrichedDisplacements()[i*ndofs] = s.getBuffer()[(i+s.getParent()->getBoundingPoints().size())*ndofs] ;
		s.getEnrichedDisplacements()[i*ndofs+offset] = s.getBuffer()[(i+s.getParent()->getBoundingPoints().size())*ndofs+offset] ;
	}
}
	
bool TwoDCohesiveForces::isActive() const 
{
	Point a = source->getBoundingPoint(0) + source->getState().getDisplacements(source->getBoundingPoint(0), false);
	Point b = source->getBoundingPoint(source->getBoundingPoints().size()/3) + source->getState().getDisplacements(source->getBoundingPoint(source->getBoundingPoints().size()/3), false);
	Point c = source->getBoundingPoint(2*source->getBoundingPoints().size()/3) + source->getState().getDisplacements(source->getBoundingPoint(2*source->getBoundingPoints().size()/3), false);
	
	double newArea = Triangle(a,b,c).area() ;
	
	if(newArea< startArea)
		return true ;
	
	return false ;
	
}





