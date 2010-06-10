
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//

#include "integrable_entity.h"
#include "../physics/fracturecriteria/fracturecriterion.h"
#include "../solvers/assembly.h"
#include "../features/boundarycondition.h"
using namespace Mu ;

const Vector & ElementState::getEnrichedDisplacements() const
{
	return this->enrichedDisplacements ;
}

Vector & ElementState::getEnrichedDisplacements()
{
	return this->enrichedDisplacements ;
}

IntegrableEntity::IntegrableEntity() : state(this), boundaryConditionCache(NULL)
{
	
}

void IntegrableEntity::applyBoundaryCondition(Assembly * a)
{
	if((getBehaviour()->changed() || boundaryConditionCache == NULL) && getBehaviour()->type != VOID_BEHAVIOUR)
	{
		if(boundaryConditionCache)
		{
			for(size_t i = 0 ; i < boundaryConditionCache->size() ; i++)
				delete (*boundaryConditionCache)[i] ;
		}
		
		delete boundaryConditionCache ;
	
		boundaryConditionCache = new std::vector<BoundaryCondition *>();
		std::valarray<Matrix> Jinv(getGaussPoints().gaussPoints.size()) ;

		for(size_t i = 0 ; i < getGaussPoints().gaussPoints.size() ;  i++)
		{
			getInverseJacobianMatrix( getGaussPoints().gaussPoints[i].first, Jinv[i] ) ;
		}

		
		for(size_t i = 0 ; i < getBoundingPoints().size() ; i++)
		{
			std::vector<BoundaryCondition *> boundaryConditionCachetmp = getBehaviour()->getBoundaryConditions(getState(),getBoundingPoint(i).id,  getShapeFunction(i), getGaussPoints(), Jinv) ;
			boundaryConditionCache->insert(boundaryConditionCache->end(), boundaryConditionCachetmp.begin(), boundaryConditionCachetmp.end()) ;
		}
		
		for(size_t i = 0 ; i < getEnrichmentFunctions().size() ; i++)
		{
			std::vector<BoundaryCondition *> boundaryConditionCachetmp = getBehaviour()->getBoundaryConditions(getState(),getEnrichmentFunction(i).getDofID(),  getEnrichmentFunction(i), getGaussPoints(), Jinv) ;
			boundaryConditionCache->insert(boundaryConditionCache->end(), boundaryConditionCachetmp.begin(), boundaryConditionCachetmp.end()) ;
		}
	}
	if(boundaryConditionCache)
	{
		for(size_t i = 0 ; i < boundaryConditionCache->size() ; i++)
		{
			if(get2DMesh())
				(*boundaryConditionCache)[i]->apply(a, get2DMesh()) ;
			else
				(*boundaryConditionCache)[i]->apply(a, get3DMesh()) ;
		}
	}
	
}

IntegrableEntity::~IntegrableEntity()
{
	if(boundaryConditionCache)
	{
		for(size_t i = 0 ; i < boundaryConditionCache->size() ; i++)
			delete (*boundaryConditionCache)[i] ;
	}
	
	delete boundaryConditionCache ;
}

const Vector & ElementState::getPreviousEnrichedDisplacements() const
{
	return this->previousEnrichedDisplacements ;
}

const ElementState & IntegrableEntity::getState() const
{
	return state ;
}

ElementState & IntegrableEntity::getState()
{
	return state ;
}

void ElementState::stepBack()
{
	if(!history.empty())
	{
		*this = history.back() ;
		history.pop_back() ;
	}
}

ElementState & ElementState::operator =(const ElementState &s)
{
	displacements.resize(s.getDisplacements().size()) ; displacements = s.getDisplacements() ;
	enrichedDisplacements.resize(s.getEnrichedDisplacements().size()) ; 
	enrichedDisplacements  = s.getEnrichedDisplacements();
	
	previousDisplacements.resize(s.getPreviousDisplacements().size()) ; 
	previousDisplacements = s.getPreviousDisplacements() ;
	previousEnrichedDisplacements.resize(s.getPreviousEnrichedDisplacements().size()) ; 
	previousEnrichedDisplacements = s. getPreviousEnrichedDisplacements();
	
	previousPreviousDisplacements.resize(s.getPreviousPreviousDisplacements().size()) ; 
	previousPreviousDisplacements = s.getPreviousPreviousDisplacements() ;
	previousPreviousEnrichedDisplacements.resize(s.getPreviousPreviousEnrichedDisplacements().size()) ; 
	previousPreviousEnrichedDisplacements = s.getPreviousPreviousEnrichedDisplacements() ;
	buffer.resize(s.getBuffer().size()) ; buffer = s.getBuffer() ;
	
	timePos = s.getTime();
	previousTimePos = s.getTime()-s.getDeltaTime();
	previousPreviousTimePos = s.getTime()-2*s.getDeltaTime();
	
	parent = s.getParent();
	return *this ;
}

ElementState::ElementState(const ElementState &s)
{
	displacements.resize(s.getDisplacements().size()) ; displacements = s.getDisplacements() ;
	enrichedDisplacements.resize(s.getEnrichedDisplacements().size()) ; 
	enrichedDisplacements  = s.getEnrichedDisplacements();
	
	previousDisplacements.resize(s.getPreviousDisplacements().size()) ; 
	previousDisplacements = s.getPreviousDisplacements() ;
	previousEnrichedDisplacements.resize(s.getPreviousEnrichedDisplacements().size()) ; 
	previousEnrichedDisplacements = s. getPreviousEnrichedDisplacements();
	
	previousPreviousDisplacements.resize(s.getPreviousPreviousDisplacements().size()) ; 
	previousPreviousDisplacements = s.getPreviousPreviousDisplacements() ;
	previousPreviousEnrichedDisplacements.resize(s.getPreviousPreviousEnrichedDisplacements().size()) ; 
	previousPreviousEnrichedDisplacements = s.getPreviousPreviousEnrichedDisplacements() ;
	buffer.resize(s.getBuffer().size()) ; buffer = s.getBuffer() ;
	
	timePos = s.getTime();
	previousTimePos = s.getTime()-s.getDeltaTime();
	previousPreviousTimePos = s.getTime()-2*s.getDeltaTime();
	
	parent = s.getParent();
}

Vector & ElementState::getPreviousEnrichedDisplacements() 
{
	return this->previousEnrichedDisplacements ;
}

const Vector & ElementState::getPreviousPreviousEnrichedDisplacements() const
{
	return this->previousPreviousEnrichedDisplacements ;
}

const Vector & ElementState::getPreviousDisplacements() const
{
	return this->previousDisplacements ;
}

Vector & ElementState::getPreviousPreviousEnrichedDisplacements() 
{
	return this->previousPreviousEnrichedDisplacements ;
}

Vector & ElementState::getPreviousDisplacements() 
{
	return this->previousDisplacements ;
}

const Vector & ElementState::getPreviousPreviousDisplacements() const
{
	return this->previousPreviousDisplacements ;
}

Vector ElementState::getAverageDisplacement() const
{
	FunctionMatrix disps = getDisplacementFunction() ;
	Matrix ret = VirtualMachine().ieval(disps, getParent()) ;
	Vector v(ret.numRows()) ;
	if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	{
		for(size_t i = 0 ;  i < v.size() ; i++)
			v[i] = ret[i][0]/parent->area() ;
	}
	else
	{
		for(size_t i = 0 ;  i < v.size() ; i++)
			v[i] = ret[i][0]/parent->volume() ;
	}
	return v ;
}

Vector & ElementState::getPreviousPreviousDisplacements() 
{
	return this->previousPreviousDisplacements ;
}

ElementState::ElementState(IntegrableEntity * s)
{
	
	parent = s ;
	this->timePos = 0 ;
	this->previousTimePos = 0 ;
	history.push_back(*this) ;
// 	size_t ndof = 2 ;
// // 	if(s->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
// // 		ndof = 3 ;
// 	this->displacements.resize(s->getBoundingPoints()->size()*ndof) ;
// 	this->previousDisplacements.resize(s->getBoundingPoints()->size()*ndof) ;
// 	this->previousPreviousDisplacements.resize(s->getBoundingPoints()->size()*ndof) ;
// 	this->enrichedDisplacements.resize(s->getEnrichmentFunctions().size()*ndof) ;
// 	this->previousEnrichedDisplacements.resize(s->getEnrichmentFunctions().size()*ndof) ;
// 	this->previousPreviousEnrichedDisplacements.resize(s->getEnrichmentFunctions().size()*ndof) ;
}

FunctionMatrix ElementState::getStressFunction(const Matrix &Jinv) const
{
	if (parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	{
		FunctionMatrix ret(2,2) ;
		FunctionMatrix temp(3,1) ;
		
		Function x_xi ;
		Function x_eta ;
		Function y_xi ;
		Function y_eta ;
		for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
		{
			x_xi += parent->getShapeFunction(j).d(XI)*displacements[j*2] ;
			x_eta += parent->getShapeFunction(j).d(ETA)*displacements[j*2] ;
			y_xi += parent->getShapeFunction(j).d(XI)*displacements[j*2+1] ;
			y_eta += parent->getShapeFunction(j).d(ETA)*displacements[j*2+1] ;
		}
		
		temp[0][0] = x_xi*Jinv[0][0] + x_eta*Jinv[0][1] ;
		temp[2][0] = (x_xi*Jinv[1][0] + x_eta*Jinv[1][1] + y_xi*Jinv[0][0] + y_eta*Jinv[0][1] )*0.5;	
		temp[1][0] = y_xi*Jinv[1][0] + y_eta*Jinv[1][1] ;
		
		Matrix cg(parent->getBehaviour()->getTensor(getParent()->getCenter())) ;
		
		temp *= cg ;
		
		ret[0][0] = temp[0][0] ; ret[0][1] = temp[2][0] ;
		ret[0][0] = temp[2][0] ; ret[0][1] = temp[1][0] ;
		
		return ret ;
	}
	else //if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
	{
		FunctionMatrix ret(3,3) ;
		FunctionMatrix temp(6,1) ;
		
		Function x_xi ;
		Function x_eta ;
		Function x_zeta ;
		Function y_xi ;
		Function y_eta ;
		Function y_zeta ;
		Function z_xi ;
		Function z_eta ;
		Function z_zeta ;
		for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
		{
			x_xi += parent->getShapeFunction(j).d(XI)*displacements[j*3] ;
			x_eta += parent->getShapeFunction(j).d(ETA)*displacements[j*3] ;
			x_zeta += parent->getShapeFunction(j).d(ZETA)*displacements[j*3] ;
			y_xi += parent->getShapeFunction(j).d(XI)*displacements[j*3+1] ;
			y_eta += parent->getShapeFunction(j).d(ETA)*displacements[j*3+1] ;
			y_zeta += parent->getShapeFunction(j).d(ZETA)*displacements[j*3+1] ;
			z_xi += parent->getShapeFunction(j).d(XI)*displacements[j*3+2] ;
			z_eta += parent->getShapeFunction(j).d(ETA)*displacements[j*3+2] ;
			z_zeta += parent->getShapeFunction(j).d(ZETA)*displacements[j*3+2] ;
		}
		
		temp[0][0] = x_xi*Jinv[0][0] + x_eta*Jinv[0][1] + x_zeta*Jinv[0][2];
		temp[1][0] = (y_xi)*Jinv[1][0] + (y_eta)*Jinv[1][1]  + (y_zeta)*Jinv[1][2];
		temp[2][0] = (z_xi)*Jinv[2][0] + (z_eta)*Jinv[2][1]  + (z_zeta)*Jinv[2][2]; 
		temp[3][0] =((y_xi ) * Jinv[2][0] + 
		             (y_eta ) * Jinv[2][1] +  
		             (y_zeta) * Jinv[2][2] + 
		             (z_xi  ) * Jinv[1][0] +
		             (z_eta ) * Jinv[1][1] + 
		             (z_zeta) * Jinv[1][2])*.5;	
		temp[4][0] =((x_xi ) * Jinv[2][0] + 
		             (x_eta ) * Jinv[2][1] +  
		             (x_zeta) * Jinv[2][2] + 
		             (z_xi  ) * Jinv[0][0] + 
		             (z_eta ) * Jinv[0][1] + 
		             (z_zeta) * Jinv[0][2])*.5;	
		temp[5][0] =((y_xi)  * Jinv[0][0] + 
		             (y_eta)  * Jinv[0][1] +  
		             (y_zeta) * Jinv[0][2] + 
		             (x_xi)   * Jinv[1][0] + 
		             (x_eta)  * Jinv[1][1] + 
		             (x_zeta) * Jinv[1][2])*.5;

		Matrix cg(parent->getBehaviour()->getTensor(getParent()->getCenter())) ;
		
		temp *= cg ;
		
		ret[0][0] = temp[0][0] ; ret[0][1] = temp[5][0] ; ret[0][2] = temp[4][0] ; 
		ret[1][0] = temp[5][0] ; ret[1][1] = temp[1][0] ; ret[1][2] = temp[3][0] ;
		ret[2][0] = temp[4][0] ; ret[2][1] = temp[3][0] ; ret[2][2] = temp[2][0] ;
		
		return ret ;
	}
	
}

double ElementState::getMaximumVonMisesStress() const 
{
	if(parent->getOrder() == LINEAR)
	{
		if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
		{
			Vector sigma = getPrincipalStresses( parent->getCenter()) ;
			return sqrt(((sigma[0]-sigma[1])*(sigma[0]-sigma[1]) + sigma[0] + sigma[1])/2.) ;
		}
		else if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
		{
			Mu::PointArray pts(4) ;
			pts[0] = &parent->getBoundingPoint(0) ;
			pts[1] = &parent->getBoundingPoint(1) ;
			pts[2] = &parent->getBoundingPoint(2) ;
			pts[3] = &parent->getBoundingPoint(3) ;
			Vector sigma  = getStress(pts) ;
			sigma[0] = sqrt((sigma[0]-sigma[1])*(sigma[0]-sigma[1]) + (sigma[0]-sigma[2])*(sigma[0]-sigma[2]) + (sigma[1]-sigma[2])*(sigma[1]-sigma[2]))/6. ;
			sigma[1] = sqrt((sigma[6]-sigma[7])*(sigma[6]-sigma[7]) + (sigma[6]-sigma[8])*(sigma[6]-sigma[8]) + (sigma[7]-sigma[8])*(sigma[7]-sigma[8]))/6. ;
			sigma[2] = sqrt((sigma[12]-sigma[13])*(sigma[12]-sigma[13]) + (sigma[12]-sigma[14])*(sigma[12]-sigma[14]) + (sigma[13]-sigma[14])*(sigma[13]-sigma[14]))/6. ;
			sigma[3] = sqrt((sigma[18]-sigma[19])*(sigma[18]-sigma[19]) + (sigma[18]-sigma[20])*(sigma[18]-sigma[20]) + (sigma[19]-sigma[20])*(sigma[19]-sigma[20]))/6. ;
			
			return std::max(std::max(std::max(sigma[0], sigma[1]), sigma[2]),sigma[3]) ;
		}
		else
		{
			return 0 ;
		}
	}
	else
	{
		if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
		{
			Vector principalStresses = getPrincipalStresses(parent->getBoundingPoints()) ;
			double maxS = 0 ;
			for(size_t i = 0 ; i < principalStresses.size()/2 ; i++)
			{
				maxS = std::max(maxS,
				                sqrt((
				                      (principalStresses[i*2+0]
				                       -principalStresses[i*2+1])
				                      *(principalStresses[i*2+0]
				                       -principalStresses[i*2+1]) + principalStresses[i*2+0] + principalStresses[i*2+1])/2.)
				                ) ;
			}
			
			return maxS ;

		}
		else if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
		{
			Mu::PointArray pts(4) ;
			pts[0] = &parent->getBoundingPoint(0) ;
			pts[1] = &parent->getBoundingPoint(2) ;
			pts[2] = &parent->getBoundingPoint(4) ;
			pts[3] = &parent->getBoundingPoint(6) ;
			Vector sigma  = getStress(pts) ;
			sigma[0] = sqrt((sigma[0]-sigma[1])*(sigma[0]-sigma[1]) + (sigma[0]-sigma[2])*(sigma[0]-sigma[2]) + (sigma[1]-sigma[2])*(sigma[1]-sigma[2]))/6 ;
			sigma[1] = sqrt((sigma[6]-sigma[7])*(sigma[6]-sigma[7]) + (sigma[6]-sigma[8])*(sigma[6]-sigma[8]) + (sigma[7]-sigma[8])*(sigma[7]-sigma[8]))/6 ;
			sigma[2] = sqrt((sigma[12]-sigma[13])*(sigma[12]-sigma[13]) + (sigma[12]-sigma[14])*(sigma[12]-sigma[14]) + (sigma[13]-sigma[14])*(sigma[13]-sigma[14]))/6 ;
			sigma[3] = sqrt((sigma[18]-sigma[19])*(sigma[18]-sigma[19]) + (sigma[18]-sigma[20])*(sigma[18]-sigma[20]) + (sigma[19]-sigma[20])*(sigma[19]-sigma[20]))/6 ;
			
			return std::max(std::max(std::max(sigma[0], sigma[1]), sigma[2]),sigma[3]) ;
		}
	}
	
	return 0 ;
}

FunctionMatrix ElementState::getStrainFunction(const Matrix &Jinv) const
{
	if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	{
		FunctionMatrix ret(2,2) ;
		
		Function x_xi ;
		Function x_eta ;
		Function y_xi ;
		Function y_eta ;
		for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
		{
			x_xi += parent->getShapeFunction(j).d(XI)*displacements[j*2] ;
			x_eta += parent->getShapeFunction(j).d(ETA)*displacements[j*2] ;
			y_xi += parent->getShapeFunction(j).d(XI)*displacements[j*2+1] ;
			y_eta += parent->getShapeFunction(j).d(ETA)*displacements[j*2+1] ;
		}
	
		ret[0][0] = x_xi*Jinv[0][0] + x_eta*Jinv[0][1] ;
		ret[0][1] = (x_xi*Jinv[1][0] + x_eta*Jinv[1][1] + y_xi*Jinv[0][0] + y_eta*Jinv[0][1] )*0.5;	
		ret[1][0] = (x_xi*Jinv[1][0] + x_eta*Jinv[1][1] + y_xi*Jinv[0][0] + y_eta*Jinv[0][1] )*0.5;	
		ret[1][1] = y_xi*Jinv[1][0] + y_eta*Jinv[1][1] ;
	
		return ret ;
	}
	else //if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
	{
		FunctionMatrix ret(3,3) ;
		FunctionMatrix temp(6,1) ;
		
		Function x_xi ;
		Function x_eta ;
		Function x_zeta ;
		Function y_xi ;
		Function y_eta ;
		Function y_zeta ;
		Function z_xi ;
		Function z_eta ;
		Function z_zeta ;
		for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
		{
			x_xi += parent->getShapeFunction(j).d(XI)*displacements[j*3] ;
			x_eta += parent->getShapeFunction(j).d(ETA)*displacements[j*3] ;
			x_zeta += parent->getShapeFunction(j).d(ZETA)*displacements[j*3] ;
			y_xi += parent->getShapeFunction(j).d(XI)*displacements[j*3+1] ;
			y_eta += parent->getShapeFunction(j).d(ETA)*displacements[j*3+1] ;
			y_zeta += parent->getShapeFunction(j).d(ZETA)*displacements[j*3+1] ;
			z_xi += parent->getShapeFunction(j).d(XI)*displacements[j*3+2] ;
			z_eta += parent->getShapeFunction(j).d(ETA)*displacements[j*3+2] ;
			z_zeta += parent->getShapeFunction(j).d(ZETA)*displacements[j*3+2] ;
		}
		
		temp[0][0] = x_xi*Jinv[0][0] + x_eta*Jinv[0][1] + x_zeta*Jinv[0][2];
		temp[1][0] = (y_xi)*Jinv[1][0] + (y_eta)*Jinv[1][1]  + (y_zeta)*Jinv[1][2];
		temp[2][0] = (z_xi)*Jinv[2][0] + (z_eta)*Jinv[2][1]  + (z_zeta)*Jinv[2][2]; 
		temp[3][0] =((y_xi ) * Jinv[2][0] + 
		             (y_eta ) * Jinv[2][1] +  
		             (y_zeta) * Jinv[2][2] + 
		             (z_xi  ) * Jinv[1][0] +
		             (z_eta ) * Jinv[1][1] + 
		             (z_zeta) * Jinv[1][2])*.5;	
		temp[4][0] =((x_xi ) * Jinv[2][0] + 
		             (x_eta ) * Jinv[2][1] +  
		             (x_zeta) * Jinv[2][2] + 
		             (z_xi  ) * Jinv[0][0] + 
		             (z_eta ) * Jinv[0][1] + 
		             (z_zeta) * Jinv[0][2])*.5;	
		temp[5][0] =((y_xi)  * Jinv[0][0] + 
		             (y_eta)  * Jinv[0][1] +  
		             (y_zeta) * Jinv[0][2] + 
		             (x_xi)   * Jinv[1][0] + 
		             (x_eta)  * Jinv[1][1] + 
		             (x_zeta) * Jinv[1][2])*.5;
		
		
		ret[0][0] = temp[0][0] ; ret[0][1] = temp[5][0] ; ret[0][2] = temp[4][0] ; 
		ret[1][0] = temp[5][0] ; ret[1][1] = temp[1][0] ; ret[1][2] = temp[3][0] ;
		ret[2][0] = temp[4][0] ; ret[2][1] = temp[3][0] ; ret[2][2] = temp[2][0] ;
		
		return ret ;
	}
	
}

FunctionMatrix ElementState::getDisplacementFunction() const
{
	if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	{
		FunctionMatrix ret(2,1) ;
	
		for(size_t j = 0 ; j < parent->getBoundingPoints().size() &&  j < displacements.size()*2; j++)
		{
			ret[0][0] += parent->getShapeFunction(j)*displacements[j*2] ;
			ret[1][0] += parent->getShapeFunction(j)*displacements[j*2+1] ;
		}
		
		for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() && j < enrichedDisplacements.size()*2; j++)
		{
			ret[0][0] += parent->getEnrichmentFunction(j) * enrichedDisplacements[j*2] ;
			ret[1][0] += parent->getEnrichmentFunction(j) * enrichedDisplacements[j*2+1] ;
		}
		
		return ret ;
	}
	else //if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
	{
		FunctionMatrix ret(3,1) ;
		
		for(size_t j = 0 ; j < parent->getBoundingPoints().size() &&  j < displacements.size()*3; j++)
		{
			ret[0][0] += parent->getShapeFunction(j)*displacements[j*3] ;
			ret[1][0] += parent->getShapeFunction(j)*displacements[j*3+1] ;
			ret[2][0] += parent->getShapeFunction(j)*displacements[j*3+2] ;
		}
		
		for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() && j < enrichedDisplacements.size()*3; j++)
		{
			ret[0][0] += parent->getEnrichmentFunction(j) * enrichedDisplacements[j*3] ;
			ret[1][0] += parent->getEnrichmentFunction(j) * enrichedDisplacements[j*3+1] ;
			ret[2][0] += parent->getEnrichmentFunction(j) * enrichedDisplacements[j*3+2] ;
		}
		
		return ret ;
	}
}


Vector ElementState::getStrain(const Point & p, bool local) const
{
	
	if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2)
	{
		Point p_ = parent->inLocalCoordinates(p) ;
		if(local)
			p_ = p ;
		double x_xi = 0;
		double x_eta = 0;
		double y_xi = 0;
		double y_eta = 0;
		VirtualMachine vm ;
		Vector lstrain(3) ;
		
		for(size_t j = 0 ; j < parent->getBoundingPoints().size(); j++)
		{
			double f_xi = vm.deval(parent->getShapeFunction(j),XI, p_) ;
			double f_eta = vm.deval(parent->getShapeFunction(j),ETA, p_) ;
			x_xi += f_xi*displacements[j*2] ;
			x_eta += f_eta*displacements[j*2] ;
			y_xi += f_xi*displacements[j*2+1] ;
			y_eta += f_eta*displacements[j*2+1] ;
		}
		
		for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() && j < enrichedDisplacements.size()*2; j++)
		{
			double f_xi = vm.deval( parent->getEnrichmentFunction(j),XI,p_ ) ;
			double f_eta = vm.deval( parent->getEnrichmentFunction(j),ETA,p_) ;
			
			x_xi += f_xi*enrichedDisplacements[j*2] ;
			x_eta += f_eta*enrichedDisplacements[j*2] ;
			y_xi += f_xi*enrichedDisplacements[j*2+1] ;
			y_eta += f_eta*enrichedDisplacements[j*2+1] ;

		}
		
		Matrix Jinv(2,2) ; parent->getInverseJacobianMatrix(p_, Jinv) ;
		lstrain[0] = (x_xi)*Jinv[0][0] + (x_eta)*Jinv[0][1] ;
		lstrain[1] = (y_xi)*Jinv[1][0] + (y_eta)*Jinv[1][1] ;
		lstrain[2] = 0.5*((x_xi)*Jinv[1][0] + (x_eta)*Jinv[1][1]  + (y_xi)*Jinv[0][0] + (y_eta)*Jinv[0][1] );	
				
		return lstrain ;
	}
	else if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3)
	{
		
		VirtualMachine vm ;
		Point p_ =  parent->inLocalCoordinates(p) ;
				
		if(local)
			p_ = p ;
		
		double x_xi = 0;
		double x_eta = 0;
		double x_zeta = 0;
		double y_xi = 0;
		double y_eta = 0;
		double y_zeta = 0;
		double z_xi = 0;
		double z_eta = 0;
		double z_zeta = 0;
		
		for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
		{
			double f_xi = vm.deval(parent->getShapeFunction(j),XI, p_) ;
			double f_eta = vm.deval(parent->getShapeFunction(j),ETA, p_) ;
			double f_zeta = vm.deval(parent->getShapeFunction(j),ZETA, p_) ;
			double x = displacements[j*3] ;
			double y = displacements[j*3+1] ;
			double z = displacements[j*3+2] ;
			
			x_xi   += f_xi   * x ;
			x_eta  += f_eta  * x ;
			x_zeta += f_zeta * x ;
			y_xi   += f_xi   * y ;
			y_eta  += f_eta  * y ;
			y_zeta += f_zeta * y ;
			z_xi   += f_xi   * z ;
			z_eta  += f_eta  * z ;
			z_zeta += f_zeta * z ;
		}
		
		Vector lstrain(6) ;
		
		for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
		{
			double f_xi = vm.deval( parent->getEnrichmentFunction(j),XI,p_ ) ;
			double f_eta = vm.deval( parent->getEnrichmentFunction(j),ETA,p_) ;
			double f_zeta = vm.deval( parent->getEnrichmentFunction(j),ZETA,p_) ;
			double x = enrichedDisplacements[j*3] ;
			double y = enrichedDisplacements[j*3+1] ;
			double z = enrichedDisplacements[j*3+2] ;
			
			x_xi += f_xi*x;
			x_eta += f_eta*x ;
			x_zeta += f_zeta*x ;
			y_xi += f_xi*y ;
			y_eta += f_eta*y ;
			y_zeta += f_zeta*y ;
			z_xi += f_xi*z ;
			z_eta += f_eta*z ;
			z_zeta += f_zeta*z ;
		}
		
		Matrix Jinv(3,3) ; parent->getInverseJacobianMatrix(p_, Jinv) ;
		lstrain[0] = (x_xi)*Jinv[0][0] + (x_eta)*Jinv[0][1]  + (x_zeta)*Jinv[0][2];
		lstrain[1] = (y_xi)*Jinv[1][0] + (y_eta)*Jinv[1][1]  + (y_zeta)*Jinv[1][2];
		lstrain[2] = (z_xi)*Jinv[2][0] + (z_eta)*Jinv[2][1]  + (z_zeta)*Jinv[2][2]; 
		
		lstrain[3] = 0.5*((y_xi  ) * Jinv[2][0] + 
		                  (y_eta ) * Jinv[2][1] +  
		                  (y_zeta) * Jinv[2][2] + 
		                  (z_xi  ) * Jinv[1][0] +
		                  (z_eta ) * Jinv[1][1] + 
		                  (z_zeta) * Jinv[1][2]);	
		
		lstrain[4] = 0.5*((x_xi  ) * Jinv[2][0] + 
		                  (x_eta ) * Jinv[2][1] +  
		                  (x_zeta) * Jinv[2][2] + 
		                  (z_xi  ) * Jinv[0][0] + 
		                  (z_eta ) * Jinv[0][1] + 
		                  (z_zeta) * Jinv[0][2]);	
		
		lstrain[5] = 0.5*((y_xi)   * Jinv[0][0] + 
		                  (y_eta)  * Jinv[0][1] +  
		                  (y_zeta) * Jinv[0][2] + 
		                  (x_xi)   * Jinv[1][0] + 
		                  (x_eta)  * Jinv[1][1] + 
		                  (x_zeta) * Jinv[1][2]);	

		
		return lstrain ;
	}
	int dofs = 2 ;
	if(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
		dofs = 3 ;
	int sz = dofs+(dofs-1)*(dofs-1) - dofs%2;
	return Vector(0., sz) ;
}

Vector ElementState::getStrain(const std::pair<Point, double> & p) const
{
	return getStrain(p.first) ;
}

Vector ElementState::getStress(const Point & p, bool local) const
{
	
	if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2)
	{
		size_t ndofs = parent->getBehaviour()->getNumberOfDegreesOfFreedom() ;
		if(ndofs!=2)
			return Vector(3) ;
		
		Point p_ = parent->inLocalCoordinates(p) ;
		if(local)
			p_ = p ;
		
		double x_xi = 0;
		double x_eta = 0;
		double y_xi = 0;
		double y_eta = 0;
		VirtualMachine vm ;
		Vector lstrain(3) ;

		for(size_t j = 0 ; j < parent->getShapeFunctions().size() ; j++)
		{
			double f_xi = vm.deval(parent->getShapeFunction(j),XI, p_) ;
			double f_eta = vm.deval(parent->getShapeFunction(j),ETA, p_) ;


			x_xi += f_xi*displacements[j*2] ;
			x_eta += f_eta*displacements[j*2] ;
			y_xi += f_xi*displacements[j*2+1] ;
			y_eta += f_eta*displacements[j*2+1] ;

		}
		
		for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
		{
			
			double f_xi = vm.deval( parent->getEnrichmentFunction(j),XI,p_ ) ;
			double f_eta = vm.deval( parent->getEnrichmentFunction(j),ETA,p_) ;
			
				x_xi += f_xi*enrichedDisplacements[j*2] ;
				x_eta += f_eta*enrichedDisplacements[j*2] ;
				y_xi += f_xi*enrichedDisplacements[j*2+1] ;
				y_eta += f_eta*enrichedDisplacements[j*2+1] ;
		}
		bool hasTimeSlices = (parent->timePlanes() > 1);
		Matrix Jinv(2+hasTimeSlices,2+hasTimeSlices) ; parent->getInverseJacobianMatrix(p_, Jinv) ;
		lstrain[0] = (x_xi)*Jinv[0][0] + (x_eta)*Jinv[0][1] ;
		lstrain[1] = (y_xi)*Jinv[1][0] + (y_eta)*Jinv[1][1] ;
		lstrain[2] = 0.5*((x_xi)*Jinv[1][0] + (x_eta)*Jinv[1][1]  + (y_xi)*Jinv[0][0] + (y_eta)*Jinv[0][1] );	
		
		Matrix cg(parent->getBehaviour()->getTensor(p_)) ;
		
		return lstrain*cg ;
	}
	else if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3)
	{
		Point p_ = parent->inLocalCoordinates(p) ;
		
		if(local)
			p_ = p ;
		
		double x_xi = 0;
		double x_eta = 0;
		double x_zeta = 0;
		double y_xi = 0;
		double y_eta = 0;
		double y_zeta = 0;
		double z_xi = 0;
		double z_eta = 0;
		double z_zeta = 0;
		VirtualMachine vm ;
		Vector lstrain(6) ;
		
		for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
		{
			double f_xi = vm.deval(parent->getShapeFunction(j),XI, p_) ;
			double f_eta = vm.deval(parent->getShapeFunction(j),ETA, p_) ;
			double f_zeta = vm.deval(parent->getShapeFunction(j),ZETA, p_) ;
			double x = displacements[j*3] ;
			double y = displacements[j*3+1] ;
			double z = displacements[j*3+2] ;
			
			x_xi   += f_xi   * x ;
			x_eta  += f_eta  * x ;
			x_zeta += f_zeta * x ;
			y_xi   += f_xi   * y ;
			y_eta  += f_eta  * y ;
			y_zeta += f_zeta * y ;
			z_xi   += f_xi   * z ;
			z_eta  += f_eta  * z ;
			z_zeta += f_zeta * z ;
		}
				
		for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
		{
			double f_xi = vm.deval( parent->getEnrichmentFunction(j),XI,p_ ) ;
			double f_eta = vm.deval( parent->getEnrichmentFunction(j),ETA,p_) ;
			double f_zeta = vm.deval( parent->getEnrichmentFunction(j),ZETA,p_) ;
			double x = enrichedDisplacements[j*3] ;
			double y = enrichedDisplacements[j*3+1] ;
			double z = enrichedDisplacements[j*3+2] ;
			
			x_xi   += f_xi   * x ;
			x_eta  += f_eta  * x ;
			x_zeta += f_zeta * x ;
			y_xi   += f_xi   * y ;
			y_eta  += f_eta  * y ;
			y_zeta += f_zeta * y ;
			z_xi   += f_xi   * z ;
			z_eta  += f_eta  * z ;
			z_zeta += f_zeta * z ;
		}
		
		Matrix Jinv(3,3) ; parent->getInverseJacobianMatrix(p_, Jinv) ;
		
		lstrain[0] = (x_xi)*Jinv[0][0] + (x_eta)*Jinv[0][1]  + (x_zeta)*Jinv[0][2];
		lstrain[1] = (y_xi)*Jinv[1][0] + (y_eta)*Jinv[1][1]  + (y_zeta)*Jinv[1][2];
		lstrain[2] = (z_xi)*Jinv[2][0] + (z_eta)*Jinv[2][1]  + (z_zeta)*Jinv[2][2]; 
		
		lstrain[3] = 0.5* ((y_xi ) * Jinv[2][0] + 
		                  (y_eta ) * Jinv[2][1] +  
		                  (y_zeta) * Jinv[2][2] + 
		                  (z_xi  ) * Jinv[1][0] +
		                  (z_eta ) * Jinv[1][1] + 
		                  (z_zeta) * Jinv[1][2]);	
		
		lstrain[4] = 0.5* ((x_xi ) * Jinv[2][0] + 
		                  (x_eta ) * Jinv[2][1] +  
		                  (x_zeta) * Jinv[2][2] + 
		                  (z_xi  ) * Jinv[0][0] + 
		                  (z_eta ) * Jinv[0][1] + 
		                  (z_zeta) * Jinv[0][2]);	
		
		lstrain[5] = 0.5* ((y_xi)  * Jinv[0][0] + 
		                  (y_eta)  * Jinv[0][1] +  
		                  (y_zeta) * Jinv[0][2] + 
		                  (x_xi)   * Jinv[1][0] + 
		                  (x_eta)  * Jinv[1][1] + 
		                  (x_zeta) * Jinv[1][2]);	
		
		Matrix cg(parent->getBehaviour()->getTensor(p_)) ;
		
		Vector ret = lstrain*cg ;
		return ret ;
	}
	
	int dofs = 2 ;
	if(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
		dofs = 3 ;
	int sz = dofs+(dofs-1)*(dofs-1) - dofs%2;
	return Vector(0., sz) ;
}

Vector ElementState::getNonEnrichedStress(const Point & p, bool local) const
{
	
	if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2)
	{
		
		Point p_ = parent->inLocalCoordinates(p) ;
		if(local)
			p_ = p ;
		
		double x_xi = 0;
		double x_eta = 0;
		double y_xi = 0;
		double y_eta = 0;
		VirtualMachine vm ;
		Vector lstrain(3) ;
		
		for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
		{
			double f_xi = vm.deval(parent->getShapeFunction(j),XI, p_) ;
			double f_eta = vm.deval(parent->getShapeFunction(j),ETA, p_) ;

			x_xi += f_xi*displacements[j*2] ;
			x_eta += f_eta*displacements[j*2] ;
			y_xi += f_xi*displacements[j*2+1] ;
			y_eta += f_eta*displacements[j*2+1] ;
			
		}

		
		Matrix Jinv(2,2) ; parent->getInverseJacobianMatrix(p_, Jinv) ;
		lstrain[0] = (x_xi)*Jinv[0][0] + (x_eta)*Jinv[0][1] ;
		lstrain[1] = (y_xi)*Jinv[1][0] + (y_eta)*Jinv[1][1] ;
		lstrain[2] = 0.5*((x_xi)*Jinv[1][0] + (x_eta)*Jinv[1][1]  + (y_xi)*Jinv[0][0] + (y_eta)*Jinv[0][1] );	
		
		Matrix cg(parent->getBehaviour()->getTensor(p_)) ;
		
		return lstrain*cg ;
	}
	else if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3)
	{
		Point p_ = parent->inLocalCoordinates(p) ;
		
		if(local)
			p_ = p ;
		
		double x_xi = 0;
		double x_eta = 0;
		double x_zeta = 0;
		double y_xi = 0;
		double y_eta = 0;
		double y_zeta = 0;
		double z_xi = 0;
		double z_eta = 0;
		double z_zeta = 0;
		VirtualMachine vm ;
		Vector lstrain(6) ;
		
		for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
		{
			double f_xi = vm.deval(parent->getShapeFunction(j),XI, p_) ;
			double f_eta = vm.deval(parent->getShapeFunction(j),ETA, p_) ;
			double f_zeta = vm.deval(parent->getShapeFunction(j),ZETA, p_) ;
			double x = displacements[j*3] ;
			double y = displacements[j*3+1] ;
			double z = displacements[j*3+2] ;
			
			x_xi   += f_xi   * x ;
			x_eta  += f_eta  * x ;
			x_zeta += f_zeta * x ;
			y_xi   += f_xi   * y ;
			y_eta  += f_eta  * y ;
			y_zeta += f_zeta * y ;
			z_xi   += f_xi   * z ;
			z_eta  += f_eta  * z ;
			z_zeta += f_zeta * z ;
		}
		
		Matrix Jinv(3,3) ; parent->getInverseJacobianMatrix(p_, Jinv) ;
		
		lstrain[0] = (x_xi)*Jinv[0][0] + (x_eta)*Jinv[0][1]  + (x_zeta)*Jinv[0][2];
		lstrain[1] = (y_xi)*Jinv[1][0] + (y_eta)*Jinv[1][1]  + (y_zeta)*Jinv[1][2];
		lstrain[2] = (z_xi)*Jinv[2][0] + (z_eta)*Jinv[2][1]  + (z_zeta)*Jinv[2][2]; 
		
		lstrain[3] = 0.5* ((y_xi ) * Jinv[2][0] + 
		                   (y_eta ) * Jinv[2][1] +  
		                   (y_zeta) * Jinv[2][2] + 
		                   (z_xi  ) * Jinv[1][0] +
		                   (z_eta ) * Jinv[1][1] + 
		                   (z_zeta) * Jinv[1][2]);	
		
		lstrain[4] = 0.5* ((x_xi ) * Jinv[2][0] + 
		                   (x_eta ) * Jinv[2][1] +  
		                   (x_zeta) * Jinv[2][2] + 
		                   (z_xi  ) * Jinv[0][0] + 
		                   (z_eta ) * Jinv[0][1] + 
		                   (z_zeta) * Jinv[0][2]);	
		
		lstrain[5] = 0.5* ((y_xi)  * Jinv[0][0] + 
		                   (y_eta)  * Jinv[0][1] +  
		                   (y_zeta) * Jinv[0][2] + 
		                   (x_xi)   * Jinv[1][0] + 
		                   (x_eta)  * Jinv[1][1] + 
		                   (x_zeta) * Jinv[1][2]);	
		
		Matrix cg(parent->getBehaviour()->getTensor(p_)) ;
		
		return lstrain*cg ;
	}
	
	int dofs = 2 ;
	if(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
		dofs = 3 ;
	int sz = dofs+(dofs-1)*(dofs-1) - dofs%2;
	return Vector(0., sz) ;
}

Vector ElementState::getNonEnrichedStrain(const Point & p, bool local) const
{
	
	if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2)
	{
		
		Point p_ = parent->inLocalCoordinates(p) ;
		if(local)
			p_ = p ;
		
		double x_xi = 0;
		double x_eta = 0;
		double y_xi = 0;
		double y_eta = 0;
		VirtualMachine vm ;
		Vector lstrain(3) ;
		
		for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
		{
			double f_xi = vm.deval(parent->getShapeFunction(j),XI, p_) ;
			double f_eta = vm.deval(parent->getShapeFunction(j),ETA, p_) ;

			x_xi += f_xi*displacements[j*2] ;
			x_eta += f_eta*displacements[j*2] ;
			y_xi += f_xi*displacements[j*2+1] ;
			y_eta += f_eta*displacements[j*2+1] ;
			
		}

		
		Matrix Jinv(2,2) ; parent->getInverseJacobianMatrix(p_, Jinv) ;
		lstrain[0] = (x_xi)*Jinv[0][0] + (x_eta)*Jinv[0][1] ;
		lstrain[1] = (y_xi)*Jinv[1][0] + (y_eta)*Jinv[1][1] ;
		lstrain[2] = 0.5*((x_xi)*Jinv[1][0] + (x_eta)*Jinv[1][1]  + (y_xi)*Jinv[0][0] + (y_eta)*Jinv[0][1] );	
		
		return lstrain ;
	}
	else if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3)
	{
		Point p_ = parent->inLocalCoordinates(p) ;
		
		if(local)
			p_ = p ;
		
		double x_xi = 0;
		double x_eta = 0;
		double x_zeta = 0;
		double y_xi = 0;
		double y_eta = 0;
		double y_zeta = 0;
		double z_xi = 0;
		double z_eta = 0;
		double z_zeta = 0;
		VirtualMachine vm ;
		Vector lstrain(6) ;
		
		for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
		{
			double f_xi = vm.deval(parent->getShapeFunction(j),XI, p_) ;
			double f_eta = vm.deval(parent->getShapeFunction(j),ETA, p_) ;
			double f_zeta = vm.deval(parent->getShapeFunction(j),ZETA, p_) ;
			double x = displacements[j*3] ;
			double y = displacements[j*3+1] ;
			double z = displacements[j*3+2] ;
			
			x_xi   += f_xi   * x ;
			x_eta  += f_eta  * x ;
			x_zeta += f_zeta * x ;
			y_xi   += f_xi   * y ;
			y_eta  += f_eta  * y ;
			y_zeta += f_zeta * y ;
			z_xi   += f_xi   * z ;
			z_eta  += f_eta  * z ;
			z_zeta += f_zeta * z ;
		}
		
		Matrix Jinv(3,3) ; parent->getInverseJacobianMatrix(p_, Jinv) ;
		lstrain[0] = (x_xi)*Jinv[0][0] + (x_eta)*Jinv[0][1]  + (x_zeta)*Jinv[0][2];
		lstrain[1] = (y_xi)*Jinv[1][0] + (y_eta)*Jinv[1][1]  + (y_zeta)*Jinv[1][2];
		lstrain[2] = (z_xi)*Jinv[2][0] + (z_eta)*Jinv[2][1]  + (z_zeta)*Jinv[2][2]; 
		
		lstrain[3] = 0.5* ((y_xi ) * Jinv[2][0] + 
		                   (y_eta ) * Jinv[2][1] +  
		                   (y_zeta) * Jinv[2][2] + 
		                   (z_xi  ) * Jinv[1][0] +
		                   (z_eta ) * Jinv[1][1] + 
		                   (z_zeta) * Jinv[1][2]);	
		
		lstrain[4] = 0.5* ((x_xi ) * Jinv[2][0] + 
		                   (x_eta ) * Jinv[2][1] +  
		                   (x_zeta) * Jinv[2][2] + 
		                   (z_xi  ) * Jinv[0][0] + 
		                   (z_eta ) * Jinv[0][1] + 
		                   (z_zeta) * Jinv[0][2]);	
		
		lstrain[5] = 0.5* ((y_xi)  * Jinv[0][0] + 
		                   (y_eta)  * Jinv[0][1] +  
		                   (y_zeta) * Jinv[0][2] + 
		                   (x_xi)   * Jinv[1][0] + 
		                   (x_eta)  * Jinv[1][1] + 
		                   (x_zeta) * Jinv[1][2]);	
		
		return lstrain ;
	}
	
	int dofs = 2 ;
	if(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
		dofs = 3 ;
	int sz = dofs+(dofs-1)*(dofs-1) - dofs%2;
	return Vector(0., sz) ;
}

Matrix ElementState::getStressMatrix(const Point & p, bool local) const
{
	if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2)
	{
		Matrix ret(2,2) ;
		Vector stress = getStress( p, local) ;
		ret[0][0] = stress[0] ;
		ret[1][1] = stress[1] ;
		ret[0][1] = ret[1][0]= .5*stress[2] ;
		return ret ;
	}
	else if(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3)
	{
		Matrix ret(3,3) ;
		Vector stress = getStress( p, local) ;
		ret[0][0] = stress[0] ;
		ret[1][1] = stress[1] ;
		ret[2][2] = stress[2] ;
		ret[0][2] = ret[2][0]= .5*stress[4] ;
		ret[1][2] = ret[2][1]= .5*stress[3] ;
		ret[0][1] = ret[1][0]= .5*stress[5] ;
		return ret ;
	}

	return Matrix( parent->getBehaviour()->getNumberOfDegreesOfFreedom()+1,  parent->getBehaviour()->getNumberOfDegreesOfFreedom()+1) ;
}

Matrix ElementState::getStrainMatrix(const Point & p, bool local) const
{
	if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2)
	{
		Matrix ret(2,2) ;
		Vector strain = getStrain( p, local) ;
		ret[0][0] = strain[0] ;
		ret[1][1] = strain[1] ;
		ret[0][1] = ret[1][0]= 0.5*strain[2] ;
		return ret ;
	}
	else if(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3)
	{
		Matrix ret(3,3) ;
		Vector strain = getStrain( p, local) ;
		ret[0][0] = strain[0] ;
		ret[1][1] = strain[1] ;
		ret[2][2] = strain[2] ;
		ret[0][2] = ret[2][0]= .5*strain[4] ;
		ret[1][2] = ret[2][1]= .5*strain[3] ;
		ret[0][1] = ret[1][0]= .5*strain[5] ;
		return ret ;
	}

	return Matrix( parent->getBehaviour()->getNumberOfDegreesOfFreedom()+1,  parent->getBehaviour()->getNumberOfDegreesOfFreedom()+1) ;
}

Vector ElementState::getStress(const std::pair<Point, double> & p) const
{
	return getStress(p.first, true) ;
}

Vector ElementState::getNonEnrichedStress(const Mu::PointArray & pts) const
{
	if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2)
	{
		if(parent->getBehaviour()->type == VOID_BEHAVIOUR)
		{
			return Vector(0., 3*pts.size()) ;
		}
		std::vector<Point> p_ ;
		for(size_t i = 0 ; i < pts.size() ; i++)
		{
			p_.push_back(parent->inLocalCoordinates(*pts[i])) ;
		}
		
	// 	
		VirtualMachine vm ;
		
		Vector lstrain(3*pts.size()) ;
		
		Matrix Jinv(2,2) ; parent->getInverseJacobianMatrix(p_[0], Jinv) ;
		
		
		for(size_t i = 0 ; i < pts.size() ; i++)
		{
			Matrix cg(parent->getBehaviour()->getTensor(p_[i])) ;
			
			double x_xi = 0;
			double x_eta = 0;
			double y_xi = 0;
			double y_eta = 0;
			
			for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				double f_xi = vm.deval(parent->getShapeFunction(j), XI, p_[i]) ;
				double f_eta = vm.deval(parent->getShapeFunction(j), ETA, p_[i]) ;
				
				if(j*2 < displacements.size())
				{
					x_xi += f_xi*displacements[j*2] ;
					x_eta += f_eta*displacements[j*2] ;
				}
				if(j*2+1 < displacements.size())
				{
					y_xi += f_xi*displacements[j*2+1] ;
					y_eta += f_eta*displacements[j*2+1] ;
				}
			}
			
			if(parent->isMoved())
				parent->getInverseJacobianMatrix(p_[i], Jinv) ;
			Vector llstrain(3) ;
			
			
			if(parent->isMoved())
				parent->getInverseJacobianMatrix(p_[i], Jinv) ;
			
			llstrain[0] = (x_xi)*Jinv[0][0] + (x_eta)*Jinv[0][1] ;
			llstrain[1] = (y_xi)*Jinv[1][0] + (y_eta)*Jinv[1][1] ;
			llstrain[2] = 0.5*((x_xi)*Jinv[1][0] + (x_eta)*Jinv[1][1]  + (y_xi)*Jinv[0][0] + (y_eta)*Jinv[0][1] );	
			
			llstrain=llstrain*cg ;
			
			lstrain[i*3 + 0] = llstrain[0] ;
			lstrain[i*3 + 1] = llstrain[1] ;
			lstrain[i*3 + 2] = llstrain[2] ;
			
		}
		return lstrain ;
	}
	else if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3)
	{
		if(parent->getBehaviour()->type == VOID_BEHAVIOUR)
		{
			return Vector(0., 6*pts.size()) ;
		}
		std::vector<Point> p_ ;
		for(size_t i = 0 ; i < pts.size() ; i++)
		{
			p_.push_back(parent->inLocalCoordinates(*pts[i])) ;
		}
		
	// 	
		VirtualMachine vm ;
		
		Vector lstrain(6*pts.size()) ;
		
		Matrix Jinv(3,3) ; parent->getInverseJacobianMatrix(p_[0], Jinv) ;
		
		
		for(size_t i = 0 ; i < pts.size() ; i++)
		{
			Matrix cg(parent->getBehaviour()->getTensor(p_[i])) ;
			
			double x_xi = 0;
			double x_eta = 0;
			double x_zeta = 0;
			double y_xi = 0;
			double y_eta = 0;
			double y_zeta = 0;
			double z_xi = 0;
			double z_eta = 0;
			double z_zeta = 0;
			
			for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				double f_xi = vm.deval(parent->getShapeFunction(j), XI, p_[i]) ;
				double f_eta = vm.deval(parent->getShapeFunction(j), ETA, p_[i]) ;
				double f_zeta = vm.deval(parent->getShapeFunction(j), ZETA, p_[i]) ;
				
				if(j*3 < displacements.size())
				{
					x_xi += f_xi*displacements[j*3] ;
					x_eta += f_eta*displacements[j*3] ;
					x_zeta += f_eta*displacements[j*3] ;
				}
				if(j*3+1 < displacements.size())
				{
					y_xi += f_xi*displacements[j*3+1] ;
					y_eta += f_eta*displacements[j*3+1] ;
					y_zeta += f_eta*displacements[j*3+1] ;
				}
				if(j*3+2 < displacements.size())
				{
					z_xi += f_xi*displacements[j*3+2] ;
					z_eta += f_eta*displacements[j*3+2] ;
					z_zeta += f_zeta*displacements[j*3+2] ;
				}
			}
			
			if(parent->isMoved())
				parent->getInverseJacobianMatrix(p_[i], Jinv) ;
			Vector llstrain(6) ;
			
			
			if(parent->isMoved())
				parent->getInverseJacobianMatrix(p_[i], Jinv) ;
			
			
			llstrain[0] = (x_xi)*Jinv[0][0] + (x_eta)*Jinv[0][1]  + (x_zeta)*Jinv[0][2];
			llstrain[1] = (y_xi)*Jinv[1][0] + (y_eta)*Jinv[1][1]  + (y_zeta)*Jinv[1][2];
			llstrain[2] = (z_xi)*Jinv[2][0] + (z_eta)*Jinv[2][1]  + (z_zeta)*Jinv[2][2]; 
			
			llstrain[3] = 0.5* ((y_xi  ) * Jinv[2][0] + 
			                   (y_eta ) * Jinv[2][1] +  
			                   (y_zeta) * Jinv[2][2] + 
			                   (z_xi  ) * Jinv[1][0] +
			                   (z_eta ) * Jinv[1][1] + 
			                   (z_zeta) * Jinv[1][2]);	
			
			llstrain[4] = 0.5* ((x_xi  ) * Jinv[2][0] + 
			                   (x_eta ) * Jinv[2][1] +  
			                   (x_zeta) * Jinv[2][2] + 
			                   (z_xi  ) * Jinv[0][0] + 
			                   (z_eta ) * Jinv[0][1] + 
			                   (z_zeta) * Jinv[0][2]);	
			
			llstrain[5] = 0.5* ((y_xi  ) * Jinv[0][0] + 
			                   (y_eta ) * Jinv[0][1] +  
			                   (y_zeta) * Jinv[0][2] + 
			                   (x_xi  ) * Jinv[1][0] + 
			                   (x_eta ) * Jinv[1][1] + 
			                   (x_zeta) * Jinv[1][2]);	
			
			llstrain=llstrain*cg ;
			
			lstrain[i*6 + 0] = llstrain[0] ;
			lstrain[i*6 + 1] = llstrain[1] ;
			lstrain[i*6 + 2] = llstrain[2] ;
			lstrain[i*6 + 3] = llstrain[3] ;
			lstrain[i*6 + 4] = llstrain[4] ;
			lstrain[i*6 + 5] = llstrain[5] ;
			
		}
		return lstrain ;
	}

	int dofs = 2 ;
	if(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
		dofs = 3 ;
	int sz = dofs+(dofs-1)*(dofs-1) - dofs%2;
	return Vector(0., sz*pts.size()) ;
}

Vector ElementState::getNonEnrichedStrain(const Mu::PointArray & pts) const
{
	if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2)
	{
		if(parent->getBehaviour()->type == VOID_BEHAVIOUR)
		{
			return Vector(0., 3*pts.size()) ;
		}
		std::vector<Point> p_ ;
		for(size_t i = 0 ; i < pts.size() ; i++)
		{
			p_.push_back(parent->inLocalCoordinates(*pts[i])) ;
		}
		
	// 	
		VirtualMachine vm ;
		
		Vector lstrain(3*pts.size()) ;
		
		Matrix Jinv(2,2) ; parent->getInverseJacobianMatrix(p_[0], Jinv) ;
		
		
		for(size_t i = 0 ; i < pts.size() ; i++)
		{
			
			double x_xi = 0;
			double x_eta = 0;
			double y_xi = 0;
			double y_eta = 0;
			
			for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				double f_xi = vm.deval(parent->getShapeFunction(j), XI, p_[i]) ;
				double f_eta = vm.deval(parent->getShapeFunction(j), ETA, p_[i]) ;
				
				if(j*2 < displacements.size())
				{
					x_xi += f_xi*displacements[j*2] ;
					x_eta += f_eta*displacements[j*2] ;
				}
				if(j*2+1 < displacements.size())
				{
					y_xi += f_xi*displacements[j*2+1] ;
					y_eta += f_eta*displacements[j*2+1] ;
				}
			}
			
			if(parent->isMoved())
				parent->getInverseJacobianMatrix(p_[i], Jinv) ;
			Vector llstrain(3) ;
			
			
			if(parent->isMoved())
				parent->getInverseJacobianMatrix(p_[i], Jinv) ;
			
			llstrain[0] = (x_xi)*Jinv[0][0] + (x_eta)*Jinv[0][1] ;
			llstrain[1] = (y_xi)*Jinv[1][0] + (y_eta)*Jinv[1][1] ;
			llstrain[2] = 0.5*((x_xi)*Jinv[1][0] + (x_eta)*Jinv[1][1]  + (y_xi)*Jinv[0][0] + (y_eta)*Jinv[0][1] );	
			
			lstrain[i*3 + 0] = llstrain[0] ;
			lstrain[i*3 + 1] = llstrain[1] ;
			lstrain[i*3 + 2] = llstrain[2] ;
			
		}
		return lstrain ;
	}
	else if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3)
	{
		if(parent->getBehaviour()->type == VOID_BEHAVIOUR)
		{
			return Vector(0., 6*pts.size()) ;
		}
		std::vector<Point> p_ ;
		for(size_t i = 0 ; i < pts.size() ; i++)
		{
			p_.push_back(parent->inLocalCoordinates(*pts[i])) ;
		}
		
	// 	
		VirtualMachine vm ;
		
		Vector lstrain(6*pts.size()) ;
		
		Matrix Jinv(3,3) ; parent->getInverseJacobianMatrix(p_[0], Jinv) ;
		
		
		for(size_t i = 0 ; i < pts.size() ; i++)
		{
			
			double x_xi = 0;
			double x_eta = 0;
			double x_zeta = 0;
			double y_xi = 0;
			double y_eta = 0;
			double y_zeta = 0;
			double z_xi = 0;
			double z_eta = 0;
			double z_zeta = 0;
			
			for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				double f_xi = vm.deval(parent->getShapeFunction(j), XI, p_[i]) ;
				double f_eta = vm.deval(parent->getShapeFunction(j), ETA, p_[i]) ;
				double f_zeta = vm.deval(parent->getShapeFunction(j), ZETA, p_[i]) ;
				
				if(j*3 < displacements.size())
				{
					x_xi += f_xi*displacements[j*3] ;
					x_eta += f_eta*displacements[j*3] ;
					x_zeta += f_eta*displacements[j*3] ;
				}
				if(j*3+1 < displacements.size())
				{
					y_xi += f_xi*displacements[j*3+1] ;
					y_eta += f_eta*displacements[j*3+1] ;
					y_zeta += f_eta*displacements[j*3+1] ;
				}
				if(j*3+2 < displacements.size())
				{
					z_xi += f_xi*displacements[j*3+2] ;
					z_eta += f_eta*displacements[j*3+2] ;
					z_zeta += f_zeta*displacements[j*3+2] ;
				}
			}
			
			if(parent->isMoved())
				parent->getInverseJacobianMatrix(p_[i], Jinv) ;
			Vector llstrain(6) ;
			
			
			if(parent->isMoved())
				parent->getInverseJacobianMatrix(p_[i], Jinv) ;
			
			
			llstrain[0] = (x_xi)*Jinv[0][0] + (x_eta)*Jinv[0][1]  + (x_zeta)*Jinv[0][2];
			llstrain[1] = (y_xi)*Jinv[1][0] + (y_eta)*Jinv[1][1]  + (y_zeta)*Jinv[1][2];
			llstrain[2] = (z_xi)*Jinv[2][0] + (z_eta)*Jinv[2][1]  + (z_zeta)*Jinv[2][2]; 
			
			llstrain[3] = 0.5* ((y_xi  ) * Jinv[2][0] + 
			                    (y_eta ) * Jinv[2][1] +  
			                    (y_zeta) * Jinv[2][2] + 
			                    (z_xi  ) * Jinv[1][0] +
			                    (z_eta ) * Jinv[1][1] + 
			                    (z_zeta) * Jinv[1][2]);	
			
			llstrain[4] = 0.5* ((x_xi  ) * Jinv[2][0] + 
			                    (x_eta ) * Jinv[2][1] +  
			                    (x_zeta) * Jinv[2][2] + 
			                    (z_xi  ) * Jinv[0][0] + 
			                    (z_eta ) * Jinv[0][1] + 
			                    (z_zeta) * Jinv[0][2]);	
			
			llstrain[5] = 0.5* ((y_xi  ) * Jinv[0][0] + 
			                    (y_eta ) * Jinv[0][1] +  
			                    (y_zeta) * Jinv[0][2] + 
			                    (x_xi  ) * Jinv[1][0] + 
			                    (x_eta ) * Jinv[1][1] + 
			                    (x_zeta) * Jinv[1][2]);	
						
			lstrain[i*6 + 0] = llstrain[0] ;
			lstrain[i*6 + 1] = llstrain[1] ;
			lstrain[i*6 + 2] = llstrain[2] ;
			lstrain[i*6 + 3] = llstrain[3] ;
			lstrain[i*6 + 4] = llstrain[4] ;
			lstrain[i*6 + 5] = llstrain[5] ;
			
		}
		return lstrain ;
	}

	int dofs = 2 ;
	if(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
		dofs = 3 ;
	int sz = dofs+(dofs-1)*(dofs-1) - dofs%2;
	return Vector(0., sz*pts.size()) ;
}

Vector ElementState::getNonEnrichedStress(const std::valarray<std::pair<Point,double> >& pts) const
{
	
	if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2)
	{
		if(parent->getBehaviour()->type == VOID_BEHAVIOUR)
		{
			return Vector(0., 3*pts.size()) ;
		}
	
		
	// 	
		VirtualMachine vm ;
		
		Vector lstrain(3*pts.size()) ;
		
		Matrix Jinv(2,2) ; parent->getInverseJacobianMatrix(pts[0].first, Jinv) ;
		
		
		for(size_t i = 0 ; i < pts.size() ; i++)
		{
			Matrix cg(parent->getBehaviour()->getTensor(pts[i].first)) ;
			
			double x_xi = 0;
			double x_eta = 0;
			double y_xi = 0;
			double y_eta = 0;
			
			for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				double f_xi = vm.deval(parent->getShapeFunction(j), XI, pts[i].first) ;
				double f_eta = vm.deval(parent->getShapeFunction(j), ETA, pts[i].first) ;
				
				if(j*2 < displacements.size())
				{
					x_xi += f_xi*displacements[j*2] ;
					x_eta += f_eta*displacements[j*2] ;
				}
				if(j*2+1 < displacements.size())
				{
					y_xi += f_xi*displacements[j*2+1] ;
					y_eta += f_eta*displacements[j*2+1] ;
				}
			}
			
			if(parent->isMoved())
				parent->getInverseJacobianMatrix(pts[i].first, Jinv) ;
			
			Vector llstrain(3) ;
			
			llstrain[0] = (x_xi)*Jinv[0][0] + (x_eta)*Jinv[0][1] ;
			llstrain[1] = (y_xi)*Jinv[1][0] + (y_eta)*Jinv[1][1] ;
			llstrain[2] = 0.5*((x_xi)*Jinv[1][0] + (x_eta)*Jinv[1][1]  + (y_xi)*Jinv[0][0] + (y_eta)*Jinv[0][1] );	
			
			llstrain=llstrain*cg ;
			
			lstrain[i*3 + 0] = llstrain[0] ;
			lstrain[i*3 + 1] = llstrain[1] ;
			lstrain[i*3 + 2] = llstrain[2] ;
			
		}
		return lstrain ;
	}
	else if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3)
	{
		if(parent->getBehaviour()->type == VOID_BEHAVIOUR)
		{
			return Vector(0., 6*pts.size()) ;
		}
		std::vector<Point> p_ ;
		for(size_t i = 0 ; i < pts.size() ; i++)
		{
			p_.push_back(parent->inLocalCoordinates(pts[i].first)) ;
		}
		
	// 	
		VirtualMachine vm ;
		
		Vector lstrain(6*pts.size()) ;
		
		Matrix Jinv(3,3) ; parent->getInverseJacobianMatrix(p_[0], Jinv) ;
		
		
		for(size_t i = 0 ; i < pts.size() ; i++)
		{
			Matrix cg(parent->getBehaviour()->getTensor(pts[i].first)) ;
			
			double x_xi = 0;
			double x_eta = 0;
			double x_zeta = 0;
			double y_xi = 0;
			double y_eta = 0;
			double y_zeta = 0;
			double z_xi = 0;
			double z_eta = 0;
			double z_zeta = 0;
			
			for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				double f_xi = vm.deval(parent->getShapeFunction(j), XI, p_[i]) ;
				double f_eta = vm.deval(parent->getShapeFunction(j), ETA, p_[i]) ;
				double f_zeta = vm.deval(parent->getShapeFunction(j), ZETA, p_[i]) ;
				
				if(j*3 < displacements.size())
				{
					x_xi += f_xi*displacements[j*3] ;
					x_eta += f_eta*displacements[j*3] ;
					x_zeta += f_eta*displacements[j*3] ;
				}
				if(j*3+1 < displacements.size())
				{
					y_xi += f_xi*displacements[j*3+1] ;
					y_eta += f_eta*displacements[j*3+1] ;
					y_zeta += f_eta*displacements[j*3+1] ;
				}
				if(j*3+2 < displacements.size())
				{
					z_xi += f_xi*displacements[j*3+2] ;
					z_eta += f_eta*displacements[j*3+2] ;
					z_zeta += f_zeta*displacements[j*3+2] ;
				}
			}
			
			if(parent->isMoved())
				parent->getInverseJacobianMatrix(p_[i], Jinv) ;
			Vector llstrain(6) ;
			
			
			if(parent->isMoved())
				parent->getInverseJacobianMatrix(p_[i], Jinv) ;
			
			
			llstrain[0] = (x_xi)*Jinv[0][0] + (x_eta)*Jinv[0][1]  + (x_zeta)*Jinv[0][2];
			llstrain[1] = (y_xi)*Jinv[1][0] + (y_eta)*Jinv[1][1]  + (y_zeta)*Jinv[1][2];
			llstrain[2] = (z_xi)*Jinv[2][0] + (z_eta)*Jinv[2][1]  + (z_zeta)*Jinv[2][2]; 
			
			llstrain[3] = 0.5* ((y_xi  ) * Jinv[2][0] + 
			                    (y_eta ) * Jinv[2][1] +  
			                    (y_zeta) * Jinv[2][2] + 
			                    (z_xi  ) * Jinv[1][0] +
			                    (z_eta ) * Jinv[1][1] + 
			                    (z_zeta) * Jinv[1][2]);	
			
			llstrain[4] = 0.5* ((x_xi  ) * Jinv[2][0] + 
			                    (x_eta ) * Jinv[2][1] +  
			                    (x_zeta) * Jinv[2][2] + 
			                    (z_xi  ) * Jinv[0][0] + 
			                    (z_eta ) * Jinv[0][1] + 
			                    (z_zeta) * Jinv[0][2]);	
			
			llstrain[5] = 0.5* ((y_xi  ) * Jinv[0][0] + 
			                    (y_eta ) * Jinv[0][1] +  
			                    (y_zeta) * Jinv[0][2] + 
			                    (x_xi  ) * Jinv[1][0] + 
			                    (x_eta ) * Jinv[1][1] + 
			                    (x_zeta) * Jinv[1][2]);	
			
			llstrain=llstrain*cg ;
			
			lstrain[i*6 + 0] = llstrain[0] ;
			lstrain[i*6 + 1] = llstrain[1] ;
			lstrain[i*6 + 2] = llstrain[2] ;
			lstrain[i*6 + 3] = llstrain[3] ;
			lstrain[i*6 + 4] = llstrain[4] ;
			lstrain[i*6 + 5] = llstrain[5] ;
			
		}
		return lstrain ;
	}

	int dofs = 2 ;
	if(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
		dofs = 3 ;
	int sz = dofs+(dofs-1)*(dofs-1) - dofs%2;
	return Vector(0., sz*pts.size()) ;
}

Vector ElementState::getNonEnrichedStrain(const std::valarray<std::pair<Point,double> >& pts) const
{
	if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2)
	{
		if(parent->getBehaviour()->type == VOID_BEHAVIOUR)
		{
			return Vector(0., 3*pts.size()) ;
		}
	
		
	// 	
		VirtualMachine vm ;
		
		Vector lstrain(3*pts.size()) ;
		
		Matrix Jinv(2,2) ; parent->getInverseJacobianMatrix(pts[0].first, Jinv) ;
		
		
		for(size_t i = 0 ; i < pts.size() ; i++)
		{
			
			double x_xi = 0;
			double x_eta = 0;
			double y_xi = 0;
			double y_eta = 0;
			
			for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				double f_xi = vm.deval(parent->getShapeFunction(j), XI, pts[i].first) ;
				double f_eta = vm.deval(parent->getShapeFunction(j), ETA, pts[i].first) ;
				
				if(j*2 < displacements.size())
				{
					x_xi += f_xi*displacements[j*2] ;
					x_eta += f_eta*displacements[j*2] ;
				}
				if(j*2+1 < displacements.size())
				{
					y_xi += f_xi*displacements[j*2+1] ;
					y_eta += f_eta*displacements[j*2+1] ;
				}
			}
			
			if(parent->isMoved())
				parent->getInverseJacobianMatrix(pts[i].first, Jinv) ;
			
			Vector llstrain(3) ;
			
			llstrain[0] = (x_xi)*Jinv[0][0] + (x_eta)*Jinv[0][1] ;
			llstrain[1] = (y_xi)*Jinv[1][0] + (y_eta)*Jinv[1][1] ;
			llstrain[2] = 0.5*((x_xi)*Jinv[1][0] + (x_eta)*Jinv[1][1]  + (y_xi)*Jinv[0][0] + (y_eta)*Jinv[0][1] );	
			
			lstrain[i*3 + 0] = llstrain[0] ;
			lstrain[i*3 + 1] = llstrain[1] ;
			lstrain[i*3 + 2] = llstrain[2] ;
			
		}
		return lstrain ;
	}
	else if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3)
	{
		if(parent->getBehaviour()->type == VOID_BEHAVIOUR)
		{
			return Vector(0., 6*pts.size()) ;
		}
		std::vector<Point> p_ ;
		for(size_t i = 0 ; i < pts.size() ; i++)
		{
			p_.push_back(parent->inLocalCoordinates(pts[i].first)) ;
		}
		
	// 	
		VirtualMachine vm ;
		
		Vector lstrain(6*pts.size()) ;
		
		Matrix Jinv(3,3) ; parent->getInverseJacobianMatrix(p_[0], Jinv) ;
		
		
		for(size_t i = 0 ; i < pts.size() ; i++)
		{
			
			double x_xi = 0;
			double x_eta = 0;
			double x_zeta = 0;
			double y_xi = 0;
			double y_eta = 0;
			double y_zeta = 0;
			double z_xi = 0;
			double z_eta = 0;
			double z_zeta = 0;
			
			for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				double f_xi = vm.deval(parent->getShapeFunction(j), XI, p_[i]) ;
				double f_eta = vm.deval(parent->getShapeFunction(j), ETA, p_[i]) ;
				double f_zeta = vm.deval(parent->getShapeFunction(j), ZETA, p_[i]) ;
				
				if(j*3 < displacements.size())
				{
					x_xi += f_xi*displacements[j*3] ;
					x_eta += f_eta*displacements[j*3] ;
					x_zeta += f_eta*displacements[j*3] ;
				}
				if(j*3+1 < displacements.size())
				{
					y_xi += f_xi*displacements[j*3+1] ;
					y_eta += f_eta*displacements[j*3+1] ;
					y_zeta += f_eta*displacements[j*3+1] ;
				}
				if(j*3+2 < displacements.size())
				{
					z_xi += f_xi*displacements[j*3+2] ;
					z_eta += f_eta*displacements[j*3+2] ;
					z_zeta += f_zeta*displacements[j*3+2] ;
				}
			}
			
			if(parent->isMoved())
				parent->getInverseJacobianMatrix(p_[i], Jinv) ;
			Vector llstrain(6) ;
			
			
			if(parent->isMoved())
				parent->getInverseJacobianMatrix(p_[i], Jinv) ;
			
			
			llstrain[0] = (x_xi)*Jinv[0][0] + (x_eta)*Jinv[0][1]  + (x_zeta)*Jinv[0][2];
			llstrain[1] = (y_xi)*Jinv[1][0] + (y_eta)*Jinv[1][1]  + (y_zeta)*Jinv[1][2];
			llstrain[2] = (z_xi)*Jinv[2][0] + (z_eta)*Jinv[2][1]  + (z_zeta)*Jinv[2][2]; 
			
			llstrain[3] = 0.5* ((y_xi  ) * Jinv[2][0] + 
			                    (y_eta ) * Jinv[2][1] +  
			                    (y_zeta) * Jinv[2][2] + 
			                    (z_xi  ) * Jinv[1][0] +
			                    (z_eta ) * Jinv[1][1] + 
			                    (z_zeta) * Jinv[1][2]);	
			
			llstrain[4] = 0.5* ((x_xi  ) * Jinv[2][0] + 
			                    (x_eta ) * Jinv[2][1] +  
			                    (x_zeta) * Jinv[2][2] + 
			                    (z_xi  ) * Jinv[0][0] + 
			                    (z_eta ) * Jinv[0][1] + 
			                    (z_zeta) * Jinv[0][2]);	
			
			llstrain[5] = 0.5* ((y_xi  ) * Jinv[0][0] + 
			                    (y_eta ) * Jinv[0][1] +  
			                    (y_zeta) * Jinv[0][2] + 
			                    (x_xi  ) * Jinv[1][0] + 
			                    (x_eta ) * Jinv[1][1] + 
			                    (x_zeta) * Jinv[1][2]);	
						
			lstrain[i*6 + 0] = llstrain[0] ;
			lstrain[i*6 + 1] = llstrain[1] ;
			lstrain[i*6 + 2] = llstrain[2] ;
			lstrain[i*6 + 3] = llstrain[3] ;
			lstrain[i*6 + 4] = llstrain[4] ;
			lstrain[i*6 + 5] = llstrain[5] ;
			
		}
		return lstrain ;
	}
	int dofs = 2 ;
	if(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
		dofs = 3 ;
	int sz = dofs+(dofs-1)*(dofs-1) - dofs%2;
	return Vector(0., sz*pts.size()) ;
}

Vector ElementState::getNonEnrichedStress(const std::valarray<std::pair<Point,double> >& pts,  const std::valarray<Matrix> &Jinv) const
{
	if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2)
	{
		if(parent->getBehaviour()->type == VOID_BEHAVIOUR)
		{
			return Vector(0., 3*pts.size()) ;
		}
		
		
	// 	
		VirtualMachine vm ;
		
		Vector lstrain(3*pts.size()) ;
		
		
		for(size_t i = 0 ; i < pts.size() ; i++)
		{
			Matrix cg(parent->getBehaviour()->getTensor(pts[i].first)) ;
			
			double x_xi = 0;
			double x_eta = 0;
			double y_xi = 0;
			double y_eta = 0;
			
			for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				double f_xi = vm.deval(parent->getShapeFunction(j), XI, pts[i].first) ;
				double f_eta = vm.deval(parent->getShapeFunction(j), ETA, pts[i].first) ;
				
				if(j*2 < displacements.size())
				{
					x_xi += f_xi*displacements[j*2] ;
					x_eta += f_eta*displacements[j*2] ;
				}
				if(j*2+1 < displacements.size())
				{
					y_xi += f_xi*displacements[j*2+1] ;
					y_eta += f_eta*displacements[j*2+1] ;
				}
			}
			
			Vector llstrain(3) ;
			
			llstrain[0] = (x_xi)*Jinv[i][0][0] + (x_eta)*Jinv[i][0][1] ;
			llstrain[1] = (y_xi)*Jinv[i][1][0] + (y_eta)*Jinv[i][1][1] ;
			llstrain[2] = 0.5*((x_xi)*Jinv[i][1][0] + (x_eta)*Jinv[i][1][1]  + (y_xi)*Jinv[i][0][0] + (y_eta)*Jinv[i][0][1] );	
			
			llstrain=llstrain*cg ;
			
			lstrain[i*3 + 0] = llstrain[0] ;
			lstrain[i*3 + 1] = llstrain[1] ;
			lstrain[i*3 + 2] = llstrain[2] ;
			
		}
		return lstrain ;
	}
	else if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3)
	{
		if(parent->getBehaviour()->type == VOID_BEHAVIOUR)
		{
			return Vector(0., 6*pts.size()) ;
		}
		std::vector<Point> p_ ;
		for(size_t i = 0 ; i < pts.size() ; i++)
		{
			p_.push_back(parent->inLocalCoordinates(pts[i].first)) ;
		}
		
	// 	
		VirtualMachine vm ;
		
		Vector lstrain(6*pts.size()) ;		
		
		for(size_t i = 0 ; i < pts.size() ; i++)
		{
			Matrix cg(parent->getBehaviour()->getTensor(pts[i].first)) ;
			
			double x_xi = 0;
			double x_eta = 0;
			double x_zeta = 0;
			double y_xi = 0;
			double y_eta = 0;
			double y_zeta = 0;
			double z_xi = 0;
			double z_eta = 0;
			double z_zeta = 0;
			
			for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				double f_xi = vm.deval(parent->getShapeFunction(j), XI, p_[i]) ;
				double f_eta = vm.deval(parent->getShapeFunction(j), ETA, p_[i]) ;
				double f_zeta = vm.deval(parent->getShapeFunction(j), ZETA, p_[i]) ;
				
				if(j*3 < displacements.size())
				{
					x_xi += f_xi*displacements[j*3] ;
					x_eta += f_eta*displacements[j*3] ;
					x_zeta += f_eta*displacements[j*3] ;
				}
				if(j*3+1 < displacements.size())
				{
					y_xi += f_xi*displacements[j*3+1] ;
					y_eta += f_eta*displacements[j*3+1] ;
					y_zeta += f_eta*displacements[j*3+1] ;
				}
				if(j*3+2 < displacements.size())
				{
					z_xi += f_xi*displacements[j*3+2] ;
					z_eta += f_eta*displacements[j*3+2] ;
					z_zeta += f_zeta*displacements[j*3+2] ;
				}
			}
			
			Vector llstrain(6) ;
			
			llstrain[0] = (x_xi)*Jinv[i][0][0] + (x_eta)*Jinv[i][0][1]  + (x_zeta)*Jinv[i][0][2];
			llstrain[1] = (y_xi)*Jinv[i][1][0] + (y_eta)*Jinv[i][1][1]  + (y_zeta)*Jinv[i][1][2];
			llstrain[2] = (z_xi)*Jinv[i][2][0] + (z_eta)*Jinv[i][2][1]  + (z_zeta)*Jinv[i][2][2]; 
			
			llstrain[3] = 0.5* ((y_xi  ) * Jinv[i][2][0] + 
			                    (y_eta ) * Jinv[i][2][1] +  
			                    (y_zeta) * Jinv[i][2][2] + 
			                    (z_xi  ) * Jinv[i][1][0] +
			                    (z_eta ) * Jinv[i][1][1] + 
			                    (z_zeta) * Jinv[i][1][2]);	
			
			llstrain[4] = 0.5* ((x_xi  ) * Jinv[i][2][0] + 
			                    (x_eta ) * Jinv[i][2][1] +  
			                    (x_zeta) * Jinv[i][2][2] + 
			                    (z_xi  ) * Jinv[i][0][0] + 
			                    (z_eta ) * Jinv[i][0][1] + 
			                    (z_zeta) * Jinv[i][0][2]);	
			
			llstrain[5] = 0.5* ((y_xi  ) * Jinv[i][0][0] + 
			                    (y_eta ) * Jinv[i][0][1] +  
			                    (y_zeta) * Jinv[i][0][2] + 
			                    (x_xi  ) * Jinv[i][1][0] + 
			                    (x_eta ) * Jinv[i][1][1] + 
			                    (x_zeta) * Jinv[i][1][2]);	
			
			
			llstrain=llstrain*cg ;
			
			lstrain[i*6 + 0] = llstrain[0] ;
			lstrain[i*6 + 1] = llstrain[1] ;
			lstrain[i*6 + 2] = llstrain[2] ;
			lstrain[i*6 + 3] = llstrain[3] ;
			lstrain[i*6 + 4] = llstrain[4] ;
			lstrain[i*6 + 5] = llstrain[5] ;
			
		}
		return lstrain ;
	}

	int dofs = 2 ;
	if(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
		dofs = 3 ;
	int sz = dofs+(dofs-1)*(dofs-1) - dofs%2;
	return Vector(0., sz*pts.size()) ;
}

Vector ElementState::getNonEnrichedStrain(const std::valarray<std::pair<Point,double> >& pts,  const std::valarray<Matrix> &Jinv) const
{
	if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2)
	{
		if(parent->getBehaviour()->type == VOID_BEHAVIOUR)
		{
			return Vector(0., 3*pts.size()) ;
		}
		
		
	// 	
		VirtualMachine vm ;
		
		Vector lstrain(3*pts.size()) ;
		
		
		for(size_t i = 0 ; i < pts.size() ; i++)
		{
			
			double x_xi = 0;
			double x_eta = 0;
			double y_xi = 0;
			double y_eta = 0;
			
			for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				double f_xi = vm.deval(parent->getShapeFunction(j), XI, pts[i].first) ;
				double f_eta = vm.deval(parent->getShapeFunction(j), ETA, pts[i].first) ;
				
				if(j*2 < displacements.size())
				{
					x_xi += f_xi*displacements[j*2] ;
					x_eta += f_eta*displacements[j*2] ;
				}
				if(j*2+1 < displacements.size())
				{
					y_xi += f_xi*displacements[j*2+1] ;
					y_eta += f_eta*displacements[j*2+1] ;
				}
			}
			
			Vector llstrain(3) ;
			
			llstrain[0] = (x_xi)*Jinv[i][0][0] + (x_eta)*Jinv[i][0][1] ;
			llstrain[1] = (y_xi)*Jinv[i][1][0] + (y_eta)*Jinv[i][1][1] ;
			llstrain[2] = 0.5*((x_xi)*Jinv[i][1][0] + (x_eta)*Jinv[i][1][1]  + (y_xi)*Jinv[i][0][0] + (y_eta)*Jinv[i][0][1] );	
	
			
			lstrain[i*3 + 0] = llstrain[0] ;
			lstrain[i*3 + 1] = llstrain[1] ;
			lstrain[i*3 + 2] = llstrain[2] ;
			
		}
		return lstrain ;
	}
	else if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3)
	{
		if(parent->getBehaviour()->type == VOID_BEHAVIOUR)
		{
			return Vector(0., 6*pts.size()) ;
		}
		std::vector<Point> p_ ;
		for(size_t i = 0 ; i < pts.size() ; i++)
		{
			p_.push_back(parent->inLocalCoordinates(pts[i].first)) ;
		}
		
	// 	
		VirtualMachine vm ;
		
		Vector lstrain(6*pts.size()) ;		
		
		for(size_t i = 0 ; i < pts.size() ; i++)
		{
			double x_xi = 0;
			double x_eta = 0;
			double x_zeta = 0;
			double y_xi = 0;
			double y_eta = 0;
			double y_zeta = 0;
			double z_xi = 0;
			double z_eta = 0;
			double z_zeta = 0;
			
			for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				double f_xi = vm.deval(parent->getShapeFunction(j), XI, p_[i]) ;
				double f_eta = vm.deval(parent->getShapeFunction(j), ETA, p_[i]) ;
				double f_zeta = vm.deval(parent->getShapeFunction(j), ZETA, p_[i]) ;
				
				if(j*3 < displacements.size())
				{
					x_xi += f_xi*displacements[j*3] ;
					x_eta += f_eta*displacements[j*3] ;
					x_zeta += f_eta*displacements[j*3] ;
				}
				if(j*3+1 < displacements.size())
				{
					y_xi += f_xi*displacements[j*3+1] ;
					y_eta += f_eta*displacements[j*3+1] ;
					y_zeta += f_eta*displacements[j*3+1] ;
				}
				if(j*3+2 < displacements.size())
				{
					z_xi += f_xi*displacements[j*3+2] ;
					z_eta += f_eta*displacements[j*3+2] ;
					z_zeta += f_zeta*displacements[j*3+2] ;
				}
			
			}
			
			Vector llstrain(6) ;
			
			
			llstrain[0] = (x_xi)*Jinv[i][0][0] + (x_eta)*Jinv[i][0][1]  + (x_zeta)*Jinv[i][0][2];
			llstrain[1] = (y_xi)*Jinv[i][1][0] + (y_eta)*Jinv[i][1][1]  + (y_zeta)*Jinv[i][1][2];
			llstrain[2] = (z_xi)*Jinv[i][2][0] + (z_eta)*Jinv[i][2][1]  + (z_zeta)*Jinv[i][2][2]; 
			
			llstrain[3] = 0.5* ((y_xi  ) * Jinv[i][2][0] + 
			                    (y_eta ) * Jinv[i][2][1] +  
			                    (y_zeta) * Jinv[i][2][2] + 
			                    (z_xi  ) * Jinv[i][1][0] +
			                    (z_eta ) * Jinv[i][1][1] + 
			                    (z_zeta) * Jinv[i][1][2]);	
			
			llstrain[4] = 0.5* ((x_xi  ) * Jinv[i][2][0] + 
			                    (x_eta ) * Jinv[i][2][1] +  
			                    (x_zeta) * Jinv[i][2][2] + 
			                    (z_xi  ) * Jinv[i][0][0] + 
			                    (z_eta ) * Jinv[i][0][1] + 
			                    (z_zeta) * Jinv[i][0][2]);	
			
			llstrain[5] = 0.5* ((y_xi  ) * Jinv[i][0][0] + 
			                    (y_eta ) * Jinv[i][0][1] +  
			                    (y_zeta) * Jinv[i][0][2] + 
			                    (x_xi  ) * Jinv[i][1][0] + 
			                    (x_eta ) * Jinv[i][1][1] + 
			                    (x_zeta) * Jinv[i][1][2]);	
			
					
			lstrain[i*6 + 0] = llstrain[0] ;
			lstrain[i*6 + 1] = llstrain[1] ;
			lstrain[i*6 + 2] = llstrain[2] ;
			lstrain[i*6 + 3] = llstrain[3] ;
			lstrain[i*6 + 4] = llstrain[4] ;
			lstrain[i*6 + 5] = llstrain[5] ;
		}
		return lstrain ;
	}
	
	int dofs = 2 ;
	if(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
		dofs = 3 ;
	int sz = dofs+(dofs-1)*(dofs-1) - dofs%2;
	return Vector(0., sz*pts.size()) ;
}


Vector ElementState::getStress(const Mu::PointArray & pts) const
{
	if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2)
	{
		if(parent->getBehaviour()->type == VOID_BEHAVIOUR)
		{
			return Vector(0., 3*pts.size()) ;
		}
		std::vector<Point> p_ ;
		for(size_t i = 0 ; i < pts.size() ; i++)
		{
			p_.push_back(parent->inLocalCoordinates(*pts[i])) ;
		}
		
	// 	
		VirtualMachine vm ;
		
		Vector lstrain(3*pts.size()) ;
		
		Matrix Jinv(2,2) ; parent->getInverseJacobianMatrix(p_[0], Jinv) ;
		
		
		for(size_t i = 0 ; i < pts.size() ; i++)
		{
			Matrix cg (parent->getBehaviour()->getTensor(p_[i])) ;
			
			double x_xi = 0;
			double x_eta = 0;
			double y_xi = 0;
			double y_eta = 0;
			
			for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				double f_xi = vm.deval(parent->getShapeFunction(j), XI, p_[i]) ;
				double f_eta = vm.deval(parent->getShapeFunction(j), ETA, p_[i]) ;
				
				if(j*2 < displacements.size())
				{
					x_xi += f_xi*displacements[j*2] ;
					x_eta += f_eta*displacements[j*2] ;
				}
				if(j*2+1 < displacements.size())
				{
					y_xi += f_xi*displacements[j*2+1] ;
					y_eta += f_eta*displacements[j*2+1] ;
				}
			}
			
			if(parent->isMoved())
				parent->getInverseJacobianMatrix(p_[i], Jinv) ;
			Vector llstrain(3) ;
			
			double delta_x_xi  = 0;
			double delta_x_eta  = 0;
			double delta_y_xi  = 0;
			double delta_y_eta  = 0;
			
			for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
			{
				double f_xi = vm.deval(parent->getEnrichmentFunction(j), XI, p_[i]) ;
				double f_eta = vm.deval(parent->getEnrichmentFunction(j), ETA, p_[i]) ;
				
				if(j*2 < enrichedDisplacements.size())
				{
					delta_x_xi += f_xi*enrichedDisplacements[j*2] ;
					delta_x_eta += f_eta*enrichedDisplacements[j*2] ;
				}
				
				if(j*2+1 < enrichedDisplacements.size())
				{
					delta_y_xi += f_xi*enrichedDisplacements[j*2+1] ;
					delta_y_eta += f_eta*enrichedDisplacements[j*2+1] ;
				}
			}
			
			if(parent->isMoved())
				parent->getInverseJacobianMatrix(p_[i], Jinv) ;
			
			llstrain[0] = (x_xi+delta_x_xi)*Jinv[0][0] + (x_eta+delta_x_eta)*Jinv[0][1] ;
			llstrain[1] = (y_xi+delta_y_xi)*Jinv[1][0] + (y_eta+delta_y_eta)*Jinv[1][1] ;
			llstrain[2] = 0.5*((x_xi+delta_x_xi)*Jinv[1][0] + (x_eta+delta_x_eta)*Jinv[1][1]  + (y_xi+delta_y_xi)*Jinv[0][0] + (y_eta+delta_y_eta)*Jinv[0][1] );	
			
			llstrain=llstrain*cg ;
			
			lstrain[i*3 + 0] = llstrain[0] ;
			lstrain[i*3 + 1] = llstrain[1] ;
			lstrain[i*3 + 2] = llstrain[2] ;
			
		}
		return lstrain ;
	}
	else if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3)
	{
		if(parent->getBehaviour()->type == VOID_BEHAVIOUR)
		{
			return Vector(0., 6*pts.size()) ;
		}
		std::vector<Point> p_ ;
		for(size_t i = 0 ; i < pts.size() ; i++)
		{
			p_.push_back(parent->inLocalCoordinates(*pts[i])) ;
		}
		
	// 	
		VirtualMachine vm ;
		
		Vector lstrain(6*pts.size()) ;
		
		Matrix Jinv(3,3) ; parent->getInverseJacobianMatrix(p_[0], Jinv) ;
		
		
		for(size_t i = 0 ; i < pts.size() ; i++)
		{
			Matrix cg(parent->getBehaviour()->getTensor(p_[i])) ;
			
			double x_xi = 0;
			double x_eta = 0;
			double x_zeta = 0;
			double y_xi = 0;
			double y_eta = 0;
			double y_zeta = 0;
			double z_xi = 0;
			double z_eta = 0;
			double z_zeta = 0;
			
			for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				double f_xi = vm.deval(parent->getShapeFunction(j), XI, p_[i]) ;
				double f_eta = vm.deval(parent->getShapeFunction(j), ETA, p_[i]) ;
				double f_zeta = vm.deval(parent->getShapeFunction(j), ZETA, p_[i]) ;
				
				if(j*3 < displacements.size())
				{
					x_xi += f_xi*displacements[j*3] ;
					x_eta += f_eta*displacements[j*3] ;
					x_zeta += f_zeta*displacements[j*3] ;
				}
				if(j*3+1 < displacements.size())
				{
					y_xi += f_xi*displacements[j*3+1] ;
					y_eta += f_eta*displacements[j*3+1] ;
					y_zeta += f_zeta*displacements[j*3+1] ;
				}
				if(j*3+2 < displacements.size())
				{
					z_xi += f_xi*displacements[j*3+2] ;
					z_eta += f_eta*displacements[j*3+2] ;
					z_zeta += f_zeta*displacements[j*3+2] ;
				}
			}
			
			if(parent->isMoved())
				parent->getInverseJacobianMatrix(p_[i], Jinv) ;
			Vector llstrain(6) ;
			
			
			for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
			{
				double f_xi = vm.deval(parent->getEnrichmentFunction(j), XI, p_[i]) ;
				double f_eta = vm.deval(parent->getEnrichmentFunction(j), ETA, p_[i]) ;
				double f_zeta = vm.deval(parent->getEnrichmentFunction(j), ZETA, p_[i]) ;
				
				if(j*3 < enrichedDisplacements.size())
				{
					x_xi += f_xi*enrichedDisplacements[j*3] ;
					x_eta += f_eta*enrichedDisplacements[j*3] ;
					x_zeta += f_zeta*enrichedDisplacements[j*3] ;
				}
				
				if(j*3+1 < enrichedDisplacements.size())
				{
					y_xi += f_xi*enrichedDisplacements[j*3+1] ;
					y_eta += f_eta*enrichedDisplacements[j*3+1] ;
					y_zeta += f_zeta*enrichedDisplacements[j*3+1] ;
				}
				
				if(j*3+2 < enrichedDisplacements.size())
				{
					z_xi += f_xi*enrichedDisplacements[j*3+2] ;
					z_eta += f_eta*enrichedDisplacements[j*3+2] ;
					z_zeta += f_zeta*enrichedDisplacements[j*3+2] ;
				}
			}
			
			if(parent->isMoved())
				parent->getInverseJacobianMatrix(p_[i], Jinv) ;
			
			
			llstrain[0] = (x_xi)*Jinv[0][0] + (x_eta)*Jinv[0][1]  + (x_zeta)*Jinv[0][2];
			llstrain[1] = (y_xi)*Jinv[1][0] + (y_eta)*Jinv[1][1]  + (y_zeta)*Jinv[1][2];
			llstrain[2] = (z_xi)*Jinv[2][0] + (z_eta)*Jinv[2][1]  + (z_zeta)*Jinv[2][2]; 
			
			llstrain[3] = 0.5* ((y_xi  ) * Jinv[2][0] + 
			                    (y_eta ) * Jinv[2][1] +  
			                    (y_zeta) * Jinv[2][2] + 
			                    (z_xi  ) * Jinv[1][0] +
			                    (z_eta ) * Jinv[1][1] + 
			                    (z_zeta) * Jinv[1][2]);	
			
			llstrain[4] = 0.5* ((x_xi  ) * Jinv[2][0] + 
			                    (x_eta ) * Jinv[2][1] +  
			                    (x_zeta) * Jinv[2][2] + 
			                    (z_xi  ) * Jinv[0][0] + 
			                    (z_eta ) * Jinv[0][1] + 
			                    (z_zeta) * Jinv[0][2]);	
			
			llstrain[5] = 0.5* ((y_xi  ) * Jinv[0][0] + 
			                    (y_eta ) * Jinv[0][1] +  
			                    (y_zeta) * Jinv[0][2] + 
			                    (x_xi  ) * Jinv[1][0] + 
			                    (x_eta ) * Jinv[1][1] + 
			                    (x_zeta) * Jinv[1][2]);	
			
			
			llstrain=llstrain*cg ;
			
			lstrain[i*6 + 0] = llstrain[0] ;
			lstrain[i*6 + 1] = llstrain[1] ;
			lstrain[i*6 + 2] = llstrain[2] ;
			lstrain[i*6 + 3] = llstrain[3] ;
			lstrain[i*6 + 4] = llstrain[4] ;
			lstrain[i*6 + 5] = llstrain[5] ;
			
		}
		return lstrain ;
	}

	int dofs = 2 ;
	if(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
		dofs = 3 ;
	int sz = dofs+(dofs-1)*(dofs-1) - dofs%2;
	return Vector(0., sz*pts.size()) ;
}

// Vector ElementState::getNonEnrichedStress(const std::valarray<std::pair<Point, double> > &p) const
// {
// 	Vector pts(3*p.size()) ;
// 	for(size_t i = 0 ; i < p.size() ; i++)
// 	{
// 		Vector str = getNonEnrichedStress(&p[i].first, true) ;
// 		for(size_t j = 0 ; j < 3 ; j++)
// 		{
// 			pts[i*3+j] = str[j] ;
// 		}
// 	}
// 	
// 	return pts ;
// }

Vector ElementState::getStress(const std::valarray<std::pair<Point, double> > & p) const
{
	if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2)
	{
		Vector pts(3*p.size()) ;
		for(size_t i = 0 ; i < p.size() ; i++)
		{
			Vector str = getStress(p[i].first, true) ;
			for(size_t j = 0 ; j < 3 ; j++)
			{
				pts[i*3+j] = str[j] ;
			}
		}
		
		return pts ;
	}
	else if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3)
	{
		Vector pts(6*p.size()) ;
		for(size_t i = 0 ; i < p.size() ; i++)
		{
			Vector str = getStress(p[i].first, true) ;
			for(size_t j = 0 ; j < 6 ; j++)
			{
				pts[i*6+j] = str[j] ;
			}
		}
		
		return pts ;
	}

	int dofs = 2 ;
	if(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
		dofs = 3 ;
	int sz = dofs+(dofs-1)*(dofs-1) - dofs%2;
	return Vector(0., sz*p.size()) ;
}

Vector ElementState::getStrain(const Mu::PointArray & pts) const
{
	if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2)
	{
		if(parent->getBehaviour()->type == VOID_BEHAVIOUR)
		{
			return Vector(0., 3*pts.size()) ;
		}
		
		std::vector<Point> p_ ;
		for(size_t i = 0 ; i < pts.size() ; i++)
		{
			p_.push_back(parent->inLocalCoordinates(*pts[i])) ;
		}
		
		VirtualMachine vm ;
		Vector lstrain(3*pts.size()) ;
		
		Matrix Jinv(2,2) ; parent->getInverseJacobianMatrix(p_[0], Jinv) ;
		
		for(size_t i = 1 ; i < pts.size() ; i++)
		{
			
			double x_xi = 0;
			double x_eta = 0;
			double y_xi = 0;
			double y_eta = 0;
			
			for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				double f_xi = vm.deval(parent->getShapeFunction(j), XI, p_[i]) ;
				double f_eta = vm.deval(parent->getShapeFunction(j), ETA, p_[i]) ;
				
				if(j*2 < displacements.size())
				{
					x_xi += f_xi*displacements[j*2] ;
					x_eta += f_eta*displacements[j*2] ;
				}
				if(j*2+1 < displacements.size())
				{
					y_xi += f_xi*displacements[j*2+1] ;
					y_eta += f_eta*displacements[j*2+1] ;
				}
			}
			
			double delta_x_xi  = 0;
			double delta_x_eta  = 0;
			double delta_y_xi  = 0;
			double delta_y_eta  = 0;
			
			for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
			{
				double f_xi = vm.deval( parent->getEnrichmentFunction(j),XI,p_[i] ) ;
				double f_eta = vm.deval( parent->getEnrichmentFunction(j),ETA,p_[i]) ;
					
				if(j*2 < enrichedDisplacements.size())
				{
					delta_x_xi += f_xi*enrichedDisplacements[j*2] ;
					delta_x_eta += f_eta*enrichedDisplacements[j*2] ;
				}
				if(j*2+1 < enrichedDisplacements.size())
				{
					delta_y_xi += f_xi*enrichedDisplacements[j*2+1] ;
					delta_y_eta += f_eta*enrichedDisplacements[j*2+1] ;
				}
			}
			
			if(parent->isMoved())
				parent->getInverseJacobianMatrix(p_[i], Jinv) ;
			
			lstrain[i*3 + 0] = (x_xi+delta_x_xi)*Jinv[0][0] + (x_eta+delta_x_eta)*Jinv[0][1] ;
			lstrain[i*3 + 1] = (y_xi+delta_y_xi)*Jinv[1][0] + (y_eta+delta_y_eta)*Jinv[1][1] ;
			lstrain[i*3 + 2] = 0.5*((x_xi+delta_x_xi)*Jinv[1][0] + (x_eta+delta_x_eta)*Jinv[1][1]  + (y_xi+delta_y_xi)*Jinv[0][0] + (y_eta+delta_y_eta)*Jinv[0][1] );	
		}
		return lstrain ;
	}
	else if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3)
	{
		if(parent->getBehaviour()->type == VOID_BEHAVIOUR)
		{
			return Vector(0., 6*pts.size()) ;
		}
		std::vector<Point> p_ ;
		for(size_t i = 0 ; i < pts.size() ; i++)
		{
			p_.push_back(parent->inLocalCoordinates(*pts[i])) ;
		}
		
	// 	
		VirtualMachine vm ;
		
		Vector lstrain(6*pts.size()) ;
		
		Matrix Jinv(3,3) ; parent->getInverseJacobianMatrix(p_[0], Jinv) ;
		
		
		for(size_t i = 0 ; i < pts.size() ; i++)
		{
			
			double x_xi = 0;
			double x_eta = 0;
			double x_zeta = 0;
			double y_xi = 0;
			double y_eta = 0;
			double y_zeta = 0;
			double z_xi = 0;
			double z_eta = 0;
			double z_zeta = 0;
			
			for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				double f_xi = vm.deval(parent->getShapeFunction(j), XI, p_[i]) ;
				double f_eta = vm.deval(parent->getShapeFunction(j), ETA, p_[i]) ;
				double f_zeta = vm.deval(parent->getShapeFunction(j), ZETA, p_[i]) ;
				
				if(j*3 < displacements.size())
				{
					x_xi += f_xi*displacements[j*3] ;
					x_eta += f_eta*displacements[j*3] ;
					x_zeta += f_zeta*displacements[j*3] ;
				}
				if(j*3+1 < displacements.size())
				{
					y_xi += f_xi*displacements[j*3+1] ;
					y_eta += f_eta*displacements[j*3+1] ;
					y_zeta += f_zeta*displacements[j*3+1] ;
				}
				if(j*3+2 < displacements.size())
				{
					z_xi += f_xi*displacements[j*3+2] ;
					z_eta += f_eta*displacements[j*3+2] ;
					z_zeta += f_zeta*displacements[j*3+2] ;
				}
			}
			
			if(parent->isMoved())
				parent->getInverseJacobianMatrix(p_[i], Jinv) ;
			Vector llstrain(6) ;
			
			
			for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
			{
				double f_xi = vm.deval(parent->getEnrichmentFunction(j), XI, p_[i]) ;
				double f_eta = vm.deval(parent->getEnrichmentFunction(j), ETA, p_[i]) ;
				double f_zeta = vm.deval(parent->getEnrichmentFunction(j), ZETA, p_[i]) ;
				
				if(j*3 < enrichedDisplacements.size())
				{
					x_xi += f_xi*enrichedDisplacements[j*3] ;
					x_eta += f_eta*enrichedDisplacements[j*3] ;
					x_zeta += f_zeta*enrichedDisplacements[j*3] ;
				}
				
				if(j*3+1 < enrichedDisplacements.size())
				{
					y_xi += f_xi*enrichedDisplacements[j*3+1] ;
					y_eta += f_eta*enrichedDisplacements[j*3+1] ;
					y_zeta += f_zeta*enrichedDisplacements[j*3+1] ;
				}
				
				if(j*3+2 < enrichedDisplacements.size())
				{
					z_xi += f_xi*enrichedDisplacements[j*3+2] ;
					z_eta += f_eta*enrichedDisplacements[j*3+2] ;
					z_zeta += f_zeta*enrichedDisplacements[j*3+2] ;
				}
			}
			
			if(parent->isMoved())
				parent->getInverseJacobianMatrix(p_[i], Jinv) ;
			
			
			llstrain[0] = (x_xi)*Jinv[0][0] + (x_eta)*Jinv[0][1]  + (x_zeta)*Jinv[0][2];
			llstrain[1] = (y_xi)*Jinv[1][0] + (y_eta)*Jinv[1][1]  + (y_zeta)*Jinv[1][2];
			llstrain[2] = (z_xi)*Jinv[2][0] + (z_eta)*Jinv[2][1]  + (z_zeta)*Jinv[2][2]; 
			
			llstrain[3] = 0.5* ((y_xi  ) * Jinv[2][0] + 
			                    (y_eta ) * Jinv[2][1] +  
			                    (y_zeta) * Jinv[2][2] + 
			                    (z_xi  ) * Jinv[1][0] +
			                    (z_eta ) * Jinv[1][1] + 
			                    (z_zeta) * Jinv[1][2]);	
			
			llstrain[4] = 0.5* ((x_xi  ) * Jinv[2][0] + 
			                    (x_eta ) * Jinv[2][1] +  
			                    (x_zeta) * Jinv[2][2] + 
			                    (z_xi  ) * Jinv[0][0] + 
			                    (z_eta ) * Jinv[0][1] + 
			                    (z_zeta) * Jinv[0][2]);	
			
			llstrain[5] = 0.5* ((y_xi  ) * Jinv[0][0] + 
			                    (y_eta ) * Jinv[0][1] +  
			                    (y_zeta) * Jinv[0][2] + 
			                    (x_xi  ) * Jinv[1][0] + 
			                    (x_eta ) * Jinv[1][1] + 
			                    (x_zeta) * Jinv[1][2]);	
			
						
			lstrain[i*6 + 0] = llstrain[0] ;
			lstrain[i*6 + 1] = llstrain[1] ;
			lstrain[i*6 + 2] = llstrain[2] ;
			lstrain[i*6 + 3] = llstrain[3] ;
			lstrain[i*6 + 4] = llstrain[4] ;
			lstrain[i*6 + 5] = llstrain[5] ;
			
		}
		return lstrain ;
	}
	int dofs = 2 ;
	if(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
		dofs = 3 ;
	int sz = dofs+(dofs-1)*(dofs-1) - dofs%2;
	return Vector(0., sz*pts.size()) ;
}

Vector ElementState::getStrain( const std::valarray<std::pair<Point, double> > &p) const
{
	if (parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2)
	{
		Vector pts(3*p.size()) ;
		for(size_t i = 0 ; i < p.size() ; i++)
		{
			Vector str = getStrain(p[i].first, true) ;
			for(size_t j = 0 ; j < 3 ; j++)
			{
				pts[i*3+j] = str[j] ;
			}
		}
		
		return pts ;
	}
	else if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3)
	{
		Vector pts(6*p.size()) ;
		for(size_t i = 0 ; i < p.size() ; i++)
		{
			Vector str = getStrain(p[i].first) ;
			for(size_t j = 0 ; j < 6 ; j++)
			{
				pts[i*6+j] = str[j] ;
			}
		}
		
		return pts ;
	}
	int dofs = 2 ;
	if(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
		dofs = 3 ;
	int sz = dofs+(dofs-1)*(dofs-1) - dofs%2;
	return Vector(0., sz*p.size()) ;
}


const Vector & ElementState::getDisplacements() const
{
	return this->displacements ;
}

Vector & ElementState::getDisplacements()
{
	return this->displacements ;
}

std::vector<double> ElementState::getEnrichedInterpolatingFactors(const Point & p, bool local) const
{

	std::vector<double> ret(parent->getEnrichmentFunctions().size()) ;
	VirtualMachine vm ;
	Point p_ = p ;
	if(!local)
		p_ = parent->inLocalCoordinates(p) ;

	for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
	{
		ret[j] = vm.eval(parent->getEnrichmentFunction(j), p_) ;
	}

	return ret;
}

std::vector<double> ElementState::getNonEnrichedInterpolatingFactors(const Point & p, bool local) const
{

	std::vector<double> ret(parent->getBoundingPoints().size()) ;
	VirtualMachine vm ;
	Point p_(p) ;
	if(!local)
		p_ = parent->inLocalCoordinates(p) ;
	for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
	{
		ret[j] = vm.eval(parent->getShapeFunction(j), p_) ;
	}
	
	return ret;
}

std::vector<double> ElementState::getInterpolatingFactors(const Point & p, bool local) const
{

	std::vector<double> ret(parent->getBoundingPoints().size()+parent->getEnrichmentFunctions().size()) ;
	VirtualMachine vm ;
	Point p_ = p ;
	if(!local)
		p_ = parent->inLocalCoordinates(p) ;

	for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
	{
		ret[j] = vm.eval(parent->getShapeFunction(j), p_) ;
	}

	for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
	{
		ret[j+parent->getBoundingPoints().size()] = vm.eval(parent->getEnrichmentFunction(j), p_) ;
	}

	return ret;
}

Vector ElementState::getDisplacements(const Point & p, bool local, bool fast, const Vector * source) const
{
	if (fast)
	{
		if (parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2)
		{
			Vector ret(0., 2) ;
			Vector d(0., parent->getBoundingPoints().size()) ;
			double rr = 2.*parent->getRadius() ;
			for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				d[j] = dist(p, parent->getBoundingPoint(j)) ;
				if(d[j] < POINT_TOLERANCE)
				{
					ret[0] += (*source)[ parent->getBoundingPoint(j).id*2] ;
					ret[1] += (*source)[ parent->getBoundingPoint(j).id*2+1] ;
					return ret ;
				}
			}
			double sd = d.sum() ;
			for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				double w = (sd-d[j])/(sd*(parent->getBoundingPoints().size()-1)) ;
				ret[0] += w*(*source)[ parent->getBoundingPoint(j).id*2] ;
				ret[1] += w*(*source)[ parent->getBoundingPoint(j).id*2+1] ;
			}
			return ret ;
		}
		else if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3)
		{
			Vector ret(0., 3) ;
			Vector d(0., parent->getBoundingPoints().size()) ;
			double rr = 2.*parent->getRadius() ;
			double t = 0 ;
			
			for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				d[j] = dist(p, parent->getBoundingPoint(j)) ;
				if(d[j] < POINT_TOLERANCE)
				{
					ret[0] += (*source)[ parent->getBoundingPoint(j).id*3] ;
					ret[1] += (*source)[ parent->getBoundingPoint(j).id*3+1] ;
					ret[2] += (*source)[ parent->getBoundingPoint(j).id*3+2] ;
					return ret ;
				}
			}
			double sd = d.sum() ;
			for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				double w = (sd-d[j])/(sd*(parent->getBoundingPoints().size()-1)) ;
				ret[0] += w*(*source)[parent->getBoundingPoint(j).id*3] ;
				ret[1] += w*(*source)[parent->getBoundingPoint(j).id*3+1] ;
				ret[2] += w*(*source)[parent->getBoundingPoint(j).id*3+2] ;
			}
			return ret ;
		}
	}
	if (parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2)
	{
		Vector ret(0., 2) ;
		VirtualMachine vm ;
		Point p_ = p ;
		
		if(!local)
			p_ = parent->inLocalCoordinates(p) ;
		
		for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
		{
			double f =  vm.eval(parent->getShapeFunction(j), p_) ;
			ret[0] += f*displacements[j*2] ;
			ret[1] += f*displacements[j*2+1] ;
		}
		
		for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
		{
			double f = vm.eval(parent->getEnrichmentFunction(j), p_) ;
			ret[0] += f*enrichedDisplacements[j*2] ;
			ret[1] += f*enrichedDisplacements[j*2+1] ;
		}
	
		return ret;
	}
	else if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3)
	{
		Vector ret(0., 3) ;
		VirtualMachine vm ;
		Point p_ = p ;
		
		if(!local)
			p_ = parent->inLocalCoordinates(p) ;
		
		for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
		{
			double f =  vm.eval(parent->getShapeFunction(j), p_) ;
			ret[0] += f*displacements[j*3] ;
			ret[1] += f*displacements[j*3+1] ;
			ret[2] += f*displacements[j*3+2] ;
		}
		
		for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
		{
			double f = vm.eval(parent->getEnrichmentFunction(j), p_) ;
			ret[0] += f*enrichedDisplacements[j*3] ;
			ret[1] += f*enrichedDisplacements[j*3+1] ;
			ret[2] += f*enrichedDisplacements[j*3+2] ;
		}
		
		return ret;
	}
	int dofs = 2 ;
	if(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
		dofs = 3 ;
	int sz = dofs+(dofs-1)*(dofs-1) - dofs%2;
	return Vector(0., sz) ;
}

Vector ElementState::getDisplacements(const std::valarray<Point> & p) const
{
	if (parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2)
	{
		Function x ;
		Function y ;
		
		std::vector<Point> p_ ;
		for(size_t i = 0 ; i < p.size() ; i++)
		{
			p_.push_back(parent->inLocalCoordinates(p[i])) ;
		}
		
		for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
		{
			x += parent->getShapeFunction(j)*displacements[j*2] ;
			y += parent->getShapeFunction(j)*displacements[j*2+1] ;
		}
	
		VirtualMachine vm ;
		Vector ret(2*p.size()) ;
		for(size_t i = 0 ; i < p.size() ; i++)
		{
			double delta_x = 0;
			double delta_y = 0;
			
			for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
			{
				delta_x += vm.eval( parent->getEnrichmentFunction(j),p_[i])*enrichedDisplacements[j*2] ;
				delta_y += vm.eval( parent->getEnrichmentFunction(j),p_[i])*enrichedDisplacements[j*2+1] ;
			}
			ret[i*2+0] = vm.eval(x+delta_x, p_[i]) ;
			ret[i*2+1] = vm.eval(y+delta_y, p_[i]) ;
		}
		return ret;
	}
	else if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3)
	{
		Function x ;
		Function y ;
		Function z ;
		
		std::vector<Point> p_ ;
		for(size_t i = 0 ; i < p.size() ; i++)
		{
			p_.push_back(parent->inLocalCoordinates(p[i])) ;
		}
		
		for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
		{
			x += parent->getShapeFunction(j)*displacements[j*3] ;
			y += parent->getShapeFunction(j)*displacements[j*3+1] ;
			z += parent->getShapeFunction(j)*displacements[j*3+2] ;
		}
		
		VirtualMachine vm ;
		Vector ret(3*p.size()) ;
		for(size_t i = 0 ; i < p.size() ; i++)
		{
			double delta_x = 0;
			double delta_y = 0;
			double delta_z = 0;
			
			for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
			{
				delta_x += vm.eval( parent->getEnrichmentFunction(j),p_[i])*enrichedDisplacements[j*3] ;
				delta_y += vm.eval( parent->getEnrichmentFunction(j),p_[i])*enrichedDisplacements[j*3+1] ;
				delta_z += vm.eval( parent->getEnrichmentFunction(j),p_[i])*enrichedDisplacements[j*3+2] ;
			}
			ret[i*3+0] = vm.eval(x+delta_x, p_[i]) ;
			ret[i*3+1] = vm.eval(y+delta_y, p_[i]) ;
			ret[i*3+2] = vm.eval(z+delta_z, p_[i]) ;
		}
		return ret;
	}
	int dofs = 2 ;
	if(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
		dofs = 3 ;
	int sz = dofs+(dofs-1)*(dofs-1) - dofs%2;
	return Vector(0., sz*p.size()) ;
}


Vector ElementState::getDisplacements(const std::vector<std::pair<Point, double> > & p ) const
{
	if (parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2)
	{
		VirtualMachine vm ;
		Vector ret(0., 2*p.size()) ;
		
		for(size_t i = 0 ; i < p.size() ; i++)
		{
			for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				double f = vm.eval(parent->getShapeFunction(j),p[i].first) ;
				ret[i*2+0] += f*displacements[j*2] ;
				ret[i*2+1] += f*displacements[j*2+1] ;
			}
				
			for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
			{
				double f = vm.eval( parent->getEnrichmentFunction(j),p[i].first) ;
				ret[i*2+0] += f*enrichedDisplacements[j*2] ;
				ret[i*2+1] += f*enrichedDisplacements[j*2+1] ;
			}
	
		}
		return ret;
	}
	else if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3)
	{
		VirtualMachine vm ;
		Vector ret(0., 3*p.size()) ;
		
		for(size_t i = 0 ; i < p.size() ; i++)
		{
			for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				double f = vm.eval(parent->getShapeFunction(j),p[i].first) ;
				ret[i*3+0] += f*displacements[j*3] ;
				ret[i*3+1] += f*displacements[j*3+1] ;
				ret[i*3+2] += f*displacements[j*3+2] ;
			}
			
			for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
			{
				double f = vm.eval( parent->getEnrichmentFunction(j),p[i].first) ;
				ret[i*3+0] += f*enrichedDisplacements[j*3] ;
				ret[i*3+1] += f*enrichedDisplacements[j*3+1] ;
				ret[i*3+2] += f*enrichedDisplacements[j*3+2] ;
			}
			
		}
		return ret;
	}
	int dofs = 2 ;
	if(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
		dofs = 3 ;
	int sz = dofs+(dofs-1)*(dofs-1) - dofs%2;
	return Vector(0., sz*p.size()) ;
}


Vector ElementState::getPreviousDisplacements(const Point & p) const
{
	if (parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2)
	{
		Function x ;
		Function y ;
		
		Point p_ = parent->inLocalCoordinates(p) ;
		
		for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
		{
			x += parent->getShapeFunction(j)*previousDisplacements[j*2] ;
			y += parent->getShapeFunction(j)*previousDisplacements[j*2+1] ;
		}
		
		double delta_x = 0;
		double delta_y = 0;
		VirtualMachine vm ;
		for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
		{
			delta_x += vm.deval( parent->getEnrichmentFunction(j),XI,p_)*previousEnrichedDisplacements[j*2] ;
			delta_y += vm.deval( parent->getEnrichmentFunction(j),ETA,p_)*previousEnrichedDisplacements[j*2+1] ;
		}
		
		Vector ret(2) ;
		ret[0] = vm.eval(x + delta_x, p) ;
		ret[1] = vm.eval(y + delta_y, p) ;
		return ret;
	}
	else if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3)
	{
		Function x ;
		Function y ;
		Function z ;
		
		Point p_ = parent->inLocalCoordinates(p) ;
		
		for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
		{
			x += parent->getShapeFunction(j)*previousDisplacements[j*3] ;
			y += parent->getShapeFunction(j)*previousDisplacements[j*3+1] ;
			z += parent->getShapeFunction(j)*previousDisplacements[j*3+2] ;
		}
		
		double delta_x = 0;
		double delta_y = 0;
		double delta_z = 0;
		VirtualMachine vm ;
		for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
		{
			delta_x += vm.deval( parent->getEnrichmentFunction(j),XI,p_)*previousEnrichedDisplacements[j*3] ;
			delta_y += vm.deval( parent->getEnrichmentFunction(j),ETA,p_)*previousEnrichedDisplacements[j*3+1] ;
			delta_z += vm.deval( parent->getEnrichmentFunction(j),ZETA,p_)*previousEnrichedDisplacements[j*3+2] ;
		}
		
		Vector ret(2) ;
		ret[0] = vm.eval(x + delta_x, p) ;
		ret[1] = vm.eval(y + delta_y, p) ;
		ret[2] = vm.eval(z + delta_z, p) ;
		return ret;
	}
	int dofs = 2 ;
	if(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
		dofs = 3 ;
	int sz = dofs+(dofs-1)*(dofs-1) - dofs%2;
	return Vector(0., sz) ;
}


Vector ElementState::getPreviousDisplacements(const std::valarray<Point> & p) const
{
	if (parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2)
	{
		Function x ;
		Function y ;
		
		std::vector<Point> p_ ;
		for(size_t i = 0 ; i < p.size() ; i++)
		{
			p_.push_back(parent->inLocalCoordinates(p[i])) ;
		}
		
		for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
		{
			x += parent->getShapeFunction(j)*previousDisplacements[j*2] ;
			y += parent->getShapeFunction(j)*previousDisplacements[j*2+1] ;
		}
		
		
		VirtualMachine vm ;
		Vector ret(2*p.size()) ;
		for(size_t i = 0 ; i < p.size() ; i++)
		{
			double delta_x = 0;
			double delta_y = 0;
			
			for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
			{
				delta_x += vm.eval( parent->getEnrichmentFunction(j),p_[i])*previousEnrichedDisplacements[j*2] ;
				delta_y += vm.eval( parent->getEnrichmentFunction(j),p_[i])*previousEnrichedDisplacements[j*2+1] ;
			}
			ret[i*2+0] = vm.eval(x+delta_x, p[i]) ;
			ret[i*2+1] = vm.eval(y+delta_y, p[i]) ;
		}
		return ret;
	}
	else if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3)
	{
		Function x ;
		Function y ;
		Function z ;
		
		std::vector<Point> p_ ;
		for(size_t i = 0 ; i < p.size() ; i++)
		{
			p_.push_back(parent->inLocalCoordinates(p[i])) ;
		}
		
		for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
		{
			x += parent->getShapeFunction(j)*previousDisplacements[j*3] ;
			y += parent->getShapeFunction(j)*previousDisplacements[j*3+1] ;
			z += parent->getShapeFunction(j)*previousDisplacements[j*3+2] ;
		}
		
		
		VirtualMachine vm ;
		Vector ret(3*p.size()) ;
		for(size_t i = 0 ; i < p.size() ; i++)
		{
			double delta_x = 0;
			double delta_y = 0;
			double delta_z = 0;
			
			for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
			{
				delta_x += vm.eval( parent->getEnrichmentFunction(j),p_[i])*previousEnrichedDisplacements[j*3] ;
				delta_y += vm.eval( parent->getEnrichmentFunction(j),p_[i])*previousEnrichedDisplacements[j*3+1] ;
				delta_z += vm.eval( parent->getEnrichmentFunction(j),p_[i])*previousEnrichedDisplacements[j*3+2] ;
			}
			ret[i*3+0] = vm.eval(x+delta_x, p[i]) ;
			ret[i*3+1] = vm.eval(y+delta_y, p[i]) ;
			ret[i*3+2] = vm.eval(z+delta_z, p[i]) ;
		}
		return ret;
	}

	int dofs = 2 ;
	if(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
		dofs = 3 ;
	int sz = dofs+(dofs-1)*(dofs-1) - dofs%2;
	return Vector(0., sz*p.size()) ;
}

Vector ElementState::getPreviousPreviousDisplacements(const Point & p) const
{
	if (parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2)
	{
		Function x ;
		Function y ;
		
		Point p_ = parent->inLocalCoordinates(p) ;
		
		for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
		{
			x += parent->getShapeFunction(j)*previousPreviousDisplacements[j*2] ;
			y += parent->getShapeFunction(j)*previousPreviousDisplacements[j*2+1] ;
		}
		
		double delta_x = 0;
		double delta_y = 0;
		VirtualMachine vm ;
		for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
		{
			delta_x += vm.deval( parent->getEnrichmentFunction(j),XI,p_)*previousPreviousEnrichedDisplacements[j*2] ;
			delta_y += vm.deval( parent->getEnrichmentFunction(j),ETA,p_)*previousPreviousEnrichedDisplacements[j*2+1] ;
		}
		
		Vector ret(2) ;
		ret[0] = vm.eval(x + delta_x, p) ;
		ret[1] = vm.eval(y + delta_y, p) ;
		return ret;
	}
	else if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3)
	{
		Function x ;
		Function y ;
		Function z ;
		
		Point p_ = parent->inLocalCoordinates(p) ;
		
		for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
		{
			x += parent->getShapeFunction(j)*previousPreviousDisplacements[j*3] ;
			y += parent->getShapeFunction(j)*previousPreviousDisplacements[j*3+1] ;
			z += parent->getShapeFunction(j)*previousPreviousDisplacements[j*3+2] ;
		}
		
		double delta_x = 0;
		double delta_y = 0;
		double delta_z = 0;
		VirtualMachine vm ;
		for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
		{
			delta_x += vm.deval( parent->getEnrichmentFunction(j),XI,p_)*previousPreviousEnrichedDisplacements[j*3] ;
			delta_y += vm.deval( parent->getEnrichmentFunction(j),ETA,p_)*previousPreviousEnrichedDisplacements[j*3+1] ;
			delta_z += vm.deval( parent->getEnrichmentFunction(j),ZETA,p_)*previousPreviousEnrichedDisplacements[j*3+2] ;
		}
		
		Vector ret(3) ;
		ret[0] = vm.eval(x + delta_x, p) ;
		ret[1] = vm.eval(y + delta_y, p) ;
		ret[2] = vm.eval(z + delta_z, p) ;
		return ret;
	}

	int dofs = 2 ;
	if(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
		dofs = 3 ;
	int sz = dofs+(dofs-1)*(dofs-1) - dofs%2;
	return Vector(0., sz) ;
}

Vector ElementState::getPreviousPreviousDisplacements(const std::valarray<Point> & p) const
{
	if (parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2)
	{
		Function x ;
		Function y ;
		
		std::vector<Point> p_ ;
		for(size_t i = 0 ; i < p.size() ; i++)
		{
			p_.push_back(parent->inLocalCoordinates(p[i])) ;
		}
		
		for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
		{
			x += parent->getShapeFunction(j)*previousPreviousDisplacements[j*2] ;
			y += parent->getShapeFunction(j)*previousPreviousDisplacements[j*2+1] ;
		}
		
		
		VirtualMachine vm ;
		Vector ret(2*p.size()) ;
		for(size_t i = 0 ; i < p.size() ; i++)
		{
			double delta_x = 0;
			double delta_y = 0;
			
			for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
			{
				delta_x += vm.eval( parent->getEnrichmentFunction(j),p_[i])*previousPreviousEnrichedDisplacements[j*2] ;
				delta_y += vm.eval( parent->getEnrichmentFunction(j),p_[i])*previousPreviousEnrichedDisplacements[j*2+1] ;
			}
			ret[i*2+0] = vm.eval(x+delta_x, p[i]) ;
			ret[i*2+1] = vm.eval(y+delta_y, p[i]) ;
		}
		return ret;
	}
	else if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3)
	{
		Function x ;
		Function y ;
		Function z ;
		
		std::vector<Point> p_ ;
		for(size_t i = 0 ; i < p.size() ; i++)
		{
			p_.push_back(parent->inLocalCoordinates(p[i])) ;
		}
		
		for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
		{
			x += parent->getShapeFunction(j)*previousPreviousDisplacements[j*3] ;
			y += parent->getShapeFunction(j)*previousPreviousDisplacements[j*3+1] ;
			z += parent->getShapeFunction(j)*previousPreviousDisplacements[j*3+2] ;
		}
		
		
		VirtualMachine vm ;
		Vector ret(3*p.size()) ;
		for(size_t i = 0 ; i < p.size() ; i++)
		{
			double delta_x = 0;
			double delta_y = 0;
			double delta_z = 0;
			
			for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
			{
				delta_x += vm.eval( parent->getEnrichmentFunction(j),p_[i])*previousPreviousEnrichedDisplacements[j*3] ;
				delta_y += vm.eval( parent->getEnrichmentFunction(j),p_[i])*previousPreviousEnrichedDisplacements[j*3+1] ;
				delta_z += vm.eval( parent->getEnrichmentFunction(j),p_[i])*previousPreviousEnrichedDisplacements[j*3+2] ;
			}
			ret[i*3+0] = vm.eval(x+delta_x, p[i]) ;
			ret[i*3+1] = vm.eval(y+delta_y, p[i]) ;
			ret[i*3+1] = vm.eval(z+delta_z, p[i]) ;
		}
		return ret;
	}

	int dofs = 2 ;
	if(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
		dofs = 3 ;
	int sz = dofs+(dofs-1)*(dofs-1) - dofs%2;
	return Vector(0., sz*p.size()) ;
}



std::pair<Vector, Vector > ElementState::getStressAndStrain(const Mu::PointArray & pts) const
{
	
	if (parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2)
	{
		if(parent->getBehaviour()->type == VOID_BEHAVIOUR)
		{
			std::cout << "void triangle" << std::endl ;
			return std::pair<Vector, Vector>(Vector(0., 3*pts.size()), Vector(0., 3*pts.size())) ;
		}
		
		std::vector<Point> p_ ;
		for(size_t i = 0 ; i < pts.size() ; i++)
		{
			p_.push_back(parent->inLocalCoordinates(*pts[i])) ;
		}
		
		VirtualMachine vm ;
		
		Vector lstrain(3*pts.size()) ;
		Vector lstress(3*pts.size()) ;
		
		size_t ndofs = parent->getBehaviour()->getNumberOfDegreesOfFreedom() ;
		int offset = ndofs-1 ;
		
		if(ndofs == 2)
		{
			for(size_t i = 0 ; i < pts.size() ; i++)
			{
				bool hasTimePlanes = (parent->timePlanes()>1) ;
				Matrix Jinv(2+hasTimePlanes,2+hasTimePlanes) ; parent->getInverseJacobianMatrix(p_[i], Jinv) ;
				Matrix cg (parent->getBehaviour()->getTensor(p_[i])) ;
				double x_xi = 0;
				double x_eta = 0;
				double y_xi = 0;
				double y_eta = 0;
	
				for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
				{
					double f_xi = vm.deval(parent->getShapeFunction(j),XI, p_[i]) ;
					double f_eta = vm.deval(parent->getShapeFunction(j),ETA, p_[i]) ;
			
					x_xi += f_xi*displacements[j*ndofs] ;
					x_eta += f_eta*displacements[j*ndofs] ;
		
					y_xi += f_xi*displacements[j*ndofs+offset] ;
					y_eta += f_eta*displacements[j*ndofs+offset] ;
				}
				
				Vector llstrain(3) ;
				
				for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
				{
					double f_xi = vm.deval( parent->getEnrichmentFunction(j),XI,p_[i] ) ;
					double f_eta = vm.deval( parent->getEnrichmentFunction(j),ETA,p_[i]) ;
					
					x_xi += f_xi*enrichedDisplacements[j*ndofs] ;
					x_eta += f_eta*enrichedDisplacements[j*ndofs] ;
					
					y_xi += f_xi*enrichedDisplacements[j*ndofs+offset] ;
					y_eta += f_eta*enrichedDisplacements[j*ndofs+offset] ;
				}
				
				llstrain[0] = (x_xi)*Jinv[0][0] + (x_eta)*Jinv[0][1] ;
				llstrain[1] = (y_xi)*Jinv[1][0] + (y_eta)*Jinv[1][1] ;
				llstrain[2] = 0.5*((x_xi)*Jinv[1][0] + (x_eta)*Jinv[1][1]  + (y_xi)*Jinv[0][0] + (y_eta)*Jinv[0][1] );	
				
				lstrain[i*3 + 0] = llstrain[0] ;
				lstrain[i*3 + 1] = llstrain[1] ;
				lstrain[i*3 + 2] = llstrain[2] ;
		
				llstrain=llstrain*cg ;
				
				lstress[i*3 + 0] = llstrain[0] ;
				lstress[i*3 + 1] = llstrain[1] ;
				lstress[i*3 + 2] = llstrain[2] ;
				
			}
		}
		std::pair<Vector, Vector> ret(lstress, lstrain) ;
		
		return ret ;
	}
	else if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3)
	{
		if(parent->getBehaviour()->type == VOID_BEHAVIOUR)
		{
			return std::pair<Vector, Vector>(Vector(0., 6*pts.size()), Vector(0., 6*pts.size())) ;
		}
		
		std::vector<Point> p_ ;
		for(size_t i = 0 ; i < pts.size() ; i++)
		{
			p_.push_back(parent->inLocalCoordinates(*pts[i])) ;
		}
		
		VirtualMachine vm ;
		
		Vector lstrain(6*pts.size()) ;
		Vector lstress(6*pts.size()) ;
		
		size_t ndofs = parent->getBehaviour()->getNumberOfDegreesOfFreedom() ;

		if(ndofs == 3)
		{
			for(size_t i = 0 ; i < pts.size() ; i++)
			{
				Matrix Jinv(3,3) ; parent->getInverseJacobianMatrix(Point(.25, .25, .25), Jinv) ;
				Matrix cg (parent->getBehaviour()->getTensor(p_[i])) ;
				Vector llstrain(0., 6) ;
				double x_xi = 0;
				double x_eta = 0;
				double x_zeta = 0;
				double y_xi = 0;
				double y_eta = 0;
				double y_zeta = 0;
				double z_xi = 0;
				double z_eta = 0;
				double z_zeta = 0;

				for(size_t j = 0 ; j < parent->getShapeFunctions().size() ; j++)
				{
					double f_xi = vm.deval(parent->getShapeFunction(j),XI, p_[i]) ;
					double f_eta = vm.deval(parent->getShapeFunction(j),ETA, p_[i]) ;
					double f_zeta = vm.deval(parent->getShapeFunction(j),ZETA, p_[i]) ;
					double x = displacements[j*3] ;
					double y = displacements[j*3+1] ;
					double z = displacements[j*3+2] ;
					
					x_xi   += f_xi   * x ;
					x_eta  += f_eta  * x ;
					x_zeta += f_zeta * x ;
					y_xi   += f_xi   * y ;
					y_eta  += f_eta  * y ;
					y_zeta += f_zeta * y ;
					z_xi   += f_xi   * z ;
					z_eta  += f_eta  * z ;
					z_zeta += f_zeta * z ;
				}
				
				
				for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
				{
					double f_xi = vm.deval( parent->getEnrichmentFunction(j),XI,p_[i] ) ;
					double f_eta = vm.deval( parent->getEnrichmentFunction(j),ETA,p_[i]) ;
					double f_zeta = vm.deval( parent->getEnrichmentFunction(j),ZETA,p_[i]) ;
					double x = enrichedDisplacements[j*3] ;
					double y = enrichedDisplacements[j*3+1] ;
					double z = enrichedDisplacements[j*3+2] ;
					
					x_xi   += f_xi   * x ;
					x_eta  += f_eta  * x ;
					x_zeta += f_zeta * x ;
					y_xi   += f_xi   * y ;
					y_eta  += f_eta  * y ;
					y_zeta += f_zeta * y ;
					z_xi   += f_xi   * z ;
					z_eta  += f_eta  * z ;
					z_zeta += f_zeta * z ;
				}
// 				std::cout << x_xi << ", " << x_eta << ", " << x_zeta << std::endl ;
// 				std::cout << y_xi << ", " << y_eta << ", " << y_zeta << std::endl ;
// 				std::cout << z_xi << ", " << z_eta << ", " << z_zeta << std::endl ;
				llstrain[0] = (x_xi)*Jinv[0][0] + (x_eta)*Jinv[0][1]  + (x_zeta)*Jinv[0][2];
				llstrain[1] = (y_xi)*Jinv[1][0] + (y_eta)*Jinv[1][1]  + (y_zeta)*Jinv[1][2];
				llstrain[2] = (z_xi)*Jinv[2][0] + (z_eta)*Jinv[2][1]  + (z_zeta)*Jinv[2][2]; 
				
				llstrain[5] = 0.5* ((y_xi ) * Jinv[2][0] + 
				                   (y_eta ) * Jinv[2][1] +  
				                   (y_zeta) * Jinv[2][2] + 
				                   (z_xi  ) * Jinv[1][0] +
				                   (z_eta ) * Jinv[1][1] + 
				                   (z_zeta) * Jinv[1][2]);	
				
				llstrain[4] = 0.5* ((x_xi ) * Jinv[2][0] + 
				                   (x_eta ) * Jinv[2][1] +  
				                   (x_zeta) * Jinv[2][2] + 
				                   (z_xi  ) * Jinv[0][0] + 
				                   (z_eta ) * Jinv[0][1] + 
				                   (z_zeta) * Jinv[0][2]);	
				
				llstrain[3] = 0.5* ((y_xi)  * Jinv[0][0] + 
				                   (y_eta)  * Jinv[0][1] +  
				                   (y_zeta) * Jinv[0][2] + 
				                   (x_xi)   * Jinv[1][0] + 
				                   (x_eta)  * Jinv[1][1] + 
				                   (x_zeta) * Jinv[1][2]);	
				
				lstrain[i*6 + 0] = llstrain[0] ;
				lstrain[i*6 + 1] = llstrain[1] ;
				lstrain[i*6 + 2] = llstrain[2] ;
				lstrain[i*6 + 3] = llstrain[3] ;
				lstrain[i*6 + 4] = llstrain[4] ;
				lstrain[i*6 + 5] = llstrain[5] ;
				
				llstrain=llstrain*cg ;
				
				lstress[i*6 + 0] = llstrain[0] ;
				lstress[i*6 + 1] = llstrain[1] ;
				lstress[i*6 + 2] = llstrain[2] ;
				lstress[i*6 + 3] = llstrain[3] ;
				lstress[i*6 + 4] = llstrain[4] ;
				lstress[i*6 + 5] = llstrain[5] ;
			}
		}
		std::pair<Vector, Vector> ret ;
		ret.first.resize(lstress.size()) ;
		ret.second.resize(lstrain.size()) ;
		ret.first = lstress ;
		ret.second = lstrain ;
		return ret ;
	}
	int dofs = 2 ;
	if(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
		dofs = 3 ;
	int sz = dofs+(dofs-1)*(dofs-1) - dofs%2;
	return std::pair<Vector, Vector>(Vector(0., sz*pts.size()), Vector(0., sz*pts.size())) ;
}

std::pair<Vector, Vector > ElementState::getStressAndStrain( std::valarray<std::pair<Point, double> > & p) const
{
	Mu::PointArray pts(p.size()) ;
	for(size_t i = 0 ; i < p.size() ; i++)
	{
		pts[i] = &p[i].first ;
	}
	
	return getStressAndStrain(pts) ;
}


Vector ElementState::getDeltaStrain(const Point & p) const
{
	if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2)
	{
		Point p_ = parent->inLocalCoordinates(p) ;
		VirtualMachine vm ;
	
		double x_xi = 0;
		double x_eta = 0;
		double y_xi = 0;
		double y_eta = 0;
		for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
		{
			double f_xi = vm.deval(parent->getShapeFunction(j), XI, p_) ;
			double f_eta = vm.deval(parent->getShapeFunction(j), ETA, p_) ;
			x_xi += f_xi*(displacements[j*2]-previousDisplacements[j*2]) ;
			x_eta += f_eta*(displacements[j*2]-previousDisplacements[j*2]) ;
			y_xi += f_xi*(displacements[j*2+1]-previousDisplacements[j*2+1]) ;
			y_eta += f_eta*(displacements[j*2+1]-previousDisplacements[j*2+1]) ;
		}
		
		
		Vector lstrain(3) ;
		
		Matrix cg (parent->getBehaviour()->getTensor(p_)) ;
		
		double delta_x_xi  = 0;
		double delta_x_eta  = 0;
		double delta_y_xi  = 0;
		double delta_y_eta  = 0;
		
		for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
		{
			delta_x_xi += vm.deval(parent->getEnrichmentFunction(j),XI,p_ )*(enrichedDisplacements[j*2] - previousEnrichedDisplacements[j*2]) ;
			delta_x_eta += vm.deval(parent->getEnrichmentFunction(j),ETA,p_)*(enrichedDisplacements[j*2] - previousEnrichedDisplacements[j*2]) ;
			delta_y_xi += vm.deval(parent->getEnrichmentFunction(j),XI,p_)*(enrichedDisplacements[j*2+1] - previousEnrichedDisplacements[j*2+1]) ;
			delta_y_eta += vm.deval(parent->getEnrichmentFunction(j),ETA,p_)*(enrichedDisplacements[j*2+1] - previousEnrichedDisplacements[j*2+1]) ;
		}
		
		Matrix Jinv(2,2) ; parent->getInverseJacobianMatrix(p_, Jinv) ;
		lstrain[0] = (x_xi+delta_x_xi)*Jinv[0][0] + (x_eta+delta_x_eta)*Jinv[0][1] ;
		lstrain[1] = (y_xi+delta_y_xi)*Jinv[1][0] + (y_eta+delta_y_eta)*Jinv[1][1] ;
		lstrain[2] = 0.5*((x_xi+delta_x_xi)*Jinv[1][0] + (x_eta+delta_x_eta)*Jinv[1][1]  + (y_xi+delta_y_xi)*Jinv[0][0] + (y_eta+delta_y_eta)*Jinv[0][1] );	
		
		return lstrain ;
	}
	else if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3)
	{
		Point p_ = parent->inLocalCoordinates(p) ;
		VirtualMachine vm ;
		
		double x_xi = 0;
		double x_eta = 0;
		double x_zeta = 0;
		double y_xi = 0;
		double y_eta = 0;
		double y_zeta = 0;
		double z_xi = 0;
		double z_eta = 0;
		double z_zeta = 0;
		for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
		{
			double f_xi = vm.deval(parent->getShapeFunction(j), XI, p_) ;
			double f_eta = vm.deval(parent->getShapeFunction(j), ETA, p_) ;
			double f_zeta = vm.deval(parent->getShapeFunction(j), ZETA, p_) ;
			
			x_xi += f_xi*(displacements[j*3]-previousDisplacements[j*3]) ;
			x_eta += f_eta*(displacements[j*3]-previousDisplacements[j*3]) ;
			x_zeta += f_zeta*(displacements[j*3]-previousDisplacements[j*3]) ;
			y_xi += f_xi*(displacements[j*3+1]-previousDisplacements[j*3+1]) ;
			y_eta += f_eta*(displacements[j*3+1]-previousDisplacements[j*3+1]) ;
			y_zeta += f_zeta*(displacements[j*3+1]-previousDisplacements[j*3+1]) ;
			z_xi += f_xi*(displacements[j*3+2]-previousDisplacements[j*3+2]) ;
			z_eta += f_eta*(displacements[j*3+2]-previousDisplacements[j*3+2]) ;
			z_zeta += f_zeta*(displacements[j*3+2]-previousDisplacements[j*3+2]) ;
		}
		
		
		Vector lstrain(6) ;
		
		Matrix cg(parent->getBehaviour()->getTensor(p_)) ;
		
		
		for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
		{
		
			double f_xi = vm.deval(parent->getEnrichmentFunction(j),XI,p_ ) ;
			double f_eta = vm.deval(parent->getEnrichmentFunction(j),ETA,p_) ;
			double f_zeta = vm.deval(parent->getEnrichmentFunction(j),ZETA,p_) ;
			
			x_xi += f_xi*(enrichedDisplacements[j*3] - previousEnrichedDisplacements[j*3]) ;
			x_eta += f_eta*(enrichedDisplacements[j*3] - previousEnrichedDisplacements[j*3]) ;
			x_zeta += f_zeta*(enrichedDisplacements[j*3] - previousEnrichedDisplacements[j*3]) ;
			y_xi += f_xi*(enrichedDisplacements[j*3+1] - previousEnrichedDisplacements[j*3+1]) ;
			y_eta += f_eta*(enrichedDisplacements[j*3+1] - previousEnrichedDisplacements[j*3+1]) ;
			y_zeta += f_zeta*(enrichedDisplacements[j*3+1] - previousEnrichedDisplacements[j*3+1]) ;
			z_xi += f_xi*(enrichedDisplacements[j*3+2] - previousEnrichedDisplacements[j*3+2]) ;
			z_eta += f_eta*(enrichedDisplacements[j*3+2] - previousEnrichedDisplacements[j*3+2]) ;
			z_zeta += f_zeta*(enrichedDisplacements[j*3+2] - previousEnrichedDisplacements[j*3+2]) ;
		}
		
		Matrix Jinv(3,3)  ; parent->getInverseJacobianMatrix(p_, Jinv) ;
		lstrain[0] = (x_xi)*Jinv[0][0] + (x_eta)*Jinv[0][1]  + (x_zeta)*Jinv[0][2];
		lstrain[1] = (y_xi)*Jinv[1][0] + (y_eta)*Jinv[1][1]  + (y_zeta)*Jinv[1][2];
		lstrain[2] = (z_xi)*Jinv[2][0] + (z_eta)*Jinv[2][1]  + (z_zeta)*Jinv[2][2];
		lstrain[3] = 0.5*((x_xi)*Jinv[1][0] + (x_eta)*Jinv[1][1]  +  (x_zeta)*Jinv[1][2]
		                  + (y_xi)*Jinv[0][0] + (y_eta)*Jinv[0][1]  + (y_zeta)*Jinv[0][2]);	
		lstrain[4] = 0.5*((x_xi)*Jinv[2][0] + (x_eta)*Jinv[2][1]  +  (x_zeta)*Jinv[2][2]
		                  + (z_xi)*Jinv[0][0] + (z_eta)*Jinv[0][1]  + (z_zeta)*Jinv[0][2]);	
		lstrain[5] = 0.5*((y_xi)*Jinv[2][0] + (y_eta)*Jinv[2][1]  +  (y_zeta)*Jinv[2][2]
		                  + (z_xi)*Jinv[1][0] + (z_eta)*Jinv[1][1]  + (z_zeta)*Jinv[1][2]);	
		
		return lstrain ;
	}
	
	return Vector(0) ;
}

Vector ElementState::getDeltaStrain(const std::pair<Point, double> & p) const
{
	return getStrain(p.first) ;
}

Vector ElementState::getDeltaStress(const Point & p) const
{
	Point p_ = parent->inLocalCoordinates(p) ;
	Function x_xi ;
	Function x_eta ;
	Function y_xi ;
	Function y_eta ;
	
	for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
	{
		x_xi += parent->getShapeFunction(j).d(XI)*(displacements[j*2]-previousDisplacements[j*2]) ;
		x_eta += parent->getShapeFunction(j).d(ETA)*(displacements[j*2]-previousDisplacements[j*2]) ;
		y_xi += parent->getShapeFunction(j).d(XI)*(displacements[j*2+1]-previousDisplacements[j*2+1]) ;
		y_eta += parent->getShapeFunction(j).d(ETA)*(displacements[j*2+1]-previousDisplacements[j*2+1]) ;
	}
	
	
	VirtualMachine vm ;
	Vector lstrain(3) ;
	
	double delta_x_xi  = 0;
	double delta_x_eta  = 0;
	double delta_y_xi  = 0;
	double delta_y_eta  = 0;
	
	for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
	{
		delta_x_xi += vm.deval( parent->getEnrichmentFunction(j),XI,p_ )*(enrichedDisplacements[j*2] - previousEnrichedDisplacements[j*2]) ;
		delta_x_eta += vm.deval( parent->getEnrichmentFunction(j),ETA,p_)*(enrichedDisplacements[j*2] - previousEnrichedDisplacements[j*2]) ;
		delta_y_xi += vm.deval( parent->getEnrichmentFunction(j),XI,p_)*(enrichedDisplacements[j*2+1] - previousEnrichedDisplacements[j*2+1]);
		delta_y_eta += vm.deval( parent->getEnrichmentFunction(j),ETA,p_)*(enrichedDisplacements[j*2+1] - previousEnrichedDisplacements[j*2+1]) ;
	}
	
	Matrix Jinv(2,2)  ; parent->getInverseJacobianMatrix(p_, Jinv) ;
	lstrain[0] = vm.eval(x_xi+delta_x_xi, p_)*Jinv[0][0] + vm.eval(x_eta+delta_x_eta, p_)*Jinv[0][1] ;
	lstrain[1] = vm.eval(y_xi+delta_y_xi, p_)*Jinv[1][0] + vm.eval(y_eta+delta_y_eta, p_)*Jinv[1][1] ;
	lstrain[2] = 0.5*(vm.eval(x_xi+delta_x_xi, p_)*Jinv[1][0] + vm.eval(x_eta+delta_x_eta, p_)*Jinv[1][1]  + vm.eval(y_xi+delta_y_xi, p_)*Jinv[0][0] + vm.eval(y_eta+delta_y_eta, p_)*Jinv[0][1] );	
	
	Matrix cg (parent->getBehaviour()->getTensor(p_)) ;
	
	return lstrain*cg ;
}

Vector ElementState::getDeltaStress(const std::pair<Point, double> & p) const
{
	return getDeltaStress(p.first) ;
}

Vector ElementState::getDeltaStress(const Mu::PointArray & pts) const
{
	if(parent->getBehaviour()->type == VOID_BEHAVIOUR)
	{
		return Vector(0., 3*pts.size()) ;
	}
	std::vector<Point> p_ ;
	for(size_t i = 0 ; i < pts.size() ; i++)
	{
		p_.push_back(parent->inLocalCoordinates(*pts[i])) ;
	}
	Function x_xi ;
	Function x_eta ;
	Function y_xi ;
	Function y_eta ;
	
	for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
	{
		x_xi += parent->getShapeFunction(j).d(XI)*(displacements[j*2] - previousDisplacements[j*2]) ;
		x_eta += parent->getShapeFunction(j).d(ETA)*(displacements[j*2] - previousDisplacements[j*2]) ;
		y_xi += parent->getShapeFunction(j).d(XI)*(displacements[j*2+1] - previousDisplacements[j*2+1]) ;
		y_eta += parent->getShapeFunction(j).d(ETA)*(displacements[j*2+1] - previousDisplacements[j*2+1]) ;
	}
	
	VirtualMachine vm ;
	
	Vector lstrain(3*pts.size()) ;
	
	Matrix Jinv(2,2)  ; parent->getInverseJacobianMatrix(p_[0], Jinv) ;
	
	for(size_t i = 0 ; i < pts.size() ; i++)
	{
		if(parent->isMoved())
			parent->getInverseJacobianMatrix(p_[i], Jinv) ;
		Vector llstrain(3) ;
		Matrix cg (parent->getBehaviour()->getTensor(p_[i])) ;
		
		double delta_x_xi  = 0;
		double delta_x_eta  = 0;
		double delta_y_xi  = 0;
		double delta_y_eta  = 0;
		
		for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
		{
			delta_x_xi += vm.deval( parent->getEnrichmentFunction(j),XI,p_[i] )*(enrichedDisplacements[j*2] - previousEnrichedDisplacements[j*2]);
			delta_x_eta += vm.deval( parent->getEnrichmentFunction(j),ETA,p_[i])*(enrichedDisplacements[j*2] - previousEnrichedDisplacements[j*2]);
			delta_y_xi += vm.deval( parent->getEnrichmentFunction(j),XI,p_[i])*(enrichedDisplacements[j*2+1] - previousEnrichedDisplacements[j*2+1]);
			delta_y_eta += vm.deval( parent->getEnrichmentFunction(j),ETA,p_[i])*(enrichedDisplacements[j*2+1] - previousEnrichedDisplacements[j*2+1]);
		}
		
		if(parent->isMoved())
			parent->getInverseJacobianMatrix(p_[i], Jinv) ;
		
		llstrain[0] = vm.eval(x_xi+delta_x_xi, p_[i])*Jinv[0][0] + vm.eval(x_eta+delta_x_eta, p_[i])*Jinv[0][1] ;
		llstrain[1] = vm.eval(y_xi+delta_y_xi, p_[i])*Jinv[1][0] + vm.eval(y_eta+delta_y_eta, p_[i])*Jinv[1][1] ;
		llstrain[2] = 0.5*(vm.eval(x_xi+delta_x_xi, p_[i])*Jinv[1][0] + vm.eval(x_eta+delta_x_eta, p_[i])*Jinv[1][1]  + vm.eval(y_xi+delta_y_xi, p_[i])*Jinv[0][0] + vm.eval(y_eta+delta_y_eta, p_[i])*Jinv[0][1] );	
		
		llstrain=llstrain*cg ;
		
		lstrain[i*3 + 0] = llstrain[0] ;
		lstrain[i*3 + 1] = llstrain[1] ;
		lstrain[i*3 + 2] = llstrain[2] ;
		
	}
	return lstrain ;
}

Vector ElementState::getDeltaStress(const std::valarray<std::pair<Point, double> > & p) const
{
	Vector pts(3*p.size()) ;
	for(size_t i = 0 ; i < p.size() ; i++)
	{
		Vector str = getDeltaStress(p[i].first) ;
		for(size_t j = 0 ; j < 3 ; j++)
		{
			pts[i*3+j] = str[j] ;
		}
	}
	
	return pts ;
}

Vector ElementState::getDeltaStrain(const Mu::PointArray & pts) const
{
	if(parent->getBehaviour()->type == VOID_BEHAVIOUR)
	{
		return Vector(0., 3*pts.size()) ;
	}
	
	std::vector<Point> p_ ;
	for(size_t i = 0 ; i < pts.size() ; i++)
	{
		p_.push_back(parent->inLocalCoordinates(*pts[i])) ;
	}
	Function x_xi ;
	Function x_eta ;
	Function y_xi ;
	Function y_eta ;
	
	for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
	{
		x_xi += parent->getShapeFunction(j).d(XI)*(displacements[j*2]- previousDisplacements[j*2]) ;
		x_eta += parent->getShapeFunction(j).d(ETA)*(displacements[j*2]- previousDisplacements[j*2]) ;
		y_xi += parent->getShapeFunction(j).d(XI)*(displacements[j*2+1]- previousDisplacements[j*2+1]) ;
		y_eta += parent->getShapeFunction(j).d(ETA)*(displacements[j*2+1]- previousDisplacements[j*2+1]) ;
	}
	
	VirtualMachine vm ;
	Vector lstrain(3*pts.size()) ;
	
	Matrix Jinv(2,2)  ; parent->getInverseJacobianMatrix(p_[0], Jinv) ;
	lstrain[0] = vm.eval(x_xi, p_[0])*Jinv[0][0] + vm.eval(x_eta, p_[0])*Jinv[0][1] ;
	lstrain[1] = vm.eval(y_xi, p_[0])*Jinv[1][0] + vm.eval(y_eta, p_[0])*Jinv[1][1] ;
	lstrain[2] = 0.5*(vm.eval(x_xi, p_[0])*Jinv[1][0] + vm.eval(x_eta, p_[0])*Jinv[1][1]  + vm.eval(y_xi, p_[0])*Jinv[0][0] + vm.eval(y_eta, p_[0])*Jinv[0][1] );	
	
	for(size_t i = 1 ; i < pts.size() ; i++)
	{
		
		double delta_x_xi  = 0;
		double delta_x_eta  = 0;
		double delta_y_xi  = 0;
		double delta_y_eta  = 0;
		
		for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
		{
			delta_x_xi += vm.deval( parent->getEnrichmentFunction(j),XI,p_[i] )*(enrichedDisplacements[j*2] - previousEnrichedDisplacements[j*2]);
				delta_x_eta += vm.deval( parent->getEnrichmentFunction(j),ETA,p_[i])*(enrichedDisplacements[j*2] - previousEnrichedDisplacements[j*2]);
				delta_y_xi += vm.deval( parent->getEnrichmentFunction(j),XI,p_[i])*(enrichedDisplacements[j*2+1] - previousEnrichedDisplacements[j*2+1]);
				delta_y_eta += vm.deval( parent->getEnrichmentFunction(j),ETA,p_[i])*(enrichedDisplacements[j*2+1] - previousEnrichedDisplacements[j*2+1]);
		}
		
		if(parent->isMoved())
			parent->getInverseJacobianMatrix(p_[i], Jinv) ;
		lstrain[i*3 + 0] = vm.eval(x_xi+delta_x_xi, p_[i])*Jinv[0][0] + vm.eval(x_eta+delta_x_eta, p_[i])*Jinv[0][1] ;
		lstrain[i*3 + 1] = vm.eval(y_xi+delta_y_xi, p_[i])*Jinv[1][0] + vm.eval(y_eta+delta_y_eta, p_[i])*Jinv[1][1] ;
		lstrain[i*3 + 2] = 0.5*(vm.eval(x_xi+delta_x_xi, p_[i])*Jinv[1][0] + vm.eval(x_eta+delta_x_eta, p_[i])*Jinv[1][1]  + vm.eval(y_xi+delta_y_xi, p_[i])*Jinv[0][0] + vm.eval(y_eta+delta_y_eta, p_[i])*Jinv[0][1] );	
	}
	return lstrain ;
}

Vector ElementState::getDeltaStrain( const std::valarray<std::pair<Point, double> > & p) const
{
	Vector pts(3*p.size()) ;
	for(size_t i = 0 ; i < p.size() ; i++)
	{
		Vector str = getStrain(p[i].first) ;
		for(size_t j = 0 ; j < 3 ; j++)
		{
			pts[i*3+j] = str[j] ;
		}
	}
	
	return pts ;
}

Vector ElementState::getDeltaDisplacements() const
{
	return this->displacements-this->previousDisplacements ;
}

Vector ElementState::getDeltaDisplacements(const Point & p) const
{
	Function x ;
	Function y ;
	
	Point p_ = parent->inLocalCoordinates(p) ;
	
	for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
	{
		x += parent->getShapeFunction(j)*(displacements[j*2]-previousDisplacements[j*2]) ;
		y += parent->getShapeFunction(j)*(displacements[j*2+1]-previousDisplacements[j*2+1]) ;
	}
	
	double delta_x = 0;
	double delta_y = 0;
	VirtualMachine vm ;
	for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
	{
		delta_x += vm.deval( parent->getEnrichmentFunction(j),XI,p_)*(enrichedDisplacements[j*2] - previousEnrichedDisplacements[j*2]);
		delta_y += vm.deval( parent->getEnrichmentFunction(j),ETA,p_)*(enrichedDisplacements[j*2+1] - previousEnrichedDisplacements[j*2+1]);
	}
	
	Vector ret(2) ;
	ret[0] = vm.eval(x + delta_x, p) ;
	ret[1] = vm.eval(y + delta_y, p) ;
	return ret;
}

Vector ElementState::getDeltaDisplacements(const std::valarray<Point> & p) const
{
	Function x ;
	Function y ;
	
	std::vector<Point> p_ ;
	for(size_t i = 0 ; i < p.size() ; i++)
	{
		p_.push_back(parent->inLocalCoordinates(p[i])) ;
	}
	
	for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
	{
		x += parent->getShapeFunction(j)*(displacements[j*2]- previousDisplacements[j*2]);
		y += parent->getShapeFunction(j)*(displacements[j*2+1]- previousDisplacements[j*2+1]);
	}
	
	for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
	{
		x += parent->getEnrichmentFunction(j).d(XI)*(enrichedDisplacements[j*2] - previousEnrichedDisplacements[j*2]);
		y += parent->getEnrichmentFunction(j).d(ETA)*(enrichedDisplacements[j*2+1] - previousEnrichedDisplacements[j*2+1]);
	}
	
	VirtualMachine vm ;
	Vector ret(2*p.size()) ;
	for(size_t i = 0 ; i < p.size() ; i++)
	{
		double delta_x = 0;
		double delta_y = 0;
		
		for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
		{
			delta_x += vm.deval( parent->getEnrichmentFunction(j),XI,p_[i])*(enrichedDisplacements[j*2] - previousEnrichedDisplacements[j*2]);
			delta_y += vm.deval( parent->getEnrichmentFunction(j),ETA,p_[i])*(enrichedDisplacements[j*2+1]- previousEnrichedDisplacements[j*2+1]) ;
		}
		ret[i*2+0] = vm.eval(x+delta_x, p[i]) ;
		ret[i*2+1] = vm.eval(y+delta_y, p[i]) ;
	}
	return ret;
}

std::pair<Vector, Vector > ElementState::getDeltaStressAndDeltaStrain(const Mu::PointArray & pts) const
{
	
	if(parent->getBehaviour()->type == VOID_BEHAVIOUR)
	{
		return std::pair<Vector, Vector>(Vector(0., 3*pts.size()), Vector(0., 3*pts.size())) ;
	}
	
	std::vector<Point> p_ ;
	for(size_t i = 0 ; i < pts.size() ; i++)
	{
		p_.push_back(parent->inLocalCoordinates(*pts[i])) ;
	}
	Function x_xi ;
	Function x_eta ;
	Function y_xi ;
	Function y_eta ;
	
	for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
	{
		x_xi += parent->getShapeFunction(j).d(XI)*(displacements[j*2]-previousDisplacements[j*2]) ;
		x_eta += parent->getShapeFunction(j).d(ETA)*(displacements[j*2]-previousDisplacements[j*2]) ;
		y_xi += parent->getShapeFunction(j).d(XI)*(displacements[j*2+1]-previousDisplacements[j*2+1]) ;
		y_eta += parent->getShapeFunction(j).d(ETA)*(displacements[j*2+1]- previousDisplacements[j*2+1]);
	}
	
	VirtualMachine vm ;
	
	Vector lstrain(3*pts.size()) ;
	Vector lstress(3*pts.size()) ;
	
	Matrix Jinv(2,2) ; parent->getInverseJacobianMatrix(p_[0], Jinv) ;
	
	for(size_t i = 0 ; i < pts.size() ; i++)
	{
		double delta_x_xi  = 0;
		double delta_x_eta  = 0;
		double delta_y_xi  = 0;
		double delta_y_eta  = 0;
		
		Vector llstrain(3) ;
		Matrix cg (parent->getBehaviour()->getTensor(p_[i])) ;
		
		for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
		{
			delta_x_xi += vm.deval( parent->getEnrichmentFunction(j),XI,p_[i] )*(enrichedDisplacements[j*2] - previousEnrichedDisplacements[j*2]);
			delta_x_eta += vm.deval( parent->getEnrichmentFunction(j),ETA,p_[i])*(enrichedDisplacements[j*2] - previousEnrichedDisplacements[j*2]);
			delta_y_xi += vm.deval( parent->getEnrichmentFunction(j),XI,p_[i])*(enrichedDisplacements[j*2+1] - previousEnrichedDisplacements[j*2+1]);
			delta_y_eta += vm.deval(parent->getEnrichmentFunction(j),ETA, p_[i])*(enrichedDisplacements[j*2+1]- previousEnrichedDisplacements[j*2+1]) ;
		}
		
		if(parent->isMoved())
			parent->getInverseJacobianMatrix(p_[i], Jinv) ;
		
		llstrain[0] = vm.eval(x_xi+delta_x_xi, p_[i])*Jinv[0][0] + vm.eval(x_eta+delta_x_eta, p_[i])*Jinv[0][1] ;
		llstrain[1] = vm.eval(y_xi+delta_y_xi, p_[i])*Jinv[1][0] + vm.eval(y_eta+delta_y_eta, p_[i])*Jinv[1][1] ;
		llstrain[2] = 0.5*(vm.eval(x_xi+delta_x_xi, p_[i])*Jinv[1][0] + vm.eval(x_eta+delta_x_eta, p_[i])*Jinv[1][1]  + vm.eval(y_xi+delta_y_xi, p_[i])*Jinv[0][0] + vm.eval(y_eta+delta_y_eta, p_[i])*Jinv[0][1] );	
		
		lstrain[i*3 + 0] = llstrain[0] ;
		lstrain[i*3 + 1] = llstrain[1] ;
		lstrain[i*3 + 2] = llstrain[2] ;
		
		llstrain=llstrain*cg ;
		
		lstress[i*3 + 0] = llstrain[0] ;
		lstress[i*3 + 1] = llstrain[1] ;
		lstress[i*3 + 2] = llstrain[2] ;
		
	}
	std::pair<Vector, Vector> ret(lstress, lstrain) ;
	
	return ret ;
}

std::pair<Vector, Vector > ElementState::getDeltaStressAndDeltaStrain( const std::valarray<std::pair<Point, double> > & p) const
{
	 Mu::PointArray pts(p.size()) ;
	for(size_t i = 0 ; i < p.size() ; i++)
	{
		pts[i] = new Point(p[i].first) ;
	}
	std::pair<Vector, Vector > ret = getDeltaStressAndDeltaStrain(pts) ;
	for(size_t i = 0 ; i < p.size() ; i++)
		delete pts[i] ;
	
	return ret ;
}

void ElementState::initialize()
{

		size_t ndofs = parent->getBehaviour()->getNumberOfDegreesOfFreedom() ;
		
		displacements.resize(parent->getBoundingPoints().size()*ndofs) ;
		displacements = 0 ;
		previousDisplacements.resize(displacements.size()) ;
		previousDisplacements = 0 ;
		previousPreviousDisplacements.resize(displacements.size()) ;
		previousPreviousDisplacements = 0 ;
		buffer.resize(displacements.size()) ;
		buffer = 0 ;
		timePos = -0.1 ;
		previousTimePos = -0.2 ;
		previousPreviousTimePos = -0.3 ;

		if(parent->getBehaviour()->getFractureCriterion())
		{
			parent->getBehaviour()->getFractureCriterion()->initialiseCache(*this) ;
		}
	
}

const Vector & ElementState::getBuffer() const
{
	return buffer ;
}

Vector & ElementState::getBuffer() 
{
	return buffer ;
}

void ElementState::step(double dt, const Vector * d)
{
	if(!history.empty())
		history.pop_back() ;
	history.push_back(ElementState(*this)) ;
	
	if(parent->getBehaviour()->type != VOID_BEHAVIOUR)
	{
		previousPreviousTimePos = previousTimePos ;
		previousTimePos = timePos ;
		timePos+=dt ;
		size_t ndofs = parent->getBehaviour()->getNumberOfDegreesOfFreedom() ;
		
		std::vector< size_t > ids = parent->getDofIds() ;
		
		if(ids.empty())
			return ;

		if(buffer.size() != parent->getBoundingPoints().size()*ndofs+parent->getEnrichmentFunctions().size()*ndofs)
			buffer.resize(parent->getBoundingPoints().size()*ndofs+parent->getEnrichmentFunctions().size()*ndofs) ;
		for(size_t i = 0 ; i < parent->getShapeFunctions().size() ; i++)
		{
			for(size_t j = 0 ; j < ndofs ; j++)
			{
				buffer[i*ndofs+j] = (*d)[ids[i]*ndofs+j] ;
			}
		}
		
		int nbp = parent->getBoundingPoints().size() ;
		for(size_t i = 0 ; i < parent->getEnrichmentFunctions().size() ; i++)
		{
			for(size_t j = 0 ; j < ndofs ; j++)
				buffer[i*ndofs+nbp*ndofs+j] = (*d)[ids[i+nbp]*ndofs+j] ;
		}
	}
}


Vector ElementState::getSpeed(const Point & p) const
{
	return (getDisplacements(p) - getPreviousDisplacements(p)) / (timePos-previousTimePos) ;
}

Vector ElementState::getSpeed(const std::valarray<Point> & p) const
{
	return (getDisplacements(p) - getPreviousDisplacements(p)) / (timePos-previousTimePos) ;
}

Vector ElementState::getSpeed() const 
{
	return (getDisplacements() - getPreviousDisplacements()) / (timePos-previousTimePos) ;
}
	
Vector ElementState::getAcceleration(const Point &) const 
{
	return ((getDisplacements() - getPreviousDisplacements()) / (timePos-previousTimePos) - (getPreviousDisplacements() - getPreviousPreviousDisplacements()) / (previousTimePos-previousPreviousTimePos))/(timePos-previousPreviousTimePos) ;
}

Vector ElementState::getAcceleration(const std::valarray<Point> & p) const 
{
	return ((getDisplacements(p) - getPreviousDisplacements(p)) / (timePos-previousTimePos) - (getPreviousDisplacements(p) - getPreviousPreviousDisplacements(p)) / (previousTimePos-previousPreviousTimePos))/(timePos-previousPreviousTimePos) ;
}

Vector ElementState::getAcceleration() const
{
	return ((getDisplacements() - getPreviousDisplacements()) / (timePos-previousTimePos) - (getPreviousDisplacements() - getPreviousPreviousDisplacements()) / (previousTimePos-previousPreviousTimePos))/(timePos-previousPreviousTimePos) ;
}

double ElementState::getTime() const
{
	return timePos ;
}

double ElementState::getDeltaTime() const
{
	return timePos - previousTimePos ;
}


double ElementState::getPrincipalAngle(const Point & p, bool local) const
{
	Vector stresses = getStress(p, local) ;
// 	double a0 = 0.5*atan2(stresses[0]-stresses[1],-stresses[2]) + M_PI ;
// 	double a1 = 0.5*atan2(stresses[0]-stresses[1],-stresses[2]) + M_PI/4. + M_PI ;
// 	double a2 = 0.5*atan2(stresses[0]-stresses[1],-stresses[2]) + M_PI/4. + M_PI ;
// 	if(std:: abs(a0) < std:: abs(a1) && std:: abs(a0) < std:: abs(a2))
// 		return a0 ;
// 	else if (std:: abs(a1) < std:: abs(a0) && std:: abs(a1) < std:: abs(a2))
// 		return a1 ;
// 	else
// 		return a2 ;
	
	return 0.5*atan2(stresses[0]-stresses[1],-stresses[2]) ;
}

Vector ElementState::getPrincipalAngle(const Mu::PointArray & v) const
{
	Vector stresses = getStress(v) ;
	Vector principal(v.size()) ;
	for(size_t i = 0 ; i < v.size() ; i++)
	{
// 		double a0 = 0.5*atan2(stresses[i*3+0]-stresses[i*3+1],-stresses[i*3+2]) + M_PI ;
// 		double a1 = 0.5*atan2(stresses[i*3+0]-stresses[i*3+1],-stresses[i*3+2]) + M_PI/4. + M_PI;
// 		double a2 = 0.5*atan2(stresses[i*3+0]-stresses[i*3+1],-stresses[i*3+2]) + M_PI/4. + M_PI ;
// 		
// 		if(std:: abs(a0) < std:: abs(a1) && std:: abs(a0) < std:: abs(a2))
// 			principal[i] = a0 ;
// 		else if (std:: abs(a1) < std:: abs(a0) && std:: abs(a1) < std:: abs(a2))
// 			principal[i] = a1 ;
// 		else
// 			principal[i] = a2 ;
		principal[i] =  0.5*atan2(stresses[i*3+0]-stresses[i*3+1],-stresses[i*3+2]) ;
	}
	return principal ;
}

Vector ElementState::getPrincipalStresses(const Point & p, bool local ) const
{
	if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	{
		Vector stresses = getStress(p,local) ;
	
		if(parent->getBehaviour()->hasInducedForces())
			stresses -= parent->getBehaviour()->getImposedStress(p) ;
		Vector lprincipal(2) ;
		lprincipal[0] = (stresses[0]+stresses[1])/2. + 
			sqrt(
				(stresses[0]-stresses[1])*(stresses[0]-stresses[1])/4. + 
				(stresses[2]*stresses[2])
				) ;
		lprincipal[1] = (stresses[0]+stresses[1])/2. - 
			sqrt(
				(stresses[0]-stresses[1])*(stresses[0]-stresses[1])/4. + 
				(stresses[2]*stresses[2])
				) ;
		

		return lprincipal ;
	}
	else
	{
		Matrix stresses = getStressMatrix(p,local) ;
		Matrix I(3,3) ; I[0][0] = 1 ; I[1][1] = 1 ; I[2][2] = 1 ;
		double m = (stresses[0][0] + stresses[1][1] + stresses[2][2])/3. ;
		Matrix Am = stresses-I*m ;
		double q = det(Am)/2. ;
		double r = std::inner_product(&Am.array()[0], &Am.array()[9],&Am.array()[0],  double(0.))/6. ;
		double phi = atan2(sqrt(r*r*r-q*q),q)/3. ;
		if(r*r*r-q*q < 1e-12)
			phi = atan(0)/3. ;
		
		if(phi < 0)
			phi+= M_PI ;
		
		Vector lprincipal(3) ;
		lprincipal[0] = m+2.*sqrt(r)*cos(phi) ;
		lprincipal[1] = m-sqrt(r)*(cos(phi)+sqrt(3)*sin(phi)) ;
		lprincipal[2] = m-sqrt(r)*(cos(phi)-sqrt(3)*sin(phi)) ;
		return lprincipal ;
	}
}

Vector ElementState::getPrincipalStrains(const Point & p, bool local ) const
{
	if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	{
		Vector strains = getStrain(p,local) ;
	
// 		if(parent->getBehaviour()->hasInducedForces())
// 			stresses -= parent->getBehaviour()->getImposedStrains(p) ;
		Vector lprincipal(2) ;
		lprincipal[0] = (strains[0]+strains[1])/2. + 
			sqrt(
				(strains[0]-strains[1])*(strains[0]-strains[1])/4. + 
				(strains[2]*strains[2])
				) ;
		lprincipal[1] = (strains[0]+strains[1])/2. - 
			sqrt(
				(strains[0]-strains[1])*(strains[0]-strains[1])/4. + 
				(strains[2]*strains[2])
				) ;
		

		return lprincipal ;
	}
	else
	{
		Matrix strains = getStrainMatrix(p,local) ;
		Matrix I(3,3) ; I[0][0] = 1 ; I[1][1] = 1 ; I[2][2] = 1 ;
		double m = (strains[0][0] + strains[1][1] + strains[2][2])/3. ;
		Matrix Am = strains-I*m ;
		double q = det(Am)/2. ;
		double r = std::inner_product(&Am.array()[0], &Am.array()[9],&Am.array()[0],  double(0.))/6. ;
		double phi = atan2(sqrt(r*r*r-q*q),q)/3. ;
		if(r*r*r-q*q < 1e-12)
			phi = atan(0)/3. ;
		
		if(phi < 0)
			phi+= M_PI ;
		
		Vector lprincipal(3) ;
		lprincipal[0] = m+2.*sqrt(r)*cos(phi) ;
		lprincipal[1] = m-sqrt(r)*(cos(phi)+sqrt(3)*sin(phi)) ;
		lprincipal[2] = m-sqrt(r)*(cos(phi)-sqrt(3)*sin(phi)) ;
		return lprincipal ;
	}
}

double ElementState::elasticEnergy() const
{
	GaussPointArray gp = parent->getGaussPoints() ;

	Vector stress = getStress(gp.gaussPoints) ;
	Vector strain = getStrain(gp.gaussPoints) ;
	size_t stresscomponents = stress.size()/gp.gaussPoints.size() ;
	double e = 0 ;
	for(size_t i = 0 ;i < gp.gaussPoints.size() ;i++)
	{
		double le = 0 ;
		for(size_t j = 0 ; j < stresscomponents ; j++)
		{
			le += stress[i*stresscomponents+j]*strain[i*stresscomponents+j] ;
		}
		e += .5*le*gp.gaussPoints[i].second ;
	}
	return e ;

}

Vector ElementState::getPrincipalStrains(const Mu::PointArray & v) const
{
	Vector strains = getStrain(v) ;
	Vector principal(2*v.size()) ;
	
	for(size_t i = 0 ; i < v.size() ; i++)
	{
		Vector lprincipal(2) ;
		lprincipal[0] = (strains[i*3]+strains[i*3+1])/2. +
			sqrt(
			      (strains[i*3]-strains[i*3+1])*(strains[i*3]-strains[i*3+1])/4. +
			      strains[i*3+2]*strains[i*3+2]*4.
			    ) ;
		lprincipal[1] = (strains[i*3]+strains[i*3+1])/2. -
			sqrt(
			      (strains[i*3]-strains[i*3+1])*(strains[i*3]-strains[i*3+1])/4. +
			      strains[i*3+2]*strains[i*3+2]*4.
			    ) ;
		principal[i*2] = lprincipal[0] ;
		principal[i*2+1] = lprincipal[1] ;
	}
	
	return principal ;
}

Vector ElementState::getPrincipalStresses(const Mu::PointArray& v) const
{
	Vector stresses = getStress(v) ;
	Vector principal(2*v.size()) ;
	
	for(size_t i = 0 ; i < v.size() ; i++)
	{
		Vector lprincipal(2) ;
		lprincipal[0] = (stresses[i*3]+stresses[i*3+1])/2. +
			sqrt(
			      (stresses[i*3]-stresses[i*3+1])*(stresses[i*3]-stresses[i*3+1])/4. +
			      stresses[i*3+2]*stresses[i*3+2]*4.
			    ) ;
		lprincipal[1] = (stresses[i*3]+stresses[i*3+1])/2. -
			sqrt(
			      (stresses[i*3]-stresses[i*3+1])*(stresses[i*3]-stresses[i*3+1])/4. +
			      stresses[i*3+2]*stresses[i*3+2]*4.
			    ) ;
		principal[i*2] = lprincipal[0] ;
		principal[i*2+1] = lprincipal[1] ;
	}
	
	return principal ;
}
