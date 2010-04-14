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

#include "homogeneised_behaviour.h"

using namespace Mu ;

HomogeneisedBehaviour::HomogeneisedBehaviour(Mesh<DelaunayTriangle, DelaunayTreeItem> * mesh2d, DelaunayTriangle * self) : LinearForm(Matrix(), false, false, 2) , mesh2d(mesh2d), self2d(self), mesh3d(NULL), self3d(NULL) 
{
	Circle c(self->getRadius()*2, self->getCircumCenter()) ;
	std::vector<DelaunayTriangle *> source = mesh2d->getConflictingElements(self->getPrimitive()) ;
	//simple averaging
	double totalArea = 0 ;
	for(size_t i = 0 ; i < source.size() ; i++)
	{
		if(source[i]->getBehaviour() && source[i]->getBehaviour()->type != VOID_BEHAVIOUR)
		{
			param = Matrix(source[i]->getBehaviour()->param.numRows(), source[i]->getBehaviour()->param.numCols()) ;
			break ;
		}
	}
	
	for(size_t i = 0 ; i < source.size() ; i++)
	{
		if(source[i]->getBehaviour())
		{
			time_d = time_d || source[i]->getBehaviour()->timeDependent() ;
			space_d = space_d || source[i]->getBehaviour()->spaceDependent() ;
			
// 			double area = source[i]->area() ;
// 			totalArea += area ;
			
			GaussPointArray gp = source[i]->getGaussPoints() ;
			if(source[i]->getBehaviour()->type != VOID_BEHAVIOUR)
			{
				for(size_t j = 0 ; j < gp.gaussPoints.size() ; j++)
				{
					param += source[i]->getBehaviour()->getTensor(gp.gaussPoints[j].first) * gp.gaussPoints[j].second  ;
					totalArea += gp.gaussPoints[j].second ;
				}
			}
			else
			{
				for(size_t j = 0 ; j < gp.gaussPoints.size() ; j++)
				{
					totalArea += gp.gaussPoints[j].second ;
				}
			}
		}
	}
	if(param.isNull())
		type = VOID_BEHAVIOUR ;
	else
		param /= totalArea ;
	v.push_back(XI);
	v.push_back(ETA);
} ;

HomogeneisedBehaviour::HomogeneisedBehaviour(Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * mesh3d, DelaunayTetrahedron * self) : LinearForm(Matrix(), false, false, 3), mesh2d(NULL), self2d(NULL), mesh3d(mesh3d), self3d(self) 
{
	std::vector<DelaunayTetrahedron *> source = mesh3d->getConflictingElements(self->getPrimitive()) ;
	//simple averaging
	double totalVolume = 0 ;
	for(size_t i = 0 ; i < source.size() ; i++)
	{
		if(source[i]->getBehaviour() && source[i]->getBehaviour()->type != VOID_BEHAVIOUR)
		{
			param = Matrix(source[i]->getBehaviour()->param.numRows(), source[i]->getBehaviour()->param.numCols()) ;
			break ;
		}
	}
	for(size_t i = 0 ; i < source.size() ; i++)
	{
		if(source[i]->getBehaviour())
		{
			time_d = time_d || source[i]->getBehaviour()->timeDependent() ;
			space_d = space_d || source[i]->getBehaviour()->spaceDependent() ;
			
			GaussPointArray gp = source[i]->getGaussPoints() ;
			if(source[i]->getBehaviour()->type != VOID_BEHAVIOUR)
			{
				for(size_t j = 0 ; j < gp.gaussPoints.size() ; j++)
				{
					param += source[i]->getBehaviour()->getTensor(gp.gaussPoints[j].first) * gp.gaussPoints[j].second ;
					totalVolume += gp.gaussPoints[j].second ;
				}
			}
			else
			{
				for(size_t j = 0 ; j < gp.gaussPoints.size() ; j++)
				{
					totalVolume += gp.gaussPoints[j].second ;
				}
			}
		}
	}
	
	if(param.isNull())
		type = VOID_BEHAVIOUR ;
	else
		param /= totalVolume ;
	
	param[3][3] /= .9 ;
	param[4][4] /= .9 ;
	param[5][5] /= .9 ;
	v.push_back(XI);
	v.push_back(ETA);
	v.push_back(ZETA);
} ;

HomogeneisedBehaviour::~HomogeneisedBehaviour() { } ;

Matrix HomogeneisedBehaviour::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{
	return VirtualMachine().ieval(Gradient(p_i) * param * Gradient(p_j, true), e,v) ;
}

void HomogeneisedBehaviour::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
// std::cout << "a--" << std::endl ;
// Jinv[0].print() ;
// std::cout << "--b" << std::endl ;
	vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v, ret) ;
// 	ret.print() ;
}

void HomogeneisedBehaviour::step(double timestep, ElementState & currentState)
{
	if(type == VOID_BEHAVIOUR)
		return ;
	if(self2d)
	{
		std::vector<DelaunayTriangle *> source = mesh2d->getConflictingElements(self2d->getPrimitive()) ;
		//simple averaging
		double totalArea = 0 ;
		for(size_t i = 0 ; i < source.size() ; i++)
		{
			if(source[i]->getBehaviour())
			{
				time_d = time_d || source[i]->getBehaviour()->timeDependent() ;
				space_d = space_d || source[i]->getBehaviour()->spaceDependent() ;
				
				
				GaussPointArray gp = source[i]->getGaussPoints() ;
				if(source[i]->getBehaviour()->type != VOID_BEHAVIOUR)
				{
					for(size_t j = 0 ; j < gp.gaussPoints.size() ; j++)
					{
						param += source[i]->getBehaviour()->getTensor(gp.gaussPoints[j].first) * gp.gaussPoints[j].second  ;
						totalArea += gp.gaussPoints[j].second ;
					}
				}
				else
				{
					for(size_t j = 0 ; j < gp.gaussPoints.size() ; j++)
					{
						totalArea += gp.gaussPoints[j].second ;
					}
				}
			}
		}
		
		param /= totalArea ;
	}
	else
	{
		std::vector<DelaunayTetrahedron *> source = mesh3d->getConflictingElements(self3d->getPrimitive()) ;
		//simple averaging
		double totalVolume = 0 ;
		for(size_t i = 0 ; i < source.size() ; i++)
		{
			if(source[i]->getBehaviour())
			{
				time_d = time_d || source[i]->getBehaviour()->timeDependent() ;
				space_d = space_d || source[i]->getBehaviour()->spaceDependent() ;
				
				GaussPointArray gp = source[i]->getGaussPoints() ;
				if(source[i]->getBehaviour()->type != VOID_BEHAVIOUR)
				{
					for(size_t j = 0 ; j < gp.gaussPoints.size() ; j++)
					{
						param += source[i]->getBehaviour()->getTensor(gp.gaussPoints[j].first) * gp.gaussPoints[j].second ;
						totalVolume += gp.gaussPoints[j].second ;
					}
				}
				else
				{
					for(size_t j = 0 ; j < gp.gaussPoints.size() ; j++)
					{
						totalVolume += gp.gaussPoints[j].second ;
					}
				}
			}
		}

		param /= totalVolume ;
	}
}

void HomogeneisedBehaviour::stepBack()
{
	if(type == VOID_BEHAVIOUR)
		return ;
	if(self2d)
	{
		std::vector<DelaunayTriangle *> source = mesh2d->getConflictingElements(self2d->getPrimitive()) ;
		//simple averaging
		double totalArea = 0 ;
		param = Matrix(source[0]->getBehaviour()->param.numRows(), source[0]->getBehaviour()->param.numCols()) ;
		for(size_t i = 0 ; i < source.size() ; i++)
		{
			time_d = time_d || source[i]->getBehaviour()->timeDependent() ;
			space_d = space_d || source[i]->getBehaviour()->spaceDependent() ;
			
			double area = source[i]->area() ;
			totalArea += area ;
			
			GaussPointArray gp = source[i]->getGaussPoints() ;
			for(size_t j = 0 ; j < gp.gaussPoints.size() ; j++)
			{
				param += source[i]->getBehaviour()->getTensor(gp.gaussPoints[j].first) * gp.gaussPoints[j].second / area ;
			}
		}
		
		param /= totalArea ;
	}
	else
	{
		std::vector<DelaunayTetrahedron *> source = mesh3d->getConflictingElements(self3d->getPrimitive()) ;
		//simple averaging
		double totalVolume = 0 ;
		param = Matrix(source[0]->getBehaviour()->param.numRows(), source[0]->getBehaviour()->param.numCols()) ;
		for(size_t i = 0 ; i < source.size() ; i++)
		{
			time_d = time_d || source[i]->getBehaviour()->timeDependent() ;
			space_d = space_d || source[i]->getBehaviour()->spaceDependent() ;
			
			double volume = source[i]->area() ;
			totalVolume += volume ;
			
			GaussPointArray gp = source[i]->getGaussPoints() ;
			for(size_t j = 0 ; j < gp.gaussPoints.size() ; j++)
			{
				param += source[i]->getBehaviour()->getTensor(gp.gaussPoints[j].first) * gp.gaussPoints[j].second / volume ;
			}
		}
		
		param /= totalVolume ;
	}
}


bool HomogeneisedBehaviour::fractured() const
{
	return false ;
}

Form * HomogeneisedBehaviour::getCopy() const 
{
	return new HomogeneisedBehaviour(*this) ;
}

void HomogeneisedBehaviour::getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector & f) const 
{
}

