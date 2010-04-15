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
#include "homogenization/elastic_homogenization.h"

using namespace Mu ;

HomogeneisedBehaviour::HomogeneisedBehaviour(Mesh<DelaunayTriangle, DelaunayTreeItem> * mesh2d, DelaunayTriangle * self) : LinearForm(Matrix(), false, false, 2) , mesh2d(mesh2d), self2d(self), mesh3d(NULL), self3d(NULL)
{
	Circle c(self->getRadius()*2, self->getCircumCenter()) ;
	source = mesh2d->getConflictingElements(self->getPrimitive()) ;
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

	homogenize() ;	

	if(param.isNull())
		type = VOID_BEHAVIOUR ;
//	else
//		param /= totalArea ;
	v.push_back(XI);
	v.push_back(ETA);
} ;

HomogeneisedBehaviour::HomogeneisedBehaviour(Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * mesh3d, DelaunayTetrahedron * self) : LinearForm(Matrix(), false, false, 3), mesh2d(NULL), self2d(NULL), mesh3d(mesh3d), self3d(self)
{
	source3d = mesh3d->getConflictingElements(self->getPrimitive()) ;
	//simple averaging
	double totalVolume = 0 ;
	for(size_t i = 0 ; i < source3d.size() ; i++)
	{
		if(source3d[i]->getBehaviour() && source3d[i]->getBehaviour()->type != VOID_BEHAVIOUR)
		{
			param = Matrix(source3d[i]->getBehaviour()->param.numRows(), source3d[i]->getBehaviour()->param.numCols()) ;
			break ;
		}
	}

	homogenize() ;
	
	if(param.isNull())
		type = VOID_BEHAVIOUR ;
//	else
//		param /= totalVolume ;
	
//	param[3][3] /= .9 ;
//	param[4][4] /= .9 ;
//	param[5][5] /= .9 ;
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
	homogenize() ;
}

void HomogeneisedBehaviour::stepBack()
{
	if(type == VOID_BEHAVIOUR)
		return ;
	homogenize() ;
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

void HomogeneisedBehaviour::homogenize()
{
	std::vector<Material> mat ;
	double totalArea ;
	if(self2d)
	{
		for(size_t i = 0 ; i < source.size() ; i++)
		{
			if(source[i]->getBehaviour())
			{
				time_d = time_d || source[i]->getBehaviour()->timeDependent() ;
				space_d = space_d || source[i]->getBehaviour()->spaceDependent() ;
				
				GaussPointArray gp = source[i]->getGaussPoints() ;
				if(source[i]->getBehaviour()->type != VOID_BEHAVIOUR)
				{
					Matrix param_gp(param.numRows(),param.numCols()) ;				
					double area_gp = 0 ;
					for(size_t j = 0 ; j < gp.gaussPoints.size() ; j++)
					{
						param_gp += source[i]->getBehaviour()->getTensor(gp.gaussPoints[j].first) * gp.gaussPoints[j].second  ;
						area_gp += gp.gaussPoints[j].second ;
					}
					mat.push_back(Material(Properties(FRACTION,source[i]->area()))) ;
					mat[i].push_back(Properties(BULK_SHEAR,param_gp / area_gp)) ;
					totalArea += source[i]->area() ;
				}
				else
				{
					mat.push_back(Material(Properties(FRACTION,source[i]->area()))) ;
					mat[i].push_back(Properties(HOOKE,std::make_pair(1e-9,0.2))) ;
					mat[i].push_back(mat[i][1].convert(BULK_SHEAR).second) ;
					totalArea += source[i]->area() ;
				}
			}
		}
	} else {
		for(size_t i = 0 ; i < source3d.size() ; i++)
		{
			if(source3d[i]->getBehaviour())
			{
				time_d = time_d || source3d[i]->getBehaviour()->timeDependent() ;
				space_d = space_d || source3d[i]->getBehaviour()->spaceDependent() ;
				
				GaussPointArray gp = source3d[i]->getGaussPoints() ;
				if(source[i]->getBehaviour()->type != VOID_BEHAVIOUR)
				{
					Matrix param_gp(param.numRows(),param.numCols()) ;				
					double area_gp = 0 ;
					for(size_t j = 0 ; j < gp.gaussPoints.size() ; j++)
					{
						param_gp += source3d[i]->getBehaviour()->getTensor(gp.gaussPoints[j].first) * gp.gaussPoints[j].second  ;
						area_gp += gp.gaussPoints[j].second ;
					}
					mat.push_back(Material(Properties(FRACTION,source3d[i]->volume()))) ;
					mat[i].push_back(Properties(BULK_SHEAR,param_gp / area_gp)) ;
					totalArea += source3d[i]->volume() ;
				}
				else
				{
					mat.push_back(Material(Properties(FRACTION,source3d[i]->volume()))) ;
					mat[i].push_back(Properties(HOOKE,std::make_pair(1e-9,0.2))) ;
					mat[i].push_back(mat[i][1].convert(BULK_SHEAR).second) ;
					totalArea += source3d[i]->volume() ;
				}
			}
		}

	}

	double f = 0 ;
	for(size_t i = 0 ; i < mat.size() ; i++)
		f += mat[i][0].getValue(0) ;
	for(size_t i = 0 ; i < mat.size() ; i++)
		mat[i][0].setValue(0,mat[i][0].getValue(0)/f) ;
	std::pair<bool,Material> hom = scheme.apply(mat) ;
	if(hom.first)
	{
		std::pair<bool,Matrix> param_hom = hom.second[0].getCauchyGreen(SPACE_TWO_DIMENSIONAL) ;
		if(param_hom.first)
			param = param_hom.second ;
	}

	
}

