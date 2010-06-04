/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010  Cyrille Dunant

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

*/

#include "mad.h"

using namespace Mu ;

MultipleAggregatingDiscontinuities::MultipleAggregatingDiscontinuities(FeatureTree * ft, DelaunayTriangle * self ,  FractureCriterion * crit) : LinearForm(Matrix(), false, false, 2), featureTree(ft) , self2d(self), self3d(NULL), leastSquares(NULL), equivalentCrack(NULL)
{
	v.push_back(XI);
	v.push_back(ETA);
	bc = new ElementDefinedBoundaryCondition(self) ;
	featureTree->addBoundaryCondition(bc);
	homogenize() ;
}

MultipleAggregatingDiscontinuities::MultipleAggregatingDiscontinuities(FeatureTree * ft, DelaunayTetrahedron * self ,  FractureCriterion * crit) : LinearForm(Matrix(), false, false, 2), featureTree(ft) , self2d(NULL), self3d(self), leastSquares(NULL), equivalentCrack(NULL)
{
	
	v.push_back(XI);
	v.push_back(ETA);
	v.push_back(ZETA);
	bc = new ElementDefinedBoundaryCondition(self) ;
	featureTree->addBoundaryCondition(bc);
	homogenize() ;
}


/** \brief Apply the law.
	*
	* The matrix is computed as: \f$ \nabla^T h_i K \nabla h_j \f$
	* @param p_i first basis polynomial.
	* @param p_j second basis polynomial.
	* @param gp Gauss Points used for the quadrature
	* @param Jinv Inverse Jacobian Matrices corresponding to the gauss points
	* @param ret Matrix to store the result
	* @param vm virtualMachine to use to compute the result
	*/
void MultipleAggregatingDiscontinuities::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v, ret) ;
}

void MultipleAggregatingDiscontinuities::step(double timestep, ElementState & currentState)
{

	if(type == VOID_BEHAVIOUR)
		return ;
	
	if(v.size() == 2)
	{
		if(equivalentCrack && currentState.getParent()->get2DMesh())
			equivalentCrack->enrich(currentState.getParent()->get2DMesh()->getLastNodeId(), currentState.getParent()->get2DMesh());
	}
	
	if(!crit->met(currentState))
		return ;
	
	featureTree->generateElements(64);
	while(!featureTree->step(timestep)) 
	{
		
	}
	
	homogenize() ;
	
	
}


void MultipleAggregatingDiscontinuities::homogenize()
{
	std::vector<Material> mat ;
	std::vector<Tag> bulkshear ;
	bulkshear.push_back(TAG_BULK_MODULUS) ;
	bulkshear.push_back(TAG_SHEAR_MODULUS) ;
	double totalArea ;
	if(self2d)
	{
		for(size_t i = 0 ; i < featureTree->getFeatures().size() ; i++)
		{
			if(featureTree->getFeature(i)->getBehaviour())
			{
				time_d = time_d || featureTree->getFeature(i)->getBehaviour()->timeDependent() ;
				space_d = space_d || featureTree->getFeature(i)->getBehaviour()->spaceDependent() ;
				
				if(featureTree->getFeature(i)->getBehaviour()->type != VOID_BEHAVIOUR)
				{
					Matrix param_gp(param.numRows(),param.numCols()) ;				

					param_gp = featureTree->getFeature(i)->getBehaviour()->getTensor(featureTree->getFeature(i)->getCenter()) ;
					Material thismat(param_gp) ;
					thismat.push_back(Properties(TAG_VOLUME,featureTree->getFeature(i)->area())) ;
					if(mat.size() == 0)
						mat.push_back(thismat) ;
					else 
					{
						bool combined = false ;
						for(size_t j = 0 ; j < mat.size() ; j++)
							combined = (combined || mat[j].combine(thismat,bulkshear,TAG_VOLUME)) ;
						if(!combined)
							mat.push_back(thismat) ;
					}
					totalArea += featureTree->getFeature(i)->area() ;
				}
				else
				{
					Material thismat(Properties(TAG_VOLUME, featureTree->getFeature(i)->area())) ;
					thismat.push_back(Properties(TAG_BULK_MODULUS, 1e-9)) ;
					thismat.push_back(Properties(TAG_SHEAR_MODULUS, 1e-9)) ;
					if(mat.size() == 0)
						mat.push_back(thismat) ;
					else 
					{
						bool combined = false ;
						for(size_t j = 0 ; j < mat.size() ; j++)
							combined = (combined || mat[j].combine(thismat,bulkshear,TAG_VOLUME)) ;
						if(!combined)
							mat.push_back(thismat) ;
					}
					totalArea += featureTree->getFeature(i)->area() ;
				}
			}
		}
	} else {
		for(size_t i = 0 ; i < featureTree->getFeatures().size() ; i++)
		{
			if(featureTree->getFeature(i)->getBehaviour())
			{
				time_d = time_d || featureTree->getFeature(i)->getBehaviour()->timeDependent() ;
				space_d = space_d || featureTree->getFeature(i)->getBehaviour()->spaceDependent() ;
				
				if(featureTree->getFeature(i)->getBehaviour()->type != VOID_BEHAVIOUR)
				{
					Matrix param_gp(param.numRows(),param.numCols()) ;				
					param_gp = featureTree->getFeature(i)->getBehaviour()->getTensor(featureTree->getFeature(i)->getCenter())  ;
					Material thismat(param_gp) ;
					thismat.push_back(Properties(TAG_VOLUME,featureTree->getFeature(i)->volume())) ;
					if(mat.size() == 0)
						mat.push_back(thismat) ;
					else 
					{
						bool combined = false ;
						for(size_t j = 0 ; j < mat.size() ; j++)
							combined = (combined || mat[j].combine(thismat,bulkshear,TAG_VOLUME)) ;
						if(!combined)
							mat.push_back(thismat) ;
					}
					totalArea += featureTree->getFeature(i)->volume() ;
				}
				else
				{
					Material thismat(Properties(TAG_VOLUME, featureTree->getFeature(i)->volume())) ;
					thismat.push_back(Properties(TAG_BULK_MODULUS, 1e-9)) ;
					thismat.push_back(Properties(TAG_SHEAR_MODULUS, 1e-9)) ;
					if(mat.size() == 0)
						mat.push_back(thismat) ;
					else 
					{
						bool combined = false ;
						for(size_t j = 0 ; j < mat.size() ; j++)
							combined = (combined || mat[j].combine(thismat,bulkshear,TAG_VOLUME)) ;
						if(!combined)
							mat.push_back(thismat) ;
					}
					totalArea += featureTree->getFeature(i)->volume() ;
				}
			}
		}

	}

	GeneralConverter fraction(TAG_VOLUME_FRACTION) ;
	AdditionConverter total(TAG_VOLUME_TOTAL) ;

	for(size_t i = 0 ; i < mat.size() ; i++)
	{
		fraction.reset() ;
		total.reset() ;
		mat[i].merge(total.homogenize(mat)) ;
		mat[i].merge(fraction.homogenize(mat[i])) ;
	}

//	std::cout << mat.size() << std::endl ;
	if(mat.size() == 1)
	{
		if(self2d)
			param = cauchyGreen(std::make_pair(mat[0].val(TAG_BULK_MODULUS,-1),mat[0].val(TAG_SHEAR_MODULUS,-1)),false,SPACE_TWO_DIMENSIONAL) ;
		else
			param = cauchyGreen(std::make_pair(mat[0].val(TAG_BULK_MODULUS,-1),mat[0].val(TAG_SHEAR_MODULUS,-1)),false,SPACE_THREE_DIMENSIONAL) ;
	} else {
		scheme.reset() ;
		Material mhom(scheme.homogenize(mat)) ;
		int nk = mhom.getIndex(TAG_BULK_MODULUS,-1) ;
		int nmu = mhom.getIndex(TAG_SHEAR_MODULUS,-1) ;
		if(scheme.isOK() && (nk+1)*(nmu+1) > 0)
		{
			if(self2d)
				param = cauchyGreen(std::make_pair(mhom[nk].val(),mhom[nmu].val()),false,SPACE_TWO_DIMENSIONAL) ;
			else
				param = cauchyGreen(std::make_pair(mhom[nk].val(),mhom[nmu].val()),false,SPACE_THREE_DIMENSIONAL) ;
		} else {
			scheme.print() ;
		}
	}
//	std::cout << std::endl ;
//	std::cout << std::endl ;
//	param.print() ;
//	std::cout << std::endl ;
//	std::cout << std::endl ;*/
	
}

std::vector<BoundaryCondition * > MultipleAggregatingDiscontinuities::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
	std::vector<BoundaryCondition * > ret ;
	if(self2d)
		ret.push_back(new ElementDefinedBoundaryCondition(self2d));
	else
		ret.push_back(new ElementDefinedBoundaryCondition(self3d));
	
	return ret ;
}
