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
#include "../geometry/geometry_base.h" 
#include "homogenization/elastic_homogenization.h"
#include "homogenization/properties_base.h"
#include "homogenization/scheme_base.h"
#include "homogenization/converter.h"
#include "../features/features.h"
#include "../features/inclusion.h"
#include "../features/sample.h"
#include "void_form.h"


using namespace Mu ;

HomogeneisedBehaviour::HomogeneisedBehaviour(Mesh<DelaunayTriangle, DelaunayTreeItem> * mesh2d, DelaunayTriangle * self) : LinearForm(Matrix(), false, false, 2) , mesh2d(mesh2d), self2d(self), mesh3d(NULL), self3d(NULL)
{
	source = mesh2d->getConflictingElements(self2d->getPrimitive()) ;
	//simple averaging
	double totalArea = 0 ;
	for(size_t i = 0 ; i < source.size() ; i++)
	{
		if(source[i]->getBehaviour() && source[i]->getBehaviour()->type != VOID_BEHAVIOUR)
		{
			Form * copy = source[i]->getBehaviour()->getCopy() ;
			source[i]->setBehaviour(copy) ;
			param = Matrix(source[i]->getBehaviour()->param.numRows(), source[i]->getBehaviour()->param.numCols()) ;
			break ;
		}
	}

	homogenize() ;	

	if(param.isNull())
		type = VOID_BEHAVIOUR ;

	v.push_back(XI);
	v.push_back(ETA);
} ;

HomogeneisedBehaviour::HomogeneisedBehaviour(std::vector<Feature *> feats, DelaunayTriangle * self) : LinearForm(Matrix(), false, false, 2), self2d(self), mesh3d(NULL), self3d(NULL)
{
	std::vector<Point> corner = self->getSamplingBoundingPoints(3) ;

	TriangularInclusion * tri = new TriangularInclusion(corner[0],corner[1],corner[2]) ;
	tri->setBehaviour(self->getBehaviour()) ;

	Point c = tri->getCircumCenter() ;
	double r = tri->getRadius()*1.1 ;

	Sample * box = new Sample(NULL,r*2.,r*2.,c.x,c.y) ;
	box->setBehaviour(new VoidForm()) ;

	subTree = new FeatureTree(box) ;
	subTree->addFeature(box, tri) ;

	for(size_t i = 0 ; i < feats.size() ; i++)
		subTree->addFeature(tri, feats[i]) ;

	subTree->sample(50) ;
	subTree->setOrder(LINEAR) ;
	subTree->generateElements() ;
	subTree->getTriangles(-1) ;
	mesh2d = subTree->get2DMesh(-1) ;
	source = mesh2d->getConflictingElements(self2d->getPrimitive()) ;
	double totalArea = 0 ;
	for(size_t i = 0 ; i < source.size() ; i++)
	{
		if(source[i]->getBehaviour() && source[i]->getBehaviour()->type != VOID_BEHAVIOUR)
		{
			Form * copy = source[i]->getBehaviour()->getCopy() ;
			source[i]->setBehaviour(copy) ;
			param = Matrix(source[i]->getBehaviour()->param.numRows(), source[i]->getBehaviour()->param.numCols()) ;
			break ;
		}
	}

	homogenize() ;	

	if(param.isNull())
		type = VOID_BEHAVIOUR ;

	v.push_back(XI);
	v.push_back(ETA);

}


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


void HomogeneisedBehaviour::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	equivalent->apply(p_i,p_j,gp,Jinv,ret,vm) ;
// std::cout << "a--" << std::endl ;
// Jinv[0].print() ;
// std::cout << "--b" << std::endl ;
//	vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v, ret) ;
// 	ret.print() ;
}

void HomogeneisedBehaviour::step(double timestep, ElementState & currentState)
{
	if(type == VOID_BEHAVIOUR)
		return ;

	std::cout << "in step()" << std::endl ;
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
	return equivalent->getCopy() ;
}

void HomogeneisedBehaviour::homogenize()
{
	Material mat ;
	std::vector<Tag> bulkshear ;
	bulkshear.push_back(TAG_BULK_MODULUS) ;
	bulkshear.push_back(TAG_SHEAR_MODULUS) ;
	double totalArea ;
	if(self2d)
	{
		for(size_t i = 0 ; i < source.size() ; i++)
		{
			if(source[i]->getBehaviour() && self2d->getPrimitive()->in(source[i]->getCenter()))
			{
				time_d = time_d || source[i]->getBehaviour()->timeDependent() ;
				space_d = space_d || source[i]->getBehaviour()->spaceDependent() ;
				
				if(source[i]->getBehaviour()->type != VOID_BEHAVIOUR)
				{
					Material thismat = source[i]->getBehaviour()->toMaterial() ;
					thismat(TAG_VOLUME,source[i]->area()) ;
					if(mat.nPhases() == 0)
						mat = mat + thismat ;
					else 
					{
						bool combined = false ;
						for(size_t j = 0 ; j < mat.nPhases() ; j++)
							combined = (combined || mat.child(j).combine(thismat,bulkshear,TAG_VOLUME)) ;
						if(!combined)
							mat = mat + thismat ;
					}
					totalArea += source[i]->area() ;
				}
				else
				{
					Material thismat(Properties(TAG_VOLUME, source[i]->area())) ;
					thismat(TAG_BULK_MODULUS, 1e-9) ;
					thismat(TAG_SHEAR_MODULUS, 1e-9) ;
					if(mat.nPhases() == 0)
						mat = mat + thismat ;
					else 
					{
						bool combined = false ;
						for(size_t j = 0 ; j < mat.nPhases() ; j++)
							combined = (combined || mat.child(j).combine(thismat,bulkshear,TAG_VOLUME)) ;
						if(!combined)
							mat = mat + thismat ;
					}
					totalArea += source[i]->area() ;
				}
			}
		}
/*	} else {
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
					Material thismat(param_gp/area_gp) ;
					thismat.push_back(Properties(TAG_VOLUME,source3d[i]->volume())) ;
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
					totalArea += source3d[i]->volume() ;
				}
				else
				{
					Material thismat(Properties(TAG_VOLUME, source3d[i]->volume())) ;
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
					totalArea += source3d[i]->volume() ;
				}
			}
		}*/

	}

	Material hom = homogenize(mat) ;
	equivalent = getEquivalentBehaviour(hom) ;

/*	GeneralConverter fraction(TAG_VOLUME_FRACTION) ;
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


Material HomogeneisedBehaviour::homogenize(Material mat)
{
	return mat ;
}

Form * HomogeneisedBehaviour::getEquivalentBehaviour(Material mat)
{
	return new VoidForm() ;
}

