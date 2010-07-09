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
#include "diffusion.h"
#include "stiffness.h"
#include "stiffness_with_imposed_deformation.h"
#include "void_form.h"
#include "homogenization/properties_base.h"
#include "homogenization/scheme_base.h"
#include "homogenization/elastic_homogenization.h"
#include "homogenization/expansion_homogenization.h"
#include "homogenization/converter.h"
#include "../features/features.h"
#include "../features/inclusion.h"
#include "../features/sample.h"
#include "../geometry/geometry_base.h" 


using namespace Mu ;

HomogeneisedBehaviour::HomogeneisedBehaviour(Mesh<DelaunayTriangle, DelaunayTreeItem> * mesh2d, DelaunayTriangle * self) : LinearForm(Matrix(), false, false, 2) , mesh2d(mesh2d), self2d(self), mesh3d(NULL), self3d(NULL)
{
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
	
	Material hom ;
	
	Material matrix = tri->getBehaviour()->toMaterial() ;
	double fmat = tri->area() ;
	for(size_t i = 0 ; i < feats.size() ; i++)
		fmat -= feats[i]->area() ;
	matrix(TAG_VOLUME,fmat) ;

	hom = hom + matrix ;
	
	for(size_t i = 0 ; i < feats.size() ; i++)
	{
		Material inc = feats[i]->getBehaviour()->toMaterial() ;
		inc(TAG_VOLUME,feats[i]->area()) ;
	}



	Material eq = homogenize(hom) ;
	equivalent = getEquivalentBehaviour(eq) ;

	v.push_back(XI);
	v.push_back(ETA);

}


HomogeneisedBehaviour::HomogeneisedBehaviour(Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * mesh3d, DelaunayTetrahedron * self) : LinearForm(Matrix(), false, false, 3), mesh2d(NULL), self2d(NULL), mesh3d(mesh3d), self3d(self)
{
	source3d = mesh3d->getConflictingElements(self->getPrimitive()) ;

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
	
	v.push_back(XI);
	v.push_back(ETA);
	v.push_back(ZETA);
} ;

HomogeneisedBehaviour::~HomogeneisedBehaviour() { } ;


void HomogeneisedBehaviour::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	equivalent->apply(p_i,p_j,gp,Jinv,ret,vm) ;
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
	double totalArea = 0. ;
	std::cout << self2d->area() << std::endl ;
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
						{
							Material tmp = mat.child(j) ;
							combined = (combined || tmp.combine(thismat,bulkshear,TAG_VOLUME)) ;
							if(combined)
							{
								mat = mat - (int) j ;
								mat = mat + tmp ;
								break ;
							}
						}
						if(!combined)
							mat = mat + thismat ;
					}
					totalArea += source[i]->area() ;
				}
				else
				{
					size_t bp = source[i]->getBoundingPoints().size() ;
					Point a = source[i]->getBoundingPoint(0) ;
					Point b = source[i]->getBoundingPoint(bp/3) ;
					Point c = source[i]->getBoundingPoint(bp*2/3) ;
					Material thismat(Properties(TAG_VOLUME, source[i]->area())) ;
					thismat(TAG_BULK_MODULUS, 1e-9) ;
					thismat(TAG_SHEAR_MODULUS, 1e-9) ;
					if(mat.nPhases() == 0)
						mat = mat + thismat ;
					else 
					{
//						thismat.print() ;
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

	std::cout << totalArea << std::endl ;

	Material hom = homogenize(mat) ;
	equivalent = getEquivalentBehaviour(hom) ;

}


Material HomogeneisedBehaviour::homogenize(Material mat)
{
	if(mat.nPhases() == 1)
		return mat.child(0) ;

	mat.makeFraction(true) ;
	if(mat.nPhases() == 2)
	{
		if(mat.child(0).val(TAG_BULK_MODULUS,-1) < 1e-6)
		{
			Material tmp = mat.child(0) ;
			mat = mat - 0 ;
			mat = mat + tmp ;
		}
		mat.build(new MoriTanaka(), false) ;

		bool imp = false ;
		for(size_t i = 0 ; i < mat.nPhases() ; i++)
			 imp = imp || (mat.child(i).getIndex(TAG_IMPOSED_STRAIN,-1) > -1) ;

		if(imp)
		{
			for(size_t i = 0 ; i < mat.nPhases() ; i++)
			{
				double def = mat.child(i).val(TAG_IMPOSED_STRAIN,-1) ;
				mat.add(TAG_EXPANSION_COEFFICIENT,def,i) ;
			}
			mat.build(new HobbsScheme(), false) ;
			mat.add(TAG_IMPOSED_STRAIN,mat(TAG_EXPANSION_COEFFICIENT)) ;
		}

		return mat ;
	}

	mat.build(new GeneralizedSelfConsistent(), false) ;

	return mat ;
}

Form * HomogeneisedBehaviour::getEquivalentBehaviour(Material mat)
{
	bool stiff = false ;
	bool diff = false ;
	bool frac = false ;
	bool imp = false ;

	if(mat.getIndex(TAG_DIFFUSION_COEFFICIENT,-1) > -1)
	{
		if(self2d)
		{
			param = Matrix(2,2) ;
			param[0][0] = mat(TAG_DIFFUSION_COEFFICIENT) ;
			param[1][1] = mat(TAG_DIFFUSION_COEFFICIENT) ;
		}
		if(self3d)
		{
			param = Matrix(3,3) ;
			param[0][0] = mat(TAG_DIFFUSION_COEFFICIENT) ;
			param[1][1] = mat(TAG_DIFFUSION_COEFFICIENT) ;
			param[2][2] = mat(TAG_DIFFUSION_COEFFICIENT) ;
		}
		diff = true ;
		return new Diffusion(param) ;
	}

	if(mat.getIndex(TAG_YOUNG_MODULUS,-1) > -1 && mat.getIndex(TAG_POISSON_RATIO,-1) > -1)
	{
		if(self2d)
			param = cauchyGreen(std::make_pair(mat(TAG_YOUNG_MODULUS),mat(TAG_POISSON_RATIO)),true,SPACE_TWO_DIMENSIONAL) ;
		if(self3d)
			param = cauchyGreen(std::make_pair(mat(TAG_YOUNG_MODULUS),mat(TAG_POISSON_RATIO)),true,SPACE_THREE_DIMENSIONAL) ;
		stiff = true ;
	} else {
		if(mat.getIndex(TAG_BULK_MODULUS,-1) > -1 && mat.getIndex(TAG_SHEAR_MODULUS,-1) > -1)
		{
			if(mat(TAG_BULK_MODULUS) < 1e-6)
				mat.print() ;
			if(self2d)
				param = cauchyGreen(std::make_pair(mat(TAG_BULK_MODULUS),mat(TAG_SHEAR_MODULUS)),false,SPACE_TWO_DIMENSIONAL) ;
			if(self3d)
				param = cauchyGreen(std::make_pair(mat(TAG_BULK_MODULUS),mat(TAG_SHEAR_MODULUS)),false,SPACE_THREE_DIMENSIONAL) ;
			stiff = true ;
		} else {
			std::cout << "warning: no mechanical behaviour detected!" << std::endl ;
		}
	}

	Vector imposed(3) ;
	if(self3d)
		imposed.resize(6) ;
	if(mat.getIndex(TAG_IMPOSED_STRAIN,-1) > -1)
	{
		imposed[0] = mat(TAG_IMPOSED_STRAIN) ;
		imposed[1] = mat(TAG_IMPOSED_STRAIN) ;
		if(imposed.size() == 6)
			imposed[2] = mat(TAG_IMPOSED_STRAIN) ;
		imp = true ;
	}

	if(stiff)
	{
		if(imp)
			return new StiffnessWithImposedDeformation(param,imposed) ;

		return new Stiffness(param) ;
	}




	return new VoidForm() ;
}

