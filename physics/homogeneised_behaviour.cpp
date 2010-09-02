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
#include "stiffness_and_fracture.h"
#include "stiffness_with_imposed_deformation.h"
#include "stiffness_with_imposed_deformation_and_fracture.h"
#include "void_form.h"
#include "fracturecriteria/fracturecriterion.h"
#include "fracturecriteria/vonmises.h"
#include "fracturecriteria/mohrcoulomb.h"
#include "fracturecriteria/ruptureenergy.h"
#include "fracturecriteria/maxstrain.h"
#include "fracturecriteria/limitstrains.h"
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

HomogeneisedBehaviour::HomogeneisedBehaviour(Mesh<DelaunayTriangle, DelaunayTreeItem> * mesh2d, DelaunayTriangle * self) : LinearForm(Matrix(), true, false, 2) , mesh2d(mesh2d), self2d(self), mesh3d(NULL), self3d(NULL)
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

HomogeneisedBehaviour::HomogeneisedBehaviour(std::vector<Feature *> feats, DelaunayTriangle * self) : LinearForm(Matrix(), true, false, 2), self2d(self), mesh3d(NULL), self3d(NULL)
{
	std::vector<Point> corner = self->getSamplingBoundingPoints(3) ;

	TriangularInclusion * tri = new TriangularInclusion(corner[0],corner[1],corner[2]) ;
	tri->setBehaviour(self->getBehaviour()) ;

	for(size_t i = 0 ; i < feats.size() ; i++)
            ft.push_back(feats[i]) ;
	
	Material hom ;
	
	Material matrix = tri->getBehaviour()->toMaterial() ;
	double fmat = tri->area() ;
	for(size_t i = 0 ; i < feats.size() ; i++)
		fmat -= feats[i]->area() ;
	matrix(TAG_VOLUME,fmat) ;

	hom = hom + matrix ;
	
	
	std::vector<Tag> bulkshear ;
	bulkshear.push_back(TAG_BULK_MODULUS) ;
	bulkshear.push_back(TAG_SHEAR_MODULUS) ;
	for(size_t i = 0 ; i < feats.size() ; i++)
	{
		Material inc = feats[i]->getBehaviour()->toMaterial() ;
		inc(TAG_VOLUME, feats[i]->area()) ;
		bool combined = false ;
		for(size_t j = 0 ; j < hom.nPhases() ; j++)
		{
			Material tmp = hom.child(j) ;
			combined = (combined || tmp.combine(inc,bulkshear,TAG_VOLUME)) ;
			if(combined)
			{
				hom = hom - (int) j ;
				hom = hom + tmp ;
				break ;
			}
		}
		if(!combined)
			hom = hom + inc ;
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

        std::vector<Point> corner = self2d->getSamplingBoundingPoints(3) ;

        TriangularInclusion * tri = new TriangularInclusion(corner[0],corner[1],corner[2]) ;
        tri->setBehaviour(self2d->getBehaviour()) ;

        Material hom ;

        Material matrix = tri->getBehaviour()->toMaterial() ;
        double fmat = tri->area() ;
        for(size_t i = 0 ; i < ft.size() ; i++)
                fmat -= ft[i]->area() ;
        matrix(TAG_VOLUME,fmat) ;

        hom = hom + matrix ;


        std::vector<Tag> bulkshear ;
        bulkshear.push_back(TAG_BULK_MODULUS) ;
        bulkshear.push_back(TAG_SHEAR_MODULUS) ;
        for(size_t i = 0 ; i < ft.size() ; i++)
        {
                Material inc = ft[i]->getBehaviour()->toMaterial() ;
                inc(TAG_VOLUME, ft[i]->area()) ;
                bool combined = false ;
                for(size_t j = 0 ; j < hom.nPhases() ; j++)
                {
                        Material tmp = hom.child(j) ;
                        combined = (combined || tmp.combine(inc,bulkshear,TAG_VOLUME)) ;
                        if(combined)
                        {
                                hom = hom - (int) j ;
                                hom = hom + tmp ;
                                break ;
                        }
                }
                if(!combined)
                        hom = hom + inc ;
        }



        Material eq = homogenize(hom) ;
        equivalent = getEquivalentBehaviour(eq) ;

        equivalent->step(timestep, currentState) ;

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

std::vector<BoundaryCondition * > HomogeneisedBehaviour::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
        return equivalent->getBoundaryConditions(s,id,p_i,gp,Jinv) ;
}

void HomogeneisedBehaviour::homogenize()
{
	Material mat ;
	std::vector<Tag> bulkshear ;
	bulkshear.push_back(TAG_BULK_MODULUS) ;
	bulkshear.push_back(TAG_SHEAR_MODULUS) ;
	double totalArea = 0. ;
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
	} else {
		for(size_t i = 0 ; i < source3d.size() ; i++)
		{
			if(source3d[i]->getBehaviour() && self3d->getPrimitive()->in(source3d[i]->getCenter()))
			{
				time_d = time_d || source3d[i]->getBehaviour()->timeDependent() ;
				space_d = space_d || source3d[i]->getBehaviour()->spaceDependent() ;
				
				if(source3d[i]->getBehaviour()->type != VOID_BEHAVIOUR)
				{
					Material thismat = source3d[i]->getBehaviour()->toMaterial() ;
					thismat(TAG_VOLUME,source3d[i]->volume()) ;
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
					totalArea += source3d[i]->volume() ;
				}
				else
				{
					size_t bp = source3d[i]->getBoundingPoints().size() ;
					Point a = source3d[i]->getBoundingPoint(0) ;
					Point b = source3d[i]->getBoundingPoint(bp/3) ;
					Point c = source3d[i]->getBoundingPoint(bp*2/3) ;
					Material thismat(Properties(TAG_VOLUME, source3d[i]->volume())) ;
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
					totalArea += source3d[i]->volume() ;
				}	
			}
		}

	}

	Material hom = homogenize(mat) ;
	equivalent = getEquivalentBehaviour(hom) ;

}


Material HomogeneisedBehaviour::homogenize(Material mat)
{
	if(mat.nPhases() == 1)
		return mat.child(0) ;

	std::vector<Tag> fracture ;
	fracture.push_back(TAG_MAX_STRAIN) ;
	fracture.push_back(TAG_MAX_COMPRESSIVE_STRAIN) ;
	fracture.push_back(TAG_MAX_TENSILE_STRAIN) ;
	fracture.push_back(TAG_MAX_STRESS) ;
	fracture.push_back(TAG_MAX_COMPRESSIVE_STRESS) ;
	fracture.push_back(TAG_MAX_TENSILE_STRESS) ;
	fracture.push_back(TAG_RUPTURE_ENERGY) ;
	
	for(size_t i = 0 ; i < fracture.size() ; i++)
	{
		int k = -1 ;
		double val = 0. ;
		for(size_t j = 0 ; j < mat.nPhases() ; j++)
		{
			if(mat.child(j).getIndex(fracture[i],-1) > -1)
			{
				double vj = std::abs(mat.child(j).val(fracture[i],-1)) ;
				if(k == -1)
				{
					val = vj ;
					k = j ;
				} else {
					if(vj > val)
					{
						val = vj ;
						k = j ;
					}
				}
			}
		}
		if(k+1 > 0)
		{
			mat(fracture[i],mat.child(k).val(fracture[i],-1)) ;
		}
	}

	
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
			if(self2d)
				param = cauchyGreen(std::make_pair(mat(TAG_BULK_MODULUS),mat(TAG_SHEAR_MODULUS)),false,SPACE_TWO_DIMENSIONAL) ;
			if(self3d)
				param = cauchyGreen(std::make_pair(mat(TAG_BULK_MODULUS),mat(TAG_SHEAR_MODULUS)),false,SPACE_THREE_DIMENSIONAL) ;
			stiff = true ;
			if(mat(TAG_BULK_MODULUS) < 1e-6)
				stiff = false ;
		} else {
			std::cout << "warning: no mechanical behaviour detected!" << std::endl ;
		}
	}
	
	FractureCriterion * crit ;
	
	if(mat.getIndex(TAG_MAX_STRAIN,-1) > -1 ||
	   mat.getIndex(TAG_MAX_TENSILE_STRAIN,-1) > -1 ||
	   mat.getIndex(TAG_MAX_COMPRESSIVE_STRAIN,-1) > -1 ||
	   mat.getIndex(TAG_MAX_STRESS,-1) > -1 ||
	   mat.getIndex(TAG_MAX_TENSILE_STRESS,-1) > -1 ||
	   mat.getIndex(TAG_MAX_COMPRESSIVE_STRESS,-1) > -1 ||
	   mat.getIndex(TAG_RUPTURE_ENERGY,-1) > -1)
	{
		frac = true ;
		if(mat.getIndex(TAG_MAX_TENSILE_STRAIN,-1) > -1 && mat.getIndex(TAG_MAX_COMPRESSIVE_STRAIN,-1) > -1)
			crit = new LimitStrains(mat(TAG_MAX_COMPRESSIVE_STRAIN),mat(TAG_MAX_TENSILE_STRAIN)) ;
		if(mat.getIndex(TAG_MAX_STRAIN,-1) > -1)
			crit = new MaximumStrain(mat(TAG_MAX_STRAIN)) ;
		if(mat.getIndex(TAG_RUPTURE_ENERGY,-1) > -1)
			crit = new RuptureEnergy(mat(TAG_RUPTURE_ENERGY)) ;
		if(mat.getIndex(TAG_MAX_TENSILE_STRESS,-1) > -1 && mat.getIndex(TAG_MAX_COMPRESSIVE_STRESS,-1) > -1)
			crit = new MohrCoulomb(mat(TAG_MAX_TENSILE_STRESS),mat(TAG_MAX_COMPRESSIVE_STRESS)) ;
		if(mat.getIndex(TAG_MAX_STRESS,-1) > -1)
			crit = new VonMises(mat(TAG_MAX_STRESS)) ;
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
/*		if(imp && frac)
                        return new StiffnessWithImposedDeformationAndFracture(param,imposed,crit) ;*/

		if(imp)
                {
                    return new StiffnessWithImposedDeformation(param,imposed) ;
                }
		
/*		if(frac)
                        return new StiffnessAndFracture(param,crit) ;*/

		return new Stiffness(param) ;
	}




	return new VoidForm() ;
}

Vector HomogeneisedBehaviour::getImposedStress(const Point & p) const
{
    return equivalent->getImposedStress(p) ;
}
