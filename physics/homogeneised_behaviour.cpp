//
// Description: 
//
//
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
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
#include "homogenization/homogenization_base.h"
#include "../features/features.h"
#include "../features/inclusion.h"
#include "../features/sample.h"
#include "../geometry/geometry_base.h" 


using namespace Mu ;

HomogeneisedBehaviour::HomogeneisedBehaviour(FeatureTree * mesh, DelaunayTriangle * self) : LinearForm(Matrix(), true, false, 2) , mesh(mesh), self2d(self), self3d(NULL), equivalent(NULL)
{
	std::vector<DelaunayTriangle * > feats = mesh->getElements2D(self->getPrimitive()) ;
	Material hom ;
	for(size_t i = 0 ; i < feats.size() ; i++)
	{
		
		Material inc = feats[i]->getBehaviour()->toMaterial(Point(1./3., 1./3.)) ;
		
		inc.setProperties(P_VOLUME, feats[i]->area()) ;
		hom.addPhase(inc) ;
		
		hom.mergePhase() ;
		
	}
	Material eq = homogenize(hom) ;
	equivalent = getEquivalentBehaviour(eq) ;
	
	v.push_back(XI);
	v.push_back(ETA);
	

	reverted = false ;
} ;

HomogeneisedBehaviour::HomogeneisedBehaviour(std::vector<Feature *> feats, DelaunayTriangle * self) : LinearForm(Matrix(), true, false, 2), self2d(self), mesh(NULL), self3d(NULL), equivalent(NULL)
{

	std::vector<Point> corner = self->getSamplingBoundingPoints(0) ;
	
	if(self->getBehaviour())
		base = self->getBehaviour()->toMaterial(Point(1./3., 1./3.)) ;
	else if(!feats.empty() && feats[0]->getBehaviour())
		base = feats[0]->getBehaviour()->toMaterial(self->inLocalCoordinates(feats[0]->getCenter())) ;
	else
	{
		std::cout << "cannot do anything." << std::endl ;
		exit(0) ;
	}
	
	
	TriangularInclusion tri(corner[0],corner[1],corner[2]) ;
	tri.setBehaviour(self->getBehaviour()->getCopy()) ;

	for(size_t i = 0 ; i < feats.size() ; i++)
            ft.push_back(feats[i]) ;
	Material hom ;
	
	Material matrix ;
	for(int i = 0 ; i < base.sizeProperties() ; i++)
	{
		matrix.setProperties(base.getProperties(i)) ;
	}
	double fmat = tri.area() ;
	for(size_t i = 0 ; i < feats.size() ; i++)
		fmat -= feats[i]->area() ;
	matrix.setProperties(P_VOLUME, fmat) ;
	hom.addPhase(matrix) ;
	for(size_t i = 0 ; i < feats.size() ; i++)
	{
		Material inc = feats[i]->getBehaviour()->toMaterial(self->inLocalCoordinates(feats[i]->getCenter())) ;
		inc.setProperties(P_VOLUME, feats[i]->area()) ;
		hom.addPhase(inc) ;
		hom.mergePhase() ;
	}
	Material eq = homogenize(hom) ;

	equivalent = getEquivalentBehaviour(eq) ;
	v.push_back(XI);
	v.push_back(ETA);
	
	reverted = false ;
}


HomogeneisedBehaviour::HomogeneisedBehaviour(FeatureTree * mesh, DelaunayTetrahedron * self) : LinearForm(Matrix(), false, false, 3), mesh(mesh), self2d(NULL), self3d(self), equivalent(NULL)
{
	source3d = mesh->getElements3D(self->getPrimitive()) ;

	double totalVolume = 0 ;
	for(size_t i = 0 ; i < source3d.size() ; i++)
	{
		if(source3d[i]->getBehaviour() && source3d[i]->getBehaviour()->type != VOID_BEHAVIOUR)
		{
			param = Matrix(source3d[i]->getBehaviour()->param.numRows(), source3d[i]->getBehaviour()->param.numCols()) ;
			break ;
		}
	}

//	homogenize() ;
	
	v.push_back(XI);
	v.push_back(ETA);
	v.push_back(ZETA);
} ;

HomogeneisedBehaviour::~HomogeneisedBehaviour() 
{ 
	delete equivalent ; 
} ;


void HomogeneisedBehaviour::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	equivalent->apply(p_i,p_j,gp,Jinv,ret,vm) ;
}

void HomogeneisedBehaviour::step(double timestep, ElementState & currentState)
{
	if(type == VOID_BEHAVIOUR)
		return ;

	if(reverted)
	{
		return ;
	}
	
	bool revert = false ;

	std::vector<Point> corner = self2d->getSamplingBoundingPoints(0) ;

	TriangularInclusion tri(corner[0],corner[1],corner[2]) ;

	Material hom ;

	Material matrix ;
	for(int i = 0 ; i < base.sizeProperties() ; i++)
	{
		matrix.setProperties(base.getProperties(i)) ;
	}
        double fmat = tri.area() ;
        for(size_t i = 0 ; i < ft.size() ; i++)
        {
        	revert |= dynamic_cast<Triangle *>(&tri)->intersects(dynamic_cast<Geometry *>(ft[i])) ;
                fmat -= ft[i]->area() ;
	}
	
	if(revert)
	{
		ft.clear() ;
		delete equivalent ;
		equivalent = getEquivalentBehaviour(base) ;
		reverted = true ;
		return ;
	}
	
	matrix.setProperties(P_VOLUME, fmat) ;
	hom.addPhase(matrix) ;

	if(self2d)
	{
		for(size_t i = 0 ; i < ft.size() ; i++)
		{
			Material inc = ft[i]->getBehaviour()->toMaterial(self2d->inLocalCoordinates(ft[i]->getCenter())) ;
			inc.setProperties(P_VOLUME, ft[i]->area()) ;
			hom.addPhase(inc) ;
			hom.mergePhase() ;
		}
	}
	else
	{
		for(size_t i = 0 ; i < ft.size() ; i++)
		{
			Material inc = ft[i]->getBehaviour()->toMaterial(self3d->inLocalCoordinates(ft[i]->getCenter())) ;
			inc.setProperties(P_VOLUME, ft[i]->area()) ;
			hom.addPhase(inc) ;
			hom.mergePhase() ;
		}
	}

	Material eq = homogenize(hom) ;
	if(equivalent)
		delete equivalent ;
	equivalent = getEquivalentBehaviour(eq) ;

	if(equivalent->timeDependent())
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
					Material thismat = source[i]->getBehaviour()->toMaterial(self2d->inLocalCoordinates(source[i]->getCenter())) ;
					thismat.setProperties(P_VOLUME,source[i]->area()) ;
					mat.addPhase(thismat) ;
					mat.mergePhase() ;
					totalArea += source[i]->area() ;
				}
				else
				{
					Material thismat ;
					thismat.setProperties(P_BULK_MODULUS, 1e-9) ;
					thismat.setProperties(P_SHEAR_MODULUS, 1e-9) ;
					thismat.setProperties(P_VOLUME,source[i]->area()) ;
					mat.addPhase(thismat) ;
					mat.mergePhase() ;
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
					Material thismat = source3d[i]->getBehaviour()->toMaterial(self3d->inLocalCoordinates(source3d[i]->getCenter())) ;
					thismat.setProperties(P_VOLUME,source3d[i]->volume()) ;
					mat.addPhase(thismat) ;
					mat.mergePhase() ;
					totalArea += source3d[i]->volume() ;
				}
				else
				{
					Material thismat ;
					thismat.setProperties(P_BULK_MODULUS, 1e-9) ;
					thismat.setProperties(P_SHEAR_MODULUS, 1e-9) ;
					thismat.setProperties(P_VOLUME,source3d[i]->volume()) ;
					mat.addPhase(thismat) ;
					mat.mergePhase() ;
					totalArea += source3d[i]->volume() ;
				}	
			}
		}

	}

	Material hom = homogenize(mat) ;
	delete equivalent ;
	FractureCriterion * frac ;
	if(self2d)
	    frac = self2d->getBehaviour()->getFractureCriterion()->getCopy() ;
	if(self3d)
	    frac = self3d->getBehaviour()->getFractureCriterion()->getCopy() ;

	equivalent = getEquivalentBehaviour(hom, frac) ;
	
}


Material HomogeneisedBehaviour::homogenize(Material mat)
{
	if(mat.sizePhase() == 1)
	{
		return mat.getPhase(0) ;
	}

	if(mat.sizePhase() == 2)
	{
		if(mat.getPhase(0).valProperties(P_BULK_MODULUS) < 1e-6)
		{
			Material tmp = mat.getPhase(0) ;
			mat.removePhase(0) ;
			mat.addPhase(tmp) ;
		}
		mat.homogenize(ELASTICITY_MORI_TANAKA) ;

		bool imp = false ;
		for(int i = 0 ; i < mat.sizePhase() ; i++)
		{
			imp |= mat.getPhase(i).hasProperties(P_EXPANSION_COEFFICIENT) ;
		}
		if(imp)
		{
			mat.homogenize(EXPANSION_HOBBS) ;
		}


		return mat ;
	}

	mat.homogenize(ELASTICITY_GENERALIZED_SELF_CONSISTENT) ;
	return mat ;
}

Form * HomogeneisedBehaviour::getEquivalentBehaviour(Material mat, FractureCriterion * frac)
{
	bool stiff = false ;
//	bool diff = false ;
	bool imp = false ;

	if(mat.hasProperties(P_YOUNG_MODULUS) && mat.hasProperties(P_POISSON_RATIO))
	{
		double E = mat.valProperties(P_YOUNG_MODULUS) ;
		double nu = mat.valProperties(P_POISSON_RATIO) ;
		if(self2d)
		{
			param = Material::cauchyGreen(std::make_pair(E, nu),true,SPACE_TWO_DIMENSIONAL) ;
		}
		if(self3d)
		{
			param = Material::cauchyGreen(std::make_pair(E, nu),true,SPACE_THREE_DIMENSIONAL) ;
		}
		stiff = true ;
	}
	else
	{
		if(mat.hasProperties(P_BULK_MODULUS) && mat.hasProperties(P_SHEAR_MODULUS))
		{
			double k = mat.valProperties(P_BULK_MODULUS) ;
			double mu = mat.valProperties(P_SHEAR_MODULUS) ;
			if(self2d)
			{
				param = Material::cauchyGreen(std::make_pair(k, mu),false,SPACE_TWO_DIMENSIONAL) ;
			}
			if(self3d)
			{
				param = Material::cauchyGreen(std::make_pair(k, mu),false,SPACE_THREE_DIMENSIONAL) ;
			}
			stiff = ( k > 1e-6) ;
		}
		else
		{
			std::cout << "warning: no mechanical behaviour detected!" << std::endl ;
		}
	}
	
	Vector def(3) ;
	if(mat.hasProperties(P_EXPANSION_COEFFICIENT))
	{
		double alpha = mat.valProperties(P_EXPANSION_COEFFICIENT) ;
		if(self2d)
		{
			def[0] = alpha ;
			def[1] = alpha ;
		}
		if(self3d)
		{
			def.resize(6) ;
			def[0] = alpha ;
			def[1] = alpha ;
			def[2] = alpha ;
		}
		imp = true ;
	}
	
	if(stiff)
	{
		if(imp)
		{
		    if(frac)
		    {
			return new StiffnessWithImposedDeformationAndFracture(param, def, frac) ;
		    }
			return new StiffnessWithImposedDeformation(param, def) ;
		}
		
		if(frac)
		    return new StiffnessAndFracture(param,frac) ;
		return new Stiffness(param) ;
	}

	return new VoidForm() ;
}

Vector HomogeneisedBehaviour::getImposedStress(const Point & p) const
{
    return equivalent->getImposedStress(p) ;
}

Matrix HomogeneisedBehaviour::getTensor(const Point & p) const
{
	return equivalent->getTensor(p) ;
}

