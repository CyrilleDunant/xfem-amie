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
#include <typeinfo>

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
#include "homogenization/phase.h"
#include "homogenization/composite.h"
#include "../features/features.h"
#include "../features/inclusion.h"
#include "../features/sample.h"
#include "../geometry/geometry_base.h"


using namespace Mu ;

HomogeneisedBehaviour::HomogeneisedBehaviour( FeatureTree *mesh, DelaunayTriangle *self ) : LinearForm( Matrix(), true, false, 2 ) , mesh( mesh ), self2d( self ), self3d( nullptr ), equivalent( nullptr )
{
	if( dynamic_cast<HomogeneisedBehaviour *>( self->getBehaviour() ) )
		original = static_cast<HomogeneisedBehaviour *>( self->getBehaviour() )->getOriginalBehaviour()->getCopy() ;
	else
		original = self->getBehaviour()->getCopy() ;

	GeneralizedSelfConsistentComposite composite( mesh->getElements2D( self->getPrimitive() ) ) ;
	equivalent = composite.getBehaviour() ;
// 	Vector alpha = static_cast<StiffnessWithImposedDeformation *>(equivalent)->imposed ;
// 	std::cout << alpha[0] << std::endl ;

	v.push_back( XI );
	v.push_back( ETA );

	reverted = false ;
} ;

HomogeneisedBehaviour::HomogeneisedBehaviour( std::vector<Feature *> feats, DelaunayTriangle *self ) : LinearForm( Matrix(), true, false, 2 ), self2d( self ), mesh( nullptr ), self3d( nullptr ), equivalent( nullptr )
{
	if( dynamic_cast<HomogeneisedBehaviour *>( self->getBehaviour() ) )
		original = static_cast<HomogeneisedBehaviour *>( self->getBehaviour() )->getOriginalBehaviour()->getCopy() ;
	else
		original = self->getBehaviour()->getCopy() ;


	VoigtMatrixMultiInclusionComposite composite( self, feats ) ;
	equivalent = composite.getBehaviour() ;
/*	equivalent = composite.inclusions[0].getBehaviour() ;
	equivalent->getTensor(Point(0,0)).print() ;*/
	
// 	if(original->getFractureCriterion())
// 	{
// 	    FractureCriterion * frac = original->getFractureCriterion() ;
// 	    Matrix C = equivalent->getTensor(Point(1./3,1./3,1./3)) ;
// 	    Vector alpha = static_cast<StiffnessWithImposedDeformation *>(equivalent)->imposed ;
// 	    equivalent = new StiffnessWithImposedDeformationAndFracture(C,alpha,frac->getCopy()) ;
// 	}
// 	Vector alpha = static_cast<StiffnessWithImposedDeformation *>(equivalent)->imposed ;
// 	Vector alpha = static_cast<StiffnessWithImposedDeformation *>(equivalent)->imposed ;
// 	std::cout << alpha[0] << std::endl ;

	v.push_back( XI );
	v.push_back( ETA );

	reverted = false ;
	
}

HomogeneisedBehaviour::HomogeneisedBehaviour( std::vector<Feature *> feats, DelaunayTetrahedron *self ) : LinearForm( Matrix(), true, false, 3 ), self3d( self ), mesh( nullptr ), self2d( nullptr ), equivalent( nullptr )
{
	
	
	if( dynamic_cast<HomogeneisedBehaviour *>( self->getBehaviour() ) )
		original = static_cast<HomogeneisedBehaviour *>( self->getBehaviour() )->getOriginalBehaviour()->getCopy() ;
	else
		original = self->getBehaviour()->getCopy() ;
	
	VoigtMatrixMultiInclusionComposite composite( self, feats ) ;
	equivalent = composite.inclusions[0].getBehaviour() ;
	equivalent->getTensor(Point(0,0)).print() ;
	
/*	if(original->getFractureCriterion())
	{
		FractureCriterion * frac = original->getFractureCriterion() ;
		Matrix C = equivalent->getTensor(Point(1./3,1./3,1./3)) ;
		Vector alpha = static_cast<StiffnessWithImposedDeformation *>(equivalent)->imposed ;
		std::cout << alpha[0] << std::endl ;
		equivalent = new StiffnessWithImposedDeformationAndFracture(C,alpha,frac->getCopy()) ;
	}*/
	
	v.push_back( XI );
	v.push_back( ETA );
	v.push_back( ZETA );
	
	reverted = false ;
	
}

void HomogeneisedBehaviour::updateEquivalentBehaviour(std::vector<Feature *> feats, DelaunayTriangle * self) 
{

	VoigtMatrixMultiInclusionComposite composite( self, feats ) ;
/*	FractureCriterion * frac = equivalent->getFractureCriterion() ;
	FractureCriterion * fracCopy = nullptr ;
	if(frac)
	{
		fracCopy = frac->getCopy() ;
	}*/
	delete equivalent ;
	equivalent = composite.getBehaviour() ;
/*	if(fracCopy)
	{
		Matrix C = equivalent->getTensor(Point(1./3,1./3,1./3)) ;
		Vector alpha = static_cast<StiffnessWithImposedDeformation *>(equivalent)->imposed ;
		delete equivalent ;
		equivalent = new StiffnessWithImposedDeformationAndFracture(C,alpha,fracCopy) ;
	}*/
	
	
}

void HomogeneisedBehaviour::updateEquivalentBehaviour(std::vector<Feature *> feats, DelaunayTetrahedron * self) 
{
	VoigtMatrixMultiInclusionComposite composite( self, feats ) ;
/*	FractureCriterion * frac = equivalent->getFractureCriterion() ;
	FractureCriterion * fracCopy = nullptr ;
	if(frac)
	{
		fracCopy = frac->getCopy() ;
	}*/
	delete equivalent ;
	equivalent = composite.getBehaviour() ;
/*	if(fracCopy)
	{
		Matrix C = equivalent->getTensor(Point(1./3,1./3,1./3)) ;
		Vector alpha = static_cast<StiffnessWithImposedDeformation *>(equivalent)->imposed ;
		delete equivalent ;
		equivalent = new StiffnessWithImposedDeformationAndFracture(C,alpha,fracCopy) ;
	}*/
	
	
}


HomogeneisedBehaviour::HomogeneisedBehaviour( FeatureTree *mesh, DelaunayTetrahedron *self ) : LinearForm( Matrix(), false, false, 3 ), mesh( mesh ), self2d( nullptr ), self3d( self ), equivalent( nullptr )
{

	original = self->getBehaviour() ;

	if( dynamic_cast<HomogeneisedBehaviour *>( original ) )
		original = static_cast<HomogeneisedBehaviour *>( original )->getOriginalBehaviour() ;

	GeneralizedSelfConsistentComposite composite( mesh->getElements3D( self->getPrimitive() ) ) ;
	equivalent = composite.getBehaviour() ;

	v.push_back( XI );
	v.push_back( ETA );
	v.push_back( ZETA );

	reverted = false ;

} ;

HomogeneisedBehaviour::~HomogeneisedBehaviour()
{
//	delete equivalent ;
} ;


void HomogeneisedBehaviour::apply( const Function &p_i, const Function &p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix &ret, VirtualMachine *vm ) const
{
// 	param.print() ;
	equivalent->apply( p_i, p_j, gp, Jinv, ret, vm ) ;
}

void HomogeneisedBehaviour::step( double timestep, ElementState &currentState, double maxscore )
{
	if( type == VOID_BEHAVIOUR )
		return ;

	if( reverted )
	{
		return ;
	}

	bool revert = false ;

	if( revert )
	{
		ft.clear() ;
		delete equivalent ;
		equivalent = original ;
		reverted = true ;
		return ;
	}


}


bool HomogeneisedBehaviour::fractured() const
{
	return false ;
}

Form * HomogeneisedBehaviour::getCopy() const
{
	Form * copy =  equivalent->getCopy() ;
	
	if(getExtra2dMeshes())
	{
		for(size_t i = 0 ; i < getExtra2dMeshes()->size() ; i++)
			copy->addMesh((*getExtra2dMeshes())[i]);
	}
	if(getExtra3dMeshes())
	{
		for(size_t i = 0 ; i < getExtra3dMeshes()->size() ; i++)
			copy->addMesh((*getExtra3dMeshes())[i]);
	}
	return copy ; 
//	return new HomogeneisedBehaviour( *this ) ;
}

std::vector<BoundaryCondition * > HomogeneisedBehaviour::getBoundaryConditions( const ElementState &s,  size_t id, const Function &p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv ) const
{
	return equivalent->getBoundaryConditions( s, id, p_i, gp, Jinv ) ;
}

Vector HomogeneisedBehaviour::getImposedStress( const Point &p, IntegrableEntity * e, int g ) const
{
	return equivalent->getImposedStress( p , e, g) ;
}

Vector HomogeneisedBehaviour::getImposedStrain( const Point &p, IntegrableEntity * e, int g ) const
{
	return equivalent->getImposedStrain( p , e, g) ;
}

Matrix HomogeneisedBehaviour::getTensor( const Point &p, IntegrableEntity * e, int g ) const
{
	return equivalent->getTensor( p, e,g  ) ;
}

