// C++ Interface: generalized_spacetime_viscoelasticity
//
// Description: Generalized visco-elastic behaviour for the Space-Time Finite Element Method
//
//
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __GENERALIZED_SPACETIME_VISCOELASTICTY_H_
#define __GENERALIZED_SPACETIME_VISCOELASTICTY_H_

#include "physics_base.h"

namespace Amie
{

typedef enum
{
	BURGER,
	KELVIN_VOIGT,
	MAXWELL,
	GENERALIZED_KELVIN_VOIGT,
	GENERALIZED_MAXWELL,
	PURE_ELASTICITY,
	PURE_VISCOSITY,
	GENERAL_VISCOELASTICITY,
} ViscoelasticModel ;

struct GeneralizedSpaceTimeViscoelasticity : public LinearForm
{
	// generalized viscosity tensor
	Matrix eta ;
	
	// density
	double rho ;
	
	// viscoelastic model represented
	ViscoelasticModel model ;
	
	// number of blocks (= 1+number of additionnal dofs)
	int blocks ;
	
	std::vector<Variable> v ;
	
	// constructor for pure elasticity or pure viscosity
	GeneralizedSpaceTimeViscoelasticity( ViscoelasticModel model, const Matrix & rig, int additionnalBlocksAfter = 0, double r = 0) ; 
	// constructor for elementary Kelvin-Voigt or Maxwell
	GeneralizedSpaceTimeViscoelasticity( ViscoelasticModel model, const Matrix & rig, const Matrix & eta, int additionnalBlocksBefore = 0, int additionnalBlocksAfter = 0, double r = 0) ; 
	// constructor for Burger (KV + Maxwell in serial)
	GeneralizedSpaceTimeViscoelasticity( ViscoelasticModel model, const Matrix & c_kv, const Matrix & e_kv, const Matrix & c_mx, const Matrix & e_mx, int additionnalBlocksBefore = 0, int additionnalBlocksAfter = 0, double r = 0) ; 
	// constructor for generalized KelvinVoigt or Maxwell
	GeneralizedSpaceTimeViscoelasticity( ViscoelasticModel model, const Matrix & c_0, std::vector<std::pair<Matrix, Matrix> > & branches, int additionnalBlocksBefore = 0, int additionnalBlocksAfter = 0, double r = 0) ;
	// constructor for generalized KelvinVoigt or Maxwell with 1 module only
	GeneralizedSpaceTimeViscoelasticity( ViscoelasticModel model, const Matrix & c_0, const Matrix & c_1, const Matrix & e_1, int additionnalBlocksBefore = 0, int additionnalBlocksAfter = 0, double r = 0) ;
	// constructor for general viscoelasticity (rig and eta are supposed symmetric)
	GeneralizedSpaceTimeViscoelasticity( const Matrix & rig, const Matrix & eta, int blocks, int additionnalBlocksAfter = 0, double r = 0) ; 

	virtual ~GeneralizedSpaceTimeViscoelasticity() ;

	virtual void apply( const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm ) const ;
	virtual void applyViscous( const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm ) const ;

	virtual bool fractured() const ;

	virtual bool isViscous() const { return true ; }

	virtual ElementState * createElementState( IntegrableEntity * e) ;

	virtual Form * getCopy() const ;

	virtual bool changed() const ;

	virtual void scale( double d )
	{
		param *= d ;
		eta *= d ;
	}
	
	virtual virtual Vector getForcesFromAppliedStress( const Vector & data, Function & shape, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<Variable> & v, bool isVolumic = false, const Vector & normal = Vector()) ;

	virtual Vector getForcesFromAppliedStress( const Function & data, size_t index, size_t externaldofs,  Function & shape, IntegrableEntity * e,const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<Variable> & v, bool isVolumic = false, const Vector & normal = Vector()) ;
	
  
} ;

struct MassAndViscosityAndStiffnessByBlock : public LinearForm
{
	std::vector<Variable> v ;
	Matrix viscosity ;
	double density ;
	int blocks ;
	
	MassAndViscosityAndStiffnessByBlock(const Matrix & rig, const Matrix & eta, double rho, int b = 1) ;
	
	virtual ~MassAndViscosityAndStiffnessByBlock() ;

	virtual XMLTree * toXML() { return new XMLTree("stiffness",param) ; } ;
	
	virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;
	
	virtual bool fractured() const { return false ; } 
	
	virtual Form * getCopy() const { return new MassAndViscosityAndStiffnessByBlock(*this) ; }

	virtual ElementState * createElementState( IntegrableEntity * e) ;
	
} ;

} ;


#endif
