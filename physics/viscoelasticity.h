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

#ifndef __VISCOELASTICTY_H_
#define __VISCOELASTICTY_H_

#include "physics_base.h"

namespace Mu
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

void placeMatrixInBlock(const Matrix & m, size_t i, size_t j, Matrix & ret) ;
void addMatrixInBlock(const Matrix & m, size_t i, size_t j, Matrix & ret) ;
void substractMatrixInBlock(const Matrix & m, size_t i, size_t j, Matrix & ret) ;
void getBlockInMatrix( const Matrix & source, size_t i, size_t j, Matrix & ret) ;


struct Viscoelasticity : public LinearForm
{
	// generalized viscosity tensor
	Matrix eta ;
	
	// density
	double rho ;
	
	// viscoelastic model represented
	ViscoelasticModel model ;
	
	// number of blocks (= 1+number of additionnal dofs)
	int blocks ;
	int effblocks ;
	
	std::vector<Variable> v ;
	
	// constructor for pure elasticity or pure viscosity
	Viscoelasticity( ViscoelasticModel model, const Matrix & rig, int additionnalBlocksAfter = 0, double r = 0) ; 
	// constructor for elementary Kelvin-Voigt or Maxwell
	Viscoelasticity( ViscoelasticModel model, const Matrix & rig, const Matrix & eta, int additionnalBlocksBefore = 0, int additionnalBlocksAfter = 0, double r = 0) ; 
	// constructor for Burger (KV + Maxwell in serial)
	Viscoelasticity( ViscoelasticModel model, const Matrix & c_kv, const Matrix & e_kv, const Matrix & c_mx, const Matrix & e_mx, int additionnalBlocksBefore = 0, int additionnalBlocksAfter = 0, double r = 0) ; 
	// constructor for generalized KelvinVoigt or Maxwell
	Viscoelasticity( ViscoelasticModel model, const Matrix & c_0, std::vector<std::pair<Matrix, Matrix> > & branches, int additionnalBlocksBefore = 0, int additionnalBlocksAfter = 0, double r = 0) ;
	// constructor for generalized KelvinVoigt or Maxwell with 1 module only
	Viscoelasticity( ViscoelasticModel model, const Matrix & c_0, const Matrix & c_1, const Matrix & e_1, int additionnalBlocksBefore = 0, int additionnalBlocksAfter = 0, double r = 0) ;
	// constructor for general viscoelasticity (rig and eta are supposed symmetric)
	Viscoelasticity( const Matrix & rig, const Matrix & eta, int blocks, int additionnalBlocksAfter = 0, double r = 0) ; 

	virtual ~Viscoelasticity() ;

	virtual void apply( const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;
	virtual void applyViscous( const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm ) const ;

	virtual bool fractured() const ;

	virtual bool isViscous() const { return true ; }

	virtual ElementState * createElementState( IntegrableEntity * e) ;

	virtual Form * getCopy() const ;

	virtual bool changed() const ;

	virtual Vector getImposedStress(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;
	virtual Vector getImposedStrain(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;
	
	virtual Matrix getTensor(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const { return param ; }	
	
	virtual Matrix getViscousTensor(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const { return eta ; }	
	
	virtual Vector getForcesFromAppliedStress( Vector & data, Function & shape, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<Variable> & v, bool isVolumic = false) ;

	virtual Vector getForcesFromAppliedStress( const Function & data, size_t index, size_t externaldofs,  Function & shape, IntegrableEntity * e,const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<Variable> & v, bool isVolumic = false) ;
	
	virtual void print() const ;
	
  
} ;

} ;


#endif
