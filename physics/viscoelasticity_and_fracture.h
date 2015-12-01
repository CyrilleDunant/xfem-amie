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

#ifndef __VISCOELASTICTY_AND_FRACTURE_H_
#define __VISCOELASTICTY_AND_FRACTURE_H_

#include "physics_base.h"
#include "viscoelasticity.h"
#include "fracturecriteria/fracturecriterion.h"
#include "damagemodels/fiberbasedisotropiclineardamage.h"

namespace Amie
{

/*PARSE ViscoelasticityAndFracture Form 
    @string<ViscoelasticModel>[model] // rheological assembly
    @value[young_modulus] // value of the Young modulus
    @value[poisson_ratio] // value of the Poisson ratio
    @object<FractureCriterion>[fracture_criterion] nullptr // fracture criterion to calculate when damage occurs
    @object<DamageModel>[damage_model] nullptr // damage algorithm
    @value[creep_characteristic_time] 1 // value of the viscosity characteristic time
    @string[file_name] // name of the file containing the modulus of the different branches
    @string<SpaceDimensionality>[dimension] SPACE_TWO_DIMENSIONAL // number of dimensions of the current simulation
    @string<planeType>[plane_type] PLANE_STRESS // 2D hypothesis (plane strain or plane stress)
 */
struct ViscoelasticityAndFracture : public Viscoelasticity
{	
	DamageModel * dfunc ;

	Matrix elasticParam ;
	Matrix viscousParam ;

	FractureCriterion * criterion ;
	
	// constructor for pure elasticity or pure viscosity
	ViscoelasticityAndFracture( ViscoelasticModel model, const Matrix & rig, FractureCriterion * c, DamageModel * d, int additionnalBlocksAfter = 0, double r = 0) ; 
	// constructor for elementary Kelvin-Voigt or Maxwell
	ViscoelasticityAndFracture( ViscoelasticModel model, const Matrix & rig, const Matrix & eta, FractureCriterion * c, DamageModel * d, int additionnalBlocksBefore = 0, int additionnalBlocksAfter = 0, double r = 0) ; 
	// constructor for Burger (KV + Maxwell in serial)
	ViscoelasticityAndFracture( ViscoelasticModel model, const Matrix & c_kv, const Matrix & e_kv, const Matrix & c_mx, const Matrix & e_mx, FractureCriterion * c, DamageModel * d, int additionnalBlocksBefore = 0, int additionnalBlocksAfter = 0, double r = 0) ; 
	// constructor for generalized KelvinVoigt or Maxwell
	ViscoelasticityAndFracture( ViscoelasticModel model, const Matrix & c_0, std::vector<std::pair<Matrix, Matrix> > & branches, FractureCriterion * c, DamageModel * d, int additionnalBlocksBefore = 0, int additionnalBlocksAfter = 0, double r = 0) ;
	// constructor for generalized KelvinVoigt or Maxwell with 1 module only
	ViscoelasticityAndFracture( ViscoelasticModel model, const Matrix & c_0, const Matrix & c_1, const Matrix & e_1, FractureCriterion * c, DamageModel * d, int additionnalBlocksBefore = 0, int additionnalBlocksAfter = 0, double r = 0) ;
	// constructor for general viscoelasticity (rig and eta are supposed symmetric)
	ViscoelasticityAndFracture( const Matrix & rig, const Matrix & eta, int blocks, FractureCriterion * c, DamageModel * d, int additionnalBlocksAfter = 0, double r = 0) ; 

	ViscoelasticityAndFracture( ViscoelasticModel model, double young, double poisson, FractureCriterion * c, DamageModel * d, double tau = 1, std::string file = std::string(), SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL, planeType pt = PLANE_STRESS, bool hooke = true, int additionnalBlocksBefore = 0, int additionnalBlocksAfter = 0) ;

	virtual ~ViscoelasticityAndFracture() ;

	virtual bool fractured() const ;

	virtual bool isViscous() const { return true ; }
	
	virtual void step(double timestep, ElementState & currentState, double maxscore) ;

	virtual FractureCriterion * getFractureCriterion() const ;

	virtual DamageModel * getDamageModel() const ;

	virtual void setFractureCriterion(FractureCriterion * frac) ;

	virtual ElementState * createElementState( IntegrableEntity * e) ;

	virtual Matrix getTensor(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;

	virtual Matrix getViscousTensor(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;

	virtual Form * getCopy() const ;

	virtual bool changed() const ;

	virtual void scale( double d )
	{
		param *= d ;
		eta *= d ;
	}
	
	virtual Vector getForcesFromAppliedStress( const Vector & data, Function & shape, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<Variable> & v, bool isVolumic = false, const Vector & normal = Vector()) ;

	virtual Vector getForcesFromAppliedStress( const Function & data, size_t index, size_t externaldofs,  Function & shape, IntegrableEntity * e,const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<Variable> & v, bool isVolumic = false, const Vector & normal = Vector()) ;
	
private:
	void setElasticAndViscousStiffnessMatrix() ;

  
} ;


} 


#endif
