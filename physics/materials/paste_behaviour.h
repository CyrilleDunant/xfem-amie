// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2013
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef PASTE_BEHAVIOUR_H
#define PASTE_BEHAVIOUR_H

#include "../weibull_distributed_stiffness.h"
#include "../../geometry/geometry_base.h"
#include "../../features/features.h"

namespace Amie
{
	struct PasteBehaviour : public WeibullDistributedStiffness
	{
	    // ultimate tensile strength: 20 MPa which is 5 * 4 MPa (instantaneous stress)
		double up ;
		double yield ;
		double c ;
		PasteBehaviour(double E = 12e9, double nu = 0.3,  double up = 0.0003, double yield = 0.0001666, double c = 12000., SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) ;
		
		virtual Form * getCopy() const ;
		
	} ;
	
	struct HydratingMechanicalCementPaste : public LinearForm
	{
		
		std::vector<Variable> v ;
		FeatureTree * diffusionTree ;
		double bulk ;
		double bulkSolid ;
		Vector imposed ;
		
		HydratingMechanicalCementPaste(FeatureTree * diffusionTree) ;
		
		virtual Form * getCopy() const ;
		
		virtual void step(double timestep, ElementState & currentState, double maxscore) ;
		
		virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;
		
		virtual bool fractured() const ;

		void makeBulkModuli(double saturation, double doh) ;
		double getCapillaryPressure(double saturation, double doh) ;

		virtual Matrix getMechanicalProperties(double saturation, double doh) ;
		virtual Vector getAutogeneousDeformation(double saturation, double doh) ;
		
		virtual ~HydratingMechanicalCementPaste();

		/** \brief return true if the damage state has been modfied*/
		virtual bool changed() const ;

		/** \brief Return the (damaged) Stifness tensor*/
		virtual Matrix getTensor(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;

		virtual bool hasInducedForces() const { return true ; }
		
		/** \brief Return the Vector of imposed Stress at the considered point. As the imposed stress is uniform, the point is ignored*/
		virtual Vector getImposedStress(const Point & p, IntegrableEntity * e, int g = -1) const ;
		virtual Vector getImposedStrain(const Point & p, IntegrableEntity * e, int g = -1) const ;

		std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;


	} ;
	
	struct HydratingDiffusionCementPaste : public LinearForm
	{
		std::vector<Variable> v ;
		Vector saturation ;
		double doh ;
		bool change ;
		
		HydratingDiffusionCementPaste() ;
		
		virtual Form * getCopy() const ;
		
		virtual bool isViscous() const { return true ; } 
		
		virtual void step(double timestep, ElementState & currentState, double maxscore) ;
		
		virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;
		virtual void applyViscous(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine *vm) const ;
		
		virtual ~HydratingDiffusionCementPaste();

		double getDegreeOfHydration() const { return doh ; }
		
		/** \brief return true if the damage state has been modfied*/
		virtual bool changed() const ;
		
		/** \brief Return the (damaged) Stifness tensor*/
		virtual Matrix getTensor(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;

		virtual double getDeltadoh(double stauration, ElementState & currentState) ;
		virtual double getDiffusionCoefficient(double stauration, ElementState & currentState) ;
		
	} ;
	
	struct ElasticOnlyPasteBehaviour : public PasteBehaviour
	{
		ElasticOnlyPasteBehaviour(double E = 12e9, double nu = 0.3, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) ;
		
		virtual Form * getCopy() const ;	
	} ;
	
	// parameters got for ./alain 1 500 500 0.29 10 200 at revision 2512
	struct ViscoElasticOnlyPasteBehaviour : public PasteBehaviour
	{
		double e_1 ;
		double e_2 ;
		int freeblocks ;
		
		ViscoElasticOnlyPasteBehaviour(double E = 12e9, double nu = 0.3, double e1 = 0.3, double e2 = 0.37, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) ;
		
		virtual Form * getCopy() const ;
	} ;

	struct ShortTermViscoElasticOnlyPasteBehaviour : public ViscoElasticOnlyPasteBehaviour
	{
		ShortTermViscoElasticOnlyPasteBehaviour(double E = 12e9, double nu = 0.3, double e1 = 0.3, double e2 = 0.3, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) : ViscoElasticOnlyPasteBehaviour(E, nu, e1,e2,dim) { } ;
		
		virtual Form * getCopy() const ;
	} ;

	typedef enum
	{
		STRAIN_CRITERION,
		STRESS_CRITERION,
		MIXED_CRITERION,
	} PasteCriterion ;

	struct ViscoDamagePasteBehaviour : public PasteBehaviour
	{
		double e_1 ;
		double e_2 ;
		int freeblocks ;
		PasteCriterion ctype ;
		double stressFraction ;
		
		ViscoDamagePasteBehaviour(double E=12e9, double nu = 0.3, double e1=0.3, double e2=0.37, double up = 0.0003, double r  = 0.00018, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) ;
		
		virtual Form * getCopy() const ;
	} ;

	struct ShortTermViscoDamagePasteBehaviour : public ViscoDamagePasteBehaviour
	{
		ShortTermViscoDamagePasteBehaviour(double E = 12e9, double nu = 0.3, double e1 = 0.3, double e2 = 0.3,double r  = 0.00018,  SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) : ViscoDamagePasteBehaviour(E, nu, e1,e2,r,dim) { } ;
		
		virtual Form * getCopy() const ;
	} ;

	struct PseudoBurgerViscoElasticOnlyPasteBehaviour : public PasteBehaviour
	{
		double e_1 ;
		double t_2 ;
		
		PseudoBurgerViscoElasticOnlyPasteBehaviour(double E=12e9, double nu = 0.3, double e1=0.3, double t2=300, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) ;
		
		virtual Form * getCopy() const ;
	} ;

	struct PseudoBurgerViscoDamagePasteBehaviour : public PasteBehaviour
	{
		double e_1 ;
		double t_2 ;
		int freeblocks ;
		PasteCriterion ctype ;
		double stressFraction ;
		
		PseudoBurgerViscoDamagePasteBehaviour(double E=12e9, double nu = 0.3, double e1=0.3, double t2=300, double up = 0.0003, double r  = 0.00018, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) ;
		
		virtual Form * getCopy() const ;
	} ;

} ;

#endif // PASTE_BEHAVIOUR
