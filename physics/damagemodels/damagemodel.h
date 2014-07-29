//
// C++ Interface: damagemodel
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __DAMAGE_MODEL_H__
#define __DAMAGE_MODEL_H__

#include "../../elements/integrable_entity.h"
#include "../../utilities/matrixops.h"

namespace Amie
{

	typedef enum{
		DISSIPATIVE,
		CONSERVATIVE
	} ConvergenceType ;

	struct PointState
	{
		bool isMet ;
		double delta ;
		double score ;
		double fraction ;
		double proximity ;
		double angleShift ;
		int mode ;
		PointState(bool met, double delta, double frac, double score, double proximity, double angleShift, double mode) : isMet(met), delta(delta), score(score), fraction(frac), proximity(proximity), angleShift(angleShift), mode(mode) {} ;
		bool operator < (const PointState & p) const {return fraction < p.fraction ; }
		void print() const
		{
			std::cout << fraction << " " << isMet << " delta = " << delta << ", score = " << score << ", proximity = " << proximity << ", angleShift = "<< angleShift<< ", mode = "<< mode<< std::endl; 
		}
	} ;
	struct RangeState
	{
		PointState up ;
		PointState down ;
		RangeState(const PointState &up, const PointState & down) : up(up), down(down) {} ;
		RangeState split(const PointState & newMid)
		{
			RangeState ret(newMid, up) ;
			up = newMid ;
			return ret ;
		}
		PointState extrapolate(double ratio = .5)
		{
			return PointState(up.isMet && down.isMet, up.delta*ratio + down.delta*(1.-ratio), up.fraction*ratio+down.fraction*(1.-ratio), up.score*ratio+down.score*(1.-ratio), up.proximity*ratio+down.proximity*(1.-ratio),up.angleShift*ratio+down.angleShift*(1.-ratio),round(up.mode*ratio+down.mode*(1.-ratio))) ;
		}

		double zeroLocation() const
		{
			double ret = down.fraction + down.delta*(up.fraction-down.fraction)/(up.delta-down.delta) ;
			if(ret >= down.fraction && ret <= up.fraction)
				return ret ;
			return -1 ;
		}

		double equilibriumLocation() const
		{
			double ret = down.fraction + down.score*(up.fraction-down.fraction)/(up.score-down.score) ;
			if(ret >= down.fraction && ret <= up.fraction)
				return ret ;
			return -1 ;

		}
		
	} ;
/** \brief Damage model interface */
class DamageModel
{
protected:
	double thresholdDamageDensity ;
	double secondaryThresholdDamageDensity ;
	double damageDensityTolerance ;
	double fraction ;
	double trialRatio ;
	
	Vector upState ;
	Vector downState ;
	double delta ;
	double effectiveDeltaFraction ;
	
	bool change ;
	bool haslimit ;
	
	Vector state ;
	Vector initalState ;

	std::vector<PointState> states ;
	
	ElementState * elementState ;
	
	ConvergenceType ctype ;
	bool alternating ;
	bool alternate ;
	bool needGlobalMaximumScore ;

	double iterationNumber ;
public:
	
	bool isNull ;
	bool converged ;
	double error ;

	DamageModel();
	
	double getThresholdDamageDensity() const;
	
	double getSecondaryThresholdDamageDensity() const;
	
	double getDamageDensityTolerance() const ;

	bool getNeedGlobalMaximumScore() const { return needGlobalMaximumScore ; }

	void setThresholdDamageDensity(double d);
	
	void setSecondaryThresholdDamageDensity(double d) ;
		
	void setDamageDensityTolerance(double d) ;
	
	void setConvergenceType(ConvergenceType ct) {ctype = ct ;}

	void setNeedGlobalMaximumScore(bool m) { needGlobalMaximumScore = m ; }
	
	double getDamageDensityTolerance() { return damageDensityTolerance ; };
	bool hasConverged() const {return converged ; }
	virtual double getDelta() const {return delta ;}
	virtual int getMode() const 
	{
		if(fractured())
			return 1 ;
		return -1 ;
	}
	virtual double getAngleShift() const { return 0. ;}
	virtual void computeDelta(const ElementState &s) = 0 ;
	
	/** \brief Return a vector of values describing the damage stage of the material
	 * 
	 * @return a Vector
	 */
	virtual const Vector & getState() const { return state ;} ;
	virtual Vector & getState(bool) ;
	
	/** \brief Increment the damage from the current state of the element considered
	 * 
	 * @param s ElementState
	 */
	virtual void step(ElementState & s, double maxscore)  ;

	virtual std::pair<Vector,Vector> computeDamageIncrement(ElementState &s) = 0 ;

	/** \brief return true is the element concerned is fractured 
		*/
	virtual bool fractured() const = 0 ;

	/** \brief Modify the rigidity matrix according to the damage model
	 * 
	 * @param m Matrix to modify
	 * @return a new Matrix
	 */
	virtual Matrix apply(const Matrix & m, const Point & p = Point(), const IntegrableEntity * e = nullptr, int g = -1) const = 0 ;
	virtual Matrix applyViscous(const Matrix & m, const Point & p = Point(), const IntegrableEntity * e = nullptr, int g = -1) { return m ; }
	
	virtual DamageModel * getCopy() const = 0 ;
	/** \brief Modify the rigidity matrix according to the damage model function (for space-time behaviour)
	 * 
	 * @param m Matrix to modify
	 * @return a new FunctionMatrix
	 */
	
	virtual bool changed() const ;

	virtual void prepare() { };
	
	virtual void postProcess();
	
	virtual bool hasInducedBoundaryConditions() const
	{
		return false ;
	} ;
	virtual bool hasInducedForces() const {return false ;}
	
	virtual std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s, size_t id,  const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
	{
		return std::vector<BoundaryCondition * >() ;
	}
	
	virtual Vector getImposedStress(const Point & p) const
	{
		return Vector(0., 0) ;
	}
	virtual Vector getImposedStrain(const Point & p) const
	{
		return Vector(0., 0) ;
	}

} ;

/** \brief Null damage model. A material with this damage model can never be damaged.*/
class NullDamage : public DamageModel
{

public:

	NullDamage() : DamageModel() { state.resize(0);} ;

	/** \brief Return the vector of variables describing the damage state
	 * 
	 * @return damage state vector
	 */

	/** \brief Do nothing
	 * 
	 * @param s ElementState
	 */
	virtual std::pair<Vector, Vector> computeDamageIncrement(ElementState & s)  { return std::make_pair(Vector(1), Vector(1)) ;} /*override*/;

	/** \brief Do nothing
	 * 
	 * @param d damage
	 */
	virtual void artificialDamageStep(double d) { } ;

	/** \brief returns 0 */
	virtual Vector getPreviousDamage() {return Vector(0) ; } ;

	/** \brief returns 0 */
	virtual Vector getPreviousPreviousDamage() {return Vector(0) ; } ;

	/** \brief Do nothing */
	virtual void artificialPreviousDamage(Vector previous, Vector previousprevious) { } ;

	/** \brief return a copy of the Matrix given as argument
	 * 
	 * @param m Matrix
	 * @return m
	 */
	virtual Matrix apply(const Matrix & m, const Point & p = Point(), const IntegrableEntity * e = nullptr, int g = -1) const {return m ;} ;
	
	virtual bool fractured() const {return false ; } ;
	
	virtual DamageModel * getCopy() const { return new NullDamage() ;}
	
	virtual void computeDelta(const ElementState & s) {} ;

// 	virtual FunctionMatrix apply(const Matrix & m) const
// 	{
// 		FunctionMatrix ret(m.numRows(), m.numCols()) ;
// 		ret += m ;
// 		return ret ;
// 	}
} ;

}

#endif
