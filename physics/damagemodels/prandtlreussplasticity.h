//
// C++ Interface: isotropiclineardamage
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MU_PRANDTL_REUSS_PLASTIC_STRAIN_H
#define MU_PRANDTL_REUSS_PLASTIC_STRAIN_H

#include "damagemodel.h"

namespace Amie {

/** \brief Isotropic linear damage model. The stifness of an affected element is scaled by a factor between 1 and 1 - .9999999
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/


/*PARSE PrandtlReussPlasticStrain DamageModel --no-suffix
*/
class PrandtlReussPlasticStrain final: public DamageModel
{
public:
    std::vector<Variable> v ;
    Matrix * param ;
    Vector imposedStrain ;
    Vector previousImposedStrain ;
    Vector imposedStrainAccumulator ;
    double plasticVariable ;

    bool broken ;
    bool inCompression ;
    bool inTension ;
    bool newtonIteration ;

    ElementState * es ;

public:
    double factor ;

    /** \brief Constructor. Set the number of degrees of freedom
     *
     * @param numDof number of degrees of freedom
     */
    PrandtlReussPlasticStrain() ;

    virtual ~PrandtlReussPlasticStrain();

    /** \brief Increment the damage
     *
     * @param s ElementState passed as a parameter
     */
    virtual std::pair<Vector, Vector> computeDamageIncrement(ElementState & s) /*override*/;
    virtual void computeDelta(ElementState & s) ;
    virtual double getAngleShift() const ;

    /** \brief compute the new stifness matrix after damage
     *
     * @param m Matrix to modify
     * @return m
     */
    virtual Matrix apply(const Matrix & m, const Point & p = Point(), const IntegrableEntity * e = nullptr, int g = -1) const;

    virtual bool hasInducedBoundaryConditions() const {
      if(previousImposedStrain.size() ==0)
	return false ;
      return true ;
    } ;
    virtual bool hasInducedForces() const {
       if(previousImposedStrain.size() ==0)
	 return false ;
       return true ;
    }


    /** \brief return true is the element concerned is fractured
    	*/
    virtual bool fractured(int direction = -1) const  ;

// 	virtual void step(ElementState & s) ;

    virtual std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s, size_t id,  const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;
    virtual Vector getImposedStress(const Point & p) const ;
    virtual Vector getImposedStrain(const Point & p) const ;
    virtual int getMode() const ;

    virtual DamageModel * getCopy() const ;

    virtual void preProcess( double timeStep, ElementState & currentState ) ;
    virtual void postProcess() ;
    double getDamage() const ;
    double getPlasticity() const;

    virtual void step( ElementState &s , double maxscore) ;
};

}

#endif
