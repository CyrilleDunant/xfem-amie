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
#ifndef MU_PLASTIC_STRAIN_H
#define MU_PLASTIC_STRAIN_H

#include "damagemodel.h"

namespace Amie {

/** \brief Isotropic linear damage model. The stifness of an affected element is scaled by a factor between 1 and 1 - .9999999
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/


/*PARSE PrandtlGrauertPlasticStrain DamageModel --no-suffix 
    @value[compressibility] 0.05 // 
    @value[limit_strain] 0.0057 //
    @value[initial_plastic_variable] 0 // 
*/
class PrandtlGrauertPlasticStrain final: public DamageModel
{
public:
    std::vector<Variable> v ;
    Matrix * param ;
    Vector imposedStrain ;
    Vector previousCompressiveImposedStrain ;
    Vector previousTensileImposedStrain ;
    Vector lastStress ;
    double c_psi ;
    double compressivePlasticVariable ;
    double tensilePlasticVariable ;
    double eps_f ;
    double kappa_0 ;
    double plasticFlowPotential(const Matrix & m) const ;

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
    PrandtlGrauertPlasticStrain(double c_psi = 0.05, double eps_f = 0.0057, double kappa_0 = 0) ;

    virtual ~PrandtlGrauertPlasticStrain();

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
        return true ;
    } ;
    virtual bool hasInducedForces() const {
        return true ;
    }


    /** \brief return true is the element concerned is fractured
    	*/
    virtual bool fractured(int direction = -1) const  ;
    
    virtual void step( ElementState &s , double maxscore) ;

// 	virtual void step(ElementState & s) ;

    virtual std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s, size_t id,  const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;
    virtual Vector getImposedStress(const Point & p) const ;
    virtual Vector getImposedStrain(const Point & p) const ;
    virtual int getMode() const ;

    virtual DamageModel * getCopy() const ;

    virtual void postProcess() ;
    double getDamage() const ;
    double getPlasticity() const;


};

}

#endif
