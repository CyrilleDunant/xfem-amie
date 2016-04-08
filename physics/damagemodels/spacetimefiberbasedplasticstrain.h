//
// C++ Interface: lineardamage
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MU_SPACETIME_FIBERPLASTICITY_H
#define MU_SPACETIME_FIBERPLASTICITY_H

#include "damagemodel.h"
#include "spacetimefiberbasedplasticstrain.h"

namespace Amie {

/*PARSE SpaceTimeFiberBasedPlasticStrain DamageModel --no-suffix
    @value[damage_increment] 1e-6 // damage increment which is applied at each step of the damage algorithm
    @value[time_tolerance] 1e-5 // time lapsed between two damage events
    @value[maximum_damage] 0.999 // damage above which an element is considered broken
    @string<bool>[fixed_orientation] false // fix the plastic direction once it has been set
*/
class SpaceTimeFiberBasedPlasticStrain : public DamageModel
{
protected:
    double fibreFraction ;
    double timeTolerance ;
    Matrix C ;
    Vector plasticStrain ;
    double orientation = 0. ;
    bool fixed ;
    std::map<std::string, double> values ;

public:
    /** \brief Constructor. Set the number of degrees of freedom
     *
     * @param numDof number of degrees of freedom
     */
    SpaceTimeFiberBasedPlasticStrain(double f = 1e-6, double tol = 0.001, double yield = 1., bool fix = false) ;

    virtual ~SpaceTimeFiberBasedPlasticStrain();

    /** \brief Increment the damage
     *
     * @param s ElementState passed as a parameter
     */
    virtual std::pair<Vector, Vector> computeDamageIncrement(ElementState & s) /*override*/;
    virtual void computeDelta(ElementState & s) ;

    /** \brief Increment the damage from the current state of the element considered
     *
     * @param s ElementState
     */
    virtual void step(ElementState & s, double maxscore)  ;

    /** \brief compute the new stifness matrix after damage
     *
     * \f$ K' = K(1-d) \f$
     * @param m Matrix to modify
     * @return the new Matrix
     */
    virtual Matrix apply(const Matrix & m, const Point & p = Point(), const IntegrableEntity * e = nullptr, int g = -1) const ;
    virtual Matrix applyViscous(const Matrix & m, const Point & p = Point(), const IntegrableEntity * e = nullptr, int g = -1) const ;

    virtual bool fractured(int direction = -1) const  ;

    virtual DamageModel * getCopy() const ;

    virtual bool hasInducedBoundaryConditions() const {
        return true ;
    } ;
    virtual bool hasInducedForces() const {
        return true ;
    }

    virtual std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s, size_t id,  const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;
    virtual Vector getImposedStress(const Point & p) const ;
    virtual Vector getImposedStrain(const Point & p) const ;


};

/*PARSE SpaceTimeFiberBasedPlasticDamage DamageModel --no-suffix
    @string[damage_function] // gets the damage variable as a function of the plastic strain
    @string[requirements] // defines the extra material variables needed to calculate the damage
    @value[damage_increment] 1e-6 // damage increment which is applied at each step of the damage algorithm
    @value[time_tolerance] 1e-5 // time lapsed between two damage events
    @value[maximum_damage] 0.999 // damage above which an element is considered broken
    @string<bool>[fixed_orientation] false // fix the plastic direction once it has been set
*/
class SpaceTimeFiberBasedPlasticDamage : public SpaceTimeFiberBasedPlasticStrain
{
protected:
    Function damage ;
    std::string yCoord ;
    std::string zCoord ;
    std::string tCoord ;
    std::string uCoord ;
    std::string vCoord ;
    std::string wCoord ;
    bool needExternalVariables ;

public:
    SpaceTimeFiberBasedPlasticDamage( std::string damage, std::string req = "",  double f = 1e-6, double tol = 1e-4, double yield = 1., bool fix = false) ;
    SpaceTimeFiberBasedPlasticDamage( Function d, std::string y ="", std::string z ="", std::string t ="", std::string u ="", std::string v ="", std::string w ="", double f = 1e-6, double tol = 1e-4, double yield = 1., bool fix = false ) ;

    virtual Matrix apply(const Matrix & m, const Point & p = Point(), const IntegrableEntity * e = nullptr, int g = -1) const ;
    virtual Matrix applyViscous(const Matrix & m, const Point & p = Point(), const IntegrableEntity * e = nullptr, int g = -1) const ;

    virtual Vector getImposedStress(const Point & p) const ;
    
    virtual void step(ElementState & s, double maxscore)  ;

    virtual DamageModel * getCopy() const ;
} ;



}

#endif
