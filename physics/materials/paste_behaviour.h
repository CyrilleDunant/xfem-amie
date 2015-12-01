// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2013
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef PASTE_BEHAVIOUR_H
#define PASTE_BEHAVIOUR_H

#include "../weibull_distributed_stiffness.h"
#include "../logarithmic_creep_with_external_parameters.h"
#include "../../geometry/geometry_base.h"
#include "../../features/features.h"

namespace Amie
{

/*PARSE PasteBehaviour Form 
    @string<bool>[elastic] FALSE // desactivates the damage if TRUE
    @string<bool>[space_time] FALSE // uses the behaviour in space-time if TRUE
    @value[young_modulus] 12e9 // value of the Young modulus
    @value[poisson_ratio] 0.3 // value of the Poisson ratio
    @value[tensile_strength] 3e6 // value of the tensile strength of the material
    @value[short_term_creep_modulus] 3.6e9 // spring of the first KV chain 
    @value[long_term_creep_modulus] 4e9 // spring of the second KV chain 
    @string<SpaceDimensionality>[dimension] SPACE_TWO_DIMENSIONAL // number of dimensions of the current simulation
    @string<planeType>[plane_type] PLANE_STRESS // 2D hypothesis (plane strain or plane stress)
    @value[material_characteristic_radius] 0.000175 // radius of the non-local damage model
    @value[variability] 0.2 // variability of the mechanical properties
    @value[blocks] 0 // additional ghost blocks for viscoelastic simulations 
 */
struct PasteBehaviour : public WeibullDistributedStiffness
{
    bool spaceTime ;
    bool elastic ;
    int freeblocks ;
    double eta10 ;
    double eta300 ;

    PasteBehaviour(bool elastic = false, bool st = false, double E = 12e9, double nu = 0.3,  double up = 3e6, double shortTermCreepFactor = 3.6e9, double longTermCreepFactor = 4e9, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL, planeType pt = PLANE_STRESS, double r = 0.000175, double var = 0.2, int blocks = 0) ;

    virtual Form * getCopy() const ;

} ;

struct HydratingMechanicalCementPaste final: public LinearForm
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

    virtual bool hasInducedForces() const {
        return true ;
    }

    /** \brief Return the Vector of imposed Stress at the considered point. As the imposed stress is uniform, the point is ignored*/
    virtual Vector getImposedStress(const Point & p, IntegrableEntity * e, int g = -1) const ;
    virtual Vector getImposedStrain(const Point & p, IntegrableEntity * e, int g = -1) const ;

    std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;


} ;

struct HydratingDiffusionCementPaste final: public LinearForm
{
    std::vector<Variable> v ;
    Vector saturation ;
    double doh ;
    bool change ;

    HydratingDiffusionCementPaste() ;

    virtual Form * getCopy() const ;

    virtual bool isViscous() const {
        return true ;
    }

    virtual void step(double timestep, ElementState & currentState, double maxscore) ;

    virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;
    virtual void applyViscous(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine *vm) const ;

    virtual ~HydratingDiffusionCementPaste();

    double getDegreeOfHydration() const {
        return doh ;
    }

    /** \brief return true if the damage state has been modfied*/
    virtual bool changed() const ;

    /** \brief Return the (damaged) Stifness tensor*/
    virtual Matrix getTensor(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;

    virtual double getDeltadoh(double stauration, ElementState & currentState) ;
    virtual double getDiffusionCoefficient(double stauration, ElementState & currentState) ;

} ;



/*PARSE LogCreepPasteBehaviour Form 
    @string<bool>[elastic] FALSE // desactivates the damage if TRUE
    @value[young_modulus] 12e9 // value of the Young modulus
    @value[poisson_ratio] 0.3 // value of the Poisson ratio
    @value[tensile_strength] 3e6 // value of the tensile strength of the material
    @value[creep_modulus] 40e9 // long-term creep viscosity
    @value[creep_characteristic_time] 2 // creep rate in the logarithmic scale
    @string<SpaceDimensionality>[dimension] SPACE_TWO_DIMENSIONAL // number of dimensions of the current simulation
    @string<planeType>[plane_type] PLANE_STRESS // 2D hypothesis (plane strain or plane stress)
    @value[material_characteristic_radius] 0.000175 // radius of the non-local damage model
    @value[variability] 0.2 // variability of the mechanical properties
 */
struct LogCreepPasteBehaviour : public LogarithmicCreepWithExternalParameters
{
    LogCreepPasteBehaviour(bool elastic = false, double E = 12e9, double nu = 0.3, double up = 3e6, double eta = 40e9, double tau = 2, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL, planeType pt = PLANE_STRESS, double r = 0.000175, double var = 0.2) ;
} ;




}

#endif // PASTE_BEHAVIOUR
