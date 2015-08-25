// C++ Interface: logarithmic_creep
//
// Description: Logarithmic visco-elastic behaviour for the Space-Time Finite Element Method
//
//
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2007-2014
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __LOGARITHMIC_CREEP_WITH_EXTERNAL_PARAMETERS_H_
#define __LOGARITHMIC_CREEP_WITH_EXTERNAL_PARAMETERS_H_

#include "logarithmic_creep.h"
#include "../elements/generalized_spacetime_viscoelastic_element_state.h"
#include "material_laws/material_laws.h"

namespace Amie
{

/* This is a generic mechanical behaviour. It is composed by a spring and a non-linear (logarithmic) dashpot placed in parallels with an imposed deformation.
 *
 * The behaviour relies on two components: a list of material parameters (mapped with their default values) and a list of material laws. The parameters and the laws must be set prior to any calculation
 *
 * The list of material parameters must include:
 *      - "young_modulus" and "poisson_ratio" for the material to be valid (or alternatively "bulk_modulus" and "shear_modulus")
 *      - "creep_modulus", "creep_poisson" and "creep_characteristic_time" for the creep law to activate (or "creep_bulk", "creep_shear" and "creep_characteristic_time") ; without the material is considered as purely elastic
 *      - "imposed_deformation" for the imposed strain to activate ; otherwise the material has an imposed strain of 0.
 *
 * During the pre- and post-processing, the material laws are applied. They are called in order of appearance in the vector of material laws.
 *
 * The local values of the parameters are stored in the element state (explicitely described as a GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables ) and not in the behaviour. The behaviour contains the default values for the parameters in case they are not found in the element state.
 */

/*PARSE LogarithmicCreepWithExternalParameters Form
    @value[] // map of strings and doubles
    @object<ExternalMaterialLawList>[relations] nullptr // vector of material laws to apply
    @object<FractureCriterion>[fracture_criterion] nullptr // fracture criterion to calculate when damage occurs
    @object<DamageModel>[damage_model] nullptr // damage algorithm
    @object<LogCreepAccumulator>[accumulator] RealTimeLogCreepAccumulator() // how the logrithmic creep is obtained
    @string<SpaceDimensionality>[dimension] SPACE_TWO_DIMENSIONAL // number of dimensions of the current simulation
    @string<planeType>[plane_type] PLANE_STRESS // 2D hypothesis (plane strain or plane stress)
*/ 
struct LogarithmicCreepWithExternalParameters : public LogarithmicCreepWithImposedDeformationAndFracture
{
public:
    planeType plane ;
protected:
    std::map<std::string, double> external ;
    std::vector< ExternalMaterialLaw * > relations ;

    virtual void makeProperties(std::map<std::string, double> & values, double kvSpringReduction = 1., double kvDashpotReduction = 1.) ;

public:
    // the argument string contains the list of external variables and their default values (see the function parseDefaultValues() in external_material_laws.h for more details)
    LogarithmicCreepWithExternalParameters(std::string args, LogCreepAccumulator * acc = new RealTimeLogCreepAccumulator(), SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL, planeType pt = PLANE_STRESS, char sep = ',') ;
    LogarithmicCreepWithExternalParameters(std::string args, std::string ptension, std::string pcompression , DamageModel * d, LogCreepAccumulator * acc = new RealTimeLogCreepAccumulator(), SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL, planeType pt = PLANE_STRESS, char sep = ',') ;
    LogarithmicCreepWithExternalParameters(std::string args, FractureCriterion * c, DamageModel * d, LogCreepAccumulator * acc = new RealTimeLogCreepAccumulator(), SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL, planeType pt = PLANE_STRESS, char sep = ',') ;
    LogarithmicCreepWithExternalParameters(std::map<std::string, double> values, ExternalMaterialLawList * mats, FractureCriterion * c, DamageModel * d, LogCreepAccumulator * acc = new RealTimeLogCreepAccumulator(), SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL, planeType pt = PLANE_STRESS) ;
    virtual ~LogarithmicCreepWithExternalParameters() ; //{ if(dfunc) {delete dfunc ; } if(criterion) {delete criterion ;} }

    virtual ElementState * createElementState( IntegrableEntity * e) ;

    virtual Form * getCopy() const ;

    virtual void step(double timestep, ElementState &s, double maxScore) ;

    virtual void preProcess( double timeStep, ElementState & currentState ) ;

    virtual bool hasInducedForces() const {
        return !fractured() && (fixCreepVariable || imposed.size()>0)  ;
    }

    void setMaterialParameters( std::map<std::string, double> & e) { external = e ; }
    void addMaterialParameter( std::string name, double defaultValue = 0) { external[name] = defaultValue ; }

    void addMaterialLaw( ExternalMaterialLaw * law) { if(law){ relations.push_back(law) ; } }
    void addMaterialLaws( std::vector<ExternalMaterialLaw *> laws) { for(size_t i = 0 ; i < laws.size() ; i++) { if(laws[i]) { relations.push_back(laws[i]);} } }
    void addMaterialLaws( ExternalMaterialLaw* laws[], size_t length) { for(size_t i = 0 ; i < length ; i++) { if(laws[i]) { relations.push_back(laws[i]);} } }
    void replaceMaterialLaw( ExternalMaterialLaw * law, int index) ;
    void removeMaterialLaw( int index ) ;
    void print() { std::cout << relations.size() << std::endl ; }

    virtual std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;


} ;


} 


#endif
