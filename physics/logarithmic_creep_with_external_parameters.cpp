#include "logarithmic_creep_with_external_parameters.h"
#include "../elements/generalized_spacetime_viscoelastic_element_state.h"

using namespace Amie ;

LogarithmicCreepWithExternalParameters::LogarithmicCreepWithExternalParameters(std::string args, FractureCriterion * c, DamageModel * d, SpaceDimensionality dim, char sep) : LogarithmicCreepWithImposedDeformationAndFracture( Matrix( 3+3*(dim == SPACE_THREE_DIMENSIONAL), 3+3*(dim == SPACE_THREE_DIMENSIONAL)), Vector(), c,d), external(parseDefaultValues(args, sep))
{
    noFracture = (c == nullptr) ;
}

LogarithmicCreepWithExternalParameters::LogarithmicCreepWithExternalParameters(std::string args, SpaceDimensionality dim, char sep) : LogarithmicCreepWithImposedDeformationAndFracture( Matrix( 3+3*(dim == SPACE_THREE_DIMENSIONAL), 3+3*(dim == SPACE_THREE_DIMENSIONAL)), Vector()), external(parseDefaultValues(args, sep))
{
    noFracture = true ;
}

ElementState * LogarithmicCreepWithExternalParameters::createElementState( IntegrableEntity * e)
{
    return new GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables(e, external) ;
}

Form * LogarithmicCreepWithExternalParameters::getCopy() const
{
    LogarithmicCreepWithExternalParameters * copy ;

    if(noFracture)
    {
        copy = new LogarithmicCreepWithExternalParameters( std::string(),  (param.numCols() == 6) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL)  ;
    }
    else
    {
        copy = new LogarithmicCreepWithExternalParameters( std::string(),criterion->getCopy(), dfunc->getCopy(),  (param.numCols() == 6) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL ) ;
        copy->dfunc->getState(true).resize(dfunc->getState().size());
        copy->dfunc->getState(true) = dfunc->getState() ;
        copy->criterion->setMaterialCharacteristicRadius(criterion->getMaterialCharacteristicRadius()) ;
        copy->dfunc->setDamageDensityTolerance(dfunc->getDamageDensityTolerance());
        copy->dfunc->setThresholdDamageDensity(dfunc->getThresholdDamageDensity());
    }

    copy->external = external ;
    for(size_t i = 0 ; i < relations.size() ; i++)
        copy->relations.push_back(relations[i]) ;

    copy->accumulator = accumulator ;

    return copy ;

}

void LogarithmicCreepWithExternalParameters::step(double timestep, ElementState &s, double maxScore)
{
    for(size_t i = 0 ; i < relations.size() ; i++)
        relations[i]->step( dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(s), timestep) ;

    LogarithmicCreepWithImposedDeformationAndFracture::step(timestep, s, maxScore) ;
}

void LogarithmicCreepWithExternalParameters::preProcess( double timeStep, ElementState & currentState )
{
    if(reducedTimeStep > POINT_TOLERANCE_2D)
        return ;

    accumulateStress(timeStep, currentState);

    GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables& state = dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(currentState) ;
    if(state.has("imposed_deformation"))
        state.set("imposed_deformation", 0.) ;


    for(size_t i = 0 ; i < relations.size() ; i++)
        relations[i]->preProcess( state, timeStep ) ;

    if(state.has("bulk_modulus"))
    {
        double k_inst = state.get("bulk_modulus", external) ;
        double mu_inst = state.get("shear_modulus", external) ;
        C = Material::cauchyGreen( k_inst, mu_inst, false, (param.numCols()==6) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL) ;
    }
    else
    {
        double E_inst = state.get("young_modulus", external) ;
        double nu_inst = state.get("poisson_ratio", external) ;
        C = Material::cauchyGreen( E_inst, nu_inst, true, (param.numCols()==6) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL) ;
    }
    
    placeMatrixInBlock( C, 0,0, param) ;
    if(state.has("creep_characteristic_time"))
    {
        isPurelyElastic = false ;
        if(state.has("creep_bulk"))
        {
            double k_visc = state.get("creep_bulk", external) ;
            double mu_visc = state.get("creep_shear", external) ;
            C = Material::cauchyGreen( k_visc, mu_visc, false, (param.numCols()==6) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL) ;
        }
        else
        {
            double E_visc = state.get("creep_modulus", external) ;
            double nu_visc = state.get("creep_poisson", external) ;
            E = Material::cauchyGreen( E_visc, nu_visc, true, (param.numCols()==6) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL) ;
        }
        tau = state.get("creep_characteristic_time", external) ;
        
        placeMatrixInBlock( C*(-1), 1,0, param) ;
        placeMatrixInBlock( C*(-1), 0,1, param) ;
        placeMatrixInBlock( C, 1,1, param) ;
        placeMatrixInBlock( E, 1,1, eta) ;

        model = MAXWELL ;
    }
    else
    {
        model = PURE_ELASTICITY ;
        isPurelyElastic = true ;
    }
    if(state.has("imposed_deformation"))
    {
        imposed.resize(C.numCols()) ;
        double a = state.get("imposed_deformation", external) ;
        for(size_t i = 0 ; i < 2+(imposed.size()==6) ; i++)
            imposed[i] = a ;
        if(std::abs(a) < POINT_TOLERANCE_2D)
            imposed.resize(0) ;

    }
    else
        imposed.resize(0) ;

    if(imposed.size() > 0 && prevImposed.size() == 0)
        prevImposed.resize(imposed.size()) ;

//    LogarithmicCreepWithImposedDeformationAndFracture::preProcess(timeStep, currentState) ;
    makeEquivalentViscosity(timeStep, currentState);

    currentState.getParent()->behaviourUpdated = true ;

}
