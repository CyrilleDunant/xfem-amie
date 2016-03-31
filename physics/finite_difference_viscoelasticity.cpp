#include "finite_difference_viscoelasticity.h"
#include "../elements/integrable_entity.h"
#include "../features/boundarycondition.h"
#include "homogenization/composite.h"

using namespace Amie ;

FiniteDifferenceViscoelasticity::FiniteDifferenceViscoelasticity(ViscoelasticModel m, const Matrix & rig, ViscoelasticFiniteDifferenceIntegration s, double t) : LinearForm(rig, false, false, (rig.numRows()/3+1)), model(m), viscoelasticInternalForces(rig.numRows()), scheme(s), theta(t)
{
    v.push_back(XI);
    v.push_back(ETA);
    if(rig.size() > 9)
        v.push_back(ZETA);

    switch(model)
    {
    case PURE_ELASTICITY:
        break ;
    case PURE_VISCOSITY:
        param = 0 ;
        tensors.push_back(rig*0) ;
        break ;
    default:
        std::cout << "warning: wrong constructor for FiniteDifferenceViscoelasticity" << std::endl ;
    }

    tensors.push_back(rig) ;

    initializeThetaCoefficient() ;

    viscoelasticInternalForces = 0 ;

}

FiniteDifferenceViscoelasticity::FiniteDifferenceViscoelasticity(ViscoelasticModel m, const Matrix & rig, const Matrix & e, ViscoelasticFiniteDifferenceIntegration s, double t) : LinearForm(rig, false, false, (rig.numRows()/3+1)), model(m), viscoelasticInternalForces(rig.numRows()), scheme(s), theta(t)
{
    v.push_back(XI);
    v.push_back(ETA);
    if(rig.size() > 9)
        v.push_back(ZETA);

    switch(model)
    {
    case KELVIN_VOIGT:
        break ;
    case MAXWELL:
        break ;
    default:
        std::cout << "warning: wrong constructor for FiniteDifferenceViscoelasticity" << std::endl ;
    }

    tensors.push_back(rig) ;
    tensors.push_back(e) ;

    initializeThetaCoefficient() ;

    viscoelasticInternalForces = 0 ;
}


FiniteDifferenceViscoelasticity::FiniteDifferenceViscoelasticity(ViscoelasticModel m, const Matrix & c_kv, const Matrix & e_kv, const Matrix & c_mx, const Matrix & e_mx, ViscoelasticFiniteDifferenceIntegration s, double t) : LinearForm(c_mx, false, false, (c_mx.numRows()/3+1)), model(m), viscoelasticInternalForces(c_mx.numRows()), scheme(s), theta(t)
{
    v.push_back(XI);
    v.push_back(ETA);
    if(c_kv.size() > 9)
        v.push_back(ZETA);

    switch(model)
    {
    case BURGER:
        break ;
    default:
        std::cout << "warning: wrong constructor for FiniteDifferenceViscoelasticity" << std::endl ;
    }

    tensors.push_back(c_mx) ;
    tensors.push_back(e_mx) ;
    tensors.push_back(c_kv) ;
    tensors.push_back(e_kv) ;

    initializeThetaCoefficient() ;

    viscoelasticInternalForces = 0 ;
}

FiniteDifferenceViscoelasticity::FiniteDifferenceViscoelasticity(ViscoelasticModel m, const Matrix & c0, std::vector<std::pair<Matrix, Matrix> > & branches, ViscoelasticFiniteDifferenceIntegration s, double t) : LinearForm(c0, false, false, (c0.numRows()/3+1)), model(m), viscoelasticInternalForces(c0.numRows()), scheme(s), theta(t)
{
    v.push_back(XI);
    v.push_back(ETA);
    if(c0.size() > 9)
        v.push_back(ZETA);

    switch(model)
    {
    case GENERALIZED_KELVIN_VOIGT:
        break ;
    case GENERALIZED_MAXWELL:
        param = c0 ;
        for(size_t i = 0 ; i < branches.size() ; i++)
            param += branches[i].first ;
        break ;
    default:
        std::cout << "warning: wrong constructor for FiniteDifferenceViscoelasticity" << std::endl ;
    }

    tensors.push_back(c0) ;
    for(size_t i = 0 ; i < branches.size() ; i++)
    {
        tensors.push_back(branches[i].first) ;
        tensors.push_back(branches[i].second) ;
    }

    initializeThetaCoefficient() ;

    viscoelasticInternalForces = 0 ;
}

FiniteDifferenceViscoelasticity::FiniteDifferenceViscoelasticity(ViscoelasticModel m, const Matrix & c0, const Matrix & c1, const Matrix & e1, ViscoelasticFiniteDifferenceIntegration s, double t) : LinearForm(c0, false, false, (c0.numRows()/3+1)), model(m), viscoelasticInternalForces(c0.numRows()), scheme(s), theta(t)
{
    v.push_back(XI);
    v.push_back(ETA);
    if(c0.size() > 9)
        v.push_back(ZETA);

    switch(model)
    {
    case GENERALIZED_KELVIN_VOIGT:
        break ;
    case GENERALIZED_MAXWELL:
        param = c0+c1 ;
        break ;
    default:
        std::cout << "warning: wrong constructor for FiniteDifferenceViscoelasticity" << std::endl ;
    }

    tensors.push_back(c0) ;
    tensors.push_back(c1) ;
    tensors.push_back(e1) ;

    initializeThetaCoefficient() ;

    viscoelasticInternalForces = 0 ;
}

FiniteDifferenceViscoelasticity::~FiniteDifferenceViscoelasticity() {}

ElementState * FiniteDifferenceViscoelasticity::createElementState( IntegrableEntity * e)
{
    int n = 0 ;
    switch(model)
    {
    case PURE_ELASTICITY:
        n = 0 ;
        break ;
    case PURE_VISCOSITY:
        n = 2 ;
        break ;
    case KELVIN_VOIGT:
        n = 2 ;
        break ;
    case MAXWELL:
        n = 2 ;
        break ;
    case BURGER:
        n = 4 ;
        break ;
    case GENERALIZED_MAXWELL:
    case GENERALIZED_KELVIN_VOIGT:
        n = (tensors.size()-1) ;
        break ;
    default:
        break ;
    }
    return new ElementStateWithInternalVariables(e, n, param.numCols()) ;
}

double getThetaForExponentialScheme( double timeStep, double charTime ) 
{
    return 1./(1.-exp(-timeStep/charTime)) - charTime/timeStep ;
}

void FiniteDifferenceViscoelasticity::preProcess( double timeStep, ElementState & currentState ) 
{
    if(timeStep < POINT_TOLERANCE)
        return ;

    VirtualMachine vm ;

    switch(model)
    {
    case PURE_ELASTICITY:
        param = tensors[0] ;
        break ;
    case PURE_VISCOSITY:
    {
        Vector strain_prev( 0., 3+3*(num_dof == 3)) ;
        Vector straindot_prev( 0., 3+3*(num_dof == 3)) ;
        if(scheme == ZIENKIEWICZ)
            theta = getThetaForExponentialScheme( timeStep, 1. ) ;
        param = tensors[1]*1./(theta*timeStep) ;
	dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, 0, strain_prev, &vm, 0) ;
	dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, 0, straindot_prev, &vm, 1) ;
        Vector tmp = strain_prev*1./(theta*timeStep) ;
        tmp += straindot_prev*(1.-theta)/theta ;
        viscoelasticInternalForces = tensors[1]*tmp ;
        break ;
    }
    case KELVIN_VOIGT:
    {
        Vector strain_prev( 0., 3+3*(num_dof == 3)) ;
        Vector straindot_prev( 0., 3+3*(num_dof == 3)) ;
        if(scheme == ZIENKIEWICZ)
            theta = getThetaForExponentialScheme( timeStep, tensors[1][0][0]/tensors[0][0][0] ) ;
        param = tensors[0] + tensors[1]*1./(theta*timeStep) ;
	dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, 0, strain_prev, &vm, 0) ;
	dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, 0, straindot_prev, &vm, 1) ;
        Vector tmp = strain_prev*1./(theta*timeStep) ;
        tmp += straindot_prev*(1.-theta)/theta ;
        viscoelasticInternalForces = tensors[1]*tmp ;
        break ;
    }
    case MAXWELL:
    {
        Vector alpha_prev( 0., 3+3*(num_dof == 3)) ;
        Vector alphadot_prev( 0., 3+3*(num_dof == 3)) ;

        if(scheme == ZIENKIEWICZ)
            theta = getThetaForExponentialScheme( timeStep, tensors[1][0][0]/tensors[0][0][0] ) ;
        Matrix S0 = tensors[0] ;
        Composite::invertTensor(S0) ;
        Matrix S1 = tensors[1]*1./(theta*timeStep) ;
        Composite::invertTensor(S1) ;
        S0 += S1 ;
        Composite::invertTensor(S0) ;
        param = S0 ;
	dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, 0, alpha_prev, &vm, 0) ;
	dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, 0, alphadot_prev, &vm, 1) ;
        Vector tmp = alpha_prev ;
        tmp += alphadot_prev*(1.-theta)*timeStep ;
        viscoelasticInternalForces = param*tmp ;
        break ;
    }
    case BURGER:
    {
        Vector alpha_prev( 0., 3+3*(num_dof == 3)) ;
        Vector alphadot_prev( 0., 3+3*(num_dof == 3)) ;

        Matrix S0 = tensors[0] ;
        Composite::invertTensor(S0) ;

        if(scheme == ZIENKIEWICZ)
            theta = getThetaForExponentialScheme( timeStep, tensors[1][0][0]/tensors[0][0][0] ) ;
        Matrix S1 = tensors[1]*1./(theta*timeStep) ;
        Composite::invertTensor(S1) ;
	dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, 0, alpha_prev, &vm, 0) ;
	dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, 0, alphadot_prev, &vm, 1) ;
        Vector strain_mx = alpha_prev+alphadot_prev*((1.-theta)*timeStep) ;


        if(scheme == ZIENKIEWICZ)
            theta = getThetaForExponentialScheme( timeStep, tensors[3][0][0]/tensors[2][0][0] ) ;
        Matrix S2 = tensors[3]+tensors[2]*(theta*timeStep) ;
        Composite::invertTensor(S2) ;
        alpha_prev = 0. ; alphadot_prev = 0. ;
	dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, 0, alpha_prev, &vm, 0) ;
	dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, 0, alphadot_prev, &vm, 1) ;
        Vector strain_kv = alpha_prev+alphadot_prev*((1.-theta)*timeStep) ;
        strain_kv = (tensors[3]*S2)*strain_kv ;

        Vector strain = strain_kv + strain_mx ;
        S0 += S1 ;
        S0 += S2 ;
        Composite::invertTensor(S0) ;
        
        param = S0 ;
        viscoelasticInternalForces = param*strain ;

        break ;
    }
    case GENERALIZED_KELVIN_VOIGT:
    {
        Vector alpha_prev( 0., 3+3*(num_dof == 3)) ;
        Vector alphadot_prev( 0., 3+3*(num_dof == 3)) ;
        Vector omega_i( 0., 3+3*(num_dof == 3)) ;
        Vector omega( 0., 3+3*(num_dof == 3)) ;
        Matrix K_i = tensors[0]*0. ;
        Matrix S = tensors[0] ;
        Composite::invertTensor(S) ;

        for(size_t i = 1 ; i < tensors.size() ; i+=2)
        {
            if(scheme == ZIENKIEWICZ)
                theta = getThetaForExponentialScheme( timeStep, tensors[i+1][0][0]/tensors[i][0][0] ) ;
            alpha_prev = 0. ; alphadot_prev = 0. ;
            dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, 0, alpha_prev, &vm, i-1) ;
            dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, 0, alphadot_prev, &vm, i) ;

            K_i = tensors[i]+tensors[i+1]*1./(theta*timeStep) ;
            Composite::invertTensor(K_i) ;

            omega_i = alpha_prev*(1./(theta*timeStep)) ;
            omega_i += alphadot_prev*(1.-theta)/theta ;

            omega += (tensors[i+1]*K_i)*omega_i ;
            S += K_i ;
        }

        Composite::invertTensor(S) ;
        param = S ;
        viscoelasticInternalForces = param*omega ;
        
        break ;
    }
    case GENERALIZED_MAXWELL:
    {
        Vector alpha_prev( 0., 3+3*(num_dof == 3)) ;
        Vector alphadot_prev( 0., 3+3*(num_dof == 3)) ;
        Vector omega_i( 0., 3+3*(num_dof == 3)) ;
        Matrix K_i = tensors[0]*0. ;

        param = tensors[0] ;
        viscoelasticInternalForces = 0 ;

        for(size_t i = 1 ; i < tensors.size() ; i+=2)
        {
            if(scheme == ZIENKIEWICZ)
                theta = getThetaForExponentialScheme( timeStep, tensors[i+1][0][0]/tensors[i][0][0] ) ;
            alpha_prev = 0. ; alphadot_prev = 0. ;
            dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, 0, alpha_prev, &vm, i-1) ;
            dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, 0, alphadot_prev, &vm, i) ;

            K_i = tensors[i]+tensors[i+1]*1./(theta*timeStep) ;
            Composite::invertTensor(K_i) ;

            omega_i = alpha_prev*(1./(theta*timeStep)) ;
            omega_i += alphadot_prev*(1.-theta)/theta ;

            param += (tensors[i] - tensors[i]*(tensors[i]*K_i)) ;
            viscoelasticInternalForces += (tensors[i]*(tensors[i+1]*K_i))*omega_i ;

        }

        break ;
    }
    case GENERAL_VISCOELASTICITY:
    {
        std::cout << "no finite difference scheme applicable for GENERAL_VISCOELASTICITY, exiting now..." << std::endl ;
        exit(0) ;
    }

    }

    currentState.getParent()->behaviourUpdated = true ;  
}


void FiniteDifferenceViscoelasticity::updateElementState(double timeStep, ElementState & currentState) const 
{
    if(timeStep < POINT_TOLERANCE)
        return ;
    VirtualMachine vm ;
    std::vector<double> coefs ;

    double localTheta = theta ;

    switch(model)
    {
    case PURE_ELASTICITY:
        break ;
    case PURE_VISCOSITY:
    {
        if(scheme == ZIENKIEWICZ)
            localTheta = getThetaForExponentialScheme( timeStep, 1. ) ;

        Vector strain_next( 0., 3+3*(num_dof == 3)) ;
        Vector strain_prev( 0., 3+3*(num_dof == 3)) ;
        Vector straindot_next( 0., 3+3*(num_dof == 3)) ;
        Vector straindot_prev( 0., 3+3*(num_dof == 3)) ;

	currentState.getAverageField( STRAIN_FIELD, strain_next, &vm, 0) ;
	dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, 0, strain_prev, &vm, 0) ;
	dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, 0, straindot_prev, &vm, 1) ;

        straindot_next = strain_next - strain_prev ;
        straindot_next /= (localTheta*timeStep) ;
        straindot_next -= straindot_prev*(1.-localTheta)/localTheta ;

	dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(strain_next, 0, 0) ;
	dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(straindot_next, 0, 1) ;

        break ;
    }
    case KELVIN_VOIGT:
    {
        Vector stress_next( 0., 3+3*(num_dof == 3)) ;
        Vector strain_next( 0., 3+3*(num_dof == 3)) ;
        Vector straindot_next( 0., 3+3*(num_dof == 3)) ;

	currentState.getAverageField( REAL_STRESS_FIELD, stress_next, &vm, 0) ;
	currentState.getAverageField( STRAIN_FIELD, strain_next, &vm, 0) ;
        stress_next -= tensors[0]*strain_next ;

        Matrix S1 = tensors[1] ;
        Composite::invertTensor(S1) ;
        straindot_next = S1*stress_next ;

	dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(strain_next, 0, 0) ;
	dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(straindot_next, 0, 1) ;

        break ;
    }
    case MAXWELL:
    {
        if(scheme == ZIENKIEWICZ)
            localTheta = getThetaForExponentialScheme( timeStep, tensors[1][0][0]/tensors[0][0][0] ) ;

        Vector alpha_prev( 0., 3+3*(num_dof == 3)) ;
        Vector alpha_next( 0., 3+3*(num_dof == 3)) ;
        Vector alphadot_prev( 0., 3+3*(num_dof == 3)) ;
        Vector alphadot_next( 0., 3+3*(num_dof == 3)) ;
        Vector stress_next( 0., 3+3*(num_dof == 3)) ;

	currentState.getAverageField( REAL_STRESS_FIELD, stress_next, &vm, 0) ;
	dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, 0, alpha_prev,&vm, 0) ;
	dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, 0, alphadot_prev,&vm, 1) ;

        Matrix S1 = tensors[1] ;
        Composite::invertTensor(S1) ;
        alphadot_next = S1*stress_next ;

        alpha_next = alpha_prev ;
        alpha_next += alphadot_next*localTheta*timeStep ;
        alpha_next += alphadot_prev*(1.-localTheta)*timeStep ;

	dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(alpha_next, 0, 0) ;
	dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(alphadot_next, 0, 1) ;

        break ;
    }
    case BURGER:
    {

        Vector alpha_prev( 0., 3+3*(num_dof == 3)) ;
        Vector alpha_next( 0., 3+3*(num_dof == 3)) ;
        Vector alphadot_prev( 0., 3+3*(num_dof == 3)) ;
        Vector alphadot_next( 0., 3+3*(num_dof == 3)) ;
        Vector stress_next( 0., 3+3*(num_dof == 3)) ;

        if(scheme == ZIENKIEWICZ)
            localTheta = getThetaForExponentialScheme( timeStep, tensors[1][0][0]/tensors[0][0][0] ) ;

	currentState.getAverageField( REAL_STRESS_FIELD, stress_next, &vm, 0) ;
	dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, 0, alpha_prev,&vm, 0) ;
	dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, 0, alphadot_prev,&vm, 1) ;

        Matrix S1 = tensors[1] ;
        Composite::invertTensor(S1) ;
        alphadot_next = S1*stress_next ;

        alpha_next = alpha_prev ;
        alpha_next += alphadot_next*localTheta*timeStep ;
        alpha_next += alphadot_prev*(1.-localTheta)*timeStep ;

	dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(alpha_next, 0, 0) ;
	dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(alphadot_next, 0, 1) ;



        if(scheme == ZIENKIEWICZ)
            localTheta = getThetaForExponentialScheme( timeStep, tensors[3][0][0]/tensors[2][0][0] ) ;

        alpha_prev = 0 ;
        alpha_next = 0 ;
        alphadot_prev = 0 ;
        alphadot_next = 0 ;
	dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, 0, alpha_prev,&vm, 2) ;
	dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, 0, alphadot_prev,&vm, 3) ;

        Matrix S2 = tensors[3]+tensors[2]*(theta*timeStep) ;
        Composite::invertTensor(S2) ;

        alpha_next = stress_next*theta*timeStep ;
        alpha_next += tensors[3]*(alpha_prev+alphadot_prev*(1.-theta)*timeStep) ;
        alpha_next = S2*alpha_next ;

        alphadot_next = (alpha_next-alpha_prev)/(theta*timeStep) ;
        alphadot_next -= alphadot_prev*(1.-theta)/theta ;

	dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(alpha_next, 0, 2) ;
	dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(alphadot_next, 0, 3) ;

        break ;
    }
    case GENERALIZED_KELVIN_VOIGT:
    {
        Vector alpha_prev( 0., 3+3*(num_dof == 3)) ;
        Vector alpha_next( 0., 3+3*(num_dof == 3)) ;
        Vector alphadot_prev( 0., 3+3*(num_dof == 3)) ;
        Vector alphadot_next( 0., 3+3*(num_dof == 3)) ;
        Vector stress_next( 0., 3+3*(num_dof == 3)) ;
        Vector omega_i( 0., 3+3*(num_dof == 3)) ;
        Matrix K_i = tensors[0]*0 ;

	currentState.getAverageField( REAL_STRESS_FIELD, stress_next, &vm, 0) ;
        for(size_t i = 1 ; i < tensors.size() ; i+=2)
        {
            if(scheme == ZIENKIEWICZ)
                localTheta = getThetaForExponentialScheme( timeStep, tensors[i+1][0][0]/tensors[i][0][0] ) ;
            alpha_prev = 0. ; alphadot_prev = 0. ;
            dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, 0, alpha_prev, &vm, i-1) ;
            dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, 0, alphadot_prev, &vm, i) ;

            K_i = tensors[i]+tensors[i+1]*1./(localTheta*timeStep) ;
            Composite::invertTensor(K_i) ;

            omega_i = alpha_prev*(1./(localTheta*timeStep)) ;
            omega_i += alphadot_prev*(1.-localTheta)/localTheta ;

            alpha_next = stress_next ;
            alpha_next += tensors[i+1]*omega_i ;
            alpha_next = K_i * alpha_next ;

            alphadot_next = (alpha_next-alpha_prev)/(localTheta*timeStep) ;
            alphadot_next -= alphadot_prev*(1.-localTheta)/localTheta ;

            dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(alpha_next, 0, i-1) ;
            dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(alphadot_next, 0, i) ;

        }
        break ;

    }
    case GENERALIZED_MAXWELL:
    {
        Vector alpha_prev( 0., 3+3*(num_dof == 3)) ;
        Vector alpha_next( 0., 3+3*(num_dof == 3)) ;
        Vector alphadot_prev( 0., 3+3*(num_dof == 3)) ;
        Vector alphadot_next( 0., 3+3*(num_dof == 3)) ;
        Vector strain_next( 0., 3+3*(num_dof == 3)) ;
        Vector omega_i( 0., 3+3*(num_dof == 3)) ;
        Matrix K_i = tensors[0]*0 ;

	currentState.getAverageField( STRAIN_FIELD, strain_next, &vm, 0) ;
        for(size_t i = 1 ; i < tensors.size() ; i+=2)
        {
            if(scheme == ZIENKIEWICZ)
                localTheta = getThetaForExponentialScheme( timeStep, tensors[i+1][0][0]/tensors[i][0][0] ) ;
            alpha_prev = 0. ; alphadot_prev = 0. ;
            dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, 0, alpha_prev, &vm, i-1) ;
            dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, 0, alphadot_prev, &vm, i) ;

            K_i = tensors[i]+tensors[i+1]*1./(localTheta*timeStep) ;
            Composite::invertTensor(K_i) ;

            omega_i = alpha_prev*(1./(localTheta*timeStep)) ;
            omega_i += alphadot_prev*(1.-localTheta)/localTheta ;

            alpha_next = tensors[i]*strain_next ;
            alpha_next += tensors[i+1]*omega_i ;
            alpha_next = K_i * alpha_next ;

            alphadot_next = (alpha_next-alpha_prev)/(localTheta*timeStep) ;
            alphadot_next -= alphadot_prev*(1.-localTheta)/localTheta ;

            dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(alpha_next, 0, i-1) ;
            dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(alphadot_next, 0, i) ;

        }
        break ;

    }
    case GENERAL_VISCOELASTICITY:
    {
        std::cout << "no finite difference scheme applicable for GENERAL_VISCOELASTICITY, exiting now..." << std::endl ;
        exit(0) ;
    }



    }
        
}


void FiniteDifferenceViscoelasticity::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
    vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v, ret) ;
}

Form * FiniteDifferenceViscoelasticity::getCopy() const
{
    switch(model)
    {
    case PURE_ELASTICITY:
        return new FiniteDifferenceViscoelasticity( PURE_ELASTICITY, tensors[0], scheme, theta) ;
    case PURE_VISCOSITY:
        return new FiniteDifferenceViscoelasticity( PURE_VISCOSITY, tensors[1], scheme, theta) ;
    case MAXWELL:
        return new FiniteDifferenceViscoelasticity( MAXWELL, tensors[0], tensors[1], scheme, theta) ;
    case KELVIN_VOIGT:
        return new FiniteDifferenceViscoelasticity( KELVIN_VOIGT, tensors[0], tensors[1], scheme, theta) ;
    case BURGER:
        return new FiniteDifferenceViscoelasticity( BURGER, tensors[2], tensors[3], tensors[0], tensors[1], scheme, theta) ;
    case GENERALIZED_MAXWELL:
    {
        std::vector<std::pair<Matrix, Matrix> > branches ;
        for(size_t i = 1 ; i < tensors.size() ; i+= 2)
        {
            branches.push_back( std::make_pair(tensors[i], tensors[i+1]) ) ;
        }
        return new FiniteDifferenceViscoelasticity( GENERALIZED_MAXWELL, tensors[0], branches, scheme, theta) ;
    }
    case GENERALIZED_KELVIN_VOIGT:
    {
        std::vector<std::pair<Matrix, Matrix> > branches ;
        for(size_t i = 1 ; i < tensors.size() ; i+= 2)
        {
            branches.push_back( std::make_pair(tensors[i], tensors[i+1]) ) ;
        }
        return new FiniteDifferenceViscoelasticity( GENERALIZED_KELVIN_VOIGT, tensors[0], branches, scheme, theta) ;
    }
    case GENERAL_VISCOELASTICITY:
        break ;
    }

    return new FiniteDifferenceViscoelasticity(PURE_ELASTICITY, param, scheme, theta) ;
}

Vector FiniteDifferenceViscoelasticity::getImposedStress(const Point & p, IntegrableEntity * e, int g) const
{
    return viscoelasticInternalForces ;
}

Vector FiniteDifferenceViscoelasticity::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const
{
    Matrix compliance = param ;
    if(compliance.array().max() == 0)
    {
         Vector ret(viscoelasticInternalForces.size()) ;
         ret = 0 ;
         return ret ;
    }

    Composite::invertTensor(compliance) ;
             
    return compliance*viscoelasticInternalForces ;
}

std::vector<BoundaryCondition * > FiniteDifferenceViscoelasticity::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{

    std::vector<BoundaryCondition * > ret ;
    if(v.size() == 2)
    {
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, viscoelasticInternalForces[0]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, viscoelasticInternalForces[1]));
// 		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, imposed[2]));
    }
    if(v.size() == 3)
    {
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_XI, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, viscoelasticInternalForces[0]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, viscoelasticInternalForces[1]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, viscoelasticInternalForces[2]));
// 		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, imposed[3]));
// 		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, imposed[4]));
// 		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_ETA_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, imposed[5]));
    }
    return ret ;
}


void FiniteDifferenceViscoelasticity::print() const
{
    std::cout << "I am a finite-difference-based viscoelastic model" ;
    switch(model)
    {
    case PURE_ELASTICITY:
        std::cout << " (elastic only)" << std::endl ;
        return ;
    case PURE_VISCOSITY:
        std::cout << " (viscous only)" << std::endl ;
        return ;
    case KELVIN_VOIGT:
        std::cout << " (kelvin-voigt)" << std::endl ;
        return ;
    case GENERALIZED_KELVIN_VOIGT:
        std::cout << " ( generalized kelvin-voigt)" << std::endl ;
        return ;
    case MAXWELL:
        std::cout << " (maxwell)" << std::endl ;
        return ;
    case GENERALIZED_MAXWELL:
        std::cout << " ( generalized maxwell)" << std::endl ;
        return ;
    case BURGER:
        std::cout << " (burger)" << std::endl ;
        return ;
    case GENERAL_VISCOELASTICITY:
        std::cout << " (general)" << std::endl ;
        return ;
    }
}

void FiniteDifferenceViscoelasticity::initializeThetaCoefficient() 
{
    switch(scheme)
    {
    case FORWARD_EULER:
        theta = 0 ;
        break ;
    case BACKWARD_EULER:
        theta = 1 ;
        break ;
    case CENTRAL_DIFFERENCE:
        theta = 0.5 ;
        break ;
    case NEWMARK:
        if(theta < 0)
             theta = 0 ;
        if(theta > 1)
             theta = 1 ;
        break ;
    case ZIENKIEWICZ:
        theta = 0.5 ;
        break ;
    }
}



