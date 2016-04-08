#include "viscoelasticity_and_fracture.h"
#include "viscoelasticity.h"
#include "fracturecriteria/maxstrain.h"
#include "../elements/generalized_spacetime_viscoelastic_element_state.h"

using namespace Amie ;

ViscoelasticityAndFracture::ViscoelasticityAndFracture(ViscoelasticModel m, const Matrix & rig, FractureCriterion * c, DamageModel * d, int n, double r) : Viscoelasticity(m, rig, n, r), dfunc(d), criterion(c)
{
    setElasticAndViscousStiffnessMatrix() ;
}

ViscoelasticityAndFracture::ViscoelasticityAndFracture(ViscoelasticModel m, const Matrix & rig, const Matrix & e, FractureCriterion * c, DamageModel * d, int b, int n, double r) : Viscoelasticity(m, rig, e, b, n, r), dfunc(d), criterion(c)
{
    setElasticAndViscousStiffnessMatrix() ;
}

ViscoelasticityAndFracture::ViscoelasticityAndFracture(const Matrix & rig, const Matrix & e, int b, FractureCriterion * c, DamageModel * d, int n, double r) : Viscoelasticity(rig, e, b, n, r), dfunc(d) , criterion(c)
{
    setElasticAndViscousStiffnessMatrix() ;
}

ViscoelasticityAndFracture::ViscoelasticityAndFracture(ViscoelasticModel m, const Matrix & c_kv, const Matrix & e_kv, const Matrix & c_mx, const Matrix & e_mx,FractureCriterion * c, DamageModel * d,  int b, int n, double r) : Viscoelasticity(m, c_kv, e_kv, c_mx, e_mx, b, n, r), dfunc(d), criterion(c)
{
    setElasticAndViscousStiffnessMatrix() ;
}

ViscoelasticityAndFracture::ViscoelasticityAndFracture(ViscoelasticModel m, const Matrix & c0, std::vector<std::pair<Matrix, Matrix> > & branches, FractureCriterion * c, DamageModel * d, int b, int n, double r) : Viscoelasticity(m, c0, branches, b, n, r), dfunc(d), criterion(c)
{
    setElasticAndViscousStiffnessMatrix() ;

}

ViscoelasticityAndFracture::ViscoelasticityAndFracture(ViscoelasticModel m, const Matrix & c0, const Matrix & c1, const Matrix & e1, FractureCriterion * c, DamageModel * d, int b, int n, double r) : Viscoelasticity(m, c0, c1, e1, b, n, r), dfunc(d), criterion(c)
{
    setElasticAndViscousStiffnessMatrix() ;
}

ViscoelasticityAndFracture::ViscoelasticityAndFracture( ViscoelasticModel model, double young, double poisson, FractureCriterion * c, DamageModel * d, double tau, std::string file, SpaceDimensionality dim, planeType pt, IsotropicMaterialParameters hooke, int b, int a) : Viscoelasticity(model, young, poisson, tau, file, dim, pt, hooke, b, a), dfunc(d), criterion(c)
{
    setElasticAndViscousStiffnessMatrix() ;
}

ViscoelasticityAndFracture::~ViscoelasticityAndFracture()
{
    delete dfunc ;
    delete criterion ;
}

ElementState * ViscoelasticityAndFracture::createElementState( IntegrableEntity * e)
{
    return new GeneralizedSpaceTimeViscoElasticElementState(e, blocks) ;
}

FractureCriterion * ViscoelasticityAndFracture::getFractureCriterion() const
{
    return criterion ;
}

DamageModel * ViscoelasticityAndFracture::getDamageModel() const
{
    return dfunc ;
}


void ViscoelasticityAndFracture::step(double timestep, ElementState & currentState, double maxscore)
{
    dfunc->step(currentState, maxscore) ;
    currentState.getParent()->behaviourUpdated = dfunc->changed() ;
    currentState.getParent()->needAssembly = dfunc->changed() ; 
}


bool ViscoelasticityAndFracture::fractured() const
{
    return dfunc->fractured() ;
}

bool ViscoelasticityAndFracture::changed() const
{
    return dfunc->changed() ;
}

Form * ViscoelasticityAndFracture::getCopy() const
{

    ViscoelasticityAndFracture * copy ;
    switch(model)
    {
    case PURE_ELASTICITY:
        copy = new ViscoelasticityAndFracture( PURE_ELASTICITY, tensors[0], criterion->getCopy(), dfunc->getCopy(), blocks-1) ;
        break ;
    case PURE_VISCOSITY:
        copy = new ViscoelasticityAndFracture( PURE_VISCOSITY, tensors[1], criterion->getCopy(), dfunc->getCopy(), blocks-1) ;
        break ;
    case MAXWELL:
        copy = new ViscoelasticityAndFracture( MAXWELL, tensors[0], tensors[1], criterion->getCopy(), dfunc->getCopy(), blocks-2) ;
        break ;
    case KELVIN_VOIGT:
        copy = new ViscoelasticityAndFracture( KELVIN_VOIGT, tensors[0], tensors[1], criterion->getCopy(), dfunc->getCopy(), blocks-1) ;
        break ;
    case GENERALIZED_MAXWELL:
    {
        std::vector<std::pair<Matrix, Matrix> > branches ;
        for(size_t i = 1 ; i < tensors.size() ; i += 2)
            branches.push_back(std::make_pair(tensors[i], tensors[i+1])) ;
        copy = new ViscoelasticityAndFracture( GENERALIZED_MAXWELL, tensors[0], branches, criterion->getCopy(), dfunc->getCopy()) ;
        break ;
    }
    case GENERALIZED_KELVIN_VOIGT:
    {
        std::vector<std::pair<Matrix, Matrix> > branches ;
        for(size_t i = 1 ; i < tensors.size() ; i += 2)
            branches.push_back(std::make_pair(tensors[i], tensors[i+1])) ;
        copy = new ViscoelasticityAndFracture( GENERALIZED_KELVIN_VOIGT, tensors[0], branches, criterion->getCopy(), dfunc->getCopy()) ;
        break ;
    }
    case BURGER:
        copy = new ViscoelasticityAndFracture( BURGER, tensors[0], tensors[1], tensors[2], tensors[3], criterion->getCopy(), dfunc->getCopy(), blocks-3) ;
        break ;
    default:
        copy = new ViscoelasticityAndFracture(  param, eta, blocks, criterion->getCopy(), dfunc->getCopy())  ;

    }
    copy->model = model ;
    copy->dfunc->getState(true).resize(dfunc->getState().size());
    copy->dfunc->getState(true) = dfunc->getState() ;
    copy->criterion->setMaterialCharacteristicRadius(criterion->getMaterialCharacteristicRadius()) ;
    copy->dfunc->setDamageDensityTolerance(dfunc->getDamageDensityTolerance());
    copy->dfunc->setThresholdDamageDensity(dfunc->getThresholdDamageDensity());

    return copy ;
}

Vector ViscoelasticityAndFracture::getForcesFromAppliedStress( const Vector & data, Function & shape, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<Variable> & v, bool isVolumic, const Vector & normal )
{
    return VirtualMachine().ieval(GradientDot( shape ) * ( data ), gp, Jinv, v) ;
}

Vector ViscoelasticityAndFracture::getForcesFromAppliedStress( const Function & data, size_t index, size_t externaldofs,  Function & shape, IntegrableEntity * e,const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<Variable> & v, bool isVolumic , const Vector & normal)
{
    VirtualMachine vm ;

    size_t n = e->getBoundingPoints().size() ;
    Vector field(0., n*externaldofs) ;
    for(size_t i = 0 ; i < n ; i++)
        field[ i*externaldofs + index ] = vm.eval( data, e->getBoundingPoint(i) ) ;

    std::vector<Vector> g(e->getGaussPoints().gaussPoints.size(), Vector(0., externaldofs)) ;
    e->getState().getExternalFieldAtGaussPoints( field, externaldofs, g) ;

    Vector f = vm.ieval( GradientDot( shape ) * g, gp, Jinv, v) ;

    field = 0 ;
    for(size_t i = 0 ; i < n ; i++)
        field[ i*externaldofs + index ] = vm.deval( data, TIME_VARIABLE, e->getBoundingPoint(i) ) ;

    e->getState().getExternalFieldAtGaussPoints( field, externaldofs, g) ;

    f += vm.ieval( Gradient( shape ) * g, gp, Jinv, v) ;

    return f ;
}

void ViscoelasticityAndFracture::setElasticAndViscousStiffnessMatrix()
{
    elasticParam.resize( param.numRows(), param.numCols()) ;
    viscousParam.resize( param.numRows(), param.numCols()) ;
    switch(model)
    {
    case PURE_ELASTICITY:
        placeMatrixInBlock( tensors[0], 0,0, elasticParam) ;
        break ;
    case PURE_VISCOSITY:
        break ;
    case KELVIN_VOIGT:
        placeMatrixInBlock( tensors[0], 0,0, elasticParam) ;
        break ;
    case MAXWELL:
        placeMatrixInBlock( tensors[0], 0,0, elasticParam) ;
        placeMatrixInBlock( tensors[0], 1,0, elasticParam) ;
        placeMatrixInBlock( -tensors[0], 0,1, elasticParam) ;
        placeMatrixInBlock( -tensors[0], 1,1, elasticParam) ;
        break ;
    case BURGER:
        placeMatrixInBlock( tensors[0], 0,0, elasticParam) ;
        placeMatrixInBlock( -tensors[0], 1,0, elasticParam) ;
        placeMatrixInBlock( -tensors[0], 2,0, elasticParam) ;
        placeMatrixInBlock( -tensors[0], 0,1, elasticParam) ;
        placeMatrixInBlock( -tensors[0], 0,2, elasticParam) ;
        placeMatrixInBlock( tensors[0], 1,1, elasticParam) ;
        placeMatrixInBlock( tensors[0], 1,2, elasticParam) ;
        placeMatrixInBlock( tensors[0], 2,1, elasticParam) ;
        placeMatrixInBlock( tensors[0], 2,2, elasticParam) ;
        addMatrixInBlock( tensors[2], 2,2, viscousParam) ;
        break ;
    case GENERALIZED_MAXWELL:
        placeMatrixInBlock( tensors[0], 0,0, elasticParam) ;
        for(size_t i = 1 ; i < tensors.size() ; i+=2)
        {
            Matrix ri = tensors[i] * (-1) ;
            addMatrixInBlock( tensors[i], 0,0, elasticParam) ;
            placeMatrixInBlock( ri, 0,i/2+1, elasticParam) ;
            placeMatrixInBlock( ri, i/2+1,0, elasticParam) ;
            placeMatrixInBlock( tensors[i], i/2+1,i/2+1, viscousParam) ;
        }
        break ;
    case GENERALIZED_KELVIN_VOIGT:
        placeMatrixInBlock( tensors[0], 0,0, elasticParam) ;
        for(size_t i = 1 ; i < tensors.size() ; i+=2)
        {
            placeMatrixInBlock( tensors[0], i/2+1,i/2+1, elasticParam) ;
            placeMatrixInBlock( -tensors[0], i/2+1,0, elasticParam) ;
            placeMatrixInBlock( -tensors[0], 0,i/2+1, elasticParam) ;
            for(size_t j = 1 ; j < i/2+1 ; j++)
            {
                placeMatrixInBlock( tensors[0], i/2+1,j, elasticParam) ;
                placeMatrixInBlock( tensors[0], j,i/2+1, elasticParam) ;
            }
            placeMatrixInBlock( tensors[i], i/2+1,i/2+1, viscousParam) ;
        }
        break ;
    default:
        elasticParam=param ;
        break ;
    }

}

Matrix ViscoelasticityAndFracture::getTensor(const Point & p, IntegrableEntity * e, int g) const
{
    if(dfunc->getState().size() == 1)
        return  dfunc->apply(elasticParam, p) + dfunc->applyViscous(viscousParam, p) ;

    Matrix ret = param ;
    Matrix tmpParam( tensors[0].numRows(), tensors[0].numCols() ) ;

    for(size_t i = 0 ; i < connectivity.size() ; i++)
    {
        getBlockInMatrix(param, connectivity[i].xplus[0], connectivity[i].yplus[0], tmpParam) ;
        tmpParam = dfunc->apply( tmpParam ) ;
        for(size_t j = 0 ; j < connectivity[i].xplus.size() ; j++)
            placeMatrixInBlock( tmpParam, connectivity[i].xplus[j], connectivity[i].yplus[j], ret ) ;
        for(size_t j = 0 ; j < connectivity[i].xminus.size() ; j++)
            placeMatrixInBlock( -tmpParam , connectivity[i].xminus[j], connectivity[i].yminus[j], ret ) ;
    }
    return ret ;
}

Matrix ViscoelasticityAndFracture::getViscousTensor(const Point & p, IntegrableEntity * e, int g) const
{
    if(dfunc->getState().size() == 1)
        return  dfunc->apply(eta, p)  ;

    Matrix ret = eta ;
    Matrix tmpParam( tensors[0].numRows(), tensors[0].numCols() ) ;

    for(size_t i = 0 ; i < connectivityViscous.size() ; i++)
    {
        getBlockInMatrix(eta, connectivityViscous[i].xplus[0], connectivityViscous[i].yplus[0], tmpParam) ;
        tmpParam = dfunc->apply( tmpParam ) ;
        for(size_t j = 0 ; j < connectivityViscous[i].xplus.size() ; j++)
            placeMatrixInBlock( tmpParam, connectivityViscous[i].xplus[j], connectivityViscous[i].yplus[j], ret ) ;
        for(size_t j = 0 ; j < connectivityViscous[i].xminus.size() ; j++)
            placeMatrixInBlock( -tmpParam , connectivityViscous[i].xminus[j], connectivityViscous[i].yminus[j], ret ) ;
    }
    return ret ;
}

void ViscoelasticityAndFracture::setFractureCriterion(FractureCriterion * frac)
{
    if(frac)
    {
        criterion = frac ;
    }

}

Vector ViscoelasticityAndFracture::getImposedStress(const Point & p, IntegrableEntity * e, int g) const
{
    if(dfunc && dfunc->hasInducedForces()) { return dfunc->getImposedStress(p) ; }
    return Vector(double(0), getTensor(p, e).numCols()/blocks) ;
}

Vector ViscoelasticityAndFracture::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const
{
    if(dfunc && dfunc->hasInducedForces()) { return dfunc->getImposedStrain(p) ; }
    return Vector(double(0), getTensor(p).numCols()/blocks) ;
}


