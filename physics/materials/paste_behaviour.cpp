// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2013
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "paste_behaviour.h"
#include "../stiffness_and_fracture.h"
#include "../stiffness.h"
#include "../viscoelasticity.h"
#include "../viscoelasticity_and_fracture.h"
#include "../fracturecriteria/mohrcoulomb.h"
#include "../damagemodels/fiberbasedisotropiclineardamage.h"
#include "../fracturecriteria/maxstrain.h"
#include "../fracturecriteria/creeprupture.h"
#include "../damagemodels/spacetimefiberbasedisotropiclineardamage.h"
#include <cstdio>

using namespace Amie ;

PasteBehaviour::PasteBehaviour(bool e, bool st, double E, double nu, double up, double e10, double e300, SpaceDimensionality dim, planeType pt, double r, double var, int b) : WeibullDistributedStiffness(E,nu, dim, -8000.*up, up, pt, var, r, YOUNG_POISSON), spaceTime(st), elastic(e), freeblocks(b), eta10(e10), eta300(e300)
{

}

Form * PasteBehaviour::getCopy() const
{
    std::default_random_engine generator(std::rand());
    std::weibull_distribution< double > distribution(5, 1);
    double weib = distribution(generator) ;
    double factor = 1. - variability + variability*weib ;

    Form * copy ;
    if(spaceTime)
    {
        Matrix C0 = param*factor ;
        Matrix C1 = C0*eta10 ;
        Matrix C2 = C0*eta300 ;
        Matrix E1 = C1*10. ;
        Matrix E2 = C2*300. ;
    
        std::vector<std::pair<Matrix, Matrix> > branches ;
        branches.push_back(std::make_pair(C1,E1));
        branches.push_back(std::make_pair(C2,E2));

        if(elastic)
        {
            copy = new Viscoelasticity(GENERALIZED_KELVIN_VOIGT, C0, branches,0,freeblocks) ;
        }
        else
        {
            SpaceTimeFiberBasedIsotropicLinearDamage * dampaste = new SpaceTimeFiberBasedIsotropicLinearDamage( 0.025, 1e-9, 0.8 ) ;
            dampaste->setLogitViscousDamageLaw(0.025, 0.3, 2.5) ;
            copy = new ViscoelasticityAndFracture(GENERALIZED_KELVIN_VOIGT, C0, branches, new SpaceTimeNonLocalMaximumStress(up*factor), dampaste, 0, freeblocks) ;
            copy->getFractureCriterion()->setMaterialCharacteristicRadius(materialRadius);
            copy->getFractureCriterion()->setScoreTolerance(1e-4) ;
            copy->getDamageModel()->setDamageDensityTolerance(0.01) ;
            copy->getDamageModel()->setThresholdDamageDensity(0.8) ;
            copy->getDamageModel()->setNeedGlobalMaximumScore(true) ;
        }

    }
    else
    {
        if(elastic)
        {
            copy = new Stiffness(param*factor) ;
        }
        else
        {
//             copy = new StiffnessAndFracture(param*factor, new /*LinearlyDecreasing*/NonLocalMohrCoulomb(up*factor,-8000.*up*factor/*,up*factor*10/E,-8000.*up*factor*10/E,E*/), new FiberBasedIsotropicLinearDamage(0.1,0.699)) ;
            copy = new StiffnessAndFracture(param*factor, new NonLocalMohrCoulomb(up*factor,-8000.*up*factor), new FiberBasedIsotropicLinearDamage(0.1,0.699)) ;
            copy->getFractureCriterion()->setMaterialCharacteristicRadius(materialRadius);
        }
    }

    return copy ;
}


LogCreepPasteBehaviour::LogCreepPasteBehaviour(bool e, double E, double nu, double up, double eta, double tau, SpaceDimensionality dim, planeType pt, double r, double var) : LogarithmicCreepWithExternalParameters("young_modulus=12e9,poisson_ration=0.2,creep_modulus=40e9,creep_poisson=0.2,creep_characteristic_time=2", new RealTimeLogCreepAccumulator(), dim, pt)
{
    external["young_modulus"] = E ;
    external["poisson_ratio"] = nu ;
    external["creep_modulus"] = eta ;
    external["creep_poisson"] = nu ;
    external["recoverable_modulus"] = eta ;
    external["recoverable_poisson"] = nu ;
    external["creep_characteristic_time"] = tau ;
    external["weibull_variability"] = var ;
    external["material_characteristic_radius"] = r ;
    if(!e && up > 0)
    {
        external["tensile_strength"] = up ;
        noFracture = false ;
        dfunc = new SpaceTimeFiberBasedIsotropicLinearDamage( 0.025, 1e-9, 0.8 ) ;
        criterion = new SpaceTimeNonLocalMaximumStress(up) ;
        criterion->setMaterialCharacteristicRadius(r);
        criterion->setScoreTolerance(1e-4) ;
        dfunc->setDamageDensityTolerance(0.01) ;
        dfunc->setThresholdDamageDensity(0.8) ;
        dfunc->setNeedGlobalMaximumScore(true) ;
    }
    else
        noFracture = true ;

    this->makeProperties( external ) ;

    this->addMaterialLaw( new WeibullDistributedMaterialLaw("young_modulus", "weibull") ) ;
}



HydratingMechanicalCementPaste::HydratingMechanicalCementPaste(FeatureTree * diffusionTree) :LinearForm(Tensor::cauchyGreen(1., 0.2,SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON), true, true, 2),  diffusionTree(diffusionTree)
{
    v.push_back(XI);
    v.push_back(ETA);

    param = getMechanicalProperties(0.999,0.3) ;
    imposed.resize(3) ;
    imposed = getAutogeneousDeformation(0.999, 0.3) ;
}

Form * HydratingMechanicalCementPaste::getCopy() const
{
    return new HydratingMechanicalCementPaste(diffusionTree) ;
}

void HydratingMechanicalCementPaste::step(double timestep, ElementState & currentState, double maxscore)
{
    DelaunayTriangle * diffusionElement = dynamic_cast<DelaunayTriangle *>(diffusionTree->get2DMesh()->getInTree(dynamic_cast<DelaunayTriangle *>(currentState.getParent())->index)) ;
    Vector saturation = diffusionElement->getState().getDisplacements() ;
    double effectiveSaturation = (saturation[0]+saturation[1]+saturation[2])/3. ;
    double doh = dynamic_cast<HydratingDiffusionCementPaste *>(diffusionElement->getBehaviour())->getDegreeOfHydration() ;

    param = getMechanicalProperties(effectiveSaturation,doh) ;
    imposed = getAutogeneousDeformation(effectiveSaturation, doh) ;

    currentState.getParent()->behaviourUpdated = true ;
}

void HydratingMechanicalCementPaste::makeBulkModuli(double effectiveSaturation, double doh)
{
    double K [10] = { 5326.7247, - 11043.576, 1084.6109, 31480.598, - 9074.1705, 84734.677, - 20941.74, 3925.6201, - 15481.504, - 66958.406  };

    bulk = K[0] + K[1]*effectiveSaturation + K[2]*doh + K[3]*effectiveSaturation*effectiveSaturation + K[4]*doh*effectiveSaturation + K[5]*doh*doh + K[6]*effectiveSaturation*effectiveSaturation*effectiveSaturation + K[7]*doh*effectiveSaturation*effectiveSaturation + K[8]*doh*doh*effectiveSaturation + K[9]*doh*doh*doh;

    if (bulk<0)
    {
        bulk=0;
    }
    bulk *= 1e6 ;

    bulkSolid = K[0] + K[2]*doh + K[5]*doh*doh + K[9]*doh*doh*doh;
    if (bulkSolid<0)
    {
        bulkSolid=0;
    }
    bulkSolid *= 1e6 ;
}

Vector HydratingMechanicalCementPaste::getImposedStress(const Point & p, IntegrableEntity * e, int g) const
{
    return (param * imposed) ;
}

Vector HydratingMechanicalCementPaste::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const
{
    return imposed ;
}

std::vector<BoundaryCondition * > HydratingMechanicalCementPaste::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
    std::vector<BoundaryCondition * > ret ;
    Vector istress = param * imposed   ;
    ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, istress[0]));
    ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, istress[1]));
    return ret ;
}

Matrix HydratingMechanicalCementPaste::getMechanicalProperties(double effectiveSaturation, double doh)
{
    double C [10]= { 6328.4269, - 5426.419, 20602.939, - 13932.055, 19700.509, 105229.09, 25708.722, - 52055.068, 47843.694, - 136323.96 };

    //in MPa
    double E = C[0] + C[1]*effectiveSaturation + C[2]*doh + C[3]*effectiveSaturation*effectiveSaturation + C[4]*doh*effectiveSaturation + C[5]*doh*doh + C[6]*effectiveSaturation*effectiveSaturation*effectiveSaturation + C[7]*doh*effectiveSaturation*effectiveSaturation + C[8]*doh*doh*effectiveSaturation + C[9]*doh*doh*doh;

    if ( E < 0 )
    {
        E = 0;
    }

    E *= 1e6 ;

    makeBulkModuli(effectiveSaturation, doh) ;

    double nu = (3.*bulk-E)/(6.*bulk) ;


    return Tensor::cauchyGreen(E, nu, SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON) ;
}

double HydratingMechanicalCementPaste::getCapillaryPressure(double saturation, double doh)
{
    double m = 0.488;
    double a = -39.623*std::pow(( 1. - std::exp(-25.057*doh) ) , 8.15) ;
    a *= 1e6 ;
    return a*std::pow((std::pow(saturation,-1./m)-1.),(1.-m));
}

Vector HydratingMechanicalCementPaste::getAutogeneousDeformation(double saturation, double doh)
{
    makeBulkModuli(saturation, doh) ;
    Vector a(0., 3) ;

    double deformation = 0 ; //pc*saturation*(1./(3*bulk) - 1./(3*bulkSolid)) ;
    if(std::isnan(deformation))
        deformation = 0 ;

    a[0] = deformation ;
    a[1] = deformation ;
    a[2] = 0. ;

    return a ;
}


void HydratingMechanicalCementPaste::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{

    vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v, ret) ;
}

bool HydratingMechanicalCementPaste::fractured() const {
    return false ;
}

HydratingMechanicalCementPaste::~HydratingMechanicalCementPaste() {}

bool HydratingMechanicalCementPaste::changed() const {
    return false ;
 
}

Matrix HydratingMechanicalCementPaste::getTensor(const Point & p, IntegrableEntity * e, int g) const
{
    return param ;
}


HydratingDiffusionCementPaste::HydratingDiffusionCementPaste(): LinearForm(Matrix(2,2), true, false, 1)
{
    v.push_back(XI);
    v.push_back(ETA);
    v.push_back(TIME_VARIABLE);
    doh = 0.3 ;
    change = false ;
    double D1_bz= (-5.412*doh*doh - 0.987*doh + 3.713)*1e-10;     //[m2/s], Duffison coefficient in function of doh

    double a = 0.04 ;
    int n = 6;
    double Abz=0.9002*doh-0.1503;    //Bolzmann function parameter in function of doh
    double Bbz=0.2815*doh+0.0297;
    double h_b=1.-1./(1.+exp((1.-Abz)/Bbz));   //Bolzmann function for inversed desorption isotherm

    param[0][0]= D1_bz*(a + (1.-a)*1./(1.+pow((1.-h_b)/0.25,n)))*24*3600;
    param[1][1]= D1_bz*(a + (1.-a)*1./(1.+pow((1.-h_b)/0.25,n)))*24*3600;
}

Form * HydratingDiffusionCementPaste::getCopy() const
{
    HydratingDiffusionCementPaste * ret = new HydratingDiffusionCementPaste() ;
    ret->doh = doh ;
    return ret ;
}

void HydratingDiffusionCementPaste::step(double timestep, ElementState & currentState, double maxscore)
{
    if(saturation.size() == 0)
        saturation.resize(6, 1.);
    else
        saturation = currentState.getDisplacements() ;
    double effectiveSaturation = (saturation[0]+saturation[1]+saturation[2])/3. ;
    effectiveSaturation = std::min(.999, effectiveSaturation) ;
    double deltaH = std::max(0.,getDeltadoh(effectiveSaturation, currentState)) ;
    doh = std::min(doh+deltaH, 1.) ;
    double D = std::abs(getDiffusionCoefficient(effectiveSaturation, currentState) ) ;
    param[0][0]=D ;
    param[1][1]=D ;
    currentState.getParent()->behaviourUpdated = true ;
}

double HydratingDiffusionCementPaste::getDeltadoh(double saturation, ElementState & currentState)
{
    double thresholdRelativeHumidity = 1.199*doh-0.1503 ;
    double maxRelativeHumidity = 0.999 ;
    double factor = (saturation-thresholdRelativeHumidity)/(maxRelativeHumidity-thresholdRelativeHumidity) ;

    if(factor < POINT_TOLERANCE)
        return 0. ;

// 	DelaunayTriangle * diffusionElement = dynamic_cast<DelaunayTriangle *>(currentState.getParent()) ;

    double deltaTime = currentState.getNodalDeltaTime() ;
// 	double currentTime = diffusionElement->getBoundingPoint(0).getT() + deltaTime*.5;

    double A1 =-0.347;
    double A2 = 278.8;
    double A3 =-0.677;
    double A4 = 1.37;

    //This is false. fix it!
    return ((-(A1/A2)*exp(-doh/A2) - (A3/A4)*exp(-doh/A4) )* deltaTime)*factor;

}

double HydratingDiffusionCementPaste::getDiffusionCoefficient(double saturation, ElementState & currentState)
{
    double D1_bz= (-5.412*doh*doh - 0.987*doh + 3.713)*1e-10;     //[m2/s], Duffison coefficient in function of doh

    double a = 0.04 ;
    int n = 6;
    double Abz=0.9002*doh-0.1503;    //Bolzmann function parameter in function of doh
    double Bbz=0.2815*doh+0.0297;
    double h_b=1.-1./(1.+exp((saturation-Abz)/Bbz));   //Bolzmann function for inversed desorption isotherm
    return D1_bz*(a + (1.-a)*1./(1.+pow((1.-h_b)/0.25,n)))*24*3600;
}


void HydratingDiffusionCementPaste::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine *vm) const
{

    ret[0][0] = (vm->ieval(VectorGradientDot(p_i) * param * VectorGradient(p_j, true),  gp, Jinv, v)
                 + vm->ieval(VectorGradient(p_i) * param * VectorGradientDot(p_j, true),  gp, Jinv, v)) ;
}

void HydratingDiffusionCementPaste::applyViscous(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine *vm) const
{
    ret[0][0] = (vm->ieval(Differential(p_i, TIME_VARIABLE) * Differential(p_j, TIME_VARIABLE),  gp, Jinv, v)
                 + vm->ieval(DoubleDifferential(p_j, TIME_VARIABLE, TIME_VARIABLE) * p_i,  gp, Jinv, v)) ;
}

HydratingDiffusionCementPaste::~HydratingDiffusionCementPaste() { }

/** \brief return true if the damage state has been modfied*/
bool HydratingDiffusionCementPaste::changed() const {
    return change ;
}

/** \brief Return the (damaged) Stifness tensor*/
Matrix HydratingDiffusionCementPaste::getTensor(const Point & p, IntegrableEntity * e , int g ) const {
    return param ;
}



