
//
// C++ Implementation: MCFT
//
// Description: implement the Modified compression field thery from Collins and Vecchio. This is a non-local, thermodynamics variant thereof
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "mcft.h"
#include "../damagemodels/damagemodel.h"
#include "../damagemodels/rotatingcrack.h"
#include "../../mesher/delaunay.h"
#include "../../mesher/delaunay_3d.h"
#include <set>

namespace Amie
{

NonLocalMCFT::NonLocalMCFT( double down, double youngModulus,  double charRad, RedistributionType r, MirrorState mirroring, double delta_x, double delta_y, double delta_z ) : FractureCriterion( mirroring, delta_x, delta_y, delta_z ),
    rtype(r), upVal( /*0.66*1e6*pow(std::abs(down*1e-6),.33)*/ .33e6*sqrt(-down*1e-6)*0.9),  downVal( down ), youngModulus(youngModulus)
{
    physicalCharacteristicRadius = charRad ;
    critStrain = -0.00163 ;//-0.0015;
    tensionCritStrain = upVal / youngModulus ;
    firstCompression = false ;
    secondCompression = false ;
    firstTension = false ;
    secondTension = false ;
    firstMet = false ;
    secondMet = false ;
    initialised = false ;
    radiusInitialised = false ;
    crackInitiated = false ;
    inRebarInfluence = false ;
    distanceToRebar = -1 ;
    effectiveInfluenceDistance = -1 ;
    smoothingType = QUARTIC_COMPACT ;
}


NonLocalMCFT::~NonLocalMCFT()
{
}

double NonLocalMCFT::getBareConcreteTensileCriterion(const ElementState & s, double pseudoYoung, double tstrain, double tstress)
{
    double maxTension = upVal;
    double maxTensionStrain = tensionCritStrain;

    if(tstrain > tensionCritStrain)
    {
        crackInitiated = true ;
        double downTestVal = tensionCritStrain ;
        double upTestVal = std::max(strain_te*100.,tstrain*100.) ;
        double factor = 1 ;
// 		double delta_tech = (strain_te-strain_ch);
        int count = 0 ;

        double sk = sqrt(k) ;
        while(std::abs(upTestVal-downTestVal) > 1e-4*tensionCritStrain && count++  <  64 )
        {
            double testVal = (upTestVal+downTestVal)*.5 ;
            double mainCurve = 1./(1.+sk*sqrt(testVal-tensionCritStrain)) ;

            if(testVal < tensionCritStrain)
                factor = 1. ;
            else //if(testVal < strain_ch)
                factor = mainCurve ;
// 			else if(testVal < strain_te)
// 				factor = mainCurve*(strain_te-testVal)/(delta_tech) ;
// 			else
// 				factor = 0 ;


            if( testVal*pseudoYoung > upVal*factor )
                upTestVal = testVal ;
            else
                downTestVal = testVal ;

        }

        maxTensionStrain = (upTestVal+downTestVal)*.5 ;
        maxTension = maxTensionStrain*pseudoYoung ;
        if(factor < POINT_TOLERANCE)
            return POINT_TOLERANCE ;

    }

    double criterion = 0 ;
    if(maxTension > POINT_TOLERANCE)
        criterion = std::abs(tstrain*pseudoYoung/maxTension) ;

    return -1.+criterion ;

    return POINT_TOLERANCE ;
}

double NonLocalMCFT::getRebarConcreteTensileCriterion(const Amie::ElementState& s, double pseudoYoung, double tstrain, double tstress, double value, double deltaCriterion)
{
    double maxTension = upVal;
//     double maxTensionAlt = upVal;
    

    if(tstrain < 0 || tstress < 0)
        return -1 ;
//     std::cout << tstrain << "  "<< tensionCritStrain << "  "<< std::endl ;
//     exit(0) ;
    if(tstrain > tensionCritStrain && tstress < upVal)
    {
        crackInitiated = true ;
        double factor = 1 ;
//         double factoralt = 1 ;

        double mainCurve = 1./(1.+sqrt(value*(tstrain-deltaCriterion))) ;
//         double mainCurveStress = 1./(1.+sqrt(value*(tstress/pseudoYoung-deltaCriterion))) ;

        if(tstrain < strain_ch*(1./(deltaCriterion/tensionCritStrain)))
        {
            factor = mainCurve ;
//             factoralt = mainCurveStress ;
        }
        else
        {
            factor = mainCurve*exp(-1.5*(tstrain-strain_ch)/(strain_te-strain_ch)) ;
//             factoralt = mainCurveStress*exp(-1.5*(tstress/pseudoYoung-strain_ch)/(strain_te-strain_ch)) ;
        }
        
        maxTension *= factor ;
//         maxTensionAlt *= factoralt ;

    }
    else
        return -1.+std::max(std::abs(tstrain/tensionCritStrain), std::abs(tstress/upVal)) ;

    double criterion0 = 0 ;
    
    if(maxTension > POINT_TOLERANCE)
        return  std::max(std::abs(tstress/maxTension), std::abs(tstrain*pseudoYoung/maxTension))-1. ;
    else
        return tstress/upVal-1. ;


    return -1.+criterion0 ; 

}

double NonLocalMCFT::getConcreteTensileCriterion(const ElementState & s, double pseudoYoung, double tstrain, double tstress)
{

// 	if(tstress <= 0 || tstrain <= 0)
// 		return -1 ;
// 	return getRebarConcreteTensileCriterion(s, pseudoYoung, tstrain, tstress) ;
// 	return getBareConcreteTensileCriterion(s, pseudoYoung, tstrain, tstress) ;

// 	if(dynamic_cast<RotatingCrack *>(s.getParent()->getBehaviour()->getDamageModel()))
// 	{
// 		if(dynamic_cast<RotatingCrack *>(s.getParent()->getBehaviour()->getDamageModel())->getT()ensionFailure)
// 			return POINT_TOLERANCE ;
// 	}

    if(!rebarLocationsAndDiameters.empty() && !inRebarInfluence && distanceToRebar < 0 && effectiveInfluenceDistance < 0)
    {
        distanceToRebar = std::abs(s.getParent()->getCenter().getY() - rebarLocationsAndDiameters[0].first) ;
        effectiveInfluenceDistance = rebarLocationsAndDiameters[0].second*7.5*2;
        inRebarInfluence = distanceToRebar < rebarLocationsAndDiameters[0].second*7.5*2 ;
// 		std::cout << rebarLocationsAndDiameters[0].first << "  "<<std::flush ;
        for(size_t i = 1 ; i < rebarLocationsAndDiameters.size() ; i++)
        {
            double distanceToRebartest = std::abs(s.getParent()->getCenter().getY() - rebarLocationsAndDiameters[i].first) ;
// 			std::cout << rebarLocationsAndDiameters[i].first << "  "<<std::flush ;
            if( distanceToRebar < rebarLocationsAndDiameters[i].second*7.5*2 && distanceToRebartest < distanceToRebar)
            {
                distanceToRebar = distanceToRebartest ;
                effectiveInfluenceDistance = rebarLocationsAndDiameters[i].second*7.5*2 ;
                inRebarInfluence = true ;
            }
        }
// 		std::cout << std::endl ;
    }

    if(!inRebarInfluence)
        return getRebarConcreteTensileCriterion(s, pseudoYoung, tstrain, tstress, k, tensionCritStrain) ;



    double f = (distanceToRebar)/(effectiveInfluenceDistance) ;
// 	double df = f ;

    double df = 3.*f*f-2.*f*f*f ;
    double factor = df*log(k)+(1.-df)*log(60.*342.) ;//40 low
    double deltaCriterion = tensionCritStrain*df ;
    if(df <= 0)
    {
        factor = log(k) ;
        deltaCriterion = tensionCritStrain ;
    }

    if(df >= 1)
    {
        factor = log(60.*342); // 40 low
        deltaCriterion = 0 ;
    }
    return getRebarConcreteTensileCriterion(s, pseudoYoung, tstrain, tstress, exp(factor), deltaCriterion) ;

// 	return df*barecrit + (1.-df)*rebcrit ;
}

double NonLocalMCFT::getConcreteCompressiveCriterion(const ElementState & s, double pseudoYoung, double cstrain, double tstress, double cstress)
{

    if(cstress > 0 || cstrain > 0)
        return -1 ;
    
    double maxCompression = downVal  ;


    if(cstrain < .1*critStrain)
    {
        double C_d = 0. ;
        double compressiveTensileRatio = -std::abs(tstress/std::min(cstress, -POINT_TOLERANCE)) ;

        if(-compressiveTensileRatio > 0.28000)
            C_d =0.35*pow(-compressiveTensileRatio-0.28, 0.8) ;

        double beta_d = 1./(1.+C_d) ;

        if(beta_d > 1)
            beta_d = 1 ;

        double f_p = beta_d*downVal ;
        double epsilon_p = beta_d*critStrain ;
        double epsratio = cstrain/epsilon_p  ;
        double n = 0.8 - f_p/17e6 ;
        double k_c = 0.67 - f_p/62e6 ;
        if(epsratio >= 1)
        {
            k_c = 1. ;
        }
        else
        {
            k_c = 0.67 - f_p/62e6 ;
        }
        maxCompression = f_p*n*(epsratio)/(n-1.+pow(epsratio,n*k_c)) ;
    }


    double criterion = std::max(std::abs(cstress/downVal), std::abs(cstrain*pseudoYoung/downVal)) ;

    if(std::abs(maxCompression) > POINT_TOLERANCE)
        criterion =  std::max(std::abs(cstress/maxCompression), std::abs(cstrain*pseudoYoung/maxCompression)) ;


    return -1.+criterion ;

}

double sqrtdecrease(double k0, double upVal, double eps_0, double strain_ch, double strain_te, double eps)
{
    if(eps < strain_ch)
        return upVal/(1.+sqrt(k0*(eps-eps_0))) ;
    else if(eps < strain_te)
        return (upVal/(1.+sqrt(k0*(eps-eps_0))))*((strain_te-eps)/(strain_te-strain_ch)) ;

    return 0 ;
}

void NonLocalMCFT::initialise( ElementState &s)
{
    double energy = /*physicalCharacteristicRadius / sqrt(s.getParent()->area())*/75. ; //N/m 32000 prev
    strain_ch = 2.*energy/(upVal) ;//*.5 energy <- // *2 energy -> 2.*energy/(1.*getMaterialCharacteristicRadius()*upVal) ;

// 	energy *= nlcorrection ;
    if(strain_ch < tensionCritStrain)
    {
        std::cout << strain_ch << " vs " << tensionCritStrain <<std::endl ;
// 		exit(0) ;
        strain_ch = 4.*tensionCritStrain ;
    }
    strain_te = 5.*strain_ch;

    double k0 = 1 ;
    double integral(0) ;

    double de = 1e-2 ;
    do
    {
        integral = 0 ;
        k0 *= 10;

        for(double i = tensionCritStrain ; i < strain_te ; i+= de*(strain_te-tensionCritStrain))
        {
            integral+= 0.5*(sqrtdecrease(k0, upVal, tensionCritStrain, strain_ch, strain_te, i)+sqrtdecrease(k0, upVal, tensionCritStrain, strain_ch, strain_te, i+de*(strain_te-tensionCritStrain)))*de*(strain_te-tensionCritStrain) ;
        }
    }
    while(integral > energy) ;


    double k_low = k0/100.  ;
    double k_high = k0 ;

    int count = 0 ;

    do
    {
        integral = 0 ;
        k = 0.5*(k_low+k_high) ;

        for(double i = tensionCritStrain ; i < strain_te ; i+= de*(strain_te-tensionCritStrain))
        {
            integral+= 0.5*(sqrtdecrease(k, upVal, tensionCritStrain, strain_ch, strain_te, i)+sqrtdecrease(k, upVal, tensionCritStrain, strain_ch, strain_te, i+de*(strain_te-tensionCritStrain)))*de*(strain_te-tensionCritStrain) ;
        }

        if(integral < energy)
            k_high = k ;
        else
            k_low = k ;
        count++ ;

    } while(std::abs(integral-energy) > 1e-3*energy && count < 64) ;



// 	k*= getMaterialCharacteristicRadius()*getMaterialCharacteristicRadius()*M_PI ;
// 		k/=1e5 ;

    if(std::abs(integral-energy) > 1e-3*energy)
    {
        std::cout << "wrong energy! " << integral << " vs "<< energy  << " k = "<< k << std::endl ;
        exit(0) ;
    }
// 	k = 1000 ;
    initialised = true ;
}


// rough crack model by bazant and gambarova
// From engineering fracture mechanics, smeared crack analysis for fracture and aggregate
// interlock in concrete, Gambarova and Valente
// std::pair<double, double> stressOnCrack(ElementState &s, double downVal, double dmax = 16)
// {
// 	double tau_0 = 0.25*downVal ;
// 	double a12 = 0.62*downVal ;
// 	double a3 = 2.45/tau_0 ;
// 	double a4 = 2.44*(1.-4./tau_0) ;
//
// 	std::pair<double, double> deltas = s.getParent()->getBehaviour()->getFractureCriterion()->getCrackOpeningAndSlip(s) ;
// 	double delta_n = deltas.first*1000. ;
// 	double delta_t = deltas.second*1000. ;
// 	if(delta_n < 0.01)
// 		return std::make_pair(0, 0) ;
//
// 	double r = delta_t/delta_n ;
// 	if(std::abs(delta_n) < POINT_TOLERANCE)
// 		r = 0 ;
// 	double f = a3+a4*std::abs(r*r*r) ;
// 	double g = 1.+a4*r*r*r*r ;
// 	double h = pow(1.+r*r, 0.25) ;
//
// 	double sigma_nt = tau_0*(1.-sqrt(2.*delta_n/dmax))*r*f/g ;
// 	double sigma_nn = a12*r*sqrt(delta_n)*sigma_nt*.25 ;
//
// 	return std::make_pair(sigma_nn, sigma_nt) ;
// }

double NonLocalMCFT::gradeAtTime(ElementState &s, double t)
{
    if(!initialised)
        initialise(s);
    
    if(s.getParent()->getBehaviour()->getDamageModel()->fractured())
        return -1 ;
        
    std::pair<Vector, Vector> sstrain = getSmoothedFields(PRINCIPAL_MECHANICAL_STRAIN_FIELD,PRINCIPAL_REAL_STRESS_FIELD , s,  t) ;
    Vector first = sstrain.second ;
    Vector second = sstrain.first ;

    double tstrain = second.max();
    double cstrain = second.min();

    double tstress = first.max();
    double cstress = first.min();


    double pseudoYoung = youngModulus*(1.-s.getParent()->getBehaviour()->getDamageModel()->getState().max());
    firstCompression = tstrain < 0 && tstress < 0;
    secondCompression = cstrain < 0 && cstress < 0;
    firstTension = tstrain > 0 && tstress > 0;
    secondTension = cstrain > 0 && cstress > 0;
    firstMet = false ;
    secondMet = false ;
// 	if(firstTension && secondTension)
// 	{
// 		tstrain -= 0.5*(tstrain-(tstrain-cstrain)) ;
// 		tstress -= 0.5*(tstress-(tstress-cstress)) ;
// 	}


    if(s.getParent()->getBehaviour()->getDamageModel()->getState().size() == 4)
    {
        double tpseudoYoung0 =  youngModulus*(1.-s.getParent()->getBehaviour()->getDamageModel()->getState()[0]) ;
        double cpseudoYoung0 =  youngModulus*(1.-s.getParent()->getBehaviour()->getDamageModel()->getState()[1]) ;
        double tpseudoYoung1 =  youngModulus*(1.-s.getParent()->getBehaviour()->getDamageModel()->getState()[2]) ;
        double cpseudoYoung1 =  youngModulus*(1.-s.getParent()->getBehaviour()->getDamageModel()->getState()[3]) ;

        double ccrit0 = -1 ;
        double tcrit0 = -1 ;
        double ccrit1 = -1 ;
        double tcrit1 = -1 ;
        firstTension = tstress >= 0 ;
        firstCompression = tstress < 0 ;
        secondTension = cstress >= 0 ;
        secondCompression = cstress < 0 ;


        if(cpseudoYoung0 > POINT_TOLERANCE*youngModulus )
        {

            ccrit0 = getConcreteCompressiveCriterion(s, cpseudoYoung0, cstrain, tstress, cstress) ;

            if(ccrit0 > 0)
                firstMet = true ;
        }

        if(tpseudoYoung0 > POINT_TOLERANCE*youngModulus )
        {
            tcrit0 = getConcreteTensileCriterion(s, tpseudoYoung0, tstrain, tstress) ;

            if(tcrit0 > 0 && tcrit0 > ccrit0)
                firstMet = true ;
        }
        double c0 = std::max(ccrit0, tcrit0) ;

        if(ccrit0 > tcrit0 && firstMet)
        {
            firstCompression = true ;
            firstTension = false ;
        }
        else
        {
            firstCompression = false ;
            firstTension = true ;
        }

        if(cpseudoYoung1 > POINT_TOLERANCE*youngModulus )
        {
            ccrit1 = getConcreteCompressiveCriterion(s, cpseudoYoung1, cstrain, tstress, cstress) ;

            if(ccrit1 > 0 && ccrit1 > ccrit0 && ccrit1 > tcrit0)
            {
                secondMet = true ;
                firstMet = false ;
            }
        }

        if(tpseudoYoung1 > POINT_TOLERANCE*youngModulus )
        {
            tcrit1 = getConcreteTensileCriterion(s, tpseudoYoung1, tstrain, tstress) ;
            if(tcrit1 > 0 && tcrit1 > ccrit0 && tcrit1 > tcrit0 && tcrit1 > ccrit1)
            {
                secondMet = true ;
                firstMet = false ;
            }
        }
        double c1 = std::max(ccrit1, tcrit1) ;
        if(ccrit1 > tcrit1 && secondMet)
        {
            secondCompression = true ;
            secondTension = false ;
        }
        else
        {
            secondCompression = false ;
            secondTension = true ;
        }

        return std::max(c0, c1) ;
    }
    else if (s.getParent()->getBehaviour()->getDamageModel()->getState().size() == 2)
    {
        double cpseudoYoung =  youngModulus*(1.-s.getParent()->getBehaviour()->getDamageModel()->getState()[1]) ;
        double tpseudoYoung =  youngModulus*(1.-s.getParent()->getBehaviour()->getDamageModel()->getState()[0]) ;
        double ccrit = -1 ;
        double tcrit = -1 ;

        if(firstTension && tpseudoYoung > POINT_TOLERANCE)
        {
            tcrit = getConcreteTensileCriterion(s, tpseudoYoung, tstrain, tstress) ;
            firstMet = true ;
        }

        if(secondCompression && cpseudoYoung > POINT_TOLERANCE)
        {
            ccrit = getConcreteCompressiveCriterion(s, cpseudoYoung, cstrain, tstress, cstress) ;
            secondMet = true ;
        }

        return std::max(ccrit,tcrit) ;
    }

        
    if(pseudoYoung < 1e-9)
    {
        firstMet = true ;
        secondMet = true ;
        return -1 ;
    }

    double ccrit = -1 ;
    double tcrit = -1 ;

    if(firstTension && pseudoYoung > POINT_TOLERANCE)
    {
        tcrit = getConcreteTensileCriterion(s, pseudoYoung, tstrain, tstress) ;
        firstMet = true ;
    }
    else if(secondCompression && pseudoYoung > POINT_TOLERANCE)
    {
        ccrit = getConcreteCompressiveCriterion(s, pseudoYoung, cstrain,tstress, cstress) ;
        secondMet = true ;
    }

    return std::max(tcrit,ccrit) ;

}

double NonLocalMCFT::grade( ElementState &s )
{
    return gradeAtTime(s, 0) ;
}

double NonLocalMCFT::getTensileLimit(const ElementState & s) const
{
    double pseudoYoung = youngModulus*(1.-std::min(s.getParent()->getBehaviour()->getDamageModel()->getState().max(), 1.-1e-12));
    std::vector<double> crits ;


    double downTestVal = 0 ;
    double upTestVal = upVal ;
    double factor = 1 ;
    double delta_tech = strain_te-strain_ch;
    while(std::abs(upTestVal-downTestVal) > 1e-4*upVal)
    {

        double testVal = (upTestVal+downTestVal)*.5/pseudoYoung ;
        double mainCurve = 1./(1.+k*sqrt((testVal-tensionCritStrain))) ;

        if(testVal < tensionCritStrain)
            factor = 1. ;
        else if(testVal < strain_ch)
            factor = mainCurve ;
        else if(testVal < strain_te)
            factor = mainCurve*(strain_te-testVal)/delta_tech ;
        else
            factor = 0 ;

        if( testVal*pseudoYoung > upVal*factor )
        {
            upTestVal = testVal*pseudoYoung ;
        }
        else
        {
            downTestVal = testVal*pseudoYoung ;
        }

    }

    return (upTestVal+downTestVal)*.5 ;

}

FractureCriterion *NonLocalMCFT::getCopy() const
{
    return new NonLocalMCFT( *this ) ;
}

NonLocalSpaceTimeMCFT::NonLocalSpaceTimeMCFT(double down, double youngModulus, double charDistance, RedistributionType r, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : NonLocalMCFT(down, youngModulus, charDistance, r , mirroring, delta_x, delta_y, delta_z) { }

NonLocalSpaceTimeMCFT::~NonLocalSpaceTimeMCFT() { }

double NonLocalSpaceTimeMCFT::grade(ElementState &s)
{    
    double gradeBefore = gradeAtTime(s, -1) ;
    double gradeAfter = gradeAtTime(s, 1) ;

    if(gradeAfter < 0)
        return gradeAfter ;
    if(gradeBefore > 0)
    {
        return .99 ;
    }

    double upTime = 1 ;
    double downTime = -1 ;
    double testTime = 0 ;
    
    while(std::abs(upTime-downTime) > 1e-7)
    {
        double gradeTest = gradeAtTime(s, testTime) ;
        if(gradeTest < 0)
            downTime = testTime ;
        else if(gradeTest > 0)
            upTime = testTime ;
        else
            return testTime ;
        
        testTime = 0.5*(downTime+upTime) ;
    }
    return 1.-(testTime*.5+.5) ;
}

FractureCriterion *NonLocalSpaceTimeMCFT::getCopy() const
{
    return new NonLocalSpaceTimeMCFT( *this ) ;
}

}
