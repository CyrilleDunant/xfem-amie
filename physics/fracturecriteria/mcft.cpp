
//
// C++ Implementation: mohrcoulomb
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007-2011
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

namespace Mu
{


NonLocalMCFT::NonLocalMCFT( double down, double youngModulus,  double charRad, RedistributionType r, MirrorState mirroring, double delta_x, double delta_y, double delta_z ) : FractureCriterion( mirroring, delta_x, delta_y, delta_z )
	, upVal( /*0.66*1e6*pow(std::abs(down*1e-6),.33)*/.25*1e6*sqrt(-down*1e-6) ), downVal( down ), youngModulus(youngModulus), rtype(r)
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
		double upTestVal = strain_te*1.1 ;
		double factor = 1 ;
		double delta_tech = (strain_te-strain_ch);
		int count = 0 ;
	
		while(std::abs(upTestVal-downTestVal) > 1e-8*tensionCritStrain && count++  <  32 )
		{
			double testVal = (upTestVal+downTestVal)*.5 ;
			double mainCurve = 1./(1.+sqrt(k)*sqrt(testVal-tensionCritStrain)) ;

			if(testVal < tensionCritStrain)
				factor = 1. ;
			else if(testVal < strain_ch)
				factor = mainCurve ;
			else if(testVal < strain_te)
				factor = mainCurve*(strain_te-testVal)/(delta_tech) ;
			else
				factor = 0 ;
			
			
			if( testVal*pseudoYoung > upVal*factor )
				upTestVal = testVal ;
			else
				downTestVal = testVal ;
			
		}

		maxTensionStrain = (upTestVal+downTestVal)*.5 ;
		maxTension = maxTensionStrain*pseudoYoung ;
		if(factor < POINT_TOLERANCE_2D)
			return POINT_TOLERANCE_2D ;

	}
	
	double criterion = 0 ;
  if(maxTensionStrain > POINT_TOLERANCE_2D)
		criterion = std::abs(tstress/(maxTension)) ;

	double strainc = -1 ;
		
	
// 	if(criterion >= 1)
// 	{
// 		strainc = 1.-1./criterion ;
// 	}
// 	else if(criterion >= 0)
// 	{
		strainc = -1.+criterion ;
// 	}
// 	else if(criterion > -1)
// 	{
// 		strainc = criterion ;
// 	}
	
	criterion = 0 ;
  if(maxTensionStrain > POINT_TOLERANCE_2D)
		criterion = std::abs(tstrain/maxTensionStrain) ;
	
	double stressc = -1 ;
		
	
	if(criterion >= 1)
	{
		stressc = 1.-1./criterion ;
	}
	else if(criterion >= 0)
	{
		stressc = -1.+criterion ;
	}
	else if(criterion > -1)
	{
		stressc = criterion ;
	}
	return strainc ; //stressc; // strain criterion misnamed
	
	return POINT_TOLERANCE_2D ;
}

double NonLocalMCFT::getRebarConcreteTensileCriterion(const Mu::ElementState& s, double pseudoYoung, double tstrain, double tstress, double value)
{

	
 	double maxTension = upVal;
 	double maxTensionStrain = tensionCritStrain;

 	if(tstrain > tensionCritStrain )
 	{
		crackInitiated = true ;
 		double downTestVal = tensionCritStrain ;
 		double upTestVal = strain_te*2. ;
 		double factor = 1 ;
 		double delta_tech = (strain_te-strain_ch)*2.;
 		double mainCurve = 0 ;
 		int count = 0 ;

		while(std::abs(upTestVal-downTestVal) > 1e-8*tensionCritStrain && count++  <  32)
		{
			double testVal = (upTestVal+downTestVal)*.5 ;
			mainCurve = 1./(1.+sqrt(value*(testVal))) ;

// 			if(testVal < strain_ch*2.)
				factor = mainCurve ;
// 			else if(testVal < strain_te*2.)
// 			{
// 				factor = mainCurve*(strain_te*2.-testVal)/(delta_tech) ;
// 			}
// 			else
// 				factor = 0 ;
			

			if( testVal*pseudoYoung > upVal*factor )
				upTestVal = testVal ;
			else
				downTestVal = testVal ;

		}

 	
 		maxTension = (upTestVal+downTestVal)*.5*pseudoYoung ;
 		maxTensionStrain = (upTestVal+downTestVal)*.5 ;
 
 		if(factor < POINT_TOLERANCE_2D)
 			return POINT_TOLERANCE_2D ;
 	}

	double criterion = 0 ;
	if(maxTensionStrain > POINT_TOLERANCE_2D)
		criterion = std::abs(tstrain/maxTensionStrain) ;

	
	double strainc = -1 ;
		
	
	if(criterion >= 1)
	{
		strainc = 1.-1./criterion ;
	}
	else if(criterion >= 0)
	{
		strainc = -1.+criterion ;
	}
	else if(criterion > -1)
	{
		strainc = criterion ;
	}
	
	
	criterion = 0 ;
	if(maxTension > POINT_TOLERANCE_2D)
		criterion = std::abs(tstress/maxTension) ;
	
	double stressc = -1 ;
		
	
// 	if(criterion >= 1)
// 	{
// 		stressc = 1.-1./criterion ;
// 	}
// 	else if(criterion >= 0)
// 	{
		stressc = -1.+criterion ;
// 	}
// 	else if(criterion > -1)
// 	{
// 		stressc = criterion ;
// 	}
	return stressc ; //strainc ;
	
	return POINT_TOLERANCE_2D ;
}

double NonLocalMCFT::getConcreteTensileCriterion(const ElementState & s, double pseudoYoung, double tstrain, double tstress)
{
	
	if(tstress <= 0 || tstrain <= 0)
		return -1 ;
// 	return getRebarConcreteTensileCriterion(s, pseudoYoung, tstrain, tstress) ;
// 	return getBareConcreteTensileCriterion(s, pseudoYoung, tstrain, tstress) ;
	
// 	if(dynamic_cast<RotatingCrack *>(s.getParent()->getBehaviour()->getDamageModel()))
// 	{
// 		if(dynamic_cast<RotatingCrack *>(s.getParent()->getBehaviour()->getDamageModel())->tensionFailure)
// 			return POINT_TOLERANCE_2D ;
// 	}

	if(!rebarLocationsAndDiameters.empty() && !inRebarInfluence && distanceToRebar < 0 && effectiveInfluenceDistance < 0)
	{
		distanceToRebar = std::abs(s.getParent()->getCenter().y - rebarLocationsAndDiameters[0].first) ;
		effectiveInfluenceDistance =  rebarLocationsAndDiameters[0].second*7.5*4;
		inRebarInfluence = distanceToRebar < rebarLocationsAndDiameters[0].second*7.5*4 ;
// 		std::cout << rebarLocationsAndDiameters[0].first << "  "<<std::flush ;
		for(size_t i = 1 ; i < rebarLocationsAndDiameters.size() ; i++)
		{
			double distanceToRebartest = std::abs(s.getParent()->getCenter().y - rebarLocationsAndDiameters[i].first) ;
// 			std::cout << rebarLocationsAndDiameters[i].first << "  "<<std::flush ;
			if( distanceToRebar < rebarLocationsAndDiameters[i].second*7.5*4 && distanceToRebartest < distanceToRebar)
			{
				distanceToRebar = distanceToRebartest ;
				effectiveInfluenceDistance = rebarLocationsAndDiameters[i].second*7.5*4 ;
				inRebarInfluence = true ;
			}
		}
// 		std::cout << std::endl ;
	}
// 	if(inRebarInfluence)
// 	  return -1 ;
/*	  std::cout << s.getParent()->getCenter().y << std::endl */;
	double barecrit = getBareConcreteTensileCriterion(s, pseudoYoung, tstrain, tstress) ;
// 	
	if(!inRebarInfluence)
		return barecrit ;

 	double rebcrit = getRebarConcreteTensileCriterion(s, pseudoYoung, tstrain, tstress) ;
	
// 	if(distanceToRebar <= 0.5*effectiveInfluenceDistance)
// 		return rebcrit ;

// 	if (!radiusInitialised)
// 	{
// 		double f = (distanceToRebar)/(effectiveInfluenceDistance) ;
// 		double df = 3.*f*f-2.*f*f*f ;
// 		setMaterialCharacteristicRadius(getMaterialCharacteristicRadius()*df+.5*getMaterialCharacteristicRadius()*(1.-df));
// 		radiusInitialised = true ;
// 	}

	double f = (distanceToRebar)/(effectiveInfluenceDistance) ;
// 	double df = f ;
	double df = 3.*f*f-2.*f*f*f ;
//
// 	return barecrit ;
	return df*barecrit + (1.-df)*rebcrit ;
}

double NonLocalMCFT::getConcreteCompressiveCriterion(const ElementState & s, double pseudoYoung, double cstrain, double tstress, double cstress)
{
	if(cstress >= 0 )
		return -1 ;
// 	if(dynamic_cast<RotatingCrack *>(s.getParent()->getBehaviour()->getDamageModel()))
// 	{
// 		if(dynamic_cast<RotatingCrack *>(s.getParent()->getBehaviour()->getDamageModel())->compressionFailure)
// 			return POINT_TOLERANCE_2D ;
// 	}
	double maxCompression = downVal  ;
	

	
	double maxCompressionStrain = downVal/pseudoYoung  ;
	
	if( cstrain < 0.1*critStrain )
	{
		double C_d = 0. ;
		double compressiveTensileRatio = -std::abs(tstress/std::min(cstress, -POINT_TOLERANCE_2D)) ;
		
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
		maxCompression = std::max(f_p*n*(epsratio)/(n-1.+pow(epsratio,n*k_c)), downVal) ;
		maxCompressionStrain = maxCompression/pseudoYoung ;

// 		if(n*f_p/(epsilon_p*pseudoYoung)-n+1. > 0)
// 		{
// 			
// 			double rcs = pow(n*f_p/(epsilon_p*pseudoYoung)-n+1., 1./(n*k_c)) ;
// 
// 			if(rcs*pseudoYoung*pseudoYoung > 1)
// 				k_c = 1. ;
// 			else
// 				k_c = 0.67 - f_p/62e6 ;
// 			
// 			rcs = pow(n*f_p/(epsilon_p*pseudoYoung)-n+1., 1./(n*k_c)) ;
// 
// 			maxCompressionStrain = rcs*epsilon_p ;
// 		}
// 		else
// 		{
// 			std::cout << "plouf " << n*f_p/(epsilon_p*pseudoYoung)-n+1. << std::endl ;
// // 			exit(0) ;
// 		}
	}
	
	double criterion = 0 ;

	if(std::abs(maxCompressionStrain) > POINT_TOLERANCE_2D)
		criterion = std::abs(cstrain/maxCompressionStrain) ;
	
// 	if(cstrain < 0.05*critStrain)
// 	{
// 		std::cout << "\n" << criterion << "\n" <<std::endl ;
// 		exit(0) ;
// 	}
	
	double strainc = -1 ;
		
	
	if(criterion >= 1)
	{
		strainc = 1.-1./criterion ;
	}
	else if(criterion >= 0)
	{
		strainc = -1.+criterion ;
	}
	else if(criterion > -1)
	{
		strainc = criterion ;
	}
	
	if(std::abs(maxCompression) > POINT_TOLERANCE_2D)
		criterion = std::abs(cstress/maxCompression) ;
	
	double stressc = -1 ;
		
	
// 	if(criterion >= 1)
// 	{
// 		stressc = 1.-1./criterion ;
// 	}
// 	else if(criterion >= 0)
// 	{
		stressc = -1.+criterion ;
// 	}
// 	else if(criterion > -1)
// 	{
// 		stressc = criterion ;
// 	}
	return stressc ;
	return POINT_TOLERANCE_2D ;
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
	double energy = 75. ; //N/m 32000 prev
	strain_ch = 2.*energy/(upVal) ;//*.5 energy <- // *2 energy -> 2.*energy/(1.*getMaterialCharacteristicRadius()*upVal) ;
	if(factors.size()==0)
		initialiseFactors(s) ;
// 	energy *= nlcorrection ;
	if(strain_ch < tensionCritStrain)
	{
		std::cout << strain_ch << " vs " << tensionCritStrain <<std::endl ;
// 		exit(0) ;
		strain_ch = tensionCritStrain ;
	}
	strain_te = 5.*strain_ch;
	
	double k0 = 1 ;
	double integral(0) ;
	
	double de = 1e-5 ;
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
		
	} while(std::abs(integral-energy) > 1e-6*energy && count < 64) ;
	
	
	
// 	k*= getMaterialCharacteristicRadius()*getMaterialCharacteristicRadius()*M_PI ;
// 		k/=1e5 ;
		
// 	if(std::abs(integral-energy) > 1e-4*energy)
// 	{
// 	  std::cout << "wrong energy! " << integral << " vs "<< energy  << " k = "<< k << std::endl ;
// 	  exit(0) ;
// 	}
// 	k = 1000 ;
	initialised = true ;
}


// rough crack model by bazant and gambarova
// From engineering fracture mechanics, smeared crack analysis for fracture and aggregate 
// interlock in concrete, Gambarova and Valente
std::pair<double, double> stressOnCrack(ElementState &s, double downVal, double dmax = 16)
{
	double tau_0 = 0.25*downVal ;
	double a12 = 0.62*downVal ;
	double a3 = 2.45/tau_0 ;
	double a4 = 2.44*(1.-4./tau_0) ;
	
	std::pair<double, double> deltas = s.getParent()->getBehaviour()->getFractureCriterion()->getCrackOpeningAndSlip(s) ;
	double delta_n = deltas.first*1000. ;
	double delta_t = deltas.second*1000. ;
	if(delta_n < 0.01)
		return std::make_pair(0, 0) ;
	
	double r = delta_t/delta_n ;
	if(std::abs(delta_n) < POINT_TOLERANCE_2D)
		r = 0 ;
	double f = a3+a4*std::abs(r*r*r) ;
	double g = 1.+a4*r*r*r*r ;
	double h = pow(1.+r*r, 0.25) ;
	
	double sigma_nt = tau_0*(1.-sqrt(2.*delta_n/dmax))*r*f/g ;
	double sigma_nn = a12*r*sqrt(delta_n)*sigma_nt*.25 ;
	
	return std::make_pair(sigma_nn, sigma_nt) ;
}

double NonLocalMCFT::grade( ElementState &s )
{
	if(!initialised)
		initialise(s);
// 	std::pair<Vector, Vector> stressStrain = smoothedPrincipalStressAndStrain(s, REAL_STRESS) ;
	std::pair<Vector, Vector> sstrain = smoothedStressAndStrain(s, REAL_STRESS) ;
	Vector first = toPrincipal(sstrain.first) ;
	Vector second = toPrincipal(sstrain.second) ;
	
	std::pair<double, double> crackstress = std::make_pair(0.,0.); 
	
	if(s.getParent()->getBehaviour()->getDamageModel()->getState().max() > .1)
		crackstress = stressOnCrack(s, downVal) ;
		
	double tstrain = second.max();//stressStrain.second.max() ; //
	double cstrain = second.min();//stressStrain.second.min();  //
	double sfactor = factors[0]/std::accumulate(&factors[0], &factors[factors.size()],double(0)) ;
// 	if(std::max(crackstress.first, 0.) > POINT_TOLERANCE_2D ||  std::min(crackstress.second, 0.) > POINT_TOLERANCE_2D)
// 		std::cout << "\n"<< std::max(crackstress.first, 0.) << "   " << std::min(crackstress.second, 0.) << std::endl ;
	double tstress = first.max(); //-std::max(std::abs(crackstress.first), 0.);   //first.max();//std::min(first.max()-2.*std::max(crackstress.first, 0.), first.max());// 
	double cstress = first.min();//+std::min(std::abs(crackstress.second), 0.);   //first.min();//std::max(first.min()+2.*std::min(crackstress.second, 0.), first.min());// 

	
	
	
	double pseudoYoung = youngModulus*(1.-s.getParent()->getBehaviour()->getDamageModel()->getState().max());
	firstCompression = tstrain < 0 ;
	secondCompression = cstrain < 0 ;
	firstTension = tstrain >=0 ;
	secondTension = cstrain >= 0 ;
	firstMet = false ;
	secondMet = false ;

	
	if(s.getParent()->getBehaviour()->getDamageModel()->getState().size() == 4)
	{
		double cpseudoYoung0 =  youngModulus*(1.-s.getParent()->getBehaviour()->getDamageModel()->getState()[1]) ;
		double tpseudoYoung0 =  youngModulus*(1.-s.getParent()->getBehaviour()->getDamageModel()->getState()[0]) ;
		double cpseudoYoung1 =  youngModulus*(1.-s.getParent()->getBehaviour()->getDamageModel()->getState()[3]) ;
		double tpseudoYoung1 =  youngModulus*(1.-s.getParent()->getBehaviour()->getDamageModel()->getState()[2]) ;
		double ccrit0 = -1 ;
		double tcrit0 = -1 ;
		double ccrit1 = -1 ;
		double tcrit1 = -1 ;
		if(firstCompression && cpseudoYoung0 > POINT_TOLERANCE_2D && tstress < 0)
		{
			ccrit0 = getConcreteCompressiveCriterion(s, cpseudoYoung0, cstrain, tstress, cstress) ;
			if(ccrit0 > 0)
				firstMet = true ;
		}
		else if(firstTension && tpseudoYoung0 > POINT_TOLERANCE_2D && tstress > 0)
		{
			tcrit0 = getConcreteTensileCriterion(s, tpseudoYoung0, tstrain, tstress) ;
			if(tcrit0 > 0)
				firstMet = true ;
		}
		double c0 = std::max(ccrit0, tcrit0) ;
		if(secondCompression && cpseudoYoung1 > POINT_TOLERANCE_2D && cstress < 0)
		{
			ccrit1 = getConcreteCompressiveCriterion(s, cpseudoYoung1, cstrain, tstress, cstress) ;
			if(ccrit1 > 0)
			{
				secondMet = true ;
			}
		}
		else if(secondTension && tpseudoYoung1 > POINT_TOLERANCE_2D && cstress > 0)
		{
			tcrit1 = getConcreteTensileCriterion(s, tpseudoYoung1, tstrain, tstress) ;
			if(tcrit1 > 0)
			{
				secondMet = true ;
			}
		}
		double c1 = std::max(ccrit1, tcrit1) ;
		if(c0 > c1)
		{
			firstMet = true ;
			secondMet = false ;
		}
		else
		{
			firstMet = false ;
			secondMet = true ;
		}
// 		std::cout << cstress << "  " << tstress << "  "<< ccrit0 << "  "<< tcrit0 << "  "<< tcrit1 << "  "<< ccrit1 << std::endl ;
		return std::max(c0, c1) ;
	}
	else if (s.getParent()->getBehaviour()->getDamageModel()->getState().size() == 2)
	{
		double cpseudoYoung =  youngModulus*(1.-s.getParent()->getBehaviour()->getDamageModel()->getState()[1]) ;
		double tpseudoYoung =  youngModulus*(1.-s.getParent()->getBehaviour()->getDamageModel()->getState()[0]) ;
		double ccrit = -1 ;
		double tcrit = -1 ;

		if(firstTension && tpseudoYoung > POINT_TOLERANCE_2D)
		{
			tcrit = getConcreteTensileCriterion(s, tpseudoYoung, tstrain, tstress) ;
			firstMet = true ;
		}

		if(secondCompression && cpseudoYoung > POINT_TOLERANCE_2D)
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
	
	double ccrit = getConcreteCompressiveCriterion(s, pseudoYoung, cstrain, tstress, cstress) ;
	double tcrit = getConcreteTensileCriterion(s, pseudoYoung, tstrain, tstress) ;
	
	if( ccrit > tcrit )
	{
		firstCompression = true ;
		secondCompression = true ;
		firstTension = false ;
		firstCompression = false ;
		secondMet = ccrit > 0 ;
		return ccrit ;
	}
	
	firstCompression = false ;
	secondCompression = false ;
	firstTension = true ;
	firstCompression = true ;
	firstMet = tcrit > 0 ;
	return tcrit ;
}

double NonLocalMCFT::getTensileLimit(const ElementState & s) const 
{
	double maxTension = upVal;
	double pseudoYoung = youngModulus*(1.-std::min(s.getParent()->getBehaviour()->getDamageModel()->getState().max(), 1.-1e-12));
	double maxTensionStrain = upVal/pseudoYoung;
	std::vector<double> crits ;
		

		double downTestVal = 0 ;
		double upTestVal = upVal ;
		double factor = 1 ;
		double delta_tech = strain_te-strain_ch;
		while(std::abs(upTestVal-downTestVal) > 1e-14*upVal)
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

} ;

FractureCriterion *NonLocalMCFT::getCopy() const
{
	return new NonLocalMCFT( *this ) ;
}

Material NonLocalMCFT::toMaterial()
{
	Material mat ;
	return mat ;
}

}
