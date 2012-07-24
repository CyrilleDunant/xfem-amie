
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
	, upVal( 1e6*.65*pow(std::abs(down*1e-6),.33) ), downVal( down ), youngModulus(youngModulus), rtype(r)
{
	physicalCharacteristicRadius = charRad ;
	critStrain = -0.0015 ;//-0.0015;
	tensionCritStrain = upVal / youngModulus ;
	firstCompression = false ;
	secondCompression = false ;
	firstTension = false ;
	secondTension = false ;
	firstMet = false ;
	secondMet = false ;
	initialised = false ;
	radiusInitialised = false ;
	scaleFactor = 1 ;
}


NonLocalMCFT::~NonLocalMCFT()
{
}

double NonLocalMCFT::getBareConcreteTensileCriterion(const ElementState & s, double pseudoYoung, double tstrain, double tstress)
{
	double maxTension = upVal*scaleFactor;
	double maxTensionStrain = tensionCritStrain;
	
	if(tstrain > tensionCritStrain )
	{
		
		double downTestVal = 0 ;
		double upTestVal = upVal ;
		double factor = 1 ;
		double delta_tech = strain_te-strain_ch;
		while(std::abs(upTestVal-downTestVal) > 1e-6*upVal)
		{
			
			double testVal = (upTestVal+downTestVal)*.5/(pseudoYoung) ;
			double mainCurve = 1./(1.+sqrt(k*(testVal-tensionCritStrain))) ;
			if(testVal < tensionCritStrain)
				factor = 1. ;
			else if(testVal < strain_ch)
				factor = mainCurve ;
			else if(testVal < strain_te)
				factor = mainCurve*(strain_te-testVal)/delta_tech ;
			else
				factor = 0 ;
		
			if( testVal*pseudoYoung > upVal*factor )
				upTestVal = testVal*pseudoYoung ;
			else
				downTestVal = testVal*pseudoYoung ;
		}
		
		
		maxTensionStrain = (upTestVal+downTestVal)*.5/(pseudoYoung) ;
		if(factor < POINT_TOLERANCE_2D)
			return POINT_TOLERANCE_2D ;
	}
	
	double criterion = 0 ;
  if(maxTensionStrain > POINT_TOLERANCE_2D)
		criterion = std::abs(tstrain/maxTensionStrain) ;
	
	if(criterion >= 1)
	{
		return 1.-1./criterion ;
	}
	else if(criterion >= 0)
	{
		return -1.+criterion ;
	}
	else if(criterion > -1)
	{
		return criterion ;
	}
	return POINT_TOLERANCE_2D ;
}

double NonLocalMCFT::getRebarConcreteTensileCriterion(const Mu::ElementState& s, double pseudoYoung, double tstrain, double tstress, double value)
{
	double maxTension = upVal*scaleFactor;
	double maxTensionStrain = tensionCritStrain;
	
	if(tstrain > tensionCritStrain )
	{
// 		for(double i = 0 ; i < 1. ; i += 0.001)
// 		{
		
			double downTestVal = 0 ;
			double upTestVal = upVal*1.2 ;
			double factor = 1 ;
			double delta_tech = strain_te-strain_ch;
			double mainCurve = 0 ;
			while(std::abs(upTestVal-downTestVal) > 1e-4*upVal)
			{
				double testVal = (upTestVal+downTestVal)*.5/(pseudoYoung) ;
				mainCurve = 1./(1.+sqrt(value*testVal)) ;

				if(testVal < tensionCritStrain)
					factor = 1. ;
				else if(testVal < strain_ch)
					factor = mainCurve ;
				else if(testVal < strain_te)
					factor = mainCurve*(strain_te-testVal)/delta_tech ;
				else
					factor = 0 ;
				

				if( testVal*pseudoYoung > upVal*factor )
					upTestVal = testVal*pseudoYoung ;
				else
					downTestVal = testVal*pseudoYoung ;

			}
		
			maxTension = (upTestVal+downTestVal)*.5*scaleFactor ;
			maxTensionStrain = (upTestVal+downTestVal)*.5/(pseudoYoung) ;
// 		}
// 		exit(0) ;
		if(factor < POINT_TOLERANCE_2D)
			return POINT_TOLERANCE_2D ;
	}
// 	std::cout << maxTensionStrain*1e6 << ", "<< maxTension*1e-6 << ",  "<< pseudoYoung*1e-9 << std::endl ;
	double criterion = 0 ;
	/*if(maxTensionStrain > POINT_TOLERANCE_2D && maxTension > POINT_TOLERANCE_2D )
	{
		if(rtype == UPPER_BOUND)
			criterion = std::min(std::abs(tstress/maxTension), std::abs(tstrain/maxTensionStrain)) ;
		else if (rtype == LOWER_BOUND)
			criterion = std::max(std::abs(tstress/maxTension), std::abs(tstrain/maxTensionStrain)) ;
		else
			criterion = (std::abs(tstress/maxTension)+ std::abs(tstrain/maxTensionStrain))*0.5 ;
	}
	else */if(maxTensionStrain > POINT_TOLERANCE_2D)
		criterion = std::abs(tstrain/maxTensionStrain) ;
// 	else if(maxTension > POINT_TOLERANCE_2D )
// 		criterion = std::abs(tstress/maxTension) ;
	
	if(criterion >= 1)
	{
		return 1.-1./criterion ;
	}
	else if(criterion >= 0)
	{
		return -1.+criterion ;
	}
	else if(criterion > -1)
	{
		return criterion ;
	}
	return POINT_TOLERANCE_2D ;
}

double NonLocalMCFT::getConcreteTensileCriterion(const ElementState & s, double pseudoYoung, double tstrain, double tstress)
{
// 	if(dynamic_cast<RotatingCrack *>(s.getParent()->getBehaviour()->getDamageModel()))
// 	{
// 		if(dynamic_cast<RotatingCrack *>(s.getParent()->getBehaviour()->getDamageModel())->tensionFailure)
// 			return POINT_TOLERANCE_2D ;
// 	}
	bool inRebarInfluence = false ;
	double distanceToRebar = -1 ;
	double effectiveInfluenceDistance = -1 ;
	if(!rebarLocationsAndDiameters.empty())
	{
		distanceToRebar = std::abs(s.getParent()->getCenter().y - rebarLocationsAndDiameters[0].first) ;
		effectiveInfluenceDistance =  rebarLocationsAndDiameters[0].second*7.5*2.;
		inRebarInfluence = distanceToRebar < rebarLocationsAndDiameters[0].second*7.5*2. ;
// 		std::cout << rebarLocationsAndDiameters[0].first << "  "<<std::flush ;
		for(size_t i = 1 ; i < rebarLocationsAndDiameters.size() ; i++)
		{
			double distanceToRebartest = std::abs(s.getParent()->getCenter().y - rebarLocationsAndDiameters[i].first) ;
// 			std::cout << rebarLocationsAndDiameters[i].first << "  "<<std::flush ;
			if( distanceToRebar < rebarLocationsAndDiameters[i].second*7.5*2. && distanceToRebartest < distanceToRebar)
			{
				distanceToRebar = distanceToRebartest ;
				effectiveInfluenceDistance = rebarLocationsAndDiameters[i].second*7.5*2. ;
				inRebarInfluence = true ;
			}
		}
// 		std::cout << std::endl ;
	}
// 	if(inRebarInfluence)
// 	  return -1 ;
/*	  std::cout << s.getParent()->getCenter().y << std::endl */;
	double barecrit = getBareConcreteTensileCriterion(s, pseudoYoung, tstrain, tstress) ;
	
	if(!inRebarInfluence)
		return barecrit ;

 	double rebcrit = getRebarConcreteTensileCriterion(s, pseudoYoung, tstrain, tstress) ;
//	return rebcrit ;


// 	if (!radiusInitialised)
// 	{
// 		double f = (distanceToRebar)/(effectiveInfluenceDistance) ;
// 		double df = 3.*f*f-2.*f*f*f ;
// 		setMaterialCharacteristicRadius(getMaterialCharacteristicRadius()*df+.5*getMaterialCharacteristicRadius()*(1.-df));
// 		radiusInitialised = true ;
// 	}

	double f = (distanceToRebar)/(effectiveInfluenceDistance) ;
	double df = f ;
//	double df = 3.*f*f-2.*f*f*f ;
//
// 	return barecrit ;
	return df*barecrit + (1.-df)*rebcrit ;
}

double NonLocalMCFT::getConcreteCompressiveCriterion(const ElementState & s, double pseudoYoung, double cstrain, double tstress, double cstress)
{
// 	if(dynamic_cast<RotatingCrack *>(s.getParent()->getBehaviour()->getDamageModel()))
// 	{
// 		if(dynamic_cast<RotatingCrack *>(s.getParent()->getBehaviour()->getDamageModel())->compressionFailure)
// 			return POINT_TOLERANCE_2D ;
// 	}
	double maxCompression = downVal*scaleFactor  ;
	double maxCompressionStrain = downVal/pseudoYoung  ;

	if(cstrain < 0.05*critStrain )
	{
		double C_d = 0. ;
		double compressiveTensileRatio = -std::abs(tstress/std::min(cstress, -POINT_TOLERANCE_2D)) ;
		
		if(-compressiveTensileRatio > 0.280001)
			C_d =0.35*pow(-compressiveTensileRatio-0.28, 0.8) ;
		
		double beta_d = 1./(1.+C_d) ;
		
		if(beta_d > 1)
			beta_d = 1 ;
		
		double f_p = beta_d*downVal ;
		double epsilon_p = beta_d*critStrain ;
		double n = 0.8 - f_p/17e6 ;
		double k_c = 0.67 - f_p/62e6 ;

		if(n*f_p/(epsilon_p*pseudoYoung)-n+1. >= 0)
		{
			
			double rcs = pow(n*f_p/(epsilon_p*pseudoYoung)-n+1., 1./(n*k_c)) ;

			if(rcs*pseudoYoung*pseudoYoung > 1)
				k_c = 1. ;
			else
				k_c = 0.67 - f_p/62e6 ;
			
			rcs = pow(n*f_p/(epsilon_p*pseudoYoung)-n+1., 1./(n*k_c)) ;

			maxCompressionStrain = rcs*epsilon_p ;
		}
	}
	
	double criterion = 0 ;

	if(std::abs(maxCompressionStrain) > POINT_TOLERANCE_2D)
		criterion = std::abs(cstrain/maxCompressionStrain) ;

	
	if(criterion >= 1)
	{
		return 1.-1./criterion ;
	}
	else if(criterion >= 0)
	{
		return -1.+criterion ;
	}
	else if(criterion > -1)
	{
		return criterion ;
	}
	return POINT_TOLERANCE_2D ;
}

void NonLocalMCFT::initialise( ElementState &s)
{
	double stressFactor = 1 ;
	double strainFactor = 1 ;
	double energy = 75. ; //N/m
	double damagedfraction = 5.;
	strain_ch = 2.*energy/(getMaterialCharacteristicRadius()*damagedfraction*upVal*stressFactor) ;//*.5 energy <- // *2 energy -> 2.*energy/(1.*getMaterialCharacteristicRadius()*upVal) ;
	
	if(strain_ch < tensionCritStrain)
	{
		std::cout << strain_ch << " vs " << tensionCritStrain <<std::endl ;
		exit(0) ;
	}
	strain_te = 5.*strain_ch;
	double del_0 = (strain_ch-tensionCritStrain*strainFactor) ;
	double del_1 = (strain_te-strain_ch) ;
	
	double k0 = 1e6 ;
	double k_low = 0 ;
	double k_high = k0 ;
	double elastic_energy = tensionCritStrain*upVal*.5*stressFactor ;
	double integral = elastic_energy ;
	do
	{
		integral = elastic_energy ;
		k = 0.5*(k_low+k_high) ;
		      
		for(double i = 0 ; i < 10000 ; i++)
		{
			integral+= (upVal*stressFactor/(1.+sqrt(k*(i)/10000.*del_0)) + upVal*stressFactor/(1.+sqrt(k*(i+1.)/10000.*del_0)))*del_0*1e-4 ;
		}
		      
		for(double i = 0 ; i < 1000 ; i++)
		{  
			integral+= (upVal*stressFactor/(1.+sqrt(k*((i)/10000.*del_1+del_0)))+upVal*stressFactor/(1.+sqrt(k*((i+1.)/10000.*del_1+del_0))))*del_1*1e-4 ;
		}
		if(integral < energy*damagedfraction)
			k_high = k ;
		else
			k_low = k ;
		
	} while(std::abs(k_low-k_high) > 1e-8*k0) ;
	
	if(std::abs(integral-energy*damagedfraction) > 1e-3*energy*damagedfraction)
	{
	  std::cout << "wrong energy! " << integral << " vs "<< energy*damagedfraction  << std::endl ;
	  exit(0) ;
	}
	initialised = true ;
}

double NonLocalMCFT::grade( ElementState &s )
{
	if(!initialised)
		initialise(s);
	
	std::pair<Vector, Vector> stressStrain = smoothedPrincipalStressAndStrain(s, REAL_STRESS) ;
	double tstrain = stressStrain.second[0];
	double cstrain = stressStrain.second[1];
	double tstress = stressStrain.first[0];
	double cstress = stressStrain.first[1];

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
		if(firstCompression && cpseudoYoung0 > POINT_TOLERANCE_2D)
		{
			ccrit0 = getConcreteCompressiveCriterion(s, cpseudoYoung0, cstrain, tstress, cstress) ;
			if(ccrit0 > 0)
				firstMet = true ;
		}
		if(firstTension && tpseudoYoung0 > POINT_TOLERANCE_2D)
		{
			tcrit0 = getConcreteTensileCriterion(s, tpseudoYoung0, tstrain, tstress) ;
			if(tcrit0 > 0)
				firstMet = true ;
		}
		double c0 = std::max(ccrit0, tcrit0) ;
		if(secondCompression && cpseudoYoung1 > POINT_TOLERANCE_2D)
		{
			ccrit1 = getConcreteCompressiveCriterion(s, cpseudoYoung1, cstrain, tstress, cstress) ;
			if(ccrit1 > 0)
			{
				secondMet = true ;
			}
		}
		if(secondTension && tpseudoYoung1 > POINT_TOLERANCE_2D)
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
