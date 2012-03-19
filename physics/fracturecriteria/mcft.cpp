
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
	scaleFactor = 1 ;
}


NonLocalMCFT::~NonLocalMCFT()
{
}

double NonLocalMCFT::getBareConcreteTensileCriterion(const ElementState & s, double pseudoYoung, double tstrain, double tstress)
{
	double stressFactor = 1 ;
	double strainFactor = 1 ;
	double modulusFactor = 1 ;
	double maxTension = upVal*scaleFactor;
	double maxTensionStrain = tensionCritStrain;
	
	if(tstrain > tensionCritStrain )
	{
		
		double downTestVal = 0 ;
		double upTestVal = upVal*stressFactor ;
		double factor = 1 ;
		double delta_tech = strain_te-strain_ch;
		while(std::abs(upTestVal-downTestVal) > 1e-8*upVal*stressFactor)
		{
			
			double testVal = (upTestVal+downTestVal)*.5/(pseudoYoung*modulusFactor) ;
			double mainCurve = 1./(1.+sqrt(k*(testVal-tensionCritStrain*strainFactor))) ;
			if(testVal < tensionCritStrain*strainFactor)
				factor = 1. ;
			else if(testVal < strain_ch)
				factor = mainCurve ;
			else if(testVal < strain_te)
				factor = mainCurve*(strain_te-testVal)/delta_tech ;
			else
				factor = 0 ;
		
			if( testVal*pseudoYoung*modulusFactor > upVal*factor*stressFactor )
				upTestVal = testVal*pseudoYoung*modulusFactor ;
			else
				downTestVal = testVal*pseudoYoung*modulusFactor ;
			
		}
		
		
		maxTension = (upTestVal+downTestVal)*.5*scaleFactor/stressFactor ;
		maxTensionStrain = (upTestVal+downTestVal)*.5/(strainFactor*pseudoYoung*modulusFactor) ;
		if(factor < POINT_TOLERANCE_2D)
			return 1. ;
	}
	
	double criterion = 0 ;
/*	if(maxTensionStrain > POINT_TOLERANCE_2D && maxTension > POINT_TOLERANCE_2D )
	{
		if(rtype == UPPER_BOUND)
			criterion = std::min(std::abs(tstress/maxTension), std::abs(tstrain/maxTensionStrain)) ;
		else if(rtype == LOWER_BOUND)
			criterion = std::max(std::abs(tstress/maxTension), std::abs(tstrain/maxTensionStrain)) ;
		else
			criterion = (std::abs(tstress/maxTension)+ std::abs(tstrain/maxTensionStrain))*.5 ;
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
	return 1 ;
}

double NonLocalMCFT::getRebarConcreteTensileCriterion(const Mu::ElementState& s, double pseudoYoung, double tstrain, double tstress, double value)
{
	double stressFactor = 1e-6 ;
	double strainFactor = 1 ;
	double modulusFactor = 1e-6 ;
	double maxTension = upVal*scaleFactor;
	double maxTensionStrain = tensionCritStrain;
	
	if(tstrain > tensionCritStrain )
	{
		double downTestVal = 0 ;
		double upTestVal = upVal*stressFactor ;
		double factor = 1 ;
		double delta_tech = (strain_te-strain_ch);

		while(std::abs(upTestVal-downTestVal) > 1e-8*upVal*stressFactor)
		{
			double testVal = (upTestVal+downTestVal)*.5/(pseudoYoung*modulusFactor) ;
			double mainCurve = 1./(1.+sqrt(value*testVal)) ;

			if(testVal < tensionCritStrain*strainFactor)
				factor = 1. ;
			else if(testVal < strain_ch)
				factor = mainCurve ;
			else if(testVal < strain_te)
			{
			  factor = mainCurve*(strain_te-testVal)/delta_tech ;
			}
			else
				factor = 0 ;

			if( testVal*pseudoYoung*modulusFactor > upVal*factor*stressFactor )
				upTestVal = testVal*pseudoYoung*modulusFactor ;
			else
				downTestVal = testVal*pseudoYoung*modulusFactor ;

		}
	
		maxTension = (upTestVal+downTestVal)*.5*scaleFactor/stressFactor ;
		maxTensionStrain = (upTestVal+downTestVal)*.5/(pseudoYoung*modulusFactor*strainFactor) ;
		
		if(factor < POINT_TOLERANCE_2D || maxTensionStrain > strain_te)
			return 1. ;
	}
	
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
	return 1 ;
}

double NonLocalMCFT::getConcreteTensileCriterion(const ElementState & s, double pseudoYoung, double tstrain, double tstress)
{
// 	if(dynamic_cast<RotatingCrack *>(s.getParent()->getBehaviour()->getDamageModel()))
// 	{
// 		if(dynamic_cast<RotatingCrack *>(s.getParent()->getBehaviour()->getDamageModel())->tensionFailure)
// 			return 0 ;
// 	}
	bool inRebarInfluence = false ;
	double distanceToRebar = -1 ;
	double effectiveInfluenceDistance = -1 ;
	if(!rebarLocationsAndDiameters.empty())
	{
		distanceToRebar = std::abs(s.getParent()->getCenter().y - rebarLocationsAndDiameters[0].first) ;
		effectiveInfluenceDistance =  rebarLocationsAndDiameters[0].second*7.5*1.5;
		inRebarInfluence = distanceToRebar < rebarLocationsAndDiameters[0].second*7.5*1.5 ;
		
		for(size_t i = 1 ; i < rebarLocationsAndDiameters.size() ; i++)
		{
			double distanceToRebartest = std::abs(s.getParent()->getCenter().y - rebarLocationsAndDiameters[i].first) ;
			if( distanceToRebar < rebarLocationsAndDiameters[i].second*7.5 && distanceToRebartest < distanceToRebar)
			{
				distanceToRebar = distanceToRebartest ;
				effectiveInfluenceDistance = rebarLocationsAndDiameters[i].second*7.5*1.5 ;
				inRebarInfluence = true ;
			}
		}
	}
	
 	if(!inRebarInfluence)
 	{
 		return getBareConcreteTensileCriterion(s, pseudoYoung, tstrain, tstress) ;
 	}
 	
  	if(distanceToRebar < effectiveInfluenceDistance*.666666)
 		return getRebarConcreteTensileCriterion(s, pseudoYoung, tstrain, tstress) ;
	

	double barecrit = getBareConcreteTensileCriterion(s, pseudoYoung, tstrain, tstress) ;
	double rebcrit = getRebarConcreteTensileCriterion(s, pseudoYoung, tstrain, tstress) ;
	
	double f = (distanceToRebar-effectiveInfluenceDistance*.666666)/(effectiveInfluenceDistance*.666666) ;
	double df = 3.*f*f-2.*f*f*f ;
	
	return df*barecrit + (1.-df)*rebcrit ;
}

double NonLocalMCFT::getConcreteCompressiveCriterion(const ElementState & s, double pseudoYoung, double cstrain, double tstress, double cstress)
{
// 	if(dynamic_cast<RotatingCrack *>(s.getParent()->getBehaviour()->getDamageModel()))
// 	{
// 		if(dynamic_cast<RotatingCrack *>(s.getParent()->getBehaviour()->getDamageModel())->compressionFailure)
// 			return 0 ;
// 	}
	double maxCompression = downVal*scaleFactor  ;
	double maxCompressionStrain = downVal/pseudoYoung  ;

	if(cstrain < 0.5*critStrain )
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
		
		double upTestVal = 2.*downVal ;
		double downTestVal = 0 ;

		if(n*f_p/(epsilon_p*pseudoYoung)-n+1. >= 0)
		{
			
			double rcs = pow(n*f_p/(epsilon_p*pseudoYoung)-n+1., 1./(n*k_c)) ;
			double testVal = rcs*epsilon_p*pseudoYoung ;

			if(testVal > epsilon_p/pseudoYoung)
				k_c = 1. ;
			else
				k_c = 0.67 - f_p/62e6 ;
			
			rcs = pow(n*f_p/(epsilon_p*pseudoYoung)-n+1., 1./(n*k_c)) ;
			testVal = rcs*epsilon_p*pseudoYoung ;

			maxCompression = testVal*scaleFactor ;
			maxCompressionStrain = testVal/pseudoYoung ;
		}
	}
	
	double criterion = 0 ;
/*	if(std::abs(maxCompressionStrain) > POINT_TOLERANCE_2D && std::abs(maxCompression) > POINT_TOLERANCE_2D )
	{
		if(rtype == UPPER_BOUND)
			criterion = std::min(std::abs(cstress/maxCompression), std::abs(cstrain/maxCompressionStrain)) ;
		else if (rtype == LOWER_BOUND)
			criterion = std::max(std::abs(cstress/maxCompression), std::abs(cstrain/maxCompressionStrain)) ;
		else
			criterion =(std::abs(cstress/maxCompression)+ std::abs(cstrain/maxCompressionStrain))*0.5 ;
	}
	else */if(std::abs(maxCompressionStrain) > POINT_TOLERANCE_2D)
		criterion = std::abs(cstrain/maxCompressionStrain) ;
// 	else if(std::abs(maxCompression) > POINT_TOLERANCE_2D )
// 		criterion = std::abs(cstress/maxCompression) ;
	
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
	return 1 ;
}

void NonLocalMCFT::initialise()
{
	double stressFactor = 1 ;
	double strainFactor = 1 ;
	double energy = 75. ; //N/m
	strain_ch = 2.*energy/(getMaterialCharacteristicRadius()*10.*upVal*stressFactor) ;//*.5 energy <- // *2 energy -> 2.*energy/(1.*getMaterialCharacteristicRadius()*upVal) ;
	
	if(strain_ch < tensionCritStrain)
	{
		std::cout << strain_ch << " vs " << tensionCritStrain <<std::endl ;
		exit(0) ;
	}
	strain_te = 5.*strain_ch;
	double del_0 = (strain_ch-tensionCritStrain*strainFactor) ;
	double del_1 = (strain_te-strain_ch) ;
	
	double k_low = 0 ;
	double k_high = 1e10 ;
	double elastic_energy = tensionCritStrain*upVal*.5*stressFactor ;
	double integral = elastic_energy ;
	do
	{
		integral = 0 ;
		k = 0.5*(k_low+k_high) ;
		      
		for(double i = 0 ; i < 10000 ; i++)
		{
			integral+= (upVal*stressFactor/(1.+sqrt(k*(i)/10000.*del_0)) + upVal*stressFactor/(1.+sqrt(k*(i+1.)/10000.*del_0)))*del_0*1e-4 ;
		}
		      
		for(double i = 0 ; i < 1000 ; i++)
		{  
			integral+= (upVal*stressFactor/(1.+sqrt(k*((i)/10000.*del_1+del_0)))+upVal*stressFactor/(1.+sqrt(k*((i+1.)/10000.*del_1+del_0))))*del_1*1e-4 ;
		}
		if(integral < energy)
			k_high = k ;
		else
			k_low = k ;
		
	} while(std::abs(k_low-k_high) > 1e-8*1e10) ;
	
	if(std::abs(integral-energy) > 1e-3*energy)
	{
	  std::cout << "wrong energy! " << integral << std::endl ;
	  exit(0) ;
	}
	initialised = true ;
}

double NonLocalMCFT::grade( ElementState &s )
{
	if(!initialised)
		initialise();
	
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
		else if(c1 > c0)
		{
			firstMet = false ;
			secondMet = true ;
		}
		else
		{
			firstMet = true ;
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
		return ccrit ;
	}
	
	firstCompression = false ;
	secondCompression = false ;
	firstTension = true ;
	firstCompression = true ;
	return tcrit ;
	return -1 ;
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
