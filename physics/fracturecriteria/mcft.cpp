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

MCFT::MCFT( double up, double down, double youngModulus, double charRad, MirrorState mirroring, double delta_x, double delta_y, double delta_z ) : FractureCriterion( mirroring, delta_x, delta_y, delta_z )
	, upVal( up ), downVal( down )
{
	critStrain = -0.0015;
	tensionCritStrain = upVal / youngModulus ;
	strainBroken = false ;
	initialised = false ;
}


MCFT::~MCFT()
{
}

double MCFT::grade( ElementState &s )
{
		if(!initialised)
	{
// 		if(s.getParent()->getRadius() > 0.5*getMaterialCharacteristicRadius())
// 		{
// 			std::cout << "too large elements!" << std::endl ;
// 			exit(0) ;
// 		}
		
		double energy = 75. ; //N/m
		strain_ch = 2.*energy/(2.*2.*getMaterialCharacteristicRadius()*upVal) ;
		
		if(strain_ch < tensionCritStrain)
		{
			std::cout << strain_ch << " vs " << tensionCritStrain <<std::endl ;
			exit(0) ;
		}
		strain_te = 5.*strain_ch;
		double del_0 = strain_ch-tensionCritStrain ;
		double del_1 = strain_te-strain_ch ;
		
		double k_low = 0 ;
		double k_high = 1e5 ;
		double elastic_energy = tensionCritStrain*upVal*.5 ;
		
		do
		{
			double integral = elastic_energy ;
			k = 0.5*(k_low+k_high) ;
			
			for(double i = 0 ; i < 20000 ; i++)
			{
				integral+= 0.5*(upVal/(1.+sqrt(k*(i)/20000.*del_0)) + upVal/(1.+sqrt(k*(i+1)/20000.*del_0)))*del_0*0.5e-4 ;
			}
			
			for(double i = 0 ; i < 20000 ; i++)
			{
				integral+= 0.5*(upVal/(1.+sqrt(k*((i)/20000.*del_1+del_0)))+upVal/(1.+sqrt(k*((i+1)/20000.*del_1+del_0))))*del_1*0.5e-4 ;
			}
			
			if(integral < energy)
			{
				k_high = k ;
			}
			else
				k_low = k ;
			
		} while(std::abs(k_low-k_high) > 1e-6) ;
		initialised = true ;
	}
	
	
	Vector pstrain = s.getPrincipalStrains( s.getParent()->getCenter() ) ;
	Vector pstress = s.getPrincipalStresses( s.getParent()->getCenter() ) ;

	double tstrain = pstrain.max();
	double cstrain = pstrain.min();
	double tstress = pstress.max();
	double cstress = pstress.min();

	double pseudoYoung = youngModulus*(1.-std::min(s.getParent()->getBehaviour()->getDamageModel()->getState().max(), 1.-1e-12));
	if(pseudoYoung < 1e-9)
		return -1 ;
	double maxCompression = downVal  ;
	double maxCompressionStrain = downVal/pseudoYoung  ;

	if(cstrain < critStrain*.45 )
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

		if(n*f_p/(epsilon_p*pseudoYoung)-n+1. < 0)
			return -1 ;
		double rcs = pow(n*f_p/(epsilon_p*pseudoYoung)-n+1., 1./(n*k_c)) ;
		double testVal = rcs*epsilon_p*pseudoYoung ;

		if(testVal > epsilon_p/pseudoYoung)
			k_c = 1. ;
		else
			k_c = 0.67 - f_p/62e6 ;
		rcs = pow(n*f_p/(epsilon_p*pseudoYoung)-n+1., 1./(n*k_c)) ;
		testVal = rcs*epsilon_p*pseudoYoung ;

		maxCompression = testVal ;
		maxCompressionStrain = testVal/pseudoYoung ;
	}
	
	
	double maxTension = upVal;
	double maxTensionStrain = upVal/pseudoYoung;
	std::vector<double> crits ;
		
	if(tstrain > tensionCritStrain*0.45 )
	{
		strainBroken = true ;
// 		pseudoYoung = std::max(tstress/tstrain, pseudoYoung) ;
// 		std::cout << 0 << "  " << 0 << std::endl ;
// 		for(double i = 0.0 ; i < 1 ; i += 0.001)
// 		{
// 		pseudoYoung = youngModulus*(1.-i) ;
		double downTestVal = 0 ;
		double upTestVal = upVal ;
		double factor = 1 ;
		double delta_tech = strain_te-strain_ch;
		while(std::abs(upTestVal-downTestVal) > 1e-6*upVal)
		{
			
			double testVal = (upTestVal+downTestVal)*.5/pseudoYoung ;
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
			{
				upTestVal = testVal*pseudoYoung ;
			}
			else
			{
				downTestVal = testVal*pseudoYoung ;
			}
			
		}
		
		
		maxTension = (upTestVal+downTestVal)*.5 ;
		maxTensionStrain = (upTestVal+downTestVal)*.5/pseudoYoung ;
// 		std::cout << maxTensionStrain<< "  " << maxTension << std::endl ;
// 		}
// 		exit(0) ;
		if(factor < POINT_TOLERANCE_2D)
			return -1. ;

	}



if(cstrain < 0 && tstrain < 0)
	metInCompression = true ;
else if(cstrain > 0 && tstrain > 0)
	metInTension = true ;
else
{
	metInCompression = std::abs(cstress/cstrain) > std::abs(tstress/tstrain) ;
	metInTension = !metInCompression;
}

if(maxCompression < 0 &&  maxCompressionStrain < 0 &&  metInCompression)
{
	crits.push_back(std::min(cstress/maxCompression, cstrain/maxCompressionStrain)) ;
}

if(maxTensionStrain > 0 && metInTension && tstrain > tensionCritStrain )
	crits.push_back(std::min(tstrain/maxTensionStrain, tstress/maxTension)) ;

if(crits.empty())
	return -1 ;

std::sort( crits.begin(), crits.end() );

if(crits.back() >= 1)
{
	return 1.-1./crits.back() ;
}
else if(crits.back() >= 0)
{
	return -1.+crits.back() ;
}
if(crits.back() > -1)
{
	return crits.back() ;
}

return -1 ;

}

FractureCriterion *MCFT::getCopy() const
{
	return new MCFT( *this ) ;
}

Material MCFT::toMaterial()
{
	Material mat ;
	return mat ;
}

double MCFT::getTensileLimit(const ElementState & s) const 
{
	double maxTension = upVal;
	double pseudoYoung = youngModulus*(1.-std::min(s.getParent()->getBehaviour()->getDamageModel()->getState().max(), 1.-1e-12));
	double maxTensionStrain = tensionCritStrain;
	std::vector<double> crits ;
		

		double downTestVal = 0 ;
		double upTestVal = upVal ;
		double factor = 1 ;
		double delta_tech = strain_te-strain_ch;
		while(std::abs(upTestVal-downTestVal) > 1e-14*upVal)
		{
			
			double testVal = (upTestVal+downTestVal)*.5/pseudoYoung ;
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


NonLocalMCFT::NonLocalMCFT( double up, double down, double youngModulus,  double charRad, RedistributionType r, MirrorState mirroring, double delta_x, double delta_y, double delta_z ) : FractureCriterion( mirroring, delta_x, delta_y, delta_z )
	, upVal( up ), downVal( down ), youngModulus(youngModulus), rtype(r)
{
	physicalCharacteristicRadius = charRad ;
	critStrain = -0.0015 ;//-0.0015;
	tensionCritStrain = upVal / youngModulus ;
	strainBroken = false ;
	initialised = false ;
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
		while(std::abs(upTestVal-downTestVal) > 1e-8*upVal)
		{
			double testVal = (upTestVal+downTestVal)*.5/pseudoYoung ;
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
		
		
		maxTension = (upTestVal+downTestVal)*.5*scaleFactor ;
		maxTensionStrain = (upTestVal+downTestVal)*.5/pseudoYoung ;
		
		if(factor < POINT_TOLERANCE_2D)
			return 1. ;
	}
	
	double criterion = 0 ;
	/*if(maxTensionStrain > POINT_TOLERANCE_2D && maxTension > POINT_TOLERANCE_2D )
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

double NonLocalMCFT::getRebarConcreteTensileCriterion(const ElementState & s, double pseudoYoung, double tstrain, double tstress)
{
	double maxTension = upVal*scaleFactor;
	double maxTensionStrain = tensionCritStrain;
	
	if(tstrain > tensionCritStrain )
	{
		double downTestVal = 0 ;
		double upTestVal = upVal ;
		double factor = 1 ;
		double delta_tech = strain_te-strain_ch;

		while(std::abs(upTestVal-downTestVal) > 1e-8*upVal)
		{
			double testVal = (upTestVal+downTestVal)*.5/pseudoYoung ;
			double mainCurve = 1./(1.+sqrt(500*testVal)) ;

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
		maxTensionStrain = (upTestVal+downTestVal)*.5/pseudoYoung ;
		
		if(factor < POINT_TOLERANCE_2D || maxTensionStrain > strain_te)
			return 1. ;
	}
	
	double criterion = 0 ;
/*	if(maxTensionStrain > POINT_TOLERANCE_2D && maxTension > POINT_TOLERANCE_2D )
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
		effectiveInfluenceDistance =  rebarLocationsAndDiameters[0].second*7.5*0.5;
		inRebarInfluence = distanceToRebar < rebarLocationsAndDiameters[0].second*7.5*0.5 ;
		
		for(size_t i = 1 ; i < rebarLocationsAndDiameters.size() ; i++)
		{
			double distanceToRebartest = std::abs(s.getParent()->getCenter().y - rebarLocationsAndDiameters[i].first) ;
			if( distanceToRebar < rebarLocationsAndDiameters[i].second*7.5*0.5 && distanceToRebartest < distanceToRebar)
			{
				distanceToRebar = distanceToRebartest ;
				effectiveInfluenceDistance = rebarLocationsAndDiameters[i].second*7.5*0.5 ;
				inRebarInfluence = true ;
			}
		}
	}
	
	if(!inRebarInfluence)
		return getBareConcreteTensileCriterion(s, pseudoYoung, tstrain, tstress) ;
	
 	if(distanceToRebar < effectiveInfluenceDistance/(0.5*7.5*.5))
 		return getRebarConcreteTensileCriterion(s, pseudoYoung, tstrain, tstress) ;
	

	double barecrit = getBareConcreteTensileCriterion(s, pseudoYoung, tstrain, tstress) ;
	double rebcrit = getRebarConcreteTensileCriterion(s, pseudoYoung, tstrain, tstress) ;
	
	double f = distanceToRebar/effectiveInfluenceDistance ;
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
	if(std::abs(maxCompressionStrain) > POINT_TOLERANCE_2D && std::abs(maxCompression) > POINT_TOLERANCE_2D )
	{
		if(rtype == UPPER_BOUND)
			criterion = std::min(std::abs(cstress/maxCompression), std::abs(cstrain/maxCompressionStrain)) ;
		else if (rtype == LOWER_BOUND)
			criterion = std::max(std::abs(cstress/maxCompression), std::abs(cstrain/maxCompressionStrain)) ;
		else
			criterion =(std::abs(cstress/maxCompression)+ std::abs(cstrain/maxCompressionStrain))*0.5 ;
	}
	else if(std::abs(maxCompressionStrain) > POINT_TOLERANCE_2D)
		criterion = std::abs(cstrain/maxCompressionStrain) ;
	else if(std::abs(maxCompression) > POINT_TOLERANCE_2D )
		criterion = std::abs(cstress/maxCompression) ;
	
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
		double energy = 75. ; //N/m
		strain_ch = 2.*energy/(2.*0.15*upVal) ;//*.5 energy <- // *2 energy -> 2.*energy/(1.*getMaterialCharacteristicRadius()*upVal) ;
		
		if(strain_ch < tensionCritStrain)
		{
			std::cout << strain_ch << " vs " << tensionCritStrain <<std::endl ;
			exit(0) ;
		}
		strain_te = 5.*strain_ch;
		double del_0 = strain_ch-tensionCritStrain ;
		double del_1 = strain_te-strain_ch ;
		
		double k_low = 0 ;
		double k_high = 1e5 ;
		double elastic_energy = tensionCritStrain*upVal*.5 ;
		
		do
		{
			double integral = elastic_energy ;
			k = 0.5*(k_low+k_high) ;
			
			for(double i = 0 ; i < 40000 ; i++)
				integral+= 0.5*(upVal/(1.+sqrt(k*(i)/40000.*del_0)) + upVal/(1.+sqrt(k*(i+1)/40000.*del_0)))*del_0*0.25e-4 ;
			
			for(double i = 0 ; i < 40000 ; i++)
				integral+= 0.5*(upVal/(1.+sqrt(k*((i)/40000.*del_1+del_0)))+upVal/(1.+sqrt(k*((i+1)/40000.*del_1+del_0))))*del_1*0.25e-4 ;
			
			if(integral < energy)
				k_high = k ;
			else
				k_low = k ;
			
		} while(std::abs(k_low-k_high) > 1e-6) ;
		
		initialised = true ;
}

double NonLocalMCFT::grade( ElementState &s )
{
	if(!initialised)
		initialise();
	
	std::pair<Vector, Vector> stressStrain = smoothedPrincipalStressAndStrain(s, REAL_STRESS) ;
	double tstrain = stressStrain.second.max();
	double cstrain = stressStrain.second.min();
	double tstress = stressStrain.first.max();
	double cstress = stressStrain.first.min();

	double pseudoYoung = youngModulus*(1.-s.getParent()->getBehaviour()->getDamageModel()->getState().max());
	
// 	if(pseudoYoung < 1e-9)
// 	{
// 		metInCompression = true ;
// 		metInTension = true ;
// 		return -1 ;
// 	}
	if(dynamic_cast<RotatingCrack *>(s.getParent()->getBehaviour()->getDamageModel()))
	{
		double cpseudoYoung =  youngModulus*(1.-s.getParent()->getBehaviour()->getDamageModel()->getState()[1]) ;
		double tpseudoYoung =  youngModulus*(1.-s.getParent()->getBehaviour()->getDamageModel()->getState()[0]) ;
		if(cpseudoYoung > 1e-6 && tpseudoYoung> 1e-6 )
		{
			double ccrit = getConcreteCompressiveCriterion(s, cpseudoYoung, cstrain, tstress, cstress) ;
			double tcrit = getConcreteTensileCriterion(s, tpseudoYoung, tstrain, tstress) ;
// 			std::cout << "tc " << ccrit << ", " << s.getParent()->getBehaviour()->getDamageModel()->getState()[1] << " : " << tcrit << ", "<< s.getParent()->getBehaviour()->getDamageModel()->getState()[0] << std::endl ;
			if( ccrit > tcrit )
			{
				metInCompression = true ;
				metInTension = false ;
				return ccrit ;
			}
			
			metInCompression = false ;
			metInTension = true ;
			return tcrit ;
		}
		else if(tpseudoYoung> 1e-6)
		{
			double tcrit = getConcreteTensileCriterion(s, tpseudoYoung, tstrain, tstress) ;
// 			std::cout << "tc " << tcrit << std::endl ;
			metInCompression = false ;
			metInTension = true ;
			return tcrit ;
		}
		else if (cpseudoYoung > 1e-6)
		{
			double ccrit = getConcreteCompressiveCriterion(s, cpseudoYoung, cstrain, tstress, cstress) ;
// 			std::cout << "c " << ccrit << "  " << std::endl ;
			metInCompression = true ;
			metInTension = false ;
			return ccrit ;
		}
// 		else if (cpseudoYoung > 1e-6 && tstress < 0)
// 		{
// 			
// 			double ccrit = getConcreteCompressiveCriterion(s, cpseudoYoung, tstrain, tstress, tstress) ;
// 			// 			std::cout << "c " << ccrit << "  " << std::endl ;
// 			metInCompression = true ;
// 			metInTension = false ;
// 			return ccrit ;
// 		}
// 		std::cout << "0 " << -1 << std::endl ;
		metInCompression = true ;
		metInTension = true ;
		return 1 ;
	}
	
	double ccrit = getConcreteCompressiveCriterion(s, pseudoYoung, cstrain, tstress, cstress) ;
	double tcrit = getConcreteTensileCriterion(s, pseudoYoung, tstrain, tstress) ;
	
	if( ccrit > tcrit )
	{
		metInCompression = true ;
		metInTension = false ;
		return ccrit ;
	}
	
	metInCompression = false ;
	metInTension = true ;
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
