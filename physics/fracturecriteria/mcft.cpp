
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
	, upVal( 0.65*1e6*.65*pow(std::abs(down*1e-6),.33) ), downVal( down ), youngModulus(youngModulus), rtype(r)
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
	
	inRebarInfluence = false ;
	distanceToRebar = -1 ;
	effectiveInfluenceDistance = -1 ;
}


NonLocalMCFT::~NonLocalMCFT()
{
}

double NonLocalMCFT::getBareConcreteTensileCriterion(const ElementState & s, double pseudoYoung, double tstrain, double tstress)
{
	double maxTension = upVal*scaleFactor;
	double maxTensionStrain = tensionCritStrain;
	if(std::abs(pseudoYoung-youngModulus) < POINT_TOLERANCE_2D)
	{
		maxTensionStrain = /*1.6**/tensionCritStrain ;
		maxTension = /*1.6**/maxTension ;
	}
	
	
	if(tstrain > /*1.6**/tensionCritStrain && std::abs(pseudoYoung-youngModulus) < POINT_TOLERANCE_2D || tstrain > tensionCritStrain && std::abs(pseudoYoung-youngModulus) >= POINT_TOLERANCE_2D)
	{
		double downTestVal = tensionCritStrain ;
		double upTestVal = strain_te*1.1 ;
		double factor = 1 ;
		double delta_tech = strain_te-strain_ch;
		int count = 0 ;
		while(std::abs(upTestVal-downTestVal) > 1e-8*tensionCritStrain && count++  <  32 )
		{
			double testVal = (upTestVal+downTestVal)*.5 ;
			double mainCurve = 1./(1.+sqrt(k)*pow(testVal-tensionCritStrain,0.8)) ;
// 
// 			factor = mainCurve ;
// 
// 			if( testVal > upVal*factor/pseudoYoung )
// 				upTestVal = testVal ;
// 			else
// 				downTestVal = testVal ;
			
			
// 			double testVal = (upTestVal+downTestVal)*.5 ;
// 			double mainCurve = 1./(1.+sqrt(50000)*pow(testVal, 0.8)) ;

			if(testVal < tensionCritStrain)
				factor = 1. ;
			else if(testVal < strain_ch)
				factor = mainCurve ;
			else if(testVal < strain_te)
// 			{
				factor = mainCurve*(strain_te-testVal)/delta_tech ;
// 				factor = exp(-175.*(testVal-tensionCritStrain)) ;
// 				factor = std::min((strain_te-testVal)/(strain_te-tensionCritStrain), mainCurve) ;
// 			}
			else
				factor = 0 ;
			
			
			if( testVal*pseudoYoung > upVal*factor )
				upTestVal = testVal ;
			else
				downTestVal = testVal ;
			
		}
		
// 		if((upTestVal+downTestVal)*.5 > strain_ch)
// 		{
// 			downTestVal = strain_ch ;
// 			upTestVal = strain_te*1.1 ;
// 			factor = 1 ;
// 			delta_tech = strain_te-strain_ch;
// 			while(std::abs(upTestVal-downTestVal) > 1e-8*tensionCritStrain)
// 			{
// 				double testVal = (upTestVal+downTestVal)*.5 ;
// 				double mainCurve = 1./(1.+sqrt(k*(testVal-tensionCritStrain))) ;
// 
// 				factor = mainCurve*(strain_te-testVal)/delta_tech ;
// 			
// 				if( testVal > upVal*factor/pseudoYoung )
// 					upTestVal = testVal ;
// 				else
// 					downTestVal = testVal ;
// 			}
// 		}
		
// 		if((upTestVal+downTestVal)*.5 > strain_te)
// 			factor = 0 ;
	
		maxTensionStrain = (upTestVal+downTestVal)*.5 ;
		
		if(factor < POINT_TOLERANCE_2D)
			return POINT_TOLERANCE_2D ;

	}
	
	double criterion = 0 ;
  if(maxTensionStrain > POINT_TOLERANCE_2D)
		criterion = std::abs(tstress/(maxTensionStrain*pseudoYoung)) ;

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
	return stressc ; //stressc; // strain criterion misnamed
	
	return POINT_TOLERANCE_2D ;
}

double NonLocalMCFT::getRebarConcreteTensileCriterion(const Mu::ElementState& s, double pseudoYoung, double tstrain, double tstress, double value)
{
	double delta = currentAngle ; // this assumes horizontal rebars
	double epsilon_0 = 2.*tensionCritStrain*cos(delta)*cos(delta) ; // strain at which the stifening activates
//	double epsilon_ys = 0.00125 ;
//	double concreteCover = .04 ; 
//	double barSpacing = .06 ;
//	double barPropCoef = 0.4 ; //plain bars (0.8 for deformed bars)
//	double strainGradientCoef = .25 ;
//	double rho_cf = 0.01; //A_steel/Aconcrete ;
//	double w = 2.*(concreteCover+barSpacing*.1)+barPropCoef*strainGradientCoef*rebarLocationsAndDiameters[0].second/rho_cf ;
//	double G = 50. ;
	
//	double epsilon_cu = 2.*cos(delta)*cos(delta)*G/(w*upVal) ;
//	double epsilon_u = epsilon_ys - upVal/(rho_cf*40e9) ;
	
//	double maxTension = upVal*scaleFactor;
//	double maxTensionStrain = tensionCritStrain;
//	if(std::abs(pseudoYoung-youngModulus) < POINT_TOLERANCE_2D)
//	{
//		maxTensionStrain = tensionCritStrain ;
//		maxTension = maxTension ;
//	}
	
//	double interactionStress = 0 ;
//	
//	if(tstrain > epsilon_0)
//	{
//		if(tstrain < epsilon_cu)
//		{
//			interactionStress = upVal*(tstrain-epsilon_0)/(epsilon_cu-epsilon_0) ;
//		}
//		else if(tstrain < epsilon_u)
//			interactionStress = upVal ;
//		else if(tstrain < epsilon_ys )
//			interactionStress = upVal*(1.-(tstrain-epsilon_u)/(epsilon_ys-epsilon_u)) ;
//		else
//			interactionStress = 0 ;
//	}
//	
	
//	if(tstrain > tensionCritStrain && std::abs(pseudoYoung-youngModulus) < POINT_TOLERANCE_2D || tstrain > tensionCritStrain && std::abs(pseudoYoung-youngModulus) >= POINT_TOLERANCE_2D)
//	{
//		double downTestVal = tensionCritStrain ;
//		double upTestVal = strain_te*1.1 ;
//		double factor = 1 ;
//		double delta_tech = strain_te-strain_ch;
//		int count = 0 ;
//		while(std::abs(upTestVal-downTestVal) > 1e-8*tensionCritStrain && count++  <  32 )
//		{
//			double testVal = (upTestVal+downTestVal)*.5 ;
//			double mainCurve = 1./(1.+sqrt(k)*pow(testVal-tensionCritStrain,0.8)) ;
//
//			if(testVal < tensionCritStrain)
//				factor = 1. ;
//			else if(testVal < strain_ch)
//				factor = mainCurve ;
//			else if(testVal < strain_te)
//				factor = mainCurve*(strain_te-testVal)/delta_tech ;
//			else
//				factor = 0 ;
//			
//			
//			if( testVal*pseudoYoung > upVal*factor+interactionStress )
//				upTestVal = testVal ;
//			else
//				downTestVal = testVal ;
//			
//		}
//
//		maxTensionStrain = (upTestVal+downTestVal)*.5 ;
//		maxTension = (upTestVal+downTestVal)*.5*scaleFactor*pseudoYoung ;
//	}
//	
	if(factors.size() == 0)
		initialiseFactors(s) ;
	double nlfactor = factors[0]/std::accumulate(&factors[0], &factors[factors.size()], double(0)) ;
	
	double corrfact = 1 ;
// 	for(size_t i = 1 ; i < 24 ; i++)
// 	{
// 		double a = 1./nlfactor * (sqrt(value*strain_ch)-log(sqrt(value*strain_ch+1)))+log(corrfact*value*strain_ch) ;
// 		a *= a ;
// 		a /= corrfact*strain_ch ;
// 		corrfact = a ;
// // 		std::cout << corrfact << std::endl ;
// 	}
	
// 	exit(0) ;
	
 	double maxTension = upVal*scaleFactor;
 	double maxTensionStrain = tensionCritStrain;

//  	value /= factors[0]/std::accumulate(&factors[0], &factors[factors.size()], double(0)) ; //getMaterialCharacteristicRadius()*getMaterialCharacteristicRadius()*M_PI ;
 	if(tstrain > epsilon_0 )
 	{
 		double downTestVal = epsilon_0 ;
 		double upTestVal = strain_te*3. ;
 		double factor = 1 ;
 		double delta_tech = strain_te*3-strain_ch*3.;
 		double mainCurve = 0 ;
 		int count = 0 ;
 		while(std::abs(upTestVal-downTestVal) > 1e-8*tensionCritStrain && count++  <  32)
 		{
 			double testVal = (upTestVal+downTestVal)*.5 ;
 			mainCurve = 1./(1.+sqrt(value*corrfact)*sqrt(testVal-epsilon_0*.5)) ;
 
 			if(testVal < strain_ch*3.)
  			factor = mainCurve ;
 			else if(testVal < strain_te*3.)
 			{
 				factor = mainCurve*(strain_te*3.-testVal)/delta_tech ;
 // 				factor = std::min((strain_te*3.-testVal)/(strain_te*3.-tensionCritStrain), mainCurve) ;
 			}
 			else
 				factor = 0 ;
 			
 
 			if( testVal*pseudoYoung > upVal*factor )
 				upTestVal = testVal ;
 			else
 				downTestVal = testVal ;
 
 		}
 	
 		maxTension = (upTestVal+downTestVal)*.5*scaleFactor*pseudoYoung ;
 		maxTensionStrain = (upTestVal+downTestVal)*.5 ;
 
 		if(factor < POINT_TOLERANCE_2D)
 			return POINT_TOLERANCE_2D ;
 	}
 	else
		return getBareConcreteTensileCriterion( s, pseudoYoung, tstrain, tstress) ;
	
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
	if(maxTensionStrain > POINT_TOLERANCE_2D)
		std::abs(tstress/(maxTensionStrain*pseudoYoung)) ;
	
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
	return strainc ; //strainc ;
	
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
// 	double df = f ;
	double df = 3.*f*f-2.*f*f*f ;
//
// 	return barecrit ;
	return df*barecrit + (1.-df)*rebcrit ;
}

double NonLocalMCFT::getConcreteCompressiveCriterion(const ElementState & s, double pseudoYoung, double cstrain, double tstress, double cstress)
{
	if(cstress >= 0 || cstrain >= 0)
		return -1 ;
// 	if(dynamic_cast<RotatingCrack *>(s.getParent()->getBehaviour()->getDamageModel()))
// 	{
// 		if(dynamic_cast<RotatingCrack *>(s.getParent()->getBehaviour()->getDamageModel())->compressionFailure)
// 			return POINT_TOLERANCE_2D ;
// 	}
	double maxCompression = downVal*scaleFactor  ;
	

	
	double maxCompressionStrain = downVal/pseudoYoung  ;
	
	if(true )
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
	
	if(std::abs(maxCompressionStrain) > POINT_TOLERANCE_2D)
		criterion = std::abs(cstress/(maxCompressionStrain*pseudoYoung)) ;
	
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
	return strainc ;
	return POINT_TOLERANCE_2D ;
}

void NonLocalMCFT::initialise( ElementState &s)
{
	double energy = 36000 ; //N/m 32000 prev
	strain_ch = 2.*energy/(8.*getMaterialCharacteristicRadius()*upVal) ;//*.5 energy <- // *2 energy -> 2.*energy/(1.*getMaterialCharacteristicRadius()*upVal) ;
	if(factors.size()==0)
		initialiseFactors(s) ;
	double nlcorrection = factors[0]/std::accumulate(&factors[0], &factors[factors.size()], double(0)) ;
// 	energy *= nlcorrection ;
	if(strain_ch < tensionCritStrain)
	{
		std::cout << strain_ch << " vs " << tensionCritStrain <<std::endl ;
		exit(0) ;
	}
	strain_te = 5.*strain_ch;
	double del_0 = (strain_ch-tensionCritStrain) ;
	double del_1 = (strain_te-strain_ch) ;
	double del = (strain_te-tensionCritStrain) ;
	
	double k0 = 1 ;
	double elastic_energy = tensionCritStrain*upVal*.5*nlcorrection ;
	double integral = elastic_energy ;
	do
	{
		integral = elastic_energy ;
		k0 *= 10;
		      
		for(double i = 0 ; i < 1e4 ; i++)
		{
// 			integral+= 0.5*(upVal/(1.+sqrt(k*(i)*1e-4*del_0)) + upVal/(1.+sqrt(k*(i+1.)*1e-4*del_0)))*del_0*1e-4 ;
			integral+= 0.5*(upVal/(1.+sqrt(k0)*pow((i)*1e-4*del,0.8)) + upVal/(1.+sqrt(k0)*pow((i+1.)*1e-4*del,0.8)))*del*1e-4 ;
		}
		      
		for(double i = 0 ; i < 10000 ; i++)
		{  
			integral+= 0.5*((del_1-(i)*1e-4*del_1)/del_1)*(upVal/(1.+sqrt(k0)*pow((i)*1e-4*del_1+del_0,0.8))+upVal/(1.+sqrt(k0)*pow((i+1.)*1e-4*del_1+del_0,0.8)))*del_1*1e-4 ;
		}
	}
	while(integral > energy) ;
	
	
	double k_low = k0/100.  ;
	double k_high = k0 ;

	int count = 0 ;
	
	do
	{
		integral = elastic_energy ;
		k = 0.5*(k_low+k_high) ;
		      
		for(double i = 0 ; i < 1e4 ; i++)
		{
// 			integral+= 0.5*(upVal/(1.+sqrt(k*(i)*1e-4*del_0)) + upVal/(1.+sqrt(k*(i+1.)*1e-4*del_0)))*del_0*1e-4 ;
			integral+= 0.5*(upVal/(1.+sqrt(k)*pow((i)*1e-4*del,0.8)) + upVal/(1.+sqrt(k)*pow((i+1.)*1e-4*del,0.8)))*del*1e-4 ;
		}
		      
		for(double i = 0 ; i < 10000 ; i++)
		{  
			integral+= 0.5*((del_1-(i)*1e-4*del_1)/del_1)*(upVal/(1.+sqrt(k)*pow((i)*1e-4*del_1+del_0,0.8))+upVal/(1.+sqrt(k)*pow((i+1.)*1e-4*del_1+del_0,0.8)))*del_1*1e-4 ;
		}
		if(integral < energy)
			k_high = k ;
		else
			k_low = k ;
		count++ ;
		
	} while(std::abs(integral-energy) > 1e-4*energy && count < 32) ;
	
// 	k*= getMaterialCharacteristicRadius()*getMaterialCharacteristicRadius()*M_PI ;
// 		k/=1e5 ;
		
	if(std::abs(integral-energy) > 1e-4*energy)
	{
	  std::cout << "wrong energy! " << integral << " vs "<< energy  << " k = "<< k << std::endl ;
	  exit(0) ;
	}
	initialised = true ;
}


// rough crack model by bazant and gambarova
// From engineering fracture mechanics, smeared crack analysis for fracture and aggregate 
// interlock in concrete, Gambarova and Valente
std::pair<double, double> stressOnCrack(ElementState &s, double downVal, double dmax = 16)
{
	double tau_0 = 0.25*downVal*1e-6 ;
	double a12 = 0.62 ;
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
	
	return std::make_pair(sigma_nn*1e6, sigma_nt*1e6) ;
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
	double tstress = first.max()-std::max(std::abs(crackstress.first)*1e6, 0.);   //first.max();//std::min(first.max()-2.*std::max(crackstress.first, 0.), first.max());// 
	double cstress = first.min()+std::min(std::abs(crackstress.second)*1e6, 0.);   //first.min();//std::max(first.min()+2.*std::min(crackstress.second, 0.), first.min());// 

	
	
	
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
// 		if(c0 > c1)
// 		{
// 			firstMet = true ;
// 			secondMet = false ;
// 		}
// 		else
// 		{
// 			firstMet = false ;
// 			secondMet = true ;
// 		}
		
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
