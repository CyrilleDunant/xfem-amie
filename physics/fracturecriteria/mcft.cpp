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
#include "../../mesher/delaunay.h"
#include "../../mesher/delaunay_3d.h"
#include <set>

namespace Mu
{

MCFT::MCFT( double up, double down, double youngModulus, double charRad, MirrorState mirroring, double delta_x, double delta_y, double delta_z ) : FractureCriterion( mirroring, delta_x, delta_y, delta_z )
	, upVal( up ), downVal( down )
{
	physicalCharacteristicRadius = charRad ;
	tensionCritStrain = up / youngModulus ;
}


MCFT::~MCFT()
{
}

double MCFT::grade( ElementState &s )
{
	Vector pstrain = s.getPrincipalStrains( s.getParent()->getCenter() ) ;
	Vector pstress = s.getPrincipalStresses( s.getParent()->getCenter() ) ;

	double tstrain = pstrain.max();
	double cstrain = pstrain.min();
	double tstress = pstress.max();
	double cstress = pstress.min();

	double critStrain = -0.002 ;
	double renormCompressionStrain = cstrain / critStrain ;
	double compressionFactor = 1 ;
	double maxCompression = -std::abs( downVal ) ;
// 	if(cstrain < critStrain)
// 	{
		compressionFactor = ( 2.*renormCompressionStrain - renormCompressionStrain * renormCompressionStrain ) / ( 0.8 - 0.34 * tstrain / critStrain ) ;
		if( compressionFactor > 1 || compressionFactor < 0 )
			compressionFactor = 1 ;
		maxCompression = -std::abs( downVal ) * compressionFactor ;
// 	}

	

	double maxTension = upVal ;

	if(tstrain > tensionCritStrain)
	{
	//Yamamoto model
//		maxTension = upVal/(1.+sqrt(2e6*(tstrain+tensionCritStrain))) ;

	//MCFT model
	maxTension = upVal / ( 1. + sqrt( 200.*tstrain ) ) ;

	//perfectly brittle
// 		maxTension = 0 ;
	}



	metInCompression = cstrain <= 0 && std::abs( cstress / maxCompression ) > std::abs( tstress / maxTension ) || tstrain <= 0;
	metInTension = tstrain >= 0 && std::abs( cstress / maxCompression ) < std::abs( tstress / maxTension ) || cstrain >= 0;


	std::vector<double> crits ;
	crits.push_back( -1 ) ;

	if( cstress <= 0 && std::abs( cstress ) >= std::abs( maxCompression ) )
	{
		metInCompression = true ;
		crits.push_back( 1. - std::abs( maxCompression / cstress ) ) ;
	}

	if( tstress >= 0 && std::abs( tstress ) >= std::abs( maxTension ) )
	{
		metInTension = true ;
		crits.push_back( 1. - std::abs( maxTension / tstress ) ) ;
	}


	if( tstress >= 0 && std::abs( tstress )  < std::abs( maxTension ) )
	{
		crits.push_back( -1. + std::abs( tstress / maxTension ) ) ;
	}

	if( cstress <= 0 && std::abs( cstress ) < std::abs( downVal ) )
	{
		crits.push_back( -1. + std::abs( cstress / downVal ) ) ;
	}

	std::sort( crits.begin(), crits.end() );
	return crits.back() ;
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


NonLocalMCFT::NonLocalMCFT( double up, double down, double youngModulus,  double charRad, MirrorState mirroring, double delta_x, double delta_y, double delta_z ) : FractureCriterion( mirroring, delta_x, delta_y, delta_z )
	, upVal( up ), downVal( down ), youngModulus(youngModulus)
{
	physicalCharacteristicRadius = charRad ;
	critStrain = -0.0015;
	tensionCritStrain = upVal / youngModulus ;
	strainBroken = false ;
	initialised = false ;
}


NonLocalMCFT::~NonLocalMCFT()
{
}

double NonLocalMCFT::grade( ElementState &s )
{
	if(!initialised)
	{
		double energy = 75. ; //N/m
		strain_ch = 2.*energy/(2.*getMaterialCharacteristicRadius()*upVal) ;
		
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
	
	std::pair<Vector, Vector> stressStrain = smoothedPrincipalStressAndStrain(s) ;

	double tstrain = stressStrain.second.max();
	double cstrain = stressStrain.second.min();
	double tstress = stressStrain.first.max();
	double cstress = stressStrain.first.min();

	double pseudoYoung = youngModulus*(1.-s.getParent()->getBehaviour()->getDamageModel()->getState().max()) ;

	double maxCompression = downVal  ;

	if(cstrain < critStrain*.5 )
	{
// 		strainBroken = true ;
		double C_d = 0. ;
		double compressiveTensileRatio = -std::abs(tstress/std::min(cstress, -POINT_TOLERANCE_2D)) ;
		
		if(-compressiveTensileRatio > 0.2801)
			C_d =0.35*pow(-compressiveTensileRatio-0.28, 0.8) ;
		
		double beta_d = 1./(1.+C_d) ;
		
		if(beta_d > 1)
			beta_d = 1 ;
		
		double f_p = beta_d*downVal ;
		double epsilon_p = beta_d*critStrain ;
		double n = 0.8 - f_p/17e6 ;
		double k_c = 0.67 - f_p/62e6 ;
		
		double upTestVal = 0 ;
		double downTestVal = 0 ;
		bool below = true ;
		int ninc = 0 ;
		while(below && ninc < 1000)
		{
			ninc++ ;
			downTestVal = upTestVal ;
			upTestVal += 1e-3* downVal ;
			if(upTestVal/pseudoYoung > epsilon_p)
				k_c = 1. ;
			else
				k_c = 0.67 - f_p/62e6 ;
			double rcs = upTestVal/(epsilon_p*pseudoYoung) ;
			double factor = n*rcs/(n-1.+pow(rcs, n*k_c)); 
			if(f_p*factor > upTestVal)
			{
				below = false ;
			}
		}
		if(ninc == 1)
		{
			upTestVal = cstress ;
			downTestVal = cstress ;
			if(upTestVal > epsilon_p/pseudoYoung)
				k_c = 1. ;
			else
				k_c = 0.67 - f_p/62e6 ;
			double rcs = upTestVal/(epsilon_p*pseudoYoung) ;
			double factor = n*rcs/(n-1.+pow(rcs, n*k_c)); 
			upTestVal = downTestVal = f_p*factor ;
		}
		if(ninc >= 1000)
		{
// 			return 1 ;

// 			upTestVal = downVal ;
// 			downTestVal = downVal ;

			upTestVal = cstress ;
			downTestVal = cstress ;
			if(upTestVal > epsilon_p/pseudoYoung)
				k_c = 1. ;
			else
				k_c = 0.67 - f_p/62e6 ;
			double rcs = upTestVal/(epsilon_p*pseudoYoung) ;
			double factor = n*rcs/(n-1.+pow(rcs, n*k_c)); 
			upTestVal = downTestVal = f_p*factor ;
			
			return 1.- std::min(cstress/upTestVal, upTestVal/cstress) ;
		}
		
		while(std::abs(upTestVal-downTestVal) > -1e-5*downVal )
		{
			double testVal = (upTestVal+downTestVal)*.5 ;

			if(testVal > epsilon_p/pseudoYoung)
				k_c = 1. ;
			else
				k_c = 0.67 - f_p/62e6 ;
			double rcs = testVal/(epsilon_p*pseudoYoung);
			double factor = n*rcs/(n-1.+pow(rcs, n*k_c)); 
			if(f_p*factor - testVal > 0)
			{
				upTestVal = testVal ;
			}
			else
			{
				downTestVal = testVal ;
			}
		}
		maxCompression = std::min((upTestVal+downTestVal)*.5, -0.) ;

	}
	
	double maxTension = upVal;
	std::vector<double> crits ;
		
	if(tstrain > tensionCritStrain )
	{
// 		std::cout << 0 << "  " << 0 << std::endl ;
// 		for(double i = 0 ; i < 1 ; i += 0.001)
// 		{
// 		pseudoYoung = youngModulus*(1.-i) ;
		double downTestVal = tensionCritStrain*pseudoYoung ;
		double upTestVal = upVal ;
		double factor = 1 ;
		while(std::abs(upTestVal-downTestVal) > 1e-6*upVal )
		{
			
			double testVal = (upTestVal+downTestVal)*.5 ;
			
			if(testVal < strain_ch*pseudoYoung)
				factor = 1./(1.+sqrt(k*(testVal/pseudoYoung-tensionCritStrain))) ;
			else if(testVal < strain_te*pseudoYoung)
				factor = 1./(1.+sqrt(k*(testVal/pseudoYoung-tensionCritStrain)))*(strain_te-testVal/pseudoYoung)/(strain_te-strain_ch) ;
			else
				factor = 0 ;
		
			if( testVal > upVal*factor )
			{
				upTestVal = testVal ;
			}
			else
			{
				downTestVal = testVal ;
			}
		}
// 		std::cout << (upTestVal+downTestVal)*.5/pseudoYoung << "  " << factor << std::endl ;
		maxTension = (upTestVal+downTestVal)*.5 ;
		
// 		}
// 		exit(0) ;
// 		if(factor < POINT_TOLERANCE_2D)
// 			return 1.- strain_te/tstrain ;

	}



metInCompression = cstress < -POINT_TOLERANCE_2D  || cstrain < -POINT_TOLERANCE_2D;
metInTension = tstrain > POINT_TOLERANCE_2D;

crits.push_back(0.) ;
if(std::abs(maxCompression) > POINT_TOLERANCE_2D)
	crits.push_back(cstress/maxCompression) ;
else
	crits.push_back(1e6) ;
// if(std::abs(maxTension) > POINT_TOLERANCE_2D)
// 	crits.push_back(tstress/maxTension) ;
// else
// 	crits.push_back(1e6) ;

std::sort( crits.begin(), crits.end() );

if(crits.back() > 1)
{
	return 1.-1./crits.back() ;
}
if(crits.back() > -POINT_TOLERANCE_2D)
{
	return -1.+crits.back() ;
}
if(crits.back() > -1)
{
	return crits.back() ;
}

return -1 ;

}

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
