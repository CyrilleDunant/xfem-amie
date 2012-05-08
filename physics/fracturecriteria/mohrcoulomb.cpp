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
#include "mohrcoulomb.h"
#include "../../mesher/delaunay.h"
#include "../../mesher/delaunay_3d.h"
#include "../damagemodels/damagemodel.h"

namespace Mu
{

MohrCoulomb::MohrCoulomb( double up, double down, MirrorState mirroring, double delta_x, double delta_y, double delta_z ) : FractureCriterion( mirroring, delta_x, delta_y, delta_z )
	, upVal( up ), downVal( down )
{
	metInTension = false ;
	metInCompression = false ;
}


MohrCoulomb::~MohrCoulomb()
{
}

double MohrCoulomb::grade( ElementState &s )
{

	if( s.getParent()->getBehaviour()->fractured() )
		return 0 ;


	Vector pstress = s.getPrincipalStresses(s.getParent()->getCenter() ) ;

	double maxStress = pstress.max() ;
	double minStress = pstress.min() ;

// 	std::cout << pstress0[0] << ", " << pstress0[1] << ", "<< pstress0[2] << std::endl ;
	metInTension = false ;
	metInCompression = false ;
	metInCompression = std::abs( minStress / downVal ) > std::abs( maxStress / upVal ) ;
	metInTension = std::abs( minStress / downVal ) < std::abs( maxStress / upVal ) ;

	if( maxStress >= upVal )
	{
		metInTension = true;

		if( minStress <= downVal )
			metInCompression = true ;

		return 1. - std::abs( upVal / maxStress ) ;
	}

	if( minStress <= downVal )
	{
		metInCompression = true ;
		return 1. - std::abs( downVal / minStress ) ;
	}

	double s0 = -1. + std::abs( maxStress / upVal );
	double s1 = -1. + std::abs( minStress / downVal ) ;

	if( minStress > 0 )
	{
		return s0 ;
	}

	if( maxStress < 0 )
	{
		return s1 ;
	}



	if( std::abs( s0 ) > std::abs( s1 ) )
		return s0 ;

	return s1;
}

FractureCriterion *MohrCoulomb::getCopy() const
{
	return new MohrCoulomb( *this ) ;
}

Material MohrCoulomb::toMaterial()
{
	Material mat ;
	return mat ;
}

NonLocalMohrCoulomb::NonLocalMohrCoulomb( double up, double down, double E, MirrorState mirroring, double delta_x, double delta_y, double delta_z ) : FractureCriterion( mirroring, delta_x, delta_y, delta_z )
	, upVal( up ), downVal( down ), stiffness(E)
{
	metInTension = false ;
	metInCompression = false ;
}


NonLocalMohrCoulomb::~NonLocalMohrCoulomb()
{

}

double NonLocalMohrCoulomb::grade( ElementState &s )
{

	if( s.getParent()->getBehaviour()->fractured() )
		return -1 ;

	std::pair<Vector, Vector> pstressStrain( smoothedPrincipalStressAndStrain(s)) ;
	Vector pstress = pstressStrain.first ;
	Vector pstrain = pstressStrain.second ;
	double maxStress = pstress.max() ;
	double minStress = pstress.min() ;
	double maxStrain = pstrain.max() ;
	double minStrain = pstrain.min() ;

// 	std::cout << pstress0[0] << ", " << pstress0[1] << ", "<< pstress0[2] << std::endl ;
	metInTension = false ;
	metInCompression = false ;
	metInCompression = std::abs( minStress / downVal ) > std::abs( maxStress / upVal ) ;
	metInTension = std::abs( minStress / downVal ) < std::abs( maxStress / upVal ) ;

	double effectiveStiffness = stiffness ;
	if(s.getParent()->getBehaviour()->getDamageModel())
		effectiveStiffness = stiffness*(1.-s.getParent()->getBehaviour()->getDamageModel()->getState().max()) ;
	
	double  upStrain = upVal/effectiveStiffness ;
	double  downStrain = downVal/effectiveStiffness ;
	std::vector<double> scores ;
	scores.push_back(-1);
	if( maxStrain >= upStrain && maxStrain > 0 )
	{
		metInTension = true;
		scores.push_back(1. - std::abs( upStrain / maxStrain ));
	}
	else if(maxStrain > 0)
			scores.push_back(-1. + std::abs( maxStrain / upStrain ));

	if( minStrain <= downStrain && minStrain < 0 )
	{
		metInCompression = true ;
		scores.push_back(1. - std::abs( downStrain / minStrain )) ;
	}
	else if(minStrain < 0 )
	{
		scores.push_back(-1. + std::abs( minStrain / downStrain )) ;
	}
	std::sort(scores.begin(), scores.end()) ;
	return scores.back() ;
}

FractureCriterion *NonLocalMohrCoulomb::getCopy() const
{
	return new NonLocalMohrCoulomb( *this ) ;
}

Material NonLocalMohrCoulomb::toMaterial()
{
	Material mat ;
	return mat ;
}



NonLocalLinearlyDecreasingMohrCoulomb::NonLocalLinearlyDecreasingMohrCoulomb( double up, double down,double limittstrain, double limitcstrain, double E, MirrorState mirroring, double delta_x, double delta_y, double delta_z ) : FractureCriterion( mirroring, delta_x, delta_y, delta_z )
	, upVal( up ), downVal( down ),limittstrain(limittstrain),limitcstrain(limitcstrain), stiffness(E)
{
	metInTension = false ;
	metInCompression = false ;
}


NonLocalLinearlyDecreasingMohrCoulomb::~NonLocalLinearlyDecreasingMohrCoulomb()
{

}

double NonLocalLinearlyDecreasingMohrCoulomb::grade( ElementState &s )
{

	if( s.getParent()->getBehaviour()->fractured() )
		return -1 ;

	std::pair<Vector, Vector> pstressStrain( smoothedPrincipalStressAndStrain(s)) ;
	Vector pstress = pstressStrain.first ;
	Vector pstrain = pstressStrain.second ;
	double maxStress = pstress.max() ;
	double minStress = pstress.min() ;
	double maxStrain = pstrain.max() ;
	double minStrain = pstrain.min() ;

// 	std::cout << pstress0[0] << ", " << pstress0[1] << ", "<< pstress0[2] << std::endl ;
	metInTension = false ;
	metInCompression = false ;
	metInCompression = std::abs( minStress / downVal ) > std::abs( maxStress / upVal ) ;
	metInTension = std::abs( minStress / downVal ) < std::abs( maxStress / upVal ) ;

	double effectiveStiffness = stiffness ;
	if(s.getParent()->getBehaviour()->getDamageModel())
		effectiveStiffness = stiffness*(1.-s.getParent()->getBehaviour()->getDamageModel()->getState().max()) ;
	
	double tfactor = 1.-(maxStrain-upVal/stiffness)/(limittstrain-upVal/stiffness) ;
	if(maxStrain > limittstrain)
		tfactor = 0 ;
	else if(maxStrain <= upVal/stiffness)
		tfactor = 1 ;
	
	double cfactor = 1.-(-minStrain+downVal/stiffness)/(-limitcstrain+downVal/stiffness) ;
	if(minStrain < limitcstrain)
		cfactor = 0 ;
	else if(minStrain > downVal/stiffness)
		cfactor = 1 ;
	
	double  upStrain = tfactor*upVal/effectiveStiffness ;
	double  downStrain = cfactor*downVal/effectiveStiffness ;
	std::vector<double> scores ;
	scores.push_back(-1);
	if( maxStrain >= upStrain && maxStrain > 0 )
	{
		metInTension = true;
		scores.push_back(1. - std::abs( upStrain / maxStrain ));
	}
	else if(maxStrain > 0 && upStrain > POINT_TOLERANCE_2D)
		scores.push_back(-1. + std::abs( maxStrain / upStrain ));
	else if(maxStrain > 0)
		return 1 ;

	if( minStrain <= downStrain && minStrain < 0 )
	{
		metInCompression = true ;
		scores.push_back(1. - std::abs( downStrain / minStrain )) ;
	}
	else if(minStrain < 0  && downStrain < -POINT_TOLERANCE_2D)
		scores.push_back(-1. + std::abs( minStrain / downStrain )) ;
	else if(minStrain < 0)
		return 1 ;
	
	std::sort(scores.begin(), scores.end()) ;

	return scores.back() ;
}

FractureCriterion *NonLocalLinearlyDecreasingMohrCoulomb::getCopy() const
{
	return new NonLocalLinearlyDecreasingMohrCoulomb( *this ) ;
}

Material NonLocalLinearlyDecreasingMohrCoulomb::toMaterial()
{
	Material mat ;
	return mat ;
}




}
