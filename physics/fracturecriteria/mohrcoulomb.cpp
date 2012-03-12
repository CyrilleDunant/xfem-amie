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
	for( size_t i = 0 ; i < testPoints.size() ; i++ )
		delete testPoints[i] ;
}

double MohrCoulomb::grade( ElementState &s )
{

	if( s.getParent()->getBehaviour()->fractured() )
		return 0 ;

	if( testPoints.size() == 0 )
	{
		if( s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
		{
			testPoints.resize( 4 );
			testPoints[0] = new Point( 0, 0, 0 ) ;
			testPoints[1] = new Point( 1, 0, 0 ) ;
			testPoints[2] = new Point( 0, 1, 0 ) ;
			testPoints[3] = new Point( 0, 0, 1 ) ;
		}
		else
		{
			testPoints.resize( 3 );
			testPoints[0] = new Point( 0, 0, 0 ) ;
			testPoints[1] = new Point( 1, 0, 0 ) ;
			testPoints[2] = new Point( 0, 1, 0 ) ;
		}
	}

	Vector pstress = s.getPrincipalStresses( testPoints, true ) ;

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
	for( size_t i = 0 ; i < testPoints.size() ; i++ )
		delete testPoints[i] ;
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
		effectiveStiffness =stiffness*(1.-s.getParent()->getBehaviour()->getDamageModel()->getState().max()) ;
	std::vector<double> scores ;
	scores.push_back(-1);
	if( maxStress >= upVal && maxStress > 0 || maxStrain > upVal/effectiveStiffness && maxStrain > 0)
	{
		metInTension = true;
		scores.push_back(std::min(1. - std::abs( upVal / maxStress ), 1. - std::abs( upVal/effectiveStiffness / maxStrain ) ));
	}
	else if(maxStress > 0)
			scores.push_back(-1. + std::abs( maxStress / upVal ));

	if( minStress <= downVal && minStress < 0 || minStrain < downVal/effectiveStiffness && minStrain < 0)
	{
		metInCompression = true ;
		scores.push_back(std::min(1. - std::abs( downVal / minStress ), 1. - std::abs( downVal/effectiveStiffness / minStrain ))) ;
	}
	else if(minStress < 0 )
	{
		scores.push_back(-1. + std::abs( minStress / downVal )) ;
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

}
