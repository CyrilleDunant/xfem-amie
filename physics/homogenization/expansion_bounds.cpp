//
// C++ Implementation: mechanical analytic homogenization
//
// Description:
//
//
// Author:  Alain Giorla, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "expansion_bounds.h"

namespace Mu
{


ExpansionBoundsScheme::ExpansionBoundsScheme(int i) : ExpansionHomogenizationScheme(i)
{
	output.push_back(TAG_EXPANSION_COEFFICIENT) ;
}

ShaperyScheme::ShaperyScheme() : ExpansionBoundsScheme(3)
{
}

Vector ShaperyScheme::process(const Matrix & data)
{
	double fmat = data[0][0] ;
	double kmat = data[0][1] ;
	double mumat= data[0][2] ;
	double amat = data[0][3] ;
	
	double finc = data[1][0] ;
	double kinc = data[1][1] ;
	double muinc= data[1][2] ;
	double ainc = data[1][3] ;

	double khom = data[2][1] ;

	double common = amat*fmat + ainc*finc ;
	double num = (khom-kmat)*(amat-ainc)/khom ;


	Vector processed(2) ;
	processed[0] = common + mumat*finc*num/(4*mumat+3*kinc) ;
	processed[1] = common - muinc*fmat*num/(4*muinc+3*kmat) ;

	return processed ;
}

RosenHashinScheme::RosenHashinScheme() : ExpansionBoundsScheme(4)
{
}

Vector RosenHashinScheme::process(const Matrix & data)
{
	double fmat = data[0][0] ;
	double kmat = data[0][1] ;
	double mumat= data[0][2] ;
	double amat = data[0][3] ;
	
	double finc = data[1][0] ;
	double kinc = data[1][1] ;
	double muinc= data[1][2] ;
	double ainc = data[1][3] ;

	double kup  = data[2][1] ;
	double klow = data[3][1] ;

	double common = amat*fmat + ainc*finc ;
	double num = fmat*finc*4*(kmat-kinc)*(amat-ainc) ;
	double denom = kmat*kinc*3 ;

	Vector processed(2) ;
	processed[0] = common + muinc*num/(denom + 4*muinc*kup) ;
	processed[1] = common + mumat*num/(denom + 4*mumat*klow) ;

	return processed ;
}










}


