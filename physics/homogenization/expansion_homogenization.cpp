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

#include "expansion_homogenization.h"

namespace Mu
{






ExpansionHomogenizationScheme::ExpansionHomogenizationScheme(size_t i) : Scheme(i)
{
	input.push_back(TAG_VOLUME_FRACTION) ;
	input.push_back(TAG_BULK_MODULUS) ;
	input.push_back(TAG_SHEAR_MODULUS) ;
	input.push_back(TAG_EXPANSION_COEFFICIENT) ;

	output.push_back(TAG_EXPANSION_COEFFICIENT) ;
}



HashinScheme::HashinScheme() : ExpansionHomogenizationScheme(2)
{
}

Vector HashinScheme::processData(const Matrix & data)
{
	double kmat = data[0][1] ;
	double amat = data[0][3] ;
	
	double finc = data[0][0] ;
	double kinc = data[0][1] ;
	double ainc = data[0][3] ;
	
	Vector processed(1) ;
	processed[0] = amat + 2*finc*kinc*(ainc-amat) / (kmat+kinc+finc*(kinc-kmat)) ;
	
	return processed ;
}






}


