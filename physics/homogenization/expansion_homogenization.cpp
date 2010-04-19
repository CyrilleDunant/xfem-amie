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






ExpansionHomogenizationScheme::ExpansionHomogenizationScheme(size_t i) : HomogenizationScheme(i, FRACTION, EXPANSION)
{
	input.push_back(BULK_SHEAR) ;
	input.push_back(EXPANSION) ;
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


