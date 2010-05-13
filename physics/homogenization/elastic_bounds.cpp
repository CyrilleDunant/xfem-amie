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

#include "elastic_bounds.h"

namespace Mu
{






ElasticBoundsScheme::ElasticBoundsScheme(int i) : Scheme(i)
{
	input.push_back(TAG_VOLUME_FRACTION) ;
	input.push_back(TAG_BULK_MODULUS) ;
	input.push_back(TAG_SHEAR_MODULUS) ;

	output.push_back(TAG_BULK_MODULUS) ;
	output.push_back(TAG_SHEAR_MODULUS) ;
	output.push_back(TAG_BULK_MODULUS) ;
	output.push_back(TAG_SHEAR_MODULUS) ;
}



HillBounds::HillBounds() : ElasticBoundsScheme(2)
{
}

Vector HillBounds::process(const Matrix & data)
{
	std::vector<Tag> avg ;
	avg.push_back(TAG_BULK_MODULUS) ;
	avg.push_back(TAG_SHEAR_MODULUS) ;
	MeanScheme p(true, true, avg) ;
	MeanScheme s(true, false, avg) ;
	Vector upper = p.process(data) ;
	Vector lower = s.process(data) ;

	Vector processed(data.numCols()*2-2) ;
	for(size_t i = 0 ; i < upper.size() ; i++)
	{
		processed[i] = upper[i] ;
		processed[upper.size()+i] = lower[i] ;
	}

	return processed ;
}






}


