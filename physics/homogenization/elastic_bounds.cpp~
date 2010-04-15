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






ElasticBoundsScheme::ElasticBoundsScheme(size_t i) : HomogenizationScheme(i, FRACTION, BULK_SHEAR)
{
	input.push_back(BULK_SHEAR) ;
	output.push_back(BULK_SHEAR) ;
}



HillBounds::HillBounds() : ElasticBoundsScheme(2)
{
}

Vector HillBounds::processData(const Matrix & data)
{
	Vector upper = MeanSeries(BULK_SHEAR).processData(data) ;
	Vector lower = MeanParallel(BULK_SHEAR).processData(data) ;

	Vector processed(data.numCols()*2-2) ;
	for(size_t i = 0 ; i < upper.size() ; i++)
	{
		processed[i] = upper[i] ;
		processed[upper.size()+i] = lower[i] ;
	}

	return processed ;
}






}


