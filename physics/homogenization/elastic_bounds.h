//
// C++ Interface: mechanical homogenization
//
// Description: Analytic elastic homogenization schemes
// In all the homogenization schemes, the first material is ALWAYS the matrix
//
//
// Author:  Alain Giorla, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef ELASTIC_BOUNDS_SCHEME
#define ELASTIC_BOUNDS_SCHEME

#include "scheme_base.h"
#include "elastic_homogenization.h"
#include "../../geometry/geometry_base.h" 

namespace Mu
{

/**
* An ElasticBoundsScheme is equivalent to an ElasticHomogenizationScheme, but it returns the
* upper bounds for the mechanical properties (index 0 and 1) and the lower bounds (index 2 and 3)
*/
class ElasticBoundsScheme : public Scheme
{
public:
	ElasticBoundsScheme(int i) ;
} ;

/**
* These bounds consist in a mean scheme in series and in parallel. They are
* wide bounds.
*/
class HillBounds : public ElasticBoundsScheme
{
public:
	HillBounds() ;
	virtual Vector process(const Matrix & data) ;

} ;



} ;

#endif
