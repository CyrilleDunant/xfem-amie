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

/* \brief void elastic scheme (returns the matrix properties)*/
class ElasticBoundsScheme : public HomogenizationScheme
{
public:
	/* \brief constructor 
	* @param i the number of phases (-1 for infinite)
	*/
	ElasticBoundsScheme(size_t i) ;
} ;

/* \brief Diluted scheme. This scheme is only valid for small fractions */
class HillBounds : public ElasticBoundsScheme
{
public:
	HillBounds() ;
	virtual Vector processData(const Matrix & data) ;

} ;



} ;

#endif
