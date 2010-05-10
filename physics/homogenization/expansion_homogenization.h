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

#ifndef EXPANSION_HOMOGENIZATION
#define EXPANSION_HOMOGENIZATION

#include "scheme_base.h"

namespace Mu
{

/* \brief void expansion homogenization scheme (returns the matrix expansion coefficient)*/
class ExpansionHomogenizationScheme : public Scheme
{
public:
	/* \brief constructor 
	* @param i the number of phases (-1 for infinite)
	*/
	ExpansionHomogenizationScheme(size_t i) ;
} ;

/* \brief Diluted scheme. This scheme is only valid for small fractions */
class HashinScheme : public ExpansionHomogenizationScheme
{
public:
	HashinScheme() ;
	virtual Vector process(const Matrix & data) ;

} ;



} ;

#endif
