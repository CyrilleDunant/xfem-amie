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

#ifndef CRACKED_HOMOGENIZATION
#define CRACKED_HOMOGENIZATION

#include "scheme_base.h"

namespace Mu
{

/* \brief void expansion homogenization scheme (returns the matrix expansion coefficient)*/
class CrackedHomogenizationScheme : public HomogenizationScheme
{
public:
	/* \brief constructor 
	* @param i the number of phases (-1 for infinite)
	*/
	CrackedHomogenizationScheme() ;
} ;

/* \brief Diluted scheme. This scheme is only valid for small fractions */
class BudianskyScheme : public CrackedHomogenizationScheme
{
public:
	BudianskyScheme() ;
	virtual Vector processData(const Matrix & data) ;

} ;

class SimplifiedBenHahaScheme : public CrackedHomogenizationScheme
{
public:
	SimplifiedBenHahaScheme() ;
	virtual Vector processData(const Matrix & data) ;

} ;



} ;

#endif
