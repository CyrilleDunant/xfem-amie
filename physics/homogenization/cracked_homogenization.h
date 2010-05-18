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

/**
* These schemes computes the mechanical properties of a single cracked material, according to the crack density
*/
class CrackedHomogenizationScheme : public Scheme
{
public:
	CrackedHomogenizationScheme() ;
} ;

/**
* This simple scheme uses an empyrical curve for the evolution of damage as a function of expansion
*/
class ExpansionDefinedCrackScheme : public CrackedHomogenizationScheme
{
protected:
	double threshold ;
	double k ;
	double m ;
public:
	ExpansionDefinedCrackScheme(double t, double k, double m) ;
	virtual Vector process(const Matrix & data) ;
} ;

class AbouChakraCrackScheme : public CrackedHomogenizationScheme
{
public:
	AbouChakraCrackScheme() ;
	virtual Vector process(const Matrix & data) ;
} ;




/* \brief Budiansky dry crack scheme */
class BudianskyDryCrackScheme : public CrackedHomogenizationScheme
{
public:
	BudianskyDryCrackScheme() ;
	virtual Vector process(const Matrix & data) ;

} ;

/**
* Simplified Budiansky scheme as defined in Ben Haha thesis
*/
class SimplifiedBenHahaScheme : public CrackedHomogenizationScheme
{
public:
	SimplifiedBenHahaScheme() ;
	virtual Vector process(const Matrix & data) ;

} ;



} ;

#endif
