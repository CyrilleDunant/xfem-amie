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

#ifndef EXPANSION_BOUNDS
#define EXPANSION_BOUNDS

#include "expansion_homogenization.h"
#include "elastic_bounds.h"

namespace Mu
{

/**
* The expansion bounds scheme give the upper and lower bounds for the expansion coefficient
* of a composite material
*/
class ExpansionBoundsScheme : public ExpansionHomogenizationScheme
{
public:
	ExpansionBoundsScheme(int i) ;

} ;


/**
* The Shapery scheme represents a 2-phasic material. However, a third input is needed, which
* is the lower bound of the bulk modulus of the composite
*/
class ShaperyScheme : public ExpansionBoundsScheme
{
public:
	ShaperyScheme() ;
	virtual Vector process(const Matrix & data) ;
} ;

/**
* The Shapery scheme represents a 2-phasic material. However, two more input are needed, 
* the upper and the lower bounds of the bulk modulus of the composite
*/
class RosenHashinScheme : public ExpansionBoundsScheme
{
public:
	RosenHashinScheme() ;
	virtual Vector process(const Matrix & data) ;
} ;


} ;

#endif
