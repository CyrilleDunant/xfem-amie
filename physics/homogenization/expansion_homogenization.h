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

/** 
* Asbtract expansion homogenization scheme. These schemes return the homogenized expansion
* coefficient of the composite material. For matric-inclusion morphologies, the matrix should
* always be the first input material
*/
class ExpansionHomogenizationScheme : public Scheme
{
public:
	/* \brief constructor 
	* @param i the number of phases (-1 for infinite)
	*/
	ExpansionHomogenizationScheme(int i) ;
} ;

/**
* Turner scheme for bi-phasic materials. 
*/
class TurnerScheme : public ExpansionHomogenizationScheme
{
public:
	TurnerScheme() ;
	virtual Vector process(const Matrix & data) ;

} ;

/**
* Kerner scheme for bi-phasic materials.
*/
class KernerScheme : public ExpansionHomogenizationScheme
{
public:
	KernerScheme() ;
	virtual Vector process(const Matrix & data) ;

} ;

/**
* Hoobs scheme for bi-phasic materials (see Ben Haha thesis)
*/
class HobbsScheme : public ExpansionHomogenizationScheme
{
public:
	HobbsScheme() ;
	virtual Vector process(const Matrix & data) ;

} ;

/**
* This scheme takes different expansion homogenization schemes and makes the mean of the results
* of each scheme, using a different ponderation for each scheme
*/
class HirschScheme : public ExpansionHomogenizationScheme
{
public:
	std::vector<ExpansionHomogenizationScheme> expansion ;
	std::vector<double> probability ;

	HirschScheme() ;
	void addScheme(ExpansionHomogenizationScheme exp, double p) ;
	virtual Vector process(const Matrix & sata) ;
} ;

/**
* Self-consistent scheme using the three-layers sphere approach, the third layer being the
* homogenized material. This scheme uses 3 phases: the matrix, the inclusion, and the homogenized
* material.
*/
class SelfConsistentSiderisScheme : public ExpansionHomogenizationScheme
{
public:
	SelfConsistentSiderisScheme() ;
	virtual Vector process(const Matrix & data) ;
} ;

/**
* This scheme returns the strain in the inclusion and in the matrix for a bi-phasic composite.
* It does not give an equivalent expansion coefficient.
*/
class SantScheme : public ExpansionHomogenizationScheme
{
public:
	SantScheme() ;
	virtual Vector process(const Matrix & data) ;
} ;

class MultiLayerExpansion : public ExpansionHomogenizationScheme
{
protected:
	double epsilon ;
public:
	MultiLayerExpansion(double eps) ;
	virtual Vector process(const Matrix & data) ;
} ;

class SelfConsistentMultiLayerExpansion : public MultiLayerExpansion
{
public:
	SelfConsistentMultiLayerExpansion(double eps) ;
	virtual Vector process(const Matrix & data) ;
} ;



} ;

#endif
