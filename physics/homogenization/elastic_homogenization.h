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

#ifndef ELASTIC_SCHEME
#define ELASTIC_SCHEME

#include "scheme_base.h"
#include "../../geometry/geometry_base.h" 

namespace Mu
{

/**
* \brief builds a CauchyGreen tensor from the YoungModulus/PoissonRatio (or BulkLModulus/ShearModulus), in 2D or in 3D
*/
Matrix cauchyGreen(std::pair<double,double> prop, bool hooke, SpaceDimensionality dim) ;

/**
* Abstract elastic scheme (returns the composite bulk (index 0) and shear moduli (index 1).
* When an homogenization scheme represents a matrix-inclusion morphology, the matrix is always 
* the first material detected.
*/
class ElasticHomogenizationScheme : public Scheme
{
public:
	ElasticHomogenizationScheme(int i) ;
} ;

/**
* Diluted scheme. This scheme is only valid for small fractions of inclusions
*/
class DilutedScheme : public ElasticHomogenizationScheme
{
public:
	DilutedScheme() ;
	virtual Vector process(const Matrix & data) ;

} ;

/**
* Generalization of the diluted scheme to any number of inclusions.
* This scheme is only valid for small fractions of inclusions
*/
class GeneralizedDilutedScheme : public ElasticHomogenizationScheme
{
public:
	GeneralizedDilutedScheme() ;
	virtual Vector process(const Matrix & data) ;
} ;


/** 
* Standard bi-phasic incremental scheme 
*/
class IncrementalScheme : public ElasticHomogenizationScheme
{
protected:
	/**
	* \brief fraction increment
	*/
	double dalpha ;
public:
	IncrementalScheme(double d) ;
	virtual Vector process(const Matrix & data) ;
} ;



/**
* standard bi-phasic Mori-Tanaka elastic homogenization scheme 
*/
class MoriTanaka : public ElasticHomogenizationScheme
{
public:
	MoriTanaka() ;
	virtual Vector process(const Matrix & data) ;

} ;

/**
* Mori-Tanaka elastic homogenization scheme generalized to n-phases material
*/
class GeneralizedMoriTanaka : public ElasticHomogenizationScheme
{
public:
	GeneralizedMoriTanaka() ;
	virtual Vector process(const Matrix & data) ;

} ;


/**
* Self-Consistent elastic homogenization scheme for bi-phasic material 
*/
class SelfConsistent : public ElasticHomogenizationScheme
{
public:
	SelfConsistent() ;
	virtual Vector process(const Matrix & data) ;

} ;

/**
* Self-Consistent elastic homogenization scheme generalized for n-phases material
*/
class GeneralizedSelfConsistent : public ElasticHomogenizationScheme
{
public:
	GeneralizedSelfConsistent() ;
	virtual Vector process(const Matrix & data) ;


} ;



} ;

#endif
