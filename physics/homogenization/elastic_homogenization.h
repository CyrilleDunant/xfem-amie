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

/* \brief builds a CauchyGreen tensor from the YoungModulus/PoissonRatio (or BulkLModulus/ShearModulus), in 2D or in 3D*/
Matrix cauchyGreen(std::pair<double,double> prop, bool hooke, SpaceDimensionality dim) ;

/* \brief void elastic scheme (returns the matrix properties)*/
class ElasticHomogenizationScheme : public Scheme
{
public:
	/* \brief constructor 
	* @param i the number of phases (-1 for infinite)
	*/
	ElasticHomogenizationScheme(int i) ;
} ;

/* \brief Diluted scheme. This scheme is only valid for small fractions */
class DilutedScheme : public ElasticHomogenizationScheme
{
public:
	DilutedScheme() ;
	virtual Vector process(const Matrix & data) ;

} ;

/* \brief Generalized diluted scheme. This scheme is only valid for small fractions of inclusions */
class GeneralizedDilutedScheme : public ElasticHomogenizationScheme
{
public:
	GeneralizedDilutedScheme() ;
	virtual Vector process(const Matrix & data) ;
} ;


/* \brief Standard bi-phasic incremental scheme */
class IncrementalScheme : public ElasticHomogenizationScheme
{
protected:
	double dalpha ;
public:
	IncrementalScheme(double d) ;
	virtual Vector process(const Matrix & data) ;
} ;



/* \brief standard bi-phasic Mori-Tanaka elastic homogenization scheme */
class MoriTanaka : public ElasticHomogenizationScheme
{
public:
	MoriTanaka() ;
	virtual Vector process(const Matrix & data) ;

} ;

/* \brief Mori-Tanaka elastic homogenization scheme generalized to n-phases material */
class GeneralizedMoriTanaka : public ElasticHomogenizationScheme
{
public:
	GeneralizedMoriTanaka() ;
	virtual Vector process(const Matrix & data) ;

} ;


/* \brief Self-Consistent elastic homogenization scheme for bi-phasic material */
class SelfConsistent : public ElasticHomogenizationScheme
{
public:
	SelfConsistent() ;
	virtual Vector process(const Matrix & data) ;

} ;

/* \brief Self-Consistent elastic homogenization scheme generalized for n-phases material */
class GeneralizedSelfConsistent : public ElasticHomogenizationScheme
{
public:
	GeneralizedSelfConsistent() ;
	virtual Vector process(const Matrix & data) ;


} ;



} ;

#endif
