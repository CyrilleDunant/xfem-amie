//
// C++ Interface: mechanical homogenization
//
// Description: 
//
//
// Author:  Alain Giorla, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef SCHEME_BASE_H
#define SCHEME_BASE_H


#include "properties_base.h"

namespace Mu
{

/* Default homogenization scheme. The actual schemes are inherited from this class*/
class HomogenizationScheme
{
protected:
	size_t nPhases ;
	std::vector<PropertiesType> input ;
	std::vector<PropertiesType> output ;

public:
	/* \brief constructor, using a single input, output is the same*/
	HomogenizationScheme(size_t n, PropertiesType in) ;
	/* \brief constructor, using a single input and a different output*/
	HomogenizationScheme(size_t n, PropertiesType in, PropertiesType out) ;
	/* \brief constructor, using a set of different input (output are the same)*/
	HomogenizationScheme(size_t n, std::vector<PropertiesType> & in) ;
	/* \brief constructor, using a set of different input, and another set of different output*/
	HomogenizationScheme(size_t n, std::vector<PropertiesType> & in, std::vector<PropertiesType> & out) ;

	/* \brief verifies that the scheme can be applied to the material*/
	virtual bool verify(const Material & mat) ;
	/* \brief returns the position of the input properties in the material*/
	virtual std::vector<std::vector<size_t> > getAllPositions(const std::vector<Material> & mat) ;
	/* \brief get the number of materials and the number of unknown*/
	virtual std::vector<size_t> getSchemeSize(const std::vector<Material> & mat) ;
	/* \brief get all data as a matrix */
	virtual Matrix getRawData(const std::vector<Material> & mat) ;
	/* \brief extracts the data, applies the scheme, and the returns the data as a new material*/
	virtual std::pair<bool, Material> apply(const std::vector<Material> & mat) ;
	/* \brief applies the model (here, returns the first material). This method is the only thing that should be overrided when creating a new model */
	virtual Vector processData(const Matrix & data) ;
} ;


/* \brief standard serial means of a set of data*/
class MeanSeries : public HomogenizationScheme
{
public:
	MeanSeries() ;
	virtual Vector processData(const Matrix & data) ;
} ;

/* \brief standard parallel means of a set of data*/
class MeanParallel : public HomogenizationScheme
{
public:
	MeanParallel() ;
	virtual Vector processData(const Matrix & data) ;
} ;




} ;

#endif // SCHEME

