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

#include "xml.h"

namespace Mu
{

typedef enum
{
	DILUTED,
	INCREMENTAL,
	MORI_TANAKA,
	SELF_CONSISTENT,
} HomogenizationScheme ;

std::pair<double,double> Enu2kmu(std::pair<double,double> E_nu) ;
std::pair<double,double> kmu2Enu(std::pair<double,double> k_mu) ;

class SimpleMaterial
{
protected:
	std::pair<double,double> Young_Poisson ;

public:
	/** \brief Simple consctructor */
	SimpleMaterial(double E, double nu) ;

	/** \brief Simple consctructor */
	SimpleMaterial(std::pair<double,double> E_nu) ;

	/** \brief Simple copy consctructor */
	SimpleMaterial(SimpleMaterial & mat) ;

	/** \brief Simple copy consctructor */
	SimpleMaterial(SimpleMaterial * mat) ;

	/** \brief Inclusion-Matrix isotropic elastic homogenization */
	SimpleMaterial(HomogenizationScheme scheme, std::pair<double,SimpleMaterial *> inclusions, SimpleMaterial * matrix) ;

	/** \brief Generalized Inclusion-Matrix isotropic elastic homogenization */
	SimpleMaterial(HomogenizationScheme scheme, std::vector<std::pair<double,SimpleMaterial *> > inclusions, SimpleMaterial * matrix) ;

	/** \brief Import from xml item*/
	SimpleMaterial(XMLTree * xml) ;

	/** \brief Get Young modulus and Poisson ratio*/
	std::pair<double,double> getEnu() {return Young_Poisson ; } ;

	/** \brief Get bulk and shear modulus*/
	std::pair<double, double> getkmu() {return Enu2kmu(Young_Poisson) ; } ;

	/** \brief Add the bulk and shear modulus of two materials*/
	void add(SimpleMaterial * mat) ;

	/** \brief Add the bulk and shear modulus of two materials*/
	void add(SimpleMaterial * mat, double f) ;

	/** \brief Multiply the bulk and shear modulus by f*/
	void multiply(double f) ;

	/** \brief Compute the relative difference between the Young modulus of the two materials*/
	double relativeDifference(SimpleMaterial * mat) ;

	/** \brief XML export*/
	XMLTree * toXML() ;
} ;



}

