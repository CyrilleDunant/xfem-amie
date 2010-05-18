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

#ifndef CONVERTER_H
#define CONVERTER_H


#include "scheme_base.h"

namespace Mu
{

/** 
* Converts any TAG_UNIVERSAL properties to a specific type
*/
class UniversalConverter : public Scheme
{
public:
	UniversalConverter(Tag t) ;

	virtual Vector process(const Matrix & data) ;
} ;

/**
* Converts a set of properties to a specific type. Each type needs a specific
* set of information. You may overwrite the following mehods to compute a specific
* property from a different formula.
*/
class GeneralConverter : public UniversalConverter
{
public:
	GeneralConverter(Tag t) ;

	virtual Vector process(const Matrix & data) ;

	/**
	* \brief get the volume from the density and the mass
	*/
	virtual double getVolume(double density, double mass) { return simpleDivision(mass,density) ; } ;
	/**
	* \brief get the volume fraction from the volume and the total volume
	*/
	virtual double getVolumeFraction(double volume, double volume_total) { return simpleDivision(volume,volume_total) ; } ;
	/**
	* \brief get the total volume from the volume and the volume fraction
	*/
	virtual double getVolumeTotal(double volume, double volume_fraction) { return simpleDivision(volume,volume_fraction) ; } ;

	/**
	* \brief get the mass from the density and the volume
	*/
	virtual double getMass(double density, double volume) { return volume*density ; } ;
	/**
	* \brief get the mass fraction from the mass and the total mass
	*/
	virtual double getMassFraction(double mass, double mass_total) { return simpleDivision(mass,mass_total) ; } ;
	/**
	* \brief get the total volume from the volume and the volume fraction
	*/
	virtual double getMassTotal(double mass, double mass_fraction) { return simpleDivision(mass,mass_fraction) ; } ;

	/**
	* \brief get the density from the volume and the mass
	*/
	virtual double getDensity(double volume, double mass) { return simpleDivision(mass,volume) ; } ;

	/**
	* \brief get the young modulus from the bulk and shear moduli
	*/
	virtual double getYoungModulus(double k, double mu) {return simpleDivision(9.*k*mu, 3.*k+mu) ; } ;
	/**
	* \brief get the poisson ratio from the bulk and shear moduli
	*/
	virtual double getPoissonRatio(double k, double mu) {return simpleDivision(3.*k-2.*mu, 6.*k+2.*mu) ; } ;
	/**
	* \brief get the bulk modulus from the young modulus and poisson ratio
	*/
	virtual double getBulkModulus(double E, double nu) {return simpleDivision(E, 3.*(1.-2.*nu)) ; } ;
	/**
	* \brief get the shear modulus from the young modulus and poisson ratio
	*/
	virtual double getShearModulus(double E, double nu) {return simpleDivision(E, 2.*(1.+nu)) ; } ;
	/**
	* \brief get the lame coefficient from the young modulus and poisson ratio
	*/
	virtual double getLameCoefficient(double E, double nu) {return simpleDivision(E*nu , (1.+nu)*(1.-2.*nu)) ; } ;

	/**
	* \brief computes the complete first elliptic integral for an ellipse defined by its larger and smaller radii
	*/
	virtual double getEllipseFirstCompleteIntegral(double a, double b) ;
	/**
	* \brief computes the second first elliptic integral for an ellipse defined by its larger and smaller radii
	*/
	virtual double getEllipseSecondCompleteIntegral(double a, double b) ;

} ;

/**
* Simple addition for a specific property type
*/
class AdditionConverter : public Scheme
{
public:
	AdditionConverter(Tag t) ;

	virtual Vector process(const Matrix & data) ;
} ;




} ;


#endif // CONVERTER_H

