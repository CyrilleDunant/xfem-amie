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

class UniversalConverter : public Scheme
{
public:
	UniversalConverter(Tag t) ;

	virtual Vector process(const Matrix & data) ;
} ;


class GeneralConverter : public UniversalConverter
{
public:
	GeneralConverter(Tag t) ;

	virtual Vector process(const Matrix & data) ;

	virtual double getVolume(double density, double mass) { return simpleDivision(mass,density) ; } ;
	virtual double getVolumeFraction(double volume, double volume_total) { return simpleDivision(volume,volume_total) ; } ;
	virtual double getVolumeTotal(double volume, double volume_fraction) { return simpleDivision(volume,volume_fraction) ; } ;

	virtual double getMass(double density, double volume) { return volume*density ; } ;
	virtual double getMassFraction(double mass, double mass_total) { return simpleDivision(mass,mass_total) ; } ;
	virtual double getMassTotal(double mass, double mass_fraction) { return simpleDivision(mass,mass_fraction) ; } ;

	virtual double getDensity(double volume, double mass) { return simpleDivision(mass,volume) ; } ;

	virtual double getYoungModulus(double k, double mu) {return simpleDivision(9.*k*mu, 3.*k+mu) ; } ;
	virtual double getPoissonRatio(double k, double mu) {return simpleDivision(3.*k-2.*mu, 6.*k+2.*mu) ; } ;
	virtual double getBulkModulus(double E, double nu) {return simpleDivision(E, 3.*(1.-2.*nu)) ; } ;
	virtual double getShearModulus(double E, double nu) {return simpleDivision(E, 2.*(1.+nu)) ; } ;

	virtual double getEllipseFirstCompleteIntegral(double a, double b) ;
	virtual double getEllipseSecondCompleteIntegral(double a, double b) ;

} ;

class AdditionConverter : public Scheme
{
public:
	AdditionConverter(Tag t) ;

	virtual Vector process(const Matrix & data) ;
} ;




} ;


#endif // CONVERTER_H

