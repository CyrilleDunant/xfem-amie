//
// C++ Implementation: mechanical analytic homogenization
//
// Description:
//
//
// Author:  Alain Giorla, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "expansion_homogenization.h"

namespace Mu
{






ExpansionHomogenizationScheme::ExpansionHomogenizationScheme(int i) : Scheme(i)
{
	input.push_back(TAG_VOLUME_FRACTION) ;
	input.push_back(TAG_BULK_MODULUS) ;
	input.push_back(TAG_SHEAR_MODULUS) ;
	input.push_back(TAG_EXPANSION_COEFFICIENT) ;

	output.push_back(TAG_EXPANSION_COEFFICIENT) ;
}



HashinScheme::HashinScheme() : ExpansionHomogenizationScheme(2)
{
}

Vector HashinScheme::process(const Matrix & data)
{
	double kmat = data[0][1] ;
	double amat = data[0][3] ;
	
	double finc = data[1][0] ;
	double kinc = data[1][1] ;
	double ainc = data[1][3] ;
	
	Vector processed(1) ;
	processed[0] = amat + 2*finc*kinc*(ainc-amat) / (kmat+kinc+finc*(kinc-kmat)) ;

//	std::cout << amat << ";" << ainc << ";" << 2*finc*kinc*(ainc-amat) << ";" << (kmat+kinc+finc*(kinc-kmat)) << std::endl ;
	
	return processed ;
}


SantScheme::SantScheme() : ExpansionHomogenizationScheme(2)
{
	input.clear() ;
	input.push_back(TAG_VOLUME_FRACTION) ;
	input.push_back(TAG_YOUNG_MODULUS) ;
	input.push_back(TAG_POISSON_RATIO) ;
	input.push_back(TAG_EXPANSION_COEFFICIENT) ;

	output.clear() ;
	output.push_back(TAG_STRAIN) ;
	output.push_back(TAG_STRAIN) ;
}

Vector SantScheme::process(const Matrix & data)
{
	double fmat = data[0][0] ;
	double Emat = data[0][1] ;
	double nmat = data[0][2] ;
	double amat = data[0][3] ;

	double finc = data[1][0] ;
	double Einc = data[1][1] ;
	double ninc = data[1][2] ;
	double ainc = data[1][3] ;

	double c = 1./finc ;

	double misfit = ainc - amat ;

	double ph = 2.*(1.-2.*ninc)*(Emat/Einc)*(c-1.) ;
	ph += c+2. + nmat*(c-4.) ;
	ph = simpleDivision(1.,ph) ;
	ph *= 2.*(c-1.)*Emat ;	
	ph *= misfit ;

	double sigmamatr = ph ;
	double sigmamatt = -ph * (2.*finc+1.)/2./(1.-finc) ;

	double epsilonmat = sigmamatt - nmat*(sigmamatr+sigmamatt) ;
	epsilonmat /= Emat ;
	epsilonmat -= amat ;

	double epsiloninc = (1.-2.*ninc)*ph ;
	epsiloninc /= Einc ;
	epsiloninc -= ainc ;

	Vector processed(2) ;
	processed[0] = epsilonmat ;
	processed[1] = epsiloninc ;

	return processed ;
}






}


