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



HobbsScheme::HobbsScheme() : ExpansionHomogenizationScheme(2)
{
}

Vector HobbsScheme::process(const Matrix & data)
{
	double kmat = data[0][1] ;
	double amat = data[0][3] ;
	
	double finc = data[1][0] ;
	double kinc = data[1][1] ;
	double ainc = data[1][3] ;
	
	Vector processed(1) ;
	processed[0] = amat + 2*finc*kinc*(ainc-amat) / (kmat+kinc+finc*(kinc-kmat)) ;

	return processed ;
}


TurnerScheme::TurnerScheme() : ExpansionHomogenizationScheme(2)
{
}

Vector TurnerScheme::process(const Matrix & data)
{
	double fmat = data[0][0] ;
	double kmat = data[0][1] ;
	double amat = data[0][3] ;
	
	double finc = data[1][0] ;
	double kinc = data[1][1] ;
	double ainc = data[1][3] ;
	
	Vector processed(1) ;
	processed[0] = (amat*kmat*fmat + ainc*kinc*finc) / (kmat*fmat + kinc*finc) ;

	return processed ;
}


KernerScheme::KernerScheme() : ExpansionHomogenizationScheme(2)
{
}

Vector KernerScheme::process(const Matrix & data)
{
	double fmat = data[0][0] ;
	double kmat = data[0][1] ;
	double mumat= data[0][2] ;
	double amat = data[0][3] ;
	
	double finc = data[1][0] ;
	double kinc = data[1][1] ;
	double ainc = data[1][3] ;
	
	Vector processed(1) ;
	processed[0] = amat*fmat + ainc*finc ;
	processed[0] += (fmat*finc)*(ainc - amat)*(kinc-kmat)/(kmat*fmat+kinc*finc+3*kmat*kinc/(4*mumat)) ;

	return processed ;
}

HirschScheme::HirschScheme() : ExpansionHomogenizationScheme(2)
{

}

void HirschScheme::addScheme(ExpansionHomogenizationScheme exp, double p)
{
	expansion.push_back(exp) ;
	probability.push_back(p) ;
	double tot = 0 ;
	for(size_t i = 0 ; i < probability.size() ; i++)
		tot += probability[i] ;
	for(size_t i = 0 ; i < probability.size() ; i++)
		probability[i] /= tot ;
}

Vector HirschScheme::process(const Matrix & data) 
{
	Vector processed(1) ;
	for(size_t i = 0 ; i < expansion.size() ; i++)
	{
		Vector exp = expansion[i].process(data) ;
		processed[0] += exp[0]*probability[i] ;
	}
	return processed ;
}








SantScheme::SantScheme() : ExpansionHomogenizationScheme(2)
{
	input.clear() ;
	input.push_back(TAG_VOLUME_FRACTION) ;
	input.push_back(TAG_YOUNG_MODULUS) ;
	input.push_back(TAG_POISSON_RATIO) ;
	input.push_back(TAG_STRAIN) ;

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


