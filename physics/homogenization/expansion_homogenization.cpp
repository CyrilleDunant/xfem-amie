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




SelfConsistentSiderisScheme::SelfConsistentSiderisScheme() : ExpansionHomogenizationScheme(3)
{
	input.clear() ;
	input.push_back(TAG_VOLUME_FRACTION) ;
	input.push_back(TAG_YOUNG_MODULUS) ;
	input.push_back(TAG_POISSON_RATIO) ;
	input.push_back(TAG_EXPANSION_COEFFICIENT) ;
}

Vector SelfConsistentSiderisScheme::process(const Matrix & data)
{
	double fmat = data[0][0] ;
	double Emat = data[0][1] ;
	double numat= data[0][2] ;
	double amat = data[0][3] ;

	double finc = data[1][0] ;
	double Einc = data[1][1] ;
	double nuinc= data[1][2] ;
	double ainc = data[1][3] ;

	double fhom = data[2][0] ;
	double Ehom = data[2][1] ;
	double nuhom= data[2][2] ;

	double M = 3.*(1.-2*numat)*Ehom*finc*fhom ;
	double K = (1.-fhom)*(1.+numat)+2.*fmat*(1.-2.*numat) ;
	K *= Einc ;
	K += Emat*2.*fmat*(1.-2.*nuinc) ;

	Vector processed(1) ;
	processed[0] = amat + Einc*M/(Ehom*fhom*K)*(ainc-amat) ;

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



MultiLayerExpansion::MultiLayerExpansion(double eps) : ExpansionHomogenizationScheme(-1) 
{
	epsilon = eps ;
	input.clear() ;
	input.push_back(TAG_VOLUME_FRACTION) ;
	input.push_back(TAG_POISSON_RATIO) ;
	input.push_back(TAG_SHEAR_MODULUS) ;

	output.clear() ;
	output.push_back(TAG_STRAIN) ;
}

Vector MultiLayerExpansion::process(const Matrix & data)
{
	int n = data.numRows() ;

	Vector f(n) ;
	f[0] = data[0][0] ;
	Vector F(n) ;
	F[0] = data[0][0] ;
	for(int i = 1 ; i < n ; i++)
	{
		f[i] = data[i][0] ;
		F[i] = F[i-1] + f[i] ;
	}

	Vector N(n) ;
	for(int i = 0 ; i < n ; i++)
		N[i] = (1. - 2.*data[i][1])/(1. + data[i][1]) ;

	Vector M(n) ;
	for(int i = 0 ; i < n ; i++)
		M[i] = data[i][2] ;

	Vector X(n) ;
	Vector Y(n) ;
	Vector Z(n) ;

	Y[1] -= N[0]/(2.*M[0]) + (F[1] + 2*F[0]*N[1])/(4 * M[1] * f[1]) ;
	Z[1] += 1./M[1] + (F[1] + 2*F[0]*N[1])/(4 * M[1] * f[1]) ;

	for(int i = 2 ; i < n ; i++)
	{
		X[i] += F[i-2] * (1+2*N[i-1]) / (4 * M[i-1] * f[i-1]) ;
		Y[i] -= 1./M[i-1] + X[i] + (F[i] + 2*F[i-1]*N[i]) / (4 * M[i] * f[i]) ;
		Z[i] += 1./M[i] + (F[i] + 2*F[i-1]*N[i]) / (4 * M[i] * f[i]) ;
	}

	Vector T(n) ;
	T[n-2] = 1. ;
	for(int i = n-3 ; i > -1 ; i--)
		T[i] = - (Z[i+2]*T[i+2] + Y[i+2]*T[i+1]) / X[i+2] ;

	for(int i = 0 ; i < n ; i++)
		std::cout << F[i] << ";" ;
	std::cout << std::endl ;

	for(int i = 0 ; i < n ; i++)
		std::cout << X[i] << ";" ;
	std::cout << std::endl ;

	for(int i = 0 ; i < n ; i++)
		std::cout << Y[i] << ";" ;
	std::cout << std::endl ;

	for(int i = 0 ; i < n ; i++)
		std::cout << Z[i] << ";" ;
	std::cout << std::endl ;

	for(int i = 0 ; i < n ; i++)
		std::cout << T[i] << ";" ;
	std::cout << std::endl ;

	double Pn1 = - epsilon / (Y[1]*T[0] + Z[1]*T[1]) ;
	std::cout << Pn1 << std::endl ;

	Vector processed(1) ;
	processed[0] = F[n-2]*(1+N[n-1])*Pn1/(2*M[n-1]*f[n-1]) ;

	return processed ;

	
}


SelfConsistentMultiLayerExpansion::SelfConsistentMultiLayerExpansion(double eps) : MultiLayerExpansion(eps) 
{

}

Vector SelfConsistentMultiLayerExpansion::process(const Matrix & data)
{
	int n = data.numRows()-1 ;

	Vector f(n+1) ;
	f[0] = data[0][0] ;
	Vector F(n+1) ;
	F[0] = data[0][0] ;
	for(int i = 1 ; i < n+1 ; i++)
	{
		f[i] = data[i][0] ;
		F[i] = F[i-1] + f[i] ;
	}

	Vector N(n+1) ;
	for(int i = 0 ; i < n+1 ; i++)
		N[i] = (1. - 2.*data[i][1])/(1. + data[i][1]) ;

	Vector M(n+1) ;
	for(int i = 0 ; i < n+1 ; i++)
		M[i] = data[i][2] ;

	Vector X(n+1) ;
	Vector Y(n+1) ;
	Vector Z(n+1) ;

	Y[1] -= N[0]/(2.*M[0]) + (F[1] + 2*F[0]*N[1])/(4 * M[1] * f[1]) ;
	Z[1] += 1./M[1] + (F[1] + 2*F[0]*N[1])/(4 * M[1] * f[1]) ;

	for(int i = 2 ; i < n+1 ; i++)
	{
		X[i] += F[i-2] * (1+2*N[i-1]) / (4 * M[i-1] * f[i-1]) ;
		Y[i] -= 1./M[i-1] + X[i] + (F[i] + 2*F[i-1]*N[i]) / (4 * M[i] * f[i]) ;
		Z[i] += 1./M[i] + (F[i] + 2*F[i-1]*N[i]) / (4 * M[i] * f[i]) ;
	}

	Vector T(n+1) ;
	T[n-1] = 1. ;
	for(int i = n-2 ; i > -1 ; i--)
		T[i] = - (Z[i+2]*T[i+2] + Y[i+2]*T[i+1]) / X[i+2] ;

	Vector S(n+1) ;
	S[n-2] = 1./X[n] ;
	for(int i = n-3 ; i > -1 ; i--)
		S[i] = - (Z[i+2]*S[i+2] + Y[i+2]*S[i+1]) / X[i+2] ;

/*	for(int i = 0 ; i < n+1 ; i++)
		std::cout << F[i] << ";" ;
	std::cout << std::endl ;

	for(int i = 0 ; i < n+1 ; i++)
		std::cout << X[i] << ";" ;
	std::cout << std::endl ;

	for(int i = 0 ; i < n+1 ; i++)
		std::cout << Y[i] << ";" ;
	std::cout << std::endl ;

	for(int i = 0 ; i < n+1 ; i++)
		std::cout << Z[i] << ";" ;
	std::cout << std::endl ;

	for(int i = 0 ; i < n+1 ; i++)
		std::cout << T[i] << ";" ;
	std::cout << std::endl ;*/

	double R = F[n-2]*(1 - N[n-1]) / (2*f[n-1]) ;

	double eS = (Y[1]*S[0] + Z[1]*S[1]) ;
	double eP = (Y[1]*T[0] + Z[1]*T[1]) ;
	eP *= (M[n-1] - R * S[n-2]) ;
	eP /= (1. + R*(T[n-2] - 1)) ;

	std::cout << eS << ";" << eP << std::endl ;

	Vector processed(1) ;
	processed[0] =  - epsilon / (eS+eP) ;

	return processed ;

	
}



}


