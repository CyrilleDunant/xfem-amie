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

#include "converter.h"
#include "../../geometry/geometry_base.h"

using namespace Mu ;


UniversalConverter::UniversalConverter(Tag t) : Scheme(1, TAG_UNIVERSAL, t)
{

}

Vector UniversalConverter::process(const Matrix & data)
{
	Vector processed(1) ;
	processed[0] = data[0][0] ;
	return processed ;
}









GeneralConverter::GeneralConverter(Tag t) : UniversalConverter(t)
{
	switch(t)
	{
	case TAG_VOLUME:
		input.clear() ;
		input.push_back(TAG_DENSITY) ;
		input.push_back(TAG_MASS) ;	
		break;
	case TAG_VOLUME_FRACTION:
		input.clear() ;
		input.push_back(TAG_VOLUME) ;
		input.push_back(TAG_VOLUME_TOTAL) ;
		break;
	case TAG_VOLUME_TOTAL:
		input.clear() ;
		input.push_back(TAG_VOLUME) ;
		input.push_back(TAG_VOLUME_FRACTION) ;
		break;
	case TAG_MASS:
		input.clear() ;
		input.push_back(TAG_DENSITY) ;
		input.push_back(TAG_VOLUME) ;
		break;
	case TAG_MASS_FRACTION:
		input.clear() ;
		input.push_back(TAG_MASS) ;
		input.push_back(TAG_MASS_TOTAL) ;
		break;
	case TAG_MASS_TOTAL:
		input.clear() ;
		input.push_back(TAG_MASS) ;
		input.push_back(TAG_MASS_FRACTION) ;
		break;
	case TAG_DENSITY:
		input.clear() ;
		input.push_back(TAG_VOLUME) ;
		input.push_back(TAG_MASS) ;
		break;
	case TAG_YOUNG_MODULUS:
		input.clear() ;
		input.push_back(TAG_BULK_MODULUS) ;
		input.push_back(TAG_SHEAR_MODULUS) ;
		break;
	case TAG_POISSON_RATIO:
		input.clear() ;
		input.push_back(TAG_BULK_MODULUS) ;
		input.push_back(TAG_SHEAR_MODULUS) ;
		break;
	case TAG_BULK_MODULUS:
		input.clear() ;
		input.push_back(TAG_YOUNG_MODULUS) ;
		input.push_back(TAG_POISSON_RATIO) ;
		break;
	case TAG_SHEAR_MODULUS:
		input.clear() ;
		input.push_back(TAG_YOUNG_MODULUS) ;
		input.push_back(TAG_POISSON_RATIO) ;
		break;
	case TAG_LAME_COEFFICIENT:
		input.clear() ;
		input.push_back(TAG_YOUNG_MODULUS) ;
		input.push_back(TAG_POISSON_RATIO) ;
		break;
	case TAG_ELLIPSE_FIRST_COMPLETE_INTEGRAL:
		input.clear() ;
		input.push_back(TAG_ELLIPSE_A) ;
		input.push_back(TAG_ELLIPSE_B) ;
		break;
	case TAG_ELLIPSE_SECOND_COMPLETE_INTEGRAL:
		input.clear() ;
		input.push_back(TAG_ELLIPSE_A) ;
		input.push_back(TAG_ELLIPSE_B) ;
		break;
	}
}

Vector GeneralConverter::process(const Matrix & data)
{
	Vector processed(1) ;
	processed[0] = data[0][0] ;
	switch(output[0])
	{
	case TAG_VOLUME:
		processed[0] = getVolume(data[0][0],data[0][1]) ;
		break;
	case TAG_VOLUME_FRACTION:
		processed[0] = getVolumeFraction(data[0][0],data[0][1]) ;
		break;
	case TAG_VOLUME_TOTAL:
		processed[0] = getVolumeTotal(data[0][0],data[0][1]) ;
		break;
	case TAG_MASS:
		processed[0] = getMass(data[0][0],data[0][1]) ;
		break;
	case TAG_MASS_FRACTION:
		processed[0] = getMassFraction(data[0][0],data[0][1]) ;
		break;
	case TAG_MASS_TOTAL:
		processed[0] = getMassTotal(data[0][0],data[0][1]) ;
		break;
	case TAG_DENSITY:
		processed[0] = getDensity(data[0][0],data[0][1]) ;
		break;
	case TAG_YOUNG_MODULUS:
		processed[0] = getYoungModulus(data[0][0],data[0][1]) ;
		break;
	case TAG_POISSON_RATIO:
		processed[0] = getPoissonRatio(data[0][0],data[0][1]) ;
		break;
	case TAG_BULK_MODULUS:
		processed[0] = getBulkModulus(data[0][0],data[0][1]) ;
		break;
	case TAG_SHEAR_MODULUS:
		processed[0] = getShearModulus(data[0][0],data[0][1]) ;
		break;
	case TAG_LAME_COEFFICIENT:
		processed[0] = getLameCoefficient(data[0][0],data[0][1]) ;
		break;
	case TAG_ELLIPSE_FIRST_COMPLETE_INTEGRAL:
		processed[0] = getEllipseFirstCompleteIntegral(data[0][0],data[0][1]) ;
		break;
	case TAG_ELLIPSE_SECOND_COMPLETE_INTEGRAL:
		processed[0] = getEllipseSecondCompleteIntegral(data[0][0],data[0][1]) ;
		break;
	}
	return processed ;
}


double GeneralConverter::getEllipseFirstCompleteIntegral(double a, double b)
{
	double k = simpleSquareRoot(1. - simpleDivision(b*b,a*a)) ;

	if(s == STATUS_BAD_HOMOGENIZATION)
		return 0.;

	if(equalsZero(k))
		return M_PI / 2. ;

	if(equalsZero(k-1.))
		return 1e9 ;

	double elliptic_first = 1. ;
	double error = 1. ;
	int n = 1 ;
	while(error > 1e-9)
	{
		double last = elliptic_first ;
		double fact = 0.5 ;
		unsigned long q = 2*n - 1 ;
		while(q > 1)
		{
			fact *= (double) q / (double) (q+1) ;
			q -= 2 ;
		}
		fact *= fact ;
		elliptic_first += fact * pow(k,(double)2*n) ;
		error = std::abs(elliptic_first - last) ;
		n++ ;
	}
	elliptic_first *= M_PI/2. ;

	return elliptic_first ;
}

double GeneralConverter::getEllipseSecondCompleteIntegral(double a, double b)
{
	double k = simpleSquareRoot(1. - simpleDivision(b*b,a*a)) ;

	if(s == STATUS_BAD_HOMOGENIZATION)
		return 0.;

	if(equalsZero(k))
		return M_PI / 2. ;

	if(equalsZero(k-1.))
		return 1. ;

	double elliptic_second = 1. ;
	double error = 1. ;
	int n = 1 ;
	while(error > 1e-6)
	{
		double last = elliptic_second ;
		double fact = 0.5 ;
		unsigned long q = 2*n - 1 ;
		while(q > 1)
		{
			fact *= (double) q / (double) (q+1) ;
			q -= 2 ;
		}
		fact *= fact ;
		elliptic_second -= fact * pow(k,(double)2*n) / (2*n-1) ;
		error = std::abs(elliptic_second - last) ;
		n++ ;
	}
	elliptic_second *= M_PI/2. ;

	return elliptic_second ;

}




AdditionConverter::AdditionConverter(Tag t) : Scheme(-1,t)
{
	if(t == TAG_MASS_TOTAL)
	{
		input.clear() ;
		input.push_back(TAG_MASS) ;
	}
	if(t == TAG_VOLUME_TOTAL)
	{
		input.clear() ;
		input.push_back(TAG_VOLUME) ;
	}
	
}

Vector AdditionConverter::process(const Matrix & data)
{
	Vector processed(1) ;
	for(size_t i = 0 ; i < data.size() ; i++)
		processed[0] += data[i][0] ;
	return processed ;
}




