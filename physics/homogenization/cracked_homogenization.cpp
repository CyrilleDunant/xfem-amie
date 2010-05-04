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

#include "cracked_homogenization.h"

namespace Mu
{






CrackedHomogenizationScheme::CrackedHomogenizationScheme() : HomogenizationScheme(1, BULK_SHEAR)
{
	input.push_back(CRACK_DENSITY) ;
	input.push_back(ELLIPSE_SHAPE) ;
	output.push_back(CRACK_DENSITY) ;
}



BudianskyScheme::BudianskyScheme() : CrackedHomogenizationScheme()
{
}

Vector BudianskyScheme::processData(const Matrix & data)
{
	double bulk = data[0][0] ;
	double shear = data[0][1] ;
	double nc = data[0][2] ;
	double a = data[0][3] ;
	double b = data[0][4] ;

	double k = sqrt(1. - (b*b)/(a*a)) ;
	double k1 = b/a ;


	// compute E(k) and K(k) 
	double elliptic_first = 1. ;
	double elliptic_second = 1. ;

	if(std::abs(k) < 1e-6)
	{
		elliptic_second = M_PI/2. ;
		elliptic_first = M_PI/2. ;
	}

	if(std::abs(k-1) > 1e-6)
	{
		double error = 1 ;
		int n = 1 ;
		while(error > 1e-6)
		{
			double last = elliptic_second ;
			int f2n = 2*n - 1 ;
			int q = f2n - 2 ;
			while(q > 1)
			{
				f2n *= q ;
				q -= 2 ;
			}
			int f2d = 2*n ;
			q = f2d - 2 ;
			while(q > 1)
			{
				f2d *= q ;
				q -= 2 ;
			}
			elliptic_first += ((double) f2n/(double) f2d)*((double) f2n/(double) f2d) * pow(k,(double)2*n) / (2*n-1) ;
			elliptic_second -= ((double) f2n/(double) f2d)*((double) f2n/(double) f2d) * pow(k,(double)2*n) / (2*n-1) ;
			error = std::abs(elliptic_second - last) ;
		}
		elliptic_first *= M_PI/2. ;
		elliptic_second *= M_PI/2. ;
	}

	double area = M_PI * a * b ;
	double perimeter = 4. * a * elliptic_second ;

	
	double epsilon = 2.*nc / M_PI * area*area / perimeter ;

	double nu = (3.*bulk-2.*shear) / (6.*bulk+2.*shear) ;
	double nux = nu + 1e-6 ;
	
	double A = (k1*k1*elliptic_first - elliptic_second) ;
	double B = k1*k1*(elliptic_second - elliptic_first) ;
	A /= (k*k*elliptic_second) ;
	B /= (k*k*elliptic_second) ;

	double epsilonx = 1 ;
	double Tx = 2. ;
	while(std::abs(epsilon - epsilonx)/epsilon > 1e-6 && nux > 1e-6)
	{
		nux -= 1e-6 ;
		Tx = 1. / (1. + A * nux) + 1. / (1. + B * nux) ;
		epsilonx = 45. * (nu - nux) ;
		epsilonx /= 8. * (1. - nux*nux) ;
		epsilonx /= (2.*(1.+3.*nu) - (1.-2.*nu)*Tx) ;
	}

	double bulkx = bulk ;
	bulkx -= bulk * 16./9. * (1. - nux*nux)/(1. - 2.*nux) * epsilon ;

	double shearx = shear ;
	shearx -= shear * 32./45.*(1. - nux)*(1. + 0.75 * Tx) * epsilon ;



/*	double kmat = data[0][1] ;
	double amat = data[0][3] ;
	
	double finc = data[0][0] ;
	double kinc = data[0][1] ;
	double ainc = data[0][3] ;*/
	
	Vector processed(3) ;
	processed[0] = bulkx ;
	processed[1] = shearx ;
	processed[2] = epsilon ;
	
	return processed ;
}

SimplifiedBenHahaScheme::SimplifiedBenHahaScheme() : CrackedHomogenizationScheme()
{
	input.pop_back() ;
	output.pop_back() ;
}

Vector SimplifiedBenHahaScheme::processData(const Matrix & data)
{
	double bulk = data[0][0] ;
	double shear = data[0][1] ;
	double nc = data[0][2] ;
	
	double nu = (3.*bulk-2.*shear) / (6.*bulk+2.*shear) ;

	double A = 16.*nc*(1.+3.*nu) ;
	double B = -45. ;
	double C = nu - A ;

	double nux = nu ;

	double delta = B*B - 4.*A*C ;
	
	if(delta > POINT_TOLERANCE)
	{
		double root1 = -B - sqrt(delta) ;
		root1 /= 2.*A ;
		double root2 = -B + sqrt(delta) ;
		root2 /= 2.*A ;

		if(root1 < 0.5)
			nux = std::max(0., root1) ;

	} else {
		if(std::abs(delta) < POINT_TOLERANCE)
		{
			nux = std::max(0.0, -B / (2.*A)) ;
			if(nux > nu)
				nux = nu ;

		} else {
			nux = nu ;

		}
	}

	double bulkx = bulk ;
	bulkx -= bulk * 16./9. * (1. - nux*nux)/(1. - 2.*nux) * nc ;

	double shearx = shear ;
	shearx -= shear * 32./45.*(1. - nux) * nc ;

	Vector processed(3) ;
	processed[0] = bulkx ;
	processed[1] = shearx ;
	processed[2] = nc ;
	
	return processed ;


}




}


