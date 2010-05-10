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
#include "manager.h"

namespace Mu
{

CrackedHomogenizationScheme::CrackedHomogenizationScheme() : Scheme(1)
{
	input.push_back(TAG_BULK_MODULUS) ;
	input.push_back(TAG_SHEAR_MODULUS) ;
	input.push_back(TAG_CRACK_DENSITY) ;

	output.push_back(TAG_BULK_MODULUS) ;
	output.push_back(TAG_SHEAR_MODULUS) ;
}



BudianskyScheme::BudianskyScheme() : CrackedHomogenizationScheme()
{
	input.push_back(TAG_ELLIPSE_A) ;
	input.push_back(TAG_ELLIPSE_B) ;

	output.push_back(TAG_CRACK_DENSITY) ;
}

Vector BudianskyScheme::process(const Matrix & data)
{
	double bulk = data[0][0] ;
	double shear = data[0][1] ;
	double nc = data[0][2] ;
	double a = data[0][3] ;
	double b = data[0][4] ;

	Matrix ellipse(1,2) ;
	ellipse[0][0] = a ;
	ellipse[0][1] = b ;

	GeometryManager man(ELLIPSE) ;

	Vector param_ellipse = man.processEllipse(ellipse) ;
	double perimeter = param_ellipse[0] ;
	double area = param_ellipse[1] ;
	double elliptic_first = param_ellipse[2] ;
	double elliptic_second = param_ellipse[3] ;

	double k1 = simpleDivision(b,a) ;
	double k = simpleSquareRoot(1. - k1*k1);

	double epsilon = 2.*nc / M_PI * area*area / perimeter ;

	double nu = (3.*bulk-2.*shear) / (6.*bulk+2.*shear) ;
	double nux = nu + 1e-6 ;
	
	double A = (k1*k1*elliptic_first - elliptic_second) ;
	double B = k1*k1*(elliptic_second - elliptic_first) ;

	double epsilonx = 1 ;
	double Tx = 2. ;
	double kEk = k*k*elliptic_second ;
	while(std::abs(epsilon - epsilonx)/epsilon > 1e-6 && nux > 1e-6)
	{
		nux -= 1e-6 ;
		Tx = kEk * ( 1. / (kEk + A * nux) + 1. / (kEk + B * nux)) ;
		epsilonx = 45. * (nu - nux) ;
		epsilonx /= 8. * (1. - nux*nux) ;
		epsilonx /= (2.*(1.+3.*nu) - (1.-2.*nu)*Tx) ;
	}

	double bulkx = bulk ;
	bulkx -= bulk * 16./9. * (1. - nux*nux)/(1. - 2.*nux) * epsilon ;

	double shearx = shear ;
	shearx -= shear * 32./45.*(1. - nux)*(1. + 0.75 * Tx) * epsilon ;



	Vector processed(3) ;
	processed[0] = bulkx ;
	processed[1] = shearx ;
	processed[2] = epsilon ;
	
	return processed ;
}

SimplifiedBenHahaScheme::SimplifiedBenHahaScheme() : CrackedHomogenizationScheme()
{

}

Vector SimplifiedBenHahaScheme::process(const Matrix & data)
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

	Vector processed(2) ;
	processed[0] = bulkx ;
	processed[1] = shearx ;
	
	return processed ;


}




}


