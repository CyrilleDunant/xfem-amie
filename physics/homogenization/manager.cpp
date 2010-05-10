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

#include "manager.h"
#include "converter.h"

using namespace Mu ;


GeometryManager::GeometryManager(GeometryType g) : Scheme(1) 
{
	gType = g ;

	switch(gType)
	{
	case CIRCLE:
		input.push_back(TAG_CIRCLE_RADIUS) ;
		output.push_back(TAG_PERIMETER) ;
		output.push_back(TAG_AREA) ;
		break;
	case ELLIPSE:	
		input.push_back(TAG_ELLIPSE_A) ;
		input.push_back(TAG_ELLIPSE_B) ;
		output.push_back(TAG_PERIMETER) ;
		output.push_back(TAG_AREA) ;
		output.push_back(TAG_ELLIPSE_FIRST_COMPLETE_INTEGRAL) ;
		output.push_back(TAG_ELLIPSE_SECOND_COMPLETE_INTEGRAL) ;
		break;
	}
}

Vector GeometryManager::process(const Matrix & data)
{
	switch(gType)
	{
		case CIRCLE:
			return processCircle(data) ;
		case ELLIPSE:
			return processEllipse(data) ;
	}
	Vector processed(output.size()) ;
	return processed ;
}

Vector GeometryManager::processCircle(const Matrix & data)
{
	Vector processed(2) ;
	double r = data[0][0] ;
	lessThanZero(r) ;
	processed[0] = 2. * M_PI * r ;
	processed[1] = M_PI * r * r ;
	return processed ;
}

Vector GeometryManager::processEllipse(const Matrix & data)
{
	GeneralConverter conv(TAG_NULL) ;

	Vector processed(4) ;
	double a = data[0][0] ;
	double b = data[0][1] ;

	lessThanZero(a) ;
	lessThanZero(b) ;
	lessThanZero(a-b) ;

	processed[2] = conv.getEllipseFirstCompleteIntegral(a,b) ;
	processed[3] = conv.getEllipseSecondCompleteIntegral(a,b) ;

	if(!conv.isOK())
		s = STATUS_BAD_HOMOGENIZATION ;

	processed[0] = 4. * a * processed[3] ;
	processed[1] = M_PI * a * b ;

	return processed ;
}










