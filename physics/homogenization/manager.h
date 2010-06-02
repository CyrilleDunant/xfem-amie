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

#ifndef MANAGER_H
#define MANAGER_H


#include "properties_base.h"
#include "scheme_base.h"
#include "../../geometry/geometry_base.h"

namespace Mu
{

/**
* The GeometryManager will give some information about a geometry (like the perimeter or the area)
* from knowing basic information (like the radius for a circle or the shape factor for an ellipse)
* Currently supported: circle, ellipse
*/
class GeometryManager : public Scheme
{
protected:
	GeometryType gType ;

public:
	GeometryManager(GeometryType g) ;

	virtual Vector process(const Matrix & data) ;
	virtual Vector processCircle(const Matrix & data) ;
	virtual Vector processEllipse(const Matrix & data) ;

} ;


} ;


#endif // MANAGER

