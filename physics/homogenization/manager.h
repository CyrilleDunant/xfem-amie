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

