//
// C++ Interface: preconditionners
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//


#ifndef __PRECONDITIONNERS_H
#define __PRECONDITIONNERS_H
#include "../matrixops.h"

namespace Mu
{
	
	struct Preconditionner
	{
		virtual ~Preconditionner() { } ;
		virtual void precondition(const Vector &, Vector &) const = 0;
	} ;
	
	struct NullPreconditionner : public Preconditionner
	{
		virtual ~NullPreconditionner() { } ;
		virtual void precondition(const Vector &v, Vector &) const ;
	} ;


} ;

#endif
