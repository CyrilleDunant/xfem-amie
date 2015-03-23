//
// C++ Interface: preconditionners
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//


#ifndef __PRECONDITIONNERS_H
#define __PRECONDITIONNERS_H
#include "../utilities/matrixops.h"

namespace Amie
{

/** \brief Abstract preconditionner interface*/
struct Preconditionner
{
    virtual ~Preconditionner() { } ;
    virtual void precondition(const Vector &, Vector &) = 0;
} ;

/** \brief Placeholder preconditionner, does nothing*/
struct NullPreconditionner : public Preconditionner
{
    virtual ~NullPreconditionner() { } ;
    virtual void precondition(const Vector &v, Vector &)  ;
} ;


}

#endif
