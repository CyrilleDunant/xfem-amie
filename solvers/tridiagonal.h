//
// C++ Interface: inversediagonal
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __TRI_DIAG_H
#define __TRI_DIAG_H

#include "preconditionners.h"
#include "../sparse/sparse_matrix.h"

namespace Amie
{

/** \brief Tridiagonal preconditionner*/
struct TriDiagonal  : public Preconditionner
{
    Vector diagonal ;
    Vector upper ;
    Vector c ;
    virtual ~TriDiagonal() {
        ;
    }
    TriDiagonal(const CoordinateIndexedSparseMatrix &A) ;
    virtual void precondition(const Vector &v,Vector &)  ;
} ;

}

#endif
