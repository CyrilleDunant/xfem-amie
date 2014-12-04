//
// C++ Interface: inversediagonal
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __SSOR_H
#define __SSOR_H

#include "preconditionners.h"
#include "../sparse/sparse_matrix.h"

namespace Amie
{

	/** \brief Preconditionner: inverse diagonal*/
	struct Ssor  : public Preconditionner
	{
		CoordinateIndexedSparseMatrix upper ;
                CoordinateIndexedSparseMatrix lower ;
                CoordinateIndexedSparseMatrix auxiliaryUpper ;
                CoordinateIndexedSparseMatrix auxiliaryLower ;
                Vector buffer ;
                double omega ;
                double omega_prev ;
                double previous_r ;
                bool converged = false ;
		virtual ~Ssor() { }
		Ssor(const CoordinateIndexedSparseMatrix &A, double omega) ;
		virtual void precondition(const Vector &v,Vector &)  ;
	} ;
	

} ;

#endif