//
// C++ Interface: multi-grid step for preconditionning purposes
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MULTIGRIDSTEP_H
#define MULTIGRIDSTEP_H

#include "preconditionners.h"
#include "multigrid.h"


namespace Mu {

/** \brief Preconditionner, perform a GS step
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
template<class MESH_T, class ETYPE>
struct MultiGridStep : public Preconditionner
{
	MultiGrid<MESH_T, ETYPE> * mg ;
	virtual ~MultiGridStep()  { };
	MultiGridStep(MultiGrid<MESH_T, ETYPE> * mg): mg(mg) { };
	
	virtual void precondition(const Vector &v,Vector & t) const 
	{
		mg->b = v ;
		mg->solve(v, NULL, 1e-8, 10, false) ;
		t = mg->x ;
	}

};

}

#endif
