// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution

#ifndef VARIABLE_H
#define VARIABLE_H

#include "../geometry/geometry_base.h"

namespace Mu
{

const double default_derivation_delta= 1e-6 ;

/** Possible variables. This is used as an index for arrays of powers, or as an argument to certain constructors. A special value of -1 is devoted to ONE as a special case of onstant variable.
 */
enum Variable
{
	ONE = -1,
	XI,
	ETA,
	ZETA,
	TIME_VARIABLE,
	U_VARIABLE,
	V_VARIABLE,
	W_VARIABLE
} ;
} ;

#endif
