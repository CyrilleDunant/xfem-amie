// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution

#ifndef VARIABLE_H
#define VARIABLE_H

#include "../geometry/geometry_base.h"

namespace Mu
{

const double default_derivation_delta= 1e-6 ;

/** \brief Possible variables. 
 * This is used as an index for arrays of powers, or as an argument to certain constructors. 
 * A special value of -1 is devoted to ONE as a special case of onstant variable.
 */
typedef enum : signed char
{
	ONE = -1,
	XI = 0,
	ETA = 1,
	ZETA = 2,
	TIME_VARIABLE = 3,
	U_VARIABLE,
	V_VARIABLE,
	W_VARIABLE
} Variable ;

} ;

#endif
