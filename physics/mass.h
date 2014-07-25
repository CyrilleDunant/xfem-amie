//
// C++ Interface: mass
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __MASS_H_
#define __MASS_H_

#include "physics_base.h"
#include "fracturecriteria/fracturecriterion.h"
#include "damagemodels/damagemodel.h"
#include "homogenization/homogenization_base.h"
#include "fracturecriteria/nonlocalvonmises.h"

namespace Mu
{

	struct Mass : public LinearForm
	{
		std::vector<Variable> v ;
		double density ;
		
		Mass(double rho, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) ;
		
		virtual ~Mass() ;
		
		virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;
		
		virtual bool fractured() const { return false ; }
		
		virtual Form * getCopy() const ;
		
	} ;
	
	
} ;

#endif
