// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef PASTE_BEHAVIOUR_H
#define PASTE_BEHAVIOUR_H

#include "../weibull_distributed_stiffness.h"
#include "../../geometry/geometry_base.h"

namespace Mu
{
	struct PasteBehaviour : public WeibullDistributedStiffness
	{
	    // ultimate tensile strength: 5.68 MPa
		double up ;
		double yield ;
		double c ;
		PasteBehaviour(double E=12e9, double nu=0.3,  double up = 0.0005, double yield = 0.00035, double c = 12000., SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) ;
		
		virtual Form * getCopy() const ;
		
	} ;
	
	struct ElasticOnlyPasteBehaviour : public PasteBehaviour
	{
		ElasticOnlyPasteBehaviour(double E=12e9, double nu=0.3, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) ;
		
		virtual Form * getCopy() const ;	
	} ;
	
	// parameters got for ./alain 1 500 500 0.29 10 200 at revision 2512
	struct ViscoElasticOnlyPasteBehaviour : public PasteBehaviour
	{
		double tau_kv ;
		double e_kv ;
		double tau_mx ;
		
		ViscoElasticOnlyPasteBehaviour(double E=12e9, double nu = 0.3, double tmx = 200, double ekv = 0.29, double tkv = 10, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) ;
		
		virtual Form * getCopy() const ;
	} ;
	
} ;

#endif // PASTE_BEHAVIOUR