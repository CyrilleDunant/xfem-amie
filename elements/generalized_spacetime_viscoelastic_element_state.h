// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution
//


#ifndef GENERALIZED_SPACETIME_VISCOELASTIC_ELEMENT_STATE
#define GENERALIZED_SPACETIME_VISCOELASTIC_ELEMENT_STATE

#include "integrable_entity.h"

namespace Mu
{

class GeneralizedSpaceTimeViscoElasticElementState : public ElementState
{
public:
	GeneralizedSpaceTimeViscoElasticElementState(IntegrableEntity * e) ;
	GeneralizedSpaceTimeViscoElasticElementState(const GeneralizedSpaceTimeViscoElasticElementState &s) ;
	GeneralizedSpaceTimeViscoElasticElementState & operator =(const GeneralizedSpaceTimeViscoElasticElementState & s) ;
	
	virtual void getField( FieldType f, const Point & p, Vector & ret, bool local, VirtualMachine * vm = nullptr, int i = 0) const ;

	virtual void getFieldAtNodes( FieldType f, Vector & ret, VirtualMachine * vm = nullptr, int i = 0) ;
	
	virtual void getField( FieldType f1, FieldType f2, const Point & p, Vector & ret1, Vector & ret2, bool local, VirtualMachine * vm = nullptr, int i = 0, int j = 0) const  ;
			
	virtual void getFieldAtNodes( FieldType f1, FieldType f2, Vector & ret1, Vector & ret2, VirtualMachine * vm = nullptr, int i = 0, int j = 0) ;
	
	virtual void getAverageField( FieldType f, Vector & ret, VirtualMachine * vm = nullptr, int i = 0, double t = 0) ;

	virtual void getAverageField( FieldType f, FieldType f_, Vector & ret, Vector & ret_, VirtualMachine * vm = nullptr, int dummy= 0, double t = 0)  ;
	
	
} ;

} ;


#endif // GENERALIZED_SPACETIME_VISCOELASTIC_ELEMENT_STATE
