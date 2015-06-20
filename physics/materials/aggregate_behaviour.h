// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2013
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef AGGREGATE_BEHAVIOUR_H
#define AGGREGATE_BEHAVIOUR_H

#include "../weibull_distributed_stiffness.h"
#include "../../geometry/geometry_base.h"

namespace Amie
{
struct AggregateBehaviour : public WeibullDistributedStiffness
{
    double up ;
    double yield ;
    double c ;
    //Yield 30 MPa
    AggregateBehaviour(double E=59e9, double nu=0.3, double up = 30e6, double yield = 0.00044, double c = 12000., SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) ;

    virtual Form * getCopy() const ;

} ;

struct ElasticOnlyAggregateBehaviour : public AggregateBehaviour
{
    ElasticOnlyAggregateBehaviour(double E=59e9, double nu=0.3, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) ;

    virtual Form * getCopy() const ;

} ;

struct ViscoElasticOnlyAggregateBehaviour : public AggregateBehaviour
{
    int freeblocks ;
    ViscoElasticOnlyAggregateBehaviour(double E=59e9, double nu=0.3, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) ;

    virtual Form * getCopy() const ;
} ;

struct ViscoDamageAggregateBehaviour : public AggregateBehaviour
{
    double rad ;
    int freeblocks ;
    ViscoDamageAggregateBehaviour(double E=59e9, double nu=0.3, double up = 0.00025, double r = 0.00025, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) ;

    virtual Form * getCopy() const ;
} ;

}

#endif // PASTE_BEHAVIOUR
