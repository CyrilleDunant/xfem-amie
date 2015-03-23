//
// C++ Implementation: weibull_distributed_stiffness
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "spatially_distributed_stiffness.h"
#include "physics_base.h"
#include "fracturecriteria/mohrcoulomb.h"
#include "fracturecriteria/maxstrain.h"
#include "fracturecriteria/ruptureenergy.h"
#include "fracturecriteria/vonmises.h"
#include "stiffness.h"
#include "stiffness_and_fracture.h"
#include "../utilities/random.h"


using namespace Amie ;

SpatiallyDistributedStiffness::SpatiallyDistributedStiffness(const Matrix & rig, const Matrix & pore, double l,double ca, double cb) : LinearForm(rig, true, true, rig.numRows()/3+1), variability(.2), pore(pore)
{
    v.push_back(XI);
    v.push_back(ETA);
    if(param.size() == 36)
        v.push_back(ZETA);

    length = l ;
    distance = l ;
    criteriona = ca ;
    criterionb = cb ;
}

SpatiallyDistributedStiffness::~SpatiallyDistributedStiffness() { }

void SpatiallyDistributedStiffness::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix &ret, VirtualMachine * vm) const
{
    vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v,ret) ;
}

bool SpatiallyDistributedStiffness::fractured() const
{
    return false ;
}

Form * SpatiallyDistributedStiffness::getCopy() const
{
    double weib = RandomNumber().weibull(1,5) ;
    double factor = 1 - variability + variability*weib ;
    Matrix newTensor(param*factor) ;
    Matrix por(pore*factor) ;
    double crita = criteriona*factor ;
    double critb = criterionb*factor ;
    newTensor = por + (newTensor - por) * distance / length ;
    crita *= (0.5+0.5*distance/length) ;
    critb *= (0.5+0.5*distance/length) ;
//	if(randomVar > 0.5)
    if(criteriona > 0)
    {
        StiffnessAndFracture* copy = new StiffnessAndFracture( newTensor, new MohrCoulomb(crita,critb) ) ;

        return copy ;
    }
    Stiffness* copy = new Stiffness( newTensor) ;

    return copy ;
//	return new Stiffness(pore) ;
//	return new Stiffness(/*pore*(1.-variability)+*/pore/**randomVar*variability*/) ;
}

void SpatiallyDistributedStiffness::setDistance(double d)
{
//  std::cout << d << std::endl ;
    distance = d;
}
