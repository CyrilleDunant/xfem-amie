//
// C++ Interface: stiffness_and_fracture
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "stiffness_with_large_strains_and_fracture.h"
#include "../mesher/delaunay.h"
#include "fracturecriteria/ruptureenergy.h"
#include "fracturecriteria/mcft.h"
#include "fracturecriteria/mohrcoulomb.h"
#include "damagemodels/nonlocalisotropiclineardamage.h"
#include "damagemodels/prandtlgrauertplasticstrain.h"
#include "damagemodels/fiberbasedisotropiclineardamage.h"

using namespace Amie ;


StiffnessWithLargeDeformationAndFracture::StiffnessWithLargeDeformationAndFracture(const Matrix & rig, FractureCriterion * crit, DamageModel * d) : LinearForm(rig, false, true, rig.numRows()/3+1)
{
    if(!d)
        dfunc = new FiberBasedIsotropicLinearDamage() ;
    else
        dfunc = d ;
    criterion = crit ;

    v.push_back(XI);
    v.push_back(ETA);
    if(param.size() == 36 )
    {
        v.push_back(ZETA);
    }

    this->time_d = false ;
    
    if(v.size() == 2)
    {
        imposed.resize(3, 0.) ;
    }
    else 
    {
        imposed.resize(6, 0.) ;
    }
    
// 	v.push_back(TIME_VARIABLE);
}

StiffnessWithLargeDeformationAndFracture::StiffnessWithLargeDeformationAndFracture(double E, double nu, FractureCriterion * crit, DamageModel * d, SpaceDimensionality dim, planeType pt, IsotropicMaterialParameters hooke) : LinearForm(Tensor::cauchyGreen(E, nu, dim, pt, hooke), false, true, dim)
{
    if(!d)
        dfunc = new FiberBasedIsotropicLinearDamage() ;
    else
        dfunc = d ;
    criterion = crit ;

    v.push_back(XI);
    v.push_back(ETA);
    if(param.size() == 36 )
    {
        v.push_back(ZETA);
    }
    if(v.size() == 2)
    {
        imposed.resize(3, 0.) ;
    }
    else 
    {
        imposed.resize(6, 0.) ;
    }}

StiffnessWithLargeDeformationAndFracture::~StiffnessWithLargeDeformationAndFracture()
{
    delete criterion ;
    delete dfunc ;
}

FractureCriterion * StiffnessWithLargeDeformationAndFracture::getFractureCriterion() const
{
    return criterion ;
}

DamageModel * StiffnessWithLargeDeformationAndFracture::getDamageModel() const
{
    return dfunc ;
}

void StiffnessWithLargeDeformationAndFracture::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix &ret, VirtualMachine * vm) const
{
    vm->ieval(Gradient(p_i) * dfunc->apply(param) * Gradient(p_j, true), gp, Jinv,v, ret) ;

}

void StiffnessWithLargeDeformationAndFracture::step(double timestep, ElementState & currentState, double maxscore)
{

    change = false ;
    if(!active)
    {
       
        for(size_t i = 0 ;  i < currentState.getParent()->getBoundingPoints().size() ; i++)
        {
            if(v.size() == 2)
            { 
                initialVolume = currentState.getParent()->area() ;
                initialPosition.push_back(currentState.getParent()->getBoundingPoint(i).getX()) ;
                initialPosition.push_back(currentState.getParent()->getBoundingPoint(i).getY()) ;
            }
            else
            {
                initialVolume = currentState.getParent()->volume() ;
                initialPosition.push_back(currentState.getParent()->getBoundingPoint(i).getX()) ;
                initialPosition.push_back(currentState.getParent()->getBoundingPoint(i).getY()) ;
                initialPosition.push_back(currentState.getParent()->getBoundingPoint(i).getZ()) ;
            }
        }
        change = true ;
        active = true ;
        
        dfunc->step(currentState, maxscore) ;
        currentState.getParent()->behaviourUpdated = dfunc->changed() ;
        currentState.getParent()->needAssembly = currentState.getParent()->behaviourUpdated ;
        return ;
    }
    
    VirtualMachine vm ;
    std::vector<double> relativeDisplacements ;
    Vector pimposed = imposed ;
    double currentVolume = 0 ;
    for(size_t i = 0 ;  i < currentState.getParent()->getBoundingPoints().size() ; i++)
    {
        if(v.size() == 2)
        {
            currentVolume  = currentState.getParent()->area() ;
            relativeDisplacements.push_back(initialPosition[i*2]-currentState.getParent()->getBoundingPoint(i).getX()) ;
            relativeDisplacements.push_back(initialPosition[i*2+1]-currentState.getParent()->getBoundingPoint(i).getY()) ;
        }
        else
        {
            currentVolume  = currentState.getParent()->volume() ;
            relativeDisplacements.push_back(initialPosition[i*3]-currentState.getParent()->getBoundingPoint(i).getX()) ;
            relativeDisplacements.push_back(initialPosition[i*3+1]-currentState.getParent()->getBoundingPoint(i).getY()) ;
            relativeDisplacements.push_back(initialPosition[i*3+2]-currentState.getParent()->getBoundingPoint(i).getZ()) ;
        }
    }

    
    if ( currentState.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
    { 
        Point p_(1./3., 1./3.) ;
        currentState.updateInverseJacobianCache(p_) ;   
        double x_xi = 0;
        double x_eta = 0;
        double y_xi = 0;
        double y_eta = 0;

        for ( size_t j = 0 ; j < currentState.getParent()->getShapeFunctions().size(); j++ )
        {
            if(j*2 >=  currentState.getDisplacements().size())
            {
                std::cerr << "displacement size mismatch" << std::endl ;
                break ;
            }
            double f_xi  = vm.deval ( currentState.getParent()->getShapeFunction ( j ), XI , p_ ) ;
            double f_eta = vm.deval ( currentState.getParent()->getShapeFunction ( j ), ETA, p_ ) ;
            
            x_xi  += f_xi * relativeDisplacements[j * 2] ;
            x_eta += f_eta * relativeDisplacements[j * 2] ;
            y_xi  += f_xi * relativeDisplacements[j * 2 + 1] ;
            y_eta += f_eta * relativeDisplacements[j * 2 + 1] ;
        }


        imposed[0] = ( x_xi ) * (*currentState.JinvCache)[0][0] + ( x_eta ) * (*currentState.JinvCache)[0][1] ;
        imposed[1] = ( y_xi ) * (*currentState.JinvCache)[1][0] + ( y_eta ) * (*currentState.JinvCache)[1][1] ;
        imposed[2] = ( ( x_xi ) * (*currentState.JinvCache)[1][0] + ( x_eta ) * (*currentState.JinvCache)[1][1]  + ( y_xi ) * (*currentState.JinvCache)[0][0] + ( y_eta ) * (*currentState.JinvCache)[0][1] );
    }
    else if ( currentState.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
    {
        Point p_(1./4., 1./4., 1./4.) ;
        currentState.updateInverseJacobianCache(p_) ;   
        double x_xi = 0;
        double x_eta = 0;
        double x_zeta = 0;
        double y_xi = 0;
        double y_eta = 0;
        double y_zeta = 0;
        double z_xi = 0;
        double z_eta = 0;
        double z_zeta = 0;

        for ( size_t j = 0 ; j < currentState.getParent()->getShapeFunctions().size() ; j++ )
        {
            double f_xi = vm.deval ( currentState.getParent()->getShapeFunction ( j ), XI, p_ ) ;
            double f_eta = vm.deval ( currentState.getParent()->getShapeFunction ( j ), ETA, p_ ) ;
            double f_zeta = vm .deval ( currentState.getParent()->getShapeFunction ( j ), ZETA, p_ ) ;
            double x = relativeDisplacements[j * 3] ;
            double y = relativeDisplacements[j * 3 + 1] ;
            double z = relativeDisplacements[j * 3 + 2] ;

            x_xi   += f_xi   * x ;
            x_eta  += f_eta  * x ;
            x_zeta += f_zeta * x ;
            y_xi   += f_xi   * y ;
            y_eta  += f_eta  * y ;
            y_zeta += f_zeta * y ;
            z_xi   += f_xi   * z ;
            z_eta  += f_eta  * z ;
            z_zeta += f_zeta * z ;
        }

        imposed[0] = ( x_xi ) * (*currentState.JinvCache)[0][0] + ( x_eta ) * (*currentState.JinvCache)[0][1]  + ( x_zeta ) * (*currentState.JinvCache)[0][2]+ (initialVolume-currentVolume)/initialVolume;
        imposed[1] = ( y_xi ) * (*currentState.JinvCache)[1][0] + ( y_eta ) * (*currentState.JinvCache)[1][1]  + ( y_zeta ) * (*currentState.JinvCache)[1][2]+ (initialVolume-currentVolume)/initialVolume;
        imposed[2] = ( z_xi ) * (*currentState.JinvCache)[2][0] + ( z_eta ) * (*currentState.JinvCache)[2][1]  + ( z_zeta ) * (*currentState.JinvCache)[2][2]+ (initialVolume-currentVolume)/initialVolume;

        imposed[3] = ( ( y_xi ) * (*currentState.JinvCache)[2][0] +
                    ( y_eta ) * (*currentState.JinvCache)[2][1] +
                    ( y_zeta ) * (*currentState.JinvCache)[2][2] +
                    ( z_xi ) * (*currentState.JinvCache)[1][0] +
                    ( z_eta ) * (*currentState.JinvCache)[1][1] +
                    ( z_zeta ) * (*currentState.JinvCache)[1][2] );

        imposed[4] = ( ( x_xi ) * (*currentState.JinvCache)[2][0] +
                    ( x_eta ) * (*currentState.JinvCache)[2][1] +
                    ( x_zeta ) * (*currentState.JinvCache)[2][2] +
                    ( z_xi ) * (*currentState.JinvCache)[0][0] +
                    ( z_eta ) * (*currentState.JinvCache)[0][1] +
                    ( z_zeta ) * (*currentState.JinvCache)[0][2] );

        imposed[5] = ( ( y_xi )   * (*currentState.JinvCache)[0][0] +
                    ( y_eta )  * (*currentState.JinvCache)[0][1] +
                    ( y_zeta ) * (*currentState.JinvCache)[0][2] +
                    ( x_xi )   * (*currentState.JinvCache)[1][0] +
                    ( x_eta )  * (*currentState.JinvCache)[1][1] +
                    ( x_zeta ) * (*currentState.JinvCache)[1][2] );
    }
    
    if(std::abs(imposed-pimposed).max() >  1e-6)
        change = true ;
    
    dfunc->step(currentState, maxscore) ;
    currentState.getParent()->behaviourUpdated = dfunc->changed() ;
    currentState.getParent()->needAssembly = currentState.getParent()->behaviourUpdated ;

}

Vector StiffnessWithLargeDeformationAndFracture::getImposedStress(const Point & p, IntegrableEntity * e, int g) const
{
    if(dfunc->hasInducedForces())
        return imposed*0. + dfunc->getImposedStress(p);
    return imposed*0. ;
}


Vector StiffnessWithLargeDeformationAndFracture::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const
{
    if(dfunc->hasInducedForces())
        return -imposed + dfunc->getImposedStrain(p);
    return -imposed ;
}

bool StiffnessWithLargeDeformationAndFracture::changed() const
{
    return dfunc->changed() || change;
}

bool StiffnessWithLargeDeformationAndFracture::fractured() const
{
    return dfunc->fractured() ;
}

Form * StiffnessWithLargeDeformationAndFracture::getCopy() const
{
    StiffnessWithLargeDeformationAndFracture * copy = new StiffnessWithLargeDeformationAndFracture(param, criterion->getCopy(), dfunc->getCopy()) ;
    copy->dfunc->getState(true).resize(dfunc->getState().size());
    copy->dfunc->getState(true) = dfunc->getState() ;
    copy->criterion->setMaterialCharacteristicRadius(criterion->getMaterialCharacteristicRadius()) ;
    copy->dfunc->setDamageDensityTolerance(dfunc->getDamageDensityTolerance());
    copy->dfunc->setThresholdDamageDensity(dfunc->getThresholdDamageDensity());

    return copy ;
}

Matrix StiffnessWithLargeDeformationAndFracture::getTensor(const Point & p, IntegrableEntity * e, int g) const
{
    return dfunc->apply(param, p, e, g) ;
}

void StiffnessWithLargeDeformationAndFracture::setFractureCriterion(FractureCriterion * frac)
{
    if(frac)
    {
        delete criterion ;
        criterion = frac ;
    }

}

std::vector<BoundaryCondition * > StiffnessWithLargeDeformationAndFracture::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
    std::vector<BoundaryCondition * > ret = dfunc->getBoundaryConditions(s,  id,  p_i, gp, Jinv);
    if(!active)
        return ret ;
    
    if(v.size() == 2)
    {
// 		Vector istress = VirtualMachine().ieval(Gradient(p_i)*(param * imposed),gp,Jinv, v)   ;
        Vector istress = getTensor(Point(1./3., 1./3.)) * imposed   ;
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, istress[0]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, istress[1]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_XI_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, istress[2]*.9));

    }
    if(v.size() == 3)
    {
// 		Vector istress = VirtualMachine().ieval(Gradient(p_i)*(param * imposed),gp,Jinv, v)   ;
        Vector istress = getTensor(Point(1./4., 1./4., 1./4.)) * imposed ;
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_XI, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, istress[0]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, istress[1]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, istress[2]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_XI_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, istress[3]*.9));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_XI_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, istress[4]*.9));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ETA_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, istress[5]*.9));
    }
    return ret ;
}


