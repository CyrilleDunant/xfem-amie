
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution
//

#include "integrable_entity.h"
#include "../physics/fracturecriteria/fracturecriterion.h"
#include "../solvers/assembly.h"
#include "../features/boundarycondition.h"
#include "../physics/damagemodels/damagemodel.h"
#include <unistd.h>

namespace Amie {


size_t fieldTypeElementarySize ( FieldType f, SpaceDimensionality dim, size_t blocks )
{
    switch ( f )
    {
    case VON_MISES_STRAIN_FIELD :
    case VON_MISES_REAL_STRESS_FIELD :
    case VON_MISES_EFFECTIVE_STRESS_FIELD :
    case PRINCIPAL_STRESS_ANGLE_FIELD :
    case PRINCIPAL_STRAIN_ANGLE_FIELD :
    case SCALAR_DAMAGE_FIELD :
    case INTERNAL_VARIABLE_FIELD:
        return 1 ;

    case DISPLACEMENT_FIELD :
    case ENRICHED_DISPLACEMENT_FIELD :
    case SPEED_FIELD :
    case FLUX_FIELD:
    case GRADIENT_FIELD:
    case PRINCIPAL_STRAIN_FIELD:
    case PRINCIPAL_MECHANICAL_STRAIN_FIELD:
    case PRINCIPAL_EFFECTIVE_STRESS_FIELD:
    case PRINCIPAL_REAL_STRESS_FIELD:
    case PRINCIPAL_IMPOSED_STRAIN_FIELD:
    case PRINCIPAL_IMPOSED_STRESS_FIELD:
    case TENSOR_DAMAGE_FIELD:
        return (dim == SPACE_THREE_DIMENSIONAL) ? 3 : 2 ;

    case TOTAL_STRAIN_FIELD:
    case STRAIN_RATE_FIELD:
    case MECHANICAL_STRAIN_FIELD:
    case EFFECTIVE_STRESS_FIELD:
    case IMPOSED_STRESS_FIELD:    
    case IMPOSED_STRAIN_FIELD:    
    case REAL_STRESS_FIELD:
    case NON_ENRICHED_STRAIN_FIELD:
    case NON_ENRICHED_STRAIN_RATE_FIELD:
    case NON_ENRICHED_EFFECTIVE_STRESS_FIELD:
    case NON_ENRICHED_REAL_STRESS_FIELD:
        return (dim == SPACE_THREE_DIMENSIONAL) ? 6 : 3 ;


    case GENERALIZED_VISCOELASTIC_DISPLACEMENT_FIELD :
    case GENERALIZED_VISCOELASTIC_ENRICHED_DISPLACEMENT_FIELD :
    case GENERALIZED_VISCOELASTIC_SPEED_FIELD :
        return (dim == SPACE_THREE_DIMENSIONAL) ? 3*blocks : 2*blocks ;

    case GENERALIZED_VISCOELASTIC_STRAIN_FIELD:
    case GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD:
    case GENERALIZED_VISCOELASTIC_EFFECTIVE_STRESS_FIELD:
    case GENERALIZED_VISCOELASTIC_REAL_STRESS_FIELD:
    case GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_FIELD:
    case GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_RATE_FIELD:
    case GENERALIZED_VISCOELASTIC_NON_ENRICHED_EFFECTIVE_STRESS_FIELD:
    case GENERALIZED_VISCOELASTIC_NON_ENRICHED_REAL_STRESS_FIELD:
        return (dim == SPACE_THREE_DIMENSIONAL) ? 6*blocks : 3*blocks ;
    default:
    {
        std::cout << f << " non-handled field !" << std::endl ;
        exit ( 0 ) ;
        return 0 ;
    }
    }
    std::cout << f << " non-existing field !" << std::endl ;
    exit ( 0 ) ;
    return 0 ;
}


ElementState * Form::createElementState ( IntegrableEntity * e )
{
    return new ElementState ( e ) ;
}

IntegrableEntity::IntegrableEntity() : boundaryConditionCache ( nullptr ), cachedGps ( nullptr )
{
    state = new ElementState ( this ) ;
    enrichmentUpdated = false ;
    behaviourUpdated = false ;
    behaviourForcesUpdated = false ;
    behaviourViscoUpdated = false ;
    needAssembly = true ;
}


Function IntegrableEntity::getZTransform() const
{
    return Function ( "1" ) ;
}

Function IntegrableEntity::getTTransform() const
{
    return Function ( "1" ) ;
}

Function IntegrableEntity::getXTransformAtCentralNodalTime() const
{
    return getXTransform() ;
}

Function IntegrableEntity::getYTransformAtCentralNodalTime() const
{
    return getYTransform() ;
}

Function IntegrableEntity::getZTransformAtCentralNodalTime() const
{
    return getZTransform() ;
}

Function IntegrableEntity::getTTransformAtCentralNodalTime() const
{
    return getTTransform() ;
}

Vector Form::getImposedStress ( const Point & p, IntegrableEntity * e, int g ) const
{
    if ( getDamageModel() && getDamageModel()->hasInducedForces() )
    {
        return getDamageModel()->getImposedStress ( p ) ;
    }

    return Vector ( double ( 0 ), getTensor ( p, e ).numCols() ) ;
}

Vector Form::getImposedStrain ( const Point & p, IntegrableEntity * e, int g ) const
{
    if ( getDamageModel() && getDamageModel()->hasInducedForces() )
    {
        return getDamageModel()->getImposedStrain ( p ) ;
    }

    return Vector ( double ( 0 ), getTensor ( p ).numCols() ) ;
}

bool Form::hasInducedForces() const
{
    return false || ( getDamageModel() && getDamageModel()->hasInducedForces() );
}

std::vector<BoundaryCondition * > Form::getBoundaryConditions ( const ElementState & s, size_t id,  const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv ) const
{
    std::vector<BoundaryCondition * > ret ;
    if ( getDamageModel() && getDamageModel()->hasInducedBoundaryConditions() )
    {
        return getDamageModel()->getBoundaryConditions ( s, id,  p_i, gp, Jinv ) ;
    }

    return  ret ;
}

// void Form::setFractureCriterion(FractureCriterion * frac)
// {
// 	if(frac)
// 		this->getFractureCriterion() = frac ;
// }

Matrix Form::getTensorDot ( const Point & p, double dt, bool calc ) const 
{
    if(!calc || dt < POINT_TOLERANCE)
        return param * 0 ;

    Point prev(p.getX(), p.getY(), p.getZ(), p.getT()-default_derivation_delta) ;
    Point next(p.getX(), p.getY(), p.getZ(), p.getT()+default_derivation_delta) ;

    return (getTensor(next)-getTensor(prev))/(dt*2.*default_derivation_delta) ;
}

void Form::getTensorDotAtGaussPoints( const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<std::pair<Matrix, Matrix> > & ret, bool calc ) const
{
    if(ret.size() != gp.gaussPoints.size())
       ret.resize( gp.gaussPoints.size(), std::make_pair(param, param*0.) ) ;

    for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
    {
        ret[g].first = getTensor( gp.gaussPoints[g].first ) ;
        ret[g].second = getTensorDot( gp.gaussPoints[g].second, 1./(Jinv[g][ Jinv[g].numRows()-1 ] [ Jinv[g].numCols()-1 ]), calc ) ;
    }
}

void Form::setSource ( const Geometry *  src ) 
{
    source = src ;
    if(getFractureCriterion())
    {
        getFractureCriterion()->setMaterialCharacteristicRadius(std::min(src->getRadius()*.5, getFractureCriterion()->getMaterialCharacteristicRadius())) ;
    }
}

Matrix Form::getViscousTensorDot ( const Point & p, double dt, bool calc ) const
{
    if(!calc || dt < POINT_TOLERANCE)
        return param * 0 ;

    Point prev(p.getX(), p.getY(), p.getZ(), p.getT()-default_derivation_delta) ;
    Point next(p.getX(), p.getY(), p.getZ(), p.getT()+default_derivation_delta) ;

    return (getViscousTensor(next)-getViscousTensor(prev))/(dt*2.*default_derivation_delta) ;
}

void Form::getViscousTensorDotAtGaussPoints( const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<std::pair<Matrix, Matrix> > & ret, bool calc ) const
{
    if(ret.size() != gp.gaussPoints.size())
       ret.resize( gp.gaussPoints.size(), std::make_pair(param*0, param*0.) ) ;

    for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
    {
        ret[g].first = getViscousTensor( gp.gaussPoints[g].first ) ;
        ret[g].second = getViscousTensorDot( gp.gaussPoints[g].second, 1./(Jinv[g][ Jinv[g].numRows()-1 ] [ Jinv[g].numCols()-1 ]), calc ) ;
    }
}

void IntegrableEntity::applyBoundaryCondition ( Assembly *a )
{
    if ( !getBehaviour() || !(getBehaviour()->getDamageModel() && (getBehaviour()->getDamageModel()->hasInducedBoundaryConditions() || !getBehaviour()->hasInducedForces())))
    {
        return ;
    }
    if ( getBehaviour()->type != VOID_BEHAVIOUR  )
    {
        if ( boundaryConditionCache && !behaviourForcesUpdated )
        {
            for ( size_t i = 0 ; i < boundaryConditionCache->size() ; i++ )
            {
                if ( ( *boundaryConditionCache ) [i] )
                {
                    if ( get2DMesh() )
                    {
                        ( *boundaryConditionCache ) [i]->apply ( a, get2DMesh() ) ;
                    }
                    else
                    {
                        ( *boundaryConditionCache ) [i]->apply ( a, get3DMesh() ) ;
                    }
                }
            }
            
            return ;
        }
        if( boundaryConditionCache && behaviourForcesUpdated ) 
            boundaryConditionCache->clear() ;
        else
            boundaryConditionCache = new std::vector<BoundaryCondition *>() ;

        int JinvSize = 3 ;
        if ( spaceDimensions() == SPACE_THREE_DIMENSIONAL && timePlanes() > 1 )
            JinvSize = 4 ;
        if ( spaceDimensions() == SPACE_TWO_DIMENSIONAL && timePlanes() == 1 )
            JinvSize = 2 ;
        std::valarray<Matrix> Jinv ( (bool) getState().JinvCache ? (*getState().JinvCache) : Matrix( JinvSize, JinvSize ),  getGaussPoints().gaussPoints.size()) ;

        if( ! getState().JinvCache )
        {
            for ( size_t i = 0 ; i < getGaussPoints().gaussPoints.size() ;  i++ )
            {
                getInverseJacobianMatrix ( getGaussPoints().gaussPoints[i].first, Jinv[i] ) ;
            }
        }


        size_t start = 0 ;
        if ( timePlanes() > 1 )
        {
            start = getBoundingPoints().size() - getBoundingPoints().size() /timePlanes() ;
        }
        for ( size_t i = start ; i < getShapeFunctions().size() ; i++ )
        {
            std::vector<BoundaryCondition *> boundaryConditionCachetmp = getBehaviour()->getBoundaryConditions ( getState(), getBoundingPoint ( i ).getId(),  getShapeFunction ( i ), getGaussPoints(), Jinv );
            boundaryConditionCache->insert(boundaryConditionCache->end(), boundaryConditionCachetmp.begin(), boundaryConditionCachetmp.end()) ;
            for ( size_t j = 0 ; j < boundaryConditionCachetmp.size() ; j++ )
            {
                if ( get2DMesh() )
                {
                    boundaryConditionCachetmp[j]->apply ( a, get2DMesh() ) ;
                }
                else
                {
                    boundaryConditionCachetmp[j]->apply ( a, get3DMesh() ) ;
                } 
            }
        }
        for ( size_t i = start ; i < getEnrichmentFunctions().size() ; i++ )
        {
            std::vector<BoundaryCondition *> boundaryConditionCachetmp = getBehaviour()->getBoundaryConditions ( getState(), getEnrichmentFunction ( i ).getDofID(),  getEnrichmentFunction ( i ), getGaussPoints(), Jinv ) ;
            boundaryConditionCache->insert(boundaryConditionCache->end(), boundaryConditionCachetmp.begin(), boundaryConditionCachetmp.end()) ;
            for ( size_t j = 0 ; j < boundaryConditionCachetmp.size() ; j++ )
            {
                if ( get2DMesh() )
                {
                    boundaryConditionCachetmp[j]->apply ( a, get2DMesh() ) ;
                }
                else
                {
                    boundaryConditionCachetmp[j]->apply ( a, get3DMesh() ) ;
                }
            }
        }

    }
}

IntegrableEntity::~IntegrableEntity()
{
    if ( boundaryConditionCache )
    {
        for ( size_t i = 0 ; i < boundaryConditionCache->size() ; i++ )
        {
            delete ( *boundaryConditionCache ) [i] ;
        }
    }
    delete state ;
    delete boundaryConditionCache ;
    delete cachedGps ;
}


void IntegrableEntity::addBoundaryCondition ( BoundaryCondition * bc )
{
    if ( !boundaryConditionCache )
    {
        boundaryConditionCache = new std::vector<BoundaryCondition *> ( 0, nullptr ) ;
    }
    boundaryConditionCache->push_back ( bc );
}

void IntegrableEntity::clearBoundaryConditions()
{
    if ( boundaryConditionCache )
    {
        for ( size_t i = 0 ; i < boundaryConditionCache->size() ; i++ )
        {
            delete ( *boundaryConditionCache ) [i] ;
        }
        delete boundaryConditionCache ;
    }
    boundaryConditionCache = nullptr ;
}


void IntegrableEntity::setState ( ElementState * s )
{
    delete state ;
    state = s ;
}

const ElementState &IntegrableEntity::getState() const
{
    return *state ;
}

ElementState &IntegrableEntity::getState()
{
    return *state ;
}


const Vector &ElementState::getDisplacements() const
{
    return this->displacements ;
}

const Vector &ElementState::getEnrichedDisplacements() const
{
    return this->enrichedDisplacements ;
}


ElementState::ElementState ( IntegrableEntity *s )
{
    JinvCache = nullptr ;
    mesh2d = nullptr ; //s->get2DMesh() ;
    mesh3d = nullptr ; //s->get3DMesh() ;
    parent = s ;
    lock = false ;
    timePos = 0 ;
    previousTimePos = 0 ;
// 	size_t ndof = 2 ;
// // 	if(s->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
// // 		ndof = 3 ;
// 	this->displacements.resize(s->getBoundingPoints()->size()*ndof) ;
// 	this->previousDisplacements.resize(s->getBoundingPoints()->size()*ndof) ;
// 	this->previousPreviousDisplacements.resize(s->getBoundingPoints()->size()*ndof) ;
// 	this->enrichedDisplacements.resize(s->getEnrichmentFunctions().size()*ndof) ;
// 	this->previousEnrichedDisplacements.resize(s->getEnrichmentFunctions().size()*ndof) ;
// 	this->previousPreviousEnrichedDisplacements.resize(s->getEnrichmentFunctions().size()*ndof) ;
}

Matrix makeStressOrStrainMatrix ( const Vector & stressOrStrain )
{
    if ( stressOrStrain.size() == 3 )
    {
        Matrix ret2 ( 2,2 ) ;
        ret2[0][0] = stressOrStrain[0] ;
        ret2[1][1] = stressOrStrain[1] ;
        ret2[1][0] = ret2[0][1] = stressOrStrain[2]  ;
        return ret2 ;
    }
    else if ( stressOrStrain.size() == 6 )
    {
        Matrix ret3 ( 3,3 ) ;
        ret3[0][0] = stressOrStrain[0] ;
        ret3[1][1] = stressOrStrain[1] ;
        ret3[2][2] = stressOrStrain[2] ;
        ret3[2][0] = ret3[0][2] = stressOrStrain[3]  ;
        ret3[2][1] = ret3[1][2] = stressOrStrain[4]  ;
        ret3[1][0] = ret3[0][1] = stressOrStrain[5]  ;
        return ret3 ;
    }
    return Matrix ( 2 + ( stressOrStrain.size() ==6 ), 2 + ( stressOrStrain.size() ==6 ) ) ;
}

bool isStrainField ( FieldType f )
{
    return f == TOTAL_STRAIN_FIELD || f == NON_ENRICHED_STRAIN_FIELD || f == MECHANICAL_STRAIN_FIELD /* || f == PRINCIPAL_STRAIN_FIELD*/ ;
}

bool isStressField ( FieldType f )
{
    return f == REAL_STRESS_FIELD || f == NON_ENRICHED_REAL_STRESS_FIELD /*|| f == PRINCIPAL_REAL_STRESS_FIELD */
           || f == EFFECTIVE_STRESS_FIELD || f == NON_ENRICHED_EFFECTIVE_STRESS_FIELD /*|| f == PRINCIPAL_EFFECTIVE_STRESS_FIELD*/  ;
}

bool isRealStressField ( FieldType f )
{
    return f == REAL_STRESS_FIELD || f == NON_ENRICHED_REAL_STRESS_FIELD || f == PRINCIPAL_REAL_STRESS_FIELD  ;
}

bool isEffectiveStressField ( FieldType f )
{
    return f == EFFECTIVE_STRESS_FIELD || f == NON_ENRICHED_EFFECTIVE_STRESS_FIELD || f == PRINCIPAL_EFFECTIVE_STRESS_FIELD  ;
}

Vector toPrincipal ( const Vector & stressOrStrain, CompositionType t )
{
    if (t ==   SINGLE_OFF_DIAGONAL_VALUES)
    {
        Vector ret ( 0., 2+ ( stressOrStrain.size() == 6 ) ) ;
        if ( ret.size() == 2 )
        {
            double trace = stressOrStrain[0] + stressOrStrain[1] ;
            double det = stressOrStrain[0]*stressOrStrain[1] - stressOrStrain[2]*stressOrStrain[2] ;
            double delta = std::sqrt(trace*trace - 4.*det) ;
            double angle =  0.5*std::atan2 ( stressOrStrain[2], stressOrStrain[0] - stressOrStrain[1] ) ;
            if(std::cos(angle) < 0)
            {
                ret[0] = (trace + delta)*.5 ;
                ret[1] = (trace - delta)*.5 ;
            }
            else
            {
                ret[0] = (trace - delta)*.5 ;
                ret[1] = (trace + delta)*.5 ;
            }
        }
        else if ( ret.size() == 3 )
        {
            Matrix mat = Amie::makeStressOrStrainMatrix ( stressOrStrain ) ;
            double trmat = -trace(mat);
            double detmat = -2.0*mat[0][1]*mat[0][2]*mat[1][2] + mat[0][0]*mat[1][2]*mat[1][2] + mat[1][1]*mat[0][2]*mat[0][2] + mat[2][2]*mat[0][1]*mat[0][1] - mat[0][0]*mat[1][1]*mat[2][2];
            double m2mat =  (mat[0][0]*mat[1][1] + mat[1][1]*mat[2][2] + mat[2][2]*mat[0][0]) - mat[0][2]*mat[0][2] - mat[0][1]*mat[0][1] - mat[1][2]*mat[1][2];
            double q =  m2mat/3. - trmat*trmat/9. ;
            double r = (trmat*m2mat - 3.*detmat)/6. - trmat*trmat*trmat/27. ;
            double d = q*q*q + r*r ;
            double r0 = std::pow(r*r - d, 1./6.);
            double phi = std::atan2(std::sqrt(-1.*d),r)/3. ;
            if ( std::abs(phi) < POINT_TOLERANCE )
            {
                phi = 0.  ;
            }
            if ( (phi) < 0. )
            {
                phi += M_PI  ;
            }
            double som = r0*std::cos(phi);
            double dif = r0*std::sin(phi) ;
            ret[0] = 2.*som - trmat/3. ;
            ret[1] = -som - trmat/3. - dif*std::sqrt(3.) ;
            ret[2] = -som - trmat/3. + dif*std::sqrt(3.) ;
        }
        return ret ;
    }
    else
    {
        Vector ret ( 0., 2+ ( stressOrStrain.size() == 6 ) ) ;
        if ( ret.size() == 2 )
        {
            double trace = stressOrStrain[0] + stressOrStrain[1] ;
            double det = stressOrStrain[0]*stressOrStrain[1] - 0.25*stressOrStrain[2]*stressOrStrain[2] ;
            double delta = std::sqrt(trace*trace - 4.*det) ;
            double angle =  0.5*std::atan2 (  0.5*stressOrStrain[2], stressOrStrain[0] - stressOrStrain[1] ) ;
            if(std::cos(angle) < 0)
            {
                ret[0] = (trace + delta)*.5 ;
                ret[1] = (trace - delta)*.5 ;
            }
            else
            {
                ret[0] = (trace - delta)*.5 ;
                ret[1] = (trace + delta)*.5 ;
            }
        }
        else if ( ret.size() == 3 )
        {
            Matrix mat = Amie::makeStressOrStrainMatrix ( stressOrStrain ) ;
            double trmat = -1.*trace(mat);
            double detmat = -2.0*mat[0][1]*mat[0][2]*mat[1][2]*0.125 + mat[0][0]*mat[1][2]*mat[1][2]*0.25 + mat[1][1]*mat[0][2]*mat[0][2]*0.25 + 0.25*mat[2][2]*mat[0][1]*mat[0][1] - mat[0][0]*mat[1][1]*mat[2][2];
            double m2mat =  (mat[0][0]*mat[1][1] + mat[1][1]*mat[2][2] + mat[2][2]*mat[0][0]) - 0.25*mat[0][2]*mat[0][2] - 0.25*mat[0][1]*mat[0][1] - 0.25*mat[1][2]*mat[1][2];
            double q =  m2mat/3. - trmat*trmat/9. ;
            double r = (trmat*m2mat - 3.*detmat)/6. - trmat*trmat*trmat/27. ;
            double d = q*q*q + r*r ;
            double r0 = std::pow(r*r - d, 1./6.);
            double phi = std::atan2(std::sqrt(-1.*d),r)/3. ;
            if ( std::abs(phi) < POINT_TOLERANCE )
            {
                phi = 0.  ;
            }
            if ( (phi) < 0. )
            {
                phi += M_PI  ;
            }
            double som = r0*std::cos(phi);
            double dif = r0*std::sin(phi) ;
            ret[0] = 2.*som - trmat/3. ;
            ret[1] = -som - trmat/3. - dif*std::sqrt(3.) ;
            ret[2] = -som - trmat/3. + dif*std::sqrt(3.) ;
        }
        return ret ;
    }

}

void ElementState::getExternalField ( Vector & nodalValues, int externaldofs, const Point & p, Vector & ret, bool local, VirtualMachine * vm ) 
{
    bool cleanup = !vm ;
    if ( !vm )
    {
        vm = new VirtualMachine() ;
    }
    Point p_ = p ;
    if ( !local )
    {
        p_ = parent->inLocalCoordinates ( p ) ;
    }

    ret = 0. ;

    for ( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
    {
        double f = vm->eval ( parent->getShapeFunction ( j ), p_ ) ;
        for ( int k = 0 ; k < externaldofs ; k++ )
        {
            ret[k] += f * nodalValues[ j*externaldofs + k ] ;
        }
    }
    if ( cleanup )
    {
        delete vm ;
    }
}

void ElementState::getExternalFieldAtGaussPoints ( Vector & nodalValues, int externaldofs, std::vector<Vector> & ret, VirtualMachine * vm ) 
{
    for ( size_t p = 0 ; p < parent->getGaussPoints().gaussPoints.size() ; p++ )
    {
        getExternalField ( nodalValues, externaldofs, parent->getGaussPoints().gaussPoints[p].first, ret[p], true, vm ) ;
    }
}

//helper function: be a bit faster
// void getStrainField ( const Point & p, Matrix * JinvCache, IntegrableEntity * parent,  const Vector & displacements, const Vector & enrichedDisplacements,  Vector & ret, VirtualMachine * vm ) 
// {
//   if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
//         {
//             double x_xi = 0;
//             double x_eta = 0;
//             double y_xi = 0;
//             double y_eta = 0;
// 
//             for ( size_t j = 0 ; j < parent->getShapeFunctions().size(); j++ )
//             {
//                 if(j*2 >= displacements.size())
//                 {
//                     std::cerr << "displacement size mismatch" << std::endl ;
//                     break ;
//                 }
//                 double f_xi  = vm ->deval ( parent->getShapeFunction ( j ), XI , p ) ;
//                 double f_eta = vm ->deval ( parent->getShapeFunction ( j ), ETA, p ) ;
//                 
//                 x_xi  += f_xi * displacements[j * 2] ;
//                 x_eta += f_eta * displacements[j * 2] ;
//                 y_xi  += f_xi * displacements[j * 2 + 1] ;
//                 y_eta += f_eta * displacements[j * 2 + 1] ;
//             }
//             for ( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
//             {
//                 double f_xi  = vm ->deval ( parent->getEnrichmentFunction ( j ), XI , p ) ;
//                 double f_eta = vm ->deval ( parent->getEnrichmentFunction ( j ), ETA, p ) ;
// 
//                 x_xi  += f_xi * enrichedDisplacements[j * 2] ;
//                 x_eta += f_eta * enrichedDisplacements[j * 2] ;
//                 y_xi  += f_xi * enrichedDisplacements[j * 2 +1] ;
//                 y_eta += f_eta * enrichedDisplacements[j * 2 + 1] ;
//             }
//             
//             ret[0] = ( x_xi ) * (*JinvCache)[0][0] + ( x_eta ) * (*JinvCache)[0][1] ;
//             ret[1] = ( y_xi ) * (*JinvCache)[1][0] + ( y_eta ) * (*JinvCache)[1][1] ;
//             ret[2] = ( ( x_xi ) * (*JinvCache)[1][0] + ( x_eta ) * (*JinvCache)[1][1]  + ( y_xi ) * (*JinvCache)[0][0] + ( y_eta ) * (*JinvCache)[0][1] );
//         }
//         else if ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3 )
//         {
//             double x_xi = 0;
//             double x_eta = 0;
//             double x_zeta = 0;
//             double y_xi = 0;
//             double y_eta = 0;
//             double y_zeta = 0;
//             double z_xi = 0;
//             double z_eta = 0;
//             double z_zeta = 0;
// 
//             for ( size_t j = 0 ; j < parent->getShapeFunctions().size() ; j++ )
//             {
//                 double f_xi = vm ->deval ( parent->getShapeFunction ( j ), XI, *p_ ) ;
//                 double f_eta = vm ->deval ( parent->getShapeFunction ( j ), ETA, *p_ ) ;
//                 double f_zeta = vm ->deval ( parent->getShapeFunction ( j ), ZETA, *p_ ) ;
//                 double x = displacements[j * 3] ;
//                 double y = displacements[j * 3 + 1] ;
//                 double z = displacements[j * 3 + 2] ;
// 
//                 x_xi   += f_xi   * x ;
//                 x_eta  += f_eta  * x ;
//                 x_zeta += f_zeta * x ;
//                 y_xi   += f_xi   * y ;
//                 y_eta  += f_eta  * y ;
//                 y_zeta += f_zeta * y ;
//                 z_xi   += f_xi   * z ;
//                 z_eta  += f_eta  * z ;
//                 z_zeta += f_zeta * z ;
//             }
// 
//             for ( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
//             {
//                 double f_xi = vm ->deval ( parent->getEnrichmentFunction ( j ), XI, *p_ ) ;
//                 double f_eta = vm ->deval ( parent->getEnrichmentFunction ( j ), ETA, *p_ ) ;
//                 double f_zeta = vm ->deval ( parent->getEnrichmentFunction ( j ), ZETA, *p_ ) ;
//                 double x = enrichedDisplacements[j * 3] ;
//                 double y = enrichedDisplacements[j * 3 + 1] ;
//                 double z = enrichedDisplacements[j * 3 + 2] ;
// 
//                 x_xi += f_xi * x;
//                 x_eta += f_eta * x ;
//                 x_zeta += f_zeta * x ;
//                 y_xi += f_xi * y ;
//                 y_eta += f_eta * y ;
//                 y_zeta += f_zeta * y ;
//                 z_xi += f_xi * z ;
//                 z_eta += f_eta * z ;
//                 z_zeta += f_zeta * z ;
//             }
// 
//             ret[0] = ( x_xi ) * (*JinvCache)[0][0] + ( x_eta ) * (*JinvCache)[0][1]  + ( x_zeta ) * (*JinvCache)[0][2];
//             ret[1] = ( y_xi ) * (*JinvCache)[1][0] + ( y_eta ) * (*JinvCache)[1][1]  + ( y_zeta ) * (*JinvCache)[1][2];
//             ret[2] = ( z_xi ) * (*JinvCache)[2][0] + ( z_eta ) * (*JinvCache)[2][1]  + ( z_zeta ) * (*JinvCache)[2][2];
// 
//             ret[3] = ( ( y_xi ) * (*JinvCache)[2][0] +
//                        ( y_eta ) * (*JinvCache)[2][1] +
//                        ( y_zeta ) * (*JinvCache)[2][2] +
//                        ( z_xi ) * (*JinvCache)[1][0] +
//                        ( z_eta ) * (*JinvCache)[1][1] +
//                        ( z_zeta ) * (*JinvCache)[1][2] );
// 
//             ret[4] = ( ( x_xi ) * (*JinvCache)[2][0] +
//                        ( x_eta ) * (*JinvCache)[2][1] +
//                        ( x_zeta ) * (*JinvCache)[2][2] +
//                        ( z_xi ) * (*JinvCache)[0][0] +
//                        ( z_eta ) * (*JinvCache)[0][1] +
//                        ( z_zeta ) * (*JinvCache)[0][2] );
// 
//             ret[5] = ( ( y_xi )   * (*JinvCache)[0][0] +
//                        ( y_eta )  * (*JinvCache)[0][1] +
//                        ( y_zeta ) * (*JinvCache)[0][2] +
//                        ( x_xi )   * (*JinvCache)[1][0] +
//                        ( x_eta )  * (*JinvCache)[1][1] +
//                        ( x_zeta ) * (*JinvCache)[1][2] );
//         }
// }

void ElementState::getField ( FieldType f, const Point & p, Vector & ret, bool local, VirtualMachine * vm, int ) 
{
    ret = 0 ;
    if(!parent->getBehaviour() || parent->getBehaviour()->type == VOID_BEHAVIOUR)
        return ;
    bool cleanup = !vm ;
    if ( !vm )
    {
        vm = new VirtualMachine() ;
    }
    int n = 0 ;
    bool cleanupp = !local ;
    const Point * p_ = &p ;
    if ( !local )
    {
        p_ = new Point(parent->inLocalCoordinates ( p ) );
    }

    switch ( f )
    {
    case SCALAR_DAMAGE_FIELD:
    {
        if(parent->getBehaviour()->getDamageModel())
            ret[0] = parent->getBehaviour()->getDamageModel()->getState().max() ;
        else
            ret[0] = 0. ;
        if ( cleanup )
        {
            delete vm ;
        }
        if(cleanupp)
            delete p_ ;
        return ;
    }
    case TENSOR_DAMAGE_FIELD:
    {
        if(parent->getBehaviour()->getDamageModel())
        {
            Vector d = parent->getBehaviour()->getDamageModel()->getState() ;
            ret = 0. ;
            for(size_t i = 0 ; i < std::min( ret.size(), d.size() ) ; i++)
               ret[i] = d[i] ;
        }
        else
            ret = 0. ;
        if ( cleanup )
        {
            delete vm ;
        }
        if(cleanupp)
            delete p_ ;
        return ;
    }
    case DISPLACEMENT_FIELD:
    {
        n =  parent->getBehaviour()->getNumberOfDegreesOfFreedom() ;
        for ( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
        {
            double f =  vm ->eval ( parent->getShapeFunction ( j ) , p_ ) ;
            for ( int k = 0 ; k < n ; k++ )
            {
                ret[k] += f * displacements[j*n+k] ;
            }
        }
        for ( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
        {
            double f =  vm ->eval ( parent->getEnrichmentFunction ( j ) , p_ ) ;
            for ( int k = 0 ; k < n ; k++ )
            {
                ret[k] += f * enrichedDisplacements[j*n+k] ;
            }
        }
        if ( cleanup )
        {
            delete vm ;
        }
        if(cleanupp)
            delete p_ ;
        return ;
    }
    case ENRICHED_DISPLACEMENT_FIELD:
    {
        n =  parent->getBehaviour()->getNumberOfDegreesOfFreedom() ;
        for ( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
        {
            double f =  vm ->eval ( parent->getEnrichmentFunction ( j ) , p_ ) ;
            for ( int k = 0 ; k < n ; k++ )
            {
                ret[k] += f * enrichedDisplacements[j*n+k] ;
            }
        }
        if(cleanupp)
            delete p_ ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case MECHANICAL_STRAIN_FIELD:
    {
        getField( TOTAL_STRAIN_FIELD, *p_, ret, true, vm) ;
        if(parent->getBehaviour() && parent->getBehaviour()->hasInducedForces())
            ret -= parent->getBehaviour()->getImposedStrain(*p_) ;
        if ( cleanup )
        {
            delete vm ;
        }
        if(cleanupp)
            delete p_ ;
        return ;
    }
    case TOTAL_STRAIN_FIELD:
    {
        if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
        {
            double x_xi = 0;
            double x_eta = 0;
            double y_xi = 0;
            double y_eta = 0;

            for ( size_t j = 0 ; j < parent->getShapeFunctions().size(); j++ )
            {
                if(j*2 >= displacements.size())
                {
                    std::cerr << "displacement size mismatch" << std::endl ;
                    break ;
                }
                double f_xi  = vm ->deval ( parent->getShapeFunction ( j ), XI , *p_ ) ;
                double f_eta = vm ->deval ( parent->getShapeFunction ( j ), ETA, *p_ ) ;
                
                x_xi  += f_xi * displacements[j * 2] ;
                x_eta += f_eta * displacements[j * 2] ;
                y_xi  += f_xi * displacements[j * 2 + 1] ;
                y_eta += f_eta * displacements[j * 2 + 1] ;
            }
            for ( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
            {
                double f_xi  = vm ->deval ( parent->getEnrichmentFunction ( j ), XI , *p_ ) ;
                double f_eta = vm ->deval ( parent->getEnrichmentFunction ( j ), ETA, *p_ ) ;

                x_xi  += f_xi * enrichedDisplacements[j * 2] ;
                x_eta += f_eta * enrichedDisplacements[j * 2] ;
                y_xi  += f_xi * enrichedDisplacements[j * 2 +1] ;
                y_eta += f_eta * enrichedDisplacements[j * 2 + 1] ;
            }

            if(!JinvCache || parent->isMoved())
            {
                if(JinvCache)
                    delete JinvCache ;
                JinvCache = new Matrix ( parent->spaceDimensions()+(parent->timePlanes()>1),parent->spaceDimensions()+(parent->timePlanes()>1)) ;
                parent->getInverseJacobianMatrix ( *p_, (*JinvCache) ) ;
            }
            
            ret[0] = ( x_xi ) * (*JinvCache)[0][0] + ( x_eta ) * (*JinvCache)[0][1] ;
            ret[1] = ( y_xi ) * (*JinvCache)[1][0] + ( y_eta ) * (*JinvCache)[1][1] ;
            ret[2] = ( ( x_xi ) * (*JinvCache)[1][0] + ( x_eta ) * (*JinvCache)[1][1]  + ( y_xi ) * (*JinvCache)[0][0] + ( y_eta ) * (*JinvCache)[0][1] );
        }
        else if ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3 )
        {
            double x_xi = 0;
            double x_eta = 0;
            double x_zeta = 0;
            double y_xi = 0;
            double y_eta = 0;
            double y_zeta = 0;
            double z_xi = 0;
            double z_eta = 0;
            double z_zeta = 0;

            for ( size_t j = 0 ; j < parent->getShapeFunctions().size() ; j++ )
            {
                double f_xi = vm ->deval ( parent->getShapeFunction ( j ), XI, *p_ ) ;
                double f_eta = vm ->deval ( parent->getShapeFunction ( j ), ETA, *p_ ) ;
                double f_zeta = vm ->deval ( parent->getShapeFunction ( j ), ZETA, *p_ ) ;
                double x = displacements[j * 3] ;
                double y = displacements[j * 3 + 1] ;
                double z = displacements[j * 3 + 2] ;

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

            for ( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
            {
                double f_xi = vm ->deval ( parent->getEnrichmentFunction ( j ), XI, *p_ ) ;
                double f_eta = vm ->deval ( parent->getEnrichmentFunction ( j ), ETA, *p_ ) ;
                double f_zeta = vm ->deval ( parent->getEnrichmentFunction ( j ), ZETA, *p_ ) ;
                double x = enrichedDisplacements[j * 3] ;
                double y = enrichedDisplacements[j * 3 + 1] ;
                double z = enrichedDisplacements[j * 3 + 2] ;

                x_xi += f_xi * x;
                x_eta += f_eta * x ;
                x_zeta += f_zeta * x ;
                y_xi += f_xi * y ;
                y_eta += f_eta * y ;
                y_zeta += f_zeta * y ;
                z_xi += f_xi * z ;
                z_eta += f_eta * z ;
                z_zeta += f_zeta * z ;
            }

            if(!JinvCache || parent->isMoved())
            {
                if(JinvCache)
                    delete JinvCache ;
                JinvCache = new Matrix ( parent->spaceDimensions()+(parent->timePlanes()>1),parent->spaceDimensions()+(parent->timePlanes()>1)) ;
                parent->getInverseJacobianMatrix ( *p_, (*JinvCache) ) ;
            }

            ret[0] = ( x_xi ) * (*JinvCache)[0][0] + ( x_eta ) * (*JinvCache)[0][1]  + ( x_zeta ) * (*JinvCache)[0][2];
            ret[1] = ( y_xi ) * (*JinvCache)[1][0] + ( y_eta ) * (*JinvCache)[1][1]  + ( y_zeta ) * (*JinvCache)[1][2];
            ret[2] = ( z_xi ) * (*JinvCache)[2][0] + ( z_eta ) * (*JinvCache)[2][1]  + ( z_zeta ) * (*JinvCache)[2][2];

            ret[3] = ( ( y_xi ) * (*JinvCache)[2][0] +
                       ( y_eta ) * (*JinvCache)[2][1] +
                       ( y_zeta ) * (*JinvCache)[2][2] +
                       ( z_xi ) * (*JinvCache)[1][0] +
                       ( z_eta ) * (*JinvCache)[1][1] +
                       ( z_zeta ) * (*JinvCache)[1][2] );

            ret[4] = ( ( x_xi ) * (*JinvCache)[2][0] +
                       ( x_eta ) * (*JinvCache)[2][1] +
                       ( x_zeta ) * (*JinvCache)[2][2] +
                       ( z_xi ) * (*JinvCache)[0][0] +
                       ( z_eta ) * (*JinvCache)[0][1] +
                       ( z_zeta ) * (*JinvCache)[0][2] );

            ret[5] = ( ( y_xi )   * (*JinvCache)[0][0] +
                       ( y_eta )  * (*JinvCache)[0][1] +
                       ( y_zeta ) * (*JinvCache)[0][2] +
                       ( x_xi )   * (*JinvCache)[1][0] +
                       ( x_eta )  * (*JinvCache)[1][1] +
                       ( x_zeta ) * (*JinvCache)[1][2] );
        }
        if ( cleanup )
        {
            delete vm ;
        }
        if(cleanupp)
            delete p_ ;
        return ;
    }
    case PRINCIPAL_STRAIN_FIELD:
    {
        Vector strains ( 0.,(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL) ? 6 : 3 ) ;
        getField ( MECHANICAL_STRAIN_FIELD, *p_, strains, true,vm ) ;

        ret = toPrincipal ( strains, DOUBLE_OFF_DIAGONAL_VALUES  ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        if(cleanupp)
            delete p_ ;
        return ;
    }
    case PRINCIPAL_MECHANICAL_STRAIN_FIELD:
    {
        Vector strains ( 0.,(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL) ? 6 : 3 ) ;
        getField ( MECHANICAL_STRAIN_FIELD, *p_, strains, true,vm ) ;

        ret = toPrincipal ( strains, DOUBLE_OFF_DIAGONAL_VALUES ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        if(cleanupp)
            delete p_ ;
        return ;
    }
    case NON_ENRICHED_STRAIN_FIELD:
    {
        if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2 )
        {
            double x_xi = 0;
            double x_eta = 0;
            double y_xi = 0;
            double y_eta = 0;

            if ( displacements.size() /2 == parent->getBoundingPoints().size() )
            {
                for ( size_t j = 0 ; j < parent->getBoundingPoints().size(); j++ )
                {
                    double f_xi = vm ->deval ( parent->getShapeFunction ( j ), XI, *p_ ) ;
                    double f_eta = vm ->deval ( parent->getShapeFunction ( j ), ETA, *p_ ) ;
                    x_xi += f_xi * displacements[j * 2] ;
                    x_eta += f_eta * displacements[j * 2] ;
                    y_xi += f_xi * displacements[j * 2 + 1] ;
                    y_eta += f_eta * displacements[j * 2 + 1] ;
                }
            }

            if(!JinvCache || parent->isMoved())
            {
                if(JinvCache)
                    delete JinvCache ;
                JinvCache = new Matrix ( parent->spaceDimensions()+(parent->timePlanes()>1),parent->spaceDimensions()+(parent->timePlanes()>1)) ;
                parent->getInverseJacobianMatrix ( *p_, (*JinvCache) ) ;
            }
            ret[0] = ( x_xi ) * (*JinvCache)[0][0] + ( x_eta ) * (*JinvCache)[0][1] ;
            ret[1] = ( y_xi ) * (*JinvCache)[1][0] + ( y_eta ) * (*JinvCache)[1][1] ;
            ret[2] = ( ( x_xi ) * (*JinvCache)[1][0] + ( x_eta ) * (*JinvCache)[1][1]  + ( y_xi ) * (*JinvCache)[0][0] + ( y_eta ) * (*JinvCache)[0][1] );
        }
        else if ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
        {
            double x_xi = 0;
            double x_eta = 0;
            double x_zeta = 0;
            double y_xi = 0;
            double y_eta = 0;
            double y_zeta = 0;
            double z_xi = 0;
            double z_eta = 0;
            double z_zeta = 0;

            for ( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
            {
                double f_xi = vm ->deval ( parent->getShapeFunction ( j ), XI, *p_ ) ;
                double f_eta = vm ->deval ( parent->getShapeFunction ( j ), ETA, *p_ ) ;
                double f_zeta = vm ->deval ( parent->getShapeFunction ( j ), ZETA, *p_ ) ;
                double x = displacements[j * 3] ;
                double y = displacements[j * 3 + 1] ;
                double z = displacements[j * 3 + 2] ;

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


            if(!JinvCache || parent->isMoved())
            {
                if(JinvCache)
                    delete JinvCache ;
                JinvCache = new Matrix ( parent->spaceDimensions()+(parent->timePlanes()>1),parent->spaceDimensions()+(parent->timePlanes()>1)) ;
                parent->getInverseJacobianMatrix ( *p_, (*JinvCache) ) ;
            }
            ret[0] = ( x_xi ) * (*JinvCache)[0][0] + ( x_eta ) * (*JinvCache)[0][1]  + ( x_zeta ) * (*JinvCache)[0][2];
            ret[1] = ( y_xi ) * (*JinvCache)[1][0] + ( y_eta ) * (*JinvCache)[1][1]  + ( y_zeta ) * (*JinvCache)[1][2];
            ret[2] = ( z_xi ) * (*JinvCache)[2][0] + ( z_eta ) * (*JinvCache)[2][1]  + ( z_zeta ) * (*JinvCache)[2][2];

            ret[3] = ( ( y_xi ) * (*JinvCache)[2][0] +
                       ( y_eta ) * (*JinvCache)[2][1] +
                       ( y_zeta ) * (*JinvCache)[2][2] +
                       ( z_xi ) * (*JinvCache)[1][0] +
                       ( z_eta ) * (*JinvCache)[1][1] +
                       ( z_zeta ) * (*JinvCache)[1][2] );

            ret[4] = ( ( x_xi ) * (*JinvCache)[2][0] +
                       ( x_eta ) * (*JinvCache)[2][1] +
                       ( x_zeta ) * (*JinvCache)[2][2] +
                       ( z_xi ) * (*JinvCache)[0][0] +
                       ( z_eta ) * (*JinvCache)[0][1] +
                       ( z_zeta ) * (*JinvCache)[0][2] );

            ret[5] = ( ( y_xi )   * (*JinvCache)[0][0] +
                       ( y_eta )  * (*JinvCache)[0][1] +
                       ( y_zeta ) * (*JinvCache)[0][2] +
                       ( x_xi )   * (*JinvCache)[1][0] +
                       ( x_eta )  * (*JinvCache)[1][1] +
                       ( x_zeta ) * (*JinvCache)[1][2] );
        }
        if ( cleanup )
        {
            delete vm ;
        }
        if(cleanupp)
            delete p_ ;
        return ;
    }
    case VON_MISES_STRAIN_FIELD:
    {
        Vector eps ( 0., ( size_t ) parent->spaceDimensions() ) ;
        getField ( PRINCIPAL_STRAIN_FIELD, *p_, eps, true,vm ) ;
        if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
        {
            ret[0] = ( 2. / 3. * ( eps[0] * eps[0] + eps[1] * eps[1] ) ) ;
        }
        else if ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
        {
            ret[0] = sqrt ( 2. / 3. * ( eps[0] * eps[0] + eps[1] * eps[1] + eps[2] * eps[2] ) ) ;
        }
        if ( cleanup )
        {
            delete vm ;
        }
        if(cleanupp)
            delete p_ ;
        return ;
    }
    case REAL_STRESS_FIELD:
    {
        getField ( MECHANICAL_STRAIN_FIELD, *p_, ret, true, vm ) ;
        if ( parent->getBehaviour()->getTensor ( *p_, parent ).numCols() != ret.size() )
        {
            ret = 0 ;
            if ( cleanup )
            {
                delete vm ;
            }
            if(cleanupp)
                delete p_ ;
            return ;
        }
        
        ret = ( Vector ) ( parent->getBehaviour()->getTensor ( *p_, parent ) * ret ) - parent->getBehaviour()->getImposedStress ( *p_, parent ) ;

        if ( cleanup )
        {
            delete vm ;
        }
        if(cleanupp)
            delete p_ ;
        return ;
    }
    case IMPOSED_STRESS_FIELD:
    {
        if ( parent->getBehaviour()->getTensor ( *p_, parent ).numCols() != ret.size() )
        {
            ret = 0 ;
            if ( cleanup )
            {
                delete vm ;
            }
            if(cleanupp)
                delete p_ ;
            return ;
        }

        ret = parent->getBehaviour()->getImposedStress ( *p_, parent ) ;

        if ( cleanup )
        {
            delete vm ;
        }
        if(cleanupp)
            delete p_ ;
        return ;
    }
    case IMPOSED_STRAIN_FIELD:
    {
        if ( parent->getBehaviour()->getTensor ( *p_, parent ).numCols() != ret.size() )
        {
            ret = 0 ;
            if ( cleanup )
            {
                delete vm ;
            }
            if(cleanupp)
                delete p_ ;
            return ;
        }

        ret = parent->getBehaviour()->getImposedStrain ( *p_, parent ) ;
	

        if ( cleanup )
        {
            delete vm ;
        }
        if(cleanupp)
            delete p_ ;
        return ;
    }
    case PRINCIPAL_IMPOSED_STRESS_FIELD:
    {

        Vector str = parent->getBehaviour()->getImposedStress ( *p_, parent ) ;
        ret = toPrincipal(str,SINGLE_OFF_DIAGONAL_VALUES) ;

        if ( cleanup )
        {
            delete vm ;
        }
        if(cleanupp)
            delete p_ ;
        return ;
    }
    case PRINCIPAL_IMPOSED_STRAIN_FIELD:
    {

        Vector str = parent->getBehaviour()->getImposedStrain ( *p_, parent ) ;
        ret = toPrincipal(str, DOUBLE_OFF_DIAGONAL_VALUES) ;

        if ( cleanup )
        {
            delete vm ;
        }
        if(cleanupp)
            delete p_ ;
        return ;
    }
    case PRINCIPAL_REAL_STRESS_FIELD:
    {
        Vector stress ( 0.,3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ;
        getField ( REAL_STRESS_FIELD, *p_, stress, true,vm ) ;
        ret = toPrincipal ( stress , SINGLE_OFF_DIAGONAL_VALUES) ;
        if ( cleanup )
        {
            delete vm ;
        }
        if(cleanupp)
            delete p_ ;
        return ;
    }
    case NON_ENRICHED_REAL_STRESS_FIELD:
    {
        getField ( NON_ENRICHED_STRAIN_FIELD, *p_, ret, true,vm ) ;
        if ( parent->getBehaviour()->getTensor ( *p_, parent ).numCols() != ret.size() )
        {
            ret = 0 ;
            if ( cleanup )
            {
                delete vm ;
            }
            if(cleanupp)
                delete p_ ;
            return ;
        }
        ret = ( Vector ) ( parent->getBehaviour()->getTensor ( *p_, parent ) * ret ) - getParent()->getBehaviour()->getImposedStress ( *p_, parent ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        if(cleanupp)
            delete p_ ;
        return ;
    }
    case VON_MISES_REAL_STRESS_FIELD:
    {
        if ( parent->getOrder() == LINEAR )
        {
            if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
            {
                Vector sigma ( 0., 2 ) ;
                Point c ( 1./3., 1./3. ) ;
                getField ( PRINCIPAL_REAL_STRESS_FIELD, c, sigma, true, vm ) ;
                ret[0] = sqrt ( ( ( sigma[0] - sigma[1] ) * ( sigma[0] - sigma[1] ) + sigma[0] * sigma[0] + sigma[1] * sigma[1] ) / 2. ) ;
                if ( cleanup )
                {
                    delete vm ;
                }
                if(cleanupp)
                    delete p_ ;
                return ;
            }
            else if ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
            {
                Amie::PointArray pts ( 4 ) ;
                pts[0] = &parent->getBoundingPoint ( 0 ) ;
                pts[1] = &parent->getBoundingPoint ( 1 ) ;
                pts[2] = &parent->getBoundingPoint ( 2 ) ;
                pts[3] = &parent->getBoundingPoint ( 3 ) ;
                Vector sigma ( 0., 24 ) ;
                getField ( PRINCIPAL_REAL_STRESS_FIELD, pts, sigma, false, vm ) ;
                sigma[0] = sqrt ( ( sigma[0] - sigma[1] ) * ( sigma[0] - sigma[1] ) + ( sigma[0] - sigma[2] ) * ( sigma[0] - sigma[2] ) + ( sigma[1] - sigma[2] ) * ( sigma[1] - sigma[2] ) ) / 6. ;
                sigma[1] = sqrt ( ( sigma[6] - sigma[7] ) * ( sigma[6] - sigma[7] ) + ( sigma[6] - sigma[8] ) * ( sigma[6] - sigma[8] ) + ( sigma[7] - sigma[8] ) * ( sigma[7] - sigma[8] ) ) / 6. ;
                sigma[2] = sqrt ( ( sigma[12] - sigma[13] ) * ( sigma[12] - sigma[13] ) + ( sigma[12] - sigma[14] ) * ( sigma[12] - sigma[14] ) + ( sigma[13] - sigma[14] ) * ( sigma[13] - sigma[14] ) ) / 6. ;
                sigma[3] = sqrt ( ( sigma[18] - sigma[19] ) * ( sigma[18] - sigma[19] ) + ( sigma[18] - sigma[20] ) * ( sigma[18] - sigma[20] ) + ( sigma[19] - sigma[20] ) * ( sigma[19] - sigma[20] ) ) / 6. ;

                ret[0] = std::max ( std::max ( std::max ( sigma[0], sigma[1] ), sigma[2] ), sigma[3] ) ;
                if ( cleanup )
                {
                    delete vm ;
                }
                if(cleanupp)
                    delete p_ ;
                return ;
            }
            return ;
        }
        else
        {
            if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
            {
                Vector principalStresses ( 0., parent->getBoundingPoints().size() *2 ) ;
                getField ( PRINCIPAL_REAL_STRESS_FIELD, parent->getBoundingPoints(), principalStresses, false, vm ) ;
                double maxS = 0 ;

                for ( size_t i = 0 ; i < principalStresses.size() / 2 ; i++ )
                {
                    maxS = std::max ( maxS,
                                      sqrt ( ( ( principalStresses[i * 2 + 0] - principalStresses[i * 2 + 1] ) * ( principalStresses[i * 2 + 0] - principalStresses[i * 2 + 1] ) + principalStresses[i * 2 + 0] * principalStresses[i * 2 + 0] + principalStresses[i * 2 + 1] * principalStresses[i * 2 + 1] ) / 2. ) ) ;
                }

                ret[0] = maxS ;
                if ( cleanup )
                {
                    delete vm ;
                }
                if(cleanupp)
                    delete p_ ;
                return ;

            }
            else if ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
            {
                Amie::PointArray pts ( 4 ) ;
                pts[0] = &parent->getBoundingPoint ( 0 ) ;
                pts[1] = &parent->getBoundingPoint ( 2 ) ;
                pts[2] = &parent->getBoundingPoint ( 4 ) ;
                pts[3] = &parent->getBoundingPoint ( 6 ) ;
                Vector sigma ( 0., 24 ) ;
                this->getField ( PRINCIPAL_REAL_STRESS_FIELD, pts, sigma, false, vm ) ;
                sigma[0] = sqrt ( ( sigma[0] - sigma[1] ) * ( sigma[0] - sigma[1] ) + ( sigma[0] - sigma[2] ) * ( sigma[0] - sigma[2] ) + ( sigma[1] - sigma[2] ) * ( sigma[1] - sigma[2] ) ) / 6 ;
                sigma[1] = sqrt ( ( sigma[6] - sigma[7] ) * ( sigma[6] - sigma[7] ) + ( sigma[6] - sigma[8] ) * ( sigma[6] - sigma[8] ) + ( sigma[7] - sigma[8] ) * ( sigma[7] - sigma[8] ) ) / 6 ;
                sigma[2] = sqrt ( ( sigma[12] - sigma[13] ) * ( sigma[12] - sigma[13] ) + ( sigma[12] - sigma[14] ) * ( sigma[12] - sigma[14] ) + ( sigma[13] - sigma[14] ) * ( sigma[13] - sigma[14] ) ) / 6 ;
                sigma[3] = sqrt ( ( sigma[18] - sigma[19] ) * ( sigma[18] - sigma[19] ) + ( sigma[18] - sigma[20] ) * ( sigma[18] - sigma[20] ) + ( sigma[19] - sigma[20] ) * ( sigma[19] - sigma[20] ) ) / 6 ;

                ret[0] = std::max ( std::max ( std::max ( sigma[0], sigma[1] ), sigma[2] ), sigma[3] ) ;
                if ( cleanup )
                {
                    delete vm ;
                }
                if(cleanupp)
                    delete p_ ;
                return ;
            }
        }
        return ;
    }
    case EFFECTIVE_STRESS_FIELD:
    {
        getField ( TOTAL_STRAIN_FIELD, *p_, ret, true,vm, 0 ) ;
        ret = ( Vector ) ( parent->getBehaviour()->param * ret ) - getParent()->getBehaviour()->getImposedStrain ( *p_, parent ) *parent->getBehaviour()->param ;
        if ( cleanup )
        {
            delete vm ;
        }
        if(cleanupp)
            delete p_ ;
        return ;
    }
    case PRINCIPAL_EFFECTIVE_STRESS_FIELD:
    {
        Vector stress ( 0.,3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ;
        getField ( EFFECTIVE_STRESS_FIELD, *p_, stress, true,vm ) ;
        ret = toPrincipal ( stress , SINGLE_OFF_DIAGONAL_VALUES) ;
        if ( cleanup )
        {
            delete vm ;
        }
        if(cleanupp)
            delete p_ ;
        return ;
    }
    case NON_ENRICHED_EFFECTIVE_STRESS_FIELD:
    {
        getField ( NON_ENRICHED_STRAIN_FIELD, *p_, ret, true,vm ) ;
        ret = ( Vector ) ( parent->getBehaviour()->param * ret ) - getParent()->getBehaviour()->getImposedStrain ( *p_, parent ) *parent->getBehaviour()->param ;
        if ( cleanup )
        {
            delete vm ;
        }
        if(cleanupp)
            delete p_ ;
        return ;
    }
    case VON_MISES_EFFECTIVE_STRESS_FIELD:
    {
        if ( parent->getOrder() == LINEAR )
        {
            if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
            {
                Vector sigma ( 0., 2 ) ;
                Point c ( 1./3., 1./3. ) ;
                getField ( PRINCIPAL_EFFECTIVE_STRESS_FIELD, c, sigma, true,vm ) ;
                ret[0] = sqrt ( ( ( sigma[0] - sigma[1] ) * ( sigma[0] - sigma[1] ) + sigma[0] * sigma[0] + sigma[1] * sigma[1] ) / 2. ) ;
                if ( cleanup )
                {
                    delete vm ;
                }
                if(cleanupp)
                    delete p_ ;
                return ;
            }
            else if ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
            {
                Amie::PointArray pts ( 4 ) ;
                pts[0] = &parent->getBoundingPoint ( 0 ) ;
                pts[1] = &parent->getBoundingPoint ( 1 ) ;
                pts[2] = &parent->getBoundingPoint ( 2 ) ;
                pts[3] = &parent->getBoundingPoint ( 3 ) ;
                Vector sigma ( 0., 24 ) ;
                getField ( PRINCIPAL_EFFECTIVE_STRESS_FIELD, pts, sigma, false, vm ) ;
                sigma[0] = sqrt ( ( sigma[0] - sigma[1] ) * ( sigma[0] - sigma[1] ) + ( sigma[0] - sigma[2] ) * ( sigma[0] - sigma[2] ) + ( sigma[1] - sigma[2] ) * ( sigma[1] - sigma[2] ) ) / 6. ;
                sigma[1] = sqrt ( ( sigma[6] - sigma[7] ) * ( sigma[6] - sigma[7] ) + ( sigma[6] - sigma[8] ) * ( sigma[6] - sigma[8] ) + ( sigma[7] - sigma[8] ) * ( sigma[7] - sigma[8] ) ) / 6. ;
                sigma[2] = sqrt ( ( sigma[12] - sigma[13] ) * ( sigma[12] - sigma[13] ) + ( sigma[12] - sigma[14] ) * ( sigma[12] - sigma[14] ) + ( sigma[13] - sigma[14] ) * ( sigma[13] - sigma[14] ) ) / 6. ;
                sigma[3] = sqrt ( ( sigma[18] - sigma[19] ) * ( sigma[18] - sigma[19] ) + ( sigma[18] - sigma[20] ) * ( sigma[18] - sigma[20] ) + ( sigma[19] - sigma[20] ) * ( sigma[19] - sigma[20] ) ) / 6. ;

                ret[0] = std::max ( std::max ( std::max ( sigma[0], sigma[1] ), sigma[2] ), sigma[3] ) ;
                if ( cleanup )
                {
                    delete vm ;
                }
                if(cleanupp)
                    delete p_ ;
                return ;
            }
            if ( cleanup )
            {
                delete vm ;
            }
            if(cleanupp)
                delete p_ ;
            return ;
        }
        else
        {
            if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
            {
                Vector principalStresses ( 0., parent->getBoundingPoints().size() *2 ) ;
                getField ( PRINCIPAL_EFFECTIVE_STRESS_FIELD, parent->getBoundingPoints(), principalStresses, false, vm ) ;
                double maxS = 0 ;

                for ( size_t i = 0 ; i < principalStresses.size() / 2 ; i++ )
                {
                    maxS = std::max ( maxS,
                                      sqrt ( ( ( principalStresses[i * 2 + 0] - principalStresses[i * 2 + 1] ) * ( principalStresses[i * 2 + 0] - principalStresses[i * 2 + 1] ) + principalStresses[i * 2 + 0] * principalStresses[i * 2 + 0] + principalStresses[i * 2 + 1] * principalStresses[i * 2 + 1] ) / 2. ) ) ;
                }

                ret[0] = maxS ;
                if ( cleanup )
                {
                    delete vm ;
                }
                if(cleanupp)
                    delete p_ ;
                return ;

            }
            else if ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
            {
                Amie::PointArray pts ( 4 ) ;
                pts[0] = &parent->getBoundingPoint ( 0 ) ;
                pts[1] = &parent->getBoundingPoint ( 2 ) ;
                pts[2] = &parent->getBoundingPoint ( 4 ) ;
                pts[3] = &parent->getBoundingPoint ( 6 ) ;
                Vector sigma ( 0., 24 ) ;
                getField ( PRINCIPAL_EFFECTIVE_STRESS_FIELD, pts, sigma, false, vm ) ;
                sigma[0] = sqrt ( ( sigma[0] - sigma[1] ) * ( sigma[0] - sigma[1] ) + ( sigma[0] - sigma[2] ) * ( sigma[0] - sigma[2] ) + ( sigma[1] - sigma[2] ) * ( sigma[1] - sigma[2] ) ) / 6 ;
                sigma[1] = sqrt ( ( sigma[6] - sigma[7] ) * ( sigma[6] - sigma[7] ) + ( sigma[6] - sigma[8] ) * ( sigma[6] - sigma[8] ) + ( sigma[7] - sigma[8] ) * ( sigma[7] - sigma[8] ) ) / 6 ;
                sigma[2] = sqrt ( ( sigma[12] - sigma[13] ) * ( sigma[12] - sigma[13] ) + ( sigma[12] - sigma[14] ) * ( sigma[12] - sigma[14] ) + ( sigma[13] - sigma[14] ) * ( sigma[13] - sigma[14] ) ) / 6 ;
                sigma[3] = sqrt ( ( sigma[18] - sigma[19] ) * ( sigma[18] - sigma[19] ) + ( sigma[18] - sigma[20] ) * ( sigma[18] - sigma[20] ) + ( sigma[19] - sigma[20] ) * ( sigma[19] - sigma[20] ) ) / 6 ;

                ret[0] = std::max ( std::max ( std::max ( sigma[0], sigma[1] ), sigma[2] ), sigma[3] ) ;
                if ( cleanup )
                {
                    delete vm ;
                }
                if(cleanupp)
                    delete p_ ;
                return ;
            }
        }
        if ( cleanup )
        {
            delete vm ;
        }
        if(cleanupp)
            delete p_ ;
        return ;
    }
    case PRINCIPAL_STRESS_ANGLE_FIELD:
    {
        Vector strains ( 0., 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ;
        getField ( REAL_STRESS_FIELD,  *p_, strains, true ) ;
        if ( std::abs ( strains ).max() < POINT_TOLERANCE )
        {
            ret = 0 ;
            if ( cleanup )
            {
                delete vm ;
            }
            if(cleanupp)
                delete p_ ;
            return ;
        }
        if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
        {
            ret[0] =  0.5*atan2 ( 0.5*strains[2], strains[0] - strains[1]  ) ;

        }
        else
        {
            ret[0] = 0.5*atan2 ( 0.5*strains[3], strains[0] - strains[1] ) ;
            ret[1] = 0.5*atan2 ( 0.5*strains[4], strains[0] - strains[2] ) ;
            ret[2] = 0.5*atan2 ( 0.5*strains[5], strains[1] - strains[2] ) ;
        }
        if ( cleanup )
        {
            delete vm ;
        }
        if(cleanupp)
            delete p_ ;
        return ;
    }
    case PRINCIPAL_STRAIN_ANGLE_FIELD:
    {
        Vector strains ( 0., 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ;
        getField ( MECHANICAL_STRAIN_FIELD,  *p_, strains, true ) ;
        if ( std::abs ( strains ).max() < POINT_TOLERANCE )
        {
            ret = 0 ;
            if ( cleanup )
            {
                delete vm ;
            }
            if(cleanupp)
                delete p_ ;
            return ;
        }
        if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
        {
            ret[0] =  0.5*atan2 ( 0.5*strains[2], strains[0] - strains[1]  ) ;
        }
        else
        {
            ret[0] = 0.5*atan2 ( 0.5*strains[3] , strains[0] - strains[1] ) ;
            ret[1] = 0.5*atan2 (  0.5*strains[4], strains[0] - strains[2] ) ;
            ret[2] = 0.5*atan2 (  0.5*strains[5], strains[1] - strains[2] ) ;
        }
        if ( cleanup )
        {
            delete vm ;
        }
        if(cleanupp)
            delete p_ ;
        return ;
    }
    case GRADIENT_FIELD:
    {
        if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 1 )
        {
            double x_xi = 0;
            double x_eta = 0;

            for ( size_t j = 0 ; j < parent->getBoundingPoints().size(); j++ )
            {
                double f_xi = vm ->deval ( parent->getShapeFunction ( j ), XI, *p_ ) ;
                double f_eta = vm ->deval ( parent->getShapeFunction ( j ), ETA, *p_ ) ;
                x_xi += f_xi * displacements[j] ;
                x_eta += f_eta * displacements[j] ;
            }

            for ( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() && j < enrichedDisplacements.size(); j++ )
            {
                double f_xi = vm ->deval ( parent->getEnrichmentFunction ( j ), XI, *p_ ) ;
                double f_eta = vm ->deval ( parent->getEnrichmentFunction ( j ), ETA, *p_ ) ;
                x_xi += f_xi * enrichedDisplacements[j] ;
                x_eta += f_eta * enrichedDisplacements[j] ;

            }

            if(!JinvCache || parent->isMoved())
            {
                if(JinvCache)
                    delete JinvCache ;
                JinvCache = new Matrix ( parent->spaceDimensions()+(parent->timePlanes()>1),parent->spaceDimensions()+(parent->timePlanes()>1)) ;
                parent->getInverseJacobianMatrix ( *p_, (*JinvCache) ) ;
            }
            
            ret[0] = ( x_xi ) * (*JinvCache)[0][0] + ( x_eta ) * (*JinvCache)[0][1] ;
            ret[1] = ( x_xi ) * (*JinvCache)[1][0] + ( x_eta ) * (*JinvCache)[1][1] ;
        }
        else if ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 1 )
        {
            double x_xi = 0;
            double x_eta = 0;
            double x_zeta = 0;

            for ( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
            {
                double f_xi = vm ->deval ( parent->getShapeFunction ( j ), XI, *p_ ) ;
                double f_eta = vm ->deval ( parent->getShapeFunction ( j ), ETA, *p_ ) ;
                double f_zeta = vm ->deval ( parent->getShapeFunction ( j ), ZETA, *p_ ) ;
                double x = displacements[j] ;

                x_xi   += f_xi   * x ;
                x_eta  += f_eta  * x ;
                x_zeta += f_zeta * x ;
            }

            for ( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
            {
                double f_xi = vm ->deval ( parent->getEnrichmentFunction ( j ), XI, *p_ ) ;
                double f_eta = vm ->deval ( parent->getEnrichmentFunction ( j ), ETA, *p_ ) ;
                double f_zeta = vm ->deval ( parent->getEnrichmentFunction ( j ), ZETA, *p_ ) ;
                double x = enrichedDisplacements[j] ;

                x_xi += f_xi * x;
                x_eta += f_eta * x ;
                x_zeta += f_zeta * x ;
            }

            if(!JinvCache || parent->isMoved())
            {
                if(JinvCache)
                    delete JinvCache ;
                JinvCache = new Matrix ( parent->spaceDimensions()+(parent->timePlanes()>1),parent->spaceDimensions()+(parent->timePlanes()>1)) ;
                parent->getInverseJacobianMatrix ( *p_, (*JinvCache) ) ;
            }
            ret[0] = ( x_xi ) * (*JinvCache)[0][0] + ( x_eta ) * (*JinvCache)[0][1]  + ( x_zeta ) * (*JinvCache)[0][2];
            ret[1] = ( x_xi ) * (*JinvCache)[1][0] + ( x_eta ) * (*JinvCache)[1][1]  + ( x_zeta ) * (*JinvCache)[1][2];
            ret[2] = ( x_xi ) * (*JinvCache)[2][0] + ( x_eta ) * (*JinvCache)[2][1]  + ( x_zeta ) * (*JinvCache)[2][2];
        }
        if ( cleanup )
        {
            delete vm ;
        }
        if(cleanupp)
            delete p_ ;
        return ;
    }
    case FLUX_FIELD:
    {
        getField ( GRADIENT_FIELD, *p_, ret, true, vm ) ;
        ret = ( Vector ) ( parent->getBehaviour()->getTensor ( *p_, parent ) * ret ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        if(cleanupp)
            delete p_ ;        
        return ;
    }
    default:
        std::cerr << "field type not supported" << std::endl ;
        exit(0) ;
    }
}

void ElementState::getField ( FieldType f, const PointArray & p, Vector & ret, bool local, VirtualMachine * vm, int ) 
{
    bool cleanup = !vm ;
    if ( !vm )
    {
        vm = new VirtualMachine() ;
    }
    Vector buffer ( 0., ret.size() /p.size() ) ;
    for ( size_t i = 0 ; i < p.size() ; i++ )
    {
        getField ( f, *p[i], buffer, local, vm ) ;
        for ( size_t j = buffer.size() *i ; j < buffer.size() * ( i+1 ) ; j++ )
        {
            ret[j] = buffer[j - buffer.size() *i] ;
        }
    }
    if ( cleanup )
    {
        delete vm ;
    }
}

void ElementState::getField ( FieldType f, const std::valarray<std::pair<Point, double> > & p, Vector & ret, bool local, VirtualMachine * vm, int ) 
{
    bool cleanup = !vm ;
    if ( !vm )
    {
        vm = new VirtualMachine() ;
    }
    Vector buffer ( 0., ret.size() /p.size() ) ;
    for ( size_t i = 0 ; i < p.size() ; i++ )
    {
        getField ( f, p[i].first, buffer, local, vm ) ;
        for ( size_t j = buffer.size() *i ; j < buffer.size() * ( i+1 ) ; j++ )
        {
            ret[j] = buffer[j - buffer.size() *i] ;
        }
    }
    if ( cleanup )
    {
        delete vm ;
    }
}

void ElementState::getFieldAtGaussPoint ( FieldType f, size_t p, Vector & ret, VirtualMachine * vm, int i )
{
    bool cleanup = !vm ;
    if ( !vm )
    {
        vm = new VirtualMachine() ;
    }
    Point p_ = parent->getGaussPoints().gaussPoints[p].first ;
    getField ( f, p_, ret, true, vm, i ) ;
    if ( cleanup )
    {
        delete vm ;
    }
}

double ElementState::getAverageField ( FieldType f, Vector & ret, VirtualMachine * vm, double time, const std::vector<double> & weights, int index )
{
    while(true) {
//         usleep(1) ;
        bool test = true;
        #pragma omp flush(lock)
        #pragma omp atomic read
        test = lock ;
        if(!test)
        {
            #pragma omp atomic write
            lock = true ;
            #pragma omp flush
            break ;
        }
    }

    bool cleanup = !vm ;
    if ( !vm )
    {
        vm = new VirtualMachine() ;
    }
    const GaussPointArray & gp = parent->getGaussPoints() ;
    size_t blocks = parent->getBehaviour()->getNumberOfDegreesOfFreedom() / parent->spaceDimensions() ;
    if( fieldTypeElementarySize ( f, parent->spaceDimensions(), blocks ) != ret.size())
        ret.resize( fieldTypeElementarySize ( f, parent->spaceDimensions(), blocks ) , 0.) ;
    else
        ret = 0 ;
//     double v = (parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? parent->area() : parent->volume() ;
    double total = 0 ;
    bool weighted = true ;
    
    if ( weights.size() != gp.gaussPoints.size() )
    {
        weighted = false ;
    }

    switch ( f )
    {
    case TOTAL_STRAIN_FIELD :

        if ( strainAtGaussPoints.size() != (parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? 3*gp.gaussPoints.size() : 6*gp.gaussPoints.size() )
        {
            strainAtGaussPoints.resize ( (parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? 3*gp.gaussPoints.size() : 6*gp.gaussPoints.size() , 0. ) ;
            strainAtGaussPointsSet = false ;
        }

        if ( !strainAtGaussPointsSet)
        {
            strainAtGaussPoints = 0 ;
            strainAtGaussPointsSet = true ;
            Vector tmp ( strainAtGaussPoints.size() /gp.gaussPoints.size() ) ;
            for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
            {
                
                getField ( f, gp.gaussPoints[i].first, tmp, true,vm, i ) ;
                for ( size_t j = 0 ; j < strainAtGaussPoints.size() /gp.gaussPoints.size() ; j++ )
                {
                    strainAtGaussPoints[i*strainAtGaussPoints.size() /gp.gaussPoints.size() +j] = tmp[j] ;
                }
                if(ret.size() != tmp.size())
                    ret.resize(tmp.size(), 0.);
                if(weighted)
                {
                    ret += tmp*gp.gaussPoints[i].second*weights[i] ;
                    total += gp.gaussPoints[i].second*weights[i] ;
                }
                else
                {
                    ret += tmp*gp.gaussPoints[i].second ;
                    total += gp.gaussPoints[i].second ;
                }

            }
            if ( total > ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? POINT_TOLERANCE*POINT_TOLERANCE : POINT_TOLERANCE*POINT_TOLERANCE*POINT_TOLERANCE))
            {
                ret /= total ;
            }
            else
            {
                ret = 0 ;
            }
            if ( cleanup )
            {
                delete vm ;
            }
        }
        else
        {
            Vector tmp ( strainAtGaussPoints.size() /gp.gaussPoints.size() ) ;
            for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
            { 
                for ( size_t j = 0 ; j < strainAtGaussPoints.size() /gp.gaussPoints.size() ; j++ )
                {
                    tmp[j] = strainAtGaussPoints[i*strainAtGaussPoints.size() /gp.gaussPoints.size() +j];
                }
                if(ret.size() != tmp.size())
                    ret.resize(tmp.size(), 0.);
                if(weighted)
                {
                    ret += tmp*gp.gaussPoints[i].second*weights[i] ;
                    total += gp.gaussPoints[i].second*weights[i] ;
                }
                else
                {
                    ret += tmp*gp.gaussPoints[i].second ;
                    total += gp.gaussPoints[i].second ;
                }
            }
            if ( total > ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? POINT_TOLERANCE*POINT_TOLERANCE : POINT_TOLERANCE*POINT_TOLERANCE*POINT_TOLERANCE) )
            {
                ret /= total ;
            }
            else
            {
                ret = 0 ;
            }
            if ( cleanup )
            {
                delete vm ;
            }
        }
        #pragma omp atomic write
        lock = false ;
        return total;
    case MECHANICAL_STRAIN_FIELD :

        if ( strainAtGaussPoints.size() != ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? 3*gp.gaussPoints.size() : 6*gp.gaussPoints.size()) )
        {
            strainAtGaussPoints.resize ( (parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? 3*gp.gaussPoints.size() : 6*gp.gaussPoints.size() , 0. ) ;
            strainAtGaussPointsSet = false ;
        }

        if ( !strainAtGaussPointsSet)
        {
            strainAtGaussPoints = 0 ;
            strainAtGaussPointsSet = true ;

            Vector tmp ( strainAtGaussPoints.size() /gp.gaussPoints.size() ) ;
            for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
            {
                
                getField ( MECHANICAL_STRAIN_FIELD, gp.gaussPoints[i].first, tmp, true,vm, i ) ;
                tmp -= getParent()->getBehaviour()->getImposedStrain(gp.gaussPoints[i].first) ;
		 
                for ( size_t j = 0 ; j < strainAtGaussPoints.size() /gp.gaussPoints.size() ; j++ )
                {
                    strainAtGaussPoints[i*strainAtGaussPoints.size() /gp.gaussPoints.size() +j] = tmp[j] ;
                }
                if(ret.size() != tmp.size())
                    ret.resize(tmp.size(), 0.);

                if(weighted)
                {
                    ret += tmp *gp.gaussPoints[i].second*weights[i] ;
                    total += gp.gaussPoints[i].second*weights[i] ;
                }
                else
                {
                    ret += tmp *gp.gaussPoints[i].second ;
                    total += gp.gaussPoints[i].second ;
                }

            }
            if ( total > ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? POINT_TOLERANCE*POINT_TOLERANCE : POINT_TOLERANCE*POINT_TOLERANCE*POINT_TOLERANCE) )
            {
                ret /= total ;
            }
            else
            {
                ret = 0 ;
            }
            if ( cleanup )
            {
                delete vm ;
            }
        }
        else
        {
            Vector tmp ( strainAtGaussPoints.size() /gp.gaussPoints.size() ) ;
            for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
            {

                for ( size_t j = 0 ; j < strainAtGaussPoints.size() /gp.gaussPoints.size() ; j++ )
                {
                    tmp[j] = strainAtGaussPoints[i*strainAtGaussPoints.size() /gp.gaussPoints.size() +j];
                }
                if(ret.size() != tmp.size())
                    ret.resize(tmp.size(), 0.);
                if(weighted)
                {
                    ret += tmp*gp.gaussPoints[i].second*weights[i] ;
                    total += gp.gaussPoints[i].second*weights[i] ;
                }
                else
                {
                    ret += tmp*gp.gaussPoints[i].second ;
                    total += gp.gaussPoints[i].second ; 
                }
            }
            if ( total > ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? POINT_TOLERANCE*POINT_TOLERANCE : POINT_TOLERANCE*POINT_TOLERANCE*POINT_TOLERANCE) )
            {
                ret /= total ;
            }
            else
            {
                ret = 0 ;
            }
            if ( cleanup )
            {
                delete vm ;
            }
        }
        #pragma omp atomic write
        lock = false ;
        return total;

    case PRINCIPAL_MECHANICAL_STRAIN_FIELD :

        if ( strainAtGaussPoints.size() != (((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? 3*gp.gaussPoints.size() : 6*gp.gaussPoints.size())) )
        {
            strainAtGaussPoints.resize ( (((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? 3*gp.gaussPoints.size() : 6*gp.gaussPoints.size())) , 0. ) ;
            strainAtGaussPointsSet = false ;
        }

        if ( !strainAtGaussPointsSet)
        {
            strainAtGaussPoints = 0 ;
            strainAtGaussPointsSet = true ;
            Vector tmp ( strainAtGaussPoints.size() /gp.gaussPoints.size() ) ;
            
                
            for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
            {
                getField (MECHANICAL_STRAIN_FIELD, gp.gaussPoints[i].first, tmp, true,vm, i ) ;
                tmp -= getParent()->getBehaviour()->getImposedStrain(gp.gaussPoints[i].first) ;
                for ( size_t j = 0 ; j < strainAtGaussPoints.size() /gp.gaussPoints.size() ; j++ )
                {
                    strainAtGaussPoints[i*strainAtGaussPoints.size() /gp.gaussPoints.size() +j] = tmp[j] ;
                }

                Vector ptmp = toPrincipal(tmp, DOUBLE_OFF_DIAGONAL_VALUES) ;
                if(ret.size() != ptmp.size())
                    ret.resize(ptmp.size(), 0.);
                if(weighted)
                {
                    ret += ptmp *gp.gaussPoints[i].second*weights[i] ;
                    total += gp.gaussPoints[i].second*weights[i] ;
                }
                else
                {
                    ret += ptmp *gp.gaussPoints[i].second ;
                    total += gp.gaussPoints[i].second ;
                }

            }
            if ( total > ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? POINT_TOLERANCE*POINT_TOLERANCE : POINT_TOLERANCE*POINT_TOLERANCE*POINT_TOLERANCE) )
            {
                ret /= total ;
            }
            else
            {
                ret = 0 ;
            }
            if ( cleanup )
            {
                delete vm ;
            }
        }
        else
        { 
            Vector tmp ( strainAtGaussPoints.size() /gp.gaussPoints.size() ) ;
            for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
            {
               
                for ( size_t j = 0 ; j < strainAtGaussPoints.size() /gp.gaussPoints.size() ; j++ )
                {
                    tmp[j] = strainAtGaussPoints[i*strainAtGaussPoints.size() /gp.gaussPoints.size() +j];
                }

                Vector ptmp = toPrincipal(tmp, DOUBLE_OFF_DIAGONAL_VALUES) ;
                if(ret.size() != ptmp.size())
                    ret.resize(ptmp.size(), 0.);
                if(weighted)
                {
                    ret += ptmp*gp.gaussPoints[i].second*weights[i] ;
                    total += gp.gaussPoints[i].second*weights[i] ;
                }
                else
                {
                    ret += ptmp*gp.gaussPoints[i].second ;
                    total += gp.gaussPoints[i].second ;
                }
            }
            if ( total > ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? POINT_TOLERANCE*POINT_TOLERANCE : POINT_TOLERANCE*POINT_TOLERANCE*POINT_TOLERANCE) )
            {
                ret /= total ;
            }
            else
            {
                ret = 0 ;
            }
            if ( cleanup )
            {
                delete vm ;
            }
        }
        #pragma omp atomic write
        lock = false ;
        return total;
    case PRINCIPAL_STRAIN_FIELD :

        if ( pstrainAtGaussPoints.size() != (((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? 2*gp.gaussPoints.size() : 3*gp.gaussPoints.size())) )
        {
            pstrainAtGaussPoints.resize ( (((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? 2*gp.gaussPoints.size() : 3*gp.gaussPoints.size())) , 0. ) ;
            pstrainAtGaussPointsSet = false ;
        }

        if ( !pstrainAtGaussPointsSet)
        {
            pstrainAtGaussPoints = 0 ;
            pstrainAtGaussPointsSet = true ;
            Vector tmp ( pstrainAtGaussPoints.size() /gp.gaussPoints.size() ) ;
            for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
            {
                getField ( f, gp.gaussPoints[i].first, tmp, true, vm, i ) ;
                for ( size_t j = 0 ; j < pstrainAtGaussPoints.size() /gp.gaussPoints.size() ; j++ )
                {
                    pstrainAtGaussPoints[i*pstrainAtGaussPoints.size() /gp.gaussPoints.size() +j] = tmp[j] ;
                }
                if(ret.size() != tmp.size())
                    ret.resize( tmp.size(), 0.);
                if(weighted)
                {
                    ret += tmp*gp.gaussPoints[i].second*weights[i] ;
                    total += gp.gaussPoints[i].second*weights[i] ;
                }
                else
                {
                    ret += tmp*gp.gaussPoints[i].second ;
                    total += gp.gaussPoints[i].second ;
                }

            }
            if ( total > ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? POINT_TOLERANCE*POINT_TOLERANCE : POINT_TOLERANCE*POINT_TOLERANCE*POINT_TOLERANCE) )
            {
                ret /= total ;
            }
            else
            {
                ret = 0 ;
            }
            if ( cleanup )
            {
                delete vm ;
            }
        }
        else
        {
            Vector tmp ( pstrainAtGaussPoints.size() /gp.gaussPoints.size() ) ;
            for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
            {
                for ( size_t j = 0 ; j < pstrainAtGaussPoints.size() /gp.gaussPoints.size() ; j++ )
                {
                    tmp[j] = pstrainAtGaussPoints[i*pstrainAtGaussPoints.size() /gp.gaussPoints.size() +j];
                }
                if(ret.size() != tmp.size())
                    ret.resize(tmp.size(), 0.);
                if(weighted)
                {
                    ret += tmp*gp.gaussPoints[i].second*weights[i] ;
                    total += gp.gaussPoints[i].second*weights[i] ;
                }
                else
                {
                    ret += tmp*gp.gaussPoints[i].second ;
                    total += gp.gaussPoints[i].second ;
                }
            }
            if ( total > ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? POINT_TOLERANCE*POINT_TOLERANCE : POINT_TOLERANCE*POINT_TOLERANCE*POINT_TOLERANCE) )
            {
                ret /= total ;
            }
            else
            {
                ret = 0 ;
            }
            if ( cleanup )
            {
                delete vm ;
            }
        }
        #pragma omp atomic write
        lock = false ;
        return total;
    case REAL_STRESS_FIELD:

        if ( stressAtGaussPoints.size() != (((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? 3*gp.gaussPoints.size() : 6*gp.gaussPoints.size())) )
        {
            stressAtGaussPoints.resize ( (((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? 3*gp.gaussPoints.size() : 6*gp.gaussPoints.size())) , 0. ) ;
            stressAtGaussPoints = false ;
        }

        if ( !stressAtGaussPointsSet)
        {
            stressAtGaussPoints = 0 ;
            stressAtGaussPointsSet = true ;
            Vector tmp ( stressAtGaussPoints.size() /gp.gaussPoints.size() ) ;

            for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
            {
                
                getField ( f, gp.gaussPoints[i].first, tmp, true, vm, i ) ;
                for ( size_t j = 0 ; j < tmp.size() ; j++ )
                {
                    stressAtGaussPoints[i*tmp.size() +j] = tmp[j] ;
                }
                if(ret.size() != tmp.size())
                    ret.resize(tmp.size(), 0.);
                if(weighted)
                {
                    ret += tmp*gp.gaussPoints[i].second*weights[i] ;
                    total += gp.gaussPoints[i].second*weights[i] ;
                }
                else
                {
                    ret += tmp*gp.gaussPoints[i].second ;
                    total += gp.gaussPoints[i].second ;
                }
            }
            if ( total > ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? POINT_TOLERANCE*POINT_TOLERANCE : POINT_TOLERANCE*POINT_TOLERANCE*POINT_TOLERANCE) )
            {
                ret /= total ;
            }
            else
            {
                ret = 0 ;
            }
            if ( cleanup )
            {
                delete vm ;
            }
        }
        else
        {
            Vector tmp ( stressAtGaussPoints.size() /gp.gaussPoints.size() ) ;
            for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
            {
                
                for ( size_t j = 0 ; j < tmp.size() ; j++ )
                {
                    tmp[j] = stressAtGaussPoints[i*tmp.size() +j];
                }
                if(ret.size() != tmp.size())
                    ret.resize(tmp.size(), 0.);
                if(weighted)
                {
                    ret += tmp*gp.gaussPoints[i].second*weights[i] ;
                    total += gp.gaussPoints[i].second*weights[i] ;
                }
                else
                {
                    ret += tmp*gp.gaussPoints[i].second ;
                    total += gp.gaussPoints[i].second ;
                }
            }
            if ( total > ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? POINT_TOLERANCE*POINT_TOLERANCE : POINT_TOLERANCE*POINT_TOLERANCE*POINT_TOLERANCE) )
            {
                ret /= total ;
            }
            else
            {
                ret = 0 ;
            }
            if ( cleanup )
            {
                delete vm ;
            }
        }
        #pragma omp atomic write
        lock = false ;
        return total;
    case PRINCIPAL_REAL_STRESS_FIELD :
        if ( pstressAtGaussPoints.size() != (((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? 2*gp.gaussPoints.size() : 3*gp.gaussPoints.size())) )
        {
            pstressAtGaussPoints.resize ( (((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? 2*gp.gaussPoints.size() : 3*gp.gaussPoints.size())) , 0. ) ;
            pstressAtGaussPoints = false ;
        }

        if ( !pstressAtGaussPointsSet)
        {
            pstressAtGaussPoints = 0 ;
            pstressAtGaussPoints = true ;
            Vector tmp ( pstressAtGaussPoints.size() /gp.gaussPoints.size() ) ;

            for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
            {
                
                getField ( f, gp.gaussPoints[i].first, tmp, true, vm,i ) ;
                for ( size_t j = 0 ; j < tmp.size() ; j++ )
                {
                    pstressAtGaussPoints[i*tmp.size() +j] = tmp[j] ;
                }
                if(ret.size() != tmp.size())
                    ret.resize(tmp.size(), 0.);
                if(weighted)
                {
                    ret += tmp*gp.gaussPoints[i].second*weights[i] ;
                    total += gp.gaussPoints[i].second*weights[i] ;
                }
                else
                {
                    ret += tmp*gp.gaussPoints[i].second ;
                    total += gp.gaussPoints[i].second ;
                }
            }
            if ( total > ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? POINT_TOLERANCE*POINT_TOLERANCE : POINT_TOLERANCE*POINT_TOLERANCE*POINT_TOLERANCE) )
            {
                ret /= total ;
            }
            else
            {
                ret = 0 ;
            }
            if ( cleanup )
            {
                delete vm ;
            }
        }
        else
        {
            Vector tmp ( pstressAtGaussPoints.size() /gp.gaussPoints.size() ) ;
            for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
            {
                
                for ( size_t j = 0 ; j < tmp.size() ; j++ )
                {
                    tmp[j] = pstressAtGaussPoints[i*tmp.size() +j];
                }
                if(ret.size() != tmp.size())
                    ret.resize(tmp.size(), 0.);
                if(weighted)
                {
                    ret += tmp*gp.gaussPoints[i].second*weights[i] ;
                    total += gp.gaussPoints[i].second*weights[i] ;
                }
                else
                {
                    ret += tmp*gp.gaussPoints[i].second ;
                    total += gp.gaussPoints[i].second ;
                }
            }
            if ( total > ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? POINT_TOLERANCE*POINT_TOLERANCE : POINT_TOLERANCE*POINT_TOLERANCE*POINT_TOLERANCE) )
            {
                ret /= total ;
            }
            else
            {
                ret = 0 ;
            }
            if ( cleanup )
            {
                delete vm ;
            }
        }
        #pragma omp atomic write
        lock = false ;
        return total;
    case EFFECTIVE_STRESS_FIELD:
        if ( stressAtGaussPoints.size() != (((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? 3*gp.gaussPoints.size() : 6*gp.gaussPoints.size())) )
        {
            stressAtGaussPoints.resize ( (((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? 3*gp.gaussPoints.size() : 6*gp.gaussPoints.size())) , 0. ) ;
            stressAtGaussPoints = false ;
        }

        if ( !stressAtGaussPointsSet)
        {
            stressAtGaussPoints = 0 ;
            stressAtGaussPointsSet = true ;
            Vector tmp ( stressAtGaussPoints.size() /gp.gaussPoints.size() ) ;

            for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
            {
                
                getField ( f, gp.gaussPoints[i].first, tmp, true, vm, i ) ;
                for ( size_t j = 0 ; j < tmp.size() ; j++ )
                {
                    stressAtGaussPoints[i*tmp.size() +j] = tmp[j] ;
                }
                if(ret.size() != tmp.size())
                    ret.resize(tmp.size(), 0.);
                if(weighted)
                {
                    ret += tmp*gp.gaussPoints[i].second*weights[i] ;
                    total += gp.gaussPoints[i].second*weights[i] ;
                }
                else
                {
                    ret += tmp*gp.gaussPoints[i].second ;
                    total += gp.gaussPoints[i].second ;
                }
            }
            if ( total > ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? POINT_TOLERANCE*POINT_TOLERANCE : POINT_TOLERANCE*POINT_TOLERANCE*POINT_TOLERANCE) )
            {
                ret /= total ;
            }
            else
            {
                ret = 0 ;
            }
            if ( cleanup )
            {
                delete vm ;
            }
        }
        else
        {
            Vector tmp ( stressAtGaussPoints.size() /gp.gaussPoints.size() ) ;
            for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
            {
                
                for ( size_t j = 0 ; j < tmp.size() ; j++ )
                {
                    tmp[j] = stressAtGaussPoints[i*tmp.size() +j];
                }
                if(ret.size() != tmp.size())
                    ret.resize(tmp.size(), 0.);
                if(weighted)
                {
                    ret += tmp*gp.gaussPoints[i].second*weights[i] ;
                    total += gp.gaussPoints[i].second*weights[i] ;
                }
                else
                {
                    ret += tmp*gp.gaussPoints[i].second ;
                    total += gp.gaussPoints[i].second ;
                }
            }
            if ( total > ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? POINT_TOLERANCE*POINT_TOLERANCE : POINT_TOLERANCE*POINT_TOLERANCE*POINT_TOLERANCE) )
            {
                ret /= total ;
            }
            else
            {
                ret = 0 ;
            }
            if ( cleanup )
            {
                delete vm ;
            }
        }
        #pragma omp atomic write
        lock = false ;
        return total;
    case PRINCIPAL_EFFECTIVE_STRESS_FIELD :
        if ( pstressAtGaussPoints.size() != (((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? 2*gp.gaussPoints.size() : 3*gp.gaussPoints.size())) )
        {
            pstressAtGaussPoints.resize ( (((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? 2*gp.gaussPoints.size() : 3*gp.gaussPoints.size())) , 0. ) ;
            pstressAtGaussPoints = false ;
        }

        if ( !pstressAtGaussPointsSet)
        {
            pstressAtGaussPoints = 0 ;
            pstressAtGaussPointsSet = true ;
            Vector tmp ( pstressAtGaussPoints.size() /gp.gaussPoints.size() ) ;

            for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
            {
                
                getField ( f, gp.gaussPoints[i].first, tmp, true, vm, i ) ;
                for ( size_t j = 0 ; j < tmp.size() ; j++ )
                {
                    pstressAtGaussPoints[i*tmp.size() +j] = tmp[j] ;
                }
                if(ret.size() != tmp.size())
                    ret.resize(tmp.size(), 0.);
                if(weighted)
                {
                    ret += tmp*gp.gaussPoints[i].second*weights[i] ;
                    total += gp.gaussPoints[i].second*weights[i] ;
                }
                else
                {
                    ret += tmp*gp.gaussPoints[i].second ;
                    total += gp.gaussPoints[i].second ;
                }
            }
            if ( total > ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? POINT_TOLERANCE*POINT_TOLERANCE : POINT_TOLERANCE*POINT_TOLERANCE*POINT_TOLERANCE) )
            {
                ret /= total ;
            }
            else
            {
                ret = 0 ;
            }
            if ( cleanup )
            {
                delete vm ;
            }
        }
        else
        {
            Vector tmp ( pstressAtGaussPoints.size() /gp.gaussPoints.size() ) ;
            for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
            {
                
                for ( size_t j = 0 ; j < tmp.size() ; j++ )
                {
                    tmp[j] = pstressAtGaussPoints[i*tmp.size() +j];
                }
                if(ret.size() != tmp.size())
                    ret.resize(tmp.size(), 0.);
                if(weighted)
                {
                    ret += tmp*gp.gaussPoints[i].second*weights[i] ;
                    total += gp.gaussPoints[i].second*weights[i] ;
                }
                else
                {
                    ret += tmp*gp.gaussPoints[i].second ;
                    total += gp.gaussPoints[i].second ;
                }
            }
            if ( total > ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? POINT_TOLERANCE*POINT_TOLERANCE : POINT_TOLERANCE*POINT_TOLERANCE*POINT_TOLERANCE) )
            {
                ret /= total ;
            }
            else
            {
                ret = 0 ;
            }
            if ( cleanup )
            {
                delete vm ;
            }
        }
        #pragma omp atomic write
        lock = false ;
        return total;

    case PRINCIPAL_STRESS_ANGLE_FIELD:
        if ( stressAtGaussPoints.size() != (((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? 3*gp.gaussPoints.size() : 6*gp.gaussPoints.size())) )
        {
            stressAtGaussPoints.resize ( (((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? 3*gp.gaussPoints.size() : 6*gp.gaussPoints.size())) , 0. ) ;
            stressAtGaussPoints = false ;
        }

        if ( !stressAtGaussPointsSet)
        {
            stressAtGaussPointsSet = true ;
            Vector tmp ( stressAtGaussPoints.size() /gp.gaussPoints.size() ) ;
            for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
            {
                
                getField ( REAL_STRESS_FIELD, gp.gaussPoints[i].first, tmp, true, vm, i ) ;
                for ( size_t j = 0 ; j < tmp.size() ; j++ )
                {
                    stressAtGaussPoints[i*tmp.size() +j] = tmp[j] ;
                }
                if(ret.size() != 1)
                    ret.resize(1, 0.);

                double w = 1 ;
                if(weighted)
                    w = weights[i] ;
                if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
                {
                    ret[0] +=  atan2 ( tmp[2], tmp[0] - tmp[1] )*gp.gaussPoints[i].second*w ;
                }
                else
                {
                    ret[0] += atan2 ( tmp[3], tmp[0] - tmp[1] )*gp.gaussPoints[i].second*w ;
                    ret[1] += atan2 ( tmp[4], tmp[0] - tmp[2] )*gp.gaussPoints[i].second*w ;
                    ret[2] += atan2 ( tmp[5], tmp[1] - tmp[2] )*gp.gaussPoints[i].second*w ;
                }
                total += gp.gaussPoints[i].second*w ;
            }
            if ( total > ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? POINT_TOLERANCE*POINT_TOLERANCE : POINT_TOLERANCE*POINT_TOLERANCE*POINT_TOLERANCE) )
            {
                ret /= total ;
            }
            else
            {
                ret = 0 ;
            }
            if ( cleanup )
            {
                delete vm ;
            }
        }
        else
        {
            Vector tmp ( stressAtGaussPoints.size() /gp.gaussPoints.size() ) ;
            for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
            {
                
                for ( size_t j = 0 ; j < tmp.size() ; j++ )
                {
                    tmp[j] = stressAtGaussPoints[i*tmp.size() +j];
                }
                if(ret.size() != 1)
                    ret.resize(1, 0.);

                double w = 1 ;
                if(weighted)
                    w = weights[i] ;
                if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
                {
                    ret[0] +=  atan2 ( tmp[2], tmp[0] - tmp[1]  )*gp.gaussPoints[i].second*w ;
                }
                else
                {
                    ret[0] += atan2 ( tmp[3], tmp[0] - tmp[1] )*gp.gaussPoints[i].second*w ;
                    ret[1] += atan2 ( tmp[4], tmp[0] - tmp[2] )*gp.gaussPoints[i].second*w ;
                    ret[2] += atan2 ( tmp[5], tmp[1] - tmp[2] )*gp.gaussPoints[i].second*w ;
                }
                total += gp.gaussPoints[i].second*w ;
            }
            if ( total > ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? POINT_TOLERANCE*POINT_TOLERANCE : POINT_TOLERANCE*POINT_TOLERANCE*POINT_TOLERANCE) )
            {
                ret /= total ;
            }
            else
            {
                ret = 0 ;
            }
            if ( cleanup )
            {
                delete vm ;
            }
        }
        #pragma omp atomic write
        lock = false ;
        return total;

    case PRINCIPAL_STRAIN_ANGLE_FIELD:
        if ( strainAtGaussPoints.size() != (((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? 3*gp.gaussPoints.size() : 6*gp.gaussPoints.size())) )
        {
            strainAtGaussPoints.resize ( (((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? 3*gp.gaussPoints.size() : 6*gp.gaussPoints.size())) , 0. ) ;
            strainAtGaussPoints = false ;
        }

        if ( !strainAtGaussPointsSet)
        {
            strainAtGaussPointsSet = true ;
            Vector tmp ( strainAtGaussPoints.size() /gp.gaussPoints.size() ) ;

            for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
            {
                
                getField ( MECHANICAL_STRAIN_FIELD, gp.gaussPoints[i].first, tmp, true, vm, i ) ;
                for ( size_t j = 0 ; j < tmp.size() ; j++ )
                {
                    strainAtGaussPoints[i*tmp.size() +j] = tmp[j] ;
                }
                if(ret.size() != 1)
                    ret.resize(1, 0.);
                
                double w = 1 ;
                if(weighted)
                    w = weights[i] ;

                if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
                {
                    ret[0] +=  atan2 ( tmp[2], tmp[0] - tmp[1] )*gp.gaussPoints[i].second*w ;
                }
                else
                {
                    ret[0] += atan2 ( tmp[3], tmp[0] - tmp[1] )*gp.gaussPoints[i].second*w ;
                    ret[1] += atan2 ( tmp[4], tmp[0] - tmp[2] )*gp.gaussPoints[i].second*w ;
                    ret[2] += atan2 ( tmp[5], tmp[1] - tmp[2] )*gp.gaussPoints[i].second*w ;
                }
                total += gp.gaussPoints[i].second*w ;
            }
            if ( total > ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? POINT_TOLERANCE*POINT_TOLERANCE : POINT_TOLERANCE*POINT_TOLERANCE*POINT_TOLERANCE) )
            {
                ret /= total ;
            }
            else
            {
                ret = 0 ;
            }
            if ( cleanup )
            {
                delete vm ;
            }
        }
        else
        {
            Vector tmp ( strainAtGaussPoints.size() /gp.gaussPoints.size() ) ;
            for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
            {
                
                for ( size_t j = 0 ; j < tmp.size() ; j++ )
                {
                    tmp[j] = strainAtGaussPoints[i*tmp.size() +j];
                }
                if(ret.size() != 1)
                    ret.resize(1, 0.);

                double w = 1 ;
                if(weighted)
                    w = weights[i] ;
                
                if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
                {
                    if(ret.size() != 1)
                        ret.resize(1, 0.);
                    ret[0] +=  atan2 ( tmp[2], tmp[0] - tmp[1] )*gp.gaussPoints[i].second*w ;
                }
                else
                {
                    if(ret.size() != 3)
                        ret.resize(3, 0.);
                    ret[0] += atan2 ( tmp[3], tmp[0] - tmp[1] )*gp.gaussPoints[i].second*w ;
                    ret[1] += atan2 ( tmp[4], tmp[0] - tmp[2] )*gp.gaussPoints[i].second*w ;
                    ret[2] += atan2 ( tmp[5], tmp[1] - tmp[2] )*gp.gaussPoints[i].second*w ;
                }
                total += gp.gaussPoints[i].second*w ;
            }
            if ( total > ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? POINT_TOLERANCE*POINT_TOLERANCE : POINT_TOLERANCE*POINT_TOLERANCE*POINT_TOLERANCE) )
            {
                ret /= total ;
            }
            else
            {
                ret = 0 ;
            }
            if ( cleanup )
            {
                delete vm ;
            }
        }
        #pragma omp atomic write
        lock = false ;
        return total;

    default :

        Vector tmp (ret) ;
        for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
        {

            getField ( f, gp.gaussPoints[i].first, tmp, true,  vm, index ) ;

            double w = 1 ;
            if(weighted)
                w = weights[i] ;
            ret += tmp*gp.gaussPoints[i].second*w ;
            total += gp.gaussPoints[i].second*w ;
        }
        if ( total > ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? POINT_TOLERANCE*POINT_TOLERANCE : POINT_TOLERANCE*POINT_TOLERANCE*POINT_TOLERANCE) )
        {
            ret /= total ;
        }
        else
        {
            ret = 0 ;
        }
        if ( cleanup )
        {
            delete vm ;
        }
        #pragma omp atomic write
        lock = false ;
        return total;
    }
}

double ElementState::getAverageField ( FieldType f, FieldType f_, Vector & ret, Vector & ret_, VirtualMachine * vm,  double time, const std::vector<double> & weights, int index )
{
    while(true) {
//         usleep(1) ;
        bool test = true;
        #pragma omp flush(lock)
        #pragma omp atomic read
        test = lock ;
        if(!test)
        {
            #pragma omp atomic write
            lock = true ;
            #pragma omp flush(lock)
            break ;
        }
    }

    bool cleanup = !vm ;
    if ( !vm )
    {
        vm = new VirtualMachine() ;
    }
    const GaussPointArray & gp = parent->getGaussPoints() ;
    ret = 0 ;
    ret_ = 0 ;
    double v = 0 ; //(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? parent->area() : parent->volume() ;
    double total = 0 ;
    bool weighted = true ;
    if ( weights.size() != gp.gaussPoints.size() )
    {
        weighted = false ;
    }

    if ( f == MECHANICAL_STRAIN_FIELD && ( f_ == EFFECTIVE_STRESS_FIELD || f_ == REAL_STRESS_FIELD ) )
    {
        if ( strainAtGaussPoints.size() != ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)?3*gp.gaussPoints.size() :6*gp.gaussPoints.size()) ||
                stressAtGaussPoints.size() != ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)?3*gp.gaussPoints.size() :6*gp.gaussPoints.size()) )
        {
            strainAtGaussPoints.resize ( ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)?3*gp.gaussPoints.size() :6*gp.gaussPoints.size()) ) ;
            stressAtGaussPoints.resize ( ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)?3*gp.gaussPoints.size() :6*gp.gaussPoints.size()) ) ;
            strainAtGaussPointsSet = false ;
            stressAtGaussPointsSet = false ;
        }


        if ( !strainAtGaussPointsSet || !stressAtGaussPointsSet)
        {
            strainAtGaussPointsSet = true ;
            stressAtGaussPointsSet = true ;
            Vector tmp ( strainAtGaussPoints.size() /gp.gaussPoints.size() ) ;
            Vector tmp_ ( stressAtGaussPoints.size() /gp.gaussPoints.size() ) ;
            for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
            {
                
                getField ( f,f_, gp.gaussPoints[i].first, tmp, tmp_, true,vm, i ) ;

                for ( size_t j = 0 ; j < strainAtGaussPoints.size() /gp.gaussPoints.size() ; j++ )
                {
                    strainAtGaussPoints[i*strainAtGaussPoints.size() /gp.gaussPoints.size() +j] = tmp[j] ;
                    stressAtGaussPoints[i*stressAtGaussPoints.size() /gp.gaussPoints.size() +j] = tmp_[j] ;
                }

                if(weighted)
                {
                    ret += tmp*gp.gaussPoints[i].second*weights[i] ;
                    ret_ += tmp_*gp.gaussPoints[i].second*weights[i] ;
                    total += gp.gaussPoints[i].second*weights[i] ;
                }
                else
                {
                    ret += tmp*gp.gaussPoints[i].second ;
                    ret_ += tmp_*gp.gaussPoints[i].second ;
                    total += gp.gaussPoints[i].second ;
                }
            }
            if ( total > ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? POINT_TOLERANCE*POINT_TOLERANCE : POINT_TOLERANCE*POINT_TOLERANCE*POINT_TOLERANCE) )
            {
                ret_ /= total ;
                ret /= total ;
            }
            else
            {
                ret = 0 ;
                ret_ = 0 ;
            }
            if ( cleanup )
            {
                delete vm ;
            }
        }
        else
        {
           Vector tmp ( strainAtGaussPoints.size() /gp.gaussPoints.size() ) ;
           Vector tmp_ ( stressAtGaussPoints.size() /gp.gaussPoints.size() ) ;
            for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
            {
                for ( size_t j = 0 ; j < strainAtGaussPoints.size() /gp.gaussPoints.size() ; j++ )
                {
                    tmp[j] = strainAtGaussPoints[i*strainAtGaussPoints.size() /gp.gaussPoints.size() +j];
                    tmp_[j] = stressAtGaussPoints[i*stressAtGaussPoints.size() /gp.gaussPoints.size() +j];
                }
                if(weighted)
                {
                    ret += tmp*gp.gaussPoints[i].second*weights[i] ;
                    ret_ += tmp_*gp.gaussPoints[i].second*weights[i] ;
                    total += gp.gaussPoints[i].second*weights[i] ;
                }
                else
                {
                    ret += tmp*gp.gaussPoints[i].second ;
                    ret_ += tmp_*gp.gaussPoints[i].second ;
                    total += gp.gaussPoints[i].second ;
                }
            }
            if ( total > ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? POINT_TOLERANCE*POINT_TOLERANCE : POINT_TOLERANCE*POINT_TOLERANCE*POINT_TOLERANCE) )
            {
                ret_ /= total ;
                ret /= total ;
            }
            else
            {
                ret = 0 ;
                ret_ = 0 ;
            }
            if ( cleanup )
            {
                delete vm ;
            }
        }
        #pragma omp atomic write
        lock = false ;
        return total ;
    }
    
    if ( f == PRINCIPAL_STRAIN_FIELD && ( f_ == PRINCIPAL_EFFECTIVE_STRESS_FIELD || f == PRINCIPAL_REAL_STRESS_FIELD ) )
    {
        if ( pstrainAtGaussPoints.size() != ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)?2*gp.gaussPoints.size() :3*gp.gaussPoints.size()) ||
                pstressAtGaussPoints.size() != ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)?2*gp.gaussPoints.size() :3*gp.gaussPoints.size()) )
        {
            pstrainAtGaussPoints.resize ( ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)?2*gp.gaussPoints.size() :3*gp.gaussPoints.size()) ) ;
            pstressAtGaussPoints.resize ( ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)?2*gp.gaussPoints.size() :3*gp.gaussPoints.size()) ) ;
            pstrainAtGaussPointsSet = false ;
            pstressAtGaussPointsSet = false ;
        }


        if ( !pstrainAtGaussPointsSet || !pstressAtGaussPointsSet)
        {

            pstrainAtGaussPointsSet = true ;
            pstressAtGaussPointsSet = true ;
            Vector tmp ( pstrainAtGaussPoints.size() /gp.gaussPoints.size() ) ;
            Vector tmp_ ( pstressAtGaussPoints.size() /gp.gaussPoints.size() ) ;
            for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
            {
                
                getField ( f,f_, gp.gaussPoints[i].first, tmp,tmp_, true, vm, i ) ;
                for ( size_t j = 0 ; j < pstrainAtGaussPoints.size() /gp.gaussPoints.size() ; j++ )
                {
                    pstrainAtGaussPoints[i*pstrainAtGaussPoints.size() /gp.gaussPoints.size() +j] = tmp[j] ;
                    pstressAtGaussPoints[i*pstressAtGaussPoints.size() /gp.gaussPoints.size() +j] = tmp_[j] ;
                }
                if(weighted)
                {
                    ret += tmp*gp.gaussPoints[i].second*weights[i] ;
                    ret_ += tmp_*gp.gaussPoints[i].second*weights[i] ;
                    total += gp.gaussPoints[i].second*weights[i] ;
                }
                else
                {
                    ret += tmp*gp.gaussPoints[i].second ;
                    ret_ += tmp_*gp.gaussPoints[i].second ;
                    total += gp.gaussPoints[i].second ;
                }

            }
            if ( total > ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? POINT_TOLERANCE*POINT_TOLERANCE : POINT_TOLERANCE*POINT_TOLERANCE*POINT_TOLERANCE) )
            {
                ret_ /= total ;
                ret /= total ;
            }
            else
            {
                ret = 0 ;
                ret_ = 0 ;
            }
            if ( cleanup )
            {
                delete vm ;
            }
        }
        else
        {
            Vector tmp ( pstrainAtGaussPoints.size() /gp.gaussPoints.size() ) ;
            Vector tmp_ ( pstressAtGaussPoints.size() /gp.gaussPoints.size() ) ;
            for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
            { 
                for ( size_t j = 0 ; j < pstrainAtGaussPoints.size() /gp.gaussPoints.size() ; j++ )
                {
                    tmp[j] = pstrainAtGaussPoints[i*pstrainAtGaussPoints.size() /gp.gaussPoints.size() +j];
                    tmp_[j] = pstressAtGaussPoints[i*pstressAtGaussPoints.size() /gp.gaussPoints.size() +j];
                }
                if(weighted)
                {
                    ret += tmp*gp.gaussPoints[i].second*weights[i] ;
                    ret_ += tmp_*gp.gaussPoints[i].second*weights[i] ;
                    total += gp.gaussPoints[i].second*weights[i] ;
                }
                else
                {
                    ret += tmp*gp.gaussPoints[i].second ;
                    ret_ += tmp_*gp.gaussPoints[i].second ;
                    total += gp.gaussPoints[i].second ;
                }
            }
            if ( total > ((parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ? POINT_TOLERANCE*POINT_TOLERANCE : POINT_TOLERANCE*POINT_TOLERANCE*POINT_TOLERANCE) )
            {
                ret_ /= total ;
                ret /= total ;
            }
            else
            {
                ret = 0 ;
                ret_ = 0 ;
            }

        }
        #pragma omp atomic write
        lock = false ;
        return total;
    }


    getAverageField ( f, ret,vm, time, weights, index ) ;
    v = getAverageField ( f_, ret_,vm, time, weights, index );
    if ( cleanup )
    {
        delete vm ;
    }
    #pragma omp atomic write
    lock = false ;
    return v ;

}


void ElementState::getField ( FieldType f1, FieldType f2, const Point & p, Vector & ret1, Vector & ret2, bool local, VirtualMachine * vm, int , int ) 
{
    bool cleanup = !vm ;
    if ( !vm )
    {
        vm = new VirtualMachine() ;
    }
    Point p_ = p ;
    if ( !local )
    {
        p_ = parent->inLocalCoordinates ( p ) ;
    }
    if ( isStrainField ( f1 ) && isStressField ( f2 ) )
    {
        getField ( f1, p, ret1, local, vm ) ;
        if ( parent->getBehaviour()->getTensor ( p_, parent ).numCols() != ret2.size() )
        {
            ret2 = 0 ;
            return ;
        }
        if ( isRealStressField ( f2 ) )
        {
            ret2 = ( Vector ) ( parent->getBehaviour()->getTensor ( p_, parent ) * ret1 ) - getParent()->getBehaviour()->getImposedStress ( p_, parent ) ;
        }
        else
        {
            ret2 = ( Vector ) ( parent->getBehaviour()->param * ret1 ) - getParent()->getBehaviour()->getImposedStress ( p_, parent ) ;
        }
        return ;
    }
    if ( f1 == PRINCIPAL_MECHANICAL_STRAIN_FIELD && ( f2 == PRINCIPAL_EFFECTIVE_STRESS_FIELD || f2 == PRINCIPAL_REAL_STRESS_FIELD ) )
    {
        Vector v1 ( 0., 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ;
        Vector v2 ( 0., v1.size() ) ;
        if ( isRealStressField ( f2 ) )
        {
            getField ( MECHANICAL_STRAIN_FIELD, REAL_STRESS_FIELD, p, v1, v2, local, vm ) ;
        }
        else
        {
            getField ( MECHANICAL_STRAIN_FIELD, EFFECTIVE_STRESS_FIELD, p, v1, v2, local, vm ) ;
        }
        ret1 = toPrincipal ( v1 , DOUBLE_OFF_DIAGONAL_VALUES) ;
        ret2 = toPrincipal ( v2 , SINGLE_OFF_DIAGONAL_VALUES ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    if ( isStrainField ( f2 ) && isStressField ( f1 ) )
    {
        getField ( f2, p, ret2, local, vm ) ;
        if ( parent->getBehaviour()->getTensor ( p_, parent ).numCols() != ret1.size() )
        {
            ret1 = 0 ;
            if ( cleanup )
            {
                delete vm ;
            }
            return ;
        }
        if ( isRealStressField ( f1 ) )
        {
            ret1 = ( Vector ) ( parent->getBehaviour()->getTensor ( p_, parent ) * ret2 ) - getParent()->getBehaviour()->getImposedStress ( p_, parent ) ;
        }
        else
        {
            ret1 = ( Vector ) ( parent->getBehaviour()->param * ret2 ) - getParent()->getBehaviour()->getImposedStress ( p_, parent ) ;
        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    if ( f2 == PRINCIPAL_MECHANICAL_STRAIN_FIELD && ( f1 == PRINCIPAL_EFFECTIVE_STRESS_FIELD || f1 == PRINCIPAL_REAL_STRESS_FIELD ) )
    {
        Vector v1 ( 0., 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ;
        Vector v2 ( 0., v1.size() ) ;
        getField ( isRealStressField ( f2 ) ? REAL_STRESS_FIELD : EFFECTIVE_STRESS_FIELD, MECHANICAL_STRAIN_FIELD, p, v1, v2, local,vm ) ;

        ret1 = toPrincipal ( v1, SINGLE_OFF_DIAGONAL_VALUES ) ;
        ret2 = toPrincipal ( v2, DOUBLE_OFF_DIAGONAL_VALUES ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    if ( f1 == GRADIENT_FIELD && f2 == FLUX_FIELD )
    {
        getField ( f1, p, ret1, local, vm ) ;
        ret2 = ( Vector ) ( parent->getBehaviour()->getTensor ( p_, parent ) * ret1 ) ;
    }
    if ( f1 == FLUX_FIELD && f2 == GRADIENT_FIELD )
    {
        getField ( f2, p, ret2, local, vm ) ;
        ret1 = ( Vector ) ( parent->getBehaviour()->getTensor ( p_, parent ) * ret2 ) ;
    }
    if ( cleanup )
    {
        delete vm ;
    }
}

void ElementState::getField ( FieldType f1, FieldType f2, const PointArray & p, Vector & ret1, Vector & ret2, bool local, VirtualMachine * vm, int , int ) 
{
    bool cleanup = !vm ;
    if ( !vm )
    {
        vm = new VirtualMachine() ;
    }
    Vector b1 ( 0., ret1.size() /p.size() ) ;
    Vector b2 ( 0., ret2.size() /p.size() ) ;
    for ( size_t i = 0 ; i < p.size() ; i++ )
    {
        getField ( f1, f2, *p[i], b1, b2, local, vm ) ;
        for ( size_t j = b1.size() *i ; j < b1.size() * ( i+1 ) ; j++ )
        {
            ret1[j] = b1[j - b1.size() *i] ;
        }
        for ( size_t j = b2.size() *i ; j < b2.size() * ( i+1 ) ; j++ )
        {
            ret2[j] = b2[j - b2.size() *i] ;
        }
    }
    if ( cleanup )
    {
        delete vm ;
    }
}

void ElementState::getField ( FieldType f1, FieldType f2, const std::valarray<std::pair<Point, double> > & p, Vector & ret1, Vector & ret2, bool local, VirtualMachine * vm, int , int ) 
{
    bool cleanup = !vm ;
    if ( !vm )
    {
        vm = new VirtualMachine() ;
    }

    Vector b1 ( 0., ret1.size() /p.size() ) ;
    Vector b2 ( 0., ret2.size() /p.size() ) ;
    for ( size_t i = 0 ; i < p.size() ; i++ )
    {
        getField ( f1, f2, p[i].first, b1, b2, local, vm ) ;
        for ( size_t j = b1.size() *i ; j < b1.size() * ( i+1 ) ; j++ )
        {
            ret1[j] = b1[j - b1.size() *i] ;
        }
        for ( size_t j = b2.size() *i ; j < b2.size() * ( i+1 ) ; j++ )
        {
            ret2[j] = b2[j - b2.size() *i] ;
        }
    }
    if ( cleanup )
    {
        delete vm ;
    }
}

void ElementState::getFieldAtGaussPoint ( FieldType f1, FieldType f2, size_t p, Vector & ret1, Vector & ret2, VirtualMachine * vm, int i, int j )
{
    bool cleanup = !vm ;
    if ( !vm )
    {
        vm = new VirtualMachine() ;
    }
    Point p_ = parent->getGaussPoints().gaussPoints[p].first ;
    getField ( f1, f2, p_, ret1, ret2, true, vm, i, j ) ;
    if ( cleanup )
    {
        delete vm ;
    }
}

std::vector<double> ElementState::getEnrichedInterpolatingFactors ( const Point &p, bool local ) const
{

    std::vector<double> ret ( parent->getEnrichmentFunctions().size() ) ;
    VirtualMachine vm ;
    Point p_ = p ;

    if ( !local )
    {
        p_ = parent->inLocalCoordinates ( p ) ;
    }

    for ( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
    {
        ret[j] = vm.eval ( parent->getEnrichmentFunction ( j ), p_ ) ;
    }

    return ret;
}


std::vector<double> ElementState::getInterpolatingFactors ( const Point &p, bool local ) const
{

    std::vector<double> ret ( parent->getBoundingPoints().size() + parent->getEnrichmentFunctions().size() ) ;
    VirtualMachine vm ;
    Point p_ = p ;

    if ( !local )
    {
        p_ = parent->inLocalCoordinates ( p ) ;
    }

    for ( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
    {
        ret[j] = vm.eval ( parent->getShapeFunction ( j ), p_ ) ;
    }

    for ( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
    {
        ret[j + parent->getBoundingPoints().size()] = vm.eval ( parent->getEnrichmentFunction ( j ), p_ ) ;
    }

    return ret;
}

void ElementState::initialize ( Mesh<DelaunayTetrahedron,DelaunayTreeItem3D> * msh)
{
    mesh3d = msh ;
    size_t ndofs = 0 ;
    if ( parent->getBehaviour() )
    {
        ndofs = parent->getBehaviour()->getNumberOfDegreesOfFreedom() ;
    }
    displacements.resize ( parent->getBoundingPoints().size() *ndofs, 0. ) ;

    if ( std::abs ( timePos - previousTimePos ) < POINT_TOLERANCE && std::abs ( timePos ) < POINT_TOLERANCE )
    {
        timePos = -0.1 ;
        previousTimePos = -0.2 ;
    }

    strainAtGaussPointsSet = false ;
    stressAtGaussPointsSet = false ;
    pstrainAtGaussPointsSet = false ;
    pstressAtGaussPointsSet = false ;
}

void ElementState::initialize ( Mesh<DelaunayTriangle,DelaunayTreeItem> * msh)
{
    mesh2d = msh ;
    size_t ndofs = 0 ;
    if ( parent->getBehaviour() )
    {
        ndofs = parent->getBehaviour()->getNumberOfDegreesOfFreedom() ;
    }
    displacements.resize ( parent->getBoundingPoints().size() *ndofs, 0. ) ;
    displacements = 0 ;

    if ( std::abs ( timePos - previousTimePos ) < POINT_TOLERANCE && std::abs ( timePos ) < POINT_TOLERANCE )
    {
        timePos = -0.1 ;
        previousTimePos = -0.2 ;
    }

    strainAtGaussPointsSet = false ;
    stressAtGaussPointsSet = false ;
    pstrainAtGaussPointsSet = false ;
    pstressAtGaussPointsSet = false ;
}

void ElementState::step ( double dt, const Vector *d )
{

    strainAtGaussPointsSet = false ;
    stressAtGaussPointsSet = false ;
    pstrainAtGaussPointsSet = false ;
    pstressAtGaussPointsSet = false ;

    if ( parent->getBehaviour() && parent->getBehaviour()->type != VOID_BEHAVIOUR )
    {
        previousTimePos = timePos ;
        timePos += dt ;
        size_t ndofs = parent->getBehaviour()->getNumberOfDegreesOfFreedom() ;

        std::vector< size_t > ids = parent->getDofIds() ;

        if ( ids.empty() )
        {
            return ;
        }

        if ( displacements.size() != parent->getBoundingPoints().size() *ndofs)
        {
            displacements.resize ( parent->getBoundingPoints().size() *ndofs) ;
        }
        
        if ( enrichedDisplacements.size() != parent->getEnrichmentFunctions().size() *ndofs )
        {
            enrichedDisplacements.resize ( parent->getEnrichmentFunctions().size() *ndofs ) ;
        }

        for ( size_t i = 0 ; i < parent->getShapeFunctions().size() ; i++ )
        {
            for ( size_t j = 0 ; j < ndofs ; j++ )
            {
                if ( ids[i] * ndofs + j < d->size() )
                {
                    displacements[i * ndofs + j] = ( *d ) [ids[i] * ndofs + j] ;
                }
                else
                {
                    displacements[i * ndofs + j] = 0 ;
                }
            }
        }

        int nbp = parent->getBoundingPoints().size() ;

        for ( size_t i = 0 ; i < parent->getEnrichmentFunctions().size() ; i++ )
        {
            for ( size_t j = 0 ; j < ndofs ; j++ )
            {
                if ( ids[i + nbp] * ndofs + j < d->size() )
                {
                    enrichedDisplacements[i * ndofs  + j] = ( *d ) [ids[i + nbp] * ndofs + j] ;
                }
                else
                {
                    enrichedDisplacements[i * ndofs + j] = 0 ;
                }
            }
        }
    }
}


double ElementState::getTime() const
{
    return timePos ;
}

double ElementState::getDeltaTime() const
{
    return timePos - previousTimePos ;
}

double ElementState::getNodalCentralTime() const
{
    if( parent->timePlanes() < 2)
        return timePos ;
    return 0.5* ( parent->getBoundingPoint ( parent->getBoundingPoints().size() -1 ).getT() + parent->getBoundingPoint ( 0 ).getT() ) ;
}

double ElementState::getNodalDeltaTime() const
{
    return parent->getBoundingPoint ( parent->getBoundingPoints().size() -1 ).getT() - parent->getBoundingPoint ( 0 ).getT() ;
}

FieldType stressFieldType ( StressCalculationMethod m )
{
    return m == REAL_STRESS ? REAL_STRESS_FIELD : EFFECTIVE_STRESS_FIELD ;
}

FieldType principalStressFieldType ( StressCalculationMethod m )
{
    return m == REAL_STRESS ? PRINCIPAL_REAL_STRESS_FIELD : PRINCIPAL_EFFECTIVE_STRESS_FIELD ;
}

std::vector<Point> ElementState::getPrincipalFrame ( const Point &p, bool local, VirtualMachine * vm, StressCalculationMethod m )
{
    bool cleanup = !vm ;
    if ( !vm )
    {
        vm = new VirtualMachine() ;
    }

    if ( getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
    {
        Vector principalStresses ( 0., 2 ) ;
        getAverageField ( principalStressFieldType ( m ), principalStresses, vm, 0. ) ;
        double principalAngle = 0.5*atan2 ( 2.*principalStresses[2], principalStresses[0]-principalStresses[1]  ) ;
        std::vector<Point> ret ;
        ret.push_back ( Point ( cos ( principalAngle ), -sin ( principalAngle ) ) );
        ret.push_back ( Point ( sin ( principalAngle ), cos ( principalAngle ) ) );
        ret.push_back ( Point ( 0,0, 1 ) );
        if ( cleanup )
        {
            delete vm ;
        }
        return ret ;
    }

    Vector stress ( 0., 6 ) ;
    Point p_ = p ;
    if ( !local )
    {
        p_ = parent->inLocalCoordinates ( p_ ) ;
    }
    this->getField ( principalStressFieldType ( m ), p_, stress, true ) ;
    Matrix stressMatrix = makeStressOrStrainMatrix ( stress ) ;

    //then, we get the principal stresses
    Vector principalStresses ( 0., 3 ) ;
    this->getField ( principalStressFieldType ( m ), p_, principalStresses, true ) ;
    std::vector<Vector> principalVectors ;
    for ( size_t i = 0 ; i <  principalStresses.size() ; ++i )
    {
        principalVectors.push_back ( Vector ( 0.,stressMatrix.numCols() ) ) ;
        Matrix m = ( stressMatrix-identity ( stressMatrix.numCols() ) *principalStresses[i] ) ;
        //Assume that the first coefficient is 1 and get a reduced system.
        Matrix m_reduced ( m ) ;
        Vector v_reduced ( 0., stressMatrix.numCols() ) ;
        v_reduced[0] = 1 ;
        m_reduced[0][0] = 1 ;
        for ( size_t j = 1 ; j < stressMatrix.numCols() ; j++ )
        {
            m_reduced[0][j] = 0 ;
        }
        for ( size_t j = 1 ; j < stressMatrix.numRows() ; j++ )
        {
            v_reduced[j] -= m_reduced[j][0] ;
            m_reduced[j][0] = 0 ;
        }

        size_t k = 1 ;
        while ( det ( m_reduced ) < 1e-6 && k < m_reduced.numCols() )
        {
            m_reduced = m ;
            v_reduced.resize ( stressMatrix.numCols(),0. ) ;
            v_reduced[k] = 1 ;
            m_reduced[k][k] = 1 ;
            for ( size_t j = 0 ; j < stressMatrix.numCols() ; j++ )
            {
                if ( j == k )
                {
                    continue ;
                }
                m_reduced[k][j] = 0 ;
            }
            for ( size_t j = 0 ; j < stressMatrix.numRows() ; j++ )
            {
                if ( j == k )
                {
                    continue ;
                }
                v_reduced[j] -= m_reduced[j][k] ;
                m_reduced[j][k] = 0 ;
            }
        }
        solveSystem ( m_reduced, v_reduced, principalVectors.back() );
    }
    std::vector<Point> ret ;


    if ( stressMatrix.numCols() == 2 )
    {
        for ( size_t i = 0 ; i <  principalVectors.size() ; ++i )
        {
            ret.push_back ( Point ( principalVectors[i][0], principalVectors[i][1] ) );
            ret.back() /= ret.back().norm() ;
        }
        ret.push_back ( Point ( 0, 0, 1 ) );
    }
    else
    {
        for ( size_t i = 0 ; i <  principalVectors.size() ; ++i )
        {
            ret.push_back ( Point ( principalVectors[i][0], principalVectors[i][1], principalVectors[i][2] ) );
            ret.back() /= ret.back().norm() ;
        }
    }
    if ( cleanup )
    {
        delete vm ;
    }
    return ret ;
}


ElementStateWithInternalVariables::ElementStateWithInternalVariables ( IntegrableEntity * e, int n_, int p_ ) : ElementState ( e ), n ( n_ ), p ( p_ )
{

}


void ElementStateWithInternalVariables::getField ( FieldType f, const Point & p, Vector & ret, bool local, VirtualMachine * vm, int i ) 
{
    if ( f != INTERNAL_VARIABLE_FIELD )
    {
        ElementState::getField ( f, p, ret, local, vm, i ) ;
    }
}

void ElementStateWithInternalVariables::getField ( FieldType f, const PointArray & p, Vector & ret, bool local, VirtualMachine * vm, int i ) 
{
    if ( f != INTERNAL_VARIABLE_FIELD )
    {
        ElementState::getField ( f, p, ret, local, vm, i ) ;
    }
}

void ElementStateWithInternalVariables::getField ( FieldType f, const std::valarray<std::pair<Point, double> > & p, Vector & ret, bool local, VirtualMachine * vm, int i ) 
{
    if ( f != INTERNAL_VARIABLE_FIELD )
    {
        ElementState::getField ( f, p, ret, local, vm, i ) ;
    }
}

// void ElementStateWithInternalVariables::getFieldAtNodes( FieldType f, Vector & ret, int i)
// {
// 	if(f != INTERNAL_VARIABLE_FIELD)
// 		ElementState::getFieldAtNodes(f, ret, i) ;
// }

void ElementStateWithInternalVariables::getFieldAtGaussPoint ( FieldType f, size_t g, Vector & ret, VirtualMachine * vm, int i )
{
    if ( f == INTERNAL_VARIABLE_FIELD )
    {
        ret = internalVariablesAtGaussPoints[g][i] ;
        return ;
    }
    else
    {
        ElementState::getFieldAtGaussPoint ( f, g, ret, vm, i ) ;
// 		Point p_ = parent->getGaussPoints().gaussPoints[g].first ;
// 		this->getField(f, p_, ret, true, i) ;
    }
}

void ElementStateWithInternalVariables::getField ( FieldType f1, FieldType f2, const Point & p, Vector & ret1, Vector & ret2, bool local, VirtualMachine *vm, int i, int j ) 
{
    if ( f1 != INTERNAL_VARIABLE_FIELD && f2 != INTERNAL_VARIABLE_FIELD )
    {
        ElementState::getField ( f1, f2, p, ret1, ret2, local, vm, i, j ) ;
    }
}

void ElementStateWithInternalVariables::getField ( FieldType f1, FieldType f2, const PointArray & p,  Vector & ret1, Vector & ret2, bool local, VirtualMachine *vm, int i, int j ) 
{
    if ( f1 != INTERNAL_VARIABLE_FIELD && f2 != INTERNAL_VARIABLE_FIELD )
    {
        ElementState::getField ( f1, f2, p, ret1, ret2, local, vm, i, j ) ;
    }
}

void ElementStateWithInternalVariables::getField ( FieldType f1, FieldType f2, const std::valarray<std::pair<Point, double> > & p,  Vector & ret1, Vector & ret2, bool local, VirtualMachine *vm, int i, int j ) 
{
    if ( f1 != INTERNAL_VARIABLE_FIELD && f2 != INTERNAL_VARIABLE_FIELD )
    {
        ElementState::getField ( f1, f2, p, ret1, ret2, local, vm, i, j ) ;
    }
}

// void ElementStateWithInternalVariables::getFieldAtNodes(FieldType f1, FieldType f2,  Vector & ret1, Vector & ret2, int i, int j)
// {
// 	if(f1 != INTERNAL_VARIABLE_FIELD && f2 != INTERNAL_VARIABLE_FIELD)
// 		ElementState::getFieldAtNodes(f1, f2, ret1, ret2, i, j) ;
// }

void ElementStateWithInternalVariables::getFieldAtGaussPoint ( FieldType f1, FieldType f2, size_t g,  Vector & ret1, Vector & ret2, VirtualMachine * vm, int i, int j )
{
    bool cleanup = !vm ;
    if ( !vm )
    {
        vm = new VirtualMachine() ;
    }
    bool done1 = false ;
    bool done2 = false ;
    if ( f1 == INTERNAL_VARIABLE_FIELD )
    {
        ret1 = internalVariablesAtGaussPoints[g][i] ;
        done1 = true ;
    }
    if ( f2 == INTERNAL_VARIABLE_FIELD )
    {
        ret2 = internalVariablesAtGaussPoints[g][j] ;
        done2 = true ;
    }
    if ( done1 && done2 )
    {
        return ;
    }
    Point p_ = parent->getGaussPoints().gaussPoints[g].first ;
    if ( done1 )
    {
        ElementState::getField ( f2, p_, ret2, true, vm, j ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    if ( done2 )
    {
        ElementState::getField ( f1, p_, ret1, true, vm, i ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    ElementState::getField ( f1, f2, p_, ret1, ret2, true, vm, i, j ) ;
    if ( cleanup )
    {
        delete vm ;
    }
}

void ElementStateWithInternalVariables::initialize (  Mesh<DelaunayTriangle,DelaunayTreeItem> * msh )
{
    ElementState::initialize ( msh) ;

    size_t ngp = parent->getGaussPoints().gaussPoints.size() ;
    std::vector<Vector> dummy ( 0 ) ;
    Vector vdummy ( 0., p ) ;

    for ( size_t g = 0 ; g < ngp ; g++ )
    {
        internalVariablesAtGaussPoints.push_back ( dummy ) ;
        for ( int k = 0 ; k < n ; k++ )
        {
            internalVariablesAtGaussPoints[g].push_back ( vdummy ) ;
        }
    }
}

void ElementStateWithInternalVariables::initialize ( Mesh<DelaunayTetrahedron,DelaunayTreeItem3D> * msh)
{
    ElementState::initialize ( msh ) ;

    size_t ngp = parent->getGaussPoints().gaussPoints.size() ;
    std::vector<Vector> dummy ( 0 ) ;
    Vector vdummy ( 0., p ) ;

    for ( size_t g = 0 ; g < ngp ; g++ )
    {
        internalVariablesAtGaussPoints.push_back ( dummy ) ;
        for ( int k = 0 ; k < n ; k++ )
        {
            internalVariablesAtGaussPoints[g].push_back ( vdummy ) ;
        }
    }
}

void ElementStateWithInternalVariables::setNumberOfGaussPoints ( int g )
{
    while ( (int)internalVariablesAtGaussPoints.size() < g )
    {
        std::vector<Vector> dummy ( 0 ) ;
        for ( int i = 0 ; i < n ; i++ )
        {
            dummy.push_back ( internalVariablesAtGaussPoints[0][i] ) ;
        }
        internalVariablesAtGaussPoints.push_back ( dummy ) ;
    }
    while ( (int)internalVariablesAtGaussPoints.size() > g )
    {
        internalVariablesAtGaussPoints.pop_back() ;
    }
}


void ElementStateWithInternalVariables::setInternalVariableAtGaussPoint ( Vector & v, size_t g, int i )
{
    internalVariablesAtGaussPoints[g][i] = v ;
}

int isGaussPoint ( const Point & p, IntegrableEntity * e )
{
    if ( e )
    {
        GaussPointArray gauss = e->getGaussPoints() ;
        for ( size_t i = 0 ; i < gauss.gaussPoints.size() ; i++ )
        {
            if ( p == gauss.gaussPoints[i].first )
            {
                return i ;
            }
        }
    }
    return -1 ;
}


Vector Form::getForcesFromAppliedStress ( const Vector & data, const Function & shape, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, const std::vector<Variable> & v, bool isVolumic, const Vector & normal, VirtualMachine *vm) const
{
    if ( isVolumic )
    {
      if(vm)
	return vm->ieval ( Gradient ( shape ) * data, gp, Jinv, v )  ;
      return VirtualMachine().ieval ( Gradient ( shape ) * data, gp, Jinv, v )  ;
    }
    Vector ret ( 0., normal.size() ) ;
    if ( normal.size() == 2 )
    {

        ret[0] = ( data[0]*normal[0]+data[2]*normal[1] ) ;
        ret[1] = ( data[2]*normal[0]+data[1]*normal[1] ) ;
    }
    else
    {
        ret[0] += ( data[0]*normal[0]+data[3]*normal[1]+data[4]*normal[2] ) ;
        ret[1] += ( data[1]*normal[1]+data[3]*normal[0]+data[5]*normal[2] ) ;
        ret[2] += ( data[2]*normal[2]+data[4]*normal[0]+data[5]*normal[1] ) ;
    }

    return ret * VirtualMachine().ieval ( shape, gp ) ;
}

Vector Form::getForcesFromAppliedStress ( const Function & data, size_t index, size_t externaldofs,  Function & shape, IntegrableEntity * e,const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<Variable> & v, bool isVolumic , const Vector& normal , VirtualMachine *vm )
{
    size_t n = e->getBoundingPoints().size() ;
    Vector field ( 0., n*externaldofs ) ;
    if(vm)
    for ( size_t i = 0 ; i < n ; i++ )
    {
      
        field[ i*externaldofs + index ] = vm->eval ( data, e->getBoundingPoint ( i ) ) ;
    }
    else
    {
      VirtualMachine vml ;
      for ( size_t i = 0 ; i < n ; i++ )
      {
	
	  field[ i*externaldofs + index ] = vml.eval ( data, e->getBoundingPoint ( i ) ) ;
      }
    }

    std::vector<Vector> g ( gp.gaussPoints.size(), Vector ( 0., externaldofs ) ) ;
    e->getState().getExternalFieldAtGaussPoints ( field, externaldofs, g ) ;

    if ( isVolumic )
    {
      if(vm)
        return vm->ieval ( Gradient ( shape ) * g, gp, Jinv, v ) ;
      else
	return VirtualMachine().ieval ( Gradient ( shape ) * g, gp, Jinv, v ) ;
    }

    Vector ret ( 0., normal.size() ) ;
    if ( normal.size() == 2 )
    {
        for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
        {
            ret[0] += ( g[0][i]*normal[0]+g[2][i]*normal[1] ) *gp.gaussPoints[i].second ;
            ret[1] += ( g[1][i]*normal[1]+g[2][i]*normal[1] ) *gp.gaussPoints[i].second ;
        }
    }
    else
    {
        for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
        {
            ret[0] += ( g[0][i]*normal[0]+g[3][i]*normal[1]+g[4][i]*normal[2] ) *gp.gaussPoints[i].second ;
            ret[1] += ( g[1][i]*normal[1]+g[3][i]*normal[0]+g[5][i]*normal[2] ) *gp.gaussPoints[i].second ;
            ret[2] += ( g[2][i]*normal[2]+g[4][i]*normal[0]+g[5][i]*normal[1] ) *gp.gaussPoints[i].second ;
        }
    }
    if(vm)
      return ret * vm->ieval ( shape, gp ) ;
    return ret * VirtualMachine().ieval ( shape, gp ) ;
}

void IntegrableEntity::setCachedGaussPoints ( GaussPointArray * gp )
{
    if ( cachedGps && cachedGps != gp )
    {
        delete cachedGps ;
    }
    cachedGps = gp ;
}

}
