#include "generalized_spacetime_viscoelastic_element_state.h"
#include "../physics/viscoelasticity.h"
#include "../physics/dual_behaviour.h"
#include "../physics/damagemodels/damagemodel.h"
#include "../mesher/delaunay.h"
#include "../utilities/random.h"
#include "../features/inclusion.h"
#include <omp.h>

using namespace Amie ;

GeneralizedSpaceTimeViscoElasticElementState::GeneralizedSpaceTimeViscoElasticElementState ( IntegrableEntity * e, int blocks ) : ElementState ( e )
{
    strainBuffer.resize( e->spaceDimensions() == SPACE_THREE_DIMENSIONAL ? 6 : 3 );
    principalBuffer.resize( e->spaceDimensions() == SPACE_THREE_DIMENSIONAL ? 3 : 2 );
    generalizedBuffer.resize( strainBuffer.size() * blocks ) ;
    generalizedBufferSecond.resize( generalizedBuffer.size() ) ;
}

GeneralizedSpaceTimeViscoElasticElementState::GeneralizedSpaceTimeViscoElasticElementState ( GeneralizedSpaceTimeViscoElasticElementState &s ) : ElementState ( s )
{

}


GeneralizedSpaceTimeViscoElasticElementState & GeneralizedSpaceTimeViscoElasticElementState::operator = ( GeneralizedSpaceTimeViscoElasticElementState & s )
{
    ElementState::operator = ( s ) ;

    return *this ;
}

GaussPointArray GeneralizedSpaceTimeViscoElasticElementState::genEquivalentGaussPointArray2D ( TriElement * trg, double time )
{
    GaussPointArray gp ;
    switch ( trg->getOrder() )
    {
    case LINEAR_TIME_LINEAR:
    case LINEAR_TIME_QUADRATIC:
        gp = gaussPointSet ( LINEAR, trg ) ;
        break ;
    case QUADRATIC_TIME_LINEAR:
    case QUADRATIC_TIME_QUADRATIC:
        gp = gaussPointSet ( QUADRATIC, trg ) ;
        break ;
    default:
        std::cout << "order not handled in switch" << std::endl ;
        exit(0) ;
        return gp ;
    }
    for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
    {
        gp.gaussPoints[i].first.getT() = time ;
    }

    if ( trg->getEnrichmentFunctions().size() > 0 )
    {
        std::vector<std::pair<Point, double> > gp_alternative ;
        for ( size_t i = 0 ; i < trg->getGaussPoints().gaussPoints.size() /3 ; i++ )
        {
            gp_alternative.push_back ( trg->getGaussPoints().gaussPoints[i*3] );
            gp_alternative.back().first.getT() = time ;
            gp_alternative.back().second *= 6. ;
        }

        gp.gaussPoints.resize ( gp_alternative.size() ) ;
        std::copy ( gp_alternative.begin(), gp_alternative.end(), &gp.gaussPoints[0] );

    }

    return gp ;
}

GaussPointArray  GeneralizedSpaceTimeViscoElasticElementState::genEquivalentGaussPointArray3D ( TetrahedralElement * tet, double time )
{
    GaussPointArray gp ;
    switch ( tet->getOrder() )
    {
    case LINEAR_TIME_LINEAR:
    case LINEAR_TIME_QUADRATIC:
        gp = gaussPointSet ( LINEAR, tet ) ;
        break ;
    case QUADRATIC_TIME_LINEAR:
    case QUADRATIC_TIME_QUADRATIC:
        gp = gaussPointSet ( QUADRATIC, tet ) ;
        break ;
    default:
        std::cout << "order not handled in switch" << std::endl ;
        exit(0) ;
        return gp ;
    }
    for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
    {
        gp.gaussPoints[i].first.getT() = time ;
    }

    if ( tet->getEnrichmentFunctions().size() > 0 )
    {
        std::vector<std::pair<Point, double> > gp_alternative ;
        for ( size_t i = 0 ; i <  tet->getGaussPoints().gaussPoints.size() /3 ; i++ )
        {
            gp_alternative.push_back ( tet->getGaussPoints().gaussPoints[i*3] );
            gp_alternative.back().first.getT() = time ;
            gp_alternative.back().second *= 6. ;
        }

        gp.gaussPoints.resize ( gp_alternative.size() ) ;
        std::copy ( gp_alternative.begin(), gp_alternative.end(), &gp.gaussPoints[0] );

    }

    return gp ;
}

double GeneralizedSpaceTimeViscoElasticElementState::getAverageField ( FieldType f, Vector & ret, VirtualMachine * vm, int dummy , double t, std::vector<double> weights )
{
//     std::cout << "plouf" << std::endl ;
    #pragma omp critical
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
    GaussPointArray gp = parent->getGaussPoints() ;
    if(ret.size() == 0)
        ret.resize( fieldTypeElementarySize( f, parent->spaceDimensions(), parent->getBehaviour()->getNumberOfDegreesOfFreedom()/parent->spaceDimensions() ) ) ;
    ret = 0 ;

    double total = 0. ;
    if ( dummy<0 )
    {
        if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
        {
            gp = GeneralizedSpaceTimeViscoElasticElementState::genEquivalentGaussPointArray2D ( dynamic_cast<TriElement *> ( parent ), t ) ;
        }
        else
        {
            gp = GeneralizedSpaceTimeViscoElasticElementState::genEquivalentGaussPointArray3D ( dynamic_cast<TetrahedralElement *> ( parent ), t ) ;
        }
    }
    if ( !vm )
    {
        vm = new VirtualMachine() ;
    }
    if ( weights.size() != gp.gaussPoints.size() )
    {
        weights.resize ( gp.gaussPoints.size(), 1. );
    }
    Vector tmp ( 0., ret.size() ) ;
    for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
    {
        Point p_ = gp.gaussPoints[i].first ;
        double w = gp.gaussPoints[i].second*weights[i] ;
        bool cached = false ;

        if(dummy < 0 && (f == GENERALIZED_VISCOELASTIC_STRAIN_FIELD || f == GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD || f == STRAIN_FIELD || f == PRINCIPAL_STRAIN_FIELD || f == STRAIN_RATE_FIELD || f == REAL_STRESS_FIELD || f == PRINCIPAL_REAL_STRESS_FIELD || f == EFFECTIVE_STRESS_FIELD || f == PRINCIPAL_EFFECTIVE_STRESS_FIELD || f == MECHANICAL_STRAIN_FIELD || f == PRINCIPAL_MECHANICAL_STRAIN_FIELD) )
        {
            if(std::abs ( t-1 ) < POINT_TOLERANCE)
            {
                tmp = getCachedFieldAtGaussPointAfter( f, gp, i, vm) ;
                cached = true ;
            }
            if(std::abs ( t+1 ) < POINT_TOLERANCE)
            {
                tmp = getCachedFieldAtGaussPointBefore( f, gp, i, vm) ;
                cached = true ;
            }
        }

        if(!cached)
        {
//             std::cout << "pliff" << std::endl ;
            getField ( f, p_, tmp, true, vm, dummy ) ;
            
        }

        ret += tmp * w ;
        total += w ;
    }
    if ( cleanup )
    {
        delete vm ;
    }

    if(total < POINT_TOLERANCE)
        ret = 0. ;
    else
        ret /= total;

    #pragma omp atomic write
    lock = false ;

    /*	if(dummy == 0)
    	{
    		int * p = nullptr ;
    		*p = 3 ;
    	}*/

    if(false) //weights.size() > 200 && abs(ret[1]) > 1e7)
    {
    	std::cout << "enrichment=" << parent->getEnrichmentFunctions().size() << std::endl ;
        std::cout << "displacements=" ;
	for(size_t i = 0 ; i < enrichedDisplacements.size() ; i++)
		std::cout << enrichedDisplacements[i] << "," ;
 	std::cout << std::endl ;
    	std::cout << "weight=" << weights.size() << std::endl ;
    	std::cout << "time=" << t << std::endl ;
    	std::cout << "field=" << f << std::endl ;
    	std::cout << "ret[0]=" << ret[0] << std::endl ;
    	std::cout << "ret[1]=" << ret[1] << std::endl ;
    	std::cout << "area=" <<(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL ? parent->area() : parent->volume()) << std::endl ;
    	std::cout << "total=" << total << std::endl ;
    }

    return total ;//parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL ? parent->area() : parent->volume() ;
}

void  GeneralizedSpaceTimeViscoElasticElementState::step ( double dt, const Vector *d )
{
    averagestressbefore.resize ( 0 );
    averagestressafter.resize ( 0 );

    averagestrainbefore.resize ( 0 );
    averagestrainafter.resize ( 0 );

    averagestrainratebefore.resize ( 0 );
    averagestrainrateafter.resize ( 0 );

    genStrainAtGaussPointBefore.resize ( 0 ) ;
    genStrainAtGaussPointAfter.resize ( 0 ) ;
    genStrainRateAtGaussPointBefore.resize ( 0 ) ;
    genStrainRateAtGaussPointAfter.resize ( 0 ) ;

    ElementState::step ( dt,d );
}

Vector GeneralizedSpaceTimeViscoElasticElementState::getCachedFieldAtGaussPointBefore( FieldType f, GaussPointArray & gp, size_t i, VirtualMachine * vm)
{
    Vector ret ;
    size_t size = gp.gaussPoints.size() * fieldTypeElementarySize( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, parent->spaceDimensions(), parent->getBehaviour()->getNumberOfDegreesOfFreedom()/parent->spaceDimensions()) ;
    size_t blocksize = size/gp.gaussPoints.size() ;
    if(f == GENERALIZED_VISCOELASTIC_STRAIN_FIELD)
    {
        // build cache if it does not exist
        if(genStrainAtGaussPointBefore.size() != size )
        {
            genStrainAtGaussPointBefore.resize( gp.gaussPoints.size() * fieldTypeElementarySize( f, parent->spaceDimensions(), parent->getBehaviour()->getNumberOfDegreesOfFreedom()/parent->spaceDimensions())) ;
            for(size_t j = 0 ; j < gp.gaussPoints.size() ; j++)
            {
                generalizedBuffer = 0. ;
                getField( f, gp.gaussPoints[j].first, generalizedBuffer, true, vm ) ;
                for(size_t k = 0 ; k < generalizedBuffer.size() ; k++)
                    genStrainAtGaussPointBefore[ j*blocksize + k] = generalizedBuffer[k] ;
            }
        }
        ret.resize( genStrainAtGaussPointBefore.size()/gp.gaussPoints.size() ) ;
        // find result in cache
        for(size_t j = 0 ; j < ret.size() ; j++)
            ret[j] = genStrainAtGaussPointBefore[ i*blocksize + j] ;
        return ret ;
    }
    if(f == STRAIN_FIELD)
    {
        // build cache if it does not exist
        if(genStrainAtGaussPointBefore.size() != size)
        {
            genStrainAtGaussPointBefore.resize( gp.gaussPoints.size() * fieldTypeElementarySize( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, parent->spaceDimensions(), parent->getBehaviour()->getNumberOfDegreesOfFreedom()/parent->spaceDimensions())) ;
            for(size_t j = 0 ; j < gp.gaussPoints.size() ; j++)
            {
                generalizedBuffer = 0. ;
                getField( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, gp.gaussPoints[j].first, generalizedBuffer, true, vm ) ;
                for(size_t k = 0 ; k < generalizedBuffer.size() ; k++)
                    genStrainAtGaussPointBefore[ j*blocksize + k] = generalizedBuffer[k] ;
            }
        }
        ret.resize( 3+(3*parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL) ) ;
        // find result in cache
        for(size_t j = 0 ; j < ret.size() ; j++)
            ret[j] = genStrainAtGaussPointBefore[ i*blocksize + j] ;
        return ret ;
    }
    if(f == MECHANICAL_STRAIN_FIELD)
    {
        generalizedBuffer = getCachedFieldAtGaussPointBefore( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, gp, i, vm) ;
        ret = getCachedFieldAtGaussPointBefore( STRAIN_FIELD, gp, i, vm) ;
        for(size_t j = 1 ; j < generalizedBuffer.size() / ret.size() ; j++)
        {
            for(size_t n = 0 ; n < ret.size() ; n++)
                ret[n] -= generalizedBuffer[j*ret.size()+n] ;
        }
        if( getParent()->getBehaviour() && getParent()->getBehaviour()->hasInducedForces() )
            ret -= getParent()->getBehaviour()->getImposedStrain( gp.gaussPoints[i].first ) ;
        return ret ;
    }
    if(f == PRINCIPAL_MECHANICAL_STRAIN_FIELD)
    {
        Vector strain = getCachedFieldAtGaussPointBefore( MECHANICAL_STRAIN_FIELD, gp, i, vm) ;
        return toPrincipal(strain, DOUBLE_OFF_DIAGONAL_VALUES) ;
    }
    if(f == GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD)
    {
        // build cache if it does not exist
        if(genStrainRateAtGaussPointBefore.size() != size)
        {
            genStrainRateAtGaussPointBefore.resize( gp.gaussPoints.size() * fieldTypeElementarySize( f, parent->spaceDimensions(), parent->getBehaviour()->getNumberOfDegreesOfFreedom()/parent->spaceDimensions())) ;
            for(size_t j = 0 ; j < gp.gaussPoints.size() ; j++)
            {
                generalizedBuffer = 0. ;
                getField( f, gp.gaussPoints[j].first, generalizedBuffer, true, vm ) ;
                for(size_t k = 0 ; k < generalizedBuffer.size() ; k++)
                    genStrainRateAtGaussPointBefore[ j*blocksize + k] = generalizedBuffer[k] ;
            }
        }
        ret.resize( genStrainRateAtGaussPointBefore.size()/gp.gaussPoints.size() ) ;
        // find result in cache
        for(size_t j = 0 ; j < ret.size() ; j++)
            ret[j] = genStrainRateAtGaussPointBefore[ i*blocksize + j] ;
        return ret ;
    }
    if(f == STRAIN_RATE_FIELD)
    {
        // build cache if it does not exist
        if(genStrainRateAtGaussPointBefore.size() != size)
        {
            genStrainRateAtGaussPointBefore.resize( gp.gaussPoints.size() * fieldTypeElementarySize( f, parent->spaceDimensions(), parent->getBehaviour()->getNumberOfDegreesOfFreedom()/parent->spaceDimensions())) ;
            for(size_t j = 0 ; j < gp.gaussPoints.size() ; j++)
            {
                generalizedBuffer = 0. ;
                getField( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, gp.gaussPoints[j].first, generalizedBuffer, true, vm ) ;
                for(size_t k = 0 ; k < generalizedBuffer.size() ; k++)
                    genStrainRateAtGaussPointBefore[ j*blocksize + k] = generalizedBuffer[k] ;
            }
        }
        ret.resize( 3+(3*parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL) ) ;
        // find result in cache
        for(size_t j = 0 ; j < ret.size() ; j++)
            ret[j] = genStrainRateAtGaussPointBefore[ i*blocksize + j] ;
        return ret ;
    }
    if(f == REAL_STRESS_FIELD)
    {
        generalizedBuffer = getCachedFieldAtGaussPointBefore( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, gp, i, vm) ;
        generalizedBufferSecond = (Vector) ((parent->getBehaviour()->getTensor( gp.gaussPoints[i].first ))*generalizedBuffer) ;
        generalizedBuffer = getCachedFieldAtGaussPointBefore( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, gp, i, vm) ;
        generalizedBufferSecond += (Vector) ((parent->getBehaviour()->getViscousTensor( gp.gaussPoints[i].first ))*generalizedBuffer) ;
        strainBuffer = parent->getBehaviour()->getImposedStress( gp.gaussPoints[i].first ) ;
        ret.resize(strainBuffer.size()) ;
        for(size_t j = 0 ; j < ret.size() ; j++)
            ret[j] = generalizedBufferSecond[j]-strainBuffer[j] ;
        return ret ;
    }
    if(f == EFFECTIVE_STRESS_FIELD)
    {
        generalizedBuffer = getCachedFieldAtGaussPointBefore( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, gp, i, vm) ;
        generalizedBufferSecond = (Vector) ((parent->getBehaviour()->param)*generalizedBuffer) ;
        generalizedBuffer = getCachedFieldAtGaussPointBefore( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, gp, i, vm) ;
        generalizedBufferSecond += (Vector) ((dynamic_cast<Viscoelasticity *>(parent->getBehaviour())->eta)*generalizedBuffer) ;
        strainBuffer = parent->getBehaviour()->getImposedStress( gp.gaussPoints[i].first ) ;
        ret.resize(strainBuffer.size()) ;
        for(size_t j = 0 ; j < ret.size() ; j++)
            ret[j] = generalizedBufferSecond[j]-strainBuffer[j] ;
        return ret ;
    }
    if(f == PRINCIPAL_STRAIN_FIELD)
    {
        strainBuffer = getCachedFieldAtGaussPointBefore( STRAIN_FIELD, gp, i, vm) ;
        return toPrincipal(strainBuffer , DOUBLE_OFF_DIAGONAL_VALUES) ;
    }
    if(f == PRINCIPAL_REAL_STRESS_FIELD)
    {
        strainBuffer = getCachedFieldAtGaussPointBefore( REAL_STRESS_FIELD, gp, i, vm) ;
        return toPrincipal(strainBuffer , SINGLE_OFF_DIAGONAL_VALUES) ;
    }
    if(f == PRINCIPAL_EFFECTIVE_STRESS_FIELD)
    {
        strainBuffer = getCachedFieldAtGaussPointBefore( EFFECTIVE_STRESS_FIELD, gp, i, vm) ;
        return toPrincipal(strainBuffer , SINGLE_OFF_DIAGONAL_VALUES) ;
    }
    return ret ;
}

void GeneralizedSpaceTimeViscoElasticElementState::getCache( FieldType f, GaussPointArray & gp, Vector & ret, VirtualMachine * vm) 
{
    size_t size = gp.gaussPoints.size() * fieldTypeElementarySize( f, parent->spaceDimensions(), parent->getBehaviour()->getNumberOfDegreesOfFreedom()/parent->spaceDimensions()) ;
    size_t blocksize = size/gp.gaussPoints.size() ;
    if(ret.size() != size)
        ret.resize( size ) ;
    for(size_t j = 0 ; j < gp.gaussPoints.size() ; j++)
    {
        generalizedBuffer = 0. ;
        getField( f, gp.gaussPoints[j].first, generalizedBuffer, true, vm ) ;
        for(size_t k = 0 ; k < generalizedBuffer.size() ; k++)
             ret[ j*blocksize + k] = generalizedBuffer[k] ;
    }
}

Vector GeneralizedSpaceTimeViscoElasticElementState::getCachedFieldAtGaussPointAfter( FieldType f, GaussPointArray & gp, size_t i, VirtualMachine * vm)
{
    Vector ret ;
    size_t size = gp.gaussPoints.size() * fieldTypeElementarySize( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, parent->spaceDimensions(), parent->getBehaviour()->getNumberOfDegreesOfFreedom()/parent->spaceDimensions()) ;
    size_t blocksize = size/gp.gaussPoints.size() ;
    if(f == GENERALIZED_VISCOELASTIC_STRAIN_FIELD)
    {
        // build cache if it does not exist
        if(genStrainAtGaussPointAfter.size() != size)
            getCache( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, gp, genStrainAtGaussPointAfter, vm) ;
        ret.resize( blocksize ) ;
        // find result in cache
        for(size_t j = 0 ; j < ret.size() ; j++)
            ret[j] = genStrainAtGaussPointAfter[ i*blocksize + j] ;
        return ret ;
    }
    if(f == STRAIN_FIELD)
    {
        // build cache if it does not exist
        if(genStrainAtGaussPointAfter.size() != size)
            getCache( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, gp, genStrainAtGaussPointAfter, vm) ;
        ret.resize( (3+(3*parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)) ) ;
        // find result in cache
        for(size_t j = 0 ; j < ret.size() ; j++)
            ret[j] = genStrainAtGaussPointAfter[ i*blocksize + j] ;
        return ret ;
    }
    if(f == MECHANICAL_STRAIN_FIELD)
    {
        generalizedBuffer = getCachedFieldAtGaussPointAfter( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, gp, i, vm) ;
        ret = getCachedFieldAtGaussPointAfter( STRAIN_FIELD, gp, i, vm) ;
        for(size_t j = 1 ; j < generalizedBuffer.size() / ret.size() ; j++)
        {
            for(size_t n = 0 ; n < ret.size() ; n++)
                ret[n] -= generalizedBuffer[j*ret.size()+n] ;
        }
        if( getParent()->getBehaviour() && getParent()->getBehaviour()->hasInducedForces() )
            ret -= getParent()->getBehaviour()->getImposedStrain( gp.gaussPoints[i].first ) ;
        return ret ;
    }
    if(f == PRINCIPAL_MECHANICAL_STRAIN_FIELD)
    {
        strainBuffer = getCachedFieldAtGaussPointAfter( MECHANICAL_STRAIN_FIELD, gp, i, vm) ;
        return toPrincipal(strainBuffer , DOUBLE_OFF_DIAGONAL_VALUES) ;
    }
    if(f == GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD)
    {
        // build cache if it does not exist
        if(genStrainRateAtGaussPointAfter.size() != size)
            getCache( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, gp, genStrainRateAtGaussPointAfter, vm) ;
        ret.resize( genStrainRateAtGaussPointAfter.size()/gp.gaussPoints.size() ) ;
        // find result in cache
        for(size_t j = 0 ; j < ret.size() ; j++)
            ret[j] = genStrainRateAtGaussPointAfter[ i*blocksize + j] ;
        return ret ;
    }
    if(f == STRAIN_RATE_FIELD)
    {
        // build cache if it does not exist
        if(genStrainRateAtGaussPointAfter.size() != size)
            getCache( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, gp, genStrainRateAtGaussPointAfter, vm) ;
        ret.resize( 3+(3*parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL) ) ;
        // find result in cache
        for(size_t j = 0 ; j < ret.size() ; j++)
            ret[j] = genStrainRateAtGaussPointAfter[ i*blocksize + j] ;
        return ret ;
    }
    if(f == REAL_STRESS_FIELD)
    {
        generalizedBuffer = getCachedFieldAtGaussPointAfter( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, gp, i, vm) ;
        generalizedBufferSecond = (Vector) ((parent->getBehaviour()->getTensor( gp.gaussPoints[i].first ))*generalizedBuffer) ;
        generalizedBuffer = getCachedFieldAtGaussPointAfter( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, gp, i, vm) ;
        generalizedBufferSecond += (Vector) ((parent->getBehaviour()->getViscousTensor( gp.gaussPoints[i].first ))*generalizedBuffer) ;
        
        strainBuffer = parent->getBehaviour()->getImposedStress( gp.gaussPoints[i].first ) ;
        ret.resize(strainBuffer.size()) ;
        for(size_t j = 0 ; j < ret.size() ; j++)
            ret[j] = generalizedBufferSecond[j]-strainBuffer[j] ;
        return ret ;
    }
    if(f == EFFECTIVE_STRESS_FIELD)
    {
        generalizedBuffer = getCachedFieldAtGaussPointAfter( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, gp, i, vm) ;
        generalizedBufferSecond = (Vector) ((parent->getBehaviour()->param)*generalizedBuffer) ;
        generalizedBuffer = getCachedFieldAtGaussPointAfter( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, gp, i, vm) ;
        generalizedBufferSecond += (Vector) ((dynamic_cast<Viscoelasticity *>(parent->getBehaviour())->eta)*generalizedBuffer) ;
        
        strainBuffer = parent->getBehaviour()->getImposedStress( gp.gaussPoints[i].first ) ;
        ret.resize(strainBuffer.size()) ;
        for(size_t j = 0 ; j < ret.size() ; j++)
            ret[j] = generalizedBufferSecond[j]-strainBuffer[j] ;
        return ret ;
    }
    if(f == PRINCIPAL_STRAIN_FIELD)
    {
        strainBuffer = getCachedFieldAtGaussPointAfter( STRAIN_FIELD, gp, i, vm) ;
        return toPrincipal(strainBuffer, DOUBLE_OFF_DIAGONAL_VALUES) ;
    }
    if(f == PRINCIPAL_REAL_STRESS_FIELD)
    {
        strainBuffer = getCachedFieldAtGaussPointAfter( REAL_STRESS_FIELD, gp, i, vm) ;
        return toPrincipal(strainBuffer, SINGLE_OFF_DIAGONAL_VALUES) ;
    }
    if(f == PRINCIPAL_EFFECTIVE_STRESS_FIELD)
    {
        strainBuffer = getCachedFieldAtGaussPointAfter( EFFECTIVE_STRESS_FIELD, gp, i, vm) ;
        return toPrincipal(strainBuffer, SINGLE_OFF_DIAGONAL_VALUES) ;
    }
    return ret ;
}


void GeneralizedSpaceTimeViscoElasticElementState::getEssentialAverageFields ( FieldType f , Vector & stress, Vector & strain, Vector & strain_rate, VirtualMachine * vm, double t )
{
    bool cleanup = !vm ;
    if ( !vm )
    {
        vm = new VirtualMachine() ;
    }

    if ( std::abs ( t-1 ) < POINT_TOLERANCE && averagestrainafter.size() > 0 )
    {
        stress = averagestressafter ;
        strain = averagestrainafter ;
        strain_rate = averagestrainrateafter ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }

    if ( std::abs ( t+1 ) < POINT_TOLERANCE && averagestrainbefore.size() > 0 )
    {

        stress = averagestressbefore ;
        strain = averagestrainbefore ;
        strain_rate = averagestrainratebefore ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }

    if ( std::abs ( t+1 ) < POINT_TOLERANCE )
    {
        averagestressbefore.resize ( stress.size() );
        averagestrainbefore.resize ( strain.size() );
        averagestrainratebefore.resize ( strain_rate.size() );
    }

    if ( std::abs ( t-1 ) < POINT_TOLERANCE )
    {
        averagestressafter.resize ( stress.size() );
        averagestrainafter.resize ( strain.size() );
        averagestrainrateafter.resize ( strain_rate.size() );
    }
    GaussPointArray gp = parent->getGaussPoints() ;
    stress = 0 ;
    strain = 0 ;
    strain_rate = 0 ;

    double total = 0 ;

    int totaldof = parent->getBehaviour()->getNumberOfDegreesOfFreedom() ;
    int realdof = parent->spaceDimensions() ;
    int blocks = totaldof / realdof ;

    Viscoelasticity * visco = dynamic_cast<Viscoelasticity *>(parent->getBehaviour()) ;

    if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
    {
        gp = GeneralizedSpaceTimeViscoElasticElementState::genEquivalentGaussPointArray2D ( static_cast<TriElement *> ( parent ), t ) ;
    }
    else
    {
        gp = GeneralizedSpaceTimeViscoElasticElementState::genEquivalentGaussPointArray3D ( static_cast<TetrahedralElement *> ( parent ), t ) ;
    }

    Vector tmpstress ( 0., strain.size() ) ;
    Vector tmpstrain ( 0., strain.size() ) ;
    Vector tmpstrainrate ( 0., strain_rate.size() ) ;
    for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ )
    {
        double w = gp.gaussPoints[i].second ;
        if(!JinvCache || parent->isMoved())
        {
            if(JinvCache)
                delete JinvCache ;
            JinvCache = new Matrix ( parent->spaceDimensions()+(parent->timePlanes()>1),parent->spaceDimensions()+(parent->timePlanes()>1)) ;
            parent->getInverseJacobianMatrix ( gp.gaussPoints[i].first, (*JinvCache) ) ;
        }

        if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
        {

            Vector dx ( 0., totaldof ) ;
            Vector dy ( 0., totaldof ) ;
            Vector dt ( 0., totaldof ) ;

            double x_xi = 0;
            double x_eta = 0;
            double y_xi = 0;
            double y_eta = 0;

            for ( size_t j = 0 ; j < parent->getShapeFunctions().size(); j++ )
            {
                double f_xi = vm->deval ( parent->getShapeFunction ( j ), XI, gp.gaussPoints[i].first ) ;
                double f_eta = vm->deval ( parent->getShapeFunction ( j ), ETA, gp.gaussPoints[i].first ) ;
                double f_tau = vm->deval ( parent->getShapeFunction ( j ), TIME_VARIABLE, gp.gaussPoints[i].first ) ;
                for ( int k = 0 ; k < totaldof ; k++ )
                {
                    dx[k] += f_xi  * displacements[j * totaldof + k] ;
                    dy[k] += f_eta * displacements[j * totaldof + k] ;
                    dt[k] += f_tau * displacements[j * totaldof + k] ;
                }
            }

            for ( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() && j < enrichedDisplacements.size() * 2; j++ )
            {
                double f_xi = vm->deval ( parent->getEnrichmentFunction ( j ), XI, gp.gaussPoints[i].first ) ;
                double f_eta = vm->deval ( parent->getEnrichmentFunction ( j ), ETA, gp.gaussPoints[i].first ) ;
                double f_tau = vm->deval ( parent->getEnrichmentFunction ( j ), TIME_VARIABLE, gp.gaussPoints[i].first ) ;
                for ( int k = 0 ; k < totaldof ; k++ )
                {
                    dx[k] += f_xi  * enrichedDisplacements[j * totaldof + k] ;
                    dy[k] += f_eta * enrichedDisplacements[j * totaldof + k] ;
                    dt[k] += f_tau * enrichedDisplacements[j * totaldof + k] ;
                }
            }


            for ( int k = 0 ; k < blocks ; k++ )
            {
                x_xi  = dx[ k * realdof + 0 ] ;
                x_eta = dy[ k * realdof + 0 ] ;
                y_xi  = dx[ k * realdof + 1 ] ;
                y_eta = dy[ k * realdof + 1 ] ;
                tmpstrain[k*3+0] = ( x_xi ) * (*JinvCache) [0][0] + ( x_eta ) * (*JinvCache) [0][1] ;//+ x_tau * Jinv[0][2];
                tmpstrain[k*3+1] = ( y_xi ) * (*JinvCache) [1][0] + ( y_eta ) * (*JinvCache) [1][1] ;//+ y_tau * Jinv[1][2] ;
                tmpstrain[k*3+2] = 0.5 * ( ( x_xi ) * (*JinvCache) [1][0] + ( x_eta ) * (*JinvCache) [1][1]  + ( y_xi ) * (*JinvCache) [0][0] + ( y_eta ) * (*JinvCache) [0][1] );//+ x_tau * Jinv[1][2]  + y_tau * Jinv[0][2]);
            }
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

            Vector dx ( 0., totaldof ) ;
            Vector dy ( 0., totaldof ) ;
            Vector dz ( 0., totaldof ) ;

            for ( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
            {
                double f_xi = vm->deval ( parent->getShapeFunction ( j ), XI, gp.gaussPoints[i].first ) ;
                double f_eta = vm->deval ( parent->getShapeFunction ( j ), ETA, gp.gaussPoints[i].first ) ;
                double f_zeta = vm->deval ( parent->getShapeFunction ( j ), ZETA, gp.gaussPoints[i].first ) ;

                for ( int k = 0 ; k < totaldof ; k++ )
                {
                    dx[k] += f_xi   * displacements[j * totaldof + k] ;
                    dy[k] += f_eta  * displacements[j * totaldof + k] ;
                    dz[k] += f_zeta * displacements[j * totaldof + k] ;
                }

            }

            for ( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
            {
                double f_xi = vm->deval ( parent->getEnrichmentFunction ( j ), XI, gp.gaussPoints[i].first ) ;
                double f_eta = vm->deval ( parent->getEnrichmentFunction ( j ), ETA, gp.gaussPoints[i].first ) ;
                double f_zeta = vm->deval ( parent->getEnrichmentFunction ( j ), ZETA, gp.gaussPoints[i].first ) ;

                for ( int k = 0 ; k < totaldof ; k++ )
                {
                    dx[k] += f_xi   * enrichedDisplacements[j * totaldof + k] ;
                    dy[k] += f_eta  * enrichedDisplacements[j * totaldof + k] ;
                    dz[k] += f_zeta * enrichedDisplacements[j * totaldof + k] ;
                }

            }

            for ( int k = 0 ; k < blocks ; k++ )
            {
                x_xi   = dx[ k * realdof + 0 ] ;
                x_eta  = dy[ k * realdof + 0 ] ;
                x_zeta = dz[ k * realdof + 0 ] ;
                y_xi   = dx[ k * realdof + 1 ] ;
                y_eta  = dy[ k * realdof + 1 ] ;
                y_zeta = dz[ k * realdof + 1 ] ;
                z_xi   = dx[ k * realdof + 2 ] ;
                z_eta  = dy[ k * realdof + 2 ] ;
                z_zeta = dz[ k * realdof + 2 ] ;
                tmpstrain[k*6+0] = ( x_xi ) * (*JinvCache)[0][0] + ( x_eta ) * (*JinvCache)[0][1]  + ( x_zeta ) * (*JinvCache)[0][2];
                tmpstrain[k*6+1] = ( y_xi ) * (*JinvCache)[1][0] + ( y_eta ) * (*JinvCache)[1][1]  + ( y_zeta ) * (*JinvCache)[1][2];
                tmpstrain[k*6+2] = ( z_xi ) * (*JinvCache)[2][0] + ( z_eta ) * (*JinvCache)[2][1]  + ( z_zeta ) * (*JinvCache)[2][2];

                tmpstrain[k*6+3] = 0.5 * ( ( y_xi ) * (*JinvCache)[2][0] +
                                           ( y_eta ) * (*JinvCache)[2][1] +
                                           ( y_zeta ) * (*JinvCache)[2][2] +
                                           ( z_xi ) * (*JinvCache)[1][0] +
                                           ( z_eta ) * (*JinvCache)[1][1] +
                                           ( z_zeta ) * (*JinvCache)[1][2] );

                tmpstrain[k*6+4] = 0.5 * ( ( x_xi ) * (*JinvCache)[2][0] +
                                           ( x_eta ) * (*JinvCache)[2][1] +
                                           ( x_zeta ) * (*JinvCache)[2][2] +
                                           ( z_xi ) * (*JinvCache)[0][0] +
                                           ( z_eta ) * (*JinvCache)[0][1] +
                                           ( z_zeta ) * (*JinvCache)[0][2] );

                tmpstrain[k*6+5] = 0.5 * ( ( y_xi )   * (*JinvCache)[0][0] +
                                           ( y_eta )  * (*JinvCache)[0][1] +
                                           ( y_zeta ) * (*JinvCache)[0][2] +
                                           ( x_xi )   * (*JinvCache)[1][0] +
                                           ( x_eta )  * (*JinvCache)[1][1] +
                                           ( x_zeta ) * (*JinvCache)[1][2] );
            }

        }

        if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
        {

            Vector dx ( 0., totaldof ) ;
            Vector dy ( 0., totaldof ) ;
            Vector dt ( 0., totaldof ) ;

            double x_xi = 0;
            double x_eta = 0;
            double y_xi = 0;
            double y_eta = 0;

            for ( size_t j = 0 ; j < parent->getBoundingPoints().size(); j++ )
            {
                double f_xi = vm->ddeval ( parent->getShapeFunction ( j ), XI, TIME_VARIABLE, gp.gaussPoints[i].first , 1e-12 ) ;
                double f_eta = vm->ddeval ( parent->getShapeFunction ( j ), ETA, TIME_VARIABLE,gp.gaussPoints[i].first , 1e-12 ) ;
                double f_tau = vm->ddeval ( parent->getShapeFunction ( j ), TIME_VARIABLE, TIME_VARIABLE,gp.gaussPoints[i].first , 1e-12 ) ;
                for ( int k = 0 ; k < totaldof ; k++ )
                {
                    dx[k] += f_xi  * displacements[j * totaldof + k] ;
                    dy[k] += f_eta * displacements[j * totaldof + k] ;
                    dt[k] += f_tau * displacements[j * totaldof + k] ;
                }

            }

            for ( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
            {
                double f_xi = vm->ddeval ( parent->getEnrichmentFunction ( j ), XI, TIME_VARIABLE,gp.gaussPoints[i].first , 1e-12 ) ;
                double f_eta = vm->ddeval ( parent->getEnrichmentFunction ( j ), ETA, TIME_VARIABLE,gp.gaussPoints[i].first , 1e-12 ) ;
                double f_tau = vm->ddeval ( parent->getEnrichmentFunction ( j ), TIME_VARIABLE,TIME_VARIABLE, gp.gaussPoints[i].first , 1e-12 ) ;
                for ( int k = 0 ; k < totaldof ; k++ )
                {
                    dx[k] += f_xi  * enrichedDisplacements[j * totaldof + k] ;
                    dy[k] += f_eta * enrichedDisplacements[j * totaldof + k] ;
                    dt[k] += f_tau * enrichedDisplacements[j * totaldof + k] ;
                }

            }


            for ( int k = 0 ; k < blocks ; k++ )
            {
                x_xi  = dx[ k * realdof + 0 ] ;
                x_eta = dy[ k * realdof + 0 ] ;
                y_xi  = dx[ k * realdof + 1 ] ;
                y_eta = dy[ k * realdof + 1 ] ;
                tmpstrainrate[k*3+0] = ( x_xi ) * (*JinvCache)[0][0] + ( x_eta ) * (*JinvCache)[0][1] ;//+ x_tau * Jinv[0][2];
                tmpstrainrate[k*3+1] = ( y_xi ) * (*JinvCache)[1][0] + ( y_eta ) * (*JinvCache)[1][1] ;//+ y_tau * Jinv[1][2] ;
                tmpstrainrate[k*3+2] = 0.5 * ( ( x_xi ) * (*JinvCache)[1][0] + ( x_eta ) * (*JinvCache)[1][1]  + ( y_xi ) * (*JinvCache)[0][0] + ( y_eta ) * (*JinvCache)[0][1] );//+ x_tau * Jinv[1][2]  + y_tau * Jinv[0][2]);
            }
            tmpstrainrate *= (*JinvCache)[2][2] ;

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

            Vector dx ( 0., totaldof ) ;
            Vector dy ( 0., totaldof ) ;
            Vector dz ( 0., totaldof ) ;

            for ( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
            {
                double f_xi = vm->ddeval ( parent->getShapeFunction ( j ), XI, TIME_VARIABLE, gp.gaussPoints[i].first , 1e-5 ) ;
                double f_eta = vm->ddeval ( parent->getShapeFunction ( j ), ETA, TIME_VARIABLE, gp.gaussPoints[i].first , 1e-5 ) ;
                double f_zeta = vm->ddeval ( parent->getShapeFunction ( j ), ZETA, TIME_VARIABLE,  gp.gaussPoints[i].first, 1e-5 ) ;

                for ( int k = 0 ; k < totaldof ; k++ )
                {
                    dx[k] += f_xi   * displacements[j * totaldof + k] ;
                    dy[k] += f_eta  * displacements[j * totaldof + k] ;
                    dz[k] += f_zeta * displacements[j * totaldof + k] ;
                }

            }

            for ( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
            {
                double f_xi = vm->ddeval ( parent->getEnrichmentFunction ( j ), XI, TIME_VARIABLE, gp.gaussPoints[i].first, 1e-5 ) ;
                double f_eta = vm->ddeval ( parent->getEnrichmentFunction ( j ), ETA, TIME_VARIABLE, gp.gaussPoints[i].first , 1e-5 ) ;
                double f_zeta = vm->ddeval ( parent->getEnrichmentFunction ( j ), ZETA, TIME_VARIABLE, gp.gaussPoints[i].first, 1e-5 ) ;

                for ( int k = 0 ; k < totaldof ; k++ )
                {
                    dx[k] += f_xi   * enrichedDisplacements[j * totaldof + k] ;
                    dy[k] += f_eta  * enrichedDisplacements[j * totaldof + k] ;
                    dz[k] += f_zeta * enrichedDisplacements[j * totaldof + k] ;
                }

            }


            for ( int k = 0 ; k < blocks ; k++ )
            {
                x_xi   = dx[ k * realdof + 0 ] ;
                x_eta  = dy[ k * realdof + 0 ] ;
                x_zeta = dz[ k * realdof + 0 ] ;
                y_xi   = dx[ k * realdof + 1 ] ;
                y_eta  = dy[ k * realdof + 1 ] ;
                y_zeta = dz[ k * realdof + 1 ] ;
                z_xi   = dx[ k * realdof + 2 ] ;
                z_eta  = dy[ k * realdof + 2 ] ;
                z_zeta = dz[ k * realdof + 2 ] ;
                tmpstrainrate[k*6+0] = ( x_xi ) * (*JinvCache)[0][0] + ( x_eta ) * (*JinvCache)[0][1]  + ( x_zeta ) * (*JinvCache)[0][2];
                tmpstrainrate[k*6+1] = ( y_xi ) * (*JinvCache)[1][0] + ( y_eta ) * (*JinvCache)[1][1]  + ( y_zeta ) * (*JinvCache)[1][2];
                tmpstrainrate[k*6+2] = ( z_xi ) * (*JinvCache)[2][0] + ( z_eta ) * (*JinvCache)[2][1]  + ( z_zeta ) * (*JinvCache)[2][2];

                tmpstrainrate[k*6+3] = 0.5 * ( ( y_xi ) * (*JinvCache)[2][0] +
                                               ( y_eta ) * (*JinvCache)[2][1] +
                                               ( y_zeta ) * (*JinvCache)[2][2] +
                                               ( z_xi ) * (*JinvCache)[1][0] +
                                               ( z_eta ) * (*JinvCache)[1][1] +
                                               ( z_zeta ) * (*JinvCache)[1][2] );

                tmpstrainrate[k*6+4] = 0.5 * ( ( x_xi ) * (*JinvCache)[2][0] +
                                               ( x_eta ) * (*JinvCache)[2][1] +
                                               ( x_zeta ) * (*JinvCache)[2][2] +
                                               ( z_xi ) * (*JinvCache)[0][0] +
                                               ( z_eta ) * (*JinvCache)[0][1] +
                                               ( z_zeta ) * (*JinvCache)[0][2] );

                tmpstrainrate[k*6+5] = 0.5 * ( ( y_xi )   * (*JinvCache)[0][0] +
                                               ( y_eta )  * (*JinvCache)[0][1] +
                                               ( y_zeta ) * (*JinvCache)[0][2] +
                                               ( x_xi )   * (*JinvCache)[1][0] +
                                               ( x_eta )  * (*JinvCache)[1][1] +
                                               ( x_zeta ) * (*JinvCache)[1][2] );
            }
            tmpstrainrate *= (*JinvCache)[3][3] ;

        }


        if ( f == REAL_STRESS_FIELD )
        {
            tmpstress = ( Vector ) ( parent->getBehaviour()->getTensor ( gp.gaussPoints[i].first, parent ) * tmpstrain )
                        + ( Vector ) ( parent->getBehaviour()->getViscousTensor ( gp.gaussPoints[i].first, parent ) * tmpstrainrate ) ;
            for ( size_t j = 0 ; j < stress.size() ; j++ )
            {
                stress[j] += tmpstress[j]*w ;
            }
        }
        else
        {
            tmpstress = ( Vector ) ( visco->param * tmpstrain )
                        + ( Vector ) ( visco->eta * tmpstrainrate ) ;
            for ( size_t j = 0 ; j < 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ; j++ )
            {
                stress[j] += tmpstress[j]*w ;
            }
        }

        strain += tmpstrain * w ;
        strain_rate += tmpstrainrate * w ;
        stress -= parent->getBehaviour()->getImposedStress ( gp.gaussPoints[i].first, parent ) *w ;

        total += w ;
    }

    stress /= total ;
    strain /= total ;
    strain_rate /= total ;

    if ( std::abs ( t+1 ) < POINT_TOLERANCE )
    {
        averagestressbefore = stress;
        averagestrainbefore = strain;
        averagestrainratebefore = strain_rate;
    }

    if ( std::abs ( t-1 ) < POINT_TOLERANCE )
    {
        averagestressafter = stress;
        averagestrainafter = strain;
        averagestrainrateafter = strain_rate;
    }

    if ( cleanup )
    {
        delete vm ;
    }
}

double GeneralizedSpaceTimeViscoElasticElementState::getAverageField ( FieldType f1, FieldType f2, Vector & r1, Vector & r2, VirtualMachine * vm , int dummy , double t , std::vector<double> weights)
{
    bool cleanup = !vm ;
    if ( !vm )
    {
        vm = new VirtualMachine() ;
    }
    double r = getAverageField ( f1, r1, vm, dummy, t , weights) ;
    getAverageField ( f2, r2, vm, dummy, t , weights) ;
    if ( cleanup )
    {
        delete vm ;
    }

    return r ;
}

void GeneralizedSpaceTimeViscoElasticElementState::getField ( FieldType f, const Point & p, Vector & ret, bool local, VirtualMachine * vm , int )  
{
    ret = 0. ;

    int totaldof = parent->getBehaviour()->getNumberOfDegreesOfFreedom() ;
    int realdof = parent->spaceDimensions() ;
    int blocks = totaldof / realdof ;
    Point p_ = p ;
    if ( !local )
    {
        p_ = parent->inLocalCoordinates ( p ) ;
    }

    Viscoelasticity * visco = dynamic_cast<Viscoelasticity *>( parent->getBehaviour() ) ;
    bool cleanup = !vm ;
    if ( !vm )
    {
        vm = new VirtualMachine() ;
    }
    switch ( f )
    {
    case SCALAR_DAMAGE_FIELD:
    {
        if(parent->getBehaviour()->getDamageModel())
        {
            Matrix trueParam = parent->getBehaviour()->getTensor( p_ ) ;
            ret[0] = 1. - trueParam[0][0]/parent->getBehaviour()->param[0][0] ;
        }
        else
            ret[0] = 0. ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case DISPLACEMENT_FIELD:
    {
        for ( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
        {
            double f =  vm->eval ( parent->getShapeFunction ( j ) , p_ ) ;
            for ( int k = 0 ; k < realdof ; k++ )
            {
                ret[k] += f * displacements[j*totaldof+k] ;
            }
        }
        for ( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
        {
            double f =  vm->eval ( parent->getEnrichmentFunction ( j ) , p_ ) ;
            for ( int k = 0 ; k < realdof ; k++ )
            {
                ret[k] += f * enrichedDisplacements[j*totaldof+k] ;
            }
        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case IMPOSED_STRESS_FIELD:
    {
        if ( parent->getBehaviour()->getTensor ( p_, parent ).numCols() != ret.size() )
        {
            ret = 0 ;
            if ( cleanup )
                delete vm ;
            return ;
        }

        ret = parent->getBehaviour()->getImposedStress ( p_, parent ) ;

        if ( cleanup )
            delete vm ;
        return ;
    }
    case GENERALIZED_VISCOELASTIC_DISPLACEMENT_FIELD:
    {

        for ( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
        {
            double f =  vm->eval ( parent->getShapeFunction ( j ) , p_ ) ;
            for ( int k = 0 ; k < totaldof ; k++ )
            {
                ret[k] += f * displacements[j*totaldof+k] ;
            }
        }
        for ( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
        {
            double f =  vm->eval ( parent->getEnrichmentFunction ( j ) , p_ ) ;
            for ( int k = 0 ; k < totaldof ; k++ )
            {
                ret[k] += f * enrichedDisplacements[j*totaldof+k] ;
            }
        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case ENRICHED_DISPLACEMENT_FIELD:
    {

        for ( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
        {
            double f =  vm->eval ( parent->getEnrichmentFunction ( j ) , p_ ) ;
            for ( int k = 0 ; k < realdof ; k++ )
            {
                ret[k] += f * enrichedDisplacements[j*totaldof+k] ;
            }
        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case GENERALIZED_VISCOELASTIC_ENRICHED_DISPLACEMENT_FIELD:
    {
        for ( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
        {
            double f =  vm->eval ( parent->getEnrichmentFunction ( j ) , p_ ) ;
            for ( int k = 0 ; k < totaldof ; k++ )
            {
                ret[k] += f * enrichedDisplacements[j*totaldof+k] ;
            }
        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case SPEED_FIELD:
    {
        for ( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
        {
            double f =  vm->deval ( parent->getShapeFunction ( j ) , TIME_VARIABLE, p_ ) ;
            for ( int k = 0 ; k < realdof ; k++ )
            {
                ret[k] += f * displacements[j*totaldof+k] ;
            }
        }
        for ( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
        {
            double f =  vm->deval ( parent->getEnrichmentFunction ( j ) , TIME_VARIABLE,  p_ ) ;
            for ( int k = 0 ; k < realdof ; k++ )
            {
                ret[k] += f * enrichedDisplacements[j*totaldof+k] ;
            }
        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case GENERALIZED_VISCOELASTIC_SPEED_FIELD:
    {
        for ( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
        {
            double f =  vm->deval ( parent->getShapeFunction ( j ) , TIME_VARIABLE, p_ ) ;
            for ( int k = 0 ; k < totaldof ; k++ )
            {
                ret[k] += f * displacements[j*totaldof+k] ;
            }
        }
        for ( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
        {
            double f =  vm->deval ( parent->getEnrichmentFunction ( j ) , TIME_VARIABLE, p_ ) ;
            for ( int k = 0 ; k < totaldof ; k++ )
            {
                ret[k] += f * enrichedDisplacements[j*totaldof+k] ;
            }
        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case STRAIN_FIELD:
    {
        if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
        {
            double x_xi = 0;
            double x_eta = 0;
            double x_tau = 0 ;
            double y_xi = 0;
            double y_eta = 0;
            double y_tau = 0 ;

            for ( size_t j = 0 ; j < parent->getBoundingPoints().size(); j++ )
            {
                double f_xi = vm->deval ( parent->getShapeFunction ( j ), XI, p_ ) ;
                double f_eta = vm->deval ( parent->getShapeFunction ( j ), ETA, p_ ) ;
                double f_tau = vm->deval ( parent->getShapeFunction ( j ), TIME_VARIABLE, p_ ) ;
                x_xi += f_xi * displacements[j * totaldof] ;
                x_eta += f_eta * displacements[j * totaldof] ;
                x_tau += f_tau * displacements[j * totaldof] ;
                y_xi += f_xi * displacements[j * totaldof + 1] ;
                y_eta += f_eta * displacements[j * totaldof + 1] ;
                y_tau += f_tau * displacements[j * totaldof + 1] ;
            }

            for ( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
            {
                double f_xi = vm->deval ( parent->getEnrichmentFunction ( j ), XI, p_ ) ;
                double f_eta = vm->deval ( parent->getEnrichmentFunction ( j ), ETA, p_ ) ;

                x_xi += f_xi * enrichedDisplacements[j * totaldof] ;
                x_eta += f_eta * enrichedDisplacements[j * totaldof] ;
                y_xi += f_xi * enrichedDisplacements[j * totaldof + 1] ;
                y_eta += f_eta * enrichedDisplacements[j * totaldof + 1] ;
            }

            if(!JinvCache || parent->isMoved())
            {
                if(JinvCache)
                    delete JinvCache ;
                JinvCache = new Matrix ( parent->spaceDimensions()+(parent->timePlanes()>1),parent->spaceDimensions()+(parent->timePlanes()>1)) ;
                parent->getInverseJacobianMatrix ( p_, (*JinvCache) ) ;
            }
            ret[0] = ( x_xi ) * (*JinvCache)[0][0] + ( x_eta ) * (*JinvCache)[0][1] ;//+ x_tau * Jinv[0][2] ;
            ret[1] = ( y_xi ) * (*JinvCache)[1][0] + ( y_eta ) * (*JinvCache)[1][1] ;//+ y_tau * Jinv[1][2] ;
            ret[2] = 0.5 * ( ( x_xi ) * (*JinvCache)[1][0] + ( x_eta ) * (*JinvCache)[1][1]  + ( y_xi ) * (*JinvCache)[0][0] + ( y_eta ) * (*JinvCache)[0][1] ) ;// + x_tau * Jinv[1][2]  + y_tau * Jinv[0][2] );

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
                double f_xi = vm->deval ( parent->getShapeFunction ( j ), XI, p_ ) ;
                double f_eta = vm->deval ( parent->getShapeFunction ( j ), ETA, p_ ) ;
                double f_zeta = vm->deval ( parent->getShapeFunction ( j ), ZETA, p_ ) ;
                double x = displacements[j * totaldof] ;
                double y = displacements[j * totaldof + 1] ;
                double z = displacements[j * totaldof + 2] ;

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
                double f_xi = vm->deval ( parent->getEnrichmentFunction ( j ), XI, p_ ) ;
                double f_eta = vm->deval ( parent->getEnrichmentFunction ( j ), ETA, p_ ) ;
                double f_zeta = vm->deval ( parent->getEnrichmentFunction ( j ), ZETA, p_ ) ;
                double x = enrichedDisplacements[j * totaldof] ;
                double y = enrichedDisplacements[j * totaldof + 1] ;
                double z = enrichedDisplacements[j * totaldof + 2] ;

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
                parent->getInverseJacobianMatrix ( p_, (*JinvCache) ) ;
            }
            ret[0] = ( x_xi ) * (*JinvCache)[0][0] + ( x_eta ) * (*JinvCache)[0][1]  + ( x_zeta ) * (*JinvCache)[0][2];
            ret[1] = ( y_xi ) * (*JinvCache)[1][0] + ( y_eta ) * (*JinvCache)[1][1]  + ( y_zeta ) * (*JinvCache)[1][2];
            ret[2] = ( z_xi ) * (*JinvCache)[2][0] + ( z_eta ) * (*JinvCache)[2][1]  + ( z_zeta ) * (*JinvCache)[2][2];

            ret[3] = 0.5 * ( ( y_xi ) * (*JinvCache)[2][0] +
                             ( y_eta ) * (*JinvCache)[2][1] +
                             ( y_zeta ) * (*JinvCache)[2][2] +
                             ( z_xi ) * (*JinvCache)[1][0] +
                             ( z_eta ) * (*JinvCache)[1][1] +
                             ( z_zeta ) * (*JinvCache)[1][2] );

            ret[4] = 0.5 * ( ( x_xi ) * (*JinvCache)[2][0] +
                             ( x_eta ) * (*JinvCache)[2][1] +
                             ( x_zeta ) * (*JinvCache)[2][2] +
                             ( z_xi ) * (*JinvCache)[0][0] +
                             ( z_eta ) * (*JinvCache)[0][1] +
                             ( z_zeta ) * (*JinvCache)[0][2] );

            ret[5] = 0.5 * ( ( y_xi )   * (*JinvCache)[0][0] +
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
        return ;
    }
    case MECHANICAL_STRAIN_FIELD:
    {
        getField( STRAIN_FIELD, p_, ret, true, vm ) ;
        getField( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, p_, generalizedBuffer, true, vm ) ;
        for(size_t j = 1 ; j < generalizedBuffer.size() / ret.size() ; j++)
        {
            for(size_t n = 0 ; n < ret.size() ; n++)
                ret[n] -= generalizedBuffer[j*ret.size()+n] ;
        }
        
        if(getParent()->getBehaviour() && getParent()->getBehaviour()->hasInducedForces())
            ret -= getParent()->getBehaviour()->getImposedStrain(p_) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case PRINCIPAL_MECHANICAL_STRAIN_FIELD:
    {
        this->getField ( MECHANICAL_STRAIN_FIELD, p_, strainBuffer, true,vm ) ;
        ret = toPrincipal ( strainBuffer , DOUBLE_OFF_DIAGONAL_VALUES) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case GENERALIZED_VISCOELASTIC_STRAIN_FIELD:
    {
//        std::cout << p_.getT() << " " ;
        if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
        {

            Vector dx ( 0., totaldof ) ;
            Vector dy ( 0., totaldof ) ;
            Vector dt ( 0., totaldof ) ;

            double x_xi = 0;
            double x_eta = 0;
            double y_xi = 0;
            double y_eta = 0;

            for ( size_t j = 0 ; j < parent->getShapeFunctions().size(); j++ )
            {
                double f_xi = vm->deval ( parent->getShapeFunction ( j ), XI, p_ ) ;
                double f_eta = vm->deval ( parent->getShapeFunction ( j ), ETA, p_ ) ;
                double f_tau = vm->deval ( parent->getShapeFunction ( j ), TIME_VARIABLE, p_ ) ;
                for ( int i = 0 ; i < totaldof ; i++ )
                {
                    dx[i] += f_xi  * displacements[j * totaldof + i] ;
                    dy[i] += f_eta * displacements[j * totaldof + i] ;
                    dt[i] += f_tau * displacements[j * totaldof + i] ;
                }

            }

            for ( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
            {
                double f_xi = vm->deval ( parent->getEnrichmentFunction ( j ), XI, p_ ) ;
                double f_eta = vm->deval ( parent->getEnrichmentFunction ( j ), ETA, p_ ) ;
                double f_tau = vm->deval ( parent->getEnrichmentFunction ( j ), TIME_VARIABLE, p_ ) ;
                for ( int i = 0 ; i < totaldof ; i++ )
                {
                    dx[i] += f_xi * enrichedDisplacements[j * totaldof + i] ;
                    dy[i] += f_eta * enrichedDisplacements[j * totaldof + i] ;
                    dt[i] += f_tau * enrichedDisplacements[j * totaldof + i] ;
                }

            }

            if(!JinvCache || parent->isMoved())
            {
                if(JinvCache)
                    delete JinvCache ;
                JinvCache = new Matrix ( parent->spaceDimensions()+(parent->timePlanes()>1),parent->spaceDimensions()+(parent->timePlanes()>1)) ;
                parent->getInverseJacobianMatrix ( p_, (*JinvCache) ) ;
            }
            for ( int i = 0 ; i < blocks ; i++ )
            {
                x_xi  = dx[ i * realdof + 0 ] ;
                x_eta = dy[ i * realdof + 0 ] ;
                y_xi  = dx[ i * realdof + 1 ] ;
                y_eta = dy[ i * realdof + 1 ] ;
                ret[i*3+0] = ( x_xi ) * (*JinvCache)[0][0] + ( x_eta ) * (*JinvCache)[0][1] ;//+ x_tau * Jinv[0][2];
                ret[i*3+1] = ( y_xi ) * (*JinvCache)[1][0] + ( y_eta ) * (*JinvCache)[1][1] ;//+ y_tau * Jinv[1][2] ;
                ret[i*3+2] = 0.5 * ( ( x_xi ) * (*JinvCache)[1][0] + ( x_eta ) * (*JinvCache)[1][1]  + ( y_xi ) * (*JinvCache)[0][0] + ( y_eta ) * (*JinvCache)[0][1] );//+ x_tau * Jinv[1][2]  + y_tau * Jinv[0][2]);
/*                ret[i*3+0] = x_xi ;
                ret[i*3+1] = y_eta ;
                ret[i*3+2] = 0 ;

                std::cout << x_xi << "/" << x_eta << "\t" ;*/

            }

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

            Vector dx ( 0., totaldof ) ;
            Vector dy ( 0., totaldof ) ;
            Vector dz ( 0., totaldof ) ;

            for ( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
            {
                double f_xi = vm->deval ( parent->getShapeFunction ( j ), XI, p_ ) ;
                double f_eta = vm->deval ( parent->getShapeFunction ( j ), ETA, p_ ) ;
                double f_zeta = vm->deval ( parent->getShapeFunction ( j ), ZETA, p_ ) ;

                for ( int i = 0 ; i < totaldof ; i++ )
                {
                    dx[i] += f_xi   * displacements[j * totaldof + i] ;
                    dy[i] += f_eta  * displacements[j * totaldof + i] ;
                    dz[i] += f_zeta * displacements[j * totaldof + i] ;
                }

                /*					x_xi   += f_xi   * x ;
                					x_eta  += f_eta  * x ;
                					x_zeta += f_zeta * x ;
                					y_xi   += f_xi   * y ;
                					y_eta  += f_eta  * y ;
                					y_zeta += f_zeta * y ;
                					z_xi   += f_xi   * z ;
                					z_eta  += f_eta  * z ;
                					z_zeta += f_zeta * z ;*/
            }

            for ( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
            {
                double f_xi = vm->deval ( parent->getEnrichmentFunction ( j ), XI, p_ ) ;
                double f_eta = vm->deval ( parent->getEnrichmentFunction ( j ), ETA, p_ ) ;
                double f_zeta = vm->deval ( parent->getEnrichmentFunction ( j ), ZETA, p_ ) ;

                for ( int i = 0 ; i < totaldof ; i++ )
                {
                    dx[i] += f_xi   * enrichedDisplacements[j * totaldof + i] ;
                    dy[i] += f_eta  * enrichedDisplacements[j * totaldof + i] ;
                    dz[i] += f_zeta * enrichedDisplacements[j * totaldof + i] ;
                }

                /*					x_xi += f_xi * x;
                					x_eta += f_eta * x ;
                					x_zeta += f_zeta * x ;
                					y_xi += f_xi * y ;
                					y_eta += f_eta * y ;
                					y_zeta += f_zeta * y ;
                					z_xi += f_xi * z ;
                					z_eta += f_eta * z ;
                					z_zeta += f_zeta * z ;*/
            }

            if(!JinvCache || parent->isMoved())
            {
                if(JinvCache)
                    delete JinvCache ;
                JinvCache = new Matrix ( parent->spaceDimensions()+(parent->timePlanes()>1),parent->spaceDimensions()+(parent->timePlanes()>1)) ;
                parent->getInverseJacobianMatrix ( p_, (*JinvCache) ) ;
            }
            for ( int i = 0 ; i < blocks ; i++ )
            {
                x_xi   = dx[ i * realdof + 0 ] ;
                x_eta  = dy[ i * realdof + 0 ] ;
                x_zeta = dz[ i * realdof + 0 ] ;
                y_xi   = dx[ i * realdof + 1 ] ;
                y_eta  = dy[ i * realdof + 1 ] ;
                y_zeta = dz[ i * realdof + 1 ] ;
                z_xi   = dx[ i * realdof + 2 ] ;
                z_eta  = dy[ i * realdof + 2 ] ;
                z_zeta = dz[ i * realdof + 2 ] ;
                ret[i*6+0] = ( x_xi ) * (*JinvCache)[0][0] + ( x_eta ) * (*JinvCache)[0][1]  + ( x_zeta ) * (*JinvCache)[0][2];
                ret[i*6+1] = ( y_xi ) * (*JinvCache)[1][0] + ( y_eta ) * (*JinvCache)[1][1]  + ( y_zeta ) * (*JinvCache)[1][2];
                ret[i*6+2] = ( z_xi ) * (*JinvCache)[2][0] + ( z_eta ) * (*JinvCache)[2][1]  + ( z_zeta ) * (*JinvCache)[2][2];

                ret[i*6+3] = 0.5 * ( ( y_xi ) * (*JinvCache)[2][0] +
                                     ( y_eta ) * (*JinvCache)[2][1] +
                                     ( y_zeta ) * (*JinvCache)[2][2] +
                                     ( z_xi ) * (*JinvCache)[1][0] +
                                     ( z_eta ) * (*JinvCache)[1][1] +
                                     ( z_zeta ) * (*JinvCache)[1][2] );

                ret[i*6+4] = 0.5 * ( ( x_xi ) * (*JinvCache)[2][0] +
                                     ( x_eta ) * (*JinvCache)[2][1] +
                                     ( x_zeta ) * (*JinvCache)[2][2] +
                                     ( z_xi ) * (*JinvCache)[0][0] +
                                     ( z_eta ) * (*JinvCache)[0][1] +
                                     ( z_zeta ) * (*JinvCache)[0][2] );

                ret[i*6+5] = 0.5 * ( ( y_xi )   * (*JinvCache)[0][0] +
                                     ( y_eta )  * (*JinvCache)[0][1] +
                                     ( y_zeta ) * (*JinvCache)[0][2] +
                                     ( x_xi )   * (*JinvCache)[1][0] +
                                     ( x_eta )  * (*JinvCache)[1][1] +
                                     ( x_zeta ) * (*JinvCache)[1][2] );
            }

        }

        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case PRINCIPAL_STRAIN_FIELD:
    {

        getField ( STRAIN_FIELD, p_, strainBuffer, true,vm ) ;
        ret = toPrincipal ( strainBuffer  , DOUBLE_OFF_DIAGONAL_VALUES) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case GENERALIZED_VISCOELASTIC_PRINCIPAL_STRAIN_FIELD:
    {

        this->getField ( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, p_, generalizedBuffer, true,vm ) ;
        for ( int i = 0 ; i < blocks ; i++ )
        {
            strainBuffer = 0 ;
            for ( size_t j = 0 ; j < strainBuffer.size() ; j++ )
            {
                strainBuffer[j] = generalizedBuffer[ i*strainBuffer.size() + j ] ;
            }
            principalBuffer = toPrincipal ( strainBuffer  , DOUBLE_OFF_DIAGONAL_VALUES) ;
            for ( size_t j = 0 ; j < principalBuffer.size() ; j++ )
            {
                ret[ i*principalBuffer.size() + j] = principalBuffer[j] ;
            }
        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case NON_ENRICHED_STRAIN_FIELD:
    {
        if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
        {
            double x_xi = 0;
            double x_eta = 0;
            double y_xi = 0;
            double y_eta = 0;

            for ( size_t j = 0 ; j < parent->getBoundingPoints().size(); j++ )
            {
                double f_xi = vm->deval ( parent->getShapeFunction ( j ), XI, p_ ) ;
                double f_eta = vm->deval ( parent->getShapeFunction ( j ), ETA, p_ ) ;
                x_xi += f_xi * displacements[j * totaldof] ;
                x_eta += f_eta * displacements[j * totaldof] ;
                y_xi += f_xi * displacements[j * totaldof + 1] ;
                y_eta += f_eta * displacements[j * totaldof + 1] ;
            }

            if(!JinvCache || parent->isMoved())
            {
                if(JinvCache)
                    delete JinvCache ;
                JinvCache = new Matrix ( parent->spaceDimensions()+(parent->timePlanes()>1),parent->spaceDimensions()+(parent->timePlanes()>1)) ;
                parent->getInverseJacobianMatrix ( p_, (*JinvCache) ) ;
            }
            ret[0] = ( x_xi ) * (*JinvCache)[0][0] + ( x_eta ) * (*JinvCache)[0][1] ;
            ret[1] = ( y_xi ) * (*JinvCache)[1][0] + ( y_eta ) * (*JinvCache)[1][1] ;
            ret[2] = 0.5 * ( ( x_xi ) * (*JinvCache)[1][0] + ( x_eta ) * (*JinvCache)[1][1]  + ( y_xi ) * (*JinvCache)[0][0] + ( y_eta ) * (*JinvCache)[0][1] );
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
                double f_xi = vm->deval ( parent->getShapeFunction ( j ), XI, p_ ) ;
                double f_eta = vm->deval ( parent->getShapeFunction ( j ), ETA, p_ ) ;
                double f_zeta = vm->deval ( parent->getShapeFunction ( j ), ZETA, p_ ) ;
                double x = displacements[j * totaldof] ;
                double y = displacements[j * totaldof + 1] ;
                double z = displacements[j * totaldof + 2] ;

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
                parent->getInverseJacobianMatrix ( p_, (*JinvCache) ) ;
            }
            ret[0] = ( x_xi ) * (*JinvCache)[0][0] + ( x_eta ) * (*JinvCache)[0][1]  + ( x_zeta ) * (*JinvCache)[0][2];
            ret[1] = ( y_xi ) * (*JinvCache)[1][0] + ( y_eta ) * (*JinvCache)[1][1]  + ( y_zeta ) * (*JinvCache)[1][2];
            ret[2] = ( z_xi ) * (*JinvCache)[2][0] + ( z_eta ) * (*JinvCache)[2][1]  + ( z_zeta ) * (*JinvCache)[2][2];

            ret[3] = 0.5 * ( ( y_xi ) * (*JinvCache)[2][0] +
                             ( y_eta ) * (*JinvCache)[2][1] +
                             ( y_zeta ) * (*JinvCache)[2][2] +
                             ( z_xi ) * (*JinvCache)[1][0] +
                             ( z_eta ) * (*JinvCache)[1][1] +
                             ( z_zeta ) * (*JinvCache)[1][2] );

            ret[4] = 0.5 * ( ( x_xi ) * (*JinvCache)[2][0] +
                             ( x_eta ) * (*JinvCache)[2][1] +
                             ( x_zeta ) * (*JinvCache)[2][2] +
                             ( z_xi ) * (*JinvCache)[0][0] +
                             ( z_eta ) * (*JinvCache)[0][1] +
                             ( z_zeta ) * (*JinvCache)[0][2] );

            ret[5] = 0.5 * ( ( y_xi )   * (*JinvCache)[0][0] +
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
        return ;
    }
    case GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_FIELD:
    {
        if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
        {
            Vector dx ( 0., totaldof ) ;
            Vector dy ( 0., totaldof ) ;

            double x_xi = 0;
            double x_eta = 0;
            double y_xi = 0;
            double y_eta = 0;

            for ( size_t j = 0 ; j < parent->getBoundingPoints().size(); j++ )
            {
                double f_xi = vm->deval ( parent->getShapeFunction ( j ), XI, p_ ) ;
                double f_eta = vm->deval ( parent->getShapeFunction ( j ), ETA, p_ ) ;
                for ( int i = 0 ; i < totaldof ; i++ )
                {
                    dx[i] += f_xi  * displacements[j * totaldof + i] ;
                    dy[i] += f_eta * displacements[j * totaldof + i] ;
                }

            }

            if(!JinvCache || parent->isMoved())
            {
                if(JinvCache)
                    delete JinvCache ;
                JinvCache = new Matrix ( parent->spaceDimensions()+(parent->timePlanes()>1),parent->spaceDimensions()+(parent->timePlanes()>1)) ;
                parent->getInverseJacobianMatrix ( p_, (*JinvCache) ) ;
            }
            for ( int i = 0 ; i < blocks ; i++ )
            {
                x_xi  = dx[ i * realdof + 0 ] ;
                x_eta = dy[ i * realdof + 0 ] ;
                y_xi  = dx[ i * realdof + 1 ] ;
                y_eta = dy[ i * realdof + 1 ] ;
                ret[i*3+0] = ( x_xi ) * (*JinvCache)[0][0] + ( x_eta ) * (*JinvCache)[0][1] ;
                ret[i*3+1] = ( y_xi ) * (*JinvCache)[1][0] + ( y_eta ) * (*JinvCache)[1][1] ;
                ret[i*3+2] = 0.5 * ( ( x_xi ) * (*JinvCache)[1][0] + ( x_eta ) * (*JinvCache)[1][1]  + ( y_xi ) * (*JinvCache)[0][0] + ( y_eta ) * (*JinvCache)[0][1] );
            }

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

            Vector dx ( 0., totaldof ) ;
            Vector dy ( 0., totaldof ) ;
            Vector dz ( 0., totaldof ) ;

            for ( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
            {
                double f_xi = vm->deval ( parent->getShapeFunction ( j ), XI, p_ ) ;
                double f_eta = vm->deval ( parent->getShapeFunction ( j ), ETA, p_ ) ;
                double f_zeta = vm->deval ( parent->getShapeFunction ( j ), ZETA, p_ ) ;

                for ( int i = 0 ; i < totaldof ; i++ )
                {
                    dx[i] += f_xi   * displacements[j * totaldof + i] ;
                    dy[i] += f_eta  * displacements[j * totaldof + i] ;
                    dz[i] += f_zeta * displacements[j * totaldof + i] ;
                }
            }


            if(!JinvCache || parent->isMoved())
            {
                if(JinvCache)
                    delete JinvCache ;
                JinvCache = new Matrix ( parent->spaceDimensions()+(parent->timePlanes()>1),parent->spaceDimensions()+(parent->timePlanes()>1)) ;
                parent->getInverseJacobianMatrix ( p_, (*JinvCache) ) ;
            }
            for ( int i = 0 ; i < blocks ; i++ )
            {
                x_xi   = dx[ i * realdof + 0 ] ;
                x_eta  = dy[ i * realdof + 0 ] ;
                x_zeta = dz[ i * realdof + 0 ] ;
                y_xi   = dx[ i * realdof + 1 ] ;
                y_eta  = dy[ i * realdof + 1 ] ;
                y_zeta = dz[ i * realdof + 1 ] ;
                z_xi   = dx[ i * realdof + 2 ] ;
                z_eta  = dy[ i * realdof + 2 ] ;
                z_zeta = dz[ i * realdof + 2 ] ;
                ret[i*6+0] = ( x_xi ) * (*JinvCache)[0][0] + ( x_eta ) * (*JinvCache)[0][1]  + ( x_zeta ) * (*JinvCache)[0][2];
                ret[i*6+1] = ( y_xi ) * (*JinvCache)[1][0] + ( y_eta ) * (*JinvCache)[1][1]  + ( y_zeta ) * (*JinvCache)[1][2];
                ret[i*6+2] = ( z_xi ) * (*JinvCache)[2][0] + ( z_eta ) * (*JinvCache)[2][1]  + ( z_zeta ) * (*JinvCache)[2][2];

                ret[i*6+3] = 0.5 * ( ( y_xi ) * (*JinvCache)[2][0] +
                                     ( y_eta ) * (*JinvCache)[2][1] +
                                     ( y_zeta ) * (*JinvCache)[2][2] +
                                     ( z_xi ) * (*JinvCache)[1][0] +
                                     ( z_eta ) * (*JinvCache)[1][1] +
                                     ( z_zeta ) * (*JinvCache)[1][2] );

                ret[i*6+4] = 0.5 * ( ( x_xi ) * (*JinvCache)[2][0] +
                                     ( x_eta ) * (*JinvCache)[2][1] +
                                     ( x_zeta ) * (*JinvCache)[2][2] +
                                     ( z_xi ) * (*JinvCache)[0][0] +
                                     ( z_eta ) * (*JinvCache)[0][1] +
                                     ( z_zeta ) * (*JinvCache)[0][2] );

                ret[i*6+5] = 0.5 * ( ( y_xi )   * (*JinvCache)[0][0] +
                                     ( y_eta )  * (*JinvCache)[0][1] +
                                     ( y_zeta ) * (*JinvCache)[0][2] +
                                     ( x_xi )   * (*JinvCache)[1][0] +
                                     ( x_eta )  * (*JinvCache)[1][1] +
                                     ( x_zeta ) * (*JinvCache)[1][2] );
            }

        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case VON_MISES_STRAIN_FIELD:
    {
        principalBuffer = 0 ;
        this->getField ( PRINCIPAL_STRAIN_FIELD, p_, principalBuffer, true ,vm ) ;
        if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
        {
            ret[0] = ( 2. / 3. * ( principalBuffer[0] * principalBuffer[0] + principalBuffer[1] * principalBuffer[1] ) ) ;
        }
        else if ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
        {
            ret[0] = sqrt ( 2. / 3. * ( principalBuffer[0] * principalBuffer[0] + principalBuffer[1] * principalBuffer[1] + principalBuffer[2] * principalBuffer[2] ) ) ;
        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case STRAIN_RATE_FIELD:
    {
        if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
        {

            Vector dx ( 0., realdof ) ;
            Vector dy ( 0., realdof ) ;
            Vector dt ( 0., realdof ) ;

            double x_xi = 0;
            double x_eta = 0;
            double y_xi = 0;
            double y_eta = 0;

            for ( size_t j = 0 ; j < parent->getBoundingPoints().size(); j++ )
            {
                double f_xi = vm->ddeval ( parent->getShapeFunction ( j ), XI, TIME_VARIABLE, p_ , 1e-12 ) ;
                double f_eta = vm->ddeval ( parent->getShapeFunction ( j ), ETA, TIME_VARIABLE,p_ , 1e-12 ) ;
                double f_tau = vm->ddeval ( parent->getShapeFunction ( j ), TIME_VARIABLE, TIME_VARIABLE,p_ , 1e-12 ) ;
                for ( int i = 0 ; i < realdof ; i++ )
                {
                    dx[i] += f_xi  * displacements[j * totaldof + i] ;
                    dy[i] += f_eta * displacements[j * totaldof + i] ;
                    dt[i] += f_tau * displacements[j * totaldof + i] ;
                }

            }

            for ( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
            {
                double f_xi = vm->ddeval ( parent->getEnrichmentFunction ( j ), XI, TIME_VARIABLE,p_ , 1e-12 ) ;
                double f_eta = vm->ddeval ( parent->getEnrichmentFunction ( j ), ETA, TIME_VARIABLE,p_ , 1e-12 ) ;
                double f_tau = vm->ddeval ( parent->getEnrichmentFunction ( j ), TIME_VARIABLE,TIME_VARIABLE, p_ , 1e-12 ) ;
                for ( int i = 0 ; i < realdof ; i++ )
                {
                    dx[i] += f_xi  * enrichedDisplacements[j * totaldof + i] ;
                    dy[i] += f_eta * enrichedDisplacements[j * totaldof + i] ;
                    dt[i] += f_tau * enrichedDisplacements[j * totaldof + i] ;
                }

            }

            if(!JinvCache || parent->isMoved())
            {
                if(JinvCache)
                    delete JinvCache ;
                JinvCache = new Matrix ( parent->spaceDimensions()+(parent->timePlanes()>1),parent->spaceDimensions()+(parent->timePlanes()>1)) ;
                parent->getInverseJacobianMatrix ( p_, (*JinvCache) ) ;
            }
            x_xi  = dx[ 0 ] ;
            x_eta = dy[ 0 ] ;
            y_xi  = dx[ 1 ] ;
            y_eta = dy[ 1 ] ;
            ret[0] = ( x_xi ) * (*JinvCache)[0][0] + ( x_eta ) * (*JinvCache)[0][1] ;//+ x_tau * Jinv[0][2];
            ret[1] = ( y_xi ) * (*JinvCache)[1][0] + ( y_eta ) * (*JinvCache)[1][1] ;//+ y_tau * Jinv[1][2] ;
            ret[2] = 0.5 * ( ( x_xi ) * (*JinvCache)[1][0] + ( x_eta ) * (*JinvCache)[1][1]  + ( y_xi ) * (*JinvCache)[0][0] + ( y_eta ) * (*JinvCache)[0][1] );//+ x_tau * Jinv[1][2]  + y_tau * Jinv[0][2]);
            ret *= (*JinvCache)[2][2] ;

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
                double f_xi = vm->ddeval ( parent->getShapeFunction ( j ), XI, TIME_VARIABLE, p_ , 1e-5 ) ;
                double f_eta = vm->ddeval ( parent->getShapeFunction ( j ), ETA, TIME_VARIABLE, p_ , 1e-5 ) ;
                double f_zeta = vm->ddeval ( parent->getShapeFunction ( j ), ZETA, TIME_VARIABLE, p_ , 1e-5 ) ;
                double x = displacements[j * totaldof] ;
                double y = displacements[j * totaldof + 1] ;
                double z = displacements[j * totaldof + 2] ;

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
                double f_xi = vm->ddeval ( parent->getEnrichmentFunction ( j ), XI, TIME_VARIABLE, p_ , 1e-5 ) ;
                double f_eta = vm->ddeval ( parent->getEnrichmentFunction ( j ), ETA, TIME_VARIABLE, p_ , 1e-5 ) ;
                double f_zeta = vm->ddeval ( parent->getEnrichmentFunction ( j ), ZETA, TIME_VARIABLE, p_ , 1e-5 ) ;
                double x = enrichedDisplacements[j * totaldof] ;
                double y = enrichedDisplacements[j * totaldof + 1] ;
                double z = enrichedDisplacements[j * totaldof + 2] ;

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
                parent->getInverseJacobianMatrix ( p_, (*JinvCache) ) ;
            }
            ret[0] = ( x_xi ) * (*JinvCache)[0][0] + ( x_eta ) * (*JinvCache)[0][1]  + ( x_zeta ) * (*JinvCache)[0][2];
            ret[1] = ( y_xi ) * (*JinvCache)[1][0] + ( y_eta ) * (*JinvCache)[1][1]  + ( y_zeta ) * (*JinvCache)[1][2];
            ret[2] = ( z_xi ) * (*JinvCache)[2][0] + ( z_eta ) * (*JinvCache)[2][1]  + ( z_zeta ) * (*JinvCache)[2][2];

            ret[3] = 0.5 * ( ( y_xi ) * (*JinvCache)[2][0] +
                             ( y_eta ) * (*JinvCache)[2][1] +
                             ( y_zeta ) * (*JinvCache)[2][2] +
                             ( z_xi ) * (*JinvCache)[1][0] +
                             ( z_eta ) * (*JinvCache)[1][1] +
                             ( z_zeta ) * (*JinvCache)[1][2] );

            ret[4] = 0.5 * ( ( x_xi ) * (*JinvCache)[2][0] +
                             ( x_eta ) * (*JinvCache)[2][1] +
                             ( x_zeta ) * (*JinvCache)[2][2] +
                             ( z_xi ) * (*JinvCache)[0][0] +
                             ( z_eta ) * (*JinvCache)[0][1] +
                             ( z_zeta ) * (*JinvCache)[0][2] );

            ret[5] = 0.5 * ( ( y_xi )   * (*JinvCache)[0][0] +
                             ( y_eta )  * (*JinvCache)[0][1] +
                             ( y_zeta ) * (*JinvCache)[0][2] +
                             ( x_xi )   * (*JinvCache)[1][0] +
                             ( x_eta )  * (*JinvCache)[1][1] +
                             ( x_zeta ) * (*JinvCache)[1][2] );

            ret *= (*JinvCache)[3][3] ;
        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD:
    {
        if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
        {

            Vector dx ( 0., totaldof ) ;
            Vector dy ( 0., totaldof ) ;
            Vector dt ( 0., totaldof ) ;

            double x_xi = 0;
            double x_eta = 0;
            double y_xi = 0;
            double y_eta = 0;

            for ( size_t j = 0 ; j < parent->getBoundingPoints().size(); j++ )
            {
                double f_xi = vm->ddeval ( parent->getShapeFunction ( j ), XI, TIME_VARIABLE, p_ , 1e-12 ) ;
                double f_eta = vm->ddeval ( parent->getShapeFunction ( j ), ETA, TIME_VARIABLE,p_ , 1e-12 ) ;
                double f_tau = vm->ddeval ( parent->getShapeFunction ( j ), TIME_VARIABLE, TIME_VARIABLE,p_ , 1e-12 ) ;
                for ( int i = 0 ; i < totaldof ; i++ )
                {
                    dx[i] += f_xi  * displacements[j * totaldof + i] ;
                    dy[i] += f_eta * displacements[j * totaldof + i] ;
                    dt[i] += f_tau * displacements[j * totaldof + i] ;
                }

            }

            for ( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
            {
                double f_xi = vm->ddeval ( parent->getEnrichmentFunction ( j ), XI, TIME_VARIABLE,p_ , 1e-12 ) ;
                double f_eta = vm->ddeval ( parent->getEnrichmentFunction ( j ), ETA, TIME_VARIABLE,p_ , 1e-12 ) ;
                double f_tau = vm->ddeval ( parent->getEnrichmentFunction ( j ), TIME_VARIABLE,TIME_VARIABLE, p_ , 1e-12 ) ;
                for ( int i = 0 ; i < totaldof ; i++ )
                {
                    dx[i] += f_xi  * enrichedDisplacements[j * totaldof + i] ;
                    dy[i] += f_eta * enrichedDisplacements[j * totaldof + i] ;
                    dt[i] += f_tau * enrichedDisplacements[j * totaldof + i] ;
                }

            }

            if(!JinvCache || parent->isMoved())
            {
                if(JinvCache)
                    delete JinvCache ;
                JinvCache = new Matrix ( parent->spaceDimensions()+(parent->timePlanes()>1),parent->spaceDimensions()+(parent->timePlanes()>1)) ;
                parent->getInverseJacobianMatrix ( p_, (*JinvCache) ) ;
            }
            for ( int i = 0 ; i < blocks ; i++ )
            {
                x_xi  = dx[ i * realdof + 0 ] ;
                x_eta = dy[ i * realdof + 0 ] ;
                y_xi  = dx[ i * realdof + 1 ] ;
                y_eta = dy[ i * realdof + 1 ] ;
                ret[i*3+0] = ( x_xi ) * (*JinvCache)[0][0] + ( x_eta ) * (*JinvCache)[0][1] ;//+ x_tau * Jinv[0][2];
                ret[i*3+1] = ( y_xi ) * (*JinvCache)[1][0] + ( y_eta ) * (*JinvCache)[1][1] ;//+ y_tau * Jinv[1][2] ;
                ret[i*3+2] = 0.5 * ( ( x_xi ) * (*JinvCache)[1][0] + ( x_eta ) * (*JinvCache)[1][1]  + ( y_xi ) * (*JinvCache)[0][0] + ( y_eta ) * (*JinvCache)[0][1] );
            }
            ret *= (*JinvCache)[2][2] ;

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

            Vector dx ( 0., totaldof ) ;
            Vector dy ( 0., totaldof ) ;
            Vector dz ( 0., totaldof ) ;

            for ( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
            {
                double f_xi = vm->ddeval ( parent->getShapeFunction ( j ), XI, TIME_VARIABLE, p_ , 1e-5 ) ;
                double f_eta = vm->ddeval ( parent->getShapeFunction ( j ), ETA, TIME_VARIABLE, p_ , 1e-5 ) ;
                double f_zeta = vm->ddeval ( parent->getShapeFunction ( j ), ZETA, TIME_VARIABLE,  p_, 1e-5 ) ;

                for ( int i = 0 ; i < totaldof ; i++ )
                {
                    dx[i] += f_xi   * displacements[j * totaldof + i] ;
                    dy[i] += f_eta  * displacements[j * totaldof + i] ;
                    dz[i] += f_zeta * displacements[j * totaldof + i] ;
                }

            }

            for ( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
            {
                double f_xi = vm->ddeval ( parent->getEnrichmentFunction ( j ), XI, TIME_VARIABLE, p_, 1e-5 ) ;
                double f_eta = vm->ddeval ( parent->getEnrichmentFunction ( j ), ETA, TIME_VARIABLE, p_ , 1e-5 ) ;
                double f_zeta = vm->ddeval ( parent->getEnrichmentFunction ( j ), ZETA, TIME_VARIABLE, p_, 1e-5 ) ;

                for ( int i = 0 ; i < totaldof ; i++ )
                {
                    dx[i] += f_xi   * enrichedDisplacements[j * totaldof + i] ;
                    dy[i] += f_eta  * enrichedDisplacements[j * totaldof + i] ;
                    dz[i] += f_zeta * enrichedDisplacements[j * totaldof + i] ;
                }

            }

            if(!JinvCache || parent->isMoved())
            {
                if(JinvCache)
                    delete JinvCache ;
                JinvCache = new Matrix ( parent->spaceDimensions()+(parent->timePlanes()>1),parent->spaceDimensions()+(parent->timePlanes()>1)) ;
                parent->getInverseJacobianMatrix ( p_, (*JinvCache) ) ;
            }

            for ( int i = 0 ; i < blocks ; i++ )
            {
                x_xi   = dx[ i * realdof + 0 ] ;
                x_eta  = dy[ i * realdof + 0 ] ;
                x_zeta = dz[ i * realdof + 0 ] ;
                y_xi   = dx[ i * realdof + 1 ] ;
                y_eta  = dy[ i * realdof + 1 ] ;
                y_zeta = dz[ i * realdof + 1 ] ;
                z_xi   = dx[ i * realdof + 2 ] ;
                z_eta  = dy[ i * realdof + 2 ] ;
                z_zeta = dz[ i * realdof + 2 ] ;
                ret[i*6+0] = ( x_xi ) * (*JinvCache)[0][0] + ( x_eta ) * (*JinvCache)[0][1]  + ( x_zeta ) * (*JinvCache)[0][2];
                ret[i*6+1] = ( y_xi ) * (*JinvCache)[1][0] + ( y_eta ) * (*JinvCache)[1][1]  + ( y_zeta ) * (*JinvCache)[1][2];
                ret[i*6+2] = ( z_xi ) * (*JinvCache)[2][0] + ( z_eta ) * (*JinvCache)[2][1]  + ( z_zeta ) * (*JinvCache)[2][2];

                ret[i*6+3] = 0.5 * ( ( y_xi ) * (*JinvCache)[2][0] +
                                     ( y_eta ) * (*JinvCache)[2][1] +
                                     ( y_zeta ) * (*JinvCache)[2][2] +
                                     ( z_xi ) * (*JinvCache)[1][0] +
                                     ( z_eta ) * (*JinvCache)[1][1] +
                                     ( z_zeta ) * (*JinvCache)[1][2] );

                ret[i*6+4] = 0.5 * ( ( x_xi ) * (*JinvCache)[2][0] +
                                     ( x_eta ) * (*JinvCache)[2][1] +
                                     ( x_zeta ) * (*JinvCache)[2][2] +
                                     ( z_xi ) * (*JinvCache)[0][0] +
                                     ( z_eta ) * (*JinvCache)[0][1] +
                                     ( z_zeta ) * (*JinvCache)[0][2] );

                ret[i*6+5] = 0.5 * ( ( y_xi )   * (*JinvCache)[0][0] +
                                     ( y_eta )  * (*JinvCache)[0][1] +
                                     ( y_zeta ) * (*JinvCache)[0][2] +
                                     ( x_xi )   * (*JinvCache)[1][0] +
                                     ( x_eta )  * (*JinvCache)[1][1] +
                                     ( x_zeta ) * (*JinvCache)[1][2] );
            }
            ret *= (*JinvCache)[3][3] ;

        }

        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case NON_ENRICHED_STRAIN_RATE_FIELD:
    {
        if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
        {
            double x_xi = 0;
            double x_eta = 0;
            double y_xi = 0;
            double y_eta = 0;

            for ( size_t j = 0 ; j < parent->getBoundingPoints().size(); j++ )
            {
                double f_xi = vm->ddeval ( parent->getShapeFunction ( j ), XI, TIME_VARIABLE, p_ , 1e-5 ) ;
                double f_eta = vm->ddeval ( parent->getShapeFunction ( j ), ETA, TIME_VARIABLE, p_ , 1e-5 ) ;
                x_xi += f_xi * displacements[j * totaldof] ;
                x_eta += f_eta * displacements[j * totaldof] ;
                y_xi += f_xi * displacements[j * totaldof + 1] ;
                y_eta += f_eta * displacements[j * totaldof + 1] ;
            }

            if(!JinvCache || parent->isMoved())
            {
                if(JinvCache)
                    delete JinvCache ;
                JinvCache = new Matrix ( parent->spaceDimensions()+(parent->timePlanes()>1),parent->spaceDimensions()+(parent->timePlanes()>1)) ;
                parent->getInverseJacobianMatrix ( p_, (*JinvCache) ) ;
            }
            
            ret[0] = ( x_xi ) * (*JinvCache)[0][0] + ( x_eta ) * (*JinvCache)[0][1] ;
            ret[1] = ( y_xi ) * (*JinvCache)[1][0] + ( y_eta ) * (*JinvCache)[1][1] ;
            ret[2] = 0.5 * ( ( x_xi ) * (*JinvCache)[1][0] + ( x_eta ) * (*JinvCache)[1][1]  + ( y_xi ) * (*JinvCache)[0][0] + ( y_eta ) * (*JinvCache)[0][1] );

            ret *= (*JinvCache)[2][2] ;
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
                double f_xi = vm->ddeval ( parent->getShapeFunction ( j ), XI, TIME_VARIABLE, p_ , 1e-5 ) ;
                double f_eta = vm->ddeval ( parent->getShapeFunction ( j ), ETA, TIME_VARIABLE, p_ , 1e-5 ) ;
                double f_zeta = vm->ddeval ( parent->getShapeFunction ( j ), ZETA, TIME_VARIABLE, p_ , 1e-5 ) ;
                double x = displacements[j * totaldof] ;
                double y = displacements[j * totaldof + 1] ;
                double z = displacements[j * totaldof + 2] ;

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
                parent->getInverseJacobianMatrix ( p_, (*JinvCache) ) ;
            }
            ret[0] = ( x_xi ) * (*JinvCache)[0][0] + ( x_eta ) * (*JinvCache)[0][1]  + ( x_zeta ) * (*JinvCache)[0][2];
            ret[1] = ( y_xi ) * (*JinvCache)[1][0] + ( y_eta ) * (*JinvCache)[1][1]  + ( y_zeta ) * (*JinvCache)[1][2];
            ret[2] = ( z_xi ) * (*JinvCache)[2][0] + ( z_eta ) * (*JinvCache)[2][1]  + ( z_zeta ) * (*JinvCache)[2][2];

            ret[3] = 0.5 * ( ( y_xi ) * (*JinvCache)[2][0] +
                             ( y_eta ) * (*JinvCache)[2][1] +
                             ( y_zeta ) * (*JinvCache)[2][2] +
                             ( z_xi ) * (*JinvCache)[1][0] +
                             ( z_eta ) * (*JinvCache)[1][1] +
                             ( z_zeta ) * (*JinvCache)[1][2] );

            ret[4] = 0.5 * ( ( x_xi ) * (*JinvCache)[2][0] +
                             ( x_eta ) * (*JinvCache)[2][1] +
                             ( x_zeta ) * (*JinvCache)[2][2] +
                             ( z_xi ) * (*JinvCache)[0][0] +
                             ( z_eta ) * (*JinvCache)[0][1] +
                             ( z_zeta ) * (*JinvCache)[0][2] );

            ret[5] = 0.5 * ( ( y_xi )   * (*JinvCache)[0][0] +
                             ( y_eta )  * (*JinvCache)[0][1] +
                             ( y_zeta ) * (*JinvCache)[0][2] +
                             ( x_xi )   * (*JinvCache)[1][0] +
                             ( x_eta )  * (*JinvCache)[1][1] +
                             ( x_zeta ) * (*JinvCache)[1][2] );

            ret *= (*JinvCache)[3][3] ;
        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_RATE_FIELD:
    {
        if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
        {
            Vector dx ( 0., totaldof ) ;
            Vector dy ( 0., totaldof ) ;

            double x_xi = 0;
            double x_eta = 0;
            double y_xi = 0;
            double y_eta = 0;

            for ( size_t j = 0 ; j < parent->getBoundingPoints().size(); j++ )
            {
                double f_xi = vm->ddeval ( parent->getShapeFunction ( j ), XI, TIME_VARIABLE, p_ , 1e-5 ) ;
                double f_eta = vm->ddeval ( parent->getShapeFunction ( j ), ETA, TIME_VARIABLE, p_ , 1e-5 ) ;
                for ( int i = 0 ; i < totaldof ; i++ )
                {
                    dx[i] += f_xi  * displacements[j * totaldof + i] ;
                    dy[i] += f_eta * displacements[j * totaldof + i] ;
                }

                /*					x_xi += f_xi * displacements[j * totaldof] ;
                					x_eta += f_eta * displacements[j * totaldof] ;
                					y_xi += f_xi * displacements[j * totaldof + 1] ;
                					y_eta += f_eta * displacements[j * totaldof + 1] ;*/
            }

            if(!JinvCache || parent->isMoved())
            {
                if(JinvCache)
                    delete JinvCache ;
                JinvCache = new Matrix ( parent->spaceDimensions()+(parent->timePlanes()>1),parent->spaceDimensions()+(parent->timePlanes()>1)) ;
                parent->getInverseJacobianMatrix ( p_, (*JinvCache) ) ;
            }
            for ( int i = 0 ; i < blocks ; i++ )
            {
                x_xi  = dx[ i * realdof + 0 ] ;
                x_eta = dy[ i * realdof + 0 ] ;
                y_xi  = dx[ i * realdof + 1 ] ;
                y_eta = dy[ i * realdof + 1 ] ;
                ret[i*3+0] = ( x_xi ) * (*JinvCache)[0][0] + ( x_eta ) * (*JinvCache)[0][1] ;
                ret[i*3+1] = ( y_xi ) * (*JinvCache)[1][0] + ( y_eta ) * (*JinvCache)[1][1] ;
                ret[i*3+2] = 0.5 * ( ( x_xi ) * (*JinvCache)[1][0] + ( x_eta ) * (*JinvCache)[1][1]  + ( y_xi ) * (*JinvCache)[0][0] + ( y_eta ) * (*JinvCache)[0][1] );
            }
            ret *= (*JinvCache)[2][2] ;

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

            Vector dx ( 0., totaldof ) ;
            Vector dy ( 0., totaldof ) ;
            Vector dz ( 0., totaldof ) ;

            for ( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
            {
                double f_xi = vm->ddeval ( parent->getShapeFunction ( j ), XI, TIME_VARIABLE, p_ , 1e-5 ) ;
                double f_eta = vm->ddeval ( parent->getShapeFunction ( j ), ETA, TIME_VARIABLE, p_ , 1e-5 ) ;
                double f_zeta = vm->ddeval ( parent->getShapeFunction ( j ), ZETA, TIME_VARIABLE,  p_, 1e-5 ) ;

                for ( int i = 0 ; i < totaldof ; i++ )
                {
                    dx[i] += f_xi   * displacements[j * totaldof + i] ;
                    dy[i] += f_eta  * displacements[j * totaldof + i] ;
                    dz[i] += f_zeta * displacements[j * totaldof + i] ;
                }

                /*					x_xi   += f_xi   * x ;
                					x_eta  += f_eta  * x ;
                					x_zeta += f_zeta * x ;
                					y_xi   += f_xi   * y ;
                					y_eta  += f_eta  * y ;
                					y_zeta += f_zeta * y ;
                					z_xi   += f_xi   * z ;
                					z_eta  += f_eta  * z ;
                					z_zeta += f_zeta * z ;*/
            }

            if(!JinvCache || parent->isMoved())
            {
                if(JinvCache)
                    delete JinvCache ;
                JinvCache = new Matrix ( parent->spaceDimensions()+(parent->timePlanes()>1),parent->spaceDimensions()+(parent->timePlanes()>1)) ;
                parent->getInverseJacobianMatrix ( p_, (*JinvCache) ) ;
            }
            for ( int i = 0 ; i < blocks ; i++ )
            {
                x_xi   = dx[ i * realdof + 0 ] ;
                x_eta  = dy[ i * realdof + 0 ] ;
                x_zeta = dz[ i * realdof + 0 ] ;
                y_xi   = dx[ i * realdof + 1 ] ;
                y_eta  = dy[ i * realdof + 1 ] ;
                y_zeta = dz[ i * realdof + 1 ] ;
                z_xi   = dx[ i * realdof + 2 ] ;
                z_eta  = dy[ i * realdof + 2 ] ;
                z_zeta = dz[ i * realdof + 2 ] ;
                ret[i*6+0] = ( x_xi ) * (*JinvCache)[0][0] + ( x_eta ) * (*JinvCache)[0][1]  + ( x_zeta ) * (*JinvCache)[0][2];
                ret[i*6+1] = ( y_xi ) * (*JinvCache)[1][0] + ( y_eta ) * (*JinvCache)[1][1]  + ( y_zeta ) * (*JinvCache)[1][2];
                ret[i*6+2] = ( z_xi ) * (*JinvCache)[2][0] + ( z_eta ) * (*JinvCache)[2][1]  + ( z_zeta ) * (*JinvCache)[2][2];

                ret[i*6+3] = 0.5 * ( ( y_xi ) * (*JinvCache)[2][0] +
                                     ( y_eta ) * (*JinvCache)[2][1] +
                                     ( y_zeta ) * (*JinvCache)[2][2] +
                                     ( z_xi ) * (*JinvCache)[1][0] +
                                     ( z_eta ) * (*JinvCache)[1][1] +
                                     ( z_zeta ) * (*JinvCache)[1][2] );

                ret[i*6+4] = 0.5 * ( ( x_xi ) * (*JinvCache)[2][0] +
                                     ( x_eta ) * (*JinvCache)[2][1] +
                                     ( x_zeta ) * (*JinvCache)[2][2] +
                                     ( z_xi ) * (*JinvCache)[0][0] +
                                     ( z_eta ) * (*JinvCache)[0][1] +
                                     ( z_zeta ) * (*JinvCache)[0][2] );

                ret[i*6+5] = 0.5 * ( ( y_xi )   * (*JinvCache)[0][0] +
                                     ( y_eta )  * (*JinvCache)[0][1] +
                                     ( y_zeta ) * (*JinvCache)[0][2] +
                                     ( x_xi )   * (*JinvCache)[1][0] +
                                     ( x_eta )  * (*JinvCache)[1][1] +
                                     ( x_zeta ) * (*JinvCache)[1][2] );
            }
            ret *= (*JinvCache)[3][3] ;

        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case REAL_STRESS_FIELD:
    {
        getField ( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, p_, generalizedBuffer, true,vm ) ;
        generalizedBufferSecond = ( Vector ) ( parent->getBehaviour()->getTensor ( p_, parent ) * generalizedBuffer ) ;
        getField ( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, p_, generalizedBuffer, true,vm ) ;
        generalizedBufferSecond += ( Vector ) ( parent->getBehaviour()->getViscousTensor ( p_, parent ) * generalizedBuffer ) ;

        for ( size_t i = 0 ; i < ret.size() ; i++ )
        {
            ret[i] = generalizedBufferSecond[i] ;
        }

        ret -= getParent()->getBehaviour()->getImposedStress ( p_, parent ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case GENERALIZED_VISCOELASTIC_REAL_STRESS_FIELD:
    {
        getField ( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, p_, generalizedBuffer, true,vm ) ;
        ret = ( Vector ) ( parent->getBehaviour()->getTensor ( p_, parent ) * generalizedBuffer ) ;
        getField ( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, p_, generalizedBuffer, true,vm ) ;
        ret += ( Vector ) ( parent->getBehaviour()->getViscousTensor ( p_, parent ) * generalizedBuffer ) ;

        strainBuffer = parent->getBehaviour()->getImposedStress ( p_, parent ) ;
        for ( size_t i = 0 ; i < strainBuffer.size() ; i++ )
        {
            ret[i] -= strainBuffer[i] ;
            for(int j = 1 ; j < blocks ; j++)
                ret[i+j*strainBuffer.size()] += strainBuffer[i] ;
        }

        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case PRINCIPAL_REAL_STRESS_FIELD:
    {
        this->getField ( REAL_STRESS_FIELD, p_, strainBuffer, true,vm ) ;
        ret = toPrincipal ( strainBuffer  , SINGLE_OFF_DIAGONAL_VALUES) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case GENERALIZED_VISCOELASTIC_PRINCIPAL_REAL_STRESS_FIELD:
    {
        this->getField ( GENERALIZED_VISCOELASTIC_REAL_STRESS_FIELD, p_, generalizedBufferSecond, true,vm ) ;
        for ( int i = 0 ; i < blocks ; i++ )
        {
            strainBuffer = 0 ;
            for ( size_t j = 0 ; j < strainBuffer.size() ; j++ )
            {
                strainBuffer[j] = generalizedBufferSecond[ i*strainBuffer.size() + j ] ;
            }
            principalBuffer = toPrincipal ( strainBuffer  , SINGLE_OFF_DIAGONAL_VALUES ) ;
            for ( size_t j = 0 ; j < principalBuffer.size() ; j++ )
            {
                ret[ i*principalBuffer.size() + j] = principalBuffer[j] ;
            }
        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case NON_ENRICHED_REAL_STRESS_FIELD:
    {
        getField ( GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_FIELD, p_, generalizedBuffer, true,vm ) ;
        generalizedBufferSecond = ( Vector ) ( parent->getBehaviour()->getTensor ( p_, parent ) * generalizedBuffer ) ;
        getField ( GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_RATE_FIELD, p_, generalizedBuffer, true,vm ) ;
        generalizedBufferSecond += ( Vector ) ( parent->getBehaviour()->getViscousTensor ( p_, parent ) * generalizedBuffer ) ;

        for ( size_t i = 0 ; i < ret.size() ; i++ )
        {
            ret[i] = generalizedBufferSecond[i] ;
        }

        ret -= getParent()->getBehaviour()->getImposedStress ( p_, parent ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case GENERALIZED_VISCOELASTIC_NON_ENRICHED_REAL_STRESS_FIELD:
    {
        getField ( GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_FIELD, p_, generalizedBuffer, true,vm ) ;
        ret = ( Vector ) ( parent->getBehaviour()->getTensor ( p_, parent ) * generalizedBuffer ) ;
        getField ( GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_RATE_FIELD, p_, generalizedBuffer, true,vm ) ;
        ret += ( Vector ) ( parent->getBehaviour()->getViscousTensor ( p_, parent ) * generalizedBuffer ) ;

        strainBuffer = parent->getBehaviour()->getImposedStress ( p_, parent ) ;
        for ( size_t i = 0 ; i < strainBuffer.size() ; i++ )
        {
            ret[i] -= strainBuffer[i] ;
            for(int j = 1 ; j < blocks ; j++)
                ret[i+j*strainBuffer.size()] += strainBuffer[i] ;
        }

        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case VON_MISES_REAL_STRESS_FIELD:
    {
        if ( parent->getOrder() == LINEAR )
        {
            if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
            {
                if ( !vm )
                {
                    vm = new VirtualMachine() ;
                }
                Point c ( 1./3., 1./3. ) ;
                this->getField ( PRINCIPAL_REAL_STRESS_FIELD, c, principalBuffer, true,vm ) ;
                ret[0] = sqrt ( ( ( principalBuffer[0] - principalBuffer[1] ) * ( principalBuffer[0] - principalBuffer[1] ) + principalBuffer[0] * principalBuffer[0] + principalBuffer[1] * principalBuffer[1] ) / 2. ) ;
                if ( cleanup )
                {
                    delete vm ;
                }
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
//					this->getField( PRINCIPAL_REAL_STRESS_FIELD, pts, sigma, false) ;
                sigma[0] = sqrt ( ( sigma[0] - sigma[1] ) * ( sigma[0] - sigma[1] ) + ( sigma[0] - sigma[2] ) * ( sigma[0] - sigma[2] ) + ( sigma[1] - sigma[2] ) * ( sigma[1] - sigma[2] ) ) / 6. ;
                sigma[1] = sqrt ( ( sigma[6] - sigma[7] ) * ( sigma[6] - sigma[7] ) + ( sigma[6] - sigma[8] ) * ( sigma[6] - sigma[8] ) + ( sigma[7] - sigma[8] ) * ( sigma[7] - sigma[8] ) ) / 6. ;
                sigma[2] = sqrt ( ( sigma[12] - sigma[13] ) * ( sigma[12] - sigma[13] ) + ( sigma[12] - sigma[14] ) * ( sigma[12] - sigma[14] ) + ( sigma[13] - sigma[14] ) * ( sigma[13] - sigma[14] ) ) / 6. ;
                sigma[3] = sqrt ( ( sigma[18] - sigma[19] ) * ( sigma[18] - sigma[19] ) + ( sigma[18] - sigma[20] ) * ( sigma[18] - sigma[20] ) + ( sigma[19] - sigma[20] ) * ( sigma[19] - sigma[20] ) ) / 6. ;

                ret[0] = std::max ( std::max ( std::max ( sigma[0], sigma[1] ), sigma[2] ), sigma[3] ) ;

                return ;
            }
            return ;
        }
        else
        {
            if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
            {
                Vector principalStresses ( 0., parent->getBoundingPoints().size() *2 ) ;
//					this->getField(PRINCIPAL_REAL_STRESS_FIELD, parent->getBoundingPoints(), principalStresses, false) ;
                double maxS = 0 ;

                for ( size_t i = 0 ; i < principalStresses.size() / 2 ; i++ )
                {
                    maxS = std::max ( maxS,
                                      sqrt ( ( ( principalStresses[i * 2 + 0] - principalStresses[i * 2 + 1] ) * ( principalStresses[i * 2 + 0] - principalStresses[i * 2 + 1] ) + principalStresses[i * 2 + 0] * principalStresses[i * 2 + 0] + principalStresses[i * 2 + 1] * principalStresses[i * 2 + 1] ) / 2. ) ) ;
                }

                ret[0] = maxS ;
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
//					this->getField( PRINCIPAL_REAL_STRESS_FIELD, pts, sigma, false) ;
                sigma[0] = sqrt ( ( sigma[0] - sigma[1] ) * ( sigma[0] - sigma[1] ) + ( sigma[0] - sigma[2] ) * ( sigma[0] - sigma[2] ) + ( sigma[1] - sigma[2] ) * ( sigma[1] - sigma[2] ) ) / 6 ;
                sigma[1] = sqrt ( ( sigma[6] - sigma[7] ) * ( sigma[6] - sigma[7] ) + ( sigma[6] - sigma[8] ) * ( sigma[6] - sigma[8] ) + ( sigma[7] - sigma[8] ) * ( sigma[7] - sigma[8] ) ) / 6 ;
                sigma[2] = sqrt ( ( sigma[12] - sigma[13] ) * ( sigma[12] - sigma[13] ) + ( sigma[12] - sigma[14] ) * ( sigma[12] - sigma[14] ) + ( sigma[13] - sigma[14] ) * ( sigma[13] - sigma[14] ) ) / 6 ;
                sigma[3] = sqrt ( ( sigma[18] - sigma[19] ) * ( sigma[18] - sigma[19] ) + ( sigma[18] - sigma[20] ) * ( sigma[18] - sigma[20] ) + ( sigma[19] - sigma[20] ) * ( sigma[19] - sigma[20] ) ) / 6 ;

                ret[0] = std::max ( std::max ( std::max ( sigma[0], sigma[1] ), sigma[2] ), sigma[3] ) ;
                return ;
            }
        }
        return ;
    }
    case EFFECTIVE_STRESS_FIELD:
    {
        getField ( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, p_, generalizedBuffer, true,vm ) ;
        generalizedBufferSecond = ( Vector ) ( visco->param * generalizedBuffer ) ;
        getField ( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, p_, generalizedBuffer, true,vm ) ;
        generalizedBufferSecond += ( Vector ) ( visco->eta * generalizedBuffer ) ;

        for ( size_t i = 0 ; i < ret.size() ; i++ )
        {
            ret[i] = generalizedBufferSecond[i] ;
        }

        ret -= getParent()->getBehaviour()->getImposedStress ( p_, parent ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case GENERALIZED_VISCOELASTIC_EFFECTIVE_STRESS_FIELD:
    {
        getField ( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, p_, generalizedBuffer, true,vm ) ;
        ret = ( Vector ) ( visco->param * generalizedBuffer )  ;
        getField ( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, p_, generalizedBuffer, true,vm ) ;
        ret += ( Vector ) ( visco->eta * generalizedBuffer ) ;

        strainBuffer = parent->getBehaviour()->getImposedStress ( p_, parent ) ;
        for ( size_t i = 0 ; i < strainBuffer.size() ; i++ )
        {
            ret[i] -= strainBuffer[i] ;
            for(int j = 1 ; j < blocks ; j++)
                ret[i+j*strainBuffer.size()] += strainBuffer[i] ;
        }

        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case PRINCIPAL_EFFECTIVE_STRESS_FIELD:
    {
        this->getField ( EFFECTIVE_STRESS_FIELD, p_, strainBuffer, true,vm ) ;
        ret = toPrincipal ( strainBuffer , SINGLE_OFF_DIAGONAL_VALUES) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case GENERALIZED_VISCOELASTIC_PRINCIPAL_EFFECTIVE_STRESS_FIELD:
    {
        this->getField ( GENERALIZED_VISCOELASTIC_EFFECTIVE_STRESS_FIELD, p_, generalizedBufferSecond, true,vm ) ;
        for ( int i = 0 ; i < blocks ; i++ )
        {
            strainBuffer = 0 ;
            for ( size_t j = 0 ; j < strainBuffer.size() ; j++ )
            {
                strainBuffer[j] = generalizedBufferSecond[ i*strainBuffer.size() + j ] ;
            }
            principalBuffer = toPrincipal ( strainBuffer  , SINGLE_OFF_DIAGONAL_VALUES ) ;
            for ( size_t j = 0 ; j < principalBuffer.size() ; j++ )
            {
                ret[ i*principalBuffer.size() + j] = principalBuffer[j] ;
            }
        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case NON_ENRICHED_EFFECTIVE_STRESS_FIELD:
    {
        getField ( GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_FIELD, p_, generalizedBuffer, true,vm ) ;
        generalizedBufferSecond = ( Vector ) ( visco->param * generalizedBuffer ) ;
        getField ( GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_RATE_FIELD, p_, generalizedBuffer, true,vm ) ;
        generalizedBufferSecond += ( Vector ) ( visco->eta * generalizedBuffer ) ;
 
        for ( size_t i = 0 ; i < ret.size() ; i++ )
        {
            ret[i] = generalizedBufferSecond[i] ;
        }

        ret -= getParent()->getBehaviour()->getImposedStress ( p_, parent ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case GENERALIZED_VISCOELASTIC_NON_ENRICHED_EFFECTIVE_STRESS_FIELD:
    {
        getField ( GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_FIELD, p_, generalizedBuffer, true,vm ) ;
        ret = ( Vector ) ( visco->param * generalizedBuffer ) ;
        getField ( GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_RATE_FIELD, p_, generalizedBuffer, true,vm ) ;
        ret += ( Vector ) ( visco->eta * generalizedBuffer ) ;

        strainBuffer = parent->getBehaviour()->getImposedStress ( p_, parent ) ;
        for ( size_t i = 0 ; i < strainBuffer.size() ; i++ )
        {
            ret[i] -= strainBuffer[i] ;
            for(int j = 1 ; j < blocks ; j++)
                ret[i+j*strainBuffer.size()] += strainBuffer[i] ;
        }

        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case VON_MISES_EFFECTIVE_STRESS_FIELD:
    {
        if ( parent->getOrder() == LINEAR )
        {
            if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
            {

                Point c ( 1./3., 1./3. ) ;
                this->getField ( PRINCIPAL_EFFECTIVE_STRESS_FIELD, c, principalBuffer, true,vm ) ;
                ret[0] = sqrt ( ( ( principalBuffer[0] - principalBuffer[1] ) * ( principalBuffer[0] - principalBuffer[1] ) + principalBuffer[0] * principalBuffer[0] + principalBuffer[1] * principalBuffer[1] ) / 2. ) ;
                if ( cleanup )
                {
                    delete vm ;
                }
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
//					this->getField( PRINCIPAL_EFFECTIVE_STRESS_FIELD, pts, sigma, false) ;
                sigma[0] = sqrt ( ( sigma[0] - sigma[1] ) * ( sigma[0] - sigma[1] ) + ( sigma[0] - sigma[2] ) * ( sigma[0] - sigma[2] ) + ( sigma[1] - sigma[2] ) * ( sigma[1] - sigma[2] ) ) / 6. ;
                sigma[1] = sqrt ( ( sigma[6] - sigma[7] ) * ( sigma[6] - sigma[7] ) + ( sigma[6] - sigma[8] ) * ( sigma[6] - sigma[8] ) + ( sigma[7] - sigma[8] ) * ( sigma[7] - sigma[8] ) ) / 6. ;
                sigma[2] = sqrt ( ( sigma[12] - sigma[13] ) * ( sigma[12] - sigma[13] ) + ( sigma[12] - sigma[14] ) * ( sigma[12] - sigma[14] ) + ( sigma[13] - sigma[14] ) * ( sigma[13] - sigma[14] ) ) / 6. ;
                sigma[3] = sqrt ( ( sigma[18] - sigma[19] ) * ( sigma[18] - sigma[19] ) + ( sigma[18] - sigma[20] ) * ( sigma[18] - sigma[20] ) + ( sigma[19] - sigma[20] ) * ( sigma[19] - sigma[20] ) ) / 6. ;

                ret[0] = std::max ( std::max ( std::max ( sigma[0], sigma[1] ), sigma[2] ), sigma[3] ) ;

                return ;
            }
            return ;
        }
        else
        {
            if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
            {
                Vector principalStresses ( 0., parent->getBoundingPoints().size() *2 ) ;
//					this->getField(PRINCIPAL_EFFECTIVE_STRESS_FIELD, parent->getBoundingPoints(), principalStresses, false) ;
                double maxS = 0 ;

                for ( size_t i = 0 ; i < principalStresses.size() / 2 ; i++ )
                {
                    maxS = std::max ( maxS,
                                      sqrt ( ( ( principalStresses[i * 2 + 0] - principalStresses[i * 2 + 1] ) * ( principalStresses[i * 2 + 0] - principalStresses[i * 2 + 1] ) + principalStresses[i * 2 + 0] * principalStresses[i * 2 + 0] + principalStresses[i * 2 + 1] * principalStresses[i * 2 + 1] ) / 2. ) ) ;
                }

                ret[0] = maxS ;
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
//					this->getField( PRINCIPAL_EFFECTIVE_STRESS_FIELD, pts, sigma, false) ;
                sigma[0] = sqrt ( ( sigma[0] - sigma[1] ) * ( sigma[0] - sigma[1] ) + ( sigma[0] - sigma[2] ) * ( sigma[0] - sigma[2] ) + ( sigma[1] - sigma[2] ) * ( sigma[1] - sigma[2] ) ) / 6 ;
                sigma[1] = sqrt ( ( sigma[6] - sigma[7] ) * ( sigma[6] - sigma[7] ) + ( sigma[6] - sigma[8] ) * ( sigma[6] - sigma[8] ) + ( sigma[7] - sigma[8] ) * ( sigma[7] - sigma[8] ) ) / 6 ;
                sigma[2] = sqrt ( ( sigma[12] - sigma[13] ) * ( sigma[12] - sigma[13] ) + ( sigma[12] - sigma[14] ) * ( sigma[12] - sigma[14] ) + ( sigma[13] - sigma[14] ) * ( sigma[13] - sigma[14] ) ) / 6 ;
                sigma[3] = sqrt ( ( sigma[18] - sigma[19] ) * ( sigma[18] - sigma[19] ) + ( sigma[18] - sigma[20] ) * ( sigma[18] - sigma[20] ) + ( sigma[19] - sigma[20] ) * ( sigma[19] - sigma[20] ) ) / 6 ;

                ret[0] = std::max ( std::max ( std::max ( sigma[0], sigma[1] ), sigma[2] ), sigma[3] ) ;
                return ;
            }
        }
        return ;
    }
    case PRINCIPAL_STRESS_ANGLE_FIELD:
    {

        this->getField ( REAL_STRESS_FIELD,  p_, strainBuffer, true, vm ) ;
        if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
        {
            ret[0] =  0.5 * atan2 ( strainBuffer[2], strainBuffer[0] - strainBuffer[1] ) ;
        }
        else
        {
            ret[0] =  0.5 * atan2 ( strainBuffer[3] , strainBuffer[0] - strainBuffer[1] ) ;
            ret[1] =  0.5 * atan2 ( strainBuffer[4] , strainBuffer[0] - strainBuffer[2] ) ;
            ret[2] =  0.5 * atan2 ( strainBuffer[5] , strainBuffer[1] - strainBuffer[2] ) ;
        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case PRINCIPAL_STRAIN_ANGLE_FIELD:
    {

        this->getField ( STRAIN_FIELD,  p_, strainBuffer, true, vm ) ;
        if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
        {
            ret[0] =  0.5 * atan2 ( strainBuffer[2], strainBuffer[0] - strainBuffer[1] ) ;
        }
        else
        {
            ret[0] =  0.5 * atan2 ( strainBuffer[3] , strainBuffer[0] - strainBuffer[1] ) ;
            ret[1] =  0.5 * atan2 ( strainBuffer[4] , strainBuffer[0] - strainBuffer[2] ) ;
            ret[2] =  0.5 * atan2 ( strainBuffer[5] , strainBuffer[1] - strainBuffer[2] ) ;
        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case GRADIENT_FIELD:
    {
        if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 1 )
        {
            double x_xi = 0;
            double x_eta = 0;

            for( size_t j = 0 ; j < parent->getBoundingPoints().size(); j++ )
            {
                double f_xi = vm->deval( parent->getShapeFunction( j ), XI, p_ ) ;
                double f_eta = vm->deval( parent->getShapeFunction( j ), ETA, p_ ) ;
                x_xi += f_xi * displacements[j] ;
                x_eta += f_eta * displacements[j] ;
            }

            for( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() && j < enrichedDisplacements.size(); j++ )
            {
                double f_xi = vm->deval( parent->getEnrichmentFunction( j ), XI, p_ ) ;
                double f_eta = vm->deval( parent->getEnrichmentFunction( j ), ETA, p_ ) ;
                x_xi += f_xi * enrichedDisplacements[j] ;
                x_eta += f_eta * enrichedDisplacements[j] ;

            }
            if(!JinvCache || parent->isMoved())
            {
                if(JinvCache)
                    delete JinvCache ;
                JinvCache = new Matrix ( parent->spaceDimensions()+(parent->timePlanes()>1),parent->spaceDimensions()+(parent->timePlanes()>1)) ;
                parent->getInverseJacobianMatrix ( p_, (*JinvCache) ) ;
            }
            parent->getInverseJacobianMatrix( p_, (*JinvCache) ) ;
            ret[0] = ( x_xi ) * (*JinvCache)[0][0] + ( x_eta ) * (*JinvCache)[0][1] ;
            ret[1] = ( x_xi ) * (*JinvCache)[1][0] + ( x_eta ) * (*JinvCache)[1][1] ;
        }
        else if( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 1 )
        {
            double x_xi = 0;
            double x_eta = 0;
            double x_zeta = 0;

            for( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
            {
                double f_xi = vm->deval( parent->getShapeFunction( j ), XI, p_ ) ;
                double f_eta = vm->deval( parent->getShapeFunction( j ), ETA, p_ ) ;
                double f_zeta = vm->deval( parent->getShapeFunction( j ), ZETA, p_ ) ;
                double x = displacements[j] ;

                x_xi   += f_xi   * x ;
                x_eta  += f_eta  * x ;
                x_zeta += f_zeta * x ;
            }

            for( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
            {
                double f_xi = vm->deval( parent->getEnrichmentFunction( j ), XI, p_ ) ;
                double f_eta = vm->deval( parent->getEnrichmentFunction( j ), ETA, p_ ) ;
                double f_zeta = vm->deval( parent->getEnrichmentFunction( j ), ZETA, p_ ) ;
                double x = enrichedDisplacements[j] ;

                x_xi += f_xi * x;
                x_eta += f_eta * x ;
                x_zeta += f_zeta * x ;
            }

            Matrix Jinv( 3, 3 ) ;
            parent->getInverseJacobianMatrix( p_, Jinv ) ;
            ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1]  + ( x_zeta ) * Jinv[0][2];
            ret[1] = ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( x_zeta ) * Jinv[1][2];
            ret[2] = ( x_xi ) * Jinv[2][0] + ( x_eta ) * Jinv[2][1]  + ( x_zeta ) * Jinv[2][2];
        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case FLUX_FIELD:
        {
            this->getField(GRADIENT_FIELD, p_, ret, true) ;
            ret = (Vector) (parent->getBehaviour()->getTensor(p_, parent) * ret) ;
            if ( cleanup )
            {
                delete vm ;
            }
            return ;
        }
    case INTERNAL_VARIABLE_FIELD:
        {
            std::cout << "field not handled" << std::endl ;
            exit(0) ;
        }
    }
}

void GeneralizedSpaceTimeViscoElasticElementState::getFieldAtNodes ( FieldType f, Vector & ret, VirtualMachine * vm, int )
{
    int totaldof = parent->getBehaviour()->getNumberOfDegreesOfFreedom() ;
    int realdof = parent->spaceDimensions() ;
    bool cleanup = !vm ;
    switch ( f )
    {
    case DISPLACEMENT_FIELD:
        for ( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
        {
            for ( int i = 0 ; i < realdof ; i++ )
            {
                ret[j*realdof+i] = displacements[ j*totaldof + i ] ;
            }
        }
        return ;
    case ENRICHED_DISPLACEMENT_FIELD:
        for ( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
        {
            for ( int i = 0 ; i < realdof ; i++ )
            {
                ret[j*realdof+i] = enrichedDisplacements[ j*totaldof + i ] ;
            }
        }
        return ;
    case GENERALIZED_VISCOELASTIC_DISPLACEMENT_FIELD:
        ret = displacements ;
        return ;
    case GENERALIZED_VISCOELASTIC_ENRICHED_DISPLACEMENT_FIELD:
        ret = enrichedDisplacements ;
        return ;
    case STRAIN_FIELD :

        ElementState::getField ( f, parent->getBoundingPoints(), ret, false,vm ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case MECHANICAL_STRAIN_FIELD :

        ElementState::getField ( f, parent->getBoundingPoints(), ret, false,vm ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case PRINCIPAL_STRAIN_FIELD :

        ElementState::getField ( f, parent->getBoundingPoints(), ret, false,vm ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case PRINCIPAL_MECHANICAL_STRAIN_FIELD :

        ElementState::getField ( f, parent->getBoundingPoints(), ret, false,vm ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case NON_ENRICHED_STRAIN_FIELD :

        ElementState::getField ( f, parent->getBoundingPoints(), ret, false,vm ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case REAL_STRESS_FIELD:

        ElementState::getField ( f,  parent->getBoundingPoints(), ret, false,vm ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case PRINCIPAL_REAL_STRESS_FIELD :

        ElementState::getField ( f, parent->getBoundingPoints(), ret, false,vm ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case NON_ENRICHED_REAL_STRESS_FIELD:

        ElementState::getField ( f,  parent->getBoundingPoints(), ret, false,vm ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case EFFECTIVE_STRESS_FIELD:
 
        ElementState::getField ( f,  parent->getBoundingPoints(), ret, false,vm ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case PRINCIPAL_EFFECTIVE_STRESS_FIELD :

        ElementState::getField ( f, parent->getBoundingPoints(), ret, false,vm ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case NON_ENRICHED_EFFECTIVE_STRESS_FIELD:

        ElementState::getField ( f,  parent->getBoundingPoints(), ret, false,vm ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    default:

        ElementState::getField ( f,  parent->getBoundingPoints(), ret, false,vm ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    ElementState::getField ( f, parent->getBoundingPoints(), ret, false,vm ) ;
    if ( cleanup )
    {
        delete vm ;
    }
}

void GeneralizedSpaceTimeViscoElasticElementState::getField ( FieldType f1, FieldType f2, const Point & p, Vector & ret1, Vector & ret2, bool local, VirtualMachine * vm, int , int ) 
{
    bool cleanup = !vm ;
    if ( !vm )
    {
        vm = new VirtualMachine() ;
    }
    getField ( f1, p, ret1, local, vm ) ;
    getField ( f2, p, ret2, local, vm ) ;
    if ( cleanup )
    {
        delete vm ;
    }
}

void GeneralizedSpaceTimeViscoElasticElementState::getFieldAtNodes ( FieldType f1, FieldType f2, Vector & ret1, Vector & ret2, VirtualMachine * vm, int , int )
{
    bool cleanup = !vm ;
    if ( !vm )
    {
        vm = new VirtualMachine() ;
    }
    getFieldAtNodes ( f1, ret1, vm ) ;
    getFieldAtNodes ( f2, ret2, vm ) ;
    if ( cleanup )
    {
        delete vm ;
    }
}



GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables::GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables ( IntegrableEntity * e, std::map<std::string, double> & v , int blocks) : GeneralizedSpaceTimeViscoElasticElementState ( e, blocks )
{
    variables = v ;

}

GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables::GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables ( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s ) : GeneralizedSpaceTimeViscoElasticElementState ( s )
{
    variables = s.getVariables() ;
}


GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables::operator = ( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s )
{
    ElementState::operator = ( s ) ;

    return *this ;
}

bool GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables::has ( std::string v ) const
{
    return ! ( variables.empty() || variables.find ( v ) == variables.end() ) ;
}

double GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables::get ( std::string v, std::map<std::string, double> & defaultValues )
{
    if ( has ( v ) )
    {
        return variables[v] ;
    }
    return defaultValues[v] ;
}

void GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables::set ( std::string v, double d )
{
    if ( has ( v ) )
    {
        variables[v] = d ;
    }
    else
    {
        variables.insert ( std::pair<std::string, double> ( v, d ) ) ;
    }
}

void GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables::add ( std::string v, double d )
{
    if ( has ( v ) )
    {
        variables[v] += d ;
    }
    else
    {
        variables.insert ( std::pair<std::string, double> ( v, d ) ) ;
    }
}

void GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables::multiply ( std::string v, double d )
{
    if ( has ( v ) )
    {
        variables[v] *= d ;
    }

}

void GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables::synchronize(std::map<std::string, double> & defaultValues)
{
    for( std::map<std::string, double>::iterator it = defaultValues.begin() ; it != defaultValues.end() ; it++)
    {
        variables[it->first] = it->second ;
    }
}

