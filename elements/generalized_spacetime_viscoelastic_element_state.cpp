#include "generalized_spacetime_viscoelastic_element_state.h"
#include "../physics/viscoelasticity.h"
#include "../physics/dual_behaviour.h"
#include "../physics/damagemodels/damagemodel.h"
#include "../mesher/delaunay.h"
#include "../utilities/random.h"
#include "../features/inclusion.h"
#include <omp.h>

using namespace Amie ;

GeneralizedSpaceTimeViscoElasticElementState::GeneralizedSpaceTimeViscoElasticElementState ( IntegrableEntity * e ) : ElementState ( e )
{


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
        GaussPointArray original = trg->getGaussPoints() ;
        for ( size_t i = 0 ; i < original.gaussPoints.size() /3 ; i++ )
        {
            gp_alternative.push_back ( original.gaussPoints[i*3] );
            gp_alternative.back().first.getT() = time ;
            gp_alternative.back().second /= 1./3. ;
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
        GaussPointArray original = tet->getGaussPoints() ;
        for ( size_t i = 0 ; i < original.gaussPoints.size() /3 ; i++ )
        {
            gp_alternative.push_back ( original.gaussPoints[i*3] );
            gp_alternative.back().first.getT() = time ;
            gp_alternative.back().second /= 1./3. ;
        }

        gp.gaussPoints.resize ( gp_alternative.size() ) ;
        std::copy ( gp_alternative.begin(), gp_alternative.end(), &gp.gaussPoints[0] );

    }

    return gp ;
}

double GeneralizedSpaceTimeViscoElasticElementState::getAverageField ( FieldType f, Vector & ret, VirtualMachine * vm, int dummy , double t, std::vector<double> weights )
{
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
	        getField ( f, p_, tmp, true, vm, dummy ) ;

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
	}

	std::cout << "weight=" << weights.size() << std::endl ;
	std::cout << "time=" << t << std::endl ;
	std::cout << "field=" << f << std::endl ;
	std::cout << "dummy=" << dummy << std::endl ;
	std::cout << "area=" <<(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL ? parent->area() : parent->volume()) << std::endl ;
	std::cout << "total=" << total << std::endl ;*/

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
	if(f == GENERALIZED_VISCOELASTIC_STRAIN_FIELD)
	{
		// build cache if it does not exist
		if(genStrainAtGaussPointBefore.size() == 0)
		{
			genStrainAtGaussPointBefore.resize( gp.gaussPoints.size() * fieldTypeElementarySize( f, parent->spaceDimensions(), parent->getBehaviour()->getNumberOfDegreesOfFreedom()/parent->spaceDimensions())) ;
			Vector tmp(genStrainAtGaussPointBefore.size()/gp.gaussPoints.size()) ;
			for(size_t j = 0 ; j < gp.gaussPoints.size() ; j++)
			{
				tmp = 0. ;
				getField( f, gp.gaussPoints[j].first, tmp, true, vm ) ;
				for(size_t k = 0 ; k < tmp.size() ; k++)
					genStrainAtGaussPointBefore[ j*ret.size() + k] = tmp[k] ;
			}
		}
		ret.resize( genStrainAtGaussPointBefore.size()/gp.gaussPoints.size() ) ;
		// find result in cache
		for(size_t j = 0 ; j < ret.size() ; j++)
			ret[j] = genStrainAtGaussPointBefore[ i*ret.size() + j] ;
		return ret ;		
	}
	if(f == STRAIN_FIELD)
	{
		// build cache if it does not exist
		if(genStrainAtGaussPointBefore.size() == 0)
		{
			genStrainAtGaussPointBefore.resize( gp.gaussPoints.size() * fieldTypeElementarySize( f, parent->spaceDimensions(), parent->getBehaviour()->getNumberOfDegreesOfFreedom()/parent->spaceDimensions())) ;
			Vector tmp(genStrainAtGaussPointBefore.size()/gp.gaussPoints.size()) ;
			for(size_t j = 0 ; j < gp.gaussPoints.size() ; j++)
			{
				tmp = 0. ;
				getField( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, gp.gaussPoints[j].first, tmp, true, vm ) ;
				for(size_t k = 0 ; k < tmp.size() ; k++)
					genStrainAtGaussPointBefore[ j*ret.size() + k] = tmp[k] ;
			}
		}
		ret.resize( 3+(3*parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL) ) ;
		// find result in cache
		for(size_t j = 0 ; j < ret.size() ; j++)
			ret[j] = genStrainAtGaussPointBefore[ i*genStrainAtGaussPointBefore.size()/gp.gaussPoints.size() + j] ;
		return ret ;		
	}
	if(f == MECHANICAL_STRAIN_FIELD)
	{
		ret = getCachedFieldAtGaussPointBefore( STRAIN_FIELD, gp, i, vm) ;
		if( getParent()->getBehaviour() && getParent()->getBehaviour()->hasInducedForces() )
			ret -= getParent()->getBehaviour()->getImposedStrain( gp.gaussPoints[i].first ) ;
		return ret ;		
	}
	if(f == PRINCIPAL_MECHANICAL_STRAIN_FIELD)
	{
		Vector strain = getCachedFieldAtGaussPointBefore( MECHANICAL_STRAIN_FIELD, gp, i, vm) ;
		return toPrincipal(strain) ;
	}
	if(f == GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD)
	{
		// build cache if it does not exist
		if(genStrainRateAtGaussPointBefore.size() == 0)
		{
			genStrainRateAtGaussPointBefore.resize( gp.gaussPoints.size() * fieldTypeElementarySize( f, parent->spaceDimensions(), parent->getBehaviour()->getNumberOfDegreesOfFreedom()/parent->spaceDimensions())) ;
			Vector tmp(genStrainRateAtGaussPointBefore.size()/gp.gaussPoints.size()) ;
			for(size_t j = 0 ; j < gp.gaussPoints.size() ; j++)
			{
				tmp = 0. ;
				getField( f, gp.gaussPoints[j].first, tmp, true, vm ) ;
				for(size_t k = 0 ; k < tmp.size() ; k++)
					genStrainRateAtGaussPointBefore[ j*ret.size() + k] = tmp[k] ;
			}
		}
		ret.resize( genStrainRateAtGaussPointBefore.size()/gp.gaussPoints.size() ) ;
		// find result in cache
		for(size_t j = 0 ; j < ret.size() ; j++)
			ret[j] = genStrainRateAtGaussPointBefore[ i*ret.size() + j] ;
		return ret ;		
	}
	if(f == STRAIN_RATE_FIELD)
	{
		// build cache if it does not exist
		if(genStrainRateAtGaussPointBefore.size() == 0)
		{
			genStrainRateAtGaussPointBefore.resize( gp.gaussPoints.size() * fieldTypeElementarySize( f, parent->spaceDimensions(), parent->getBehaviour()->getNumberOfDegreesOfFreedom()/parent->spaceDimensions())) ;
			Vector tmp(genStrainRateAtGaussPointBefore.size()/gp.gaussPoints.size()) ;
			for(size_t j = 0 ; j < gp.gaussPoints.size() ; j++)
			{
				tmp = 0. ;
				getField( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, gp.gaussPoints[j].first, tmp, true, vm ) ;
				for(size_t k = 0 ; k < tmp.size() ; k++)
					genStrainRateAtGaussPointBefore[ j*ret.size() + k] = tmp[k] ;
			}
		}
		ret.resize( 3+(3*parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL) ) ;
		// find result in cache
		for(size_t j = 0 ; j < ret.size() ; j++)
			ret[j] = genStrainRateAtGaussPointBefore[ i*genStrainRateAtGaussPointBefore.size()/gp.gaussPoints.size() + j] ;
		return ret ;		
	}
	if(f == REAL_STRESS_FIELD)
	{
		Vector strain = getCachedFieldAtGaussPointBefore( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, gp, i, vm) ;
		Vector strainrate = getCachedFieldAtGaussPointBefore( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, gp, i, vm) ;
		Vector tmp = (Vector) ((parent->getBehaviour()->getTensor( gp.gaussPoints[i].first ))*strain) + (Vector) ((parent->getBehaviour()->getViscousTensor( gp.gaussPoints[i].first ))*strainrate) ;
		Vector imp = parent->getBehaviour()->getImposedStress( gp.gaussPoints[i].first ) ;
		ret.resize(imp.size()) ;
		for(size_t j = 0 ; j < ret.size() ; j++)
			ret[j] = tmp[j]-imp[j] ;
		return ret ;
	}
	if(f == EFFECTIVE_STRESS_FIELD)
	{
		Vector strain = getCachedFieldAtGaussPointBefore( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, gp, i, vm) ;
		Vector strainrate = getCachedFieldAtGaussPointBefore( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, gp, i, vm) ;
		Vector tmp = (Vector) ((parent->getBehaviour()->param)*strain) + (Vector) (((dynamic_cast<Viscoelasticity *>(parent->getBehaviour()))->eta)*strainrate) ;
		Vector imp = parent->getBehaviour()->getImposedStress( gp.gaussPoints[i].first ) ;
		ret.resize(imp.size()) ;
		for(size_t j = 0 ; j < ret.size() ; j++)
			ret[j] = tmp[j]-imp[j] ;
		return ret ;
	}
	if(f == PRINCIPAL_STRAIN_FIELD)
	{
		Vector strain = getCachedFieldAtGaussPointBefore( STRAIN_FIELD, gp, i, vm) ;
		return toPrincipal(strain) ;
	}
	if(f == PRINCIPAL_REAL_STRESS_FIELD)
	{
		Vector stress = getCachedFieldAtGaussPointBefore( REAL_STRESS_FIELD, gp, i, vm) ;
		return toPrincipal(stress) ;
	}
	if(f == PRINCIPAL_EFFECTIVE_STRESS_FIELD)
	{
		Vector stress = getCachedFieldAtGaussPointBefore( EFFECTIVE_STRESS_FIELD, gp, i, vm) ;
		return toPrincipal(stress) ;
	}
	return ret ;
}

Vector GeneralizedSpaceTimeViscoElasticElementState::getCachedFieldAtGaussPointAfter( FieldType f, GaussPointArray & gp, size_t i, VirtualMachine * vm) 
{
	Vector ret ;
	if(f == GENERALIZED_VISCOELASTIC_STRAIN_FIELD)
	{
		// build cache if it does not exist
		if(genStrainAtGaussPointAfter.size() == 0)
		{
			genStrainAtGaussPointAfter.resize( gp.gaussPoints.size() * fieldTypeElementarySize( f, parent->spaceDimensions(), parent->getBehaviour()->getNumberOfDegreesOfFreedom()/parent->spaceDimensions())) ;
			Vector tmp(genStrainAtGaussPointAfter.size()/gp.gaussPoints.size()) ;
			for(size_t j = 0 ; j < gp.gaussPoints.size() ; j++)
			{
				tmp = 0. ;
				getField( f, gp.gaussPoints[j].first, tmp, true, vm ) ;
				for(size_t k = 0 ; k < tmp.size() ; k++)
					genStrainAtGaussPointAfter[ j*ret.size() + k] = tmp[k] ;
			}
		}
		ret.resize( genStrainAtGaussPointAfter.size()/gp.gaussPoints.size() ) ;
		// find result in cache
		for(size_t j = 0 ; j < ret.size() ; j++)
			ret[j] = genStrainAtGaussPointAfter[ i*ret.size() + j] ;
		return ret ;		
	}
	if(f == STRAIN_FIELD)
	{
		// build cache if it does not exist
		if(genStrainAtGaussPointAfter.size() == 0)
		{
			genStrainAtGaussPointAfter.resize( gp.gaussPoints.size() * fieldTypeElementarySize( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, parent->spaceDimensions(), parent->getBehaviour()->getNumberOfDegreesOfFreedom()/parent->spaceDimensions())) ;
			Vector tmp(genStrainAtGaussPointAfter.size()/gp.gaussPoints.size()) ;
			for(size_t j = 0; j < gp.gaussPoints.size() ; j++)
			{
				tmp = 0. ;
				getField( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, gp.gaussPoints[j].first, tmp, true, vm ) ;
				for(size_t k = 0 ; k < tmp.size() ; k++)
					genStrainAtGaussPointAfter[ j*ret.size() + k] = tmp[k] ;
			}
		}
		ret.resize( 3+(3*parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL) ) ;
		// find result in cache
		for(size_t j = 0 ; j < ret.size() ; j++)
			ret[j] = genStrainAtGaussPointAfter[ i*genStrainAtGaussPointAfter.size()/gp.gaussPoints.size() + j] ;
		return ret ;		
	}
	if(f == MECHANICAL_STRAIN_FIELD)
	{
		ret = getCachedFieldAtGaussPointAfter( STRAIN_FIELD, gp, i, vm) ;
		if( getParent()->getBehaviour() && getParent()->getBehaviour()->hasInducedForces() )
			ret -= getParent()->getBehaviour()->getImposedStrain( gp.gaussPoints[i].first ) ;
		return ret ;		
	}
	if(f == PRINCIPAL_MECHANICAL_STRAIN_FIELD)
	{
		Vector strain = getCachedFieldAtGaussPointAfter( MECHANICAL_STRAIN_FIELD, gp, i, vm) ;
		return toPrincipal(strain) ;
	}
	if(f == GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD)
	{
		// build cache if it does not exist
		if(genStrainRateAtGaussPointAfter.size() == 0)
		{
			genStrainRateAtGaussPointAfter.resize( gp.gaussPoints.size() * fieldTypeElementarySize( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, parent->spaceDimensions(), parent->getBehaviour()->getNumberOfDegreesOfFreedom()/parent->spaceDimensions())) ;
			Vector tmp(genStrainRateAtGaussPointAfter.size()/gp.gaussPoints.size()) ;
			for(size_t j = 0 ; j < gp.gaussPoints.size() ; j++)
			{
				tmp = 0. ;
				getField( f, gp.gaussPoints[j].first, tmp, true, vm ) ;
				for(size_t k = 0 ; k < tmp.size() ; k++)
					genStrainRateAtGaussPointAfter[ j*ret.size() + k] = tmp[k] ;
			}
		}
		ret.resize( genStrainRateAtGaussPointAfter.size()/gp.gaussPoints.size() ) ;
		// find result in cache
		for(size_t j = 0 ; j < ret.size() ; j++)
			ret[j] = genStrainRateAtGaussPointAfter[ i*ret.size() + j] ;
		return ret ;		
	}
	if(f == STRAIN_RATE_FIELD)
	{
		// build cache if it does not exist
		if(genStrainRateAtGaussPointAfter.size() == 0)
		{
			genStrainRateAtGaussPointAfter.resize( gp.gaussPoints.size() * fieldTypeElementarySize( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, parent->spaceDimensions(), parent->getBehaviour()->getNumberOfDegreesOfFreedom()/parent->spaceDimensions())) ;
			Vector tmp(genStrainRateAtGaussPointAfter.size()/gp.gaussPoints.size()) ;
			for(size_t j = 0 ; j < gp.gaussPoints.size() ; j++)
			{
				tmp = 0. ;
				getField( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, gp.gaussPoints[j].first, tmp, true, vm ) ;
				for(size_t k = 0 ; k < tmp.size() ; k++)
					genStrainRateAtGaussPointAfter[ j*ret.size() + k] = tmp[k] ;
			}
		}
		ret.resize( 3+(3*parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL) ) ;
		// find result in cache
		for(size_t j = 0 ; j < ret.size() ; j++)
			ret[j] = genStrainRateAtGaussPointAfter[ i*genStrainRateAtGaussPointAfter.size()/gp.gaussPoints.size() + j] ;
		return ret ;		
	}
	if(f == REAL_STRESS_FIELD)
	{
		Vector strain = getCachedFieldAtGaussPointAfter( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, gp, i, vm) ;
		Vector strainrate = getCachedFieldAtGaussPointAfter( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, gp, i, vm) ;
		Vector tmp = (Vector) ((parent->getBehaviour()->getTensor( gp.gaussPoints[i].first ))*strain) + (Vector) ((parent->getBehaviour()->getViscousTensor( gp.gaussPoints[i].first ))*strainrate) ;
		Vector imp = parent->getBehaviour()->getImposedStress( gp.gaussPoints[i].first ) ;
		ret.resize(imp.size()) ;
		for(size_t j = 0 ; j < ret.size() ; j++)
			ret[j] = tmp[j]-imp[j] ;
		return ret ;
	}
	if(f == EFFECTIVE_STRESS_FIELD)
	{
		Vector strain = getCachedFieldAtGaussPointAfter( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, gp, i, vm) ;
		Vector strainrate = getCachedFieldAtGaussPointAfter( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, gp, i, vm) ;
		Vector tmp = (Vector) ((parent->getBehaviour()->param)*strain) + (Vector) (((dynamic_cast<Viscoelasticity *>(parent->getBehaviour()))->eta)*strainrate) ;
		Vector imp = parent->getBehaviour()->getImposedStress( gp.gaussPoints[i].first ) ;
		ret.resize(imp.size()) ;
		for(size_t j = 0 ; j < ret.size() ; j++)
			ret[j] = tmp[j]-imp[j] ;
		return ret ;
	}
	if(f == PRINCIPAL_STRAIN_FIELD)
	{
		Vector strain = getCachedFieldAtGaussPointAfter( STRAIN_FIELD, gp, i, vm) ;
		return toPrincipal(strain) ;
	}
	if(f == PRINCIPAL_REAL_STRESS_FIELD)
	{
		Vector stress = getCachedFieldAtGaussPointAfter( REAL_STRESS_FIELD, gp, i, vm) ;
		return toPrincipal(stress) ;
	}
	if(f == PRINCIPAL_EFFECTIVE_STRESS_FIELD)
	{
		Vector stress = getCachedFieldAtGaussPointAfter( EFFECTIVE_STRESS_FIELD, gp, i, vm) ;
		return toPrincipal(stress) ;
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

    Form * visco = parent->getBehaviour() ;

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
        Matrix Jinv ( 3,3 ) ;
        parent->getInverseJacobianMatrix ( gp.gaussPoints[i].first, Jinv ) ;

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
                tmpstrain[k*3+0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;//+ x_tau * Jinv[0][2];
                tmpstrain[k*3+1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;//+ y_tau * Jinv[1][2] ;
                tmpstrain[k*3+2] = 0.5 * ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] );//+ x_tau * Jinv[1][2]  + y_tau * Jinv[0][2]);
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
                tmpstrain[k*6+0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1]  + ( x_zeta ) * Jinv[0][2];
                tmpstrain[k*6+1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1]  + ( y_zeta ) * Jinv[1][2];
                tmpstrain[k*6+2] = ( z_xi ) * Jinv[2][0] + ( z_eta ) * Jinv[2][1]  + ( z_zeta ) * Jinv[2][2];

                tmpstrain[k*6+3] = 0.5 * ( ( y_xi ) * Jinv[2][0] +
                                           ( y_eta ) * Jinv[2][1] +
                                           ( y_zeta ) * Jinv[2][2] +
                                           ( z_xi ) * Jinv[1][0] +
                                           ( z_eta ) * Jinv[1][1] +
                                           ( z_zeta ) * Jinv[1][2] );

                tmpstrain[k*6+4] = 0.5 * ( ( x_xi ) * Jinv[2][0] +
                                           ( x_eta ) * Jinv[2][1] +
                                           ( x_zeta ) * Jinv[2][2] +
                                           ( z_xi ) * Jinv[0][0] +
                                           ( z_eta ) * Jinv[0][1] +
                                           ( z_zeta ) * Jinv[0][2] );

                tmpstrain[k*6+5] = 0.5 * ( ( y_xi )   * Jinv[0][0] +
                                           ( y_eta )  * Jinv[0][1] +
                                           ( y_zeta ) * Jinv[0][2] +
                                           ( x_xi )   * Jinv[1][0] +
                                           ( x_eta )  * Jinv[1][1] +
                                           ( x_zeta ) * Jinv[1][2] );
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
                tmpstrainrate[k*3+0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;//+ x_tau * Jinv[0][2];
                tmpstrainrate[k*3+1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;//+ y_tau * Jinv[1][2] ;
                tmpstrainrate[k*3+2] = 0.5 * ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] );//+ x_tau * Jinv[1][2]  + y_tau * Jinv[0][2]);
            }
            tmpstrainrate *= Jinv[2][2] ;

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
                tmpstrainrate[k*6+0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1]  + ( x_zeta ) * Jinv[0][2];
                tmpstrainrate[k*6+1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1]  + ( y_zeta ) * Jinv[1][2];
                tmpstrainrate[k*6+2] = ( z_xi ) * Jinv[2][0] + ( z_eta ) * Jinv[2][1]  + ( z_zeta ) * Jinv[2][2];

                tmpstrainrate[k*6+3] = 0.5 * ( ( y_xi ) * Jinv[2][0] +
                                               ( y_eta ) * Jinv[2][1] +
                                               ( y_zeta ) * Jinv[2][2] +
                                               ( z_xi ) * Jinv[1][0] +
                                               ( z_eta ) * Jinv[1][1] +
                                               ( z_zeta ) * Jinv[1][2] );

                tmpstrainrate[k*6+4] = 0.5 * ( ( x_xi ) * Jinv[2][0] +
                                               ( x_eta ) * Jinv[2][1] +
                                               ( x_zeta ) * Jinv[2][2] +
                                               ( z_xi ) * Jinv[0][0] +
                                               ( z_eta ) * Jinv[0][1] +
                                               ( z_zeta ) * Jinv[0][2] );

                tmpstrainrate[k*6+5] = 0.5 * ( ( y_xi )   * Jinv[0][0] +
                                               ( y_eta )  * Jinv[0][1] +
                                               ( y_zeta ) * Jinv[0][2] +
                                               ( x_xi )   * Jinv[1][0] +
                                               ( x_eta )  * Jinv[1][1] +
                                               ( x_zeta ) * Jinv[1][2] );
            }
            tmpstrainrate *= Jinv[3][3] ;

        }


        if ( f == REAL_STRESS_FIELD )
        {
            tmpstress = ( Vector ) ( visco->getTensor ( gp.gaussPoints[i].first, parent ) * tmpstrain )
                        + ( Vector ) ( visco->getViscousTensor ( gp.gaussPoints[i].first, parent ) * tmpstrainrate ) ;
            for ( size_t j = 0 ; j < stress.size() ; j++ )
            {
                stress[j] += tmpstress[j]*w ;
            }
        }
        else
        {
            tmpstress = ( Vector ) ( visco->param * tmpstrain )
                        + ( Vector ) ( visco->getViscousTensor ( gp.gaussPoints[i].first, parent ) * tmpstrainrate ) ;
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

void GeneralizedSpaceTimeViscoElasticElementState::getField ( FieldType f, const Point & p, Vector & ret, bool local, VirtualMachine * vm , int )  const
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

    Form * visco = ( parent->getBehaviour() ) ;
    bool cleanup = !vm ;
    if ( !vm )
    {
        vm = new VirtualMachine() ;
    }
    switch ( f )
    {
    case SCALAR_DAMAGE_FIELD:
	if(parent->getBehaviour()->getDamageModel())
		ret[0] = parent->getBehaviour()->getDamageModel()->getState().max() ;
	else
		ret[0] = 0. ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case DISPLACEMENT_FIELD:
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
    case GENERALIZED_VISCOELASTIC_DISPLACEMENT_FIELD:
//			std::cout << ret.size() << "\t" << totaldof << std::endl ;
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
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
    case ENRICHED_DISPLACEMENT_FIELD:
        if ( !vm )
        {
            vm = new VirtualMachine() ;
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
    case GENERALIZED_VISCOELASTIC_ENRICHED_DISPLACEMENT_FIELD:
        if ( !vm )
        {
            vm = new VirtualMachine() ;
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
    case SPEED_FIELD:
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
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
    case GENERALIZED_VISCOELASTIC_SPEED_FIELD:
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
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
    case STRAIN_FIELD:
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
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

            for ( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() && j < enrichedDisplacements.size() * 2; j++ )
            {
                double f_xi = vm->deval ( parent->getEnrichmentFunction ( j ), XI, p_ ) ;
                double f_eta = vm->deval ( parent->getEnrichmentFunction ( j ), ETA, p_ ) ;

                x_xi += f_xi * enrichedDisplacements[j * totaldof] ;
                x_eta += f_eta * enrichedDisplacements[j * totaldof] ;
                y_xi += f_xi * enrichedDisplacements[j * totaldof + 1] ;
                y_eta += f_eta * enrichedDisplacements[j * totaldof + 1] ;
            }

            Matrix Jinv ( 3,3 ) ;
            parent->getInverseJacobianMatrix ( p_, Jinv ) ;
            ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;//+ x_tau * Jinv[0][2] ;
            ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;//+ y_tau * Jinv[1][2] ;
            ret[2] = 0.5 * ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] ) ;// + x_tau * Jinv[1][2]  + y_tau * Jinv[0][2] );

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

            Matrix Jinv ( 4,4 );
            parent->getInverseJacobianMatrix ( p_, Jinv ) ;
            ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1]  + ( x_zeta ) * Jinv[0][2];
            ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1]  + ( y_zeta ) * Jinv[1][2];
            ret[2] = ( z_xi ) * Jinv[2][0] + ( z_eta ) * Jinv[2][1]  + ( z_zeta ) * Jinv[2][2];

            ret[3] = 0.5 * ( ( y_xi ) * Jinv[2][0] +
                             ( y_eta ) * Jinv[2][1] +
                             ( y_zeta ) * Jinv[2][2] +
                             ( z_xi ) * Jinv[1][0] +
                             ( z_eta ) * Jinv[1][1] +
                             ( z_zeta ) * Jinv[1][2] );

            ret[4] = 0.5 * ( ( x_xi ) * Jinv[2][0] +
                             ( x_eta ) * Jinv[2][1] +
                             ( x_zeta ) * Jinv[2][2] +
                             ( z_xi ) * Jinv[0][0] +
                             ( z_eta ) * Jinv[0][1] +
                             ( z_zeta ) * Jinv[0][2] );

            ret[5] = 0.5 * ( ( y_xi )   * Jinv[0][0] +
                             ( y_eta )  * Jinv[0][1] +
                             ( y_zeta ) * Jinv[0][2] +
                             ( x_xi )   * Jinv[1][0] +
                             ( x_eta )  * Jinv[1][1] +
                             ( x_zeta ) * Jinv[1][2] );
        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case MECHANICAL_STRAIN_FIELD:
	getField( STRAIN_FIELD, p_, ret, true, vm ) ;
	if(getParent()->getBehaviour() && getParent()->getBehaviour()->hasInducedForces())
		ret -= getParent()->getBehaviour()->getImposedStrain(p_) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case PRINCIPAL_MECHANICAL_STRAIN_FIELD:
    {
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        Vector strains ( 0.,3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ;
        this->getField ( MECHANICAL_STRAIN_FIELD, p_, strains, true,vm ) ;
        ret = toPrincipal ( strains ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case GENERALIZED_VISCOELASTIC_STRAIN_FIELD:
        if ( !vm )
        {
            vm = new VirtualMachine() ;
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
                double f_xi = vm->deval ( parent->getShapeFunction ( j ), XI, p_ ) ;
                double f_eta = vm->deval ( parent->getShapeFunction ( j ), ETA, p_ ) ;
                double f_tau = vm->deval ( parent->getShapeFunction ( j ), TIME_VARIABLE, p_ ) ;
                for ( int i = 0 ; i < totaldof ; i++ )
                {
                    dx[i] += f_xi  * displacements[j * totaldof + i] ;
                    dy[i] += f_eta * displacements[j * totaldof + i] ;
                    dt[i] += f_tau * displacements[j * totaldof + i] ;
                }

                /*					x_xi += f_xi * displacements[j * totaldof] ;
                					x_eta += f_eta * displacements[j * totaldof] ;
                					y_xi += f_xi * displacements[j * totaldof + 1] ;
                					y_eta += f_eta * displacements[j * totaldof + 1] ;*/
            }

            for ( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() && j < enrichedDisplacements.size() * 2; j++ )
            {
                double f_xi = vm->deval ( parent->getEnrichmentFunction ( j ), XI, p_ ) ;
                double f_eta = vm->deval ( parent->getEnrichmentFunction ( j ), ETA, p_ ) ;
                double f_tau = vm->deval ( parent->getEnrichmentFunction ( j ), TIME_VARIABLE, p_ ) ;
                for ( int i = 0 ; i < totaldof ; i++ )
                {
                    dx[i] += f_xi  * enrichedDisplacements[j * totaldof + i] ;
                    dy[i] += f_eta * enrichedDisplacements[j * totaldof + i] ;
                    dt[i] += f_tau * enrichedDisplacements[j * totaldof + i] ;
                }

                /*					x_xi += f_xi * enrichedDisplacements[j * totaldof] ;
                					x_eta += f_eta * enrichedDisplacements[j * totaldof] ;
                					y_xi += f_xi * enrichedDisplacements[j * totaldof + 1] ;
                					y_eta += f_eta * enrichedDisplacements[j * totaldof + 1] ;*/
            }

            Matrix Jinv ( 3,3 ) ;
            parent->getInverseJacobianMatrix ( p_, Jinv ) ;
            for ( int i = 0 ; i < blocks ; i++ )
            {
                x_xi  = dx[ i * realdof + 0 ] ;
                x_eta = dy[ i * realdof + 0 ] ;
                y_xi  = dx[ i * realdof + 1 ] ;
                y_eta = dy[ i * realdof + 1 ] ;
                ret[i*3+0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;//+ x_tau * Jinv[0][2];
                ret[i*3+1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;//+ y_tau * Jinv[1][2] ;
                ret[i*3+2] = 0.5 * ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] );//+ x_tau * Jinv[1][2]  + y_tau * Jinv[0][2]);
//					std::cout << ret.size() << std::endl ;
            }
// 				Vector strainns(3) ;
// 				this->getField(STRAIN_FIELD, p_, strainns, true);
//
// 				std::cout << parent->getEnrichmentFunctions().size()<<"-" << ret[0] - strainns[0] << " ; " ;

// 				ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;
// 				ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;
// 				ret[2] = 0.5 * ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] );
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

            Matrix Jinv ( 4, 4 ) ;
            parent->getInverseJacobianMatrix ( p_, Jinv ) ;
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
                ret[i*6+0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1]  + ( x_zeta ) * Jinv[0][2];
                ret[i*6+1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1]  + ( y_zeta ) * Jinv[1][2];
                ret[i*6+2] = ( z_xi ) * Jinv[2][0] + ( z_eta ) * Jinv[2][1]  + ( z_zeta ) * Jinv[2][2];

                ret[i*6+3] = 0.5 * ( ( y_xi ) * Jinv[2][0] +
                                     ( y_eta ) * Jinv[2][1] +
                                     ( y_zeta ) * Jinv[2][2] +
                                     ( z_xi ) * Jinv[1][0] +
                                     ( z_eta ) * Jinv[1][1] +
                                     ( z_zeta ) * Jinv[1][2] );

                ret[i*6+4] = 0.5 * ( ( x_xi ) * Jinv[2][0] +
                                     ( x_eta ) * Jinv[2][1] +
                                     ( x_zeta ) * Jinv[2][2] +
                                     ( z_xi ) * Jinv[0][0] +
                                     ( z_eta ) * Jinv[0][1] +
                                     ( z_zeta ) * Jinv[0][2] );

                ret[i*6+5] = 0.5 * ( ( y_xi )   * Jinv[0][0] +
                                     ( y_eta )  * Jinv[0][1] +
                                     ( y_zeta ) * Jinv[0][2] +
                                     ( x_xi )   * Jinv[1][0] +
                                     ( x_eta )  * Jinv[1][1] +
                                     ( x_zeta ) * Jinv[1][2] );
            }

        }

        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case PRINCIPAL_STRAIN_FIELD:
    {
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        Vector strains ( 0.,3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ;
        this->getField ( STRAIN_FIELD, p_, strains, true,vm ) ;
        ret = toPrincipal ( strains ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case GENERALIZED_VISCOELASTIC_PRINCIPAL_STRAIN_FIELD:
    {
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        Vector strains ( 0.,blocks* ( 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ) ;
        this->getField ( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, p_, strains, true,vm ) ;
        Vector tmp ( 0., 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ;
        Vector ptmp ( 0., 2+ ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ;
        for ( int i = 0 ; i < blocks ; i++ )
        {
            for ( size_t j = 0 ; j < tmp.size() ; j++ )
            {
                tmp[j] = strains[ i*tmp.size() + j ] ;
            }
            ptmp = toPrincipal ( tmp ) ;
            for ( size_t j = 0 ; j < tmp.size() ; j++ )
            {
                ret[ i*ptmp.size() + j] = ptmp[j] ;
            }
        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case NON_ENRICHED_STRAIN_FIELD:
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
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

            Matrix Jinv ( 3,3 ) ;
            parent->getInverseJacobianMatrix ( p_, Jinv ) ;
            ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;
            ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;
            ret[2] = 0.5 * ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] );
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

            Matrix Jinv ( 4, 4 ) ;
            parent->getInverseJacobianMatrix ( p_, Jinv ) ;
            ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1]  + ( x_zeta ) * Jinv[0][2];
            ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1]  + ( y_zeta ) * Jinv[1][2];
            ret[2] = ( z_xi ) * Jinv[2][0] + ( z_eta ) * Jinv[2][1]  + ( z_zeta ) * Jinv[2][2];

            ret[3] = 0.5 * ( ( y_xi ) * Jinv[2][0] +
                             ( y_eta ) * Jinv[2][1] +
                             ( y_zeta ) * Jinv[2][2] +
                             ( z_xi ) * Jinv[1][0] +
                             ( z_eta ) * Jinv[1][1] +
                             ( z_zeta ) * Jinv[1][2] );

            ret[4] = 0.5 * ( ( x_xi ) * Jinv[2][0] +
                             ( x_eta ) * Jinv[2][1] +
                             ( x_zeta ) * Jinv[2][2] +
                             ( z_xi ) * Jinv[0][0] +
                             ( z_eta ) * Jinv[0][1] +
                             ( z_zeta ) * Jinv[0][2] );

            ret[5] = 0.5 * ( ( y_xi )   * Jinv[0][0] +
                             ( y_eta )  * Jinv[0][1] +
                             ( y_zeta ) * Jinv[0][2] +
                             ( x_xi )   * Jinv[1][0] +
                             ( x_eta )  * Jinv[1][1] +
                             ( x_zeta ) * Jinv[1][2] );
        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_FIELD:
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
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

                /*					x_xi += f_xi * displacements[j * totaldof] ;
                					x_eta += f_eta * displacements[j * totaldof] ;
                					y_xi += f_xi * displacements[j * totaldof + 1] ;
                					y_eta += f_eta * displacements[j * totaldof + 1] ;*/
            }

            Matrix Jinv ( 3,3 ) ;
            parent->getInverseJacobianMatrix ( p_, Jinv ) ;
            for ( int i = 0 ; i < blocks ; i++ )
            {
                x_xi  = dx[ i * realdof + 0 ] ;
                x_eta = dy[ i * realdof + 0 ] ;
                y_xi  = dx[ i * realdof + 1 ] ;
                y_eta = dy[ i * realdof + 1 ] ;
                ret[i*3+0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;
                ret[i*3+1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;
                ret[i*3+2] = 0.5 * ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] );
            }
// 				ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;
// 				ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;
// 				ret[2] = 0.5 * ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] );
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


            Matrix Jinv ( 4, 4 ) ;
            parent->getInverseJacobianMatrix ( p_, Jinv ) ;
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
                ret[i*6+0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1]  + ( x_zeta ) * Jinv[0][2];
                ret[i*6+1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1]  + ( y_zeta ) * Jinv[1][2];
                ret[i*6+2] = ( z_xi ) * Jinv[2][0] + ( z_eta ) * Jinv[2][1]  + ( z_zeta ) * Jinv[2][2];

                ret[i*6+3] = 0.5 * ( ( y_xi ) * Jinv[2][0] +
                                     ( y_eta ) * Jinv[2][1] +
                                     ( y_zeta ) * Jinv[2][2] +
                                     ( z_xi ) * Jinv[1][0] +
                                     ( z_eta ) * Jinv[1][1] +
                                     ( z_zeta ) * Jinv[1][2] );

                ret[i*6+4] = 0.5 * ( ( x_xi ) * Jinv[2][0] +
                                     ( x_eta ) * Jinv[2][1] +
                                     ( x_zeta ) * Jinv[2][2] +
                                     ( z_xi ) * Jinv[0][0] +
                                     ( z_eta ) * Jinv[0][1] +
                                     ( z_zeta ) * Jinv[0][2] );

                ret[i*6+5] = 0.5 * ( ( y_xi )   * Jinv[0][0] +
                                     ( y_eta )  * Jinv[0][1] +
                                     ( y_zeta ) * Jinv[0][2] +
                                     ( x_xi )   * Jinv[1][0] +
                                     ( x_eta )  * Jinv[1][1] +
                                     ( x_zeta ) * Jinv[1][2] );
            }

        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case VON_MISES_STRAIN_FIELD:
    {
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        Vector eps ( 0., ( size_t ) parent->spaceDimensions() ) ;
        this->getField ( PRINCIPAL_STRAIN_FIELD, p_, eps, true ,vm ) ;
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
        return ;
    }
    case STRAIN_RATE_FIELD:
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
        {
		// 				if(totaldof == 2)
// 				{
// 					this->getField( STRAIN_RATE_FIELD, p_, ret, true) ;
// 					return ;
// 				}


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

                /*					x_xi += f_xi * displacements[j * totaldof] ;
                					x_eta += f_eta * displacements[j * totaldof] ;
                					y_xi += f_xi * displacements[j * totaldof + 1] ;
                					y_eta += f_eta * displacements[j * totaldof + 1] ;*/
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

                /*					x_xi += f_xi * enrichedDisplacements[j * totaldof] ;
                					x_eta += f_eta * enrichedDisplacements[j * totaldof] ;
                					y_xi += f_xi * enrichedDisplacements[j * totaldof + 1] ;
                					y_eta += f_eta * enrichedDisplacements[j * totaldof + 1] ;*/
            }

            Matrix Jinv ( 3,3 ) ;
            parent->getInverseJacobianMatrix ( p_, Jinv ) ;
                x_xi  = dx[ 0 ] ;
                x_eta = dy[ 0 ] ;
                y_xi  = dx[ 1 ] ;
                y_eta = dy[ 1 ] ;
                ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;//+ x_tau * Jinv[0][2];
                ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;//+ y_tau * Jinv[1][2] ;
                ret[2] = 0.5 * ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] );//+ x_tau * Jinv[1][2]  + y_tau * Jinv[0][2]);
                //					std::cout << ret.size() << std::endl ;
            ret *= Jinv[2][2] ;

// 				ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;
// 				ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;
// 				ret[2] = 0.5 * ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] );
//
// 				ret *= Jinv[2][2] ;
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

            Matrix Jinv ( 4, 4 ) ;
            parent->getInverseJacobianMatrix ( p_, Jinv ) ;
            ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1]  + ( x_zeta ) * Jinv[0][2];
            ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1]  + ( y_zeta ) * Jinv[1][2];
            ret[2] = ( z_xi ) * Jinv[2][0] + ( z_eta ) * Jinv[2][1]  + ( z_zeta ) * Jinv[2][2];

            ret[3] = 0.5 * ( ( y_xi ) * Jinv[2][0] +
                             ( y_eta ) * Jinv[2][1] +
                             ( y_zeta ) * Jinv[2][2] +
                             ( z_xi ) * Jinv[1][0] +
                             ( z_eta ) * Jinv[1][1] +
                             ( z_zeta ) * Jinv[1][2] );

            ret[4] = 0.5 * ( ( x_xi ) * Jinv[2][0] +
                             ( x_eta ) * Jinv[2][1] +
                             ( x_zeta ) * Jinv[2][2] +
                             ( z_xi ) * Jinv[0][0] +
                             ( z_eta ) * Jinv[0][1] +
                             ( z_zeta ) * Jinv[0][2] );

            ret[5] = 0.5 * ( ( y_xi )   * Jinv[0][0] +
                             ( y_eta )  * Jinv[0][1] +
                             ( y_zeta ) * Jinv[0][2] +
                             ( x_xi )   * Jinv[1][0] +
                             ( x_eta )  * Jinv[1][1] +
                             ( x_zeta ) * Jinv[1][2] );

            ret *= Jinv[3][3] ;
        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD:
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
        {
// 				if(totaldof == 2)
// 				{
// 					this->getField( STRAIN_RATE_FIELD, p_, ret, true) ;
// 					return ;
// 				}


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

                /*					x_xi += f_xi * displacements[j * totaldof] ;
                					x_eta += f_eta * displacements[j * totaldof] ;
                					y_xi += f_xi * displacements[j * totaldof + 1] ;
                					y_eta += f_eta * displacements[j * totaldof + 1] ;*/
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

                /*					x_xi += f_xi * enrichedDisplacements[j * totaldof] ;
                					x_eta += f_eta * enrichedDisplacements[j * totaldof] ;
                					y_xi += f_xi * enrichedDisplacements[j * totaldof + 1] ;
                					y_eta += f_eta * enrichedDisplacements[j * totaldof + 1] ;*/
            }

            Matrix Jinv ( 3,3 ) ;
            parent->getInverseJacobianMatrix ( p_, Jinv ) ;
            for ( int i = 0 ; i < blocks ; i++ )
            {
                x_xi  = dx[ i * realdof + 0 ] ;
                x_eta = dy[ i * realdof + 0 ] ;
                y_xi  = dx[ i * realdof + 1 ] ;
                y_eta = dy[ i * realdof + 1 ] ;
                ret[i*3+0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;//+ x_tau * Jinv[0][2];
                ret[i*3+1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;//+ y_tau * Jinv[1][2] ;
                ret[i*3+2] = 0.5 * ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] );//+ x_tau * Jinv[1][2]  + y_tau * Jinv[0][2]);
                //					std::cout << ret.size() << std::endl ;
            }
            ret *= Jinv[2][2] ;
//	    std::cout << 2./Jinv[2][2] << "    ";
// 				ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;
// 				ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;
// 				ret[2] = 0.5 * ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] );
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

            Matrix Jinv ( 4, 4 ) ;
            parent->getInverseJacobianMatrix ( p_, Jinv ) ;
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
                ret[i*6+0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1]  + ( x_zeta ) * Jinv[0][2];
                ret[i*6+1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1]  + ( y_zeta ) * Jinv[1][2];
                ret[i*6+2] = ( z_xi ) * Jinv[2][0] + ( z_eta ) * Jinv[2][1]  + ( z_zeta ) * Jinv[2][2];

                ret[i*6+3] = 0.5 * ( ( y_xi ) * Jinv[2][0] +
                                     ( y_eta ) * Jinv[2][1] +
                                     ( y_zeta ) * Jinv[2][2] +
                                     ( z_xi ) * Jinv[1][0] +
                                     ( z_eta ) * Jinv[1][1] +
                                     ( z_zeta ) * Jinv[1][2] );

                ret[i*6+4] = 0.5 * ( ( x_xi ) * Jinv[2][0] +
                                     ( x_eta ) * Jinv[2][1] +
                                     ( x_zeta ) * Jinv[2][2] +
                                     ( z_xi ) * Jinv[0][0] +
                                     ( z_eta ) * Jinv[0][1] +
                                     ( z_zeta ) * Jinv[0][2] );

                ret[i*6+5] = 0.5 * ( ( y_xi )   * Jinv[0][0] +
                                     ( y_eta )  * Jinv[0][1] +
                                     ( y_zeta ) * Jinv[0][2] +
                                     ( x_xi )   * Jinv[1][0] +
                                     ( x_eta )  * Jinv[1][1] +
                                     ( x_zeta ) * Jinv[1][2] );
            }
            ret *= Jinv[3][3] ;

        }

        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case NON_ENRICHED_STRAIN_RATE_FIELD:
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
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

            Matrix Jinv ( 3,3 ) ;
            parent->getInverseJacobianMatrix ( p_, Jinv ) ;
            ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;
            ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;
            ret[2] = 0.5 * ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] );

            ret *= Jinv[2][2] ;
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

            Matrix Jinv ( 4, 4 ) ;
            parent->getInverseJacobianMatrix ( p_, Jinv ) ;
            ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1]  + ( x_zeta ) * Jinv[0][2];
            ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1]  + ( y_zeta ) * Jinv[1][2];
            ret[2] = ( z_xi ) * Jinv[2][0] + ( z_eta ) * Jinv[2][1]  + ( z_zeta ) * Jinv[2][2];

            ret[3] = 0.5 * ( ( y_xi ) * Jinv[2][0] +
                             ( y_eta ) * Jinv[2][1] +
                             ( y_zeta ) * Jinv[2][2] +
                             ( z_xi ) * Jinv[1][0] +
                             ( z_eta ) * Jinv[1][1] +
                             ( z_zeta ) * Jinv[1][2] );

            ret[4] = 0.5 * ( ( x_xi ) * Jinv[2][0] +
                             ( x_eta ) * Jinv[2][1] +
                             ( x_zeta ) * Jinv[2][2] +
                             ( z_xi ) * Jinv[0][0] +
                             ( z_eta ) * Jinv[0][1] +
                             ( z_zeta ) * Jinv[0][2] );

            ret[5] = 0.5 * ( ( y_xi )   * Jinv[0][0] +
                             ( y_eta )  * Jinv[0][1] +
                             ( y_zeta ) * Jinv[0][2] +
                             ( x_xi )   * Jinv[1][0] +
                             ( x_eta )  * Jinv[1][1] +
                             ( x_zeta ) * Jinv[1][2] );

            ret *= Jinv[3][3] ;
        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_RATE_FIELD:
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
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

            Matrix Jinv ( 3,3 ) ;
            parent->getInverseJacobianMatrix ( p_, Jinv ) ;
            for ( int i = 0 ; i < blocks ; i++ )
            {
                x_xi  = dx[ i * realdof + 0 ] ;
                x_eta = dy[ i * realdof + 0 ] ;
                y_xi  = dx[ i * realdof + 1 ] ;
                y_eta = dy[ i * realdof + 1 ] ;
                ret[i*3+0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;
                ret[i*3+1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;
                ret[i*3+2] = 0.5 * ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] );
            }
            ret *= Jinv[2][2] ;
// 				ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;
// 				ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;
// 				ret[2] = 0.5 * ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] );
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

            Matrix Jinv ( 4, 4 ) ;
            parent->getInverseJacobianMatrix ( p_, Jinv ) ;
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
                ret[i*6+0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1]  + ( x_zeta ) * Jinv[0][2];
                ret[i*6+1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1]  + ( y_zeta ) * Jinv[1][2];
                ret[i*6+2] = ( z_xi ) * Jinv[2][0] + ( z_eta ) * Jinv[2][1]  + ( z_zeta ) * Jinv[2][2];

                ret[i*6+3] = 0.5 * ( ( y_xi ) * Jinv[2][0] +
                                     ( y_eta ) * Jinv[2][1] +
                                     ( y_zeta ) * Jinv[2][2] +
                                     ( z_xi ) * Jinv[1][0] +
                                     ( z_eta ) * Jinv[1][1] +
                                     ( z_zeta ) * Jinv[1][2] );

                ret[i*6+4] = 0.5 * ( ( x_xi ) * Jinv[2][0] +
                                     ( x_eta ) * Jinv[2][1] +
                                     ( x_zeta ) * Jinv[2][2] +
                                     ( z_xi ) * Jinv[0][0] +
                                     ( z_eta ) * Jinv[0][1] +
                                     ( z_zeta ) * Jinv[0][2] );

                ret[i*6+5] = 0.5 * ( ( y_xi )   * Jinv[0][0] +
                                     ( y_eta )  * Jinv[0][1] +
                                     ( y_zeta ) * Jinv[0][2] +
                                     ( x_xi )   * Jinv[1][0] +
                                     ( x_eta )  * Jinv[1][1] +
                                     ( x_zeta ) * Jinv[1][2] );
            }
            ret *= Jinv[3][3] ;

        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case REAL_STRESS_FIELD:
    {
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        Vector strains ( 0., blocks* ( 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ) ;
        Vector speeds ( 0., blocks* ( 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ) ;
        getField ( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, p_, strains, true,vm ) ;
        getField ( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, p_, speeds, true,vm ) ;
        Vector stresses = ( Vector ) ( visco->getTensor ( p_, parent ) * strains )
                          + ( Vector ) ( visco->getViscousTensor ( p_, parent ) * speeds ) ;
        for ( size_t i = 0 ; i < ret.size() ; i++ )
        {
            ret[i] = stresses[i] ;
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
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
// 			visco->getTensor(p_,parent).print() ;
// 			visco->getViscousTensor(p_, parent).print() ;
// 			exit(0) ;
        Vector strains ( 0., blocks* ( 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ) ;
        Vector speeds ( 0., blocks* ( 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ) ;
        getField ( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, p_, strains, true,vm ) ;
        getField ( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, p_, speeds, true,vm ) ;
        ret = ( Vector ) ( visco->getTensor ( p_, parent ) * strains )
              + ( Vector ) ( visco->getViscousTensor ( p_, parent ) * speeds ) ;
        Vector stresses = visco->getImposedStress ( p_, parent ) ;
        for ( size_t i = 0 ; i < stresses.size() ; i++ )
        {
            ret[i] -= stresses[i] ;
	    for(int j = 1 ; j < blocks ; j++)
		ret[i+j*stresses.size()] += stresses[i] ;
        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case PRINCIPAL_REAL_STRESS_FIELD:
    {
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        Vector stress ( 0.,3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ;
        this->getField ( REAL_STRESS_FIELD, p_, stress, true,vm ) ;
        ret = toPrincipal ( stress ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case GENERALIZED_VISCOELASTIC_PRINCIPAL_REAL_STRESS_FIELD:
    {
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        Vector stress ( 0.,blocks* ( 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ) ;
        this->getField ( GENERALIZED_VISCOELASTIC_REAL_STRESS_FIELD, p_, stress, true,vm ) ;
        Vector tmp ( 0., 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ;
        Vector ptmp ( 0., 2+ ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ;
        for ( int i = 0 ; i < blocks ; i++ )
        {
            for ( size_t j = 0 ; j < tmp.size() ; j++ )
            {
                tmp[j] = stress[ i*tmp.size() + j ] ;
            }
            ptmp = toPrincipal ( tmp ) ;
            for ( size_t j = 0 ; j < tmp.size() ; j++ )
            {
                ret[ i*ptmp.size() + j] = ptmp[j] ;
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
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        Vector strains ( 0., blocks* ( 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ) ;
        Vector speeds ( 0., blocks* ( 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ) ;
        this->getField ( GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_FIELD, p_, strains, true,vm ) ;
        this->getField ( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, p_, speeds, true,vm ) ;
        Vector imposed = visco->getImposedStrain ( p_, parent ) ;
        for ( size_t i = 0 ; i < strains.size() ; i++ )
        {
            strains[i] -= imposed[i] ;
        }
        Vector stresses = ( Vector ) ( visco->getTensor ( p_, parent ) * strains )
                          + ( Vector ) ( visco->getViscousTensor ( p_, parent ) * speeds ) ;
        for ( size_t i = 0 ; i < 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ; i++ )
        {
            ret[i] = stresses[i] ;
        }
//			ret -= getParent()->getBehaviour()->getImposedStress(p_, parent) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case GENERALIZED_VISCOELASTIC_NON_ENRICHED_REAL_STRESS_FIELD:
    {
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        Vector strains ( 0., blocks* ( 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ) ;
        Vector speeds ( 0., blocks* ( 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ) ;
        this->getField ( GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_FIELD, p_, strains, true,vm ) ;
        this->getField ( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, p_, speeds, true,vm ) ;
        ret = ( Vector ) ( visco->getTensor ( p_, parent ) * strains )
              + ( Vector ) ( visco->getViscousTensor ( p_, parent ) * speeds ) ;
        Vector stresses = visco->getImposedStress ( p_, parent ) ;
        for ( size_t i = 0 ; i < stresses.size() ; i++ )
        {
            ret[i] -= stresses[i] ;
        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case VON_MISES_REAL_STRESS_FIELD:

        if ( parent->getOrder() == LINEAR )
        {
            if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
            {
                if ( !vm )
                {
                    vm = new VirtualMachine() ;
                }
                Vector sigma ( 0., 2 ) ;
                Point c ( 1./3., 1./3. ) ;
                this->getField ( PRINCIPAL_REAL_STRESS_FIELD, c, sigma, true,vm ) ;
                ret[0] = sqrt ( ( ( sigma[0] - sigma[1] ) * ( sigma[0] - sigma[1] ) + sigma[0] * sigma[0] + sigma[1] * sigma[1] ) / 2. ) ;
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
    case EFFECTIVE_STRESS_FIELD:
    {
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        Vector strains ( 0., blocks* ( 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ) ;
        Vector speeds ( 0., blocks* ( 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ) ;
        this->getField ( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, p_, strains, true,vm ) ;
        this->getField ( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, p_, speeds, true,vm ) ;
        Vector stresses = ( Vector ) ( visco->param * strains )
                          + ( Vector ) ( visco->getViscousTensor ( p_, parent ) * speeds ) ;
        for ( size_t i = 0 ; i < 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ; i++ )
        {
            ret[i] = stresses[i] ;
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
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        Vector strains ( 0., blocks* ( 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ) ;
        Vector speeds ( 0., blocks* ( 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ) ;
        this->getField ( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, p_, strains, true,vm ) ;
        this->getField ( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, p_, speeds, true,vm ) ;
        ret = ( Vector ) ( visco->param * strains )
              + ( Vector ) ( visco->getViscousTensor ( p_, parent ) * speeds ) ;
        Vector stresses = visco->getImposedStress ( p_, parent ) ;
        for ( size_t i = 0 ; i < stresses.size() ; i++ )
        {
            ret[i] -= stresses[i] ;
        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case PRINCIPAL_EFFECTIVE_STRESS_FIELD:
    {
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        Vector stress ( 0.,3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ;
        this->getField ( EFFECTIVE_STRESS_FIELD, p_, stress, true,vm ) ;
        ret = toPrincipal ( stress ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case GENERALIZED_VISCOELASTIC_PRINCIPAL_EFFECTIVE_STRESS_FIELD:
    {
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        Vector stress ( 0.,blocks* ( 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ) ;
        this->getField ( GENERALIZED_VISCOELASTIC_EFFECTIVE_STRESS_FIELD, p_, stress, true,vm ) ;
        Vector tmp ( 0., 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ;
        Vector ptmp ( 0., 2+ ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ;
        for ( int i = 0 ; i < blocks ; i++ )
        {
            for ( size_t j = 0 ; j < tmp.size() ; j++ )
            {
                tmp[j] = stress[ i*tmp.size() + j ] ;
            }
            ptmp = toPrincipal ( tmp ) ;
            for ( size_t j = 0 ; j < tmp.size() ; j++ )
            {
                ret[ i*ptmp.size() + j] = ptmp[j] ;
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
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        Vector strains ( 0., blocks* ( 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ) ;
        Vector speeds ( 0., blocks* ( 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ) ;
        this->getField ( GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_FIELD, p_, strains, true,vm ) ;
        this->getField ( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, p_, speeds, true,vm ) ;
        Vector stresses = ( Vector ) ( visco->param * strains )
                          + ( Vector ) ( visco->getViscousTensor ( p_, parent ) * speeds ) ;
        for ( size_t i = 0 ; i < 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ; i++ )
        {
            ret[i] = stresses[i] ;
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
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        Vector strains ( 0., blocks* ( 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ) ;
        Vector speeds ( 0., blocks* ( 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ) ;
        this->getField ( GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_FIELD, p_, strains, true,vm ) ;
        this->getField ( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, p_, speeds, true,vm ) ;
        ret = ( Vector ) ( visco->param * strains )
              + ( Vector ) ( visco->getViscousTensor ( p_, parent ) * speeds ) ;
        Vector stresses = visco->getImposedStress ( p_, parent ) ;
        for ( size_t i = 0 ; i < stresses.size() ; i++ )
        {
            ret[i] -= stresses[i] ;
        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case VON_MISES_EFFECTIVE_STRESS_FIELD:
        if ( parent->getOrder() == LINEAR )
        {
            if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
            {
                if ( !vm )
                {
                    vm = new VirtualMachine() ;
                }
                Vector sigma ( 0., 2 ) ;
                Point c ( 1./3., 1./3. ) ;
                this->getField ( PRINCIPAL_EFFECTIVE_STRESS_FIELD, c, sigma, true,vm ) ;
                ret[0] = sqrt ( ( ( sigma[0] - sigma[1] ) * ( sigma[0] - sigma[1] ) + sigma[0] * sigma[0] + sigma[1] * sigma[1] ) / 2. ) ;
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
    case PRINCIPAL_STRESS_ANGLE_FIELD:
    {
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        Vector strains ( 0., 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ;
        this->getField ( REAL_STRESS_FIELD,  p_, strains, true, vm ) ;
        if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
        {
            ret[0] =  0.5 * atan2 ( strains[2], strains[0] - strains[1] ) ;
        }
        else
        {
            ret[0] =  0.5 * atan2 ( strains[3] , strains[0] - strains[1] ) ;
            ret[1] =  0.5 * atan2 ( strains[4] , strains[0] - strains[2] ) ;
            ret[2] =  0.5 * atan2 ( strains[5] , strains[1] - strains[2] ) ;
        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case PRINCIPAL_STRAIN_ANGLE_FIELD:
    {
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        Vector strains ( 0., 3+3* ( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) ) ;
        this->getField ( STRAIN_FIELD,  p_, strains, true, vm ) ;
        if ( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
        {
            ret[0] =  0.5 * atan2 ( strains[2], strains[0] - strains[1] ) ;
        }
        else
        {
            ret[0] =  0.5 * atan2 ( strains[3] , strains[0] - strains[1] ) ;
            ret[1] =  0.5 * atan2 ( strains[4] , strains[0] - strains[2] ) ;
            ret[2] =  0.5 * atan2 ( strains[5] , strains[1] - strains[2] ) ;
        }
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    }
    case GRADIENT_FIELD:
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

            Matrix Jinv( 2, 2 ) ;
            parent->getInverseJacobianMatrix( p_, Jinv ) ;
            ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;
            ret[1] = ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1] ;
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
        return ;
    case FLUX_FIELD:
        this->getField(GRADIENT_FIELD, p_, ret, true) ;
        ret = (Vector) (parent->getBehaviour()->getTensor(p_, parent) * ret) ;
        return ;
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
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        ElementState::getField ( f, parent->getBoundingPoints(), ret, false,vm ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case MECHANICAL_STRAIN_FIELD :
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        ElementState::getField ( f, parent->getBoundingPoints(), ret, false,vm ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case PRINCIPAL_STRAIN_FIELD :
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        ElementState::getField ( f, parent->getBoundingPoints(), ret, false,vm ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case PRINCIPAL_MECHANICAL_STRAIN_FIELD :
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        ElementState::getField ( f, parent->getBoundingPoints(), ret, false,vm ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case NON_ENRICHED_STRAIN_FIELD :
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        ElementState::getField ( f, parent->getBoundingPoints(), ret, false,vm ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case REAL_STRESS_FIELD:
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        ElementState::getField ( f,  parent->getBoundingPoints(), ret, false,vm ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case PRINCIPAL_REAL_STRESS_FIELD :
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        ElementState::getField ( f, parent->getBoundingPoints(), ret, false,vm ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case NON_ENRICHED_REAL_STRESS_FIELD:
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        ElementState::getField ( f,  parent->getBoundingPoints(), ret, false,vm ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case EFFECTIVE_STRESS_FIELD:
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        ElementState::getField ( f,  parent->getBoundingPoints(), ret, false,vm ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case PRINCIPAL_EFFECTIVE_STRESS_FIELD :
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        ElementState::getField ( f, parent->getBoundingPoints(), ret, false,vm ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    case NON_ENRICHED_EFFECTIVE_STRESS_FIELD:
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
        ElementState::getField ( f,  parent->getBoundingPoints(), ret, false,vm ) ;
        if ( cleanup )
        {
            delete vm ;
        }
        return ;
    default:
        if ( !vm )
        {
            vm = new VirtualMachine() ;
        }
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

void GeneralizedSpaceTimeViscoElasticElementState::getField ( FieldType f1, FieldType f2, const Point & p, Vector & ret1, Vector & ret2, bool local, VirtualMachine * vm, int , int )  const
{
    bool cleanup = !vm ;
    if ( !vm )
    {
        vm = new VirtualMachine() ;
    }
    this->getField ( f1, p, ret1, local, vm ) ;
    this->getField ( f2, p, ret2, local, vm ) ;
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
    this->getFieldAtNodes ( f1, ret1, vm ) ;
    this->getFieldAtNodes ( f2, ret2, vm ) ;
    if ( cleanup )
    {
        delete vm ;
    }
}



GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables::GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables ( IntegrableEntity * e, std::map<std::string, double> & v ) : GeneralizedSpaceTimeViscoElasticElementState ( e )
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
//    else
//    {
//        variables.insert ( std::pair<std::string, double> ( v, 0. ) ) ;
//    }
}

void GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables::synchronize(std::map<std::string, double> & defaultValues) 
{
	for( std::map<std::string, double>::iterator it = defaultValues.begin() ; it != defaultValues.end() ; it++)
	{
		variables[it->first] = it->second ;
	}
//	variables = defaultValues ;
}

