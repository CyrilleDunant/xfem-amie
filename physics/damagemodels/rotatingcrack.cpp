//
// C++ Implementation: isotropiclineardamage
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "rotatingcrack.h"
#include "damagemodel.h"
#include "../../elements/integrable_entity.h"
#include "../fracturecriteria/fracturecriterion.h"
#include "../orthotropicstiffness.h"
#include "../stiffness.h"

namespace Amie
{

RotatingCrack::RotatingCrack ( double E, double nu ) :  E ( E ), nu ( nu )
{
    getState ( true ).resize ( 4, 0. );
    getState(true)[0] = 1e-4 ;
    isNull = false ;
    currentAngle = M_PI*.25 ;
    initialAngle = M_PI*.25 ;
    factor = 1 ;
    es = nullptr ;
    needGlobalMaximumScore = true ;
    firstTension = false ;
    secondTension = false ;
    firstTensionFailure = false ;
    secondTensionFailure = false ;
    firstCompressionFailure = false ;
    secondCompressionFailure = false ;
    firstMet = false ;
    secondMet = false ;
    alternating = true ;
    alternate = true ;
    postprocheck = false ;
// 	ctype = CONSERVATIVE ;
    stiff = new OrthotropicStiffness ( E,E,E/ ( 1.-nu*nu ),nu, 0. ) ;
    roughsampling = true ;
    iterationcount = 0 ;
    iterationNumber = 48 ;
    damageDensityTolerance =  std::max(0.25/pow(2.,iterationNumber/4), 1e-4) ; //1e-8 ;//1. / pow( 2., 14 );
}


double damageAtAngle ( const std::vector<std::pair<double, double> > & increments , double angle )
{
    double ret = 0 ;
    double rettest = 0 ;

    for ( size_t i = 0 ; i <  increments.size() ; i++ )
    {
        rettest += increments[i].second ;
        if ( cos ( increments[i].first - angle ) > POINT_TOLERANCE )
        {
            double a = cos ( increments[i].first - angle ) * cos ( increments[i].first - angle ) ;
            ret +=  a * increments[i].second ;
        }
    }

    return ret ;
}

int RotatingCrack::getMode() const
{
    if ( es && es->getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet() &&
            ( (!firstTension && es->getParent()->getBehaviour()->getFractureCriterion()->directionInTension ( 0 ))
              || (firstTension && es->getParent()->getBehaviour()->getFractureCriterion()->directionInCompression ( 0 ))
              || (!secondTension && es->getParent()->getBehaviour()->getFractureCriterion()->directionInTension ( 1 ))
              || (secondTension && es->getParent()->getBehaviour()->getFractureCriterion()->directionInCompression ( 1 ))
              || firstMet != es->getParent()->getBehaviour()->getFractureCriterion()->directionMet ( 0 )
              || secondMet != es->getParent()->getBehaviour()->getFractureCriterion()->directionMet ( 1 ) )
       )
    {
        return 1 ;
    }
    return -1 ;
}

double RotatingCrack::getAngleShift() const
{
    return 0 ;
}


void RotatingCrack::step(ElementState & s, double maxscore)
{

    if(!alternate)
    {
        DamageModel::step(s, maxscore) ;
        
        if ( change )
        {

            double E_0 = E ;
            double E_1 = E ;
            double fs = firstTension  ? getState() [1] : getState() [1] ;
            double ss = secondTension ? getState() [2] : getState() [2] ;

            E_0 *=  1. - fs  ;
            E_1 *=  1. - ss  ;

            double nunu = nu ;
            if ( getState().max() > POINT_TOLERANCE )
            {
                nunu = nu*exp ( -1./ ( 1.-std::max ( fs,ss ) ) ) ;
                E_0 /= 1.-nunu*nunu ;
                E_1 /= 1.-nunu*nunu ;
            }

            double G = E_0*E_1/ ( E_0+E_1 ) ;
            if ( E_0+E_1 < POINT_TOLERANCE*E )
            {
                G = 0 ;
            }
            
            stiff->setStiffness ( E_0, E_1, G, nunu ) ;        
            
        } 
    }
    else
    {

        double max = maxscore ;
        if(!needGlobalMaximumScore)
            max = s.getParent()->getBehaviour()->getFractureCriterion()->getMaxScoreInNeighbourhood(s)  ;

        s.getParent()->getBehaviour()->getFractureCriterion()->setChange( s , max) ;
        change = false ;
        if(!s.getParent()->getBehaviour()->getFractureCriterion())
        {
            converged = true ;
            return ;
        }
        
        if( !s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet() )
        {
            s.getParent()->getBehaviour()->getFractureCriterion()->setCheckpoint( false );
            alternateCheckpoint = false ;
            // this is necessary because we want to trigger side-effects
            //for example, plasticstrain gets a pointer to s
            computeDamageIncrement( s ) ;
            converged = true ;
            return ;
        }
        
        if( s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint() ) // initiate iteration
        {
            angles_scores.clear();
            s.getParent()->getBehaviour()->getFractureCriterion()->setCheckpoint( false );
            alternateCheckpoint = false ;

            if(!fractured())
            {
                converged = false ;
                change = true ;
            }
            else
            {
                converged = true ;
                change = false ;
                return ;
            }

            
            initialAngle = currentAngle ;               
            roughsampling = true ;
            iterationcount = 0 ;         
            
//             if(firstTension)
                currentAngle = std::max(initialAngle-M_PI*.25, 0.) ;

            stiff->setAngle ( currentAngle ) ;   

            return ;

        }
        else
        {
           iterationcount++ ;
            angles_scores.push_back( std::make_pair(currentAngle, s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState()));
            std::sort(angles_scores.begin(),angles_scores.end());
            
            //first produce a rough sampling           
            double minAngle = std::max(initialAngle-M_PI*.25, 0.) ;
            double maxAngle = std::min(initialAngle+M_PI*.25, M_PI*.5) ;
            double deltaAngle = std::max(maxAngle-minAngle, 0.25*M_PI) ;
            if(angles_scores.back().first < maxAngle - .2*deltaAngle && roughsampling)
            {  
                currentAngle += .2*deltaAngle ;
//                 if(firstTension)
                    stiff->setAngle ( currentAngle ) ;     
                change = true ;
                return ;
            }
            else if(roughsampling)
            { 
                //find a triplet of points around the minimum               
                double phi = (1.+sqrt(5))*.5 ;
                std::vector<std::pair<double, double>> newset ;
                double minscore = angles_scores[0].second ;
                size_t minidex = 0 ;
                for(size_t i = 1 ; i < angles_scores.size() ; i++)
                {
                    if(angles_scores[i].second < minscore)
                    {
                        minscore = angles_scores[i].second ;
                        minidex = i ;
                    }
                }
                if(minidex == 0)
                {
                    newset.push_back(angles_scores.back()) ;
                    currentAngle = phi/(1.+phi) * angles_scores.back().first + 1./(1.+phi) * angles_scores[minidex+1].first;
                    newset.push_back(angles_scores[minidex+1]);
                }
                else if(minidex == (angles_scores.size()-1))
                {
                    newset.push_back(angles_scores[minidex-1]) ;
                    currentAngle = phi/(1.+phi) * angles_scores[minidex-1].first + 1./(1.+phi) * angles_scores.front().first;
                    newset.push_back(angles_scores.front()) ; 
                }
                else
                {
                    newset.push_back(angles_scores[minidex-1]) ;
                    currentAngle = phi/(1.+phi) * angles_scores[minidex-1].first + 1./(1.+phi) * angles_scores[minidex+1].first;
                    newset.push_back(angles_scores[minidex+1]) ; 
                }
                angles_scores = newset ;
//                 if(firstTension)
                    stiff->setAngle ( currentAngle ) ;     
                roughsampling = false ;
                change = true ;
                return ;
            }
           
           

            double minscore = angles_scores[0].second ;
            size_t minidex = 0 ;
            std::vector<std::pair<double, double>> newset ;
            for(size_t i = 1 ; i < angles_scores.size() ; i++)
            {
                if(angles_scores[i].second < minscore)
                {
                    minscore = angles_scores[i].second ;
                    minidex = i ;
                }
            }
            if(minidex == 0)
            {
                currentAngle = angles_scores.back().first-M_PI*.5 + angles_scores[minidex+1].first - angles_scores[minidex].first;
                newset.push_back(angles_scores.back()); 
                newset.push_back(angles_scores[minidex]); 
                newset.push_back(angles_scores[minidex+1]); 
               
            }
            else if(minidex == (angles_scores.size()-1))
            {
                currentAngle = angles_scores[minidex-1].first + angles_scores.front().first+M_PI*.5 - angles_scores[minidex].first;
                newset.push_back(angles_scores[minidex-1]);
                newset.push_back(angles_scores[minidex]); 
                newset.push_back(angles_scores.front()); 
                
            }
            else
            {
                currentAngle = angles_scores[minidex-1].first + angles_scores[minidex+1].first - angles_scores[minidex].first;
                newset.push_back(angles_scores[minidex-1]);
                newset.push_back(angles_scores[minidex]); 
                newset.push_back(angles_scores[minidex+1]); 
            }
            angles_scores = newset ;
            if(currentAngle > M_PI*.5)
                currentAngle -= M_PI*.5 ;
            if(currentAngle < 0)
                currentAngle += M_PI*.5 ;
            stiff->setAngle ( currentAngle ) ; 
            change = true ;
            if(iterationcount > iterationNumber/2)
            {

                initialAngle = currentAngle ;
                alternate = false ;
                alternateCheckpoint = true ;
            }    
        }
    }
    
}

std::pair< Vector, Vector > RotatingCrack::computeDamageIncrement ( ElementState &s )
{

    Vector range (1., 4 ) ;
    
    es = &s ;

    firstTension  = s.getParent()->getBehaviour()->getFractureCriterion()->directionInTension ( 0 ) ;
    secondTension = s.getParent()->getBehaviour()->getFractureCriterion()->directionInTension ( 1 ) ;


    if ( s.getParent()->getBehaviour()->getFractureCriterion()->directionMet ( 0 ) )
    {
        firstMet = true ;

            range[0] = getState() [0] ;

    }
    else
    {
        range[1] = getState() [1] ;
        range[0] = getState() [0] ;
        firstMet = false ;
    }

    if ( s.getParent()->getBehaviour()->getFractureCriterion()->directionMet ( 1 ) )
    {
        secondMet = true ;
//         if ( secondTension )
//         {
            range[3] = getState() [3] ;
//         }
//         else
//         {
//             range[2] = getState() [2] ;
//         }

    }
    else
    {
        range[3] = getState() [3] ;
        range[2] = getState() [2] ;
        secondMet = false ;
    }

    return std::make_pair ( getState(),  range ) ;
}

Matrix RotatingCrack::apply ( const Matrix &m, const Point & p , const IntegrableEntity * e , int g ) const
{
    return stiff->getTensor ( Point() ) *factor ;
}


void  RotatingCrack::computeDelta ( ElementState &s )
{

    Vector range ( 1., 4 ) ;
    delta = 1. ;
    return ;
    if ( s.getParent()->getBehaviour()->getFractureCriterion()->directionInTension ( 0 ) )
    {
		firstTension = true ;
        range[1] = getState() [1] ;
    }

    if ( s.getParent()->getBehaviour()->getFractureCriterion()->directionInCompression ( 0 ) )
    {
		firstTension = false ;
        range[0] = getState() [0] ;
    }

    if ( s.getParent()->getBehaviour()->getFractureCriterion()->directionInTension ( 1 ) )
    {
		secondTension = true ;
        range[3] = getState() [3] ;
    }

    if ( s.getParent()->getBehaviour()->getFractureCriterion()->directionInCompression ( 1 ) )
    {
		secondTension = false ;
        range[2] = getState() [2] ;
    }

//     range[1] = std::max ( range[1], range[3] ) ;
//     range[3] = range[1] ;

    delta = ( range-state ).max() ;
}


bool RotatingCrack::fractured() const
{
//     return state.max() >= .95 ;
// 	if ( fraction < 0 )
    return false ;

// 	return (firstTension && firstTensionFailure || !firstTension && firstCompressionFailure) || ( secondTension && secondTensionFailure || !secondTension && secondCompressionFailure ) ;

}

void addAndConsolidate ( std::vector<std::pair<double, double> > & target, std::vector<double> & weights, double a, double v, double tol = 1e-2 )
{
    for ( size_t i = 0 ; i < target.size() ; i++ )
    {
        if ( std::abs ( target[i].first - a ) < tol )
        {
            double newa = ( target[i].first * weights[i] + a ) / ( weights[i] + 1 ) ;
            double newv = target[i].second + v ;
            target[i] = std::make_pair ( newa, newv ) ;
            weights[i]++ ;
            return ;
        }
        if ( target[i].first > a )
        {
            target.insert ( target.begin() +i, std::make_pair ( a, v ) ) ;
            weights.insert ( weights.begin() +i, 1 ) ;
            return  ;
        }
    }
    target.push_back ( std::make_pair ( a, v ) ) ;
    weights.push_back ( 1 ) ;
}

void RotatingCrack::postProcess()
{

}

RotatingCrack::~RotatingCrack()
{
    delete stiff ;
}



FixedCrack::FixedCrack ( double E, double nu ) :  E ( E ), nu ( nu )
{
    getState ( true ).resize ( 4, 0. );
    getState( true )[0] = 0.998 ;
    isNull = false ;
    currentAngle = 0 ;
    factor = 1 ;
    es = nullptr ;
    firstTension = false ;
    secondTension = false ;
    firstTensionFailure = false ;
    secondTensionFailure = false ;
    firstCompressionFailure = false ;
    secondCompressionFailure = false ;
    angleset = false ;
}


int FixedCrack::getMode() const
{
    if ( es && es->getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet() &&
            ( (!firstTension && es->getParent()->getBehaviour()->getFractureCriterion()->directionInTension ( 0 ))
              || (firstTension && es->getParent()->getBehaviour()->getFractureCriterion()->directionInCompression ( 0 ))
              || (!secondTension && es->getParent()->getBehaviour()->getFractureCriterion()->directionInTension ( 1 ))
              || (secondTension && es->getParent()->getBehaviour()->getFractureCriterion()->directionInCompression ( 1 ) ))
       )
    {
        return 1 ;
    }
    return -1 ;
}

double FixedCrack::getAngleShift() const
{
    return 0 ;
}

std::pair< Vector, Vector > FixedCrack::computeDamageIncrement ( ElementState &s )
{
    Vector range ( 1., 4 ) ;

    if ( s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint() &&
            s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet() )
    {
        es = &s ;

        firstTension =  s.getParent()->getBehaviour()->getFractureCriterion()->directionInTension ( 0 ) ? true : false ;
        secondTension = s.getParent()->getBehaviour()->getFractureCriterion()->directionInTension ( 1 ) ? true : false ;


        if ( s.getParent()->getBehaviour()->getFractureCriterion()->directionMet ( 0 ) )
        {
            if ( firstTension  )
            {
                range[1] = getState() [1] ;
            }
            else if ( !firstTension )
            {
				range[0] = getState()[0] ;
            }
            else
            {
                range[1] = getState() [1] ;
                range[0] = getState() [0] ;
            }
        }
        else
        {
            range[1] = getState() [1] ;
            range[0] = getState() [0] ;
        }

        if ( s.getParent()->getBehaviour()->getFractureCriterion()->directionMet ( 1 ) )
        {
            if ( secondTension )
            {
                range[3] = getState() [3] ;
            }
            else if ( !secondTension )
            {
				range[2] = getState()[2] ;
            }
            else
            {
                range[1] = getState() [1] ;
                range[0] = getState() [0] ;
            }
        }
        else
        {
            range[3] = getState() [3] ;
            range[2] = getState() [2] ;
        }

    }
    else if ( s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint() )
    {
        firstTension =  s.getParent()->getBehaviour()->getFractureCriterion()->directionInTension ( 0 ) ? true : false ;
        secondTension = s.getParent()->getBehaviour()->getFractureCriterion()->directionInTension ( 1 ) ? true : false ;
    }

    return std::make_pair ( getState(),  range ) ;
}

Matrix FixedCrack::apply ( const Matrix & m, const Point & p, const IntegrableEntity * e, int g ) const
{

    if ( getState().max() < POINT_TOLERANCE )
    {
        return m ;
    }

// 	if(fractured())
// 		return m *0 ;

    double E_0 = E ;
    double E_1 = E ;
    double fs = getState() [0] ;
    double ss = getState() [2] ;
    if ( !firstTension )
    {
        fs = getState() [1] ;
        E_0 *= 0.5 ;
    }
    if ( !secondTension )
    {
        ss = getState() [3] ;
        E_1 *= 0.5 ;
    }

    E_0 *= ( 1. - fs ) ;
    E_1 *= ( 1. - ss ) ;

// 	double maxE = std::max(E_0, E_1) ;
// 	if(E_0 < E_1)
// 		E_0 = std::max(E_0, E_1*1e-4) ;
// 	if(E_1 < E_0)
// 		E_1 = std::max(E_1, E_0*1e-4) ;


    return OrthotropicStiffness ( E_0,
                                  E_1,
                                  E * ( 1.-std::max ( fs, ss ) ) * ( 1. - nu ) * .5,
                                  nu,
                                  currentAngle ).getTensor ( Point() ) *factor ;


}


void  FixedCrack::computeDelta ( ElementState &s )
{
    Vector range ( 1., 4 ) ;

    if ( s.getParent()->getBehaviour()->getFractureCriterion()->directionInTension ( 0 ) )
    {
        firstTension = true ;
        range[1] = getState() [1] ;
    }

    if ( s.getParent()->getBehaviour()->getFractureCriterion()->directionInCompression ( 0 ) )
    {
        firstTension = false ;
        range[0] = getState() [0] ;
    }

    if ( s.getParent()->getBehaviour()->getFractureCriterion()->directionInTension ( 1 ) )
    {
        secondTension = true ;
        range[3] = getState() [3] ;
    }

    if ( s.getParent()->getBehaviour()->getFractureCriterion()->directionInCompression ( 1 ) )
    {
        secondTension = false ;
        range[2] = getState() [2] ;
    }

    delta = ( range-state ).max() ;
}

bool FixedCrack::fractured() const
{
// 	if ( fraction < 0 )
    return false ;

    return ( (firstTension && firstTensionFailure) || (!firstTension && firstCompressionFailure) ) || ( (secondTension && secondTensionFailure) || (!secondTension && secondCompressionFailure) ) ;

}

void FixedCrack::postProcess()
{
    if ( converged )
    {
        getState ( true ) [0] = std::max ( getState ( true ) [0],getState ( true ) [2] ) ;
        getState ( true ) [2] = std::max ( getState ( true ) [0],getState ( true ) [2] ) ;
        getState ( true ) [1] = std::max ( getState ( true ) [1],getState ( true ) [3] ) ;
        getState ( true ) [3] = std::max ( getState ( true ) [1],getState ( true ) [3] ) ;
    }
    if ( converged && getState() [0] >= thresholdDamageDensity )
    {
        firstTensionFailure = true ;
        getState ( true ) [0] = 1. ;
    }
    if ( converged && getState() [1] >= thresholdDamageDensity )
    {
        firstCompressionFailure = true ;
        getState ( true ) [1] = 1. ;
    }
    if ( converged && getState() [2] >= thresholdDamageDensity )
    {
        secondTensionFailure = true ;
        getState ( true ) [2] = 1. ;
    }
    if ( converged && getState() [3] >= thresholdDamageDensity )
    {
        secondCompressionFailure = true ;
        getState ( true ) [3] = 1. ;
    }
}

FixedCrack::~FixedCrack()
{
}

DamageModel * RotatingCrack::getCopy() const
{
    RotatingCrack * ret = new RotatingCrack(E,nu) ;
    ret->factor = factor ;
    ret->copyEssentialParameters( this ) ;
    return ret ;
}

DamageModel * FixedCrack::getCopy() const
{
    FixedCrack * ret = new FixedCrack(E,nu) ;
    ret->factor = factor ;
    ret->copyEssentialParameters( this ) ;
    return ret ;
}

}
