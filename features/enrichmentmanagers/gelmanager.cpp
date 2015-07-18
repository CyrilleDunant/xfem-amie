#include "gelmanager.h"
#include "../../physics/materials/gel_behaviour.h"

namespace Amie
{

GelManager::GelManager(FeatureTree * ftree, double zonedensity, const std::vector<Feature *> & aggregates, double dr, double initialRadius) :deltaRadius(dr), reactedArea(0),aggregateArea(0), ftree(ftree)
{
    if(aggregates.empty())
        return ;

    int tries = 6000 ;
    double radiusFraction = 500 ;


    double xmin = aggregates[0]->getCenter().getX()-aggregates[0]->getRadius() ;
    double xmax = aggregates[0]->getCenter().getX()+aggregates[0]->getRadius() ;
    double ymin = aggregates[0]->getCenter().getY()-aggregates[0]->getRadius() ;
    double ymax = aggregates[0]->getCenter().getY()+aggregates[0]->getRadius() ;

    for(auto agg : aggregates)
    {
        xmin = std::min(xmin, agg->getCenter().getX()-agg->getRadius()) ;
        xmax = std::max(xmax, agg->getCenter().getX()+agg->getRadius()) ;
        ymin = std::min(ymin, agg->getCenter().getY()-agg->getRadius()) ;
        ymax = std::max(ymax, agg->getCenter().getY()+agg->getRadius()) ;
    }

    double area = (xmax-xmin)*(ymax-ymin) ;
    Rectangle baseGeometry(xmax-xmin, ymax-ymin, 0.5*(xmax+xmin), 0.5*(ymax+ymin)) ;
    int nzones = round(area*zonedensity) ;

    std::vector<ExpansiveZone *> zonesToPlace ;

    int trycount = 0 ;
    for( int i = 0 ; i < nzones && trycount < tries*nzones; i++ )
    {
        trycount++ ;
        Point pos( ( ( double )rand() / RAND_MAX - .5 ) * ( (xmax-xmin) - initialRadius * radiusFraction ), ( ( double )rand() / RAND_MAX - .5 ) * ( (ymax-ymin) - initialRadius * radiusFraction ) ) ;
        pos += baseGeometry.getCenter() ;
        bool alone  = true ;

        for( size_t j = 0 ; j < zonesToPlace.size() ; j++ )
        {
            if( dist( pos, zonesToPlace[j]->Circle::getCenter() ) < 0.0015 )
            {
                alone = false ;
                break ;
            }
        }

        if( alone )
            zonesToPlace.push_back( new ExpansiveZone( nullptr, initialRadius/15., pos.getX(), pos.getY(), new GelBehaviour() ) ) ;
        else
            i-- ;
    }

    std::map<Feature *, int> zonesPerIncs ;

    for( size_t i = 0 ; i < zonesToPlace.size() ; i++ )
    {
        bool placed = false ;

        for( size_t j = 0 ; j < aggregates.size() ; j++ )
        {
            if( dist( zonesToPlace[i]->getCenter(), aggregates[j]->getCenter() ) < aggregates[j]->getRadius() - 0.0002 /*&& incs[j]->getRadius() <= 0.008 && incs[j]->getRadius() > 0.004*/ && baseGeometry.in( zonesToPlace[i]->getCenter() ) )
            {
                zonesPerIncs[aggregates[j]]++ ; ;
                ftree->addFeature( aggregates[j], zonesToPlace[i] ) ;
                zones.push_back( std::make_pair( zonesToPlace[i], aggregates[j] ) ) ;
                reactedArea += zonesToPlace[i]->area() ;
                placed = true ;
                break ;
            }
        }

        if( !placed )
            delete zonesToPlace[i] ;
    }

    int count = 0 ;

    for( const auto & i : zonesPerIncs )
    {
        aggregateArea += i.first->area() ;
        count += i.second ;
//      std::cout << aggregateArea << "  " << count << std::endl ;
    }
    if(deltaRadius < 0)
        deltaRadius = sqrt( aggregateArea * 0.2 / ( ( double )zones.size() * M_PI ) ) / 800. ;
}

double GelManager::getReactedFraction() const
{
    return reactedArea/aggregateArea ;
}

double GelManager::getAggregateArea() const 
{
    return aggregateArea ;
}

bool GelManager::step(double dt, Vector * v, Mesh< DelaunayTriangle, DelaunayTreeItem >* dtree)
{
    if(dt < POINT_TOLERANCE)
        return false ;

    Feature *current = nullptr ;

    if( !zones.empty() )
        current = zones[0].second ;

    double currentArea = 0 ;
    int currentNumber = 0 ;
    int stoppedReaction = 0 ;
    reactedArea = 0 ;
    ftree->forceEnrichmentChange();

    for( size_t z = 0 ; z < zones.size() ; z++ )
    {
        zones[z].first->setRadius( zones[z].first->getRadius() + deltaRadius ) ;

//         std::cout << zones[z].first->getRadius() + deltaRadius << std::endl ;
        if( zones[z].second == current )
        {
            currentArea += zones[z].first->area() ;
            currentNumber++ ;
        }
        else
        {
            if( currentArea / zones[z - 1].second->area() > reactiveFraction )
            {
                stoppedReaction++ ;

                for( int m = 0 ; m < currentNumber ; m++ )
                {
                    reactedArea -= zones[z - 1 - m].first->area() ;
                    zones[z - 1 - m].first->setRadius( zones[z].first->getRadius() - deltaRadius ) ;
                    reactedArea += zones[z - 1 - m].first->area() ;
                }
            }

            currentArea = zones[z].first->area() ;
            currentNumber = 1 ;
            current = zones[z].second ;
        }

        reactedArea += zones[z].first->area() ;
    }
    
    return true ;
}

double GelManager::getReactiveFraction() const
{
    return reactiveFraction ;
}
void GelManager::setReactiveFraction(double r)
{
    reactiveFraction = r ;
}

bool GelManager::step(double dt, Vector * v, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * dtree)
{
    return true ;
}

void GelManager::setDeltaRadius(double dr)
{
    deltaRadius = dr ;
}
}
