#include "gelmanager.h"
#include "../../physics/materials/gel_behaviour.h"

namespace Amie
{

GelManager::GelManager(FeatureTree * f) : deltaRadius(0), reactedArea(0), aggregateArea(0), reactiveFraction(0), ftree(f)
{

}

GelManager::GelManager(FeatureTree * ftree, double zonedensity, const std::vector<Feature *> & aggregates,double reactiveFraction, double dr, double initialRadius) :deltaRadius(dr), reactedArea(0),aggregateArea(0), reactiveFraction(reactiveFraction), ftree(ftree)
{
    iterationCounter = 0 ;
    if(aggregates.empty())
        return ;
    double xmin = aggregates[0]->getCenter().getX()-aggregates[0]->getRadius() ;
    double xmax = aggregates[0]->getCenter().getX()+aggregates[0]->getRadius() ;
    double ymin = aggregates[0]->getCenter().getY()-aggregates[0]->getRadius() ;
    double ymax = aggregates[0]->getCenter().getY()+aggregates[0]->getRadius() ;

    for(auto agg : aggregates)
    {
        aggregateArea += agg->area() ; 
        xmin = std::min(xmin, agg->getCenter().getX()-agg->getRadius()) ;
        xmax = std::max(xmax, agg->getCenter().getX()+agg->getRadius()) ;
        ymin = std::min(ymin, agg->getCenter().getY()-agg->getRadius()) ;
        ymax = std::max(ymax, agg->getCenter().getY()+agg->getRadius()) ;
    }
    int tries = 6000 ;
    double area = (xmax-xmin)*(ymax-ymin) ;
    Rectangle baseGeometry(xmax-xmin, ymax-ymin, 0.5*(xmax+xmin), 0.5*(ymax+ymin)) ;
    int nzones = round(area*zonedensity) ;

    if(deltaRadius < 0)
    {
        deltaRadius = (sqrt((aggregateArea*reactiveFraction)/(nzones*M_PI))-initialRadius)/800. ;
    }
    std::vector<ExpansiveZone *> zonesToPlace ;

    int trycount = 0 ;
    for( int i = 0 ; i < nzones && trycount < tries*nzones; i++ )
    {
        trycount++ ;
        Point pos( ( ( double )rand() / RAND_MAX - .5 ) * ( xmax-xmin), ( ( double )rand() / RAND_MAX - .5 ) * ( ymax-ymin) ) ;
        pos += baseGeometry.getCenter() ;
        bool alone  = true ;

        for( size_t j = 0 ; j < zonesToPlace.size() ; j++ )
        {
            if( dist( pos, zonesToPlace[j]->Circle::getCenter() ) < deltaRadius*100. )
            {
                alone = false ;
                break ;
            }
        }

        if( alone )
            zonesToPlace.push_back( new ExpansiveZone( nullptr, initialRadius, pos.getX(), pos.getY(), new GelBehaviour() ) ) ;
        else
            i-- ;
    }

    std::map<Feature *, int> zonesPerIncs ;

    for( size_t i = 0 ; i < zonesToPlace.size() ; i++ )
    {
        bool placed = false ;

        for( size_t j = 0 ; j < aggregates.size() ; j++ )
        {
            if( dist( zonesToPlace[i]->getCenter(), aggregates[j]->getCenter() ) < aggregates[j]->getRadius() - deltaRadius*50. && baseGeometry.in( zonesToPlace[i]->getCenter() ) )
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

    aggregateArea = 0 ;
    for( const auto & i : zonesPerIncs )
    {
        aggregateArea += i.first->area() ;
        count += i.second ;
    }
    deltaRadius *= .5 ;
}

GelManager::GelManager( FeatureTree * f, InclusionFamily * inc, int gel, double rf, double dr)  : deltaRadius(dr), reactedArea(0), aggregateArea(0), reactiveFraction(rf), ftree(f)
{
    aggregateArea = 0 ;
    reactedArea = 0 ;

    std::vector<Feature *> found ;
    for(size_t i = 0 ; i < inc->features[gel].size() ; i++)
    {
        ExpansiveZone * test = dynamic_cast<ExpansiveZone *>( inc->features[gel][i] ) ;
        Feature * agg = inc->features[gel][i]->getFather() ;
        if( test != nullptr && agg != nullptr )
        {
            zones.push_back( std::make_pair( test, agg ) ) ;
            reactedArea += test->area() ;
            if(std::find( found.begin(), found.end(), agg ) == found.end())
            {
                found.push_back( agg ) ;
                aggregateArea += agg->area() ;
            }
        }
    }
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
//         if( zones[z].first->intersects( zones[z].second))
//             continue ;
        
        zones[z].first->setRadius( zones[z].first->getRadius() + deltaRadius ) ;

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

FunctionBasedGelManager::FunctionBasedGelManager( FeatureTree * f, InclusionFamily * inc, int gel, std::string function, double rf) : GelManager(f, inc, gel, rf, -1), radius(function)
{
    reactiveFraction = rf ;
    aggregateArea = 0 ;
    reactedArea = 0 ;
    deltaRadius = 0 ;

    std::vector<Feature *> found ;
    for(size_t i = 0 ; i < inc->features[gel].size() ; i++)
    {
        ExpansiveZone * test = dynamic_cast<ExpansiveZone *>( inc->features[gel][i] ) ;
        Feature * agg = inc->features[gel][i]->getFather() ;
        if( test != nullptr && agg != nullptr )
        {
            zones.push_back( std::make_pair( test, agg ) ) ;
            reactedArea += test->area() ;
            if(std::find( found.begin(), found.end(), agg ) == found.end())
            {
                found.push_back( agg ) ;
                aggregateArea += agg->area() ;
            }
        }
    }
}

bool FunctionBasedGelManager::step(double dt, Vector * v, Mesh< DelaunayTriangle, DelaunayTreeItem >* dtree)
{
    if(dt < POINT_TOLERANCE)
        return false ;

    setDeltaRadius( VirtualMachine().eval( radius, 0,0,0, ftree->getCurrentTime() ) - VirtualMachine().eval( radius, 0,0,0, std::max(0., ftree->getCurrentTime()-dt) ) ) ;

    return GelManager::step( dt, v, dtree ) ;
}

SpaceTimeGelManager::SpaceTimeGelManager( FeatureTree * f, InclusionFamily * inc, int gel, std::string r, double rf)  : GelManager(f, inc, gel, rf, -1)
{
    Function radius(r) ;

    aggregateArea = 0 ;
    reactedArea = 0 ;

    std::vector<Feature *> found ;
    for(size_t i = 0 ; i < inc->features[gel].size() ; i++)
    {
        GrowingExpansiveZone * test = dynamic_cast<GrowingExpansiveZone *>( inc->features[gel][i] ) ;
        Feature * agg = inc->features[gel][i]->getFather() ;
        if( test != nullptr && agg != nullptr )
        {
            stzones.push_back( std::make_pair( test, agg ) ) ;
            test->setRadiusFunction( radius ) ;
            if(std::find( found.begin(), found.end(), agg ) == found.end())
            {
                found.push_back( agg ) ;
                aggregateArea += agg->area() ;
            }
        }
    }
}

bool SpaceTimeGelManager::step(double dt, Vector * v, Mesh< DelaunayTriangle, DelaunayTreeItem >* dtree)
{
    if(dt < POINT_TOLERANCE)
        return false ;

    Feature *current = nullptr ;

    if( !stzones.empty() )
        current = stzones[0].second ;

    double currentArea = 0 ;
    int currentNumber = 0 ;
    int stoppedReaction = 0 ;
    reactedArea = 0 ;
    ftree->forceEnrichmentChange();

    std::vector<Circle> tests ;

    for( size_t z = 0 ; z < stzones.size() ; z++ )
    {
        tests.push_back( stzones[z].first->circleAtTime( Point(0,0,0, ftree->getCurrentTime() ) ) ) ;

        if( stzones[z].second == current )
        {
            currentArea += tests[z].area() ;
            currentNumber++ ;
        }
        else
        {
            if( currentArea / stzones[z - 1].second->area() > reactiveFraction )
            {
                stoppedReaction++ ;

                for( int m = 0 ; m < currentNumber ; m++ )
                {
                    Function finalRadius( tests[z - 1 - m].getRadius() ) ;
                    stzones[z - 1 - m].first->setRadiusFunction( finalRadius ) ;
                }
            }

            currentArea = tests[z].area() ;
            currentNumber = 1 ;
            current = stzones[z].second ;
        }

        reactedArea += tests[z].area() ;
    }
    
    return true ;
}

}
