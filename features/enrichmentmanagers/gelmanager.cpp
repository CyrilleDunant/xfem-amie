#include "gelmanager.h"
#include "../../physics/materials/gel_behaviour.h"

namespace Amie
{

GelManager::GelManager(FeatureTree * f) : deltaRadius(0), reactedArea(0), aggregateArea(0), reactiveFraction(0), ftree(f)
{

}

GelManager::GelManager(FeatureTree * ftree, double zonedensity, const std::vector<Feature *> & aggregates,double reactiveFraction) :deltaRadius(1e-5), reactedArea(0),aggregateArea(0), reactiveFraction(reactiveFraction), ftree(ftree)
{
    double initialRadius = 1e-14 ;
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
    int tries = 600000 ;
    double area = (xmax-xmin)*(ymax-ymin) ;
    Rectangle baseGeometry(xmax-xmin, ymax-ymin, 0.5*(xmax+xmin), 0.5*(ymax+ymin)) ;
    int nzones = round(area*zonedensity) ;

    if(deltaRadius < 0)
    {
        deltaRadius = (sqrt((aggregateArea*reactiveFraction)/(nzones*M_PI))-initialRadius)/800. ;
    }
    std::vector<ExpansiveZone *> zonesToPlace ;

    int trycount = 0 ;
    std::srand(6) ;
    for( int i = 0 ; i < nzones && trycount < tries*nzones; i++ )
    {
        trycount++ ;
        Point pos( ( ( double )std::rand() / RAND_MAX - .5 ) * ( xmax-xmin), ( ( double )std::rand() / RAND_MAX - .5 ) * ( ymax-ymin) ) ;
        pos += baseGeometry.getCenter() ;
        bool alone  = true ;

        for( size_t j = 0 ; j < zonesToPlace.size() ; j++ )
        {
            if( dist( pos, zonesToPlace[j]->Circle::getCenter() ) < deltaRadius*200. )
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
	    if( dist( zonesToPlace[i]->getCenter(), aggregates[j]->getCenter() ) < aggregates[j]->getRadius() && baseGeometry.in( zonesToPlace[i]->getCenter() ) && 
	        dist( zonesToPlace[i]->getCenter(), aggregates[j]->getCenter() ) > aggregates[j]->getRadius() - deltaRadius*50.
	    )
	    {
	      int count = 0 ;
	      Point vec = zonesToPlace[i]->getCenter()-aggregates[j]->getCenter() ;
	      vec /= vec.norm() ;
	      while(dist( zonesToPlace[i]->getCenter(), aggregates[j]->getCenter() ) > aggregates[j]->getRadius() - deltaRadius*50. && count++ < 50)
	      {
		 zonesToPlace[i]->getCenter() += vec*deltaRadius ;
	      }
	    }
            if( dist( zonesToPlace[i]->getCenter(), aggregates[j]->getCenter() ) < aggregates[j]->getRadius() - deltaRadius*50. && baseGeometry.in( zonesToPlace[i]->getCenter() ) )
            {
                zonesPerIncs[aggregates[j]]++ ;
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
    for( size_t z = 0 ; z < zones.size() ; z++ )
        zones[z].first->setRadius( deltaRadius*4 ) ;
}

GelManager::GelManager( FeatureTree * f, InclusionFamily * inc, double rf)  : deltaRadius(1e-5), reactedArea(0), aggregateArea(0), reactiveFraction(rf), ftree(f)
{
    aggregateArea = 0 ;
    reactedArea = 0 ;

    std::vector<Feature *> found ;
    for(size_t i = 0 ; i < inc->features.size() ; i++)
    {
        ExpansiveZone * test = dynamic_cast<ExpansiveZone *>( inc->features[i] ) ;
        Feature * agg = inc->features[i]->getFather() ;
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


double GelManager::areaForDr(double dr) const
{
    std::valarray<bool> mask(true, zones.size()) ;
    size_t iterator = 0;
    Feature *current = nullptr ;

    if( !zones.empty() )
        current = zones[0].second ;

    double currentArea = 0 ;
    int currentNumber = 0 ;
    for( size_t z = 0 ; z < zones.size() ; z++ )
    {
//         if( zones[z].first->intersects( zones[z].second))
//             continue ;
        
        zones[z].first->setRadius( zones[z].first->getRadius() + dr ) ;

        if( zones[z].second == current )
        {
            currentArea += (zones[z].first->getRadius() + dr)*(zones[z].first->getRadius() + dr)*M_PI ;
            currentNumber++ ;
            
        }
        else
        {
            if( currentArea / zones[z - 1].second->area() > reactiveFraction )
            {
                for( int m = 0 ; m < currentNumber ; m++ )
                {
                    mask[iterator+m] = false ;
                }
            }
            iterator+=currentNumber ;

            currentArea = (zones[z].first->getRadius() + dr)*(zones[z].first->getRadius() + dr)*M_PI ;
            currentNumber = 1 ;
            current = zones[z].second ;
        }
    }
    double ret = 0 ;
    for( size_t z = 0 ; z < zones.size() ; z++ )
        if(mask[z])
            ret += (zones[z].first->getRadius() + dr)*(zones[z].first->getRadius() + dr)*M_PI ;
    return ret ;
}

bool GelManager::step(double dt, Vector * v, Mesh< DelaunayTriangle, DelaunayTreeItem >* dtree)
{
    if(dt < POINT_TOLERANCE)
        return false ;
    
    if(steps == 0)
    {
        initialArea = 0 ;
        for( size_t z = 0 ; z < zones.size() ; z++ )
        {
            initialArea += zones[z].first->area() ;
        }
        deltaRadius = initialArea*initialArea/(totalSteps*totalSteps) ;
    }
    steps++ ;
    
    double targetArea = initialArea + aggregateArea*reactiveFraction*steps/totalSteps ;
//     std::cout << "target Area: " <<targetArea << " = "<< aggregateArea << "  "<< reactiveFraction << "  " << steps/totalSteps <<std::endl ;
    
    double minArea = 0 ;
    for( size_t z = 0 ; z < zones.size() ; z++ )
    {
        minArea += zones[z].first->area() ;
    }
    double mindr = 0 ;
    double maxdr = deltaRadius ;
    
    double maxArea = areaForDr(maxdr) ;

    
    do {
        maxdr *= 2. ;
        maxArea = areaForDr(maxdr) ;
//         std::cout << "max " <<maxArea << " vs " << targetArea << std::endl ;
    } while(maxArea < targetArea) ;
    double dr = 0 ;
    
    do {
        dr = (maxdr+mindr)*.5 ;
        double currentArea = areaForDr(dr) ;
//         std::cout << currentArea << " vs " << targetArea << std::endl ;
        
        if(currentArea > targetArea)
            maxdr = dr ;
        else
            mindr = dr ;
    } while(maxdr-mindr > 1e-32) ;
    deltaRadius = (maxdr+mindr)*.5 ;
//     std::cout << "increasing Delta: " << (maxdr+mindr)*.5 << std::endl ;
    
    
    Feature *current = nullptr ;

    if( !zones.empty() )
        current = zones[0].second ;

    double currentArea = 0 ;
    int currentNumber = 0 ;
    int stoppedReaction = 0 ;
    
    
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
                    zones[z - 1 - m].first->setRadius( std::max(zones[z - 1 - m].first->getRadius() - deltaRadius, sqrt((aggregateArea*reactiveFraction)/(zones.size()*totalSteps*M_PI)) ) ) ;
                }
            }

            currentArea = zones[z].first->area() ;
            currentNumber = 1 ;
            current = zones[z].second ;
        }
    }
    reactedArea = 0 ;
    for( size_t z = 0 ; z < zones.size() ; z++ )
        reactedArea += zones[z].first->area() ;
    
//     std::cout << "target area: " <<targetArea << "  resulting area: "<< reactedArea<< std::endl ;
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

FunctionBasedGelManager::FunctionBasedGelManager( FeatureTree * f, InclusionFamily * inc, std::string function, double rf) : GelManager(f, inc, rf), radius(function)
{
    reactiveFraction = rf ;
    aggregateArea = 0 ;
    reactedArea = 0 ;
    deltaRadius = 0 ;

    std::vector<Feature *> found ;
    for(size_t i = 0 ; i < inc->features.size() ; i++)
    {
        ExpansiveZone * test = dynamic_cast<ExpansiveZone *>( inc->features[i] ) ;
        Feature * agg = inc->features[i]->getFather() ;
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

SpaceTimeGelManager::SpaceTimeGelManager( FeatureTree * f, InclusionFamily * inc, std::string r, double rf)  : GelManager(f, inc, rf)
{
    Function radius(r) ;

    aggregateArea = 0 ;
    reactedArea = 0 ;

    std::vector<Feature *> found ;
    for(size_t i = 0 ; i < inc->features.size() ; i++)
    {
        GrowingExpansiveZone * test = dynamic_cast<GrowingExpansiveZone *>( inc->features[i] ) ;
        Feature * agg = inc->features[i]->getFather() ;
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
