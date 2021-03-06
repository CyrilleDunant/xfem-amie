// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2006-2011
//
// Copyright: See COPYING file that comes with this distribution
//

#include "crack.h"
#include "../physics/fracturecriteria/fracturecriterion.h"
#include "../elements/integrable_entity.h"
#include "../polynomial/vm_function_extra.h"

#include <fstream>

namespace Amie {

BranchedCrack::BranchedCrack(Feature *father, Point * a, Point * b) : EnrichmentFeature(father), SegmentedLine(std::valarray<Point * >(2))
{
	boundingPoints[0] = a ;
	boundingPoints[1] = b ;
	this->SegmentedLine::center = *a ;
	branches.push_back ( this ) ;
	scorePropagation = false ;
	
	Segment seg1(*a, *b) ;
	Segment seg0(*b, *a) ;
	double angle1 = seg1.normal().angle() ;
	double angle0 = seg0.normal().angle() ;
	
	if(father->in(*a))
		tips.push_back(std::make_pair(boundingPoints[0], angle0)) ;
	if(father->in(*b))
		tips.push_back(std::make_pair(boundingPoints[1], angle1)) ;

	changed = false ;
	isUpdated = false ;

	enrichementRadius = 0.00075 ;
	
}

bool operator ==(const std::pair<Amie::Point*, double> & a, const Amie::Point* b)
{
	return a.first == b ;
}

BranchedCrack::BranchedCrack(Point * a, Point * b) : EnrichmentFeature(nullptr), SegmentedLine(std::valarray<Point * >(2))
{
	boundingPoints[0] = a ;
	boundingPoints[1] = b ;
	this->SegmentedLine::center = *a ;
	branches.push_back ( this ) ;
	scorePropagation = false ;
	
	Segment seg1(*a, *b) ;
	Segment seg0(*b, *a) ;
	
	double angle1 = seg1.normal().angle() ;
	double angle0 = seg0.normal().angle() ;
	
	tips.push_back(std::make_pair(boundingPoints[0], angle0)) ;
	tips.push_back(std::make_pair(boundingPoints[1], angle1)) ;

	changed = false ;
	isUpdated = false ;

	enrichementRadius = 0.00075 ;
}

void BranchedCrack::branch(Point* fromTip, const std::vector<Point *> & newTip)
{
	for(auto i = tips.begin() ; i != tips.end() ; ++i)
	{
		if(i->first == fromTip)
		{
			tips.erase(i) ;
			break ;
		}
	}
	
	std::map<double, Point *> sortedNewTips ;
	for(std::vector<Point * >::const_iterator j = newTip.begin() ; j !=newTip.end() ; ++j)
	{
		double angle = atan2((*j)->getY()-fromTip->getY(), (*j)->getX()-fromTip->getX()) ;
		tips.push_back ( std::make_pair(*j, angle) ) ;
		sortedNewTips[angle] = *j ;
	}
	
	std::vector<Point *> newTipSorted ;
	for(std::map<double, Point *>::const_iterator j = sortedNewTips.begin() ; j !=sortedNewTips.end() ; ++j)
		newTipSorted.push_back(j->second) ;
	
	SegmentedLine * branchToExtend = nullptr ;
	bool fromHead = false ;
	for ( size_t i = 0 ; i < branches.size() ; i++ )
	{
		if ( branches[i]->getHead() == fromTip )
		{
			branchToExtend = branches[i] ;
			fromHead = true ;
			break ;
		}

		if ( branches[i]->getTail() == fromTip )
		{
			branchToExtend = branches[i] ;
			break ;
		}
	}
	
	for(size_t i = 1 ; i < newTipSorted.size() ; ++i)
	{
		std::valarray<Point * > newBranch ( 2 ) ;
		newBranch[0] = fromTip ;
		newBranch[1] = newTipSorted[i] ;
		branches.push_back ( new SegmentedLine ( newBranch ) ) ;
	}
	
	for(size_t i = 0 ; i < newTipSorted.size()-1 ; ++i)
	{
		std::valarray<Point *> forkComponent0 ( 3 ) ;
		forkComponent0[0] = newTipSorted[i] ;
		forkComponent0[1] = fromTip ;
		forkComponent0[2] = newTipSorted[i+1] ;
		forks.push_back ( new SegmentedLine ( forkComponent0 ) ) ;
	}
	


	std::valarray<Point *> forkComponent1 ( 3 ) ;
	std::valarray<Point *> forkComponent2 ( 3 ) ;

	std::valarray<Point *> newBP( branchToExtend->getBoundingPoints().size() +1 ) ;

	if ( fromHead )
	{
		std::copy ( &branchToExtend->getBoundingPoints() [0],
					&branchToExtend->getBoundingPoints() [newBP.size()-1],
					&newBP[1] ) ;
		newBP[0] = newTipSorted.front() ;
		branchToExtend->setBoundingPoints(newBP) ;
		forkComponent1[0] = newBP[0] ;
		forkComponent1[1] = newBP[1] ;
		forkComponent1[2] = newBP[2] ;
		forkComponent2[0] = newTipSorted.back() ;
		forkComponent2[1] = newBP[1] ;
		forkComponent2[2] = newBP[2] ;

	}
	else
	{
		std::copy ( &branchToExtend->getBoundingPoints() [0],
					&branchToExtend->getBoundingPoints() [newBP.size()-1],
					&newBP[0] ) ;
		newBP[newBP.size()-1] = newTipSorted.front() ;
		branchToExtend->setBoundingPoints(newBP) ;
		forkComponent1[0] = newBP[newBP.size()-3] ;
		forkComponent1[1] = newBP[newBP.size()-2] ;
		forkComponent1[2] = newBP[newBP.size()-1] ;
		forkComponent2[0] = newBP[newBP.size()-3] ;
		forkComponent2[1] = newBP[newBP.size()-2] ;
		forkComponent2[2] = newTipSorted.back() ;
	}

	forks.push_back ( new SegmentedLine ( forkComponent1 ) ) ;
	forks.push_back ( new SegmentedLine ( forkComponent2 ) ) ;
}

void BranchedCrack::branch ( Point* fromTip, Point * newTip0, Point * newTip1 )
{
	for(auto i = tips.begin() ; i !=tips.end() ; ++i)
	{
		if(i->first == fromTip)
		{
			tips.erase(i) ;
			break ;
		}
	}
	tips.push_back ( std::make_pair(newTip0, atan2(newTip0->getY()-fromTip->getY(), newTip0->getX()-fromTip->getX())) ) ;
	tips.push_back ( std::make_pair(newTip1, atan2(newTip1->getY()-fromTip->getY(), newTip1->getX()-fromTip->getX())) ) ;
	std::valarray<Point * > newBranch ( 2 ) ;
	newBranch[0] = fromTip ;
	newBranch[1] = newTip0 ;

	SegmentedLine * branchToExtend = nullptr ;
	bool fromHead = false ;
	for ( size_t i = 0 ; i < branches.size() ; i++ )
	{
		if ( branches[i]->getHead() == fromTip )
		{
			branchToExtend = branches[i] ;
			fromHead = true ;
			break ;
		}

		if ( branches[i]->getTail() == fromTip )
		{
			branchToExtend = branches[i] ;
			break ;
		}
	}

	branches.push_back ( new SegmentedLine ( newBranch ) ) ;
	std::valarray<Point *> forkComponent0 ( 3 ) ;
	forkComponent0[0] = newTip0 ;
	forkComponent0[1] = fromTip ;
	forkComponent0[2] = newTip1 ;

	std::valarray<Point *> forkComponent1 ( 3 ) ;
	std::valarray<Point *> forkComponent2 ( 3 ) ;

	std::valarray<Point *> newBP ( branchToExtend->getBoundingPoints().size() +1 ) ;
	if ( fromHead )
	{
		std::copy ( &branchToExtend->getBoundingPoints() [0],
		            &branchToExtend->getBoundingPoints() [newBP.size()-1],
		            &newBP[1] ) ;
		newBP[0] = newTip0 ;
		branchToExtend->setBoundingPoints(newBP) ;
		forkComponent1[0] = newBP[0] ;
		forkComponent1[1] = newBP[1] ;
		forkComponent1[2] = newBP[2] ;
		forkComponent2[0] = newTip1 ;
		forkComponent2[1] = newBP[1] ;
		forkComponent2[2] = newBP[2] ;

	}
	else
	{
		std::copy ( &branchToExtend->getBoundingPoints() [0],
		            &branchToExtend->getBoundingPoints() [newBP.size()-1],
		            &newBP[0] ) ;
		newBP[newBP.size()-1] = newTip0 ;
		branchToExtend->setBoundingPoints(newBP) ;
		forkComponent1[0] = newBP[newBP.size()-3] ;
		forkComponent1[1] = newBP[newBP.size()-2] ;
		forkComponent1[2] = newBP[newBP.size()-1] ;
		forkComponent2[0] = newBP[newBP.size()-3] ;
		forkComponent2[1] = newBP[newBP.size()-2] ;
		forkComponent2[2] = newTip1 ;
	}

	forks.push_back ( new SegmentedLine ( forkComponent0 ) ) ;
	forks.push_back ( new SegmentedLine ( forkComponent1 ) ) ;
	forks.push_back ( new SegmentedLine ( forkComponent2 ) ) ;

}

double BranchedCrack::propagationAngleFromTip(const std::pair<Point *, double> & tip, Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree)
{
	double acount = 0 ;
	double aangle = 0 ;
		
	Point lastDir( sin(tip.second) , cos(tip.second)) ;
	
	Circle c(enrichementRadius, *(tip.first)) ;
	std::vector<DelaunayTriangle*> disk = dtree->getConflictingElements (&c) ;
	std::vector<DelaunayTriangle*> tris ;
	while( disk.size() < 16)
	{
		c.setRadius(c.getRadius()*2.) ;
		disk = dtree->getConflictingElements ( &c ) ;
	}


	for ( size_t i = 0 ; i < disk.size() ; i++ )
	{
		
		if ( !disk[i]->in ( *getTail() )
		     && (!Line(*tip.first, lastDir).intersects(disk[i]->getPrimitive())
				|| c.in(disk[i]->getCenter()))
			)
		{
			tris.push_back ( disk[i] ) ;
		}
	}
	
	if ( !tris.empty() )
	{
		
		for ( size_t i = 0 ; i < tris.size() ; i++ )
		{
			
			DelaunayTriangle * current = tris[i] ;
			Vector a(0.,1) ;
			Vector principalStresses(0., 2) ;
			
			for ( size_t j = 0 ; j < current->getBoundingPoints().size() ; j++ )
			{
				
				Point currentDir ( current->getBoundingPoint ( j ).getX()-tip.first->getX(),
				                   current->getBoundingPoint ( j ).getY()-tip.first->getY() ) ;
				
				if ( ( currentDir*lastDir ) > 0 )
				{
					current->getState().getField( PRINCIPAL_REAL_STRESS_FIELD, current->getBoundingPoint(j), principalStresses, false) ;
					
					if ( currentDir.norm() > 1e-8 )
					{
						current->getState().getField( PRINCIPAL_STRESS_ANGLE_FIELD, current->getBoundingPoint(j), a, false) ;
						aangle += a[0]*2. ;
						acount++ ;
					}
				}
			}
		}
	}
	
	if ( acount )
		aangle /= acount ;

	return aangle ;
}

std::pair<double, double> BranchedCrack::computeJIntegralAtTip ( std::pair<Point *, double> & tip, Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree )
{
	Point direction ( cos(tip.second), sin(tip.second)) ;
	Segment tipSegment ( *(tip.first) , direction) ;
	Circle c ( enrichementRadius*.5, *(tip.first) ) ;
	
	
	std::vector<DelaunayTriangle *> disk = dtree->getConflictingElements ( &c ) ;
	std::vector<DelaunayTriangle *> ring ;
	if(disk.size() == 0)
		return std::make_pair(0., 0.) ;
	
	if(disk.size() == 1)
	{
		influenceRadius = disk[0]->getRadius() ;
		c.setRadius(disk[0]->getRadius()) ;
		disk = dtree->getConflictingElements ( &c ) ;
	}
	
	std::vector<std::pair<Segment *, DelaunayTriangle *> > gamma ;
	
	
	for ( size_t j = 0 ; j < disk.size() ; j++ )
	{
		Triangle * t = dynamic_cast<Triangle *> ( disk[j] ) ;
		if ( t->intersects ( &c ) )
		{
			ring.push_back ( disk[j] ) ;
			
			Segment A ( *disk[j]->first, *disk[j]->second ) ;
			Segment B ( *disk[j]->third, *disk[j]->first ) ;
			Segment C ( *disk[j]->second, *disk[j]->third ) ;
			
			if ( tipSegment.intersects ( t ) )
			{
				
				if ( A.intersects ( &c ) && B.intersects ( &c ) && !c.in ( C.midPoint() ) )
				{
					Point ipoint = C.intersection ( tipSegment ) ;
					gamma.push_back ( std::make_pair ( new Segment ( *disk[j]->third, ipoint ),disk[j] ) ) ;
					gamma.push_back ( std::make_pair ( new Segment ( *disk[j]->second, ipoint ),disk[j] ) ) ;
				}
				if ( B.intersects ( &c ) && C.intersects ( &c ) && !c.in ( A.midPoint() ) )
				{
					Point ipoint = A.intersection ( tipSegment ) ;
					gamma.push_back ( std::make_pair ( new Segment ( *disk[j]->first, ipoint ),disk[j] ) ) ;
					gamma.push_back ( std::make_pair ( new Segment ( *disk[j]->second, ipoint ),disk[j] ) ) ;
				}
				if ( C.intersects ( &c ) && A.intersects ( &c ) && !c.in ( B.midPoint() ) )
				{
					Point ipoint = B.intersection ( tipSegment ) ;
					gamma.push_back ( std::make_pair ( new Segment ( *disk[j]->third, ipoint ), disk[j] ) ) ;
					gamma.push_back ( std::make_pair ( new Segment ( *disk[j]->first, ipoint ),disk[j] ) ) ;
				}
			}
			else
			{
				if ( A.intersects ( &c ) && B.intersects ( &c ) && !c.in ( C.midPoint() ) )
					gamma.push_back ( std::make_pair ( new Segment ( C.first(), C.second() ), disk[j] ) ) ;
				if ( B.intersects ( &c ) && C.intersects ( &c ) && !c.in ( A.midPoint() ) )
					gamma.push_back ( std::make_pair ( new Segment ( A.first(), A.second() ), disk[j] ) ) ;
				if ( C.intersects ( &c ) && A.intersects ( &c ) && !c.in ( B.midPoint() ) )
					gamma.push_back ( std::make_pair ( new Segment ( B.first(), B.second() ), disk[j] ) ) ;
			}
		}
	}
	
	double freeEnergy0 = 0 ;
	double freeEnergy1 = 0 ;
	
	for ( size_t j = 0 ; j < gamma.size() ; j++ )
	{
		if ( gamma[j].second->getBehaviour()->type != VOID_BEHAVIOUR )
		{
			// 		FunctionMatrix displacements = gamma[j].second->getState().getDisplacementFunction() ;
			std::valarray<std::pair<Point, double> > gaussPointsi = gamma[j].first->getGaussPoints() ;
			std::valarray<std::pair<Point, double> > gaussPoints(3) ;
			gaussPoints[0] = gaussPointsi[0] ; 
			gaussPoints[1] = gaussPointsi[1] ; 
			gaussPoints[2] = std::make_pair ( gamma[j].first->midPoint(), ( double ) 1 ) ;
			Point n = gamma[j].first->normal ( gamma[j].second->getCenter() ) ;
			
			
			Vector stepLengthal ( 2 ) ; stepLengthal[0] = n.getX() ; stepLengthal[1] = n.getY() ;
			Vector d ( 2 ) ; d[0] = -n.getY() ;   d[1] = n.getX() ;
			double ilocal0 = 0 ;
			double ilocal1 = 0 ;
			
			Point localVector = gamma[j].second->inLocalCoordinates ( tipSegment.vector() ) ;
			Point localNormal ( -localVector.getY(), localVector.getX() ) ;
			Vector stress(0., 3) ;
			Vector strain(0., 3) ;
			
			for ( size_t k = 0 ; k < gaussPoints.size() ; k++ )
			{
				Point localGaussPoint = gamma[j].second->inLocalCoordinates ( gaussPoints[k].first ) ;
				gamma[j].second->getState().getField(TOTAL_STRAIN_FIELD, REAL_STRESS_FIELD, localGaussPoint, strain, stress, true) ;
				Matrix sigma = makeStressOrStrainMatrix(stress) ;
				Matrix epsilon = makeStressOrStrainMatrix(strain) ;
				double sigma_epsilon = sigma[0][0]*epsilon[0][0] + sigma[0][1]*epsilon[0][1] + sigma[1][0]*epsilon[1][0] + sigma[1][1]*epsilon[1][1];//0.5*std::inner_product(&sigma.array()[0], &sigma.array()[sigma.size()], &epsilon.array()[0], 0. );
				
				Vector T = sigma * stepLengthal ;
				
				
				ilocal0 += ( 0.5*sigma_epsilon*stepLengthal[0] - ( T[0]*epsilon[0][0] + T[1]*epsilon[0][1] ) ) *gaussPoints[k].second ;
				ilocal1 += ( 0.5*sigma_epsilon*stepLengthal[1] - ( T[0]*epsilon[0][1] + T[1]*epsilon[1][1] ) ) *gaussPoints[k].second ;
				
			}
			freeEnergy0 += ilocal0*gamma[j].first->norm() ;
			freeEnergy1 += ilocal1*gamma[j].first->norm() ;
		}
	}
	return std::make_pair ( freeEnergy0, freeEnergy1 ) ;
}

void BranchedCrack::grow( Point* fromTip, Point* newTip)
{
	for(auto i = tips.begin() ; i !=tips.end() ; ++i)
	{
		if(i->first == fromTip)
		{
			tips.erase(i) ;
			break ;
		}
	}

	SegmentedLine * branchToExtend = nullptr ;
	bool fromHead = false ;
	for ( size_t i = 0 ; i < branches.size() ; i++ )
	{
		if ( branches[i]->getHead() == fromTip )
		{
			branchToExtend = branches[i] ;
			fromHead = true ;
			break ;
		}
		
		if ( branches[i]->getTail() == fromTip )
		{
			branchToExtend = branches[i] ;
			break ;
		}
	}

	tips.push_back ( std::make_pair(newTip, Segment(*fromTip,*newTip).normal().angle()) ) ;
	
	std::valarray<Point *> newBP ( branchToExtend->getBoundingPoints().size() +1 ) ;
	if ( fromHead )
	{
		std::copy ( &branchToExtend->getBoundingPoints() [0],
		            &branchToExtend->getBoundingPoints() [newBP.size()-1],
		            &newBP[1] ) ;
		newBP[0] = newTip ;
		
	}
	else
	{
		std::copy ( &branchToExtend->getBoundingPoints() [0],
		            &branchToExtend->getBoundingPoints() [newBP.size()-1],
		            &newBP[0] ) ;
		newBP[newBP.size()-1] = newTip ;
	}
	branchToExtend->setBoundingPoints(newBP) ;
}

void BranchedCrack::move(Point * toMove, const Point & target) 
{    
    for(auto i = tips.begin() ; i !=tips.end() ; ++i)
    {
        if(i->first == toMove)
        {
            for ( size_t j = 0 ; j < branches.size() ; j++ )
            {
                if(branches[j]->getBoundingPoints()[ branches[j]->getBoundingPoints().size()-1] == toMove)
                {
                    i->second = Segment(*branches[j]->getBoundingPoints()[ branches[j]->getBoundingPoints().size()-2],target).normal().angle() ;
                    *(i->first) = target ;
                    return ;
                }
                if(branches[j]->getBoundingPoints()[ 0 ] == toMove)
                {
                    i->second = Segment(*branches[j]->getBoundingPoints()[1],target).normal().angle() ;
                    *(i->first) = target ;
                    return ;
                }
            }
        }
    }
    
    for ( size_t j = 0 ; j < branches.size() ; j++ )
    {
        if(branches[j]->getBoundingPoints()[ branches[j]->getBoundingPoints().size()-2] == toMove)
        {
            for(auto i = tips.begin() ; i !=tips.end() ; ++i)
            {
                if(i->first == branches[j]->getBoundingPoints()[ branches[j]->getBoundingPoints().size()-1])
                {
                    i->second = Segment(target, *branches[j]->getBoundingPoints()[ branches[j]->getBoundingPoints().size()-1]).normal().angle() ;
                    *branches[j]->getBoundingPoints()[ branches[j]->getBoundingPoints().size()-2] = target ;
                    return ;
                }
            }
        }
        if(branches[j]->getBoundingPoints()[ 1 ] == toMove)
        {
            for(auto i = tips.begin() ; i !=tips.end() ; ++i)
            {
                if(i->first == branches[j]->getBoundingPoints()[ 0])
                {
                    i->second = Segment(target, *branches[j]->getBoundingPoints()[0]).normal().angle() ;
                    *branches[j]->getBoundingPoints()[ 1 ] = target ;
                    return ;
                }
            }
        }
    }
    
    for ( size_t j = 0 ; j < branches.size() ; j++ )
    {
     for(size_t i = 0 ; i < branches[j]->getBoundingPoints().size() ; i++)
     {
         if(toMove == branches[j]->getBoundingPoints()[i])
             *branches[j]->getBoundingPoints()[i] = target ;
     }
    }
    
    
}

void BranchedCrack::merge ( BranchedCrack & newSet)
{
	//first, we find the intersection point.
	if(newSet.getTips().empty())
	{
		std::cerr << "Refusing merge with tipless crack" << std::endl ;
		return ;
	}
	
	std::pair<Point *, double> tipForMerge = newSet.getTips()[0] ;
	double dist = 0 ;
	for(size_t i = 0 ; i < tips.size() ; i++)
	{
		dist = std::max(squareDist(tips[i].first, tipForMerge.first), dist) ;
	}
	double maxdist = dist*2. ;
	std::pair<Point*, Point*> targetSegment ;
	SegmentedLine * toBranch = nullptr;
	
	
	for(size_t i = 0 ; i < newSet.getTips().size() ; i++)
	{
		std::pair<Point *, double> tipForMergePotential = newSet.getTips()[i] ;
		double distPotential = maxdist ;
		std::pair<Point*, Point*> targetSegmentPotential ;
		SegmentedLine * toBranchPotential = nullptr;
		for(size_t j = 0 ; j < branches.size() ; j++)
		{
			for(size_t k = 0 ; k < branches[j]->getBoundingPoints().size()-1 ; k++)
			{
				Point mid = (branches[j]->getBoundingPoint(k)
				             +branches[j]->getBoundingPoint(k+1))*.5 ;
				double dtest = squareDist(mid, *tipForMergePotential.first) ;
				
				if(dtest < distPotential)
				{
					distPotential = dtest ;
					targetSegmentPotential = std::make_pair(&branches[j]->getBoundingPoint(k),
					                          &branches[j]->getBoundingPoint(k+1)) ;
					
					toBranchPotential = branches[j] ;
				}
				
			}
			
		}
		
		if(distPotential< dist)
		{
			dist = distPotential ;
			tipForMerge = tipForMergePotential ;
			targetSegment = targetSegmentPotential ;
			toBranch = toBranchPotential ;
		}
	}
	
	//a new merge point is instantiated which is the intersection between 
	// the nearest segment and the prolonged branch
	SegmentedLine * fromBranch = nullptr ;
	bool isTail = false ;
	
	for(size_t i = 0 ; i < newSet.getBranches().size() ; i++)
	{
		if(newSet.getBranches()[i]->getHead() == tipForMerge.first)
		{
			fromBranch = newSet.getBranches()[i] ;
			break ;
		}
		
		if(newSet.getBranches()[i]->getTail() == tipForMerge.first)
		{
			fromBranch = newSet.getBranches()[i] ;
			isTail = true ;
			break ;
		}
	}
	
	Point intersection ;
	
	if(isTail)
	{
		Line to(*targetSegment.first, *targetSegment.first-*targetSegment.second) ;

		intersection = to.projection(fromBranch->getBoundingPoint(fromBranch->getBoundingPoints().size()-1));
		
		if(!Segment(*targetSegment.first, *targetSegment.second).on(intersection))
			intersection = (*targetSegment.first-*targetSegment.second)*.5 ;
	}
	else
	{
		Line to(*targetSegment.first, *targetSegment.first-*targetSegment.second) ;
		intersection = to.projection( fromBranch->getBoundingPoint(0)) ;
		
		if(!Segment(*targetSegment.first, *targetSegment.second).on(intersection))
			intersection = (*targetSegment.first-*targetSegment.second)*.5 ;
	}
	
	//then 1) the tip is removed as a tip, 
	newSet.getTips().erase(std::find(newSet.getTips().begin(), newSet.getTips().end(), tipForMerge)) ;
	
	//     2) the corresponding branch is extended
	
	Point * intersectionPoint ;
	if(isTail)
	{
		std::valarray<Point *> points = fromBranch->getBoundingPoints() ;
		fromBranch->getBoundingPoints().resize(fromBranch->getBoundingPoints().size()+1) ;
		std::copy(&points[0], &points[points.size()], &fromBranch->getBoundingPoints()[0]) ;
		fromBranch->getBoundingPoints()[points.size()] = new Point(intersection) ;
		intersectionPoint = fromBranch->getBoundingPoints()[points.size()] ;
	}
	else
	{
		std::valarray<Point *> points = fromBranch->getBoundingPoints() ;
		fromBranch->getBoundingPoints().resize(fromBranch->getBoundingPoints().size()+1) ;
		std::copy(&points[0], &points[points.size()], &fromBranch->getBoundingPoints()[1]) ;
		fromBranch->getBoundingPoints()[0] = new Point(intersection) ;
		intersectionPoint = fromBranch->getBoundingPoints()[0] ;
	}
	
	//     3) the original branch is cut
	for(size_t k = 0 ; k < toBranch->getBoundingPoints().size()-1 ; k++)
	{
// 		toBranch->getBoundingPoint(k).print() ;
// 		toBranch->getBoundingPoint(k+1).print() ;
// 		targetSegment.first->print() ;
// 		targetSegment.second->print() ;
		if(toBranch->getBoundingPoint(k) == *(targetSegment.first) && toBranch->getBoundingPoint(k+1) == *(targetSegment.second))
		{
			std::valarray<Point *> points = toBranch->getBoundingPoints() ;
			toBranch->getBoundingPoints().resize(points.size()+1) ;
			std::copy(&points[0], &points[k+1], &toBranch->getBoundingPoints()[0]) ;
			
			if(isTail)
				toBranch->getBoundingPoints()[k+1] = fromBranch->getTail() ;
			else
				toBranch->getBoundingPoints()[k+1] = fromBranch->getHead() ;
			
			std::copy(&points[k+1], &points[points.size()], &toBranch->getBoundingPoints()[k+2]) ;

			break ;
		}
	}
	
	//     4) a new fork is added
	std::valarray<Point*> b0(3) ;
	b0[0] = targetSegment.first ;
	b0[1] = intersectionPoint ;
	b0[2] = targetSegment.second ;
	
	std::valarray<Point*> b1(3) ;
	b1[0] = targetSegment.first ;
	b1[1] = intersectionPoint ;
	b1[2] = tipForMerge.first ;
	
	std::valarray<Point*> b2(3) ;
	b2[0] = targetSegment.second ;
	b2[1] = intersectionPoint ;
	b2[2] = tipForMerge.first ;
	forks.push_back ( new SegmentedLine ( b0 ) ) ;
	forks.push_back ( new SegmentedLine ( b1 ) ) ;
	forks.push_back ( new SegmentedLine ( b2 ) ) ;
	
	//     5) the elements are migrated
	for(size_t i =  0 ; i < newSet.getTips().size() ; i++)
		tips.push_back(newSet.getTips()[i]) ;
	for(size_t i =  0 ; i < newSet.getBranches().size() ; i++)
		branches.push_back(newSet.getBranches()[i]) ;
	for(size_t i =  0 ; i < newSet.getForks().size() ; i++)
		forks.push_back(newSet.getForks()[i]) ;
	
	//     6) cleanup
	newSet.getTips().clear() ;
	newSet.getBranches().clear() ;
	newSet.getForks().clear() ;
}

const std::vector<SegmentedLine *> & BranchedCrack::getBranches() const
{
	return branches ;
}

void BranchedCrack::enrichTips(size_t & startid, Mesh<DelaunayTriangle,DelaunayTreeItem> * dt)
{
	for(size_t i =  0 ; i < tips.size(); i++)
	{
		enrichTip(startid, dt, tips[i]) ;
	}
}

Function getBlendingFunction(const std::map<const Point *, int> & dofIds, const DelaunayTriangle * t)
{
// 	return Function("1") ;

// if(t->getOrder() == QUADRATIC)
// {
// 	TriElement father(QUADRATIC) ;
// 	if(dofIds.find(t->first) != dofIds.end() && dofIds.find(t->second) == dofIds.end() && dofIds.find(t->third) == dofIds.end())
// 	{
// 		return father.getShapeFunction(0) + 0.25*father.getShapeFunction(1)+ 0.25*father.getShapeFunction(5);
// 	}
// 	
// 	if(dofIds.find(t->first) == dofIds.end() && dofIds.find(t->second) != dofIds.end() && dofIds.find(t->third) == dofIds.end())
// 	{
// 		return father.getShapeFunction(2) + 0.25*father.getShapeFunction(1)+ 0.25*father.getShapeFunction(3);
// 	}
// 	
// 	if(dofIds.find(t->first) == dofIds.end() && dofIds.find(t->second) == dofIds.end() && dofIds.find(t->third) != dofIds.end())
// 	{
// 		return father.getShapeFunction(4) + 0.25*father.getShapeFunction(3)+ 0.25*father.getShapeFunction(5);
// 	}
// 	
// 	if(dofIds.find(t->first) == dofIds.end() && dofIds.find(t->second) != dofIds.end() && dofIds.find(t->third) != dofIds.end())
// 	{
// 		return father.getShapeFunction(2)+father.getShapeFunction(3)+father.getShapeFunction(4) + 0.25*father.getShapeFunction(1)+ 0.25*father.getShapeFunction(5);
// 	}
// 	
// 	if(dofIds.find(t->first) != dofIds.end() && dofIds.find(t->second) == dofIds.end() && dofIds.find(t->third) != dofIds.end())
// 	{
// 		return father.getShapeFunction(0) + father.getShapeFunction(5) + father.getShapeFunction(4) + 0.25*father.getShapeFunction(1) +0.25*father.getShapeFunction(3);
// 	}
// 	
// 	if(dofIds.find(t->first) != dofIds.end() && dofIds.find(t->second) != dofIds.end() && dofIds.find(t->third) == dofIds.end())
// 	{
// 		return father.getShapeFunction(1)+father.getShapeFunction(0)+father.getShapeFunction(2) + 0.25*father.getShapeFunction(3) + 0.25*father.getShapeFunction(5);
// 	}
// }


	TriElement father(LINEAR) ;
// 	Function f ;
// 	for(size_t i = 0 ; i < t->getBoundingPoints().size() ; i++)
// 	{
// 		if(dofIds.find(&(t->getBoundingPoint(i))) != dofIds.end())
// 			f += father.getShapeFunction(i) ;
// 	}
// 	return f ;
	
	if(dofIds.find(t->first) != dofIds.end() && dofIds.find(t->second) == dofIds.end() && dofIds.find(t->third) == dofIds.end())
	{
		return father.getShapeFunction(0) ;
	}
	
	if(dofIds.find(t->first) == dofIds.end() && dofIds.find(t->second) != dofIds.end() && dofIds.find(t->third) == dofIds.end())
	{
		return father.getShapeFunction(1) ;
	}
	
	if(dofIds.find(t->first) == dofIds.end() && dofIds.find(t->second) == dofIds.end() && dofIds.find(t->third) != dofIds.end())
	{
		return father.getShapeFunction(2) ;
	}
	
	if(dofIds.find(t->first) == dofIds.end() && dofIds.find(t->second) != dofIds.end() && dofIds.find(t->third) != dofIds.end())
	{
		return 1-father.getShapeFunction(0) ;
	}
	
	if(dofIds.find(t->first) != dofIds.end() && dofIds.find(t->second) == dofIds.end() && dofIds.find(t->third) != dofIds.end())
	{
		return 1-father.getShapeFunction(1) ;
	}
	
	if(dofIds.find(t->first) != dofIds.end() && dofIds.find(t->second) != dofIds.end() && dofIds.find(t->third) == dofIds.end())
	{
		return 1-father.getShapeFunction(2) ;
	}
	
	return Function("1") ;
}

Function getBlendingFunction(const Circle * eps, const DelaunayTriangle * t)
{
	TriElement father(LINEAR) ;
	
	if(eps->in(*t->first) && !eps->in(*t->second) && !eps->in(*t->third) )
	{
		return father.getShapeFunction(0) ;
	}
	
	if(!eps->in(*t->first)  && eps->in(*t->second)  && !eps->in(*t->third) )
	{
		return father.getShapeFunction(1) ;
	}
	
	if(!eps->in(*t->first)  && !eps->in(*t->second)  && eps->in(*t->third) )
	{
		return father.getShapeFunction(2) ;
	}
	
	if(!eps->in(*t->first)  && eps->in(*t->second)  && eps->in(*t->third) )
	{
		return 1.-father.getShapeFunction(0) ;
	}
	
	if(eps->in(*t->first)  && !eps->in(*t->second)  && eps->in(*t->third) )
	{
		return 1.-father.getShapeFunction(1) ;
	}
	
	if(eps->in(*t->first)  && eps->in(*t->second)  && !eps->in(*t->third) )
	{
		return 1.-father.getShapeFunction(2) ;
	}
	
	return Function() ;
}

void BranchedCrack::enrichTip(size_t & lastId, Mesh<DelaunayTriangle,DelaunayTreeItem> * dt, const std::pair<Point *, double> & tip)
{
	Circle epsilon(enrichementRadius, Point(*(tip.first))) ;
	std::vector<DelaunayTriangle *> triangles = dt->getConflictingElements(&epsilon) ;
	std::vector<DelaunayTriangle *> test = dt->getConflictingElements(tip.first) ;
	if(test.empty())
		return ;
	std::map<Point *, size_t> done ;
	VirtualMachine vm ;
	double angle = tip.second;
	TriElement father (triangles[0]->getOrder()) ;
	std::valarray<Function> shapefunc = TriElement ( LINEAR ).getShapeFunctions() ;
	
	for(size_t  i = 0 ; i < triangles.size() ; i++)
	{
		
// 		if(!triangles[i]->in(*(tip.first)))
// 			continue ;
		enrichmentMap.insert(triangles[i]) ;
		
		if(!triangles[i]->enrichmentUpdated)
			triangles[i]->clearEnrichment(getPrimitive()) ;
		triangles[i]->enrichmentUpdated = true ;
		std::vector<Point> hint ;
		Line crossing ( *tip.first, Point(cos(angle), sin(angle)) ) ;
		if(crossing.intersects(static_cast<Triangle *>(triangles[i])))
		{
			std::vector<Point > completeIntersection = crossing.intersection ( static_cast<Triangle *>(triangles[i]) ) ;
			
			if(!completeIntersection.empty())
			{
				Point intersectionTransformed = triangles[i]->inLocalCoordinates ( completeIntersection[0] ) ;
				hint.push_back ( intersectionTransformed ) ;
			}
			
			if(triangles[i]->in(*tip.first))
			{
				Point singularityTransformed = triangles[i]->inLocalCoordinates ( *tip.first );
				hint.push_back ( singularityTransformed ) ;
			}
			
			if(completeIntersection.size() == 2)
			{
				Point intersectionBisTransformed = triangles[i]->inLocalCoordinates ( completeIntersection[1] ) ;
				hint.push_back ( intersectionBisTransformed ) ;
			}
		}
		
		
		Function x = triangles[i]->getXTransform() ;
		Function y = triangles[i]->getYTransform() ;

		double rotatedSingularityX = tip.first->getX()*cos ( angle ) + tip.first->getY()*sin ( angle ) ;
		double rotatedSingularityY = tip.first->getY()*cos ( angle ) - tip.first->getX()*sin ( angle ) ;
		Function rotatedX = x*cos ( angle ) + y*sin ( angle ) ;
		Function rotatedY = y*cos ( angle ) - x*sin ( angle );
		Function x_alt = x - tip.first->getX() ;
		Function y_alt = y - tip.first->getY() ;

		Function theta_alt = f_atan2 ( rotatedY-rotatedSingularityY, rotatedX-rotatedSingularityX );
		Function r_alt = Function(*tip.first, triangles[i]);//f_sqrt ( (x_alt^2)  + (y_alt^2) );
		

		Function x_ = x_alt ; //f_curvilinear_x(getPrimitive(), (tip.first == getHead()), x, y) ; 
		Function y_ = y_alt ; //f_curvilinear_y(getPrimitive(), (tip.first == getHead()), x, y) ; 


		Function theta = f_atan2 ( y_, x_ );
		
		Function f0 = f_sqrt ( r_alt ) *f_sin ( theta_alt*.5 ) ;
		Function f1 = f_sqrt ( r_alt ) *f_cos ( theta_alt*.5 ) ;
		Function f2 = f_sqrt ( r_alt ) *f_sin ( theta_alt*.5 ) *f_sin ( theta_alt ) ;
		Function f3 = f_sqrt ( r_alt ) *f_cos ( theta_alt*.5 ) *f_sin ( theta_alt ) ;
		f0.setNumberOfDerivatives(0);
		f1.setNumberOfDerivatives(0);
		f2.setNumberOfDerivatives(0);
		f3.setNumberOfDerivatives(0);
		bool hinted = false ;
		int pcount = 0 ;
		

			Function blend = Function("1") ;
		

		if(epsilon.in(*triangles[i]->first))
		{
			pcount++ ;
			int usedId = 0 ;
			if(done.find(triangles[i]->first) == done.end())
			{
				done[triangles[i]->first] = lastId ;
				usedId = lastId ;
				lastId += 4 ;
			}
			else
			{
				usedId = done[triangles[i]->first] ;
			}
			Point pt = triangles[i]->inLocalCoordinates( *triangles[i]->first) ;
			Function f = shapefunc[0]* ( f0 - vm.eval ( f0,  pt) ) * blend ;
			if(!hinted)
			{
				f.setIntegrationHint ( hint ) ;
				hinted = true ;
			}
			
			f.setPoint ( triangles[i]->first) ;
			f.setDofID ( usedId ) ;
			triangles[i]->setEnrichment (  f , getPrimitive()  ) ;
			
			f = shapefunc[0]* ( f1 - vm.eval ( f1, pt ) ) * blend ;
			f.setPoint ( triangles[i]->first) ;
			f.setDofID ( usedId+1 ) ;
			triangles[i]->setEnrichment ( f , getPrimitive() ) ;
			
			f = shapefunc[0]* ( f2 - vm.eval ( f2, pt ) ) * blend ;
			f.setPoint ( triangles[i]->first) ;
			f.setDofID ( usedId+2 ) ;
			triangles[i]->setEnrichment ( f , getPrimitive() ) ;
			
			f = shapefunc[0]* ( f3 - vm.eval ( f3, pt ) ) * blend ;
			f.setPoint (triangles[i]->first ) ;
			f.setDofID ( usedId+3 ) ;
			triangles[i]->setEnrichment ( f , getPrimitive() ) ;
		}
		
		if(epsilon.in(*triangles[i]->second))
		{
			pcount++ ;
			int usedId = 0 ;
			if(done.find(triangles[i]->second) == done.end())
			{
				done[triangles[i]->second] = lastId ;
				usedId = lastId ;
				lastId += 4 ;
			}
			else
			{
				usedId = done[triangles[i]->second] ;
			}
			Point pt = triangles[i]->inLocalCoordinates( *triangles[i]->second) ;
			Function f = shapefunc[1]* ( f0 - vm.eval ( f0,  pt) ) * blend ;
			if(!hinted)
			{
				f.setIntegrationHint ( hint ) ;
				hinted = true ;
			}
			
			f.setPoint ( triangles[i]->second) ;
			f.setDofID ( usedId ) ;
			triangles[i]->setEnrichment (  f , getPrimitive()  ) ;
			
			f = shapefunc[1]* ( f1 - vm.eval ( f1, pt ) ) * blend ;
			f.setPoint ( triangles[i]->second) ;
			f.setDofID ( usedId+1 ) ;
			triangles[i]->setEnrichment ( f , getPrimitive() ) ;
			
			f = shapefunc[1]* ( f2 - vm.eval ( f2, pt ) ) * blend ;
			f.setPoint ( triangles[i]->second) ;
			f.setDofID ( usedId+2 ) ;
			triangles[i]->setEnrichment ( f , getPrimitive() ) ;
			
			f = shapefunc[1]* ( f3 - vm.eval ( f3, pt ) ) * blend ;
			f.setPoint (triangles[i]->second ) ;
			f.setDofID ( usedId+3 ) ;
			triangles[i]->setEnrichment ( f , getPrimitive() ) ;
		}
		
		if(epsilon.in(*triangles[i]->third))
		{
			pcount++ ;
			int usedId = 0 ;
			if(done.find(triangles[i]->third) == done.end())
			{
				done[triangles[i]->third] = lastId ;
				usedId = lastId ;
				lastId += 4 ;
			}
			else
			{
				usedId = done[triangles[i]->third] ;
			}
			Point pt = triangles[i]->inLocalCoordinates( *triangles[i]->third) ;
			Function f = shapefunc[2]* ( f0 - vm.eval ( f0,  pt) ) * blend ;
			if(!hinted)
			{
				f.setIntegrationHint ( hint ) ;
				hinted = true ;
			}
			
			f.setPoint ( triangles[i]->third) ;
			f.setDofID ( usedId ) ;
			triangles[i]->setEnrichment (  f , getPrimitive()  ) ;
			
			f = shapefunc[2]* ( f1 - vm.eval ( f1, pt ) ) * blend ;
			f.setPoint ( triangles[i]->third) ;
			f.setDofID ( usedId+1 ) ;
			triangles[i]->setEnrichment ( f , getPrimitive() ) ;
			
			f = shapefunc[2]* ( f2 - vm.eval ( f2, pt ) ) * blend ;
			f.setPoint ( triangles[i]->third) ;
			f.setDofID ( usedId+2 ) ;
			triangles[i]->setEnrichment ( f , getPrimitive() ) ;
			
			f = shapefunc[2]* ( f3 - vm.eval ( f3, pt ) ) * blend ;
			f.setPoint (triangles[i]->third ) ;
			f.setDofID ( usedId+3 ) ;
			triangles[i]->setEnrichment ( f , getPrimitive() ) ;
		}
		
		if(pcount == (int)triangles[i]->getBoundingPoints().size())
			tipEnrichmentMap.insert(triangles[i]) ;
		
	}

	for(std::map<Point *, size_t>::const_iterator i = done.begin() ; i != done.end() ; ++i)
	{
		donePoints.insert(i->first) ;
	}
}

void BranchedCrack::enrichForks(size_t & startid, Mesh<DelaunayTriangle,DelaunayTreeItem> * dt)
{
	for(size_t i =  0 ; i < forks.size() ; i++)
	{
		enrichSegmentedLine(startid, dt, forks[i], &forks[i]->getBoundingPoint(1)) ;
	}
}

void BranchedCrack::enrichBranches(size_t & startid, Mesh<DelaunayTriangle,DelaunayTreeItem> * dt)
{
	
	for(size_t i =  0 ; i < branches.size() ; i++)
	{
		enrichSegmentedLine(startid, dt, branches[i]) ;
	}
}

void BranchedCrack::enrichSegmentedLine(size_t & lastId, Mesh<DelaunayTriangle,DelaunayTreeItem> * dt, const SegmentedLine * line)
{
	std::vector<Segment> intersectingSegments ;
	for ( size_t j = 1 ; j < line->getBoundingPoints().size() ; j++ )
	{
		intersectingSegments.push_back ( Segment ( line->getBoundingPoint ( j-1 ), line->getBoundingPoint ( j ) ) ) ;
	}
		
	std::vector<DelaunayTriangle *> tris = dt->getConflictingElements(line) ;
	
	std::valarray<Function> shapefunc = TriElement ( LINEAR ).getShapeFunctions() ;
	VirtualMachine vm ;
	std::map<Point *, size_t> done ;
	
	for ( size_t i = 0 ; i < tris.size() ; i++ )
	{
		if(tipEnrichmentMap.find(tris[i]) != tipEnrichmentMap.end())
		{
			continue ;
		}
		
		if(forkEnrichmentMap.find(tris[i]) != forkEnrichmentMap.end())
		{
			continue ;
		}
		
		if(!tris[i]->enrichmentUpdated)
			tris[i]->clearEnrichment(line) ;
		
		tris[i]->enrichmentUpdated = true ;
		DelaunayTriangle *e = tris[i] ;
		enrichmentMap.insert(e) ;
		std::vector<Point> intersection = line->intersection ( e->getPrimitive() ) ;
		
		for ( size_t j = 0 ; j < line->getBoundingPoints().size() ; j++ )
		{
			if(e->in(line->getBoundingPoint(j)))
			{
				intersection.push_back(line->getBoundingPoint(j));
			}
		}

		if(intersection.empty())
		{
			continue ;
		}
// 			e->setNonLinearBehaviour( new TwoDCohesiveForces(e, dynamic_cast<SegmentedLine *>(this)) ) ;


		std::vector<Point> hint ;
		std::vector<Point> transformed ;

		for ( size_t k = 0 ; k < intersection.size() ; k++ )
		{
			transformed.push_back ( e->inLocalCoordinates ( intersection[k] ) ) ;
			hint.push_back ( transformed.back() ) ;
		}
		

		Function s ( intersectingSegments, e ) ;
		Function dx = e->getXTransform() ;
		Function dy = e->getYTransform() ;
		s.setVariableTransform(XI, dx);
		s.setVariableTransform(ETA, dy);
		s.setNumberOfDerivatives(0);
// 		Function s = f_sign(e->getYTransform()) ;

		int usedId = 0 ;
		if(done.find(e->first) == done.end())
		{
			done[e->first] = lastId ;
			usedId = lastId ;
			lastId++ ;
		}
		else
		{
			usedId = done[e->first] ;
		}

		double val = vm.eval ( s, Point ( 0,1 ) ) ;
		Function f = shapefunc[0]*(s - val) ;
		f.setIntegrationHint ( hint ) ;
		f.setPoint ( e->first ) ;
		f.setDofID ( usedId ) ;
		e->setEnrichment ( f , getPrimitive() ) ;
		
		if(done.find(e->second) == done.end())
		{
				done[e->second] = lastId ;
				usedId = lastId ;
				lastId++ ;
		}
		else
		{
			usedId = done[e->second] ;
		}
		f = shapefunc[1]*(s - vm.eval(s, Point(0,0))) ;
		f.setPoint ( e->second) ;
		f.setDofID ( usedId ) ;
		e->setEnrichment ( f , getPrimitive() ) ;

		if(done.find(e->third) == done.end())
		{
			done[e->third] = lastId ;
			usedId = lastId ;
			lastId++ ;
		}
		else
		{
			usedId = done[e->third] ;
		}
		
		f = shapefunc[2] *(s - vm.eval(s, Point(1,0)))  ;

		f.setPoint (e->third ) ;
		f.setDofID ( usedId ) ;
		e->setEnrichment ( f , getPrimitive() ) ;

	}
	
	for(std::map<Point *, size_t>::const_iterator i = done.begin() ; i != done.end() ; ++i)
	{
		donePoints.insert(i->first) ;
	}
	
}

void BranchedCrack::enrichSegmentedLine(size_t & lastId, Mesh<DelaunayTriangle,DelaunayTreeItem> * dt, const SegmentedLine * line, const Point * target)
{
	std::vector<DelaunayTriangle *> tris = dt->getConflictingElements(line) ;
	std::vector<Segment> intersectingSegments ;
	for ( size_t j = 1 ; j < line->getBoundingPoints().size() ; j++ )
	{
		intersectingSegments.push_back ( Segment ( line->getBoundingPoint ( j-1 ), line->getBoundingPoint ( j ) ) ) ;
	}
	
	std::valarray<Function> shapefunc = TriElement ( LINEAR ).getShapeFunctions() ;
	VirtualMachine vm ;
	std::map<Point *, size_t> done ;
	DelaunayTriangle * centralTriangle = nullptr;
	
	for ( size_t i = 0 ; i < tris.size() ; i++ )
	{
		if(tris[i]->in(*target))
		{
			centralTriangle = tris[i] ;
			break ;
		}
	}

	for ( size_t i = 0 ; i < tris.size() ; i++ )
	{
		bool performEnrichment = tris[i]->in(*target) ;
		if(centralTriangle)
		{
            std::vector<DelaunayTriangle * > neighbourhood = dt->getNeighbourhood(centralTriangle) ;
			for ( size_t j = 0 ; j < neighbourhood.size() ; j++ )
			{
				if ( neighbourhood[j]->isAlive() && neighbourhood[j] == tris[i] )
				{
					performEnrichment = true ;
					break ;
				}
			}
		}
		else
			performEnrichment = false ;
		
		if(performEnrichment)
		{
			
// 			if(tipEnrichmentMap.find(tris[i]) != tipEnrichmentMap.end())
// 			{
// 				continue ;
// 			}
			
			
			if(!tris[i]->enrichmentUpdated)
				tris[i]->clearEnrichment(line) ;
			
			tris[i]->enrichmentUpdated = true ;
			DelaunayTriangle *e = tris[i] ;
			enrichmentMap.insert(e) ;
			forkEnrichmentMap.insert(e) ;
			std::vector<Point> intersection = line->intersection ( e->getPrimitive() ) ;
			
			for ( size_t j = 0 ; j < line->getBoundingPoints().size() ; j++ )
			{
				if(e->in(line->getBoundingPoint(j)))
				{
					intersection.push_back(line->getBoundingPoint(j));
				}
			}

			if(intersection.empty())
				continue ;

	// 			e->setNonLinearBehaviour( new TwoDCohesiveForces(e, dynamic_cast<SegmentedLine *>(this)) ) ;


			std::vector<Point> hint ;
			std::vector<Point> transformed ;

			for ( size_t k = 0 ; k < intersection.size() ; k++ )
			{
				transformed.push_back ( e->inLocalCoordinates ( intersection[k] ) ) ;
			}


			for ( size_t k = 0 ; k < transformed.size() ; k++ )
			{
				hint.push_back ( transformed[k] ) ;
			}
			
			Function s ( intersectingSegments,e) ;

			int usedId = 0 ;
			if(done.find(e->first) == done.end())
			{
				done[e->first] = lastId ;
				usedId = lastId ;
				lastId++ ;
			}
			else
			{
				usedId = done[e->first] ;
			}

			Function f = shapefunc[0]* ( s - vm.eval ( s, Point ( 0,1 ) ) ) ;
			f.setIntegrationHint ( hint ) ;
			f.setPoint ( e->first ) ;
			f.setDofID ( usedId ) ;
			e->setEnrichment ( f , getPrimitive() ) ;
			
			if(done.find(e->second) == done.end())
			{
				done[e->second] = lastId ;
				usedId = lastId ;
				lastId++ ;
			}
			else
			{
				usedId = done[e->second] ;
			}
			f = shapefunc[1]* ( s - vm.eval ( s, Point ( 0,0 ) ) ) ;
			f.setPoint ( e->second ) ;
			f.setDofID ( usedId ) ;
			e->setEnrichment ( f , getPrimitive() ) ;

			if(done.find(e->third) == done.end())
			{
				done[e->third] = lastId ;
				usedId = lastId ;
				lastId++ ;
			}
			else
			{
				usedId = done[e->third] ;
			}
			f = shapefunc[2] * ( s - vm.eval ( s, Point ( 1,0 ) ) ) ;

			f.setPoint ( e->third ) ;
			f.setDofID ( usedId ) ;
			e->setEnrichment ( f , getPrimitive() ) ;

			std::vector<DelaunayTriangle *> toEnrichAlso ;
            std::vector<DelaunayTriangle *> neighbourhood = dt->getNeighbourhood(e) ;
			for (auto & n : neighbourhood )
			{
				if ( n->isAlive() 
					&& !line->intersects(n->getPrimitive()))
					toEnrichAlso.push_back ( n ) ;
			}


			for ( size_t j = 0 ; j < toEnrichAlso.size() ; j++ )
			{
				if(!toEnrichAlso[j]->enrichmentUpdated)
					toEnrichAlso[j]->clearEnrichment( line) ;
				toEnrichAlso[j]->enrichmentUpdated = true ;
				
				DelaunayTriangle * elem = toEnrichAlso[j] ;
				enrichmentMap.insert(elem) ;

				transformed.clear() ;

				for ( size_t k = 0 ; k < intersection.size() ; k++ )
				{
					transformed.push_back ( elem->inLocalCoordinates ( intersection[k] ) ) ;
				}

				Function s_ ( intersectingSegments, elem) ;
				hint.clear() ;

				if ( done.find(elem->first) != done.end())
				{
					Function f = shapefunc[0]* ( s_ - vm.eval ( s_, Point ( 0,1 ) ) ) ;
					f.setPoint ( elem->first ) ;
					f.setDofID ( done[elem->first] ) ;
					elem->setEnrichment ( f , getPrimitive()) ;
					
				}

				if (done.find(elem->second) != done.end())
				{
					Function f = shapefunc[1]* ( s_ - vm.eval ( s_, Point ( 0,0 ) ) ) ;
					f.setPoint ( elem->second ) ;
					f.setDofID ( done[elem->second] ) ;
					elem->setEnrichment ( f , getPrimitive()) ;
				}
				if (done.find(elem->third) != done.end())
				{
					Function f = shapefunc[2]* ( s_ - vm.eval ( s_, Point ( 1,0 ) ) ) ;
					f.setPoint ( elem->third ) ;
					f.setDofID ( done[elem->third] ) ;
					elem->setEnrichment ( f , getPrimitive()) ;
				}
			}
		}
	}
	
	for(std::map<Point *, size_t>::const_iterator i = done.begin() ; i != done.end() ; ++i)
	{
		donePoints.insert(i->first) ;
	}
	
}

double BranchedCrack::getEnrichementRadius() const
{
	return enrichementRadius ;
}

void BranchedCrack::setEnrichementRadius(double newRadius)
{
	enrichementRadius = newRadius ;
}

const std::vector<std::pair<Point *, double> > & BranchedCrack::getTips() const
{
	return tips ;
}

const std::vector<SegmentedLine * > & BranchedCrack::getForks() const
{
	return forks ;
}

std::vector<SegmentedLine *> & BranchedCrack::getBranches()
{
	return branches ;
}

std::vector<std::pair<Point *, double> > & BranchedCrack::getTips()
{
	return tips ;
}

std::vector<SegmentedLine * > & BranchedCrack::getForks()
{
	return forks ;
}

bool BranchedCrack::isEmpty() const
{
	return branches.empty() ;
}

void BranchedCrack::enrich(size_t & counter,  Mesh<DelaunayTriangle,DelaunayTreeItem> * dtree)
{
	freeIds.clear();
	std::vector<DelaunayTriangle *> toClear ;
	for(auto i = enrichmentMap.begin() ; i != enrichmentMap.end() ; ++i)
	{
		toClear.push_back(*i) ;
	}
	
	for(auto i = toClear.begin() ; i != toClear.end() ; ++i)
	{
		std::vector<size_t> idpts = (*i)->clearEnrichment(getPrimitive()) ;
		(*i)->enrichmentUpdated = true;
	}
	
	enrichmentMap.clear() ;
	tipEnrichmentMap.clear() ;
	forkEnrichmentMap.clear() ;
	enrichTips(counter, dtree) ;
	enrichForks(counter, dtree) ;
	enrichBranches(counter, dtree) ;
}

std::vector<DelaunayTriangle*> BranchedCrack::getElements2D(FeatureTree* dt)
{
	
	std::vector<DelaunayTriangle*> ret ;

	for(size_t i = 0 ; i < branches.size() ; i++)
	{
		std::vector<DelaunayTriangle*> tris = dt->get2DMesh()->getConflictingElements(branches[i]) ;
		ret.insert(ret.end(), tris.begin(), tris.end()) ;
	}
	
	std::sort(ret.begin(), ret.end()) ;
	auto e = std::unique(ret.begin(), ret.end()) ;
	ret.erase(e, ret.end()) ;
	
	return ret ;
}

std::vector<DelaunayTetrahedron*> BranchedCrack::getElements3D(FeatureTree* dt)
{
	return std::vector<DelaunayTetrahedron*>() ;
}

bool BranchedCrack::interacts(Feature* f, double d) const
{
	bool ret = false ;

	for(size_t i = 0 ; i < branches.size() ; i++)
	{
		ret = ret || f->intersects(branches[i]) ;

		if (ret)
			return true ;
	}

	return false ;
}

std::vector<Amie::Geometry*> BranchedCrack::getRefinementZones(size_t level ) const
{
	std::vector<Geometry *> ret ;
	
	if ( level > 0 )
	{
		for ( size_t i = 0 ; i < tips.size() ; i++ )
		{
			ret.push_back ( new Circle ( 0.2, *tips[i].first ) ) ;
		}
	}
	if ( level > 1 )
	{
		for ( size_t i = 0 ; i < tips.size() ; i++ )
		{
			ret.push_back ( new Circle ( 0.15, *tips[i].first ) ) ;
		}
	}
	if ( level > 3 )
	{
		for ( size_t i = 0 ; i < tips.size() ; i++ )
		{
			ret.push_back ( new Circle ( 0.1, *tips[i].first ) ) ;
		}
	}
	return ret ;
}

void BranchedCrack::print() const
{
	std::cout << " == Branching crack == " << std::endl ;
	std::cout << " branches : " << branches.size() << std::endl;
	for(size_t i = 0 ; i < branches.size() ; i++)
	{
		std::cout << "branch " << i << std::endl ;
		
		for(size_t j = 0 ; j < branches[i]->getBoundingPoints().size() ; j++)
		{
			std::cout << "("<< branches[i]->getBoundingPoint(j).getX() << ", " << branches[i]->getBoundingPoint(j).getY() << ")" << std::endl ;
		}
	}

	
	std::cout << " tips : " << tips.size() << std::endl;
	for(size_t i = 0 ; i < tips.size() ; i++)
	{
		std::cout << "tip " << i << " : "<< "("<< tips[i].first->getX() << ", " << tips[i].first->getY() << ")"<< std::endl ;
	}
	std::cout << " forks : " << forks.size() << std::endl;
	for(size_t i = 0 ; i < forks.size() ; i++)
	{
		std::cout << "fork " << i << std::endl ;
		
		for(size_t j = 0 ; j < forks[i]->getBoundingPoints().size() ; j++)
		{
			std::cout << "("<< forks[i]->getBoundingPoint(j).getX() << ", " << forks[i]->getBoundingPoint(j).getY() << ")" << std::endl ;
		}
		
	}
	std::cout << std::endl ;

}

void BranchedCrack::sample(double linearDensity, double surfaceDensityFactor, Sampler * sampler)
{
	return ;
}

BranchedCrack::~BranchedCrack() { }

bool BranchedCrack::isVoid(const Point&) const
{
	return false ;
}

std::vector<Point*> BranchedCrack::getSamplingPoints() const
{
	return std::vector<Point*>() ;
}

bool BranchedCrack::enrichmentTarget(DelaunayTriangle* tri)
{
	return enrichmentMap.find(tri) != enrichmentMap.end() ;
}

void BranchedCrack::setScorePropagationMethod()
{
	scorePropagation = true ;
}
	
void BranchedCrack::setEnergyPropagationMethod()
{
	scorePropagation = false ;
}

void BranchedCrack::step(double dt, Vector* v, Mesh< DelaunayTriangle, DelaunayTreeItem >* dtree)
{
	changed = false ;
	isUpdated = false ;
// 	return ;
	std::vector<Point *> tipsToGrow ; 
	std::vector<double> angles ; 
	double pdistance = enrichementRadius * .5 ;
	std::vector<double> originalAngles ;
	if(!scorePropagation)
	{
		for(size_t i = 0 ; i < tips.size() ; i++)
		{
			std::pair<double, double> energy  = computeJIntegralAtTip(tips[i], dtree);
			if(energy.first*energy.first + energy.second*energy.second > 0)
			{
				tipsToGrow.push_back(tips[i].first) ;
				originalAngles.push_back(tips[i].second);
				changed = true ;
				isUpdated = true ;
				angles.push_back(propagationAngleFromTip(tips[i], dtree)) ;
			}
		}
	}
	else
	{
		pdistance = -1 ;
		for(size_t i = 0 ; i < tips.size() ; i++)
		{
			
			DelaunayTriangle * tri = dtree->getUniqueConflictingElement(tips[i].first) ;
			pdistance = std::max(tri->getRadius()*2.,pdistance) ;
			double score = -1 ;
			if(tri->getBehaviour()->getFractureCriterion())
				score = tri->getBehaviour()->getFractureCriterion()->getScoreAtState() ;
			if(score > 0 || true)
			{
				tipsToGrow.push_back(tips[i].first) ;
				changed = true ;
				isUpdated = true ;
				angles.push_back(tri->getBehaviour()->getFractureCriterion()->getSmoothedField(PRINCIPAL_STRESS_ANGLE_FIELD, tri->getState())[0]+M_PI*.5) ;
			}
		}
	}
	
	for(size_t i = 0 ; i < tipsToGrow.size() ; i++)
	{
		double sign = -1 ;
		if(originalAngles[i] < 0)
			sign = 1 ;
				
		grow(tipsToGrow[i], new Point(tipsToGrow[i]->getX() + sign * pdistance * cos(angles[i]),
																	tipsToGrow[i]->getY() + sign * pdistance * sin(angles[i]))) ;
	}
}

bool BranchedCrack::moved() const
{
	return changed ;
}
}
