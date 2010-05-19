// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2006-2007
//
// Copyright: See COPYING file that comes with this distribution
//

#include "crack.h"

#include <fstream>

using namespace Mu ;

BranchedCrack::BranchedCrack(Feature *father, Point * a, Point * b) : EnrichmentFeature(father), SegmentedLine(std::valarray<Point * >(2))
{
	boundingPoints[0] = a ;
	boundingPoints[1] = b ;
	this->SegmentedLine::center = *a ;

	if(father->in(*a))
		tips.push_back(std::make_pair(a, atan2(a->y-b->y, a->x-b->x))) ;
	if(father->in(*b))
		tips.push_back(std::make_pair(b, atan2(b->y-a->y, b->x-a->x))) ;

	changed = false ;

	enrichementRadius = 0.00075 ;
	
}

bool operator ==(const std::pair<Mu::Point*, double> & a, const Mu::Point* b)
{
	return a.first == b ;
}

BranchedCrack::BranchedCrack(Point * a, Point * b) : EnrichmentFeature(NULL), SegmentedLine(std::valarray<Point * >(2))
{
	if(a->x < b->x)
	{
		boundingPoints[0] = a ;
		boundingPoints[1] = b ;
		this->SegmentedLine::center = *a ;
	}
	else
	{
		boundingPoints[0] = a ;
		boundingPoints[1] = b ;
		this->SegmentedLine::center = *b ;
	}
	branches.push_back ( this ) ;

	tips.push_back(std::make_pair(boundingPoints[0], atan2(boundingPoints[0]->y-boundingPoints[1]->y, boundingPoints[0]->x-boundingPoints[1]->x))) ;
// 	if(std::abs(tips.back().second - M_PI) < std::numeric_limits<double>::epsilon())
// 		tips.back().second = - M_PI ; 
	tips.push_back(std::make_pair(boundingPoints[1], atan2(boundingPoints[1]->y-boundingPoints[0]->y, boundingPoints[1]->x-boundingPoints[0]->x))) ;
// 	if(std::abs(tips.back().second - M_PI) < std::numeric_limits<double>::epsilon())
// 		tips.back().second = - M_PI ; 

	changed = true ;

	enrichementRadius = 0.00075 ;
}

void BranchedCrack::branch(Point* fromTip, const std::vector<Point *> & newTip)
{
	for(std::vector<std::pair<Point *, double> >::iterator i = tips.begin() ; i != tips.end() ; ++i)
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
		double angle = atan2((*j)->y-fromTip->y, (*j)->x-fromTip->x) ;
		tips.push_back ( std::make_pair(*j, angle) ) ;
		sortedNewTips[angle] = *j ;
	}
	
	std::vector<Point *> newTipSorted ;
	for(std::map<double, Point *>::const_iterator j = sortedNewTips.begin() ; j !=sortedNewTips.end() ; ++j)
		newTipSorted.push_back(j->second) ;
	
	SegmentedLine * branchToExtend = NULL ;
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
	for(std::vector<std::pair<Point *, double> >::iterator i = tips.begin() ; i !=tips.end() ; ++i)
	{
		if(i->first == fromTip)
		{
			tips.erase(i) ;
			break ;
		}
	}
	tips.push_back ( std::make_pair(newTip0, atan2(newTip0->y-fromTip->y, newTip0->x-fromTip->x)) ) ;
	tips.push_back ( std::make_pair(newTip1, atan2(newTip1->y-fromTip->y, newTip1->x-fromTip->x)) ) ;
	std::valarray<Point * > newBranch ( 2 ) ;
	newBranch[0] = fromTip ;
	newBranch[1] = newTip0 ;

	SegmentedLine * branchToExtend = NULL ;
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
		
	Point lastDir( cos(tip.second) , sin(tip.second)) ;
	
	Circle c(enrichementRadius, *(tip.first)) ;
	std::vector<DelaunayTriangle*> disk = dtree->getConflictingElements (&c) ;
	std::vector<DelaunayTriangle*> tris ;
	while( disk.size() < 16)
	{
		c.setRadius(c.getRadius()*2.) ;
		disk = dtree->getConflictingElements ( &c ) ;
	}


	DelaunayTriangle * tailElem = NULL;
	for ( size_t i = 0 ; i < disk.size() ; i++ )
	{
		
		if ( !disk[i]->in ( *getTail() )
		     && (!Line(*tip.first, lastDir).intersects(disk[i]->getPrimitive())
				|| c.in(disk[i]->getCenter()))
			)
		{
			tris.push_back ( disk[i] ) ;
		}
		else if ( disk[i]->in ( *tip.first ) )
		{
			tailElem = disk[i] ;
		}
	}
	
	if ( !tris.empty() )
	{
		
		for ( size_t i = 0 ; i < tris.size() ; i++ )
		{
			
			DelaunayTriangle * current = tris[i] ;
			
			for ( size_t j = 0 ; j < current->getBoundingPoints().size() ; j++ )
			{
				
				Point currentDir ( current->getBoundingPoint ( j ).x-tip.first->x,
				                   current->getBoundingPoint ( j ).y-tip.first->y ) ;
				
				if ( ( currentDir*lastDir ) > 0 )
				{
					Vector principalStresses = current->getState().getPrincipalStresses (
						current->getBoundingPoint ( j ) ) ;
					double maxPrincipalStressCurrent = std::abs ( principalStresses[0] );
					
					if ( current->getBehaviour()->type == VOID_BEHAVIOUR )
						maxPrincipalStressCurrent = 1000 ;
					
					if ( currentDir.norm() > 1e-8 )
					{
						aangle += current->getState().getPrincipalAngle (
							current->getBoundingPoint ( j ) ) ;
						acount++ ;
					}
				}
			}
		}
	}
	
	if ( acount )
		aangle /= acount ;
	aangle += M_PI*.25 ;
	if(std::abs(aangle-tip.second) > M_PI*.5)
		aangle += M_PI ;
	return aangle ;
}

std::pair<double, double> BranchedCrack::computeJIntegralAtTip ( std::pair<Point *, double> & tip, Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree )
{
	Point direction ( cos(tip.second), sin(tip.second)) ;
	Segment tipSegment ( *(tip.first) , direction) ;
	Circle c ( enrichementRadius, *(tip.first) ) ;
	
	
	std::vector<DelaunayTriangle *> disk = dtree->getConflictingElements ( &c ) ;
	std::vector<DelaunayTriangle *> ring ;
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
			std::vector<std::pair<Point, double> > gaussPoints = gamma[j].first->getGaussPoints() ;
			gaussPoints.push_back ( std::make_pair ( gamma[j].first->midPoint(), ( double ) 1 ) );
			Point n = gamma[j].first->normal ( gamma[j].second->getCenter() ) ;
			
			
			Point ln = gamma[j].second->inLocalCoordinates ( n ) ;
			Vector stepLengthal ( 2 ) ; stepLengthal[0] = n.x ; stepLengthal[1] = n.y ;
			Vector d ( 2 ) ; d[0] = -n.y ;   d[1] = n.x ;
			double ilocal0 = 0 ;
			double ilocal1 = 0 ;
			
			Point localVector = gamma[j].second->inLocalCoordinates ( tipSegment.vector() ) ;
			Point localNormal ( -localVector.y, localVector.x ) ;
			
			for ( size_t k = 0 ; k < gaussPoints.size() ; k++ )
			{
				Point localGaussPoint = gamma[j].second->inLocalCoordinates ( gaussPoints[k].first ) ;
				Matrix sigma (gamma[j].second->getState().getStressMatrix ( localGaussPoint,true )) ;
				Matrix epsilon (gamma[j].second->getState().getStrainMatrix ( localGaussPoint,true )) ;
				double sigma_epsilon = sigma[0][0]*epsilon[0][0] + sigma[0][1]*epsilon[0][1] + sigma[1][0]*epsilon[1][0] + sigma[1][1]*epsilon[1][1];//0.5*std::inner_product(&sigma.array()[0], &sigma.array()[sigma.size()], &epsilon.array()[0], 0. );
				
				Vector T = sigma * stepLengthal ;
				
				
				ilocal0 += ( 0.5*sigma_epsilon*stepLengthal[0] - ( T[0]*epsilon[0][0] + T[1]*epsilon[0][1] ) ) *gaussPoints[k].second  ;
				ilocal1 += ( 0.5*sigma_epsilon*stepLengthal[1] - ( T[0]*epsilon[0][1] + T[1]*epsilon[1][1] ) ) *gaussPoints[k].second;
				
			}
			freeEnergy0 += ilocal0*gamma[j].first->norm() ;
			freeEnergy1 += ilocal1*gamma[j].first->norm() ;
		}
	}
	return std::make_pair ( freeEnergy0, freeEnergy1 ) ;
}

void BranchedCrack::grow( Point* fromTip, Point* newTip)
{
	for(std::vector<std::pair<Point *, double> >::iterator i = tips.begin() ; i !=tips.end() ; ++i)
	{
		if(i->first == fromTip)
		{
			tips.erase(i) ;
			break ;
		}
	}

	SegmentedLine * branchToExtend = NULL ;
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

	std::vector<Point> inter = Circle(enrichementRadius, *newTip).intersection(branchToExtend) ;
	tips.push_back ( std::make_pair(newTip, atan2(newTip->y-inter[0].y, newTip->x-inter[0].x)) ) ;

	
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
	SegmentedLine * toBranch = NULL;
	
	
	for(size_t i = 0 ; i < newSet.getTips().size() ; i++)
	{
		std::pair<Point *, double> tipForMergePotential = newSet.getTips()[i] ;
		double distPotential = maxdist ;
		std::pair<Point*, Point*> targetSegmentPotential ;
		SegmentedLine * toBranchPotential = NULL;
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
	SegmentedLine * fromBranch = NULL ;
	bool isHead = false ;
	bool isTail = false ;
	
	for(size_t i = 0 ; i < newSet.getBranches().size() ; i++)
	{
		if(newSet.getBranches()[i]->getHead() == tipForMerge.first)
		{
			fromBranch = newSet.getBranches()[i] ;
			isHead = true ;
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
		toBranch->getBoundingPoint(k).print() ;
		toBranch->getBoundingPoint(k+1).print() ;
		targetSegment.first->print() ;
		targetSegment.second->print() ;
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
	for(size_t i =  0 ; i < tips.size() ; i++)
	{
		enrichTip(startid, dt, tips[i]) ;
	}
}

void BranchedCrack::enrichTip(size_t & startid, Mesh<DelaunayTriangle,DelaunayTreeItem> * dt, const std::pair<Point *, double> & tip)
{
	Circle epsilon(enrichementRadius, Point(*(tip.first))) ;
	std::vector<DelaunayTriangle *> triangles = dt->getConflictingElements(&epsilon) ;
	std::map<Point *, size_t> done ;
	VirtualMachine vm ;
	double angle = tip.second ;
	std::valarray<Function> shapefunc(3) ;
	Matrix xi(2,2) ; xi[1][0] = 1 ;
	Matrix eta(2,2) ; eta[0][1] = 1 ;
	Matrix one(2,2) ; one[0][0] = 1 ;

	shapefunc[0] = Function(eta) ;
	shapefunc[1] = Function(one-xi-eta) ;
	shapefunc[2] = Function(xi) ;
	
	std::valarray <Point> pointLocal(3) ;
	pointLocal[0] = Point ( 0,1 ) ;
	pointLocal[1] = Point ( 0,0 ) ;
	pointLocal[2] = Point ( 1,0 )  ;
	
	for(size_t  i = 0 ; i < triangles.size() ; i++)
	{
		enrichmentMap.insert(triangles[i]) ;
		
		if(!triangles[i]->enrichmentUpdated)
			triangles[i]->clearEnrichment(static_cast<SegmentedLine *>(this)) ;
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
			
			if(!completeIntersection.size() == 2)
			{
				Point intersectionBisTransformed = triangles[i]->inLocalCoordinates ( completeIntersection[1] ) ;
				hint.push_back ( intersectionBisTransformed ) ;
			}
		}
		
		
		Function x = triangles[i]->getXTransform() ;
		Function y = triangles[i]->getYTransform() ;

		double rotatedSingularityX = tip.first->x*cos ( angle ) + tip.first->y*sin ( angle ) ;
		double rotatedSingularityY = tip.first->y*cos ( angle ) - tip.first->x*sin ( angle ) ;
		Function rotatedX = x*cos ( angle ) + y*sin ( angle ) ;
		Function rotatedY = y*cos ( angle ) - x*sin ( angle );
		Function x_alt = x - tip.first->x ;
		Function y_alt = y - tip.first->y ;

		Function theta_alt = f_atan2 ( rotatedY-rotatedSingularityY, rotatedX-rotatedSingularityX );
		Function r_alt = f_sqrt ( (x_alt^2)  + (y_alt^2) );


		Function x_ = f_curvilinear_x(getPrimitive(), (tip.first == getHead()), x, y) ; 
		Function y_ = f_curvilinear_y(getPrimitive(), (tip.first == getHead()), x, y) ; 


		Function theta = f_atan2 ( y_, x_ );
		Function r = f_sqrt ( (x_^2)  + (y_^2) );
		
		Function f0 = f_sqrt ( r_alt ) *f_sin ( theta_alt/2 );
		Function f1 = f_sqrt ( r_alt ) *f_cos ( theta_alt/2 );
		Function f2 = f_sqrt ( r_alt ) *f_sin ( theta_alt/2 ) *f_sin ( theta_alt );
		Function f3 = f_sqrt ( r_alt ) *f_cos ( theta_alt/2 ) *f_sin ( theta_alt );
// 		Function f4 = f_sqrt ( r ) *f_sin ( theta/2 ) *f_sin ( theta * 2);
// 		Function f5 = f_sqrt ( r ) *f_cos ( theta/2 ) *f_sin ( theta * 2);
		

// 		if(!triangles[i]->in(*tip.first) && getBoundingPoints().size() == 3)
// 		{
// 			for(double k = 0  ; k < 100 ; k++)
// 			{
// 				for(double l = 0  ; l < 100 ; l++)
// 				{
// 					if( k+l < 100 )
// 						std::cout << vm.eval(theta, k/100., l/100.) << "  "<< std::flush ;
// 					else
// 						std::cout << "0  " << std::flush ;
// 				}
// 				std::cout << std::endl ;
// 			}
// // 			exit(0) ;
// 		}
// 
// 		if(!triangles[i]->in(*tip.first) && getBoundingPoints().size() == 3)
// 		{
// 			for(double k = 0  ; k < 100 ; k++)
// 			{
// 				for(double l = 0  ; l < 100 ; l++)
// 				{
// 					if( k+l < 100 )
// 						std::cout << vm.eval(r, k/100., l/100.) << "  " << std::flush ;
// 					else
// 						std::cout << "0  " << std::flush ;
// 				}
// 				std::cout << std::endl ;
// 			}
// 			
// 			exit(0) ;
// 		}
		
		std::vector<Point *> currentPoint ;
		currentPoint.push_back(triangles[i]->first) ;
		currentPoint.push_back(triangles[i]->second) ;
		currentPoint.push_back(triangles[i]->third) ;
		
		bool hinted = false ;
		int pcount = 0 ;
		for(size_t p = 0 ; p < 3 ;p++)
		{
			if(epsilon.in(*currentPoint[p]))
			{
				pcount++ ;
				int usedId = 0 ;
				if(done.find(currentPoint[p]) == done.end())
				{
					if(freeIds.size() < 4)
					{
						done[currentPoint[p]] = startid ;
						usedId = startid ;
						startid+=4 ;
					}
					else
					{
						done[currentPoint[p]] = *freeIds.begin() ;
						usedId = *freeIds.begin() ;
						std::set<size_t>::iterator iter = freeIds.begin() ;
						++iter ; ++iter ; ++iter ; ++iter ;
						freeIds.erase(freeIds.begin(), iter);
					}
				}
				else
				{
					usedId = done[currentPoint[p]] ;
				}

				Function f = shapefunc[p]* ( f0 - vm.eval ( f0, pointLocal[p] ) ) ;
				if(!hinted)
				{
					f.setIntegrationHint ( hint ) ;
					hinted = true ;
				}
				f.setPoint ( currentPoint[p]) ;
				f.setDofID ( usedId ) ;
				triangles[i]->setEnrichment (  f , static_cast<SegmentedLine *>(this)  ) ;
				
				f = shapefunc[p]* ( f1 - vm.eval ( f1, pointLocal[p] ) ) ;
				f.setPoint ( currentPoint[p]) ;
				f.setDofID ( usedId+1 ) ;
				triangles[i]->setEnrichment ( f , static_cast<SegmentedLine *>(this) ) ;
				
				f = shapefunc[p]* ( f2 - vm.eval ( f2, pointLocal[p] ) ) ;
				f.setPoint ( currentPoint[p]) ;
				f.setDofID ( usedId+2 ) ;
				triangles[i]->setEnrichment ( f , static_cast<SegmentedLine *>(this) ) ;
				
				f = shapefunc[p]* ( f3 - vm.eval ( f3, pointLocal[p] ) ) ;
				f.setPoint (currentPoint[p] ) ;
				f.setDofID ( usedId+3 ) ;
				triangles[i]->setEnrichment ( f , static_cast<SegmentedLine *>(this) ) ;
				
// 				f = shapefunc[p]* ( f4 - vm.eval ( f4, pointLocal[p] ) ) ;
// 				f.setPoint ( currentPoint[p]) ;
// 				f.setDofID ( usedId+4 ) ;
// 				triangles[i]->setEnrichment ( f , static_cast<SegmentedLine *>(this) ) ;
// 				
// 				f = shapefunc[p]* ( f5 - vm.eval ( f5, pointLocal[p] ) ) ;
// 				f.setPoint (currentPoint[p] ) ;
// 				f.setDofID ( usedId+5 ) ;
// 				triangles[i]->setEnrichment ( f , static_cast<SegmentedLine *>(this) ) ;
			}
		}
		if(pcount == 3)
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

void BranchedCrack::enrichSegmentedLine(size_t & startid, Mesh<DelaunayTriangle,DelaunayTreeItem> * dt, const SegmentedLine * line)
{
	std::vector<DelaunayTriangle *> tris = dt->getConflictingElements(line) ;
	
	std::valarray<Function> shapefunc = TriElement ( LINEAR ).getShapeFunctions() ;
	VirtualMachine vm ;
	std::map<Point *, size_t> done ;
	
	for ( size_t i = 0 ; i < tris.size() ; i++ )
	{
		while(tipEnrichmentMap.find(tris[i]) != tipEnrichmentMap.end() && i < tris.size())
		{
			i++ ;
		}
		if(i >= tris.size())
			break ;
		if(tipEnrichmentMap.find(tris[i]) != tipEnrichmentMap.end())
			break ;
		
		
		if(!tris[i]->enrichmentUpdated)
			tris[i]->clearEnrichment(static_cast<SegmentedLine *>(this)) ;
		
		tris[i]->enrichmentUpdated = true ;
		DelaunayTriangle *e = tris[i] ;
		enrichmentMap.insert(e) ;
		std::vector<Point> intersection ;
		

		for ( size_t j = 1 ; j < line->getBoundingPoints().size() ; j++ )
		{
			Segment test ( getBoundingPoint ( j-1 ), getBoundingPoint ( j ) ) ;

			if ( test.intersects ( static_cast<Triangle *> ( e ) ) )
			{

				//there is either one or two intersection points
				std::vector<Point> temp_intersection = test.intersection ( static_cast<Triangle *> ( e ) ) ;


				if ( temp_intersection.size() == 2 ) //we simply cross the element
				{
					if ( Segment ( temp_intersection[0], temp_intersection[1] ).vector() * test.vector() < 0 )
					{
						std::swap ( temp_intersection[0], temp_intersection[1] ) ;
					}
					intersection = temp_intersection ;

				}
				else //then there are kinks
				{
					intersection.push_back ( temp_intersection[0] ) ;
					intersection.push_back ( getBoundingPoint ( j ) ) ;
					for ( size_t k = j ; k < this->getBoundingPoints().size()-1 ; k++ )
					{
						Segment test_for_kink ( getBoundingPoint ( k ), getBoundingPoint ( k+1 ) ) ;

						if ( !test_for_kink.intersects ( static_cast<Triangle *> ( e ) ) )
						{
							intersection.push_back ( getBoundingPoint ( k+1 ) ) ;
						}
						else
						{
							intersection.push_back ( test_for_kink.intersection ( static_cast<Triangle *> ( e ) ) [0] ) ;
						}
					}
				}
			}
		}


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
		
		std::vector<Segment> intersectingSegments ;
		for ( size_t j = 1 ; j < this->getBoundingPoints().size() ; j++ )
		{
			intersectingSegments.push_back ( Segment ( getBoundingPoint ( j-1 ), getBoundingPoint ( j ) ) ) ;
		}
		Function s ( intersectingSegments, e->getXTransform(), e->getYTransform() ) ;

		int usedId = 0 ;
		if(done.find(e->first) == done.end())
		{
			if(freeIds.empty())
			{
				done[e->first] = startid ;
				usedId = startid ;
				startid++ ;
			}
			else
			{
				done[e->first] = *freeIds.begin() ;
				usedId = *freeIds.begin() ;
				freeIds.erase(freeIds.begin());
			}
		}
		else
		{
			usedId = done[e->first] ;
		}

		Function f = shapefunc[0]* ( s - vm.eval ( s, Point ( 0,1 ) ) ) ;
		f.setIntegrationHint ( hint ) ;
		f.setPoint ( e->first ) ;
		f.setDofID ( usedId ) ;
		e->setEnrichment ( f , static_cast<SegmentedLine *>(this) ) ;
		
		if(done.find(e->second) == done.end())
		{
			if(freeIds.empty())
			{
				done[e->second] = startid ;
				usedId = startid ;
				startid++ ;
			}
			else
			{
				done[e->second] = *freeIds.begin() ;
				usedId = *freeIds.begin() ;
				freeIds.erase(freeIds.begin());
			}
		}
		else
		{
			usedId = done[e->second] ;
		}
		f = shapefunc[1]* ( s - vm.eval ( s, Point ( 0,0 ) ) ) ;
		f.setPoint ( e->second ) ;
		f.setDofID ( usedId ) ;
		e->setEnrichment ( f , static_cast<SegmentedLine *>(this) ) ;

		if(done.find(e->third) == done.end())
		{
			if(freeIds.empty())
			{
				done[e->third] = startid ;
				usedId = startid ;
				startid++ ;
			}
			else
			{
				done[e->third] = *freeIds.begin() ;
				usedId = *freeIds.begin() ;
				freeIds.erase(freeIds.begin());
			}
		}
		else
		{
			usedId = done[e->third] ;
		}
		f = shapefunc[2] * ( s - vm.eval ( s, Point ( 1,0 ) ) ) ;

		f.setPoint ( e->third ) ;
		f.setDofID ( usedId ) ;
		e->setEnrichment ( f , static_cast<SegmentedLine *>(this) ) ;

		std::vector<DelaunayTriangle *> toEnrichAlso ;
		for ( size_t j = 0 ; j < e->neighbourhood.size() ; j++ )
		{
			if ( e->getNeighbourhood(j)->isAlive() 
				&& !line->intersects(static_cast<Triangle *>(e->getNeighbourhood(j))))
				toEnrichAlso.push_back ( e->getNeighbourhood(j) ) ;
		}


		for ( size_t j = 0 ; j < toEnrichAlso.size() ; j++ )
		{
			if(!toEnrichAlso[j]->enrichmentUpdated)
				toEnrichAlso[j]->clearEnrichment( static_cast<SegmentedLine *>(this)) ;
			toEnrichAlso[j]->enrichmentUpdated = true ;
			
			DelaunayTriangle * elem = toEnrichAlso[j] ;
			enrichmentMap.insert(elem) ;

			transformed.clear() ;

			for ( size_t k = 0 ; k < intersection.size() ; k++ )
			{
				transformed.push_back ( elem->inLocalCoordinates ( intersection[k] ) ) ;
			}

			Function s_ ( intersectingSegments, elem->getXTransform(), elem->getYTransform() ) ;
			hint.clear() ;

			if ( done.find(elem->first) != done.end())
			{
				Function f = shapefunc[0]* ( s_ - vm.eval ( s_, Point ( 0,1 ) ) ) ;
				f.setPoint ( elem->first ) ;
				f.setDofID ( done[elem->first] ) ;
				elem->setEnrichment ( f , static_cast<SegmentedLine *>(this)) ;
				
			}

			if (done.find(elem->second) != done.end())
			{
				Function f = shapefunc[1]* ( s_ - vm.eval ( s_, Point ( 0,0 ) ) ) ;
				f.setPoint ( elem->second ) ;
				f.setDofID ( done[elem->second] ) ;
				elem->setEnrichment ( f , static_cast<SegmentedLine *>(this)) ;
			}
			if (done.find(elem->third) != done.end())
			{
				Function f = shapefunc[2]* ( s_ - vm.eval ( s_, Point ( 1,0 ) ) ) ;
				f.setPoint ( elem->third ) ;
				f.setDofID ( done[elem->third] ) ;
				elem->setEnrichment ( f , static_cast<SegmentedLine *>(this)) ;
			}
		}
	}
	
	for(std::map<Point *, size_t>::const_iterator i = done.begin() ; i != done.end() ; ++i)
	{
		donePoints.insert(i->first) ;
	}
	
}

void BranchedCrack::enrichSegmentedLine(size_t & startid, Mesh<DelaunayTriangle,DelaunayTreeItem> * dt, const SegmentedLine * line, const Point * target)
{
	std::vector<DelaunayTriangle *> tris = dt->getConflictingElements(line) ;
	
	std::valarray<Function> shapefunc = TriElement ( LINEAR ).getShapeFunctions() ;
	VirtualMachine vm ;
	std::map<Point *, size_t> done ;
	
	for ( size_t i = 0 ; i < tris.size() ; i++ )
	{
		if(tris[i]->in(*target))
		{
			while(tipEnrichmentMap.find(tris[i]) != tipEnrichmentMap.end() && i < tris.size())
			{
				i++ ;
			}
			if(i >= tris.size())
				break ;
			if(tipEnrichmentMap.find(tris[i]) != tipEnrichmentMap.end())
				break ;
			
			
			if(!tris[i]->enrichmentUpdated)
				tris[i]->clearEnrichment(static_cast<SegmentedLine *>(this)) ;
			
			tris[i]->enrichmentUpdated = true ;
			DelaunayTriangle *e = tris[i] ;
			enrichmentMap.insert(e) ;
			std::vector<Point> intersection ;
			

			for ( size_t j = 1 ; j < line->getBoundingPoints().size() ; j++ )
			{
				Segment test ( getBoundingPoint ( j-1 ), getBoundingPoint ( j ) ) ;

				if ( test.intersects ( static_cast<Triangle *> ( e ) ) )
				{

					//there is either one or two intersection points
					std::vector<Point> temp_intersection = test.intersection ( static_cast<Triangle *> ( e ) ) ;


					if ( temp_intersection.size() == 2 ) //we simply cross the element
					{
						if ( Segment ( temp_intersection[0], temp_intersection[1] ).vector() * test.vector() < 0 )
						{
							std::swap ( temp_intersection[0], temp_intersection[1] ) ;
						}
						intersection = temp_intersection ;

					}
					else //then there are kinks
					{
						intersection.push_back ( temp_intersection[0] ) ;
						intersection.push_back ( getBoundingPoint ( j ) ) ;
						for ( size_t k = j ; k < this->getBoundingPoints().size()-1 ; k++ )
						{
							Segment test_for_kink ( getBoundingPoint ( k ), getBoundingPoint ( k+1 ) ) ;

							if ( !test_for_kink.intersects ( static_cast<Triangle *> ( e ) ) )
							{
								intersection.push_back ( getBoundingPoint ( k+1 ) ) ;
							}
							else
							{
								intersection.push_back ( test_for_kink.intersection ( static_cast<Triangle *> ( e ) ) [0] ) ;
							}
						}
					}
				}
			}


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
			
			std::vector<Segment> intersectingSegments ;
			for ( size_t j = 1 ; j < this->getBoundingPoints().size() ; j++ )
			{
				intersectingSegments.push_back ( Segment ( getBoundingPoint ( j-1 ), getBoundingPoint ( j ) ) ) ;
			}
			Function s ( intersectingSegments, e->getXTransform(), e->getYTransform() ) ;

			int usedId = 0 ;
			if(done.find(e->first) == done.end())
			{
				if(freeIds.empty())
				{
					done[e->first] = startid ;
					usedId = startid ;
					startid++ ;
				}
				else
				{
					done[e->first] = *freeIds.begin() ;
					usedId = *freeIds.begin() ;
					freeIds.erase(freeIds.begin());
				}
			}
			else
			{
				usedId = done[e->first] ;
			}

			Function f = shapefunc[0]* ( s - vm.eval ( s, Point ( 0,1 ) ) ) ;
			f.setIntegrationHint ( hint ) ;
			f.setPoint ( e->first ) ;
			f.setDofID ( usedId ) ;
			e->setEnrichment ( f , static_cast<SegmentedLine *>(this) ) ;
			
			if(done.find(e->second) == done.end())
			{
				if(freeIds.empty())
				{
					done[e->second] = startid ;
					usedId = startid ;
					startid++ ;
				}
				else
				{
					done[e->second] = *freeIds.begin() ;
					usedId = *freeIds.begin() ;
					freeIds.erase(freeIds.begin());
				}
			}
			else
			{
				usedId = done[e->second] ;
			}
			f = shapefunc[1]* ( s - vm.eval ( s, Point ( 0,0 ) ) ) ;
			f.setPoint ( e->second ) ;
			f.setDofID ( usedId ) ;
			e->setEnrichment ( f , static_cast<SegmentedLine *>(this) ) ;

			if(done.find(e->third) == done.end())
			{
				if(freeIds.empty())
				{
					done[e->third] = startid ;
					usedId = startid ;
					startid++ ;
				}
				else
				{
					done[e->third] = *freeIds.begin() ;
					usedId = *freeIds.begin() ;
					freeIds.erase(freeIds.begin());
				}
			}
			else
			{
				usedId = done[e->third] ;
			}
			f = shapefunc[2] * ( s - vm.eval ( s, Point ( 1,0 ) ) ) ;

			f.setPoint ( e->third ) ;
			f.setDofID ( usedId ) ;
			e->setEnrichment ( f , static_cast<SegmentedLine *>(this) ) ;

			std::vector<DelaunayTriangle *> toEnrichAlso ;
			for ( size_t j = 0 ; j < e->neighbourhood.size() ; j++ )
			{
				if ( e->getNeighbourhood(j)->isAlive() 
					&& !line->intersects(static_cast<Triangle *>(e->getNeighbourhood(j))))
					toEnrichAlso.push_back ( e->getNeighbourhood(j) ) ;
			}


			for ( size_t j = 0 ; j < toEnrichAlso.size() ; j++ )
			{
				if(!toEnrichAlso[j]->enrichmentUpdated)
					toEnrichAlso[j]->clearEnrichment( static_cast<SegmentedLine *>(this)) ;
				toEnrichAlso[j]->enrichmentUpdated = true ;
				
				DelaunayTriangle * elem = toEnrichAlso[j] ;
				enrichmentMap.insert(elem) ;

				transformed.clear() ;

				for ( size_t k = 0 ; k < intersection.size() ; k++ )
				{
					transformed.push_back ( elem->inLocalCoordinates ( intersection[k] ) ) ;
				}

				Function s_ ( intersectingSegments, elem->getXTransform(), elem->getYTransform() ) ;
				hint.clear() ;

				if ( done.find(elem->first) != done.end())
				{
					Function f = shapefunc[0]* ( s_ - vm.eval ( s_, Point ( 0,1 ) ) ) ;
					f.setPoint ( elem->first ) ;
					f.setDofID ( done[elem->first] ) ;
					elem->setEnrichment ( f , static_cast<SegmentedLine *>(this)) ;
					
				}

				if (done.find(elem->second) != done.end())
				{
					Function f = shapefunc[1]* ( s_ - vm.eval ( s_, Point ( 0,0 ) ) ) ;
					f.setPoint ( elem->second ) ;
					f.setDofID ( done[elem->second] ) ;
					elem->setEnrichment ( f , static_cast<SegmentedLine *>(this)) ;
				}
				if (done.find(elem->third) != done.end())
				{
					Function f = shapefunc[2]* ( s_ - vm.eval ( s_, Point ( 1,0 ) ) ) ;
					f.setPoint ( elem->third ) ;
					f.setDofID ( done[elem->third] ) ;
					elem->setEnrichment ( f , static_cast<SegmentedLine *>(this)) ;
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
	for(std::set<DelaunayTriangle *>::iterator i = enrichmentMap.begin() ; i != enrichmentMap.end() ; ++i)
	{
		toClear.push_back(*i) ;
	}
	
	for(std::vector<DelaunayTriangle *>::iterator i = toClear.begin() ; i != toClear.end() ; ++i)
	{
		std::vector<size_t> idpts = (*i)->clearEnrichment(getPrimitive()) ;
		for(size_t j = 0 ; j < idpts.size() ; j++)
			freeIds.insert(idpts[j]) ;
		(*i)->enrichmentUpdated = true;
	}
	
	enrichmentMap.clear() ;
	tipEnrichmentMap.clear() ;
	enrichTips(counter, dtree) ;
	enrichForks(counter, dtree) ;
	enrichBranches(counter, dtree) ;
}

void BranchedCrack::computeCenter()
{
	for(size_t i = 0 ; i < tips.size() ; i++)
	{
		SegmentedLine::center += *tips[i].first*(1./tips.size()) ;
	}
}

std::vector<DelaunayTriangle*> BranchedCrack::getElements(Mesh<DelaunayTriangle,DelaunayTreeItem>* dt)
{
	
	std::vector<DelaunayTriangle*> ret ;

	for(size_t i = 0 ; i < branches.size() ; i++)
	{
		std::vector<DelaunayTriangle*> tris = dt->getConflictingElements(branches[i]) ;
		ret.insert(ret.end(), tris.begin(), tris.end()) ;
	}
	
	std::sort(ret.begin(), ret.end()) ;
	ret.erase(std::unique(ret.begin(), ret.end()), ret.end()) ;
	
	return ret ;
}

std::vector<DelaunayTetrahedron*> BranchedCrack::getElements(Mesh<DelaunayTetrahedron, DelaunayTreeItem3D>*)
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

Point* BranchedCrack::pointAfter(size_t)
{
	return NULL ;
}

std::vector<Mu::Geometry*> BranchedCrack::getRefinementZones(size_t level ) const
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
			std::cout << "("<< branches[i]->getBoundingPoint(j).x << ", " << branches[i]->getBoundingPoint(j).y << ")" << std::endl ;
		}
	}

	
	std::cout << " tips : " << tips.size() << std::endl;
	for(size_t i = 0 ; i < tips.size() ; i++)
	{
		std::cout << "tip " << i << " : "<< "("<< tips[i].first->x << ", " << tips[i].first->y << ")"<< std::endl ;
	}
	std::cout << " forks : " << forks.size() << std::endl;
	for(size_t i = 0 ; i < forks.size() ; i++)
	{
		std::cout << "fork " << i << std::endl ;
		
		for(size_t j = 0 ; j < forks[i]->getBoundingPoints().size() ; j++)
		{
			std::cout << "("<< forks[i]->getBoundingPoint(j).x << ", " << forks[i]->getBoundingPoint(j).y << ")" << std::endl ;
		}
		
	}
	std::cout << std::endl ;

}

void BranchedCrack::sample(size_t)
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

void BranchedCrack::step(double dt, Vector*, Mesh<DelaunayTriangle,DelaunayTreeItem> * dtree)
{
	changed = false ;
	std::vector<Point *> tipsToGrow ; 
	std::vector<double> angles ; 
	for(size_t i = 0 ; i < tips.size() ; i++)
	{
		std::pair<double, double> energy  = computeJIntegralAtTip(tips[i], dtree);
		if(energy.first*energy.first + energy.second*energy.second > 0)
		{
			tipsToGrow.push_back(tips[i].first) ;
			changed = true ;
			angles.push_back(propagationAngleFromTip(tips[i], dtree)) ;
		}
	}
	
	for(size_t i = 0 ; i < tipsToGrow.size() ; i++)
	{
		grow(tipsToGrow[i], new Point(tipsToGrow[i]->x+enrichementRadius*.5*cos(angles[i]),
		                              tipsToGrow[i]->y+enrichementRadius*.5* sin(angles[i]))) ;
	}
}

void BranchedCrack::snap(Mesh<DelaunayTriangle,DelaunayTreeItem>*)
{
	
}

bool BranchedCrack::moved() const
{
	return changed ;
}
