// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2006-2007
//
// Copyright: See COPYING file that comes with this distribution
//

#include "crack.h"

#include <fstream>

using namespace Mu ;

BranchedCrack::BranchedCrack(Feature *father, Point * a, Point * b) : EnrichmentFeature(father), SegmentedLine(std::valarray<Point * >(0))
{
	std::valarray<Point * > newBranch ( 2 ) ;
	newBranch[0] = a ;
	newBranch[1] = b ;
	branches.push_back ( new SegmentedLine ( newBranch ) ) ;

	if(father->in(*a))
		tips.push_back(std::make_pair(a, atan2(a->y-b->y, a->x-b->x))) ;
	if(father->in(*b))
		tips.push_back(std::make_pair(b, atan2(b->y-a->y, b->x-a->x))) ;

	changed = false ;
	
}

bool operator ==(const std::pair<Mu::Point*, double> & a, const Mu::Point* b)
{
	return a.first == b ;
}

BranchedCrack::BranchedCrack(Point * a, Point * b) : EnrichmentFeature(NULL), SegmentedLine(std::valarray<Point * >(0))
{
	std::valarray<Point * > newBranch ( 2 ) ;
	newBranch[0] = a ;
	newBranch[1] = b ;
	branches.push_back ( new SegmentedLine ( newBranch ) ) ;

	tips.push_back(std::make_pair(a, atan2(a->y-b->y, a->x-b->x))) ;
	tips.push_back(std::make_pair(b, atan2(b->y-a->y, b->x-a->x))) ;

	changed = false ;
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
	forkComponent0[1] = newTip1 ;

	std::valarray<Point *> forkComponent1 ( 3 ) ;
	std::valarray<Point *> forkComponent2 ( 3 ) ;

	std::valarray<Point *> newBP ( branchToExtend->getBoundingPoints().size() +1 ) ;
	if ( fromHead )
	{
		std::copy ( &branchToExtend->getBoundingPoints() [0],
		            &branchToExtend->getBoundingPoints() [newBP.size()-1],
		            &newBP[1] ) ;
		newBP[0] = newTip0 ;
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

void BranchedCrack::enrichTips(size_t & startid, DelaunayTree * dt)
{
	for(size_t i =  0 ; i < tips.size() ; i++)
	{
		enrichTip(startid, dt, tips[i], 0.) ;
	}
}

void BranchedCrack::enrichTip(size_t & startid, DelaunayTree * dt, const std::pair<Point *, double> & tip, double angle)
{
	Circle epsilon(enrichementRadius, *tip.first) ;
	std::vector<DelaunayTriangle *> triangles = dt->conflicts(&epsilon) ;
	std::map<Point *, size_t> done ;
	VirtualMachine vm ;
	
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
		std::vector<Point> hint ;
		Line crossing ( *tip.first, Point(cos(angle), sin(angle)) ) ;
		std::vector<Point > completeIntersection = crossing.intersection ( static_cast<Triangle *>(triangles[i]) ) ;
		
		if(completeIntersection.size() == 2)
		{
			Point intersectionTransformed = triangles[i]->inLocalCoordinates ( completeIntersection[0] ) ;
			Point intersectionBisTransformed = triangles[i]->inLocalCoordinates ( completeIntersection[1] ) ;
			Point singularityTransformed = triangles[i]->inLocalCoordinates ( *tip.first );
			
			hint.push_back ( intersectionTransformed ) ;
			hint.push_back ( intersectionBisTransformed ) ;
			hint.push_back ( singularityTransformed ) ;
		}
		
		
		Function x = triangles[i]->getXTransform() ;
		Function y = triangles[i]->getYTransform() ;
		double rotatedSingularityX = tip.first->x*cos ( angle ) + tip.first->y*sin ( angle ) ;
		double rotatedSingularityY = -tip.first->x*sin ( angle ) + tip.first->y*cos ( angle ) ;
		Function rotatedX = x*cos ( angle ) + y*sin ( angle ) ;
		Function rotatedY = x*sin ( -angle ) + y*cos ( angle ) ;
		Function x_ = x - tip.first->x ;
		Function y_ = y - tip.first->y ;
		Function theta = f_atan2 ( rotatedY-rotatedSingularityY, rotatedX-rotatedSingularityX );
		Function r = f_sqrt ( ( x_^2 ) + ( y_^2 ) );
		
		Function f0 = f_sqrt ( r ) *f_sin ( theta/2 );
		Function f1 = f_sqrt ( r ) *f_cos ( theta/2 );
		Function f2 = f_sqrt ( r ) *f_sin ( theta/2 ) *f_cos ( theta );
		Function f3 = f_sqrt ( r ) *f_cos ( theta/2 ) *f_cos ( theta );
		
		std::vector<Point *> currentPoint ;
		currentPoint.push_back(triangles[i]->first) ;
		currentPoint.push_back(triangles[i]->second) ;
		currentPoint.push_back(triangles[i]->third) ;
		
		for(size_t p = 0 ; p < 3 ;p++)
		{
			if(epsilon.in(*currentPoint[p]))
			{
				int usedId = 0 ;
				if(done.find(currentPoint[p]) == done.end())
				{
					done[currentPoint[p]] = startid ;
					usedId = startid ;
					startid+=4 ;
				}
				else
				{
					usedId = done[currentPoint[p]] ;
				}

				Function f = shapefunc[p]* ( f0 - vm.eval ( f0, pointLocal[p] ) ) ;
				f.setIntegrationHint ( hint ) ;
				f.setPointID ( currentPoint[p]->id ) ;
				f.setDofID ( usedId ) ;
				triangles[i]->setEnrichment (  f  ) ;
				
				f = shapefunc[p]* ( f1 - vm.eval ( f1, pointLocal[p] ) ) ;
				f.setIntegrationHint ( hint ) ;
				f.setPointID ( currentPoint[p]->id ) ;
				f.setDofID ( usedId+1 ) ;
				triangles[i]->setEnrichment ( f  ) ;
				
				f = shapefunc[p]* ( f2 - vm.eval ( f2, pointLocal[p] ) ) ;
				f.setIntegrationHint ( hint ) ;
				f.setPointID ( currentPoint[p]->id ) ;
				f.setDofID ( usedId+2 ) ;
				triangles[i]->setEnrichment ( f  ) ;
				
				f = shapefunc[p]* ( f3 - vm.eval ( f3, pointLocal[p] ) ) ;
				f.setIntegrationHint ( hint ) ;
				f.setPointID (currentPoint[p]->id ) ;
				f.setDofID ( usedId+3 ) ;
				triangles[i]->setEnrichment ( f  ) ;
			}
		}
		
	}

	for(std::map<Point *, size_t>::const_iterator i = done.begin() ; i != done.end() ; ++i)
	{
		donePoints.insert(i->first) ;
	}
}

void BranchedCrack::enrichForks(size_t & startid, DelaunayTree * dt)
{
	for(size_t i =  0 ; i < forks.size() ; i++)
	{
		enrichSegmentedLine(startid, dt, forks[i]) ;
	}
}

void BranchedCrack::enrichBranches(size_t & startid, DelaunayTree * dt)
{
	for(size_t i =  0 ; i < branches.size() ; i++)
	{
		enrichSegmentedLine(startid, dt, branches[i]) ;
	}
}

void BranchedCrack::enrichSegmentedLine(size_t & startid, DelaunayTree * dt, const SegmentedLine * line)
{
	std::vector<DelaunayTriangle *> tris = dt->conflicts(line) ;
	std::valarray<Function> shapefunc = TriElement ( LINEAR ).getShapeFunctions() ;
	VirtualMachine vm ;
	std::map<Point *, size_t> done ;
	
	for ( size_t i = 0 ; i < tris.size() ; i++ )
	{
		DelaunayTriangle *e = tris[i] ;

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

		if ( transformed.size() >2 )
		{
			hint.push_back ( transformed[0] ) ;
			for ( size_t k = 1 ; k < transformed.size()-1 ; k++ )
			{
				hint.push_back ( transformed[k] ) ;
			}
			hint.push_back ( transformed[transformed.size()-1] ) ;
		}
		else
		{
			for ( size_t k = 0 ; k < transformed.size() ; k++ )
			{
				hint.push_back ( transformed[k] ) ;
			}
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
			done[e->first] = startid ;
			usedId = startid ;
			startid++ ;
		}
		else
		{
			usedId = done[e->first] ;
		}

		Function f = shapefunc[0]* ( s - vm.eval ( s, Point ( 0,1 ) ) ) ;
		f.setIntegrationHint ( hint ) ;
		f.setPointID ( e->first->id ) ;
		f.setDofID ( usedId ) ;
		e->setEnrichment ( f  ) ;

		if(done.find(e->first) == done.end())
		{
			done[e->first] = startid ;
			usedId = startid ;
			startid++ ;
		}
		else
		{
			usedId = done[e->first] ;
		}
		f = shapefunc[1]* ( s - vm.eval ( s, Point ( 0,0 ) ) ) ;
		f.setIntegrationHint ( hint ) ;
		f.setPointID ( e->second->id ) ;
		f.setDofID ( usedId ) ;
		e->setEnrichment ( f  ) ;

		if(done.find(e->first) == done.end())
		{
			done[e->first] = startid ;
			usedId = startid ;
			startid++ ;
		}
		else
		{
			usedId = done[e->first] ;
		}
		f = shapefunc[2]* ( s - vm.eval ( s, Point ( 1,0 ) ) ) ;

		f.setIntegrationHint ( hint ) ;
		f.setPointID ( e->third->id ) ;
		f.setDofID ( usedId ) ;
		e->setEnrichment ( f  ) ;

		std::vector<DelaunayTriangle *> toEnrichAlso ;
		for ( size_t j = 0 ; j < e->neighbourhood.size() ; j++ )
		{
			if ( e->getNeighbourhood(j)->isAlive() && !line->intersects(static_cast<Triangle *>(e)))
				toEnrichAlso.push_back ( e->getNeighbourhood(j) ) ;
		}


		for ( size_t j = 0 ; j < toEnrichAlso.size() ; j++ )
		{
			DelaunayTriangle * elem = toEnrichAlso[j] ;

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
				f.setIntegrationHint ( hint ) ;
				f.setPointID ( elem->first->id ) ;
				f.setDofID ( done[elem->first] ) ;
				elem->setEnrichment ( f ) ;

			}

			if (done.find(elem->second) != done.end())
			{
				Function f = shapefunc[1]* ( s_ - vm.eval ( s_, Point ( 0,0 ) ) ) ;
				f.setIntegrationHint ( hint ) ;
				f.setPointID ( elem->second->id ) ;
				f.setDofID ( done[elem->second] ) ;
				elem->setEnrichment ( f ) ;

			}
			if (done.find(elem->third) != done.end())
			{
				Function f = shapefunc[2]* ( s_ - vm.eval ( s_, Point ( 1,0 ) ) ) ;
				f.setIntegrationHint ( hint ) ;
				f.setPointID ( elem->third->id ) ;
				f.setDofID ( done[elem->third] ) ;
				elem->setEnrichment ( f ) ;
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

void BranchedCrack::enrich(size_t & counter,  DelaunayTree * dtree)
{
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

std::vector<DelaunayTriangle*> BranchedCrack::getTriangles(DelaunayTree*)
{
	return std::vector<DelaunayTriangle*>() ;
}

std::vector<DelaunayTetrahedron*> BranchedCrack::getTetrahedrons(DelaunayTree_3D*)
{
	return std::vector<DelaunayTetrahedron*>() ;
}

bool BranchedCrack::interacts(Feature* f) const
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
	std::cout << " forks : " << forks.size()/3 << std::endl;
	for(size_t i = 0 ; i < forks.size()/3 ; i++)
	{
		std::cout << "fork " << i << std::endl ;
		
		for(size_t j = 0 ; j < forks[i*3]->getBoundingPoints().size() ; j++)
		{
			std::cout << "("<< forks[i*3]->getBoundingPoint(j).x << ", " << forks[i*3]->getBoundingPoint(j).y << ")" << std::endl ;
		}
		for(size_t j = 0 ; j < forks[i*3+1]->getBoundingPoints().size() ; j++)
		{
			std::cout << "("<< forks[i*3+1]->getBoundingPoint(j).x << ", " << forks[i*3+1]->getBoundingPoint(j).y << ")" << std::endl ;
		}
		for(size_t j = 0 ; j < forks[i*3+2]->getBoundingPoints().size() ; j++)
		{
			std::cout << "("<< forks[i*3+2]->getBoundingPoint(j).x << ", " << forks[i*3+2]->getBoundingPoint(j).y << ")" << std::endl ;
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
	return enrichmentMap.find(std::make_pair(tri, tri->first)) != enrichmentMap.end()
		|| enrichmentMap.find(std::make_pair(tri, tri->second)) != enrichmentMap.end()
		|| enrichmentMap.find(std::make_pair(tri, tri->third)) != enrichmentMap.end() ;
}

void BranchedCrack::step(double, Vector*, const DelaunayTree*)
{
	
}

void BranchedCrack::snap(DelaunayTree*)
{
	
}

bool BranchedCrack::moved() const
{
	return changed ;
}


std::vector<Geometry *> Crack::getRefinementZones ( size_t level ) const
{
	std::vector<Geometry *> ret ;

	if ( level > 0 )
	{
		for ( size_t i = 0 ; i < boundingPoints.size() ; i++ )
		{
			ret.push_back ( new Circle ( 0.2, getBoundingPoint ( i ) ) ) ;
		}
	}
	if ( level > 1 )
	{
		for ( size_t i = 0 ; i < boundingPoints.size() ; i++ )
		{
			ret.push_back ( new Circle ( 0.15, getBoundingPoint ( i ) ) ) ;
		}
	}
	if ( level > 3 )
	{
		for ( size_t i = 0 ; i < boundingPoints.size() ; i++ )
		{
			ret.push_back ( new Circle ( 0.1, getBoundingPoint ( i ) ) ) ;
		}
	}
	return ret ;
}

std::vector<DelaunayTriangle *> Crack::getTriangles ( DelaunayTree * dt )
{
	std::vector<DelaunayTriangle *> vec  ;
	std::vector<DelaunayTriangle *> visited ;

	vec = dt->conflicts(this->getPrimitive()) ;

	std::vector<DelaunayTriangle *> inhead  = dt->conflicts ( boundary ) ;
	for ( size_t i = 0 ; i< inhead.size() ; i++ )
	{
		vec.push_back ( inhead[i] ) ;
	}

	inhead  = dt->conflicts ( boundary2 ) ;

	for ( size_t i = 0 ; i< inhead.size() ; i++ )
	{
		vec.push_back ( inhead[i] ) ;
	}

	std::sort ( vec.begin(), vec.end() ) ;
	std::vector<DelaunayTriangle *>::iterator e = std::unique ( vec.begin(), vec.end() );
	vec.erase ( e, vec.end() ) ;
	return vec ;
}

std::vector<DelaunayTriangle *> Crack::getIntersectingTriangles ( DelaunayTree * dt )
{
	return dt->conflicts(this->getPrimitive()) ;
}

Crack::Crack ( Feature * father, const std::valarray<Point *> & points, double radius ) : EnrichmentFeature ( father ), SegmentedLine ( points )
{
  infRad = radius ;// enrichment radius
  
  
  this->isEnrichmentFeature = true ;
  
  double maxr = 0.25*dist ( getHead(),getTail() ) ;
  if ( infRad > maxr )
    infRad = maxr ;
  
  this->boundary = new Circle ( infRad, getHead() ) ;
  this->boundary2 = new Circle ( infRad, getTail() ) ;
  changed = true ;
  // set increment in crack growth for teach step
  stepLength = 0.25*infRad ;
  criticalJ = 0.0 ;
}

Crack::Crack ( const std::valarray<Point *> & points, double radius ) : EnrichmentFeature ( NULL ), SegmentedLine ( points )
{
	infRad = radius ;
	this->isEnrichmentFeature = true ;

	double maxr = 0.25*dist ( getHead(),getTail() ) ;
	if ( infRad > maxr )
		infRad = maxr ;

	this->boundary = new Circle ( infRad, getHead() ) ;
	this->boundary2 = new Circle ( infRad, getTail() ) ;
	changed = true ;
	stepLength = .25*infRad ;
	criticalJ = 0.0 ;
}

/**
SB New constructor 
@param radius enrichment radius for the crack tips
@param gfac number such that the growth increment is gfac*enrichment radius
@param critJ critical value of the J integral for propagation (Fracture toughness, in fact this is a material param...)
**/
Crack::Crack ( Feature * father, const std::valarray<Point *> & points, double radius, double gfac, double critJ ) : EnrichmentFeature ( father ), SegmentedLine ( points )
{
  infRad = radius ;// enrichment radius
  
  
  this->isEnrichmentFeature = true ;
  
  double maxr = 0.25*dist ( getHead(),getTail() ) ;
  if ( infRad > maxr )
    infRad = maxr ;
  
  this->boundary = new Circle ( infRad, getHead() ) ;
  this->boundary2 = new Circle ( infRad, getTail() ) ;
  changed = true ;
  // set increment in crack growth for teach step
  stepLength = gfac*infRad ;
  criticalJ = critJ ;
}

/**
SB New constructor 
@param radius enrichment radius for the crack tips
@param gfac number such that the growth increment is gfac*enrichment radius
@param critJ critical value of the J integral for propagation (Fracture toughness, in fact this is a material param...)
**/
Crack::Crack ( const std::valarray<Point *> & points, double radius, double gfac, double critJ ) : EnrichmentFeature ( NULL ), SegmentedLine ( points )
{
	infRad = radius ;
	this->isEnrichmentFeature = true ;

	double maxr = 0.25*dist ( getHead(),getTail() ) ;
	if ( infRad > maxr )
		infRad = maxr ;

	this->boundary = new Circle ( infRad, getHead() ) ;
	this->boundary2 = new Circle ( infRad, getTail() ) ;
	changed = true ;
	stepLength = gfac*infRad ;
	criticalJ = critJ ;
}


void Crack::setInfluenceRadius ( double r )
{
	
	delete this->boundary ;
	delete this->boundary2 ;
	double maxr = 0.25*dist ( getHead(),getTail() ) ;
	if ( r > maxr )
		r = maxr ;
	this->infRad = r ;
	stepLength = .25*infRad ;
	this->boundary = new Circle ( r, getHead() ) ;
	this->boundary2 = new Circle ( r, getTail() ) ;
}

/**
SB to set the critical J, influence radius and crack growth increment per step
@param radius enrichment radius for the crack tips
@param gfac number such that the growth increment is gfac*enrichment radius
@param critJ critical value of the J integral for propagation (Fracture toughness, in fact this is a material param...)
**/
void Crack::setParams ( double r, double gfac, double critJ)
{
	
	delete this->boundary ;
	delete this->boundary2 ;
	double maxr = 0.25*dist ( getHead(),getTail() ) ;
	if ( r > maxr )
		r = maxr ;
	this->infRad = r ;
	stepLength = gfac*infRad ;
	this->boundary = new Circle ( r, getHead() ) ;
	this->boundary2 = new Circle ( r, getTail() ) ;
	criticalJ = critJ;
}

bool Crack::enrichmentTarget ( DelaunayTriangle * t )
{
	if ( t==NULL )
		return false ;
	for ( size_t i = 0 ; i < boundingPoints.size()-1 ; i++ )
	{
		Segment s ( getBoundingPoint ( i ), getBoundingPoint ( i+1 ) ) ;
		if ( s.intersects ( dynamic_cast<Triangle *> ( t ) ) )
			return true ;
	}

	return false ;
}

void Crack::enrich ( size_t & counter, DelaunayTree * dtree )
{
	std::vector<DelaunayTriangle *> tris = this->getTriangles ( dtree ) ;
	if ( tris.size() < 2 )
	{
		std::cout << "crack too small, cowardly discarding" << std::endl ;
		return ;
	}

	map.update ( &tris, this, counter ) ;

	std::valarray<Function> shapefunc = TriElement ( LINEAR ).getShapeFunctions() ;
	VirtualMachine vm ;

	//neg head -- wrong
	//no neg -- wrong
	//swaptail -- correct
	//atan != atan2 -- wrong
	double headAngle = atan2 ( getHead()->y-getBoundingPoint ( 1 ).y, 
	                           getHead()->x-getBoundingPoint ( 1 ).x ) ;
	double tailAngle = atan2 ( getTail()->y-getBoundingPoint ( getBoundingPoints().size()-2 ).y,
	                           getTail()->x-getBoundingPoint ( getBoundingPoints().size()-2 ).x) ;

	std::vector<Segment> intersectingSegments ;
	for ( size_t j = 1 ; j < this->getBoundingPoints().size() ; j++ )
	{
		intersectingSegments.push_back ( Segment ( getBoundingPoint ( j-1 ), getBoundingPoint ( j ) ) ) ;
	}

	for ( size_t i = 0 ; i < tris.size() ; i++ )
	{
		DelaunayTriangle *e = tris[i] ;

		std::vector<Point> intersection ;

		for ( size_t j = 1 ; j < this->getBoundingPoints().size() ; j++ )
		{
			Segment test ( getBoundingPoint ( j-1 ), getBoundingPoint ( j ) ) ;

			if ( test.intersects ( dynamic_cast<Triangle *> ( e ) ) )
			{

				//there is either one or two intersection points
				std::vector<Point> temp_intersection = test.intersection ( dynamic_cast<Triangle *> ( e ) ) ;


				if ( temp_intersection.size() == 2 ) //we simply cross the element
				{
					if ( Segment ( temp_intersection[0], temp_intersection[1] ).vector() * test.vector() < 0 )
					{
						std::swap ( temp_intersection[0], temp_intersection[1] ) ;
					}
					intersection = temp_intersection ;

				}
				else if ( e->in ( *getHead() ) || e->in ( *getTail() ) ) // First : check for singularity.
				{
					intersection.push_back ( temp_intersection[0] ) ;
				}
				else //then there are kinks
				{
					intersection.push_back ( temp_intersection[0] ) ;
					intersection.push_back ( getBoundingPoint ( j ) ) ;
					for ( size_t k = j ; k < this->getBoundingPoints().size()-1 ; k++ )
					{
						Segment test_for_kink ( getBoundingPoint ( k ), getBoundingPoint ( k+1 ) ) ;

						if ( !test_for_kink.intersects ( dynamic_cast<Triangle *> ( e ) ) )
						{
							intersection.push_back ( getBoundingPoint ( k+1 ) ) ;
						}
						else
						{
							intersection.push_back ( test_for_kink.intersection ( dynamic_cast<Triangle *> ( e ) ) [0] ) ;
						}
					}
				}
			}
		}

		if ( intersection.size() == 1 )
		{
// 			e->setNonLinearBehaviour( new TwoDCohesiveForces(e, dynamic_cast<SegmentedLine *>(this)) ) ;

			Point singularity ;

			double angle = 0 ;

			if ( e->in ( *getHead() ) )
			{
				singularity = *getHead() ;
				intersection.push_back ( singularity ) ;
				angle = headAngle ;
			}
			else
			{
				singularity = *getTail() ;
				intersection.push_back ( singularity ) ;
				angle = tailAngle ;
			}

			Line crossing ( intersection[0], intersection[0]-singularity ) ;
			std::vector<Point > completeIntersection = crossing.intersection ( dynamic_cast<Triangle *> ( e ) ) ;
			Point intersectionTransformed = e->inLocalCoordinates ( completeIntersection[0] ) ;
			Point intersectionBisTransformed = e->inLocalCoordinates ( completeIntersection[1] ) ;

			Point singularityTransformed = e->inLocalCoordinates ( singularity );
			std::vector<Point> hint ;


			hint.push_back ( intersectionTransformed ) ;
			hint.push_back ( intersectionBisTransformed ) ;
			hint.push_back ( singularityTransformed ) ;

			Function x = e->getXTransform() ;
			Function y = e->getYTransform() ;
			double rotatedSingularityX = singularity.x*cos ( angle ) + singularity.y*sin ( angle ) ;
			double rotatedSingularityY = -singularity.x*sin ( angle ) + singularity.y*cos ( angle ) ;
			Function rotatedX = x*cos ( angle ) + y*sin ( angle ) ;
			Function rotatedY = x*sin ( -angle ) + y*cos ( angle ) ;
			Function x_ = x - singularity.x ;
			Function y_ = y - singularity.y ;
			Function theta = f_atan2 ( rotatedY-rotatedSingularityY, rotatedX-rotatedSingularityX );
			Function r = f_sqrt ( ( x_^2 ) + ( y_^2 ) );

			Function f0 = f_sqrt ( r ) *f_sin ( theta/2 );
			Function f1 = f_sqrt ( r ) *f_cos ( theta/2 );
			Function f2 = f_sqrt ( r ) *f_sin ( theta/2 ) *f_cos ( theta );
			Function f3 = f_sqrt ( r ) *f_cos ( theta/2 ) *f_cos ( theta );

			std::vector <Point> pointLocal ;
			pointLocal.push_back ( Point ( 0,1 ) ) ;
			pointLocal.push_back ( Point ( 0,0 ) ) ;
			pointLocal.push_back ( Point ( 1,0 ) ) ;
			std::vector <Point *> pointGlobal ;
			pointGlobal.push_back ( e->first ) ;
			pointGlobal.push_back ( e->second ) ;
			pointGlobal.push_back ( e->third ) ;

// 			Triangle parent ;
// 			for(size_t k  = 0 ; k < 100 ;  k++)
// 			{
// 				for(size_t l  = 0 ; l < 100 ;  l++)
// 				{
// 					if(parent.in(Point((double)k/99., (double)l/99.)))
// 						std::cout << vm.deval(f1, XI, (double)k/99., (double)l/99.) << "   "<< std::flush ;
// 					else
// 						std::cout << "0   "<< std::flush ;
// 				}
// 				std::cout << std::endl ;
// 			}
// 			Function f = shapefunc[0]*(f0 - vm.eval(f0, pointLocal[0])) ;

// 			for(size_t i = 0 ; i < 80 ; i++)
// 			{
// 				for(size_t j = 0 ; j < 80 ; j++)
// 				{
// 					std::cout << vm.deval(f0,ETA, (double)i/80,(double)j/80) << "   " << std::flush ;
// 				}
// 				std::cout << std::endl ;
// 			}



			for ( size_t j = 0 ; j < 3 ; j++ )
			{
				if ( map.getEnrichment ( pointGlobal[j]->id ).getType() == SINGULAR_ENRICHMENT )
				{
					std::vector<size_t> ids = map.getEnrichment ( pointGlobal[j]->id ).getID() ;

					Function f = shapefunc[j]* ( f0 - vm.eval ( f0, pointLocal[j] ) ) ;
					f.setIntegrationHint ( hint ) ;
					f.setPointID ( pointGlobal[j]->id ) ;
					f.setDofID ( ids[0] ) ;
					e->setEnrichment (f) ;

					f = shapefunc[j]* ( f1 - vm.eval ( f1, pointLocal[j] ) ) ;
					f.setIntegrationHint ( hint ) ;
					f.setPointID ( pointGlobal[j]->id ) ;
					f.setDofID ( ids[1] ) ;
					e->setEnrichment (f) ;

					f = shapefunc[j]* ( f2 - vm.eval ( f2, pointLocal[j] ) ) ;
					f.setIntegrationHint ( hint ) ;
					f.setPointID ( pointGlobal[j]->id ) ;
					f.setDofID ( ids[2] ) ;
					e->setEnrichment (f ) ;

					f = shapefunc[j]* ( f3 - vm.eval ( f3, pointLocal[j] ) ) ;
					f.setIntegrationHint ( hint ) ;
					f.setPointID ( pointGlobal[j]->id ) ;
					f.setDofID ( ids[3] ) ;
					e->setEnrichment (f ) ;

				}
			}

			std::vector<DelaunayTriangle *> toEnrichAlso ;
			for ( size_t j = 0 ; j < e->neighbourhood.size() ; j++ )
			{
				if ( e->getNeighbourhood(j)->isAlive() && !enrichmentTarget ( e->getNeighbourhood(j) ) )
					toEnrichAlso.push_back ( e->getNeighbourhood(j) ) ;
			}


			for ( size_t j = 0 ; j < toEnrichAlso.size() ; j++ )
			{
				DelaunayTriangle * elem = toEnrichAlso[j] ;
				if ( !toEnrichAlso[j]->in ( singularity ) )
				{

					hint.clear() ;
					intersectionTransformed = elem->inLocalCoordinates ( intersection[0] ) ;
					singularityTransformed =  elem->inLocalCoordinates ( singularity ) ;

					x = elem->getXTransform() ;
					y = elem->getYTransform() ;
					rotatedX = x*cos ( angle ) + y*sin ( angle ) ;
					rotatedY = x*sin ( -angle ) + y*cos ( angle ) ;
					x_ = x - singularity.x ;
					y_ = y - singularity.y ;
					theta = f_atan2 ( rotatedY-rotatedSingularityY, rotatedX-rotatedSingularityX );
					r = f_sqrt ( ( x_^2 ) + ( y_^2 ) );

					f0 = f_sqrt ( r ) *f_sin ( theta/2 );
					f1 = f_sqrt ( r ) *f_cos ( theta/2 );
					f2 = f_sqrt ( r ) *f_sin ( theta/2 ) *f_cos ( theta );
					f3 = f_sqrt ( r ) *f_cos ( theta/2 ) *f_cos ( theta );

					pointGlobal.clear() ;
					pointGlobal.push_back ( elem->first ) ;
					pointGlobal.push_back ( elem->second ) ;
					pointGlobal.push_back ( elem->third ) ;

					for ( size_t j = 0 ; j < 3 ; j++ )
					{
						if ( map.getEnrichment ( pointGlobal[j]->id ).getType() == SINGULAR_ENRICHMENT )
						{

							std::vector<size_t> ids = map.getEnrichment ( pointGlobal[j]->id ).getID() ;

							Function f = shapefunc[j]* ( f0 - vm.eval ( f0, pointLocal[j] ) ) ;
							f.setIntegrationHint ( hint ) ;
							f.setPointID ( pointGlobal[j]->id ) ;
							f.setDofID ( ids[0] ) ;
							elem->setEnrichment (f ) ;

							f = shapefunc[j]* ( f1 - vm.eval ( f1, pointLocal[j] ) ) ;
							f.setIntegrationHint ( hint ) ;
							f.setPointID ( pointGlobal[j]->id ) ;
							f.setDofID ( ids[1] ) ;
							elem->setEnrichment (f ) ;

							f = shapefunc[j]* ( f2 - vm.eval ( f2, pointLocal[j] ) ) ;
							f.setIntegrationHint ( hint ) ;
							f.setPointID ( pointGlobal[j]->id ) ;
							f.setDofID ( ids[2] ) ;
							elem->setEnrichment (f ) ;

							f = shapefunc[j]* ( f3 - vm.eval ( f3, pointLocal[j] ) ) ;
							f.setIntegrationHint ( hint ) ;
							f.setPointID ( pointGlobal[j]->id ) ;
							f.setDofID ( ids[3] ) ;
							elem->setEnrichment (f ) ;

						}
					}
				}
			}

		}


		if ( map.getEnrichment ( e->first->id ).getType() == SPLIT_ENRICHMENT &&
		        map.getEnrichment ( e->second->id ).getType() == SPLIT_ENRICHMENT &&
		        map.getEnrichment ( e->third->id ).getType() == SPLIT_ENRICHMENT
		   )
		{
// 			e->setNonLinearBehaviour( new TwoDCohesiveForces(e, dynamic_cast<SegmentedLine *>(this)) ) ;


			std::vector<Point> hint ;
			std::vector<Point> transformed ;

			for ( int k = 0 ; k < (int)intersection.size()-1 ; k++ )
			{
				transformed.push_back ( e->inLocalCoordinates ( intersection[k] ) ) ;
				transformed.push_back ( e->inLocalCoordinates ( intersection[k]*.5+ intersection[k+1]*.5) ) ;
				Triangle().project(&transformed.back()) ;
				
			}
			if(!intersection.empty())
				transformed.push_back ( e->inLocalCoordinates ( intersection[intersection.size()-1] ) ) ;

			for ( size_t k = 0 ; k < transformed.size() ; k++ )
			{
				hint.push_back ( transformed[k] ) ;
			}

			Function s ( intersectingSegments, e->getXTransform(), e->getYTransform() ) ;

			std::vector<size_t> ids = map.getEnrichment ( e->first->id ).getID() ;
			Function f = shapefunc[0]* ( s - vm.eval ( s, Point ( 0,1 ) ) ) ;
			f.setIntegrationHint ( hint ) ;
			f.setPointID ( e->first->id ) ;
			f.setDofID ( ids[0] ) ;
			e->setEnrichment ( f  ) ;

			ids = map.getEnrichment ( e->second->id ).getID() ;
			f = shapefunc[1]* ( s - vm.eval ( s, Point ( 0,0 ) ) ) ;
			f.setIntegrationHint ( hint ) ;
			f.setPointID ( e->second->id ) ;
			f.setDofID ( ids[0] ) ;
			e->setEnrichment ( f  ) ;

			ids = map.getEnrichment ( e->third->id ).getID() ;
			f = shapefunc[2]* ( s - vm.eval ( s, Point ( 1,0 ) ) ) ;
			f.setIntegrationHint ( hint ) ;
			f.setPointID ( e->third->id ) ;
			f.setDofID ( ids[0] ) ;
			e->setEnrichment ( f  ) ;

			std::vector<DelaunayTriangle *> toEnrichAlso ;
			for ( size_t j = 0 ; j < e->neighbourhood.size() ; j++ )
			{
				toEnrichAlso.push_back ( e->getNeighbourhood(j) ) ;
			}


			for ( size_t j = 0 ; j < toEnrichAlso.size() ; j++ )
			{
				DelaunayTriangle * elem = toEnrichAlso[j] ;

				Function s_ ( intersectingSegments, elem->getXTransform(), elem->getYTransform() ) ;
				hint.clear() ;

				if ( e->isVertexByID ( elem->first ) &&
				        map.getEnrichment ( elem->first->id ).getType() == SPLIT_ENRICHMENT )
				{
					std::vector<size_t> ids = map.getEnrichment ( elem->first->id ).getID() ;
					Function f = shapefunc[0]* ( s_ - vm.eval ( s_, Point ( 0,1 ) ) ) ;
					f.setIntegrationHint ( hint ) ;
					f.setPointID ( elem->first->id ) ;
					f.setDofID ( ids[0] ) ;
					elem->setEnrichment ( f ) ;

				}

				if ( e->isVertexByID ( elem->second ) &&
				        map.getEnrichment ( elem->second->id ).getType() == SPLIT_ENRICHMENT )
				{
					std::vector<size_t> ids =map.getEnrichment ( elem->second->id ).getID() ;
					Function f = shapefunc[1]* ( s_ - vm.eval ( s_, Point ( 0,0 ) ) ) ;

					f.setIntegrationHint ( hint ) ;
					f.setPointID ( elem->second->id ) ;
					f.setDofID ( ids[0] ) ;
					elem->setEnrichment ( f  ) ;

				}
				if ( e->isVertexByID ( elem->third ) &&
				        map.getEnrichment ( elem->third->id ).getType() == SPLIT_ENRICHMENT )
				{
					std::vector<size_t> ids = map.getEnrichment ( elem->third->id ).getID() ;
					Function f = shapefunc[2]* ( s_ - vm.eval ( s_, Point ( 1,0 ) ) ) ;
					f.setIntegrationHint ( hint ) ;
					f.setPointID ( elem->third->id ) ;
					f.setDofID ( ids[0] ) ;
					elem->setEnrichment ( f  ) ;

				}
			}
		}

		if ( !( e->in ( *getHead() ) || e->in ( *getTail() ) ) )
		{
			Point singularity ;
			double angle = 0
			               ;
			if ( squareDist ( e->getCenter(), *getHead() ) < squareDist ( e->getCenter(), *getTail() ) )
			{
				singularity = *getHead() ;
				angle= headAngle ;
			}
			else
			{
				singularity = *getTail() ;
				angle= tailAngle ;
			}

			std::vector<Point> hint ;
			Line cutting(singularity, Point(cos(angle), sin(angle))) ;

			std::vector<Point> intersec = cutting.intersection ( dynamic_cast<Triangle *> ( e ) ) ;


			if ( intersec.size() == 2 )
			{

				Point ilocal0 = e->inLocalCoordinates ( intersection[0] ) ;
				Point ilocal1 = e->inLocalCoordinates ( intersection[1] ) ;

				hint.push_back ( ilocal0 ) ;
				hint.push_back ( ilocal1 ) ;
			}

			Function x = e->getXTransform() ;
			Function y = e->getYTransform() ;
			double rotatedSingularityX = singularity.x*cos ( angle ) + singularity.y*sin ( angle ) ;
			double rotatedSingularityY = -singularity.x*sin ( angle ) + singularity.y*cos ( angle ) ;
			Function rotatedX = x*cos ( angle ) + y*sin ( angle ) ;
			Function rotatedY = x*sin ( -angle ) + y*cos ( angle ) ;
			Function x_ = x - singularity.x ;
			Function y_ = y - singularity.y ;
			Function theta = f_atan2 ( rotatedY-rotatedSingularityY, rotatedX-rotatedSingularityX );
			Function r = f_sqrt ( ( x_^2 ) + ( y_^2 ) );

			Function f0 = f_sqrt ( r ) *f_sin ( theta/2 );
			Function f1 = f_sqrt ( r ) *f_cos ( theta/2 );
			Function f2 = f_sqrt ( r ) *f_sin ( theta/2 ) *f_cos ( theta );
			Function f3 = f_sqrt ( r ) *f_cos ( theta/2 ) *f_cos ( theta );

			std::vector <Point> pointsLocal ;
			pointsLocal.push_back ( Point ( 0,1 ) ) ;
			pointsLocal.push_back ( Point ( 0,0 ) ) ;
			pointsLocal.push_back ( Point ( 1,0 ) ) ;
			std::vector <Point *> pointsGlobal ;
			pointsGlobal.push_back ( e->first ) ;
			pointsGlobal.push_back ( e->second ) ;
			pointsGlobal.push_back ( e->third ) ;
			for ( size_t j = 0 ; j <3 ; j++ )
			{
				if ( map.getEnrichment ( pointsGlobal[j]->id ).getType() == SINGULAR_ENRICHMENT )
				{
					std::vector<size_t> ids = map.getEnrichment ( pointsGlobal[j]->id ).getID() ;
					Function f = shapefunc[j]* ( f0 - vm.eval ( f0, pointsLocal[j] ) ) ;
					f.setIntegrationHint ( hint ) ;
					f.setPointID ( pointsGlobal[j]->id ) ;
					f.setDofID ( ids[0] ) ;
					e->setEnrichment ( f  ) ;

					f = shapefunc[j]* ( f1 - vm.eval ( f1, pointsLocal[j] ) ) ;
					f.setIntegrationHint ( hint ) ;
					f.setPointID ( pointsGlobal[j]->id ) ;
					f.setDofID ( ids[1] ) ;
					e->setEnrichment (  f  ) ;

					f = shapefunc[j]* ( f2 - vm.eval ( f2, pointsLocal[j] ) ) ;
					f.setIntegrationHint ( hint ) ;
					f.setPointID ( pointsGlobal[j]->id ) ;
					f.setDofID ( ids[2] ) ;
					e->setEnrichment (  f  ) ;

					f = shapefunc[j]* ( f3 - vm.eval ( f3, pointsLocal[j] ) ) ;
					f.setIntegrationHint ( hint ) ;
					f.setPointID ( pointsGlobal[j]->id ) ;
					f.setDofID ( ids[3] ) ;
					e->setEnrichment ( f  ) ;

				}
			}

			std::vector<DelaunayTriangle *> toEnrichAlso ;
			for ( size_t j = 0 ; j < e->neighbourhood.size() ; j++ )
			{
				if ( e->getNeighbourhood(j)->isAlive() && !enrichmentTarget ( e->getNeighbourhood(j) ) )
					toEnrichAlso.push_back ( e->getNeighbourhood(j) ) ;
			}
			for ( size_t j = 0 ; j < toEnrichAlso.size() ; j++ )
			{
				DelaunayTriangle * elem = toEnrichAlso[j] ;

				hint.clear() ;

				x = elem->getXTransform() ;
				y = elem->getYTransform() ;
				rotatedX = x*cos ( angle ) + y*sin ( angle ) ;
				rotatedY = x*sin ( -angle ) + y*cos ( angle ) ;
				x_ = x - singularity.x ;
				y_ = y - singularity.y ;
				theta = f_atan2 ( rotatedY-rotatedSingularityY, rotatedX-rotatedSingularityX );
				r = f_sqrt ( ( x_^2 ) + ( y_^2 ) );

				f0 = f_sqrt ( r ) *f_sin ( theta/2 );
				f1 = f_sqrt ( r ) *f_cos ( theta/2 );
				f2 = f_sqrt ( r ) *f_sin ( theta/2 ) *f_cos ( theta );
				f3 = f_sqrt ( r ) *f_cos ( theta/2 ) *f_cos ( theta );


				pointsGlobal.clear() ;
				pointsGlobal.push_back ( elem->first ) ;
				pointsGlobal.push_back ( elem->second ) ;
				pointsGlobal.push_back ( elem->third ) ;

				for ( size_t j = 0 ; j < 3 ; j++ )
				{
					if ( map.getEnrichment ( pointsGlobal[j]->id ).getType() == SINGULAR_ENRICHMENT )
					{

						std::vector<size_t> ids = map.getEnrichment ( pointsGlobal[j]->id ).getID() ;

						Function f = shapefunc[j]* ( f0 - vm.eval ( f0, pointsLocal[j] ) ) ;
						f.setIntegrationHint ( hint ) ;
						f.setPointID ( pointsGlobal[j]->id ) ;
						f.setDofID ( ids[0] ) ;
						elem->setEnrichment (  f  ) ;

						f = shapefunc[j]* ( f1 - vm.eval ( f1, pointsLocal[j] ) ) ;
						f.setIntegrationHint ( hint ) ;
						f.setPointID ( pointsGlobal[j]->id ) ;
						f.setDofID ( ids[1] ) ;
						elem->setEnrichment (  f  ) ;

						f = shapefunc[j]* ( f2 - vm.eval ( f2, pointsLocal[j] ) ) ;
						f.setIntegrationHint ( hint ) ;
						f.setPointID ( pointsGlobal[j]->id ) ;
						f.setDofID ( ids[2] ) ;
						elem->setEnrichment (  f  ) ;

						f = shapefunc[j]* ( f3 - vm.eval ( f3, pointsLocal[j] ) ) ;
						f.setIntegrationHint ( hint ) ;
						f.setPointID ( pointsGlobal[j]->id ) ;
						f.setDofID ( ids[3] ) ;
						elem->setEnrichment (  f  ) ;
					}
				}
			}
		}
	}

}

bool Crack::interacts ( Feature * f ) const
{
	for ( PointSet::const_iterator i = this->begin() ; i < this->end() ; i++ )
		if ( f->inBoundary ( ( *i ) ) )
			return true ;
	return false ;
}

Point * Crack::pointAfter ( size_t i )
{
	Point * to_insert = new Point ( *boundingPoints[i]*0.5 + *boundingPoints[ ( i+1 ) %boundingPoints.size() ]*0.5 ) ;
	std::valarray<Point *> temp ( boundingPoints.size() +1 ) ;
	std::copy ( &boundingPoints[0], &boundingPoints[i], &temp[0] ) ;
	temp[i+1] = to_insert ;
	std::copy ( &boundingPoints[i+1], &boundingPoints[boundingPoints.size() ], &temp[i+2] ) ;
	boundingPoints.resize ( temp.size() ) ;
	std::copy ( &temp[0],&temp[temp.size() ] , &boundingPoints[0] ) ;
	return to_insert ;
}

bool Crack::inBoundary ( const Point & v ) const
{
	return false ;
	Point p ( v.x, v.y );

	this->project ( &p ) ;

	return squareDist2D(v, p) < infRad*infRad ;
}

std::vector<Point *> Crack::getSamplingPoints() const
{
	std::vector<Point *> ret ;

	Point v0 = getBoundingPoint ( 0 ) - getBoundingPoint ( 1 ) ;
	v0 = v0/v0.norm() ;
	Point n0_0 ( v0.y, -v0.x ) ;
	Point n0_1 ( -v0.y, v0.x ) ;
	Point v1 = getBoundingPoint ( boundingPoints.size()-1 ) - getBoundingPoint ( boundingPoints.size()-2 ) ;
	Point n1_0 ( v1.y, -v1.x ) ;
	Point n1_1 ( -v1.y, v1.x ) ;
	v1 = v1/v1.norm() ;

	ret.push_back ( new Point ( getBoundingPoint ( 0 ) + v0*0.9*infRad ) ) ;
	ret.push_back ( new Point ( getBoundingPoint ( 0 ) - v0*0.45*infRad + n0_0*0.45*infRad ) ) ;
	ret.push_back ( new Point ( getBoundingPoint ( 0 ) - v0*0.45*infRad + n0_1*0.45*infRad ) ) ;

	ret.push_back ( new Point ( getBoundingPoint ( boundingPoints.size()-1 ) + v1*0.9*infRad ) ) ;
	ret.push_back ( new Point ( getBoundingPoint ( boundingPoints.size()-1 ) - v1*0.45*infRad + n1_0*0.45*infRad ) ) ;
	ret.push_back ( new Point ( getBoundingPoint ( boundingPoints.size()-1 ) - v1*0.45*infRad + n1_1*0.45*infRad ) ) ;

	return ret ;

}

bool Crack::inBoundary ( const Point *v ) const
{
	return false ;
	Point p ( v->x, v->y );

	this->project ( &p ) ;

	return squareDist2D(v, &p) < infRad*infRad ;

}

Crack::~Crack()
{
}

Crack::EnrichmentData::EnrichmentData ( std::vector<size_t> id, EnrichmentType t ) : iD ( id ), type ( t ) { };

Crack::EnrichmentData::EnrichmentData ( size_t id, EnrichmentType t ) : type ( t )
{
	iD.push_back ( id ) ;
}

Crack::EnrichmentData::EnrichmentData() :  iD ( 0 ), type ( VOID_ENRICHMENT )
{
}

std::vector<size_t> Crack::EnrichmentData::getID() const
{
	return iD ;
}

void Crack::EnrichmentData::setID ( const std::vector<size_t> i )
{
	iD =i ;
}

Crack::EnrichmentType Crack::EnrichmentData::getType() const
{
	return type ;
}

void Crack::EnrichmentData::setType ( EnrichmentType t )
{
	type = t ;
}

std::pair<double, double> Crack::computeJIntegralAtHead ( double dt, const DelaunayTree * dtree )
{
	Point direction ( *getHead()- getBoundingPoint ( 1 ) ) ;
	Segment tip ( *getHead() , getBoundingPoint ( 1 ) ) ;
	Circle c ( boundary->getRadius()*1.5, boundary->getCenter() ) ;

	std::vector<DelaunayTriangle *> disk = dtree->conflicts ( &c ) ;
	std::vector<DelaunayTriangle *> ring ;

	if(disk.size() == 1)
	{
		setInfluenceRadius(disk[0]->getRadius()) ;
		c.setRadius(disk[0]->getRadius()) ;
		disk = dtree->conflicts ( &c ) ;
	}
	
	std::vector<std::pair<Segment *, DelaunayTriangle *> > gamma ;


	for ( size_t j = 0 ; j < disk.size() ; j++ )
	{
		Triangle * t = dynamic_cast<Triangle *> ( disk[j] ) ;
		if ( t->intersects ( &c ) )
		{
			ring.push_back ( disk[j] ) ;

			Segment A ( *disk[j]->first, *disk[j]->second ) ;
			Segment B ( *disk[j]->first, *disk[j]->third ) ;
			Segment C ( *disk[j]->second, *disk[j]->third ) ;

			if ( tip.intersects ( t ) )
			{
				if ( A.intersects ( &c ) && B.intersects ( &c ) && !c.in ( C.midPoint() ) )
				{
					Point ipoint = C.intersection ( tip ) ;
					gamma.push_back ( std::make_pair ( new Segment ( *disk[j]->third, ipoint ),disk[j] ) ) ;
					gamma.push_back ( std::make_pair ( new Segment ( *disk[j]->second, ipoint ),disk[j] ) ) ;
				}
				if ( B.intersects ( &c ) && C.intersects ( &c ) && !c.in ( A.midPoint() ) )
				{
					Point ipoint = A.intersection ( tip ) ;
					gamma.push_back ( std::make_pair ( new Segment ( *disk[j]->first, ipoint ),disk[j] ) ) ;
					gamma.push_back ( std::make_pair ( new Segment ( *disk[j]->second, ipoint ),disk[j] ) ) ;
				}
				if ( C.intersects ( &c ) && A.intersects ( &c ) && !c.in ( B.midPoint() ) )
				{
					Point ipoint = B.intersection ( tip ) ;
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

			Point localVector = gamma[j].second->inLocalCoordinates ( tip.vector() ) ;
			Point localNormal ( -localVector.y, localVector.x ) ;

			for ( size_t k = 0 ; k < gaussPoints.size() ; k++ )
			{
				Point localGaussPoint = gamma[j].second->inLocalCoordinates ( gaussPoints[k].first ) ;
				Matrix sigma = gamma[j].second->getState().getStressMatrix ( localGaussPoint,true ) ;
				Matrix epsilon = gamma[j].second->getState().getStrainMatrix ( localGaussPoint,true ) ;
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

std::pair<double, double> Crack::computeJIntegralAtTail ( double dt, const DelaunayTree * dtree )
{
	Point direction ( *getTail()- getBoundingPoint ( getBoundingPoints().size()-2 ) ) ;
	Segment tip ( *getTail() , getBoundingPoint ( getBoundingPoints().size()-2 ) ) ;
	Circle c ( boundary2->getRadius()*1.5, boundary2->getCenter() ) ;


	std::vector<DelaunayTriangle *> disk = dtree->conflicts ( &c ) ;
	std::vector<DelaunayTriangle *> ring ;
	if(disk.size() == 1)
	{
		setInfluenceRadius(disk[0]->getRadius()) ;
		c.setRadius(disk[0]->getRadius()) ;
		disk = dtree->conflicts ( &c ) ;
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

			if ( tip.intersects ( t ) )
			{

				if ( A.intersects ( &c ) && B.intersects ( &c ) && !c.in ( C.midPoint() ) )
				{
					Point ipoint = C.intersection ( tip ) ;
					gamma.push_back ( std::make_pair ( new Segment ( *disk[j]->third, ipoint ),disk[j] ) ) ;
					gamma.push_back ( std::make_pair ( new Segment ( *disk[j]->second, ipoint ),disk[j] ) ) ;
				}
				if ( B.intersects ( &c ) && C.intersects ( &c ) && !c.in ( A.midPoint() ) )
				{
					Point ipoint = A.intersection ( tip ) ;
					gamma.push_back ( std::make_pair ( new Segment ( *disk[j]->first, ipoint ),disk[j] ) ) ;
					gamma.push_back ( std::make_pair ( new Segment ( *disk[j]->second, ipoint ),disk[j] ) ) ;
				}
				if ( C.intersects ( &c ) && A.intersects ( &c ) && !c.in ( B.midPoint() ) )
				{
					Point ipoint = B.intersection ( tip ) ;
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

			Point localVector = gamma[j].second->inLocalCoordinates ( tip.vector() ) ;
			Point localNormal ( -localVector.y, localVector.x ) ;

			for ( size_t k = 0 ; k < gaussPoints.size() ; k++ )
			{
				Point localGaussPoint = gamma[j].second->inLocalCoordinates ( gaussPoints[k].first ) ;
				Matrix sigma = gamma[j].second->getState().getStressMatrix ( localGaussPoint,true ) ;
				Matrix epsilon = gamma[j].second->getState().getStrainMatrix ( localGaussPoint,true ) ;
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


void Crack::setCriticalJ(double newJ)
{
	criticalJ = newJ ;
}

double Crack::getCriticalJ() const
{
	return criticalJ ;
}

void Crack::step( double dt, std::valarray<double> *, const DelaunayTree * dtree )
{

	changed = false ;

// 	return ;

	double norm = .00025 ;

	std::pair<double, double> headJ = computeJIntegralAtHead ( dt, dtree ) ;
	Vector J ( 2 ) ; J[0] = headJ.first ; J[1] = headJ.second ;
	std::cout << "at head : " << J[0] << ", " << J[1] << std::endl ;
	std::cout << "J angle is " << atan2 ( J[1], J[0] ) << std::endl ;

// 	Circle atHead ( infRad, this->boundary->getCenter() ) ;
	std::vector<DelaunayTriangle *> disk = dtree->conflicts ( this->boundary ) ;
	Circle c0(*static_cast<Circle *>(boundary)) ;
	disk = dtree->conflicts (&c0) ;
	while( disk.size() < 16)
	{
		c0.setRadius(c0.getRadius()*2.) ;
		disk = dtree->conflicts ( &c0 ) ;
	}
	std::vector<DelaunayTriangle *> tris ;
	Point direction ;
	double count = 0 ;
	double acount = 0 ;
	double aangle = 0 ;
	double currentAngle = atan2 ( getHead()->y - boundingPoints[1]->y,
					getHead()->x -boundingPoints[1]->x ) ;
	Point lastDir ( getHead()->x-getBoundingPoint ( 1 ).x, getHead()->y-getBoundingPoint ( 1 ).y ) ;


	criticalJ = 0.;

	if(sqrt(J[0]*J[0] + J[1]*J[1]) >= criticalJ)
	{
		DelaunayTriangle * headElem = NULL ;
		for ( size_t i = 0 ; i < disk.size() ; i++ )
		{
			if ( !disk[i]->in ( *getHead() )
				&& (!disk[i]->intersects ( dynamic_cast<SegmentedLine *> ( this ) )
				    || boundary->in(disk[i]->getCenter()))
			)
				tris.push_back ( disk[i] ) ;
			else if ( disk[i]->in ( *getHead() ) )
				headElem = disk[i] ;
		}
		
		if ( !tris.empty() )
		{
	
			for ( size_t i = 0 ; i < tris.size() ; i++ )
			{
	
				DelaunayTriangle * current = tris[i] ;
				for ( size_t j = 0 ; j < current->getBoundingPoints().size() ; j++ )
				{
	
					Point currentDir ( current->getBoundingPoint ( j ).x-getHead()->x, current->getBoundingPoint ( j ).y-getHead()->y ) ;
					if ( ( currentDir*lastDir ) > 0 )
					{
						Vector principalStresses = current->getState().getPrincipalStresses( 
							current->getBoundingPoint ( j ) ) ;
						double maxPrincipalStressCurrent = 1./std::abs(
							principalStresses[0] );
						if ( current->getBehaviour()->type == VOID_BEHAVIOUR )
							maxPrincipalStressCurrent = 1000 ;
	
						if ( currentDir.norm() > 1e-8 )
						{
							aangle += current->getState().getPrincipalAngle (
								current->getBoundingPoint ( j ) ) ;
							acount++ ;
							direction += ( currentDir/currentDir.norm() )
								*maxPrincipalStressCurrent ;
							count += maxPrincipalStressCurrent ;
						}
					}
				}
			}
		}
		if ( count> 1e-8 )
			direction/=direction.norm() ;
		else
		{
			direction = Point ( getHead()->x-boundingPoints[1]->x,
					getHead()->y-boundingPoints[1]->y ) ;
			direction /= direction.norm() ;
	
		}
	
		direction*= stepLength ;
		if ( acount )
			aangle /= acount ;
	
		if ( headElem )
		{
			changed = true ;

			if ( headElem->getBehaviour()->type == VOID_BEHAVIOUR )
			{
				setInfluenceRadius ( 0 ) ;
			}
			else
			{
				direction.set ( stepLength*cos ( aangle + M_PI*.25), 
				                stepLength*sin ( aangle + M_PI*.25) ) ;
	
				if ( ( direction*lastDir ) < 0 )
					direction.set ( stepLength*cos ( aangle+M_PI + M_PI*.25), 
					                stepLength*sin ( aangle+M_PI + M_PI*.25) ) ;
	
				std::valarray<Point *> newPoints ( boundingPoints.size() +1 ) ;
				newPoints[0] = new Point ( getHead()->x+direction.x, getHead()->y+direction.y ) ;
				for ( size_t i = 1 ; i < newPoints.size() ; i++ )
				{
					newPoints[i] = boundingPoints[i-1] ;
				}
				setBoundingPoints ( newPoints ) ;
				setInfluenceRadius ( infRad ) ;
			}
		}
	
	}
	
	std::pair<double, double> tailJ = computeJIntegralAtTail ( dt, dtree ) ;
	J[0] = tailJ.first ; J[1] = tailJ.second ;
	std::cout << "at tail : " << J[0] << ", " << J[1] << std::endl ;
	std::cout << "J angle is " << atan2 ( J[1], J[0] ) << std::endl ;


// 	Circle atTail ( infRad, this->boundary2->getCenter() ) ;
	Circle c(*static_cast<Circle *>(boundary2)) ;
	disk = dtree->conflicts (&c) ;
	while( disk.size() < 16)
	{
		c.setRadius(c.getRadius()*2.) ;
		disk = dtree->conflicts ( &c ) ;
	}
	tris.clear() ;
	if(sqrt(J[0]*J[0] + J[1]*J[1]) >= criticalJ)
	{
		DelaunayTriangle * tailElem = NULL;
		for ( size_t i = 0 ; i < disk.size() ; i++ )
		{
			
			if ( !disk[i]->in ( *getTail() )
				&& (!disk[i]->intersects ( getPrimitive() )
				    || boundary2->in(disk[i]->getCenter()))
			)
			{
				tris.push_back ( disk[i] ) ;
			}
			else if ( disk[i]->in ( *getTail() ) )
			{
				tailElem = disk[i] ;
			}
		}
	
		count = 0 ;
		acount = 0 ;
		aangle = 0 ;
		direction = Point() ;
	
		currentAngle = atan2 ( getTail()->y-boundingPoints[boundingPoints.size()-2]->y,
				getTail()->x-boundingPoints[boundingPoints.size()-2]->x ) ;
	
		lastDir = Point ( getTail()->x-getBoundingPoint ( getBoundingPoints().size()-2 ).x, getTail()->y-getBoundingPoint ( getBoundingPoints().size()-2 ).y ) ;
		if ( !tris.empty() )
		{
	
			for ( size_t i = 0 ; i < tris.size() ; i++ )
			{
	
				DelaunayTriangle * current = tris[i] ;
	
				for ( size_t j = 0 ; j < current->getBoundingPoints().size() ; j++ )
				{
	
					Point currentDir ( current->getBoundingPoint ( j ).x-getTail()->x,
					                   current->getBoundingPoint ( j ).y-getTail()->y ) ;
	
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
							direction += ( currentDir/currentDir.norm() )
								*maxPrincipalStressCurrent ;
							count += maxPrincipalStressCurrent ;
						}
					}
				}
			}
		}
	
		if ( count> 1e-8 )
			direction/=direction.norm() ;
		else
		{
			direction = Point ( getTail()->x-boundingPoints[boundingPoints.size()-2]->x,
					getTail()->y-boundingPoints[boundingPoints.size()-2]->y ) ;
			direction /= direction.norm() ;
	
		}
		direction*= stepLength ;
		if ( acount )
			aangle /= acount ;
	
		if ( tailElem )
		{
			changed = true ;

			if ( tailElem->getBehaviour()->type == VOID_BEHAVIOUR )
			{
				setInfluenceRadius ( 0 ) ;
			}
			else
			{
	
				direction.set ( stepLength*cos ( aangle+ M_PI*.25), 
				                stepLength*sin ( aangle + M_PI*.25) ) ;
	
				if ( ( direction*lastDir ) < 0 )
					direction.set ( stepLength*cos ( aangle+M_PI + M_PI*.25), 
					                stepLength*sin ( aangle+M_PI + M_PI*.25) ) ;

	
				std::valarray<Point *> newPoints ( boundingPoints.size() +1 ) ;
				for ( size_t i = 0 ; i < boundingPoints.size() ; i++ )
				{
					newPoints[i] = boundingPoints[i] ;
				}
				newPoints[boundingPoints.size() ] = new Point ( getTail()->x+direction.x,
					getTail()->y+direction.y ) ;
				setBoundingPoints ( newPoints ) ;
				setInfluenceRadius ( infRad ) ;
			}
		}
	}
	
	print() ;
	
}

void Crack::print() const
{
	for(size_t i= 0 ; i < getBoundingPoints().size() ;i++)
	{
		getBoundingPoint(i).print();
	}
}

bool Crack::EnrichmentData::enriched ( DelaunayTriangle * t )
{
	std::pair<std::pair<int, int>, int > t_ ( std::make_pair ( t->first->id, t->second->id ),  t->third->id ) ;
	if ( state.find ( t_ ) == state.end() )
	{
		state[t_] = false ;
		return false ;
	}

	return state[t_] ;
}

void Crack::EnrichmentData::setEnriched ( DelaunayTriangle * t, bool b )
{
	std::pair<std::pair<int, int>, int > t_ ( std::make_pair ( t->first->id, t->second->id ),  t->third->id ) ;
	state[t_] = b ;
}

bool Crack::EnrichmentMap::inMap ( size_t p ) const
{
	return props.find ( p ) != props.end() ;
}

Crack::EnrichmentMap::EnrichmentMap()
{
}

void Crack::EnrichmentMap::EnrichmentMap::update ( std::vector<DelaunayTriangle *> * my_triangles, Crack * myself, size_t &start )
{
	props.clear() ;
	for ( size_t i = 0 ; i < my_triangles->size() ; i++ )
	{
		props[ ( *my_triangles ) [i]->first->id] = EnrichmentData ( 0, VOID_ENRICHMENT ) ;
		props[ ( *my_triangles ) [i]->second->id] = EnrichmentData ( 0, VOID_ENRICHMENT ) ;
		props[ ( *my_triangles ) [i]->third->id] = EnrichmentData ( 0, VOID_ENRICHMENT ) ;
	}

	for ( size_t i = 0 ; i < my_triangles->size() ; i++ )
	{

		if ( ( *my_triangles ) [i]->getBehaviour()->type != VOID_BEHAVIOUR )
		{
			bool intersection = ( myself->boundary->in ( * ( *my_triangles ) [i]->first )   ||
			                      myself->boundary->in ( * ( *my_triangles ) [i]->second )  ||
			                      myself->boundary->in ( * ( *my_triangles ) [i]->third ) )  ||
			                    ( myself->boundary2->in ( * ( *my_triangles ) [i]->first ) ||
			                      myself->boundary2->in ( * ( *my_triangles ) [i]->second ) ||
			                      myself->boundary2->in ( * ( *my_triangles ) [i]->third ) ) ||
			                    ( *my_triangles ) [i]->in ( *myself->getHead() ) ||
			                    ( *my_triangles ) [i]->in ( *myself->getTail() );


			if ( intersection )
			{

				if ( props[ ( *my_triangles ) [i]->first->id].getType() == VOID_ENRICHMENT )
				{
					std::vector<size_t> id0 ;
					id0.push_back ( start++ ) ;
					id0.push_back ( start++ ) ;
					id0.push_back ( start++ ) ;
					id0.push_back ( start++ ) ;
					props[ ( *my_triangles ) [i]->first->id] = EnrichmentData ( id0, SINGULAR_ENRICHMENT ) ;
				}
				if ( props[ ( *my_triangles ) [i]->second->id].getType() == VOID_ENRICHMENT )
				{
					std::vector<size_t> id0 ;
					id0.push_back ( start++ ) ;
					id0.push_back ( start++ ) ;
					id0.push_back ( start++ ) ;
					id0.push_back ( start++ ) ;
					props[ ( *my_triangles ) [i]->second->id] = EnrichmentData ( id0, SINGULAR_ENRICHMENT ) ;
				}
				if ( props[ ( *my_triangles ) [i]->third->id].getType() == VOID_ENRICHMENT )
				{
					std::vector<size_t> id0 ;
					id0.push_back ( start++ ) ;
					id0.push_back ( start++ ) ;
					id0.push_back ( start++ ) ;
					id0.push_back ( start++ ) ;
					props[ ( *my_triangles ) [i]->third->id] = EnrichmentData ( id0, SINGULAR_ENRICHMENT ) ;
				}
			}
		}
	}

	for ( size_t i = 0 ; i < my_triangles->size() ; i++ )
	{
		if ( ( *my_triangles ) [i]->getBehaviour()->type != VOID_BEHAVIOUR )
		{
			if ( props[ ( *my_triangles ) [i]->first->id].getType() == VOID_ENRICHMENT )
			{
				props[ ( *my_triangles ) [i]->first->id] = EnrichmentData ( start++, SPLIT_ENRICHMENT ) ;
			}

			if ( props[ ( *my_triangles ) [i]->second->id].getType() == VOID_ENRICHMENT )
			{
				props[ ( *my_triangles ) [i]->second->id] = EnrichmentData ( start++, SPLIT_ENRICHMENT ) ;
			}

			if ( props[ ( *my_triangles ) [i]->third->id].getType() == VOID_ENRICHMENT )
			{
				props[ ( *my_triangles ) [i]->third->id] = EnrichmentData ( start++, SPLIT_ENRICHMENT ) ;
			}
		}
	}

}

Crack::EnrichmentData Crack::EnrichmentMap::getEnrichment ( size_t p )
{
	if ( props.find ( p ) == props.end() )
		return EnrichmentData ( 0, VOID_ENRICHMENT ) ;
	return props[p] ;
}


void Crack::snap ( DelaunayTree * dtree )
{
	std::vector<DelaunayTreeItem *> targets = dtree->conflicts ( this->getHead() ) ;

	std::map<double, Point *> positions ;

	for ( size_t i = 0 ; i< targets.size() ; i++ )
	{
		if ( targets[i]->isTriangle )
		{
			DelaunayTriangle * t =  dynamic_cast<DelaunayTriangle *> ( targets[i] ) ;
			positions[dist ( *getHead(), t->getCenter() ) ] = &t->getCenter() ;
		}
	}
	if ( !positions.empty() )
		this->getHead()->set ( positions.begin()->second->x, positions.begin()->second->y ) ;

	targets = dtree->conflicts ( this->getTail() ) ;

	positions.clear() ;

	for ( size_t i = 0 ; i< targets.size() ; i++ )
	{
		if ( targets[i]->isTriangle )
		{
			DelaunayTriangle * t =  dynamic_cast<DelaunayTriangle *> ( targets[i] ) ;
			positions[dist ( *getHead(), t->getCenter() ) ] = &t->getCenter() ;
		}
	}
	if ( !positions.empty() )
		this->getTail()->set ( positions.begin()->second->x, positions.begin()->second->y ) ;

}

bool Crack::moved() const
{
	return changed ;
}

/**
  Prints the geometry of the crack to a file
**/
void Crack::printFile(const std::string& filename) const // SB
{
	std::fstream filestr;
	filestr.open (filename.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
  
  for(size_t i = 0 ; i < this->getBoundingPoints().size() ; i++)
    {
      filestr << this->getBoundingPoint(i).x << "  " << this->getBoundingPoint(i).y << std::endl ;
    }
  filestr.close();
}
