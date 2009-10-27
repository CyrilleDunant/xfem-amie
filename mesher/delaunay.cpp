//
// C++ Implementation: delaunay
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//


#include "delaunay.h"
#include <limits>

#define DEBUG
#undef DEBUG

using namespace Mu ;

DelaunayTreeItem::DelaunayTreeItem( DelaunayTree * t, DelaunayTreeItem * father,  const Point * c) : stepson(0), neighbour(0), son(0)
{
	tree = t ;
	this->stepfather = NULL ;
	this->father = father;
	this->m_c  = c ;
	this->dead = false ;
	this->visited = false ;
	this->erased = false ;
	
	this->isPlane =false ;
	this->isTriangle =false ;
	this->isDeadTriangle =false ;
	this->first = NULL ;
	this->second = NULL ;
	this->third = NULL ;

	index = 0 ;
	if(tree)
	{
		index = tree->tree.size() ;
		tree->tree.push_back(this) ;
	}
}
	
bool DelaunayTriangle::isConflicting(const Geometry * g) const
{
// 	Circle c(getRadius(), getCircumCenter()) ;
	return g->in(*first) 
		|| g->in(*second) 
		|| g->in(*third) 
		|| in(g->getCenter()) 
		|| g->intersects(getPrimitive()) ;
	
}

bool DelaunayDeadTriangle::isConflicting(const Geometry * g) const
{
	Triangle t(*first, *second, *third) ;
	return g->in(*first) 
		|| g->in(*second) 
		|| g->in(*third) 
		|| inCircumCircle(g->getCenter()) 
		|| g->intersects(&t) ;
	
}

bool DelaunayDemiPlane::isConflicting(const Geometry * g) const
{
	if(g->in(*first))
		return true ;
	if(g->in(*second))
		return true ;
	if(inCircumCircle(g->getCenter()))
		return true ;
	
	return false ;
}

bool DelaunayTreeItem::isConflicting(const Geometry * g) const
{
	return inCircumCircle(g->getCenter()) || g->in(*first) || g->in(*second)  ;
}



DelaunayTreeItem::~DelaunayTreeItem()
{
}
	
const Point * DelaunayTreeItem::killer() const 
{
	return m_k ;
}

const Point * DelaunayTreeItem::creator() const 
{
	return m_c ;
}

void  DelaunayTreeItem::setCreator(const Point * p)  
{
	m_c = p ;
}

size_t DelaunayTreeItem::numberOfCommonVertices(const DelaunayTreeItem * s) const
{
	if(this->isTriangle && s->isTriangle)
	{
		size_t ret = 0 ;
		
		if(s->first == this->first || s->first == this->second || s->first == this->third)
			ret++ ;
		if(s->second == this->first || s->second == this->second || s->second == this->third)
			ret++ ;
		if(s->third == this->first || s->third == this->second || s->third == this->third)
			ret++ ;
	
		assert(ret < 4) ;
			
		return ret ;
	}
	if(this->isTriangle && s->isPlane)
	{
		size_t ret = 0 ;
		
		if(isAligned(s->first, s->second, first))
			ret++ ;
		if(isAligned(s->first, s->second, second))
			ret++ ;
		if(isAligned(s->first, s->second, third))
			ret++ ;
		
		assert(ret < 4) ;
		
		return ret ;
	}
	if(this->isPlane && s->isPlane)
	{
		size_t ret = 0 ;
		
		if(isAligned(s->first, s->second, first))
			ret++ ;
		if(isAligned(s->first, s->second, second))
			ret++ ;
		
		assert(ret < 4) ;
		
		return ret ;
	}
	
	if(this->isPlane && s->isTriangle)
	{
		size_t ret = 0 ;
		
		if(isAligned(this->first, this->second, s->first))
			ret++ ;
		if(isAligned(this->first, this->second, s->second))
			ret++ ;
		if(isAligned(this->first, this->second, s->third))
			ret++ ;
		
		assert(ret < 4) ;
		
		return ret ;
	}
	
	return 0 ;
}

void Star::updateNeighbourhood()
{
	std::vector<DelaunayTreeItem *> items ;

	for(std::vector<DelaunayTreeItem *>::iterator i =  treeitem.begin() ; i != treeitem.end() ;++i)
	{	
		for(size_t j = 0 ; j < (*i)->son.size() ; j++)
			items.push_back((*i)->getSon(j)) ;
		for(size_t j = 0 ; j < (*i)->stepson.size() ; j++)
			items.push_back((*i)->getStepson(j)) ;
		for(size_t j = 0 ; j < (*i)->neighbour.size() ; j++)
			items.push_back((*i)->getNeighbour(j)) ;
// 		items.insert(items.end() , &treeitem[i]->son[0] , &treeitem[i]->son[treeitem[i]->son.size()]) ;
// 		items.insert(items.end() , treeitem[i]->neighbour.begin() , treeitem[i]->neighbour.end()) ;
	}
	std::sort(items.begin(), items.end()) ;
	std::vector<DelaunayTreeItem *>::iterator e = std::unique(items.begin(), items.end()) ;
	items.erase(e, items.end()) ;
	
	bool goOn = true ;
	while(goOn)
	{
		goOn = false ;
		for(std::vector<DelaunayTreeItem *>::iterator i = items.begin() ; i != items.end() ; ++i)
		{
			for(std::vector<DelaunayTreeItem *>::iterator j = i+1 ; j != items.end() ; ++j)
			{
				if((*i)->isDuplicate((*j)) && (*i)->isAlive() && (*j)->isAlive())
				{
					for(size_t k = 0 ; k < (*j)->neighbour.size() ; k++)
					{
						makeNeighbours((*j)->getNeighbour(k), (*i)) ;
					}
					
					(*j)->kill((*j)->creator()) ;
					(*j)->erased = true ;
					std::vector<unsigned int> newSons;
					newSons.insert(newSons.end(), 
						&(*j)->father->son[0] , 
						&(*j)->father->son[(*j)->father->son.size()]) ;
					newSons.erase(std::find(
								newSons.begin(), 
								newSons.end(), 
								(*j)->index
								)
							) ;
					(*j)->father->son.resize((*j)->father->son.size()-1);
					std::copy(newSons.begin(), newSons.end(), &(*j)->father->son[0] ) ;

					std::vector<unsigned int> newStepsons;
					newStepsons.insert(newStepsons.end(), 
							&(*j)->stepfather->stepson[0],
							&(*j)->stepfather->stepson[(*j)->stepfather->stepson.size()]) ;
					newStepsons.erase(std::find(
								newStepsons.begin(), 
								newStepsons.end(), 
								(*j)->index
								)
							) ;
					(*j)->stepfather->stepson.resize((*j)->stepfather->stepson.size()-1);
					std::copy( newStepsons.begin(), newStepsons.end(),&(*j)->stepfather->stepson[0]) ;
					goOn = true ;
					items.erase(j) ;
					break ;
				}
			}
			
			if(goOn)
				break ;
		}
	}

	for(std::vector<DelaunayTreeItem *>::iterator i = items.begin() ; i != items.end() ; ++i)
	{
		for(std::vector<DelaunayTreeItem *>::iterator j = i+1 ; j != items.end() ; ++j)
		{
			if(!(*i)->erased && !(*j)->erased )
				makeNeighbours((*i), (*j)) ;
		}
	}
}

bool DelaunayTreeItem::isDuplicate(const DelaunayTreeItem * t) const
{
	return t!=this && this->isTriangle && t->isTriangle && this->numberOfCommonVertices(t) == 3 ;
}

void DelaunayTree::addSharedNodes(size_t nodes_per_side, size_t time_planes, double timestep, const TriElement * father)
{
	std::vector<DelaunayTriangle *> tri = getTriangles() ;

	for(std::vector<DelaunayTriangle *>::iterator i = tri.begin() ; i != tri.end() ; ++i)
	{
		
		(*i)->visited = true ;
			
		size_t nodes_per_plane = nodes_per_side*3+3 ;
		
		std::valarray<Point *> newPoints(nodes_per_plane*time_planes) ;
		std::valarray<bool> done(false, nodes_per_plane*time_planes) ;
		
		for(size_t plane = 0 ; plane < time_planes ; plane++)
		{
			for(size_t side = 0 ; side < 3 ; side++)
			{
				Point a((*i)->getBoundingPoint(side)) ;
				Point b((*i)->getBoundingPoint((side+1)%3)) ;
				
				if(time_planes> 1)
				{
					a.t = (double)plane*(timestep/(double)(time_planes-1))-timestep/2.;
					b.t = (double)plane*(timestep/(double)(time_planes-1))-timestep/2.;
				}
				for(size_t node = 0 ; node < nodes_per_side+1 ; node++)
				{
					double fraction = (double)(node)/((double)nodes_per_side+1) ;
					Point proto = a*(1.-fraction) + b*fraction ;
					Point * foundPoint = NULL ;
					
					for(size_t j = 0 ; j< (*i)->getBoundingPoints().size() ; j++)
					{
						if((*i)->getBoundingPoint(j) == proto)
						{
							foundPoint = &(*i)->getBoundingPoint(j) ;
							break ;
						}
					}
					
					if(!foundPoint)
					{
						for(size_t j = 0 ; j < (*i)->neighbourhood.size() ; j++)
						{
							if((*i)->getNeighbourhood(j)->visited)
							{
								DelaunayTriangle * n = (*i)->getNeighbourhood(j) ;
								for(size_t k = 0 ; k < n->getBoundingPoints().size();k++)
								{
									if(n->getBoundingPoint(k) == proto)
									{
										foundPoint = &n->getBoundingPoint(k) ;
										break ;
									}
								}
								
								if(foundPoint)
									break ;
							}
						}
					}
					
					if(!done[nodes_per_plane*plane+side*(nodes_per_side+1)+node])
					{
						if(foundPoint)
						{
							newPoints[nodes_per_plane*plane+side*(nodes_per_side+1)+node]  = foundPoint ;
						}
						else
						{
							additionalPoints.push_back(new Point(proto) ) ;
							newPoints[nodes_per_plane*plane+side*(nodes_per_side+1)+node]  = additionalPoints.back();
							newPoints[nodes_per_plane*plane+side*(nodes_per_side+1)+node]->id = global_counter++ ;
						}
						
						done[nodes_per_plane*plane+side*(nodes_per_side+1)+node] = true ;
					}
				}
			}
		}
		
		(*i)->setBoundingPoints(newPoints) ;
	}
			
	
	for(std::vector<DelaunayTriangle *>::iterator i = tri.begin() ; i != tri.end() ; ++i)
	{
		(*i)->clearVisited() ;
	}

}

void DelaunayTree::refresh(TriElement *father)
{
	
	for(std::vector<DelaunayTreeItem *>::iterator i = tree.begin() ; i != tree.end() ; ++i)
	{
		if((*i)->isAlive() && (*i)->isTriangle)
		{
			(static_cast<DelaunayTriangle *>(*i))->refresh(father) ;
		}
	}
}

size_t DelaunayTree::numPoints() const
{
	return this->global_counter ;
}

void DelaunayTreeItem::conflicts(std::pair<std::vector<DelaunayTriangle *>, std::vector<DelaunayTreeItem *> > & ret, const Geometry *g)
{
	if(visited)
	{
		return ;
	}
	visited = true ;
	
	ret.second.push_back(this) ;
	
	
// 	if(!isConflicting(g) )
// 	{
// 		return ;
// 	}
	
	for (size_t i  = 0 ;  i < stepson.size() ; i++)
	{
		bool limit = false ;
// 		if(!getStepson(i)->visited && getStepson(i)->isTriangle && !getStepson(i)->isDeadTriangle)
// 		{
// 			DelaunayTriangle * t = static_cast<DelaunayTriangle *>(getStepson(i)) ;
// 			limit = std::abs(
// 			                  squareDist2D(
// 			                                t->getCircumCenter(),g->getCenter())
// 			                  -(t->getRadius()+g->getRadius())*(t->getRadius()+g->getRadius())
// 			                ) 
// 				< .1*t->getRadius()*t->getRadius() ;
// 		}
// 		if(!getStepson(i)->visited && !getStepson(i)->isDeadTriangle)
// 		{
// 			DelaunayDeadTriangle * t = static_cast<DelaunayDeadTriangle *>(getStepson(i)) ;
// 			limit = std::abs(
// 			                  squareDist2D(
// 			                                t->getCircumCenter(),&g->getCenter())
// 			                  -(t->getRadius()+g->getRadius())*(t->getRadius()+g->getRadius()))
// 				< .1*t->getRadius()*t->getRadius() ;
// 		}
		
		if( (!getStepson(i)->visited && getStepson(i)->isConflicting(g)) || limit)
		{
			getStepson(i)->conflicts(ret, g) ;
		}
	}
	
	for (size_t i  = 0 ;  i < son.size() ; i++)
	{
		bool limit = false ;

// 		if(!getSon(i)->visited && getSon(i)->isTriangle && !getSon(i)->isDeadTriangle)
// 		{
// 			DelaunayTriangle * t = static_cast<DelaunayTriangle *>(getSon(i)) ;
// 			limit = std::abs(
// 			                  squareDist2D(t->getCircumCenter(),g->getCenter())
// 			                  -(t->getRadius()+g->getRadius())*(t->getRadius()+g->getRadius())) 
// 				< .01*t->getRadius()*t->getRadius() ;
// 		}
// 		if(!getSon(i)->visited && getSon(i)->isDeadTriangle)
// 		{
// 			DelaunayDeadTriangle * t = static_cast<DelaunayDeadTriangle *>(getSon(i)) ;
// 			limit = std::abs(
// 			                  squareDist2D(t->getCircumCenter(),&g->getCenter())
// 			                  -(t->getRadius()+g->getRadius())*(t->getRadius()+g->getRadius())) 
// 				< .01*t->getRadius()*t->getRadius() ;
// 		}

		if( (!getSon(i)->visited && getSon(i)->isConflicting(g)) || limit)
		{
			getSon(i)->conflicts(ret,g) ;
		}
	}
	
	if(isAlive() && isTriangle && !isDeadTriangle && isConflicting(g))
	{
		ret.first.push_back(static_cast<DelaunayTriangle *>(this)) ;
	}
	
	for (size_t i  = 0 ;  i < neighbour.size() ; i++)
	{
		bool limit = false ;

// 		if(!getNeighbour(i)->visited && getNeighbour(i)->isTriangle&& !getNeighbour(i)->isDeadTriangle)
// 		{
// 			DelaunayTriangle * t = static_cast<DelaunayTriangle *>(getNeighbour(i)) ;
// 			limit = std::abs(
// 			                  squareDist2D(t->getCircumCenter(),g->getCenter())
// 			                  -(t->getRadius()+g->getRadius())*(t->getRadius()+g->getRadius())) 
// 				< .01*t->getRadius()*t->getRadius() ;
// 		}
// 		if(!getNeighbour(i)->visited && getNeighbour(i)->isDeadTriangle)
// 		{
// 			DelaunayDeadTriangle * t = static_cast<DelaunayDeadTriangle *>(getNeighbour(i)) ;
// 			limit = std::abs(squareDist2D(
// 			                               t->getCircumCenter(),&g->getCenter())
// 			                 -(t->getRadius()+g->getRadius())*(t->getRadius()+g->getRadius())) 
// 				< .01*t->getRadius()*t->getRadius() ;
// 		}

		if( (!getNeighbour(i)->visited && getNeighbour(i)->isConflicting(g)) || limit)
		{
			getNeighbour(i)->conflicts(ret, g) ;
		}
	}
}

void DelaunayTreeItem::conflicts(std::pair<std::vector<Point *>, std::vector<DelaunayTreeItem *> > & ret, const Segment *s)  
{
	std::pair<std::vector<Point *>, std::vector<DelaunayTreeItem *> > temp ;
	
	if(visited)
		return  ;
	visited = true ;
	ret.second.push_back(this) ;
	
	if(isTriangle)
	{
		std::vector<Point *> intersections ;
		Segment s0((*first), (*second)) ;
		Segment s1((*second), (*third)) ;
		Segment s2((*third), (*first)) ;
		
		if (s->intersects(s0) )
		{
			Point * i = new Point(s->intersection(s0)) ;
			intersections.push_back(i) ;
		}
		if (s->intersects(s1) )
		{
			Point * i = new Point(s->intersection(s1)) ;
			intersections.push_back(i) ;
		}
		if (s->intersects(s2) )
		{
			Point * i = new Point(s->intersection(s2)) ;

			intersections.push_back(i) ;
		}
			
		if(intersections.size() == 2)
		{
			Point *i = new Point( (*intersections[0])*0.5 + (*intersections[1])*0.5 ) ;
			ret.first.push_back(intersections[0]) ;
			ret.first.push_back(i) ;
			ret.first.push_back(intersections[1]) ;
			
		}
		if(intersections.size() == 1)
			ret.first.push_back(intersections[0]) ;
		
	}
	else if(isPlane)
	{
		Segment s0((*first), (*second)) ;
		if(s->intersects(s0))
		{
			Point * i = new Point(s->intersection(s0)) ;

			ret.first.push_back(i) ;
		}
	}
	
	if(!isAlive())
	{
		for (size_t i  = 0 ;  i < ret.first.size() ; i++)
			delete ret.first[i] ;
		
		ret.first.clear() ;
	}
	
	for (size_t i  = 0 ;  i < stepson.size() ; i++)
	{
		if( !getStepson(i)->visited) 
		{
			temp.first.clear() ;
			temp.second.clear() ;
			getStepson(i)->conflicts(temp,s) ;
			for(size_t j  = 0 ; j < temp.first.size() ;j++)
			{
				bool unique = true ;
				for(size_t k  = 0 ; k < ret.first.size() && unique == true ;k++)
				{
					if(*temp.first[j] == *ret.first[k])
						unique = false ; 
				}
				
				if(unique)
					ret.first.push_back(temp.first[j]) ;
			}
			ret.second.insert(ret.second.end(), temp.second.begin(), temp.second.end()) ;
		}
	}
	
	for (size_t i  = 0 ;  i < son.size() ; i++)
	{
		if( !getSon(i)->visited )
		{
			temp.first.clear() ;
			temp.second.clear() ;
			getSon(i)->conflicts(temp,s) ;
			for(size_t j  = 0 ; j < temp.first.size() ;j++)
			{
				bool unique = true ;
				for(size_t k  = 0 ; k < ret.first.size() && unique == true ;k++)
				{
					if(*temp.first[j] == *ret.first[k])
						unique = false ; 
				}
				
				if(unique)
					ret.first.push_back(temp.first[j]) ;
			}
			ret.second.insert(ret.second.end(), temp.second.begin(), temp.second.end()) ;
		}
	}
	
	for (size_t i  = 0 ;  i < neighbour.size() ; i++)
	{
		if( !getNeighbour(i)->visited )
		{
			temp.first.clear() ;
			temp.second.clear() ;
			getNeighbour(i)->conflicts(temp,s) ;
			for(size_t j  = 0 ; j < temp.first.size() ;j++)
			{
				bool unique = true ;
				for(size_t k  = 0 ; k < ret.first.size() && unique == true ;k++)
				{
					if(*temp.first[j] == *ret.first[k])
						unique = false ; 
				}
				
				if(unique)
					ret.first.push_back(temp.first[j]) ;
			}
			ret.second.insert(ret.second.end(), temp.second.begin(), temp.second.end()) ;
		}
	}
	
	return;
}

void DelaunayTreeItem::conflicts(std::pair<std::vector<DelaunayTreeItem *>, std::vector<DelaunayTreeItem *> >& ret,const Point *p)
{
// 	print() ;
	if(visited)
		return  ;
	visited = true ;
	ret.second.push_back(this) ;
	if(!inCircumCircle(*p))
		return  ;
	
// 	std::cerr << "stepsons" << std::endl ;
// 	Point center = (*first+*second+*third)/3. ;
	for (size_t i  = 0 ;  i < stepson.size() ; i++)
	{
		bool limit = false ;
		if(!getStepson(i)->visited && getStepson(i)->isTriangle && !getStepson(i)->isDeadTriangle)
		{
			DelaunayTriangle * t = static_cast<DelaunayTriangle *>(getStepson(i)) ;
			limit = std::abs(squareDist2D(t->getCircumCenter(),*p)-t->getRadius()*t->getRadius()) 
				< 20000.*POINT_TOLERANCE*POINT_TOLERANCE ;
		}
		if(!getStepson(i)->visited && getStepson(i)->isDeadTriangle)
		{
			DelaunayDeadTriangle * t = static_cast<DelaunayDeadTriangle *>(getStepson(i)) ;
			limit = std::abs(squareDist2D(t->getCircumCenter(),p)-t->getRadius()*t->getRadius()) 
				< 20000.*POINT_TOLERANCE*POINT_TOLERANCE ;
		}
		
		if( (!getStepson(i)->visited && getStepson(i)->inCircumCircle(*p)) || limit) 
		{
			getStepson(i)->conflicts(ret,p) ;
		}
	}
	
	
	for (size_t i  = 0 ;  i < son.size() ; i++)
	{
		bool limit = false ;
		if(!getSon(i)->visited && getSon(i)->isTriangle && !getSon(i)->isDeadTriangle)
		{
			DelaunayTriangle * t = static_cast<DelaunayTriangle *>(getSon(i)) ;
			limit = std::abs(squareDist2D(t->getCircumCenter(),*p)-t->getRadius()*t->getRadius()) 
				< 20000.*POINT_TOLERANCE*POINT_TOLERANCE ;
		}
		if(!getSon(i)->visited && getSon(i)->isDeadTriangle)
		{
			DelaunayDeadTriangle * t = static_cast<DelaunayDeadTriangle *>(getSon(i)) ;
			limit = std::abs(squareDist2D(t->getCircumCenter(),p)-t->getRadius()*t->getRadius()) 
				< 20000.*POINT_TOLERANCE*POINT_TOLERANCE ;
		}
		
		if( (!getSon(i)->visited && getSon(i)->inCircumCircle(*p)) || limit)
		{
			getSon(i)->conflicts(ret,p) ;
		}
		
	}
	
	if(isAlive())
	{
		ret.first.push_back(this) ;
	}
	
	for (size_t i  = 0 ;  i < neighbour.size() ; i++)
	{
		
		bool limit = false ;
		if(!getNeighbour(i)->visited && getNeighbour(i)->isTriangle && !getNeighbour(i)->isDeadTriangle)
		{
			DelaunayTriangle * t = static_cast<DelaunayTriangle *>(getNeighbour(i)) ;
			limit = std::abs(squareDist2D(t->getCircumCenter(),*p)-t->getRadius()*t->getRadius()) 
				< 20000.*POINT_TOLERANCE*POINT_TOLERANCE ;
		}
		if(!getNeighbour(i)->visited && getNeighbour(i)->isDeadTriangle)
		{
			DelaunayDeadTriangle * t = static_cast<DelaunayDeadTriangle *>(getNeighbour(i)) ;
			limit = std::abs(squareDist2D(t->getCircumCenter(),p)-t->getRadius()*t->getRadius()) 
				< 20000.*POINT_TOLERANCE*POINT_TOLERANCE ;
		}
// 		limit = true ;
		if( (!getNeighbour(i)->visited && getNeighbour(i)->inCircumCircle(*p)) || limit)
		{
			getNeighbour(i)->conflicts(ret,p) ;
		}
	}
}

void DelaunayTreeItem::removeNeighbour(DelaunayTreeItem * t)
{
	unsigned int* e = std::find(&neighbour[0], &neighbour[neighbour.size()], t->index) ;

	if(e !=  &neighbour[neighbour.size()])
	{
		std::valarray<unsigned int>  newneighbours(neighbour.size()-1) ;
		std::copy(&neighbour[0], e, &newneighbours[0]) ;
		std::copy(e+1, &neighbour[neighbour.size()], &newneighbours[e-&neighbour[0]]) ;
		neighbour.resize(neighbour.size()-1) ;
		neighbour = newneighbours ;
	}
}

void DelaunayTriangle::removeNeighbourhood(DelaunayTriangle * t)
{
	unsigned int* e = std::find(&neighbourhood[0], &neighbourhood[neighbourhood.size()], t->index) ;

	if(e !=  &neighbour[neighbour.size()])
	{
		std::valarray<unsigned int>  newneighbours(neighbourhood.size()-1) ;
		std::copy(&neighbourhood[0], e, &newneighbours[0]) ;
		std::copy(e+1, &neighbourhood[neighbourhood.size()], &newneighbours[e-&neighbourhood[0]]) ;
		neighbourhood.resize(neighbourhood.size()-1) ;
		neighbourhood = newneighbours ;
	}
}
	
void DelaunayTreeItem::addNeighbour(DelaunayTreeItem * t)
{
	assert(t != this) ;

	if(std::find(&neighbour[0], &neighbour[neighbour.size()],t->index) !=  &neighbour[neighbour.size()])
	{
		return ;
	}
	
	if(t->isAlive())
	{
		std::valarray<unsigned int>  newneighbours(neighbour.size()+1) ;
		std::copy(&neighbour[0], &neighbour[neighbour.size()], &newneighbours[0]) ;
		newneighbours[neighbour.size()] = t->index ;
		neighbour.resize(neighbour.size()+1) ;
		neighbour = newneighbours ;
	}
}

void DelaunayTriangle::addNeighbourhood(DelaunayTriangle * t)
{
	assert(t != this) ;

	if(std::find(&neighbourhood[0], &neighbourhood[neighbourhood.size()],t->index) !=  &neighbourhood[neighbourhood.size()])
	{
		return ;
	}
	
	if(t->isAlive())
	{
		std::valarray<unsigned int>  newneighbours(neighbourhood.size()+1) ;
		std::copy(&neighbourhood[0], &neighbourhood[neighbourhood.size()], &newneighbours[0]) ;
		newneighbours[neighbourhood.size()] = t->index ;
		neighbourhood.resize(neighbourhood.size()+1) ;
		neighbourhood = newneighbours ;
	}
}

DelaunayTriangle * DelaunayTriangle::getNeighbourhood(size_t i) const
{
	return static_cast<DelaunayTriangle *>(tree->tree[neighbourhood[i]]) ;
}

DelaunayTreeItem * DelaunayTreeItem::getNeighbour(size_t i) const
{
	return tree->tree[neighbour[i]] ;
}

DelaunayTreeItem * DelaunayTreeItem::getSon(size_t i) const
{
	return tree->tree[son[i]] ;
}

DelaunayTreeItem * DelaunayTreeItem::getStepson(size_t i) const
{
	return tree->tree[stepson[i]] ;
}

void DelaunayTreeItem::addStepson(DelaunayTreeItem * s)
{
	if(s == this)
		return ;
	
	std::valarray<unsigned int>  newstepson(stepson.size()+1) ;
	std::copy(&stepson[0], &stepson[stepson.size()], &newstepson[0]) ;
	newstepson[stepson.size()] = s->index ;
	stepson.resize(stepson.size()+1) ;
	stepson = newstepson ;
	s->setStepfather(this) ;
	addNeighbour(s) ;
}

void DelaunayTreeItem::addSon(DelaunayTreeItem * s)
{
	std::valarray<unsigned int>  newson(son.size()+1) ;
	std::copy(&son[0], &son[son.size()], &newson[0]) ;
	newson[son.size()] = s->index ;
	son.resize(son.size()+1) ;
	son = newson ;
}
	
void DelaunayTreeItem::removeSon(DelaunayTreeItem * t)
{
	if(!son.size())
		return ;
	unsigned int* e = std::find(&son[0], &son[son.size()], t->index) ;
	
	if(e !=  &son[son.size()])
	{
		std::valarray<unsigned int>  newson(son.size()-1) ;
		std::copy(&son[0], e, &newson[0]) ;
		std::copy(e+1, &son[son.size()], &newson[e-&son[0]]) ;
		son.resize(son.size()-1) ;
		son = newson ;
	}
}

void DelaunayTreeItem::removeStepson(DelaunayTreeItem * t)
{
	if(stepson.size() == 0)
		return ;
	unsigned int* e = std::find(&stepson[0], &stepson[stepson.size()], t->index) ;
	if(e !=  &stepson[stepson.size()])
	{
		std::valarray<unsigned int>  newstepson(stepson.size()-1) ;
		std::copy(&stepson[0], e, &newstepson[0]) ;
		std::copy(e+1, &stepson[stepson.size()], &newstepson[e-&stepson[0]]) ;
		stepson.resize(stepson.size()-1) ;
		stepson = newstepson ;
	}
}


void DelaunayTreeItem::erase(const Point * p)
{
	dead = true ;
	m_k  = p ;
}

void DelaunayTreeItem::kill(const Point * p)
{
	dead = true ;
	m_k  = p ;
	
	for(size_t i = 0 ; i < this->neighbour.size() ; i++)
	{
		this->getNeighbour(i)->removeNeighbour(this) ;
	}
}
	
bool DelaunayTreeItem::isAlive() const
{
	return !this->dead ;
}
	


void DelaunayTreeItem::setStepfather(DelaunayTreeItem * s)
{
	stepfather = s ;
	addNeighbour(s) ;
}
	
void DelaunayTreeItem::clearVisited()
{
	visited = false ;
}

DelaunayTriangle::DelaunayTriangle(DelaunayTree *t, DelaunayTreeItem * father,  Point *p0,  Point *p1, Point *p2,  Point * c) : TriElement(p0, p1, p2), DelaunayTreeItem(t, father, c)
{
	first = &getBoundingPoint(0) ;
	second = &getBoundingPoint(1) ;
	third = &getBoundingPoint(2) ;
	
	assert(in(this->getCenter())) ;
	
	isPlane = false ;
	isTriangle = true ;
	assert(first->id > -1) ;
	assert(second->id > -1) ;
	assert(third->id > -1) ;
	cachedGPs = NULL ;
}

DelaunayTriangle::DelaunayTriangle() : DelaunayTreeItem(NULL, NULL, NULL)
{
	first = &getBoundingPoint(0) ;
	second = &getBoundingPoint(1) ;
	third = &getBoundingPoint(2) ;
	
	isPlane = false ;
	isTriangle = true ;
	cachedGPs = NULL ;
}

DelaunayDeadTriangle::DelaunayDeadTriangle( DelaunayTriangle * parent) : DelaunayTreeItem(*parent),
center(parent->getCircumCenter()), 
radius(parent->getRadius()), sqradius(radius*radius)
{
	index = parent->index ;
	tree = parent->tree ;
	
	stepson.resize(parent->stepson.size()) ;
	stepson = parent->stepson ;
	
	neighbour.resize(parent->neighbour.size()) ;
	neighbour = parent->neighbour ;
	
	son.resize(parent->son.size()) ;
	son = parent->son ;
	
	father = parent->father ;
	stepfather = parent->stepfather ;

	isTriangle = true ;
	isPlane = false ;
	isDeadTriangle = true ;
	visited =false ;
	
	first = parent->first ;
	second = parent->second ;
	third = parent->third ;
	dead = true ;
}
	

double DelaunayDeadTriangle::getRadius() const
{
	return radius ;
}

const Point * DelaunayDeadTriangle::getCircumCenter() const
{
	return &center ;
}

bool DelaunayDeadTriangle::inCircumCircle(const Point & p) const
{


	if(p.x > center.x+1.01*radius)
		return false ;
	if(p.x < center.x-1.01*radius)
		return false ;
	if(p.y > center.y+1.01*radius)
		return false ;
	if(p.y < center.y-1.01*radius)
		return false ;

	if(squareDist2D(center, p) < .99*sqradius)
		return true ;

	double delta = POINT_TOLERANCE*radius ;
	Point a(p) ; a.x += delta ; a.y += delta ;
	Point c(p) ; c.x += delta ; c.y -= delta ; 
	Point e(p) ; e.x -= delta ; e.y += delta ;
	Point h(p) ; h.x -= delta ; h.y -= delta ; 
	return  squareDist2D(center, a) < sqradius 
		&&  squareDist2D(center, c) < sqradius
		&&  squareDist2D(center, e) < sqradius
		&&  squareDist2D(center, h) < sqradius;	
	
}

bool DelaunayDeadTriangle::isNeighbour( const DelaunayTreeItem * t) const
{
	size_t cv = this->numberOfCommonVertices(t) ;
	return (cv == 2);
}

bool DelaunayDeadTriangle::isVertex(const Point *p) const
{
	return (*p == *first) || (*p == *second) || (*p == *third)  ;
}

bool DelaunayDeadTriangle::isVertexByID(const Point *p) const
{
	return p == first || p == second || p == third  ;
}

bool DelaunayDeadTriangle::in( const Point & p) const
{
	return Triangle(third, first, second).in(p) ;
}

void DelaunayDeadTriangle::print() const
{

	std::cout << "(" << first->x << ", " << first->y <<  ") " ;
	std::cout << "(" << second->x << ", " << second->y <<  ") " ;
	std::cout << "(" << third->x << ", " << third->y <<  ") " ;
	std::cout <<  ":: "<< isAlive() << std::endl ;
}
	
inline std::pair<Point*, Point*> DelaunayDeadTriangle::commonEdge(const DelaunayTreeItem * t) const
{
	if(t == this)
	{
		return std::pair<Point*, Point*>(NULL, NULL) ;
	}
	
	if(t->isTriangle)
	{
		if(this->isVertex(t->first))
		{
			if(this->isVertex(t->second))
				return std::pair< Point*,  Point*>(t->first , t->second) ;
			if(this->isVertex(t->third))
				return std::pair< Point*,  Point*>(t->first , t->third) ;
			
			return std::pair<Point*, Point*>(NULL, NULL) ;
		}
		
		if(this->isVertex(t->second))
		{
			if(this->isVertex(t->third))
				return std::pair< Point*,  Point*>(t->third , t->second) ;
			
			return std::pair<Point*, Point*>(NULL, NULL) ;
		}
		
	}
	else
	{
		if(isAligned(t->first, first, second) && isAligned(t->second, first, second))
			return std::pair< Point*,  Point*>(first, second ) ;
		if(isAligned(t->first, third, second) && isAligned(t->second, third, second))
			return std::pair< Point*,  Point*>(third, second ) ;
		if(isAligned(t->first, third, first) && isAligned(t->second, third, first))
			return std::pair< Point*,  Point*>(first, third ) ;
	}
	
	return std::pair< Point*,  Point*>(NULL, NULL) ;
}


DelaunayTriangle::~DelaunayTriangle()
{
	delete cachedGPs ;
	for(size_t i = 0 ; i <  cachedElementaryMatrix.size() ; i++)
	{
		cachedElementaryMatrix[i].clear() ;
	}
}

DelaunayDemiPlane::~DelaunayDemiPlane()
{
// 	first = NULL ;
// 	second = NULL ;
// 	third = NULL ;
}
	
inline bool DelaunayTriangle::isVertex(const Point * p) const
{
	return (*p == *first || *p == *second || *p == *third) ;
}

bool DelaunayTriangle::isVertexByID(const Point * p) const
{
	return p == first || p == second || p == third ;
}
	

inline bool DelaunayTriangle::hasVertexByID(const std::valarray<Point *> * p) const
{
	for(size_t i = 0 ; i < p->size() ; i++)
	{
		if((*p)[i]->id == first->id || (*p)[i]->id == second->id || (*p)[i]->id == third->id)
			return true ;
	}
	return false ;
}

std::pair<Point*, Point*> DelaunayTriangle::commonEdge(const DelaunayTreeItem * t) const
{
	if(t == this)
	{
		return std::pair<Point*, Point*>(NULL, NULL) ;
	}
	
	if(t->isTriangle || t->isDeadTriangle)
	{
		if(this->isVertex(t->first))
		{
			if(this->isVertex(t->second))
				return std::pair< Point*,  Point*>(t->first , t->second) ;
			if(this->isVertex(t->third))
				return std::pair< Point*,  Point*>(t->first , t->third) ;
			
			return std::pair<Point*, Point*>(NULL, NULL) ;
		}
		
		if(this->isVertex(t->second))
		{
			if(this->isVertex(t->third))
				return std::pair< Point*,  Point*>(t->third , t->second) ;
			
			return std::pair<Point*, Point*>(NULL, NULL) ;
		}
		
	}
	else if (numberOfCommonVertices(t) == 2)
	{
		if(!isAligned(first, t->first, t->second))
			return std::pair< Point*,  Point*>(third, second ) ;
		if(!isAligned(second, t->first, t->second))
			return std::pair< Point*,  Point*>(third, first ) ;
		if(!isAligned(t->third, t->first, t->second))
			return std::pair< Point*,  Point*>(first, second ) ;
		
		//now, this means they are all "aligned"
		Segment seg(*t->first,*t->second) ;
		double d0 = squareDist2D(*first, seg.project(*first)) ; 
		double d1 = squareDist2D(*second, seg.project(*second)) ; 
		double d2 = squareDist2D(*third, seg.project(*third)) ; 
		if(d0 >= d1 && d0 >=d2)
			return std::pair< Point*,  Point*>(third, second ) ;
		else if(d1 >= d0 && d1 >=d2)
			return std::pair< Point*,  Point*>(third, first ) ;
		
		return std::pair< Point*,  Point*>(first, second ) ;
	}
	
	return std::pair< Point*,  Point*>(NULL, NULL) ;
}

std::pair< Point*,  Point*> DelaunayDemiPlane::commonEdge(const DelaunayTreeItem * t) const
{
	if(t->isTriangle)
	{
		Point test = (*first)*0.5 + (*second)*0.5 ;
		
		if(isAligned(test, (*t->first), (*t->second)))
		{
			return std::pair< Point*,  Point*>(t->first, t->second ) ;
		}
		if(isAligned(test, (*t->third), (*t->second)))
		{
			return std::pair< Point*,  Point*>(t->third, t->second ) ;
		}
		if(isAligned(test, (*t->first), (*t->third)))
		{
			return std::pair< Point*,  Point*>( t->first, t->third ) ;
		}
	}
	else
	{
		if((t->first == first && t->second == second) ||
		   (t->first == second &&  t->second == first))
			return std::pair< Point*,  Point*>(t->first , t->second) ;
	}
	
	return std::pair< Point*,  Point*>(NULL, NULL) ;
}
	

void DelaunayDemiPlane::merge(DelaunayDemiPlane * p)
{
	if(isAlive() && p != this && p->isAlive())
	{

		if( (first == p->first && second == p->second ) || (first == p->second && second == p->first) )
		{
			for(size_t i = 0 ; i <  p->neighbour.size() ; i++)
			{
				for(size_t j = 0 ; j <  neighbour.size() ; j++)
				{
					if(p->getNeighbour(i)->isNeighbour(getNeighbour(j)))
						makeNeighbours(getNeighbour(j), p->getNeighbour(i)) ;
				}
			}
			p->kill(first) ;
			kill(first) ;
			return ;
		}
		
		if(isAligned(first, second, p->first) && isAligned(first, second, p->second))
		{

			for(size_t i = 0 ; i <  p->neighbour.size() ; i++)
			{
				makeNeighbours(this, p->getNeighbour(i)) ;
			}
			for(size_t i = 0 ; i <  p->son.size() ; i++)
			{
				p->getSon(i)->father = this ;
				addSon(p->getSon(i)) ;
			}
			for(size_t i = 0 ; i <  p->stepson.size() ; i++)
			{
				p->getStepson(i)->stepfather = this ;
				addStepson(p->getStepson(i)) ;
				makeNeighbours(this, p->getStepson(i)) ;
			}
			p->father->addSon(this) ;
			p->addSon(this) ;
			p->kill(first) ;
		}
	}
}


std::pair< Point*,  Point*> DelaunayDeadTriangle::nearestEdge(const Point & p) const
{
	std::map<double, Point> cen ;
	Point c0(((*first) + (*second))/2.) ;
	cen[squareDist2D(c0, p)] = c0 ;
	Point c1(((*third) + (*second))/2.) ;
	cen[squareDist2D(c1, p)] = c1 ;
	Point c2(((*third) + (*first))/2.) ;
	cen[squareDist2D(c2, p)] = c2 ;
	
	if(cen.begin()->second == c0)
		return std::pair< Point*,  Point*>(first, second) ;
	if(cen.begin()->second == c1)
		return std::pair< Point*,  Point*>(third, second) ;
	if(cen.begin()->second == c2)
		return std::pair< Point*,  Point*>(third, first) ;
	
	return std::pair< Point*,  Point*>(first, second) ;
}

std::pair< Point*,  Point*> DelaunayTriangle::nearestEdge(const Point & p) const
{
	std::map<double, Point> cen ;
	Point c0(((*first) + (*second))/2.) ;
	cen[squareDist2D(c0, p)] = c0 ;
	Point c1(((*third) + (*second))/2.) ;
	cen[squareDist2D(c1, p)] = c1 ;
	Point c2(((*third) + (*first))/2.) ;
	cen[squareDist2D(c2, p)] = c2 ;
	
	if(cen.begin()->second == c0)
		return std::pair< Point*,  Point*>(first, second) ;
	if(cen.begin()->second == c1)
		return std::pair< Point*,  Point*>(third, second) ;
	if(cen.begin()->second == c2)
		return std::pair< Point*,  Point*>(third, first) ;
	
	return std::pair< Point*,  Point*>(first, second) ;
}
	
bool DelaunayTriangle::inCircumCircle(const Point &p) const 
{
// 	return  (circumCenter.x -p.x)*(circumCenter.x -p.x) + (circumCenter.y -p.y)*(circumCenter.y -p.y)< radius*radius -POINT_TOLERANCE  ;
	return  Triangle::inCircumCircle(p) ;
}
	
void DelaunayTriangle::insert(std::vector<DelaunayTreeItem *> &ret, Point *p,  Star* s)
{
	if (visited)
		return ;

	visited = true ;
	for (size_t i = 0 ; i < neighbour.size() ; i++)
	{
		if(this->numberOfCommonVertices(getNeighbour(i)) == 2)
		{
			std::pair< Point*,  Point*> pp = this->commonEdge(getNeighbour(i)) ;

			if (!getNeighbour(i)->visited && !getNeighbour(i)->inCircumCircle(*p))
			{

				if(Triangle(*p, *pp.first, *pp.second).area() > getRadius()*std::numeric_limits<double>::epsilon())
				{
					DelaunayTriangle *ss = new DelaunayTriangle(this->tree, this, p, pp.first, pp.second, p) ;
					addSon(ss) ;
	
					getNeighbour(i)->addStepson(ss) ;
					ret.push_back(ss) ;
				}

			}
		}

	}
	this->kill(p) ;
	s->updateNeighbourhood() ;

}
	

bool DelaunayTriangle::isNeighbour(const DelaunayTreeItem * t)  const
{
// 	if(t->isPlane)
// 		return ((int)isAligned(t->first, t->second, first) + (int)isAligned(t->first, t->second, second) + (int)isAligned(t->first, t->second, third)) == 2;
	
	return this->numberOfCommonVertices(t) == 2 ;
}

bool DelaunayTriangle::isInNeighbourhood(const DelaunayTriangle * t) const
{
	return this->numberOfCommonVertices(t) > 0 ;
}

void DelaunayTriangle::print() const
{
	for(size_t i = 0 ; i < this->getBoundingPoints().size() ; i++)
	{
		std::cerr << "(" << getBoundingPoint(i).x <<  ";" << getBoundingPoint(i).y <<") " ;
	}
	std::cerr <<  ":: "<< isAlive() << std::endl ;
}


// void DelaunayTriangle::displace(Vector * eps)
// {
// 	(*eps)[first->id*2]+=(*eps)[first->id*2] ;
// 	(*eps)[first->id*2+1]+=(*eps)[first->id*2+1] ;
// }

DelaunayDemiPlane::DelaunayDemiPlane(DelaunayTree *t, DelaunayTreeItem * father,  Point  * _begin,  Point  * _end,  Point  * p,  Point * c) : DelaunayTreeItem(t, father, c)
{
	second  = _end ;
	first = _begin ;
	third = p ;
	dead = false ;
	vector =(*second)- (*first) ;
	Point pseudonormal = (*third) - (*first);
	direction = (vector.x*pseudonormal.y - vector.y*pseudonormal.x) ;
	isPlane =true ;
	isTriangle = false ;
}

std::pair< Point*,  Point*> DelaunayDemiPlane::nearestEdge(const Point & p) const
{
	return std::pair< Point*,  Point*>(first, second) ;
}
	
bool DelaunayDemiPlane::inCircumCircle(const Point &p) const
{
	return (fma(vector.x,(p.y - first->y), - vector.y*(p.x - first->x)) * direction < -10*std::numeric_limits<double>::epsilon()) ;
}
	
bool DelaunayDemiPlane::isVertex(const Point *p) const
{
	return ( (*p) == (*first) || (*p) == (*second) ) || isAligned(p, first, second) ;
}
	

bool  DelaunayDemiPlane::isNeighbour(const DelaunayTreeItem * t) const
{
	
	if(t->isTriangle)
	{
		return t->isNeighbour(this) ;
	}
	else
	{
		return (first ==  t->first || first == t->first ||
		        second == t->first || second == t->second) ;
	}
return false ;
}

void DelaunayDemiPlane::insert(std::vector<DelaunayTreeItem *> & ret, Point *p, Star *s)
{
	if (visited)
		return ;
	visited = true ;
	
	std::vector<Point*> lims ;
	
	for(size_t i = 0 ; i < neighbour.size() ; i++)
	{
		if(numberOfCommonVertices(getNeighbour(i)) == 2)
		{
			std::pair< Point*,  Point*> pp = getNeighbour(i)->commonEdge(this) ;
			if(std::find(lims.begin(), lims.end(), pp.first) != lims.end())
				lims.erase(std::find(lims.begin(), lims.end(), pp.first)) ;
			else
				lims.push_back(pp.first) ;
			
			if(std::find(lims.begin(), lims.end(), pp.second) != lims.end())
				lims.erase(std::find(lims.begin(), lims.end(), pp.second)) ;
			else
				lims.push_back(pp.second) ;
			
			if (!getNeighbour(i)->visited && !getNeighbour(i)->inCircumCircle(*p))
			{

				assert(getNeighbour(i)->isNeighbour(this) );
				
				if(Triangle(*p, *pp.first, *pp.second).area() > static_cast<DelaunayTriangle *>(getNeighbour(i))->getRadius()*std::numeric_limits<double>::epsilon())
				{
					DelaunayTriangle *ss = new DelaunayTriangle(this->tree, this, p, pp.first, pp.second , p) ;
					ss->visited = true ;
					addSon(ss) ;
					
					getNeighbour(i)->addStepson(ss) ;
		
					ret.push_back(ss) ;
				}
			}
		}
	}
	

	if(lims.empty())
	{
		lims.push_back(first) ;
		lims.push_back(second) ;
	}
	
	if(!isAligned(lims[0], p, lims[1]))
	{
		this->kill(p) ;
		DelaunayDemiPlane *p0 = new DelaunayDemiPlane(this->tree,this, lims[0], p, lims[1], p) ;

		DelaunayDemiPlane *p1 = new DelaunayDemiPlane(this->tree,this, lims[1], p, lims[0], p) ;

		addSon(p0) ;
		addSon(p1) ;
		
		ret.push_back(p0) ;
		ret.push_back(p1) ;
		
		for(size_t i = 0 ; i < son.size()-2 ; i++)
		{
			if(getSon(i)->isNeighbour(p0))
				getSon(i)->addStepson(p0) ;
			
			if(getSon(i)->isNeighbour(p1))
				getSon(i)->addStepson(p1) ;
			
			for(size_t j = i ; j < son.size()-2 ; j++)
				if(getSon(i)->isNeighbour(getSon(j)))
					makeNeighbours(getSon(i), getSon(j)) ;
		}
	}
	else
	{
		this->kill(p) ;
		DelaunayDemiPlane *p0 = new DelaunayDemiPlane(this->tree,this, lims[0], third, lims[1], p) ;
		
		DelaunayDemiPlane *p1 = new DelaunayDemiPlane(this->tree,this, lims[1], third, lims[0], p) ;
		
		addSon(p0) ;
		addSon(p1) ;
		
		ret.push_back(p0) ;
		ret.push_back(p1) ;
		
		for(size_t i = 0 ; i < son.size()-2 ; i++)
		{
			if(getSon(i)->isNeighbour(p0))
				getSon(i)->addStepson(p0) ;
			
			if(getSon(i)->isNeighbour(p1))
				getSon(i)->addStepson(p1) ;
			
			for(size_t j = i ; j < son.size()-2 ; j++)
				if(getSon(i)->isNeighbour(getSon(j)))
					makeNeighbours(getSon(i), getSon(j)) ;
		}
	}
	
	s->updateNeighbourhood() ;

	return  ;
}

void DelaunayDemiPlane::print() const 
{
	std::cerr << "###############(" << first->x << ", " << first->y << ") (" <<
		second->x << ", " << second->y << ")" << "X (" << third->x << ", " << third->y << ") :: " << isAlive()  << std::endl ;
}

void makeNeighbours(DelaunayTreeItem *t0, DelaunayTreeItem *t1 ) 
{
	if(t0 == t1 || !t0->isAlive() || !t1->isAlive())
		return ;
	if(t0->isPlane && t1->isPlane)
	   return ;
	if(!t0->isNeighbour(t1))
		return ;
	t0->addNeighbour(t1) ;
	t1->addNeighbour(t0) ;
} 

void updateNeighbours(std::vector<DelaunayTreeItem *> * t)
{
	for(size_t i = 0 ; i < t->size() ; i++)
	{
		for(size_t j = i ; j < t->size() ; j++)
		{
			if((*t)[i]->isNeighbour((*t)[j]))
			{
				makeNeighbours( (*t)[i],(*t)[j] ) ;
			}
		}
	}
} 

DelaunayRoot::DelaunayRoot(DelaunayTree *tt, Point * p0, Point * p1, Point * p2) : DelaunayTreeItem(tt, NULL, NULL)
{
	isPlane = false ;
	isTriangle =false ;
	isDeadTriangle =false ;
	this->father = NULL ;
	DelaunayTriangle *t = new DelaunayTriangle(tt, this, p0, p1, p2, NULL) ;
	DelaunayDemiPlane * pl0 = new DelaunayDemiPlane(tt,this, p0, p1, p2, NULL);
	DelaunayDemiPlane * pl1 = new DelaunayDemiPlane(tt,this, p0, p2, p1, NULL);
	DelaunayDemiPlane * pl2 = new DelaunayDemiPlane(tt,this, p1, p2, p0, NULL);
	
	makeNeighbours(t,pl0 ) ;
	makeNeighbours(t,pl1) ;
	makeNeighbours(t,pl2) ;

	addSon(t) ;
	addSon(pl0) ;
	addSon(pl1) ;
	addSon(pl2) ;
	kill(p0) ;
}
	
void DelaunayRoot::print() const
{
	std::cerr << "I am root !" << std::endl ;
}
	
bool DelaunayRoot::isVertex(const Point *p) const
{
	return false ;
}
	
bool DelaunayRoot::inCircumCircle(const Point & p) const 
{
	return true ;
}
	
std::pair< Point*,  Point*> DelaunayRoot::nearestEdge(const Point & p) const
{
	assert(false) ;
	return std::pair< Point*,  Point*>(NULL, NULL) ;
}
	
void DelaunayRoot::insert(std::vector<DelaunayTreeItem *> & ret,Point *p, Star *s)
{
	
	for (size_t i  = 0 ;  i < son.size() ; i++)
	{
		
		std::vector<DelaunayTreeItem *> temp ;
		getSon(i)->insert(temp, p, s) ;
		
		for(size_t j = 0 ; j< temp.size() ; j++)
		{
			ret.push_back(temp[j]) ;
		}
		
	}
	
	updateNeighbours(&ret) ;
	return  ;
}



 void DelaunayRoot::conflicts(std::pair<std::vector<DelaunayTriangle *>, std::vector<DelaunayTreeItem *> > & ret, const Geometry *g)
{

	visited = true ;
	
	for (size_t i  = 0 ;  i < son.size() ; i++)
	{
		std::pair<std::vector<DelaunayTriangle *>, std::vector<DelaunayTreeItem *> > temp  ;
		getSon(i)->conflicts(temp,g) ;
		ret.first.insert(ret.first.end(),temp.first.begin(), temp.first.end()) ;
		ret.second.insert(ret.second.end(),temp.second.begin(), temp.second.end()) ;
	}
	
	return ;
}


 void DelaunayRoot::conflicts(std::pair<std::vector<DelaunayTreeItem *>, std::vector<DelaunayTreeItem *> > & ret,const Point *p)
{

	visited = true ;
	
	for (size_t i  = 0 ;  i < son.size() ; i++)
	{
		std::pair<std::vector<DelaunayTreeItem *>, std::vector<DelaunayTreeItem *> > temp  ;
		getSon(i)->conflicts(temp,p) ;
		ret.first.insert(ret.first.end(),temp.first.begin(), temp.first.end()) ;
		ret.second.insert(ret.second.end(),temp.second.begin(), temp.second.end()) ;
	}
	
	return  ;
}

Star::Star(std::vector<DelaunayTreeItem *> *t, const Point *p)
{
	treeitem.insert(treeitem.end(), t->begin(), t->end()) ;
	for(size_t i = 0 ; i < t->size() ; i++)
	{
		if((*t)[i]->isTriangle)
		{
			this->edge.push_back((*t)[i]->first) ;
			this->edge.push_back((*t)[i]->second) ;
			this->edge.push_back((*t)[i]->third) ;
		}
		else if((*t)[i]->isPlane)
		{
			this->edge.push_back((*t)[i]->first) ;
			this->edge.push_back((*t)[i]->second) ;
		}
	}
	std::sort(edge.begin(), edge.end()) ;
	edge.erase(unique(edge.begin(), edge.end()), edge.end()) ;
}
	
size_t Star::size()
{
	return edge.size() ;
}
	
const Point * Star::getEdge(size_t i) const
{
	assert(i < edge.size()) ;
	return this->edge[i] ;
}


DelaunayTree::DelaunayTree(Point * p0, Point *p1, Point *p2)
{
	neighbourhood = false ;
	this->global_counter = 3;
	p0->id = 0 ; p1->id = 1 ; p2->id = 2 ;
	DelaunayRoot *root = new DelaunayRoot( this, p0, p1, p2) ;
	plane.push_back(static_cast<DelaunayDemiPlane *>(root->getSon(1))) ;
	plane.push_back(static_cast<DelaunayDemiPlane *>(root->getSon(2))) ;
	plane.push_back(static_cast<DelaunayDemiPlane *>(root->getSon(3))) ;
}
	
DelaunayTree::~DelaunayTree() 
{ 	
	std::valarray<Point *> nularray(0) ;
	for(size_t i = 0 ;  i < this->tree.size() ; i++)
	{

		if(this->tree[i] && this->tree[i]->isTriangle && !this->tree[i]->isDeadTriangle )
		{
			DelaunayTriangle * t = static_cast<DelaunayTriangle *>(tree[i]) ;
			
			t->setBoundingPoints(nularray) ;
		}

	}
	
	for(size_t i = 0 ;  i < this->tree.size() ; i++)
	{
		delete this->tree[i] ;
	}

	for(size_t i = 0 ; i < additionalPoints.size() ; i++)
		delete additionalPoints[i] ;
	
} ;

void DelaunayTree::insertIf( Point *p, std::vector<SamplingCriterion *> v, double minScore )
{
	std::vector<DelaunayTreeItem *> cons = this->conflicts(p) ;
	
	//We store all the pointers to the affected elements of the tree, and make copies of those
	std::vector<DelaunayTreeItem *> backup = tree ;
	
	
	Star * s = new Star(&cons, p) ;
	
	std::vector<DelaunayTreeItem *> ret ;
	
	for(size_t i = 0 ; i < cons.size() ; i++)
	{
		std::vector<DelaunayTreeItem *> temp ;
		cons[i]->insert(temp,p, s) ;
		ret.insert(ret.end(), temp.begin(), temp.end()) ;
	}
	
	for(size_t i =  0 ; i < ret.size() ; i++)
	{
		if(ret[i]->isTriangle)
		{
			double score = 0;
			for(size_t j = 0 ; j < v.size() ; j++)
				if (v[j]->meetsCriterion((DelaunayTriangle *)(ret[i])))
					score += 1/v.size() ;
			
			if (score < minScore)
			{
				tree = backup ;	
				delete s ;
				for(size_t j =  0 ; j < ret.size() ; j++)
				{
					delete ret[j] ;
				}

				for(size_t j = 0 ; j < cons.size() ; j++)
					cons[j]->clearVisited() ;
				
				return ;
			}
			else
			{
				p->id =this->global_counter++ ;
			}
		}
	}
	
	bool weGotPlanes = false ;
	
	for(size_t j = 0 ; j< ret.size() ; j++)
	{
		if(ret[j]->isAlive() && ret[j]->isPlane )
		{
			weGotPlanes = true ;
			
			if(ret[j]->neighbour.size() > 1 || ret[j]->neighbour.size() == 0)
			{
				ret[j]->kill(p) ;
			}
			else
			{
				plane.push_back((DelaunayDemiPlane*)(ret[j])) ;
			}
		}
	}
	
	if(weGotPlanes)
	{
		for(int k = 0 ; k< (int)plane.size()-1 ; k++)
		{
			for(int j = k ; j< (int)plane.size() ; j++)
			{
				plane[j]->merge(plane[k]) ;
			}
		}
	}
	
	for(size_t i = 0 ; i < ret.size() ; i++)
		ret[i]->clearVisited() ;

	for(size_t i = 0 ; i < cons.size() ; i++)
		cons[i]->clearVisited() ;
	
	for(size_t i = 0 ; i < ret.size() ; i++)
	{
		
		assert(ret[i]->isPlane  || ret[i]->neighbour.size() == 3) ;
		
	}
	
	std::vector<DelaunayDemiPlane *> * hull = this->getConvexHull() ;
	plane.clear() ;
	plane.insert(plane.end(), hull->begin(), hull->end()) ;
	delete hull ;
	delete s ;
}

void DelaunayTree::insert(Segment *s)
{
	std::vector<Point *> cons = this->conflicts(s);
	
	this->insert(new Point(s->first())) ;
	
	for(size_t i = 0 ;  i < cons.size() ; i++)
			this->insert(cons[i]) ;
		
	this->insert(new Point(s->second())) ;
}

void DelaunayTree::insert(Point *p)
{
	std::vector<DelaunayTreeItem *> cons = this->conflicts(p) ;
		
	if(cons.empty())
	{
// 		std::cerr << "Failed insertion : in nothing !" << std::endl ;
		return ;
	}
	
	neighbourhood = false ;
	
	for(size_t i = 0 ; i < cons.size() ; i++)
	{
		if(cons[i]->isVertex(p)) 
		{
// 			std::cerr << "Vertex collision !" << std::endl ;
			return ;
		}
	}
	
	p->id = this->global_counter++ ;

	Star * s = new Star(&cons, p) ;
	
	std::vector<DelaunayTreeItem *> ret ;
	
	for(size_t i = 0 ; i < cons.size() ; i++)
	{
		std::vector<DelaunayTreeItem *> temp ;
		cons[i]->insert(temp,p, s) ;
		ret.insert(ret.end(), temp.begin(), temp.end()) ;
	}
	
	s->updateNeighbourhood() ;

	bool weGotPlanes = false ;
	
	for(size_t j = 0 ; j< ret.size() ; j++)
	{
		if(ret[j]->isAlive() && ret[j]->isPlane )
		{
			weGotPlanes = true ;
			plane.push_back((DelaunayDemiPlane*)(ret[j])) ;
		}
	}
	
	if(weGotPlanes)
	{
		for(int k = 0 ; k< (int)plane.size()-1 ; k++)
		{
			for(int j = k ; j< (int)plane.size() ; j++)
			{
				plane[j]->merge(plane[k]) ;
			}
		}
	}

	for(size_t i = 0 ; i < ret.size() ; i++)
	{
		ret[i]->clearVisited() ;

	}
	for(size_t i = 0 ; i < cons.size() ; i++)
	{
		cons[i]->clearVisited() ;
	}
	

	std::vector<DelaunayDemiPlane *> * hull = this->getConvexHull() ;
	plane.clear() ;
	plane.insert(plane.end(), hull->begin(), hull->end()) ;

	for(size_t i = 0 ; i < cons.size() ; i++)
	{
		if(!cons[i]->isAlive() && cons[i]->isTriangle && !cons[i]->isDeadTriangle)
		{
			DelaunayDeadTriangle* dt = new DelaunayDeadTriangle(static_cast<DelaunayTriangle *>(cons[i])) ;
			dt->clearVisited() ;
			std::valarray<Point *> nularray(0) ;
			static_cast<DelaunayTriangle *>(cons[i])->setBoundingPoints(nularray) ;
			tree[cons[i]->index] = dt ;
			delete cons[i] ;
		}
	}
	delete hull ;
	delete s ;
	
}



std::vector<DelaunayTreeItem *> DelaunayTree::conflicts( const Point *p) const
{
	std::pair< std::vector<DelaunayTreeItem *>,std::vector<DelaunayTreeItem *> > cons ;
	this->tree[0]->conflicts(cons,p) ;
	
	for(size_t i = 0 ; i < plane.size() ; i++)
	{
	
		if(!plane[i]->visited)
		{
			std::pair< std::vector<DelaunayTreeItem *>,std::vector<DelaunayTreeItem *> > temp ;
			plane[i]->conflicts(temp,p) ;
			
			cons.first.insert(cons.first.end(), temp.first.begin(),temp.first.end()) ;
			cons.second.insert(cons.second.end(), temp.second.begin(),temp.second.end()) ;
			
			
			for(size_t j = 0 ; j < plane[i]->neighbour.size() ; j++)
			{
				if(!plane[i]->getNeighbour(j)->visited)
				{
					std::pair< std::vector<DelaunayTreeItem *>, std::vector<DelaunayTreeItem *> > temp ;
					plane[i]->getNeighbour(j)->conflicts(temp,p) ;
					
					cons.first.insert(cons.first.end(), temp.first.begin(),temp.first.end()) ;
					cons.second.insert(cons.second.end(), temp.second.begin(),temp.second.end()) ;
				}
			}
			
		}
		
	}

	std::stable_sort(cons.first.begin(), cons.first.end()) ;
	std::vector<DelaunayTreeItem *>::iterator e = std::unique(cons.first.begin(), cons.first.end()) ;
	cons.first.erase(e, cons.first.end()) ;
	
	
	for(size_t i = 0 ; i < cons.second.size() ; i++)
	{
		cons.second[i]->clearVisited() ;
	}
	
	std::vector<DelaunayTreeItem *> ret ;
	ret.insert(ret.end(), cons.first.begin(), cons.first.end()) ;	
	return ret ;
}

std::vector<Point *> DelaunayTree::conflicts( const Segment *s) const
{
	std::pair< std::vector<Point *>,std::vector<DelaunayTreeItem *> > cons;
	
	this->tree[0]->conflicts(cons,s) ;
	
	for(size_t i = 0 ; i < plane.size() ; i++)
	{
		
		if(!plane[i]->visited)
		{
			std::pair< std::vector<Point *>,std::vector<DelaunayTreeItem *> > temp ;
			plane[i]->conflicts(temp,s) ;
			
			cons.first.insert(cons.first.end(), temp.first.begin(),temp.first.end()) ;
			cons.second.insert(cons.second.end(), temp.second.begin(),temp.second.end()) ;			
			
			for(size_t j = 0 ; j < plane[i]->neighbour.size() ; j++)
			{
				std::pair< std::vector<Point *>, std::vector<DelaunayTreeItem *> > temp ;
				plane[i]->getNeighbour(j)->conflicts(temp,s) ;
				
				cons.first.insert(cons.first.end(), temp.first.begin(),temp.first.end()) ;
				cons.second.insert(cons.second.end(), temp.second.begin(),temp.second.end()) ;
			}
			
		}
		
	}
	
	for(size_t i = 0 ; i < cons.second.size() ; i++)
	{
		cons.second[i]->clearVisited() ;
	}
	
	std::vector<Point *> ret ;
	ret.insert(ret.end(), cons.first.begin(), cons.first.end()) ;
	
	return ret ;
}


std::vector<DelaunayTriangle *> DelaunayTree::conflicts(const Geometry *g) const
{
// 	Circle neighbourhood(1e-8, g->getCenter()) ; 
//magic value : basically, we don't want to have the center a vertex of the mesh.
	
// 	std::vector<DelaunayTriangle *> ret ;
	std::pair< std::vector<DelaunayTriangle *>,std::vector<DelaunayTreeItem *> > cons ;
	if(!tree.empty())
		this->tree[0]->conflicts(cons, g) ;
	
	for(size_t i = 0 ; i < plane.size() ; i++)
	{
		
		if(!plane[i]->visited)
		{
			std::pair< std::vector<DelaunayTriangle *>,std::vector<DelaunayTreeItem *> > temp ;
			plane[i]->conflicts(temp,g) ;
			
			cons.first.insert(cons.first.end(), temp.first.begin(),temp.first.end()) ;
			cons.second.insert(cons.second.end(), temp.second.begin(),temp.second.end()) ;			
			
			for(size_t j = 0 ; j < plane[i]->neighbour.size() ; j++)
			{
				std::pair< std::vector<DelaunayTriangle *>, std::vector<DelaunayTreeItem *> > temp ;
				plane[i]->getNeighbour(j)->conflicts(temp,g) ;
				
				cons.first.insert(cons.first.end(), temp.first.begin(),temp.first.end()) ;
				cons.second.insert(cons.second.end(), temp.second.begin(),temp.second.end()) ;
			}
			
		}
		
	}

	for(size_t i = 0 ; i < cons.second.size() ; i++)
		cons.second[i]->clearVisited() ;
	
	return cons.first ;
}

std::vector<DelaunayDemiPlane *> * DelaunayTree::getConvexHull()
{
	std::vector<DelaunayDemiPlane *> * ret = new std::vector<DelaunayDemiPlane *> ;
	
	for(size_t i = 0 ; i < plane.size() ; i++)
	{
		if(plane[i]->isAlive())
			ret->push_back(plane[i]) ;
	}
	
	return ret ;
}

std::vector<DelaunayTriangle *> DelaunayTree::getTriangles(bool buildNeighbourhood)
{

	std::vector<DelaunayTriangle *> ret ;
	
	for(size_t i = 0 ; i < tree.size() ; i++)
	{
		if(tree[i]->isTriangle && !tree[i]->isDeadTriangle)
		{
			ret.push_back((DelaunayTriangle *)(tree[i])) ;
		}
	}
	
	if(!neighbourhood && buildNeighbourhood)
	{
		std::cerr << "\r building neighbourhood... element 0/" << ret.size() << std::flush ;
		for( size_t i = 0 ; i < ret.size() ;i++)
		{
			ret[i]->neighbourhood.resize(0) ;
			ret[i]->clearVisited() ;
		}
		
		
		
		for( size_t i = 0 ; i < ret.size() ;i++)
		{
			if(i%100 == 0)
				std::cerr << "\r building neighbourhood... element " << i <<"/" << ret.size() << std::flush ;
			
			std::vector<DelaunayTriangle *> tocheck ;
			std::vector<DelaunayTriangle *> toclean ;
			for(size_t j = 0 ; j < ret[i]->neighbour.size() ; j++)
			{
				if(ret[i]->getNeighbour(j)->isTriangle && !ret[i]->getNeighbour(j)->visited)
				{
					tocheck.push_back(static_cast<DelaunayTriangle *>(ret[i]->getNeighbour(j)));
					ret[i]->getNeighbour(j)->visited = true ;
					toclean.push_back(*tocheck.rbegin()) ;
					ret[i]->addNeighbourhood(*tocheck.rbegin()) ;
				}
			}
			
			while(!tocheck.empty())
			{
				std::vector<DelaunayTriangle *> tocheck_temp ;
				for(size_t k = 0 ; k < tocheck.size() ; k++)
				{
					for(size_t j = 0 ; j < tocheck[k]->neighbour.size() ; j++)
					{
						if(
						    tocheck[k]->getNeighbour(j)->isTriangle 
						    && !tocheck[k]->getNeighbour(j)->visited
						    && tocheck[k]->getNeighbour(j) != ret[i]
						    && static_cast<DelaunayTriangle *>(tocheck[k]->getNeighbour(j))->isInNeighbourhood(ret[i])
						  )
						{
							tocheck_temp.push_back(static_cast<DelaunayTriangle *>(tocheck[k]->getNeighbour(j)));
							tocheck[k]->getNeighbour(j)->visited = true ;
							toclean.push_back(*tocheck_temp.rbegin()) ;
							ret[i]->addNeighbourhood(*tocheck_temp.rbegin()) ;
						}
					}
				}
				
				tocheck = tocheck_temp ;
			}
			for(size_t j = 0 ; j < toclean.size() ; j++)
			{
				toclean[j]->clearVisited() ;
			}
			
// 			for(size_t j = i+1 ; j <ret.size() ;j++)
// 			{
// 				if(ret[j]->isInNeighbourhood(ret[i]))
// 				{
// 					ret[j]->neighbourhood.push_back(ret[i]) ;
// 					ret[i]->neighbourhood.push_back(ret[j]) ;
// 				}
// 			}
		}
		
		std::cerr << " ...done" << std::endl ;
		neighbourhood = true ;
	}
	
	return ret ;
}

void DelaunayTree::print() const
{
	size_t alive = 0 ;
	std::cerr << "we have a total of " << tree.size() << " elements" << std::endl ;
	for(size_t i = 0 ; i < tree.size() ; i++)
	{
		if (tree[i]->isAlive())
		{
			alive++ ;
		}
	}
	std::cerr << "of which " << alive << " are alive" << std::endl ;
// 	#ifdef DEBUG
	for(size_t i = 0 ; i < tree.size() ; i++)
	{
		if (tree[i]->isAlive())
		{
			tree[i]->print() ;
			for(size_t j = 0 ; j< tree[i]->neighbour.size() ; j++)
			{
				std::cerr << "\t --> " ;
				tree[i]->getNeighbour(j)->print() ;
			}
		}
	}
// 	#endif
}

void DelaunayTriangle::refresh(const TriElement * father)
{
	this->TriElement::refresh(father) ;
	
// 	if(!cachedGPs)
// 		cachedGPs = new GaussPointArray(getSubTriangulatedGaussPoints()) ; 
	
// 	this->computeCenter() ;
}

std::vector<std::vector<Matrix> > & DelaunayTriangle::getElementaryMatrix() 
{
	if(!behaviourUpdated && !enrichmentUpdated)
	{
		return cachedElementaryMatrix ;
	}
	
	if(enrichmentUpdated)
	{
		for(size_t i = 0 ; i <  cachedElementaryMatrix.size() ; i++)
		{
			cachedElementaryMatrix[i].clear() ;
		}
		cachedElementaryMatrix.clear() ;
		cachedGaussPoints = getSubTriangulatedGaussPoints() ;
	}
	
	std::valarray<Matrix> Jinv ;
	std::vector<size_t > dofs = getDofIds() ;
	if(moved)
	{
		Jinv.resize(cachedGaussPoints.gaussPoints.size()) ;
		for(size_t i = 0 ; i < cachedGaussPoints.gaussPoints.size() ;  i++)
		{
			getInverseJacobianMatrix( cachedGaussPoints.gaussPoints[i].first, Jinv[i]) ;
		}
	}
	else
	{
		Matrix J ;
		getInverseJacobianMatrix(Point( 1./3.,1./3.), J ) ;
		Jinv.resize(cachedGaussPoints.gaussPoints.size(),J) ;
	}
	int size = getBehaviour()->getNumberOfDegreesOfFreedom() ;
	
	if(enrichmentUpdated)
	{
		for(size_t i = 0 ; i < dofs.size() ; i++)
		{
			std::vector< Matrix > v_j ;
			
			for(size_t j = 0 ; j < dofs.size() ; j++)
			{
				v_j.push_back(Matrix(size,size)) ;
			}
			
			cachedElementaryMatrix.push_back(v_j) ;
		}
	}
	
	VirtualMachine vm ;
	
	for(size_t i = 0 ; i < getShapeFunctions().size() ; i++)
	{
		 behaviour->apply(getShapeFunction(i), getShapeFunction(i),cachedGaussPoints, Jinv, cachedElementaryMatrix[i][i], &vm) ;
		
		for(size_t j = i+1 ; j < getShapeFunctions().size() ; j++)
		{
			 behaviour->apply(getShapeFunction(i), getShapeFunction(j),cachedGaussPoints, Jinv,cachedElementaryMatrix[i][j], &vm) ;
			 behaviour->apply(getShapeFunction(j), getShapeFunction(i),cachedGaussPoints, Jinv,cachedElementaryMatrix[j][i], &vm) ;
		}
		for(size_t j = 0 ; j < getEnrichmentFunctions().size() ; j++)
		{
			 behaviour->apply(getShapeFunction(i), getEnrichmentFunction(j),cachedGaussPoints, Jinv,cachedElementaryMatrix[i][j+getShapeFunctions().size()], &vm) ;
			 behaviour->apply(getEnrichmentFunction(j), getShapeFunction(i),cachedGaussPoints, Jinv,cachedElementaryMatrix[j+getShapeFunctions().size()][i], &vm) ;
		}
	}
	for(size_t i = 0 ; i < getEnrichmentFunctions().size() ; i++)
	{
		 behaviour->apply(getEnrichmentFunction(i), getEnrichmentFunction(i),cachedGaussPoints, Jinv,cachedElementaryMatrix[i+getShapeFunctions().size()][i+getShapeFunctions().size()], &vm) ;
		
		for(size_t j = 0 ; j < getEnrichmentFunctions().size() ; j++)
		{
			 behaviour->apply(getEnrichmentFunction(i), getEnrichmentFunction(j),cachedGaussPoints, Jinv,cachedElementaryMatrix[i+getShapeFunctions().size()][j+getShapeFunctions().size()], &vm) ;
			 behaviour->apply(getEnrichmentFunction(j), getEnrichmentFunction(i),cachedGaussPoints, Jinv,cachedElementaryMatrix[j+getShapeFunctions().size()][i+getShapeFunctions().size()], &vm) ;
		}
	}

	enrichmentUpdated = false ;
	behaviourUpdated = false ;
	if(behaviour->hasInducedForces())
		cachedForces.resize(0) ;
// 	exit(0) ;
	return cachedElementaryMatrix ;
}

std::vector<std::vector<Matrix> > DelaunayTriangle::getNonLinearElementaryMatrix() 
{
	std::vector<size_t > dofs = getDofIds() ;
	std::vector<std::vector<Matrix> > mother ;
	
	if(!this->getNonLinearBehaviour()->isActive())
	{
		for(size_t i = 0 ; i < dofs.size() ; i++)
		{
			std::vector< Matrix > v_j ;
			
			for(size_t j = 0 ; j < dofs.size() ; j++)
			{
				v_j.push_back(Matrix()) ;
			}
			
			mother.push_back(v_j) ;
		}
		
		return mother ;
	}
	
	std::valarray<Matrix> Jinv ;
	GaussPointArray gp = getSubTriangulatedGaussPoints() ;
	
	if(moved)
	{
		Jinv.resize(gp.gaussPoints.size(), Matrix()) ;
		for(size_t i = 0 ; i < gp.gaussPoints.size() ;  i++)
		{
			getInverseJacobianMatrix( gp.gaussPoints[i].first, Jinv[i] ) ;
		}
	}
	else
	{
		Matrix J ;
		getInverseJacobianMatrix(Point( 1./3.,1./3.) , J);
		Jinv.resize(gp.gaussPoints.size(), J) ;
	}

	int size = nonlinbehaviour->getNumberOfDegreesOfFreedom() ;
	
	
	for(size_t i = 0 ; i < dofs.size() ; i++)
	{
		std::vector< Matrix > v_j ;
		
		for(size_t j = 0 ; j < dofs.size() ; j++)
		{
			v_j.push_back(Matrix(size,size)) ;
		}
		
		mother.push_back(v_j) ;
	}
	
	VirtualMachine vm ;
	
	for(size_t i = 0 ; i < getShapeFunctions().size() ; i++)
	{
		 nonlinbehaviour->apply(getShapeFunction(i), getShapeFunction(i),gp, Jinv, mother[i][i], &vm) ;
		
		for(size_t j = i+1 ; j < getShapeFunctions().size() ; j++)
		{
			nonlinbehaviour->apply(getShapeFunction(i), getShapeFunction(j),gp, Jinv, mother[i][j], &vm) ;
			nonlinbehaviour->apply(getShapeFunction(j), getShapeFunction(i),gp, Jinv, mother[j][i], &vm) ;
		}
		for(size_t j = 0 ; j < getEnrichmentFunctions().size() ; j++)
		{
			nonlinbehaviour->apply(getShapeFunction(i), getEnrichmentFunction(j),gp, Jinv, mother[i][j+getShapeFunctions().size()], &vm) ;
			nonlinbehaviour->apply(getEnrichmentFunction(j), getShapeFunction(i),gp, Jinv, mother[j+getShapeFunctions().size()][i], &vm) ;
		}
	}
	for(size_t i = 0 ; i < getEnrichmentFunctions().size() ; i++)
	{
		nonlinbehaviour->apply(getEnrichmentFunction(i), getEnrichmentFunction(i),gp, Jinv, mother[i+getShapeFunctions().size()][i+getShapeFunctions().size()], &vm) ;
		
		for(size_t j = 0 ; j < getEnrichmentFunctions().size() ; j++)
		{
			nonlinbehaviour->apply(getEnrichmentFunction(i), getEnrichmentFunction(j),gp, Jinv, mother[i+getShapeFunctions().size()][j+getShapeFunctions().size()], &vm) ;
			nonlinbehaviour->apply(getEnrichmentFunction(j), getEnrichmentFunction(i),gp, Jinv, mother[j+getShapeFunctions().size()][i+getShapeFunctions().size()], &vm) ;
		}
	}
	
	return mother ;
}

std::vector<Point *> DelaunayTriangle::getIntegrationHints() const
{
	std::vector<Point *> to_add ;
	to_add.push_back(new Point(0,1)) ;
	to_add.push_back(new Point(0,0)) ;
	to_add.push_back(new Point(1,0)) ;
	

// 	std::vector<std::vector<Segment> > segments ;
// 	
// 	for(int i = 0 ; i <  (int)getEnrichmentFunctions().size() ; i++)
// 	{
// 		std::vector<Segment> temp ;
// 		for(int j = 0 ; j < (int)getEnrichmentFunction(i).getIntegrationHint().size()-1 ; j++)
// 		{
// 			temp.push_back(Segment(getEnrichmentFunction(i).getIntegrationHint(j), getEnrichmentFunction(i).getIntegrationHint(j+1))) ;
// 		}
// 		
// 		segments.push_back(temp) ;
// 	}
// 	for(size_t i = 0 ; i < segments.size() ; i++) // per shape function
// 	{
// 		for(size_t j = i+1 ; j < segments.size() ; j++) // per other shape function
// 		{
// 			for(size_t k = 0 ; k < segments[i].size() ; k++) //per segment in first SF
// 			{
// 				for(size_t l = 0 ; l < segments[j].size() ; l++) //per segment in second SF
// 				{
// 					if(segments[i][k].intersects(segments[j][l]))
// 					{
// 						bool go = true ;
// 						Point test = segments[i][k].intersection(segments[j][l]) ;
// 						for(int k = 0 ; k < to_add.size()  ; k++ )
// 						{
// 							if(squareDist2D(&test, to_add[k]) 
// 							   < POINT_TOLERANCE)
// 							{
// 								go = false ;
// 								break ;
// 							}
// 						}
// 						if(go)
// 						{
// 							to_add.push_back(new Point(test)) ;
// 							if(to_add.back()->x< 0)
// 								to_add.back()->x = 0 ;
// 							if(to_add.back()->y < 0)
// 								to_add.back()->y = 0 ;
// 						}
// 					}
// 				}
// 			}
// 		}
// 	}
	
	for(size_t i = 0 ; i <  getEnrichmentFunctions().size() ; i++)
	{
		
		for(size_t j = 0 ; j < getEnrichmentFunction(i).getIntegrationHint().size() ; j++)
		{
			bool go = true ;
			for(int k = 0 ; k < to_add.size()  ; k++ )
			{
				if(squareDist2D(getEnrichmentFunction(i).getIntegrationHint(j), *to_add[k]) 
					< POINT_TOLERANCE*POINT_TOLERANCE)
				{
					go = false ;
					break ;
				}
			}
			
			if(go /*&& f.in(getEnrichmentFunction(i).getIntegrationHint(j))*/
			  )
			{
				to_add.push_back(new Point(getEnrichmentFunction(i).getIntegrationHint(j))) ;
				if(to_add.back()->x < 0)
					to_add.back()->x = 0 ;
				if(to_add.back()->y < 0)
					to_add.back()->y = 0 ;
			}
		}
	}
	
	return to_add ;
}

GaussPointArray DelaunayTriangle::getSubTriangulatedGaussPoints() const
{
	if(!enrichmentUpdated)
		return cachedGaussPoints ;

	GaussPointArray gp = getGaussPoints() ; 
	size_t numberOfRefinements = 8;
	
	double tol = 1e-6 ;
	double position_tol = 4.*POINT_TOLERANCE ;
	double infinity = .15 ;
	VirtualMachine vm ;
	if(getEnrichmentFunctions().size() > 0)
	{
		std::vector<std::pair<Point, double> > gp_alternative ;
		VirtualMachine vm ;
		std::vector<Point *> to_add = getIntegrationHints();
		std::vector<Point *> pointsToCleanup = to_add;
		std::vector<DelaunayTriangle *> triangleToCleanup;
		std::vector<DelaunayTriangle *> tri ;
		std::vector<bool> pass ;
		int passNum = 0;
		double J = jacobianAtPoint(Point(1./3., 1./3.)) ;
		double lastError = 10 ;
		size_t maxGradientIndex = 0 ;
		std::vector<double> grads(getEnrichmentFunctions().size(), 0.) ;
		
		double ndivs = 30 ;
		for(double k = 0  ; k < ndivs ; k++)
		{
			for(double l = 0  ; l < ndivs ; l++)
			{
				if( k+l < ndivs )
					gp_alternative.push_back(std::make_pair(Point(k/ndivs, l/ndivs), 1.)) ;
			}
		}
		double a = area() ; 
		for(size_t i =0 ; i < gp_alternative.size() ; i++)
		{
			gp_alternative[i].second = 0.25*a/gp_alternative.size() ;
		}
		
		if(gp.gaussPoints.size() < gp_alternative.size())
		{
			
			gp.gaussPoints.resize(gp_alternative.size()) ;
			std::copy(gp_alternative.begin(), gp_alternative.end(), &gp.gaussPoints[0]);
			gp.id = -1 ;
		}
		return gp ;

		DelaunayTree * dt = new DelaunayTree(to_add[0], to_add[1], to_add[2]) ;
		TriElement f(LINEAR) ;
		if(to_add.size() > 4)
			std::random_shuffle(to_add.begin()+3, to_add.end()) ;

		for(size_t i = 3 ; i < to_add.size() ; i++)
		{
			dt->insert(to_add[i]) ;
		}

		tri = dt->getTriangles(false) ;
		// original set of triangles
		
// 		
		pass.resize(tri.size(), false) ;
		
		for(size_t j = 0 ; j < tri.size() ; j++)
		{
			for(size_t k = 0 ; k <  getEnrichmentFunctions().size() ; k++)
			{
				
				double dx = vm.deval(getEnrichmentFunction(k),XI, tri[j]->getCenter()) ;
				double dy = vm.deval(getEnrichmentFunction(k),ETA, tri[j]->getCenter()) ;
				grads[k] = std::max((dx*dx+dy*dy), grads[k]) ;
			}
		}		
		
		for(size_t i = 0 ; i < grads.size() ; i++)
			if(grads[i] > grads[maxGradientIndex])
				maxGradientIndex = i ;

		std::vector<double> triIntegral ;
		double integral = 0 ;
		for(size_t j = 0 ; j < tri.size() ; j++)
		{
			double dx = vm.deval(getEnrichmentFunction(maxGradientIndex),XI, tri[j]->getCenter()) ;
			double dy = vm.deval(getEnrichmentFunction(maxGradientIndex),ETA, tri[j]->getCenter()) ;
			double t_integral = (dx*dx+dy*dy)*tri[j]->area()*J*.5 ;

			triIntegral.push_back(t_integral) ;
			integral += t_integral ;
			
		}

// 		vm.print(getEnrichmentFunction(maxGradientIndex)) ;// 		if(to_add.size() < 5)
// 			numberOfRefinements = 0 ;
		for(size_t i = 0 ; i < numberOfRefinements ; i++)
		{

			double newIntegral = 0;

			for(size_t j = 0 ; j < tri.size() ; j++)
			{
				if(!pass[j] )
				{
					Function x = XTransform(tri[j]->getBoundingPoints(), f.getShapeFunctions()) ;
					Function y = YTransform(tri[j]->getBoundingPoints(), f.getShapeFunctions()) ;
					
					double error = 0 ;

					tri[j]->setOrder(QUADTREE_REFINED) ;
					GaussPointArray gpquad = tri[j]->getGaussPoints() ;
					tri[j]->setOrder(LINEAR) ;
// 					GaussPointArray gplin = tri[j]->getGaussPoints() ;
					
					
					for(size_t m = 0 ; m < gpquad.gaussPoints.size() ; m++)
					{
						gpquad.gaussPoints[m].first.set(vm.eval(x, gpquad.gaussPoints[m].first), vm.eval(y, gpquad.gaussPoints[m].first)) ;
					}
// 					for(size_t m = 0 ; m < gplin.gaussPoints.size() ; m++)
// 					{
// 						gplin.gaussPoints[m].first.set(vm.eval(x, gplin.gaussPoints[m].first), vm.eval(y, gplin.gaussPoints[m].first)) ;
// 					}
					
					double newTriIntegral = 0 ;

					for(size_t m = 0 ; m < gpquad.gaussPoints.size() ; m++)
					{
						double dx = vm.deval(getEnrichmentFunction(maxGradientIndex),XI, gpquad.gaussPoints[m].first) ;
						double dy = vm.deval(getEnrichmentFunction(maxGradientIndex),ETA, gpquad.gaussPoints[m].first) ;
						newTriIntegral += (dx*dx+dy*dy)*gpquad.gaussPoints[m].second*J ;
					}
					
					newIntegral += newTriIntegral ;
					error = std::abs((triIntegral[j]-newTriIntegral)) ;;
					if(error < tol )
					{
						pass[j] = true ;
					}
				}
			}

			
			std::vector<DelaunayTriangle *> newTris ;
			std::vector<bool> newPass ;
			triIntegral.clear() ;
			for(size_t j = 0 ; j < tri.size() ; j++)
			{
			
				if(!pass[j] && tri[j]->area() > sqrt(std::numeric_limits<double>::epsilon()))
				{

					std::pair<std::vector<DelaunayTriangle *>, std::vector<Point *> > q =
						quad(tri[j]) ;
					newTris.insert(newTris.end(),q.first.begin(), q.first.end()) ;
			
					newPass.push_back(false) ;
					newPass.push_back(false) ;
					newPass.push_back(false) ;
					newPass.push_back(false) ;
					pointsToCleanup.insert(pointsToCleanup.end(),
							q.second.begin(), 
							q.second.end()) ;
					triangleToCleanup.insert(triangleToCleanup.end(), 
								q.first.begin(), 
								q.first.end()) ;

				}
				else if(pass[j])
				{
					newTris.push_back(tri[j]) ;
					newPass.push_back(true) ;
					passNum++ ;
				}
				else
				{
					newTris.push_back(tri[j]) ;
					newPass.push_back(false) ;
				}
			}

			tri = newTris ;
			pass = newPass ;
			if(std::abs(integral-newIntegral) < 1e-10)
				break ;
			integral = newIntegral ;
			
			for(size_t j = 0 ; j < tri.size() ; j++)
			{
				double dx = vm.deval(getEnrichmentFunction(maxGradientIndex),XI, tri[j]->getCenter()) ;
				double dy = vm.deval(getEnrichmentFunction(maxGradientIndex),ETA, tri[j]->getCenter()) ;
				double t_integral = (dx*dx+dy*dy)*tri[j]->area()*J*.5 ;

				triIntegral.push_back(t_integral) ;
				
			}
		}
// 		
		for(size_t i = 0 ; i < tri.size() ; i++)
		{

			Function x = XTransform(tri[i]->getBoundingPoints(), f.getShapeFunctions()) ;
			Function y = YTransform(tri[i]->getBoundingPoints(), f.getShapeFunctions()) ;
			
			tri[i]->setOrder(CUBIC) ;
			GaussPointArray gp_temp = tri[i]->getGaussPoints() ;
			tri[i]->setOrder(LINEAR) ;
			
			for(size_t j = 0 ; j < gp_temp.gaussPoints.size() ; j++)
			{

				gp_temp.gaussPoints[j].first.set(vm.eval(x, gp_temp.gaussPoints[j].first), vm.eval(y, gp_temp.gaussPoints[j].first)) ;
				if(moved)
					gp_temp.gaussPoints[j].second *= this->jacobianAtPoint(gp_temp.gaussPoints[j].first) ;
				else
					gp_temp.gaussPoints[j].second *= J; 
				
				gp_alternative.push_back(gp_temp.gaussPoints[j]) ;
			}
		}

		delete dt ;
		if(numberOfRefinements)
		{

			std::valarray<Point *> nularray(0) ;
	
			for(size_t i = 0 ; i < triangleToCleanup.size() ; i++)
			{
				triangleToCleanup[i]->setBoundingPoints(nularray) ;
				delete triangleToCleanup[i];
			}
		
		}
		for(size_t i = 0 ; i < pointsToCleanup.size() ; i++)
			delete pointsToCleanup[i] ;
		
		
		if(gp.gaussPoints.size() < gp_alternative.size())
		{
			
			gp.gaussPoints.resize(gp_alternative.size()) ;
			std::copy(gp_alternative.begin(), gp_alternative.end(), &gp.gaussPoints[0]);
			gp.id = -1 ;
		}
	}
	return gp ;
}



Vector DelaunayTriangle::getNonLinearForces()
{
	std::vector<size_t> dofs = getDofIds() ;
	Vector forces(dofs.size()*2) ;
	
	if(!this->getNonLinearBehaviour()->isActive())
	{
		forces = 0 ;
		return forces ;
	}
	
	std::valarray<Matrix> Jinv ;
	
	GaussPointArray gp  = getSubTriangulatedGaussPoints() ;
	
	if(moved)
	{
		Jinv.resize(gp.gaussPoints.size(), Matrix()) ;
		for(size_t i = 0 ; i < gp.gaussPoints.size() ;  i++)
		{
			getInverseJacobianMatrix( gp.gaussPoints[i].first, Jinv[i]) ;
		}
	}
	else
	{
		Matrix J ;
		getInverseJacobianMatrix(Point( 1./3.,1./3.) , J) ;
		Jinv.resize(gp.gaussPoints.size(), J) ;
	}
	
	size_t numdof = getBehaviour()->getNumberOfDegreesOfFreedom() ;
	Vector f(0., numdof) ;
	for(size_t i = 0 ; i < getShapeFunctions().size() ; i++)
	{

		behaviour->getForces(this->getState(), getShapeFunction(i) ,gp, Jinv, f) ;
			
		for(size_t j = 0 ; j < numdof ; j++)
			forces[i*numdof+j]+=f[j];

	}
		
	for(size_t i = 0 ; i < getEnrichmentFunctions().size() ; i++)
	{
		behaviour->getForces(this->getState(), getEnrichmentFunction(i) ,gp, Jinv, f) ;
		
		for(size_t j = 0 ; j < numdof ; j++)
			forces[(i+getShapeFunctions().size())*numdof+j]+=f[j];
	}

	
	return forces ;
}

Vector DelaunayTriangle::getForces()
{
	if(!behaviourUpdated && !enrichmentUpdated && cachedForces.size() > 0)
	{
		return cachedForces ;
	}
	
	if(enrichmentUpdated)
		cachedGaussPoints = getSubTriangulatedGaussPoints() ;

	size_t numdof = getBehaviour()->getNumberOfDegreesOfFreedom() ;
	std::vector<size_t> dofs = getDofIds() ;
	cachedForces.resize(dofs.size()*numdof, double(0)) ;
	std::valarray<Matrix> Jinv ;

	if(moved)
	{
		Jinv.resize( cachedGaussPoints.gaussPoints.size(), Matrix()) ;
		for(size_t i = 0 ; i < cachedGaussPoints.gaussPoints.size() ;  i++)
		{
			getInverseJacobianMatrix( cachedGaussPoints.gaussPoints[i].first, Jinv[i] ) ;
		}
	}
	else
	{
		Matrix J ;
		getInverseJacobianMatrix(Point( 1./3.,1./3.), J ) ;

		Jinv.resize(cachedGaussPoints.gaussPoints.size(), J) ;
	}
// 	
	
	Vector f(0., numdof) ;
	for(size_t i = 0 ; i < getShapeFunctions().size() ; i++)
	{
		behaviour->getForces(this->getState(), getShapeFunction(i),cachedGaussPoints, Jinv, f) ;
		for(size_t j = 0 ; j < numdof ; j++)
			cachedForces[i*numdof+j]+=f[j];

	}
		
	for(size_t i = 0 ; i < getEnrichmentFunctions().size() ; i++)
	{
		behaviour->getForces(this->getState(), getEnrichmentFunction(i) ,cachedGaussPoints, Jinv, f) ;
		for(size_t j = 0 ; j < numdof ; j++)
			cachedForces[(i+getShapeFunctions().size())*numdof+j]+=f[j];
	}
	enrichmentUpdated = false ;
	behaviourUpdated = false ;
	return cachedForces ;
}

DelaunayDeadTriangle::~DelaunayDeadTriangle() { }

std::pair<std::vector<DelaunayTriangle *>, std::vector<Point *> > Mu::quad(const DelaunayTriangle * t)
{
	std::vector<DelaunayTriangle* > tris ;
	std::vector<Point *> points ;

	points.push_back(new Point(*t->first + (*t->second - *t->first)*.5) ) ;
	points.push_back(new Point(*t->first + (*t->third - *t->first)*.5 )) ;
	points.push_back(new Point(*t->second + (*t->third - *t->second)*.5 )) ;
	
	tris.push_back(new DelaunayTriangle(NULL,NULL, points[0], points[1], t->first, NULL)) ;
	tris.push_back(new DelaunayTriangle(NULL,NULL,points[0], points[2], t->second, NULL)) ;
	tris.push_back(new DelaunayTriangle(NULL,NULL,points[1], points[2], t->third, NULL)) ;
	tris.push_back(new DelaunayTriangle(NULL,NULL,points[0], points[1], points[2], NULL)) ;
	return std::make_pair(tris, points) ;
}
