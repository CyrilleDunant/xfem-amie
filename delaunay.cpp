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

#define DEBUG
#undef DEBUG

using namespace Mu ;

DelaunayTreeItem::DelaunayTreeItem( DelaunayTreeItem * father,  const Point * c) 
{
	this->stepfather = NULL ;
	this->father = father;
	this->m_c  = c ;
	this->dead = false ;
	this->visited = false ;
	
	this->isPlane =false ;
	this->isTriangle =false ;
	this->first = NULL ;
	this->second = NULL ;
	this->third = NULL ;
}
	
bool DelaunayTreeItem::isConflicting(const Geometry * g) const
{
	return inCircumCircle(g->getCenter()) || g->in(*first) || g->in(*second) || ( isTriangle && g->in(*third)) ;
}


DelaunayTreeItem::~DelaunayTreeItem()
{
	this->neighbour.clear() ;
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
	
	for(size_t i = 0 ; i < treeitem.size() ;i++)
	{
		items.insert(items.end() , treeitem[i]->son.begin() , treeitem[i]->son.end()) ;
		items.insert(items.end() , treeitem[i]->stepson.begin() , treeitem[i]->stepson.end()) ;
		items.insert(items.end() , treeitem[i]->neighbour.begin() , treeitem[i]->neighbour.end()) ;
	}
	
	for( size_t i = 0 ; i < items.size() ;i++)
	{
		for(size_t j = i+1 ; j <items.size() ;j++)
		{
			makeNeighbours(items[i], items[j]) ;
		}
	}

}

void DelaunayTree::addSharedNodes(size_t nodes_per_side, size_t time_planes, double timestep, const TriElement * father)
{
	std::vector<DelaunayTriangle *> tri = getTriangles() ;
	for(size_t i = 0 ; i < tri.size() ; i++)
	{
		
		tri[i]->visited = true ;
			
		size_t nodes_per_plane = nodes_per_side*3+3 ;
		
		std::valarray<Point *> newPoints(nodes_per_plane*time_planes) ;
		std::valarray<bool> done(false, nodes_per_plane*time_planes) ;
		
		for(size_t plane = 0 ; plane < time_planes ; plane++)
		{
			for(size_t side = 0 ; side < 3 ; side++)
			{
				Point a(tri[i]->getBoundingPoint(side)) ;
				Point b(tri[i]->getBoundingPoint((side+1)%3)) ;
				
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
					
					for(size_t j = 0 ; j< tri[i]->getBoundingPoints().size() ; j++)
					{
						if(tri[i]->getBoundingPoint(j) == proto)
						{
							foundPoint = &tri[i]->getBoundingPoint(j) ;
							break ;
						}
					}
					
					if(!foundPoint)
					{
						for(size_t j = 0 ; j < tri[i]->neighbour.size() ; j++)
						{
							if(tri[i]->neighbour[j]->isTriangle && tri[i]->neighbour[j]->visited)
							{
								DelaunayTriangle * n = static_cast<DelaunayTriangle *>(tri[i]->neighbour[j]) ;
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
							newPoints[nodes_per_plane*plane+side*(nodes_per_side+1)+node]  = new Point(proto) ;
							newPoints[nodes_per_plane*plane+side*(nodes_per_side+1)+node]->id = global_counter++ ;
						}
						
						done[nodes_per_plane*plane+side*(nodes_per_side+1)+node] = true ;
					}
				}
			}
		}
		
		tri[i]->setBoundingPoints(newPoints) ;
	}
			
	
	for(size_t i = 0 ; i < tri.size() ; i++)
	{
		tri[i]->clearVisited() ;
	}

}

void DelaunayTree::refresh(TriElement *father, bool compile)
{
	if(compile)
	{
		for(size_t k = 0 ; k < father->getShapeFunctions().size() ; k++)
		{
			father->getShapeFunction(k).compile() ;
		}
	}
	
	for(size_t i = 0 ; i < this->tree.size() ; i++)
	{
		if(this->tree[i]->isAlive() && this->tree[i]->isTriangle)
		{
			((DelaunayTriangle *)(this->tree[i]))->refresh(father) ;
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
	
	
	if(!inCircumCircle(g->getCenter()) && (!g->in(*first) && !g->in(*second) && ( !g->in(*third) && isTriangle ) || ( g->in(*third) && isPlane )) )
	{
		return ;
	}
	
	for (size_t i  = 0 ;  i < stepson.size() ; i++)
	{
		if( !stepson[i]->visited && !(!stepson[i]->inCircumCircle(g->getCenter()) && (!g->in(*stepson[i]->first) && !g->in(*stepson[i]->second) && ( !g->in(*stepson[i]->third) && stepson[i]->isTriangle ) || ( g->in(*stepson[i]->third) && stepson[i]->isPlane ))))
		{
			std::pair<std::vector<DelaunayTriangle *>, std::vector<DelaunayTreeItem *> > temp ;
			stepson[i]->conflicts(temp, g) ;
			ret.first.insert(ret.first.end(), temp.first.begin(), temp.first.end()) ;
			ret.second.insert(ret.second.end(), temp.second.begin(), temp.second.end()) ;
		}
	}
	
	for (size_t i  = 0 ;  i < son.size() ; i++)
	{
		if( !son[i]->visited && !(!son[i]->inCircumCircle(g->getCenter()) && (!g->in(*son[i]->first) && !g->in(*son[i]->second) && ( !g->in(*son[i]->third) && son[i]->isTriangle ) || ( g->in(*son[i]->third) && son[i]->isPlane ))))
		{
			std::pair<std::vector<DelaunayTriangle *>, std::vector<DelaunayTreeItem *> > temp ;
			son[i]->conflicts(temp,g) ;
			ret.first.insert(ret.first.end(), temp.first.begin(), temp.first.end()) ;
			ret.second.insert(ret.second.end(), temp.second.begin(), temp.second.end()) ;
		}
	}
	
	if(isAlive() && isTriangle)
	{
		ret.first.push_back((DelaunayTriangle *)(this)) ;
	}
	
	for (size_t i  = 0 ;  i < neighbour.size() ; i++)
	{
		if( !neighbour[i]->visited && !(!neighbour[i]->inCircumCircle(g->getCenter()) && (!g->in(*neighbour[i]->first) && !g->in(*neighbour[i]->second) && ( !g->in(*neighbour[i]->third) && neighbour[i]->isTriangle ) || ( g->in(*neighbour[i]->third) && neighbour[i]->isPlane ))))
		{
			std::pair<std::vector<DelaunayTriangle *>, std::vector<DelaunayTreeItem *> > temp  ;
			neighbour[i]->conflicts(temp, g) ;
			ret.first.insert(ret.first.end(), temp.first.begin(), temp.first.end()) ;
			ret.second.insert(ret.second.end(), temp.second.begin(), temp.second.end()) ;
			
		}
	}
	return ;
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
		if( !stepson[i]->visited) 
		{
			temp.first.clear() ;
			temp.second.clear() ;
			stepson[i]->conflicts(temp,s) ;
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
		if( !son[i]->visited )
		{
			temp.first.clear() ;
			temp.second.clear() ;
			son[i]->conflicts(temp,s) ;
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
		if( !neighbour[i]->visited )
		{
			temp.first.clear() ;
			temp.second.clear() ;
			neighbour[i]->conflicts(temp,s) ;
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
	for (size_t i  = 0 ;  i < stepson.size() ; i++)
	{
		bool limit = false ;
		if(!stepson[i]->visited && stepson[i]->isTriangle)
		{
			DelaunayTriangle * t = static_cast<DelaunayTriangle *>(stepson[i]) ;
			limit = std::abs(dist(t->getCircumCenter(),*p)-t->getRadius()) < 2.*POINT_TOLERANCE ;
		}
		
		if( (!stepson[i]->visited && stepson[i]->inCircumCircle(*p)) || limit) 
		{
			std::pair<std::vector<DelaunayTreeItem *>, std::vector<DelaunayTreeItem *> > temp ;
			stepson[i]->conflicts(temp,p) ;
			ret.first.insert(ret.first.end(), temp.first.begin(), temp.first.end()) ;
			ret.second.insert(ret.second.end(), temp.second.begin(), temp.second.end()) ;
		}
	}
// 	std::cerr << "sons" << std::endl ;
	for (size_t i  = 0 ;  i < son.size() ; i++)
	{
		bool limit = false ;
		if(!son[i]->visited && son[i]->isTriangle)
		{
			DelaunayTriangle * t = static_cast<DelaunayTriangle *>(son[i]) ;
			limit = std::abs(dist(t->getCircumCenter(),*p)-t->getRadius()) < 2.*POINT_TOLERANCE ;
		}
		
		if( (!son[i]->visited && son[i]->inCircumCircle(*p)) || limit)
		{
			std::pair<std::vector<DelaunayTreeItem *>, std::vector<DelaunayTreeItem *> > temp ;
			son[i]->conflicts(temp,p) ;
			ret.first.insert(ret.first.end(), temp.first.begin(), temp.first.end()) ;
			ret.second.insert(ret.second.end(), temp.second.begin(), temp.second.end()) ;
		}
	}
	
	if(isAlive())
	{
		ret.first.push_back(this) ;
	}
// 	std::cerr << "neighbour" << std::endl ;
	for (size_t i  = 0 ;  i < neighbour.size() ; i++)
	{
		
		bool limit = false ;
		if(!neighbour[i]->visited && neighbour[i]->isTriangle)
		{
			DelaunayTriangle * t = static_cast<DelaunayTriangle *>(neighbour[i]) ;
			limit = std::abs(dist(t->getCircumCenter(),*p)-t->getRadius()) < 2.*POINT_TOLERANCE ;
		}
		
		if( (!neighbour[i]->visited && neighbour[i]->inCircumCircle(*p)) || limit)
		{
			std::pair<std::vector<DelaunayTreeItem *>, std::vector<DelaunayTreeItem *> > temp  ;
			neighbour[i]->conflicts(temp,p) ;
			ret.first.insert(ret.first.end(), temp.first.begin(), temp.first.end()) ;
			ret.second.insert(ret.second.end(), temp.second.begin(), temp.second.end()) ;
		}
	}
	
	return ;
}

void DelaunayTreeItem::removeNeighbour(DelaunayTreeItem * t)
{
	for(std::vector<DelaunayTreeItem *>::iterator i = neighbour.begin() ; i!=neighbour.end() ;++i)
	{
		if(*i == t)
		{
			neighbour.erase(i) ;
			return ;
		}
	}
}

void DelaunayTreeItem::removeDeadNeighbour(DelaunayTreeItem * t)
{
	if(std::find(deadneighbour.begin(), deadneighbour.end(),t) !=  neighbour.end())
	{
		deadneighbour.erase(std::find(deadneighbour.begin(), deadneighbour.end(), t)) ;
	}
}
	
void DelaunayTreeItem::addNeighbour(DelaunayTreeItem * t)
{
	assert(t != this) ;
	if(std::find(neighbour.begin(), neighbour.end(), t) != neighbour.end())
	{
		return ;
	}
	
	if(t->isAlive())
	{
		neighbour.push_back(t) ;
	}
}
	

void DelaunayTreeItem::addDeadNeighbour(DelaunayTreeItem * t)
{
	assert(t != this) ;
	if(std::find(deadneighbour.begin(), deadneighbour.end(), t) != deadneighbour.end())
	{
		return ;
	}
	
	deadneighbour.push_back(t) ;
}

DelaunayTreeItem * DelaunayTreeItem::getNeighbour(size_t i)
{
	return this->neighbour[i] ;
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
		this->neighbour[i]->removeNeighbour(this) ;
	}
}
	
bool DelaunayTreeItem::isAlive() const
{
	return !this->dead ;
}
	
void DelaunayTreeItem::addStepson(DelaunayTreeItem * s)
{
	if(s == this)
		return ;
	stepson.push_back(s) ;
	s->setStepfather(this) ;
	addNeighbour(s) ;
}

void DelaunayTreeItem::addSon(DelaunayTreeItem * s)
{
	son.push_back(s) ;
}
	
void DelaunayTreeItem::removeSon(DelaunayTreeItem * t)
{
	if(std::find(son.begin(), son.end(),t) !=  son.end())
	{
		son.erase(std::find(son.begin(), son.end(), t)) ;
	}
}

void DelaunayTreeItem::removeStepson(DelaunayTreeItem * t)
{
	if(std::find(stepson.begin(), stepson.end(),t) !=  stepson.end())
	{
		stepson.erase(std::find(stepson.begin(), stepson.end(), t)) ;
	}
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

DelaunayTriangle::DelaunayTriangle(DelaunayTreeItem * father,  Point *p0,  Point *p1, Point *p2,  Point * c) : TriElement(p0, p1, p2), DelaunayTreeItem(father, c)
{
	first = &getBoundingPoint(0) ;
	second = &getBoundingPoint(1) ;
	third = &getBoundingPoint(2) ;
	
	assert(in(this->getCenter())) ;
	
	neighbour.clear() ;
	isPlane = false ;
	isTriangle = true ;
	assert(first->id > -1) ;
	assert(second->id > -1) ;
	assert(third->id > -1) ;
	cachedGPs = NULL ;
}

DelaunayTriangle::DelaunayTriangle() : DelaunayTreeItem(NULL, NULL)
{
	first = &getBoundingPoint(0) ;
	second = &getBoundingPoint(1) ;
	third = &getBoundingPoint(2) ;
	
	assert(in(this->getCenter())) ;
	
	neighbour.clear() ;
	isPlane = false ;
	isTriangle = true ;
	cachedGPs = NULL ;
}
	
DelaunayTriangle::~DelaunayTriangle()
{
	delete cachedGPs ;
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
	return (p->id == first->id || p->id == second->id || p->id == third->id) ;
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

std::pair<Point*, Point*> DelaunayTriangle::commonEdge(const DelaunayTreeItem * t) 
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
		if(isAligned((*t->first), (*first), (*second)) && isAligned((*t->second), (*first), (*second)))
			return std::pair< Point*,  Point*>(first, second ) ;
		if(isAligned((*t->first), (*third), (*second)) && isAligned((*t->second), (*third), (*second)))
			return std::pair< Point*,  Point*>(third, second ) ;
		if(isAligned((*t->first), (*third), (*first)) && isAligned((*t->second), (*third), (*first)))
			return std::pair< Point*,  Point*>(first, third ) ;
	}
	
	return std::pair< Point*,  Point*>(NULL, NULL) ;
}

std::pair< Point*,  Point*> DelaunayDemiPlane::commonEdge(const DelaunayTreeItem * t) 
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
					if(p->neighbour[i]->isNeighbour(neighbour[j]))
						makeNeighbours(neighbour[j], p->neighbour[i]) ;
				}
			}
			p->kill(first) ;
			kill(first) ;
			return ;
		}
		if(isAligned((*first), (*second), (*p->first)) && isAligned((*first), (*second), (*p->second)) &&
		   isAligned((*first), (*p->first) ,(*p->second)) && isAligned((*second), (*p->first) ,(*p->second)))
		{

			for(size_t i = 0 ; i <  p->neighbour.size() ; i++)
			{
				makeNeighbours(this, p->neighbour[i]) ;
			}
			for(size_t i = 0 ; i <  p->son.size() ; i++)
			{
				p->son[i]->father = this ;
				son.push_back(p->son[i]) ;
			}
			for(size_t i = 0 ; i <  p->stepson.size() ; i++)
			{
				p->stepson[i]->stepfather = this ;
				stepson.push_back(p->stepson[i]) ;
				makeNeighbours(this, p->stepson[i]) ;
			}
			p->father->son.push_back(this) ;
			p->son.push_back(this) ;
			p->kill(first) ;
		}
	}
}


std::pair< Point*,  Point*> DelaunayTriangle::nearestEdge(const Point p)
{
	std::map<double, Point> cen ;
	Point c0(((*first) + (*second))/2.) ;
	cen[squareDist(c0, p)] = c0 ;
	Point c1(((*third) + (*second))/2.) ;
	cen[squareDist(c1, p)] = c1 ;
	Point c2(((*third) + (*first))/2.) ;
	cen[squareDist(c2, p)] = c2 ;
	
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
	return  this->Triangle::inCircumCircle(p) ;

}
	
void DelaunayTriangle::insert(std::vector<DelaunayTreeItem *> &ret, Point *p,  Star* s)
{

	for (size_t i = 0 ; i < neighbour.size() ; i++)
	{
		if(this->numberOfCommonVertices(neighbour[i]) == 2)
		{
			std::pair< Point*,  Point*> pp = this->commonEdge(neighbour[i]) ;
			if (!neighbour[i]->inCircumCircle(*p))
			{

				if(!isAligned((*p), (*pp.first), (*pp.second) ))
				{
					DelaunayTriangle *ss = new DelaunayTriangle(this, p, pp.first, pp.second, p) ;
					son.push_back(ss) ;

					neighbour[i]->addStepson(ss) ;
					ret.push_back(ss) ;
					
				}

			}
		}

	}
	this->kill(p) ;
	s->updateNeighbourhood() ;

}
	

bool DelaunayTriangle::isNeighbour(DelaunayTreeItem * t)  
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
		std::cerr << "(" << getBoundingPoint(i).x << ";" << getBoundingPoint(i).y <<  ";" << getBoundingPoint(i).t <<") " ;
	}
	std::cerr <<  ":: "<< isAlive() << std::endl ;
}


// void DelaunayTriangle::displace(Vector * eps)
// {
// 	(*eps)[first->id*2]+=(*eps)[first->id*2] ;
// 	(*eps)[first->id*2+1]+=(*eps)[first->id*2+1] ;
// }

DelaunayDemiPlane::DelaunayDemiPlane(DelaunayTreeItem * father,  Point  * _begin,  Point  * _end,  Point  * p,  Point * c) : DelaunayTreeItem(father, c)
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

std::pair< Point*,  Point*> DelaunayDemiPlane::nearestEdge(const Point p) 
{
	return std::pair< Point*,  Point*>(first, second) ;
}
	
bool DelaunayDemiPlane::inCircumCircle(const Point &p) const
{
	return ((vector.x*(p.y - first->y) - vector.y*(p.x - first->x)) * direction < -10*std::numeric_limits<double>::epsilon()) ;
}
	
bool DelaunayDemiPlane::isVertex(const Point *p) const
{
	return ( (*p) == (*first) || (*p) == (*second) ) || isAligned((*p), (*first), (*second)) ;
}
	

bool  DelaunayDemiPlane::isNeighbour(DelaunayTreeItem * t)
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

	std::vector<Point*> lims ;
	
	for(size_t i = 0 ; i < neighbour.size() ; i++)
	{
		if(numberOfCommonVertices(neighbour[i]) == 2)
		{
			std::pair< Point*,  Point*> pp = neighbour[i]->commonEdge(this) ;
			if(std::find(lims.begin(), lims.end(), pp.first) != lims.end())
				lims.erase(std::find(lims.begin(), lims.end(), pp.first)) ;
			else
				lims.push_back(pp.first) ;
			
			if(std::find(lims.begin(), lims.end(), pp.second) != lims.end())
				lims.erase(std::find(lims.begin(), lims.end(), pp.second)) ;
			else
				lims.push_back(pp.second) ;
			
			if (!neighbour[i]->inCircumCircle(*p))
			{

				assert(neighbour[i]->isNeighbour(this) );
				
				if(!isAligned((*p), (*first), (*second)))
				{

					DelaunayTriangle *ss = new DelaunayTriangle(this, p, pp.first, pp.second , p) ;
					ss->visited = true ;

					son.push_back(ss) ;
					neighbour[i]->addStepson(ss) ;
		
					ret.push_back(ss) ;
				}
			}
		}
	}
	this->kill(p) ;

	for(size_t i = 0 ; i < deadneighbour.size() ; i++)
	{
		std::pair< Point*, Point*> pp = deadneighbour[i]->commonEdge(this) ;
		if(std::find(lims.begin(), lims.end(), pp.first) != lims.end())
			lims.erase(std::find(lims.begin(), lims.end(), pp.first)) ;
		else
			lims.push_back(pp.first) ;
		
		if(std::find(lims.begin(), lims.end(), pp.second) != lims.end())
			lims.erase(std::find(lims.begin(), lims.end(), pp.second)) ;
		else
			lims.push_back(pp.second) ;
	}
	if(lims.empty())
	{
		lims.push_back(first) ;
		lims.push_back(second) ;
	}
	
	assert(lims.size() < 3) ;
	
	assert(!isAligned((*lims[0]), (*p), (*lims[1]))) ;
	assert((*lims[0]) != (*lims[1])) ;
	
	DelaunayDemiPlane *p0 = new DelaunayDemiPlane(this, lims[0], p, lims[1], p) ;

	DelaunayDemiPlane *p1 = new DelaunayDemiPlane(this, lims[1], p, lims[0], p) ;

	son.push_back(p0) ;
	son.push_back(p1) ;
	
	ret.push_back(p0) ;
	ret.push_back(p1) ;
	
	for(size_t i = 0 ; i < son.size()-2 ; i++)
	{
		if(son[i]->isNeighbour(p0))
			son[i]->addStepson(p0) ;
		
		if(son[i]->isNeighbour(p1))
			son[i]->addStepson(p1) ;
		
		for(size_t j = i ; j < son.size()-2 ; j++)
			if(son[i]->isNeighbour(son[j]))
				makeNeighbours(son[i], son[j]) ;
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

DelaunayRoot::DelaunayRoot(Point * p0, Point * p1, Point * p2) : DelaunayTreeItem(NULL, NULL)
{
	isPlane = false ;
	isTriangle =false ;
	this->father = NULL ;
	DelaunayTriangle *t = new DelaunayTriangle(this, p0, p1, p2, NULL) ;
	DelaunayDemiPlane * pl0 = new DelaunayDemiPlane(this, p0, p1, p2, NULL);
	DelaunayDemiPlane * pl1 = new DelaunayDemiPlane(this, p0, p2, p1, NULL);
	DelaunayDemiPlane * pl2 = new DelaunayDemiPlane(this, p1, p2, p0, NULL);
	
	makeNeighbours(t,pl0 ) ;
	makeNeighbours(t,pl1) ;
	makeNeighbours(t,pl2) ;
	makeNeighbours(pl1,pl0 ) ;
	makeNeighbours(pl2,pl1) ;
	makeNeighbours(pl0,pl2) ;

	son.push_back(t) ;
	son.push_back(pl0) ;
	son.push_back(pl1) ;
	son.push_back(pl2) ;
	kill(p0) ;
}
	
void DelaunayRoot::print() const
{
	std::cerr << "I am root !" << std::endl ;
}


DelaunayTreeItem * DelaunayRoot::getSon(size_t i)
{
	return son[i] ;
}
	
bool DelaunayRoot::isVertex(const Point *p) const
{
	return false ;
}
	
bool DelaunayRoot::inCircumCircle(const Point & p) const 
{
	return true ;
}
	
std::pair< Point*,  Point*> DelaunayRoot::nearestEdge(const Point p)
{
	assert(false) ;
	return std::pair< Point*,  Point*>(NULL, NULL) ;
}
	
void DelaunayRoot::insert(std::vector<DelaunayTreeItem *> & ret,Point *p, Star *s)
{
	
	for (size_t i  = 0 ;  i < son.size() ; i++)
	{
		
		std::vector<DelaunayTreeItem *> temp ;
		son[i]->insert(temp, p, s) ;
		
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
		son[i]->conflicts(temp,g) ;
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
		son[i]->conflicts(temp,p) ;
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
	DelaunayRoot *root = new DelaunayRoot( p0, p1, p2) ;
	tree.push_back(root) ;
	tree.push_back(root->getSon(0)) ;
	tree.push_back(root->getSon(1)) ;
	tree.push_back(root->getSon(2)) ;
	tree.push_back(root->getSon(3)) ;
	plane.push_back(dynamic_cast<DelaunayDemiPlane *>(root->getSon(1))) ;
	plane.push_back(dynamic_cast<DelaunayDemiPlane *>(root->getSon(2))) ;
	plane.push_back(dynamic_cast<DelaunayDemiPlane *>(root->getSon(3))) ;
}
	
DelaunayTree::~DelaunayTree() 
{ 	
	std::vector<Point *> pts ;
	std::valarray<Point *> nularray(0) ;
	for(size_t i = 0 ;  i < this->tree.size() ; i++)
	{

		if(this->tree[i]->isTriangle)
		{
			DelaunayTriangle * t = dynamic_cast<DelaunayTriangle *>(tree[i]) ;
			
			for(size_t j = 0 ; j < t->getBoundingPoints().size() ;j++)
				pts.push_back(&t->getBoundingPoint(j)) ;
			
			t->setBoundingPoints(nularray) ;
		}

	}
	
	for(size_t i = 0 ;  i < this->tree.size() ; i++)
	{
		delete this->tree[i] ;
	}
	
	std::sort(pts.begin(), pts.end()) ;
	std::vector<Point *>::iterator e = std::unique(pts.begin(), pts.end()) ;
	pts.erase(e, pts.end()) ;
	
	for(size_t i = 0 ;  i < pts.size() ; i++)
	{
		delete pts[i] ;
	}
	
// 	for(size_t i = 0 ;  i < arr.size() ; i++)
// 	{
// 		delete arr[i] ;
// 	}
	
// 	for(size_t i = 0 ;  i < this->plane.size() ; i++)
// 	{
// // 		for(size_t j = 0 ; j < this->tree[i]->size() ; j++)
// // 			delete this->tree[i]->getPoint(j) ;
// // 		
// 		delete this->plane[i] ;
// 	}
} ;

void DelaunayTree::insertIf( Point *p, std::vector<SamplingCriterion *> v, double minScore )
{
	std::vector<DelaunayTreeItem *> cons = this->conflicts(p) ;
	
	
	for(size_t i = 0 ; i < cons.size() ; i++)
	{
		cons[i]->deadneighbour.clear() ;
	}
	
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
			
			if(ret[j]->neighbour.size() > 1 || ret[j]->neighbour.empty())
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
		
		if(ret[i]->isAlive())
			tree.push_back(ret[i]) ;
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
	
	for(size_t i = 0 ; i < cons.size() ; i++)
	{
		cons[i]->deadneighbour.clear() ;
	}
	

	Star * s = new Star(&cons, p) ;
	
	std::vector<DelaunayTreeItem *> ret ;
	
	for(size_t i = 0 ; i < cons.size() ; i++)
	{
		std::vector<DelaunayTreeItem *> temp ;
		cons[i]->insert(temp,p, s) ;
		ret.insert(ret.end(), temp.begin(), temp.end()) ;
	}
	
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
	
	
	for(size_t i = 0 ; i < ret.size() ; i++)
	{		
		assert((ret[i]->isTriangle && ret[i]->neighbour.size() ==3) || ret[i]->isPlane) ;
		tree.push_back(ret[i]) ;
	}

	std::vector<DelaunayDemiPlane *> * hull = this->getConvexHull() ;
	plane.clear() ;
	plane.insert(plane.end(), hull->begin(), hull->end()) ;
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
				std::pair< std::vector<DelaunayTreeItem *>, std::vector<DelaunayTreeItem *> > temp ;
				plane[i]->neighbour[j]->conflicts(temp,p) ;
				
				cons.first.insert(cons.first.end(), temp.first.begin(),temp.first.end()) ;
				cons.second.insert(cons.second.end(), temp.second.begin(),temp.second.end()) ;
			}
			
		}
		
	}
	
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
				plane[i]->neighbour[j]->conflicts(temp,s) ;
				
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
	Circle neighbourhood(0.01, g->getCenter()) ; 
//magic value : basically, we don't want to have the center a vertex of the mesh.
	
	std::vector<DelaunayTriangle *> ret ;
	std::pair< std::vector<DelaunayTreeItem *>,std::vector<DelaunayTreeItem *> > cons ;
	if(!tree.empty())
		this->tree[0]->conflicts(cons, &g->getCenter()) ;
	
// 	for(size_t i = 0 ; i < plane.size() ; i++)
// 	{
// 	
// 		if(!plane[i]->visited)
// 		{
// 			
// 			for(size_t j = 0 ; j < plane[i]->neighbour.size() ; j++)
// 			{
// 				std::pair< std::vector<DelaunayTriangle *>, std::vector<DelaunayTreeItem *> >  temp = plane[i]->neighbour[j]->conflicts(g) ;
// 				
// 				cons.first.insert(cons.first.end(), temp.first.begin(),temp.first.end()) ;
// 				cons.second.insert(cons.second.end(), temp.second.begin(),temp.second.end()) ;
// 			}
// 			
// 		}
// 		
// 	}
	
	std::set<DelaunayTriangle *> toCheck ;
	
	for(size_t i = 0 ; i < cons.second.size() ; i++)
		cons.second[i]->clearVisited() ;
	
	
	for(size_t i = 0 ; i < cons.first.size() ; i++)
	{
		if(cons.first[i]->isTriangle && cons.first[i]->isConflicting(g))
		{
			DelaunayTriangle * t = static_cast<DelaunayTriangle *>(cons.first[i]) ;
			t->visited = true ;
			cons.second.push_back(t) ;
			ret.push_back(t) ;
			for(size_t j = 0 ; j < t->neighbourhood.size() ; j++)
			{
				if(!t->neighbourhood[j]->visited && !t->neighbourhood[j]->visited)
				{
					t->neighbourhood[j]->visited = true ;
					toCheck.insert(t->neighbourhood[j]) ;
				}
			}
		}
	}
	
	cons.second.clear() ;
	for(size_t i = 0 ; i < cons.first.size() ; ++i)
		cons.first[i]->clearVisited() ;
	
	while(!toCheck.empty())
	{
		std::set<DelaunayTriangle *> temp ;
		
		for(std::set<DelaunayTriangle *>::iterator i = toCheck.begin() ; i != toCheck.end() ; ++i)
		{
			(*i)->visited = true ;
			cons.second.push_back(*i) ;
			
			if((*i)->isConflicting(g))
			{
				ret.push_back(*i) ;
			}
		}
		
		for(std::set<DelaunayTriangle *>::iterator i = toCheck.begin() ; i != toCheck.end() ; ++i)
		{

			for(size_t j = 0 ; j< (*i)->neighbourhood.size() ; j++)
			{
				if(!(*i)->neighbourhood[j]->visited )
				{
					(*i)->neighbourhood[j]->visited = true ;
					temp.insert((*i)->neighbourhood[j]) ;
				}
			}
		}
		
		toCheck = temp ;
		
		for(std::set<DelaunayTriangle *>::iterator i = toCheck.begin() ; i != toCheck.end() ; ++i)
			(*i)->clearVisited() ;
		
	}
		
	std::sort(ret.begin(), ret.end()) ;
	std::vector<DelaunayTriangle *>::iterator en = std::unique(ret.begin(), ret.end()) ;
	ret.erase(en, ret.end()) ;
	
	for(size_t i = 0 ; i < cons.second.size() ; i++)
		cons.second[i]->clearVisited() ;
	
	return ret ;
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
		if(tree[i]->isAlive() && tree[i]->isTriangle)
		{
			ret.push_back((DelaunayTriangle *)(tree[i])) ;
		}
	}

	if(!neighbourhood && buildNeighbourhood)
	{
		std::cerr << "\r building neighbourhood... element 0/" << ret.size() << std::flush ;
		for( size_t i = 0 ; i < ret.size() ;i++)
		{
			ret[i]->neighbourhood.clear() ;
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
				if(ret[i]->neighbour[j]->isTriangle && !ret[i]->neighbour[j]->visited)
				{
					tocheck.push_back(static_cast<DelaunayTriangle *>(ret[i]->neighbour[j]));
					ret[i]->neighbour[j]->visited = true ;
					toclean.push_back(*tocheck.rbegin()) ;
					ret[i]->neighbourhood.push_back(*tocheck.rbegin()) ;
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
						    tocheck[k]->neighbour[j]->isTriangle 
						    && !tocheck[k]->neighbour[j]->visited
						    && tocheck[k]->neighbour[j] != ret[i]
						    && static_cast<DelaunayTriangle *>(tocheck[k]->neighbour[j])->isInNeighbourhood(ret[i])
						  )
						{
							tocheck_temp.push_back(static_cast<DelaunayTriangle *>(tocheck[k]->neighbour[j]));
							tocheck[k]->neighbour[j]->visited = true ;
							toclean.push_back(*tocheck_temp.rbegin()) ;
							ret[i]->neighbourhood.push_back(*tocheck_temp.rbegin()) ;
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
	#ifdef DEBUG
	for(size_t i = 0 ; i < tree.size() ; i++)
	{
		if (tree[i]->isAlive() && tree[i]->neighbour.size() > 3 )
		{
			tree[i]->print() ;
			for(size_t j = 0 ; j< tree[i]->neighbour.size() ; j++)
			{
				std::cerr << "\t --> " ;
				tree[i]->neighbour[j]->print() ;
			}
		}
	}
	#endif
}

void DelaunayTriangle::refresh(const TriElement * father)
{
	this->TriElement::refresh(father) ;
	
	if(!cachedGPs)
		cachedGPs = new std::valarray<std::pair<Point, double> >(getSubTriangulatedGaussPoints()) ; 
	
// 	this->computeCenter() ;
}

std::vector<std::vector<Matrix> > DelaunayTriangle::getElementaryMatrix() const 
{

	std::vector<std::vector<Matrix > > mother ;
	std::valarray<std::pair<Point, double> > gp =getSubTriangulatedGaussPoints() ;
	std::valarray<Matrix> Jinv ;
	std::vector<std::pair<size_t, Function> > dofs = getDofs() ;

// 	if(moved)
// 	{
		Jinv.resize(gp.size()) ;
		for(size_t i = 0 ; i < gp.size() ;  i++)
		{
			Jinv[i] = getInverseJacobianMatrix( gp[i].first ) ;
		}
// 	}
// 	else
// 	{
// 		Matrix J = this->getInverseJacobianMatrix(Point( 1./3.,1./3.) ) ;
// 		Jinv.resize(gp.size()) ;
// 		for(size_t i = 0 ; i < gp.size() ;  i++)
// 		{
// 			Jinv[i] = J ;
// 		}
// 	}
	
	for(size_t i = 0 ; i < dofs.size() ; i++)
	{
		std::vector< Matrix > v_j ;
		
		for(size_t j = 0 ; j < dofs.size() ; j++)
		{
			v_j.push_back(Matrix()) ;
		}
		
		mother.push_back(v_j) ;
	}
	
	
	for(size_t i = 0 ; i < dofs.size() ; i++)
	{
		mother[i][i] = behaviour->apply(dofs[i].second, dofs[i].second,gp, Jinv) ;
		
		for(size_t j = i+1 ; j < dofs.size() ; j++)
		{
			mother[i][j] = behaviour->apply(dofs[i].second, dofs[j].second,gp, Jinv) ;
			mother[j][i] = mother[i][j].transpose() ;
		}
	}

	return mother ;
}

std::vector<std::vector<Matrix> > DelaunayTriangle::getNonLinearElementaryMatrix() const 
{
	std::vector<std::pair<size_t, Function> > dofs = getDofs() ;
	std::vector<std::vector<Matrix> > mother ;
	
	if(state == NULL)
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
	std::valarray<std::pair<Point, double> > gp = getSubTriangulatedGaussPoints() ;
	this->nonlinbehaviour->setState(this->getState()) ;
	
// 	if(moved)
// 	{
		Jinv.resize(gp.size()) ;
		for(size_t i = 0 ; i < gp.size() ;  i++)
		{
			Jinv[i] = getInverseJacobianMatrix( gp[i].first ) ;
		}
// 	}
// 	else
// 	{
// 		Matrix J = this->getInverseJacobianMatrix(Point( 1./3.,1./3.) ) ;
// 		Jinv.resize(gp.size()) ;
// 		for(size_t i = 0 ; i < gp.size() ;  i++)
// 		{
// 			Jinv[i] = J ;
// 		}
// 	}
	
	for(size_t i = 0 ; i < dofs.size() ; i++)
	{
		std::vector< Matrix > v_j ;
		
		for(size_t j = 0 ; j < dofs.size() ; j++)
		{
			v_j.push_back(Matrix()) ;
		}
		
		mother.push_back(v_j) ;
	}
	
	
	for(size_t i = 0 ; i < dofs.size() ; i++)
	{
		mother[i][i] =nonlinbehaviour->apply(dofs[i].second, dofs[i].second,gp, Jinv) ;
		
		for(size_t j = i+1 ; j < dofs.size() ; j++)
		{
			mother[i][j] = nonlinbehaviour->apply(dofs[i].second, dofs[j].second,gp, Jinv) ;
			mother[j][i] = mother[i][j].transpose() ;
		}
	}
	return mother ;
}



std::valarray<std::pair<Point, double> > DelaunayTriangle::getSubTriangulatedGaussPoints() const
{
	std::valarray<std::pair<Point, double> > gp = getGaussPoints() ; 
	
	VirtualMachine vm ;
	
	if(getEnrichmentFunctions().size() > 0 )
	{
		std::vector<std::pair<Point, double> > gp_alternative ;
		std::vector<Point> to_add ;
		std::vector<Point> to_add_extra ;
		to_add.push_back(Point(0,1)) ;
		to_add.push_back(Point(0,0)) ;
		to_add.push_back(Point(1,0)) ;
		TriElement father(LINEAR) ;
		
		for(size_t i = 0 ; i <  getEnrichmentFunctions().size() ; i++)
		{
			for(size_t j = 0 ; j < getEnrichmentFunction(i).second.getIntegrationHint().size() ; j++)
			{
				if(squareDist(getEnrichmentFunction(i).second.getIntegrationHint(j), to_add[0]) > 1e-8 && 
				   squareDist(getEnrichmentFunction(i).second.getIntegrationHint(j), to_add[1]) > 1e-8 && 
				   squareDist(getEnrichmentFunction(i).second.getIntegrationHint(j), to_add[2]) > 1e-8 &&
				   father.in(getEnrichmentFunction(i).second.getIntegrationHint(j)) )
				{
					bool ok = true ;
					for(size_t k = 0 ; k < to_add_extra.size() ; k++)
					{
						if(squareDist(getEnrichmentFunction(i).second.getIntegrationHint(j), to_add_extra[k]) < 1e-8)
						{
							ok = false ;
							break ;
						}
					}
					if(ok && father.in(getEnrichmentFunction(i).second.getIntegrationHint(j)))
						to_add_extra.push_back(getEnrichmentFunction(i).second.getIntegrationHint(j)) ;
				}
			}
		}
		
		to_add.insert(to_add.end(), to_add_extra.begin(),to_add_extra.end() ) ;
		
		DelaunayTree dt(new Point(to_add[0]), new Point(to_add[1]), new Point(to_add[2])) ;
		for(size_t i = 3 ; i < to_add.size() ; i++)
		{
			dt.insert(new Point(to_add[i])) ;
		}
		
		std::vector<DelaunayTriangle *> tri = dt.getTriangles(false) ;

		size_t numberOfRefinements = 0 ;
		
		for(size_t i = 0 ; i < numberOfRefinements ; i++)
		{
			tri = dt.getTriangles(false) ;
			std::vector<Point> quadtree ;
			for(size_t i = 0 ; i < tri.size() ; i++)
			{
				quadtree.push_back((*tri[i]->first+*tri[i]->second)*.5) ;
				quadtree.push_back((*tri[i]->first+*tri[i]->third)*.5) ;
				quadtree.push_back((*tri[i]->third+*tri[i]->second)*.5) ;
			}
			std::stable_sort(quadtree.begin(), quadtree.end()) ;
			std::vector<Point>::iterator e = std::unique(quadtree.begin(), quadtree.end()) ;
			quadtree.erase(e, quadtree.end()) ;
			
			for(size_t i = 0 ; i < quadtree.size() ; i++)
			{
				dt.insert(new Point(quadtree[i])) ;
			}
		}
		
// 		dt.addSharedNodes(1) ;
		tri = dt.getTriangles(false) ;
		dt.refresh( &father, false) ;

// 		if(moved)
// 		{
			for(size_t i = 0 ; i < tri.size() ; i++)
			{

// 				double jmin =  (*tri)[i]->jacobianAtPoint(Point(1./3.,1./3.)) ;
				Function x = tri[i]->getXTransform() ;
				Function y = tri[i]->getYTransform() ;
				tri[i]->setOrder(QUADRATIC) ;

				std::valarray<std::pair<Point, double> > gp_temp = tri[i]->getGaussPoints() ;
				
				for(size_t j = 0 ; j < gp_temp.size() ; j++)
				{
// 					gp_temp[j].second /= jmin ;
					gp_temp[j].first.set(vm.eval(x, gp_temp[j].first), vm.eval(y, gp_temp[j].first)) ;
					gp_temp[j].second *= this->jacobianAtPoint(gp_temp[j].first) ;
					gp_alternative.push_back(gp_temp[j]) ;
				}
			}

// 		}
// 		else
// 		{
// 			double ja = this->jacobianAtPoint(Point(1./3.,1./3.)) ;
// 			for(size_t i = 0 ; i < tri.size() ; i++)
// 			{
// 
// // 				double jmin =  (*tri)[i]->jacobianAtPoint(Point(1./3.,1./3.)) ;
// 				Function x = tri[i]->getXTransform() ;
// 				Function y = tri[i]->getYTransform() ;
// 				tri[i]->setOrder(QUADRATIC) ;
// 
// 				std::valarray<std::pair<Point, double> > gp_temp = tri[i]->getGaussPoints() ;
// 				
// 				for(size_t j = 0 ; j < gp_temp.size() ; j++)
// 				{
// // 					gp_temp[j].second /= jmin ;
// 					gp_temp[j].second *= ja ;
// 					gp_temp[j].first.set(vm.eval(x, gp_temp[j].first), vm.eval(y, gp_temp[j].first)) ;
// 					gp_alternative.push_back(gp_temp[j]) ;
// 				}
// 			}
// 		}
		
		
		gp.resize(gp_alternative.size()) ;
		std::copy(gp_alternative.begin(), gp_alternative.end(), &gp[0]);
		
	}
	
	return gp ;
}



Vector DelaunayTriangle::getNonLinearForces() const
{
	std::vector<std::pair<size_t, Function> > dofs = getDofs() ;
	Vector forces(dofs.size()*2) ;
	
	if(state == NULL)
	{
		std::cerr << "no state " << std::endl ;
		forces = 0 ;
		return forces ;
	}
	
	if(!this->getNonLinearBehaviour()->isActive())
	{
// 		std::cerr << "not active" << std::endl ;
		forces = 0 ;
		return forces ;
	}
	
	std::valarray<Matrix> Jinv ;
	
	std::valarray<std::pair<Point, double> > gp  = getSubTriangulatedGaussPoints() ;
	
// 	if(moved)
// 	{
		Jinv.resize(gp.size()) ;
		for(size_t i = 0 ; i < gp.size() ;  i++)
		{
			Jinv[i] = getInverseJacobianMatrix( gp[i].first ) ;
		}
// 	}
// 	else
// 	{
// 		Matrix J = this->getInverseJacobianMatrix(Point( 1./3.,1./3.) ) ;
// 		Jinv.resize(gp.size()) ;
// 		for(size_t i = 0 ; i < gp.size() ;  i++)
// 		{
// 			Jinv[i] = J ;
// 		}
// 	}
	
	
	for(size_t i = 0 ; i < dofs.size() ; i++)
	{
		for(size_t j = 0 ; j < dofs.size() ; i++)
		{
			Vector f = getNonLinearBehaviour()->getForces(this->getState(), dofs[i].second, dofs[j].second,gp, Jinv) ;
			forces[i*2] += f[0] ;
			forces[i*2+1] += f[1] ;
		}
	}
	return forces ;
}

Vector DelaunayTriangle::getForces() const
{
	std::vector<std::pair<size_t, Function> > dofs = getDofs() ;
	Vector forces(dofs.size()*2) ;
	
	std::valarray<Matrix> Jinv ;
	std::valarray<std::pair<Point, double> > gp = getSubTriangulatedGaussPoints() ;
	
// 	if(moved)
// 	{
		Jinv.resize(gp.size()) ;
		for(size_t i = 0 ; i < gp.size() ;  i++)
		{
			Jinv[i] = getInverseJacobianMatrix( gp[i].first ) ;
		}
// 	}
// 	else
// 	{
// 		Matrix J = this->getInverseJacobianMatrix(Point( 1./3.,1./3.) ) ;
// 
// 		Jinv.resize(gp.size()) ;
// 		for(size_t i = 0 ; i < gp.size() ;  i++)
// 		{
// 			Jinv[i] = J ;
// 		}
// 	}
// 	
	size_t numdof = getBehaviour()->getNumberOfDegreesOfFreedom() ;
	int offset = numdof-1 ;
	
	for(size_t i = 0 ; i < dofs.size() ; i++)
	{
		for(size_t j = 0 ; j < dofs.size() ; j++)
		{
			Vector f = behaviour->getForces( this->getState(), dofs[i].second, dofs[j].second,gp, Jinv) ;
			forces[i*numdof] += f[0] ;
			if(offset)
				forces[i*numdof+offset] += f[offset] ;
		}
	}
	
	return forces ;
}

