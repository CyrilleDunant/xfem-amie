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
	
bool DelaunayTriangle::isConflicting(const Geometry * g) const
{

	return g->in(*first) || g->in(*second) || g->in(*third) || inCircumCircle(g->getCenter()) ;
	
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
						for(size_t j = 0 ; j < tri[i]->neighbourhood.size() ; j++)
						{
							if(tri[i]->neighbourhood[j]->visited)
							{
								DelaunayTriangle * n = tri[i]->neighbourhood[j] ;
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
			(static_cast<DelaunayTriangle *>(this->tree[i]))->refresh(father) ;
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
	
	
	if(!isConflicting(g) )
	{
		return ;
	}
	
	for (size_t i  = 0 ;  i < stepson.size() ; i++)
	{
		bool limit = false ;
		if(!stepson[i]->visited && stepson[i]->isTriangle)
		{
			DelaunayTriangle * t = static_cast<DelaunayTriangle *>(stepson[i]) ;
			limit = std::abs(squareDist2D(t->getCircumCenter(),g->getCenter())-(t->getRadius()+g->getRadius())*(t->getRadius()+g->getRadius())) < POINT_TOLERANCE ;
		}
		
		if( !stepson[i]->visited && stepson[i]->isConflicting(g) || limit)
		{
			std::pair<std::vector<DelaunayTriangle *>, std::vector<DelaunayTreeItem *> > temp ;
			stepson[i]->conflicts(temp, g) ;
			ret.first.insert(ret.first.end(), temp.first.begin(), temp.first.end()) ;
			ret.second.insert(ret.second.end(), temp.second.begin(), temp.second.end()) ;
		}
	}
	
	for (size_t i  = 0 ;  i < son.size() ; i++)
	{
		bool limit = false ;

		if(!son[i]->visited && son[i]->isTriangle)
		{
			DelaunayTriangle * t = static_cast<DelaunayTriangle *>(son[i]) ;
			limit = std::abs(squareDist2D(t->getCircumCenter(),g->getCenter())-(t->getRadius()+g->getRadius())*(t->getRadius()+g->getRadius())) < POINT_TOLERANCE ;
		}

		if( !son[i]->visited && son[i]->isConflicting(g) || limit)
		{
			std::pair<std::vector<DelaunayTriangle *>, std::vector<DelaunayTreeItem *> > temp ;
			son[i]->conflicts(temp,g) ;
			ret.first.insert(ret.first.end(), temp.first.begin(), temp.first.end()) ;
			ret.second.insert(ret.second.end(), temp.second.begin(), temp.second.end()) ;
		}
	}
	
	if(isAlive() && isTriangle)
	{
		ret.first.push_back(static_cast<DelaunayTriangle *>(this)) ;
	}
	
	for (size_t i  = 0 ;  i < neighbour.size() ; i++)
	{
		bool limit = false ;

		if(!neighbour[i]->visited && neighbour[i]->isTriangle)
		{
			DelaunayTriangle * t = static_cast<DelaunayTriangle *>(neighbour[i]) ;
			limit = std::abs(squareDist2D(t->getCircumCenter(),g->getCenter())-(t->getRadius()+g->getRadius())*(t->getRadius()+g->getRadius())) < POINT_TOLERANCE ;
		}

		if( !neighbour[i]->visited && neighbour[i]->isConflicting(g) || limit)
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
			limit = std::abs(squareDist2D(t->getCircumCenter(),*p)-t->getRadius()*t->getRadius()) < 2.*POINT_TOLERANCE*POINT_TOLERANCE ;
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
			limit = std::abs(squareDist2D(t->getCircumCenter(),*p)-t->getRadius()*t->getRadius()) < 2.*POINT_TOLERANCE*POINT_TOLERANCE ;
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
			limit = std::abs(squareDist2D(t->getCircumCenter(),*p)-t->getRadius()*t->getRadius()) < 2.*POINT_TOLERANCE*POINT_TOLERANCE ;
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
	
	if(std::find_if(neighbour.begin(), neighbour.end(), EqItems(t, 1e-11)) != neighbour.end())
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

std::pair<Point*, Point*> DelaunayTriangle::commonEdge(const DelaunayTreeItem * t) const
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
					if(p->neighbour[i]->isNeighbour(neighbour[j]))
						makeNeighbours(neighbour[j], p->neighbour[i]) ;
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
		if(this->numberOfCommonVertices(neighbour[i]) == 2)
		{
			std::pair< Point*,  Point*> pp = this->commonEdge(neighbour[i]) ;
			if (!neighbour[i]->inCircumCircle(*p))
			{
				if(!isAligned(p, pp.first, pp.second ))
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
		std::cerr << "(" << getBoundingPoint(i).id <<  ";" << getBoundingPoint(i).t <<") " ;
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
	return (fma(vector.x,(p.y - first->y), - vector.y*(p.x - first->x)) * direction < -10*std::numeric_limits<double>::epsilon()) ;
}
	
bool DelaunayDemiPlane::isVertex(const Point *p) const
{
	return ( (*p) == (*first) || (*p) == (*second) ) || isAligned(p, first, second) ;
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
	if (visited)
		return ;
	visited = true ;
	
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
				
				if(!isAligned(p, first, second))
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
	
	assert(!isAligned(lims[0], p, lims[1])) ;
	assert((*lims[0]) != (*lims[1])) ;

	if(!isAligned(lims[0], p, lims[1]))
	{
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
	}
	else
	{
		DelaunayDemiPlane *p0 = new DelaunayDemiPlane(this, lims[0], third, lims[1], p) ;
		
		DelaunayDemiPlane *p1 = new DelaunayDemiPlane(this, lims[1], third, lims[0], p) ;
		
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
/*	std::vector<Point *> pts ;*/
	std::valarray<Point *> nularray(0) ;
	for(size_t i = 0 ;  i < this->tree.size() ; i++)
	{

		if(this->tree[i]->isTriangle)
		{
			DelaunayTriangle * t = dynamic_cast<DelaunayTriangle *>(tree[i]) ;
			
// 			for(size_t j = 0 ; j < t->getBoundingPoints().size() ;j++)
// 				pts.push_back(&t->getBoundingPoint(j)) ;
			
			t->setBoundingPoints(nularray) ;
		}

	}
	
	for(size_t i = 0 ;  i < this->tree.size() ; i++)
	{
		delete this->tree[i] ;
	}
	
// 	std::sort(pts.begin(), pts.end()) ;
// 	std::vector<Point *>::iterator e = std::unique(pts.begin(), pts.end()) ;
// 	pts.erase(e, pts.end()) ;
	
// 	for(size_t i = 0 ;  i < pts.size() ; i++)
// 	{
// 		delete pts[i] ;
// 		pts[i] = 0;
// 	}
	
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
				if(!plane[i]->neighbour[j]->visited)
				{
					std::pair< std::vector<DelaunayTreeItem *>, std::vector<DelaunayTreeItem *> > temp ;
					plane[i]->neighbour[j]->conflicts(temp,p) ;
					
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
				plane[i]->neighbour[j]->conflicts(temp,g) ;
				
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
// 	#ifdef DEBUG
	for(size_t i = 0 ; i < tree.size() ; i++)
	{
		if (tree[i]->isAlive())
		{
			tree[i]->print() ;
			for(size_t j = 0 ; j< tree[i]->neighbour.size() ; j++)
			{
				std::cerr << "\t --> " ;
				tree[i]->neighbour[j]->print() ;
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

std::vector<std::vector<Matrix> > DelaunayTriangle::getElementaryMatrix() const 
{

	std::vector<std::vector<Matrix > > mother ;
	GaussPointArray gp = getSubTriangulatedGaussPoints() ;
	std::valarray<Matrix> Jinv ;
	std::vector<size_t > dofs = getDofIds() ;

	if(moved)
	{
		Jinv.resize(gp.gaussPoints.size()) ;
		for(size_t i = 0 ; i < gp.gaussPoints.size() ;  i++)
		{
			Jinv[i] = getInverseJacobianMatrix( gp.gaussPoints[i].first ) ;
		}
	}
	else
	{
		Matrix J = this->getInverseJacobianMatrix(Point( 1./3.,1./3.) ) ;
		Jinv.resize(gp.gaussPoints.size()) ;
		for(size_t i = 0 ; i < gp.gaussPoints.size() ;  i++)
		{
			Jinv[i] = J ;
		}
	}

	int size = getBehaviour()->getNumberOfDegreesOfFreedom() ;
	
	
	for(size_t i = 0 ; i < dofs.size() ; i++)
	{
		std::vector< Matrix > v_j ;
		
		for(size_t j = 0 ; j < dofs.size() ; j++)
		{
			v_j.push_back(Matrix(size,size)) ;
		}
		
		mother.push_back(v_j) ;
	}
	
	
	for(size_t i = 0 ; i < getShapeFunctions().size() ; i++)
	{
		mother[i][i] = behaviour->apply(getShapeFunction(i), getShapeFunction(i),gp, Jinv) ;
		
		for(size_t j = i+1 ; j < getShapeFunctions().size() ; j++)
		{
			mother[i][j] = behaviour->apply(getShapeFunction(i), getShapeFunction(j),gp, Jinv) ;
			mother[j][i] = behaviour->apply(getShapeFunction(j), getShapeFunction(i),gp, Jinv) ;
		}
		for(size_t j = 0 ; j < getEnrichmentFunctions().size() ; j++)
		{
			mother[i][j+getShapeFunctions().size()] = behaviour->apply(getShapeFunction(i), getEnrichmentFunction(j),gp, Jinv) ;
			mother[j+getShapeFunctions().size()][i] = behaviour->apply(getEnrichmentFunction(j), getShapeFunction(i),gp, Jinv) ;
		}
	}
	for(size_t i = 0 ; i < getEnrichmentFunctions().size() ; i++)
	{
		mother[i+getShapeFunctions().size()][i+getShapeFunctions().size()] = behaviour->apply(getEnrichmentFunction(i), getEnrichmentFunction(i),gp, Jinv) ;
		
		for(size_t j = 0 ; j < getEnrichmentFunctions().size() ; j++)
		{
			mother[i+getShapeFunctions().size()][j+getShapeFunctions().size()] = behaviour->apply(getEnrichmentFunction(i), getEnrichmentFunction(j),gp, Jinv) ;
			mother[j+getShapeFunctions().size()][i+getShapeFunctions().size()] = behaviour->apply(getEnrichmentFunction(j), getEnrichmentFunction(i),gp, Jinv) ;
		}
	}

	return mother ;
}

std::vector<std::vector<Matrix> > DelaunayTriangle::getNonLinearElementaryMatrix() const 
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
	
// 	if(moved)
// 	{
		Jinv.resize(gp.gaussPoints.size()) ;
		for(size_t i = 0 ; i < gp.gaussPoints.size() ;  i++)
		{
			Jinv[i] = getInverseJacobianMatrix( gp.gaussPoints[i].first ) ;
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
	
	
	for(size_t i = 0 ; i < getShapeFunctions().size() ; i++)
	{
		mother[i][i] = behaviour->apply(getShapeFunction(i), getShapeFunction(i),gp, Jinv) ;
		
		for(size_t j = i+1 ; j < getShapeFunctions().size() ; j++)
		{
			mother[i][j] = behaviour->apply(getShapeFunction(i), getShapeFunction(j),gp, Jinv) ;
			mother[j][i] = behaviour->apply(getShapeFunction(j), getShapeFunction(i),gp, Jinv) ;
		}
		for(size_t j = 0 ; j < getEnrichmentFunctions().size() ; j++)
		{
			mother[i][j+getShapeFunctions().size()] = behaviour->apply(getShapeFunction(i), getEnrichmentFunction(j),gp, Jinv) ;
			mother[j+getShapeFunctions().size()][i] = behaviour->apply(getEnrichmentFunction(j), getShapeFunction(i),gp, Jinv) ;
		}
	}
	for(size_t i = 0 ; i < getEnrichmentFunctions().size() ; i++)
	{
		mother[i+getShapeFunctions().size()][i+getShapeFunctions().size()] = behaviour->apply(getEnrichmentFunction(i), getEnrichmentFunction(i),gp, Jinv) ;
		
		for(size_t j = 0 ; j < getEnrichmentFunctions().size() ; j++)
		{
			mother[i+getShapeFunctions().size()][j+getShapeFunctions().size()] = behaviour->apply(getEnrichmentFunction(i), getEnrichmentFunction(j),gp, Jinv) ;
			mother[j+getShapeFunctions().size()][i+getShapeFunctions().size()] = behaviour->apply(getEnrichmentFunction(j), getEnrichmentFunction(i),gp, Jinv) ;
		}
	}
	
	return mother ;
}



GaussPointArray DelaunayTriangle::getSubTriangulatedGaussPoints() const
{
	
	GaussPointArray gp = getGaussPoints() ; 
	
	VirtualMachine vm ;
	
	if(getEnrichmentFunctions().size() > 0)
	{
		std::vector<std::pair<Point, double> > gp_alternative ;
		std::vector<Point> to_add ;
		to_add.push_back(Point(0,1)) ;
		to_add.push_back(Point(0,0)) ;
		to_add.push_back(Point(1,0)) ;
		TriElement father(LINEAR) ;
		
		for(size_t i = 0 ; i <  getEnrichmentFunctions().size() ; i++)
		{
			for(size_t j = 0 ; j < getEnrichmentFunction(i).getIntegrationHint().size() ; j++)
			{
				bool go = true ;
				
				for(size_t k = 0 ; k < 3  ; k++ )
				{
					if(squareDist2D(getEnrichmentFunction(i).getIntegrationHint(j),to_add[k]) < 1e-12)
					{
						go = false ;
						break ;
					}
				}
				if(go)
					to_add.push_back(getEnrichmentFunction(i).getIntegrationHint(j)) ;
			}
		}

		std::sort(to_add.begin()+3, to_add.end()) ;
		std::vector<Point>::iterator e = std::unique(to_add.begin()+3, to_add.end(), PointEqTol(1e-12)) ;
		to_add.erase(e, to_add.end()) ;
		
		DelaunayTree dt(&to_add[0], &to_add[1], &to_add[2]) ;
		for(size_t i = 3 ; i < to_add.size() ; i++)
		{
			dt.insert(&to_add[i]) ;
		}
		
		std::vector<DelaunayTriangle *> tri = dt.getTriangles(false) ;
		std::vector<Point *> pointsToCleanup ;
		std::vector<DelaunayTriangle *> triangleToCleanup;
		size_t numberOfRefinements =  3;
		
		for(size_t i = 0 ; i < numberOfRefinements ; i++)
		{
			std::vector<DelaunayTriangle *> newTris ;
			for(size_t j = 0 ; j < tri.size() ; j++)
			{
				std::pair<std::vector<DelaunayTriangle *>, std::vector<Point *> > q = quad(tri[j]) ;
				newTris.insert(newTris.end(),q.first.begin(), q.first.end()) ;
				pointsToCleanup.insert(pointsToCleanup.end(),q.second.begin(), q.second.end()) ;
				if(i)
					triangleToCleanup.push_back(tri[j]) ;
			}
			tri = newTris ;
		}

		
		for(size_t i = 0 ; i < tri.size() ; i++)
			tri[i]->refresh(&father) ;
		
		double J = this->jacobianAtPoint(Point(1./3., 1./3.)) ;
		
		for(size_t i = 0 ; i < tri.size() ; i++)
		{

			Function x = tri[i]->getXTransform() ;
			Function y = tri[i]->getYTransform() ;
			tri[i]->setOrder(LINEAR) ;

			GaussPointArray gp_temp = tri[i]->getGaussPoints() ;
			
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

		if(numberOfRefinements)
		{
			std::valarray<Point *> nularray(0) ;
	
			for(size_t i = 0 ; i < triangleToCleanup.size() ; i++)
			{
				triangleToCleanup[i]->setBoundingPoints(nularray) ;
				delete triangleToCleanup[i];
			}
		
			for(size_t i = 0 ; i < tri.size() ; i++)
			{
				tri[i]->setBoundingPoints(nularray) ;
				delete tri[i] ;
			}
		}
		for(size_t i = 0 ; i < pointsToCleanup.size() ; i++)
			delete pointsToCleanup[i] ;

		gp.gaussPoints.resize(gp_alternative.size()) ;
		std::copy(gp_alternative.begin(), gp_alternative.end(), &gp.gaussPoints[0]);
		gp.id = -1 ;
	}
	
	return gp ;
}



Vector DelaunayTriangle::getNonLinearForces() const
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
		Jinv.resize(gp.gaussPoints.size()) ;
		for(size_t i = 0 ; i < gp.gaussPoints.size() ;  i++)
		{
			Jinv[i] = getInverseJacobianMatrix( gp.gaussPoints[i].first ) ;
		}
	}
	else
	{
		Matrix J = this->getInverseJacobianMatrix(Point( 1./3.,1./3.) ) ;
		Jinv.resize(gp.gaussPoints.size()) ;
		for(size_t i = 0 ; i < gp.gaussPoints.size() ;  i++)
		{
			Jinv[i] = J ;
		}
	}
	
	size_t numdof = getBehaviour()->getNumberOfDegreesOfFreedom() ;
	
	for(size_t i = 0 ; i < getShapeFunctions().size() ; i++)
	{

			Vector f = behaviour->getForces(this->getState(), getShapeFunction(i) ,gp, Jinv) ;
			
		for(size_t j = 0 ; j < numdof ; j++)
			forces[i*numdof+j]+=f[j];

	}
		
	for(size_t i = 0 ; i < getEnrichmentFunctions().size() ; i++)
	{
		Vector f = behaviour->getForces(this->getState(), getEnrichmentFunction(i) ,gp, Jinv) ;
		
		for(size_t j = 0 ; j < numdof ; j++)
			forces[(i+getShapeFunctions().size())*numdof+j]+=f[j];
	}

	
	return forces ;
}

Vector DelaunayTriangle::getForces() const
{
	std::vector<size_t> dofs = getDofIds() ;
	Vector forces(dofs.size()*2) ;
	
	std::valarray<Matrix> Jinv ;
	GaussPointArray gp = getSubTriangulatedGaussPoints() ;
	
	if(moved)
	{
		Jinv.resize(gp.gaussPoints.size()) ;
		for(size_t i = 0 ; i < gp.gaussPoints.size() ;  i++)
		{
			Jinv[i] = getInverseJacobianMatrix( gp.gaussPoints[i].first ) ;
		}
	}
	else
	{
		Matrix J = this->getInverseJacobianMatrix(Point( 1./3.,1./3.) ) ;

		Jinv.resize(gp.gaussPoints.size()) ;
		for(size_t i = 0 ; i < gp.gaussPoints.size() ;  i++)
		{
			Jinv[i] = J ;
		}
	}
// 	
	size_t numdof = getBehaviour()->getNumberOfDegreesOfFreedom() ;
	

	for(size_t i = 0 ; i < getShapeFunctions().size() ; i++)
	{

		Vector f = behaviour->getForces(this->getState(), getShapeFunction(i),gp, Jinv) ;
			
		for(size_t j = 0 ; j < numdof ; j++)
			forces[i*numdof+j]+=f[j];

	}
		
	for(size_t i = 0 ; i < getEnrichmentFunctions().size() ; i++)
	{
		Vector f = behaviour->getForces(this->getState(), getEnrichmentFunction(i) ,gp, Jinv) ;
		
		for(size_t j = 0 ; j < numdof ; j++)
			forces[(i+getShapeFunctions().size())*numdof+j]+=f[j];
	}
	
	return forces ;
}

std::pair<std::vector<DelaunayTriangle *>, std::vector<Point *> > Mu::quad(const DelaunayTriangle * t)
{
	std::vector<DelaunayTriangle* > tris ;
	std::vector<Point *> points ;

	points.push_back(new Point(*t->first + (*t->second - *t->first)*.5) ) ;
	points.push_back(new Point(*t->first + (*t->third - *t->first)*.5 )) ;
	points.push_back(new Point(*t->second + (*t->third - *t->second)*.5 )) ;
	
	tris.push_back(new DelaunayTriangle(NULL, points[0], points[1], t->first, NULL)) ;
	tris.push_back(new DelaunayTriangle(NULL,points[0], points[2], t->second, NULL)) ;
	tris.push_back(new DelaunayTriangle(NULL,points[1], points[2], t->third, NULL)) ;
	tris.push_back(new DelaunayTriangle(NULL,points[0], points[1], points[2], NULL)) ;
	return std::make_pair(tris, points) ;
}

