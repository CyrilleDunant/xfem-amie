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

#include "delaunay_3d.h"

#define DEBUG
// #undef DEBUG

using namespace Mu ;

DelaunayTreeItem_3D::DelaunayTreeItem_3D( DelaunayTreeItem_3D * father,  const Point * c) 
{
	this->stepfather = NULL ;
	this->father = father;
	this->m_c  = c ;
	this->dead = false ;
	this->erased = false ;
	this->visited = false ;
	
	this->isSpace =false ;
	this->isTetrahedron =false ;
	this->first = NULL ;
	this->second = NULL ;
	this->third = NULL ;
	this->fourth = NULL ;
}
	
DelaunayTreeItem_3D::~DelaunayTreeItem_3D()
{
// 	for(size_t i = 0 ; i < neighbour.size() ; i++)
// 	{
// 		neighbour[i]->removeNeighbour(this) ;
// 	}
}
	
const Point * DelaunayTreeItem_3D::killer() const 
{
	return m_k ;
}

const Point * DelaunayTreeItem_3D::creator() const 
{
	return m_c ;
}

void  DelaunayTreeItem_3D::setCreator(const Point * p)  
{
	m_c = p ;
}

size_t DelaunayTreeItem_3D::numberOfCommonVertices(const DelaunayTreeItem_3D * s) const
{
	if(this->isTetrahedron && s->isTetrahedron)//changed to 3D
	{
		size_t ret = 0 ;
		
		ret +=  (s->first == this->first || s->first  == this->second  || s->first  == this->third  || s->first  == this->fourth);
		ret += (s->second == this->first || s->second == this->second  || s->second == this->third  || s->second == this->fourth);
		ret +=  (s->third == this->first || s->third  == this->second  || s->third  == this->third  || s->third  == this->fourth) ;
		ret += (s->fourth == this->first || s->fourth == this->second  || s->fourth == this->third  || s->fourth == this->fourth) ;
		
		assert(ret < 5) ;
			
		return ret ;
	}
	if(this->isTetrahedron && s->isSpace)
	{
		size_t ret = 0 ;
		
		if(isCoplanar(s->first, s->second, s->third, first))
			ret++ ;
		if(isCoplanar(s->first, s->second, s->third, second))
			ret++ ;
		if(isCoplanar(s->first, s->second, s->third, third))
			ret++ ;
		if(isCoplanar(s->first, s->second, s->third, fourth))
			ret++ ;
		
		assert(ret < 4) ;
		
		return ret ;
	}
	if(this->isSpace && s->isSpace)
	{
		size_t ret = 0 ;
		
		if(isCoplanar(this->first, this->second, this->third ,s->first))
			ret++ ;
		if(isCoplanar(this->first, this->second, this->third,s->second))
			ret++ ;
		if(isCoplanar(this->first, this->second, this->third,s->third))
			ret++ ;
		
		assert(ret < 4) ;
		
		return ret ;
	}
	
	if(this->isSpace && s->isTetrahedron)
	{
		size_t ret = 0 ;
		
		if(isCoplanar(this->first, this->second, this->third ,s->first))
			ret++ ;
		if(isCoplanar(this->first, this->second, this->third,s->second))
			ret++ ;
		if(isCoplanar(this->first, this->second, this->third,s->third))
			ret++ ;
		if(isCoplanar(this->first, this->second, this->third,s->fourth))
			ret++ ;
		
		assert(ret < 4) ;
		
		return ret ;
	}
	
	return 0 ;
}


void Star_3D::updateNeighbourhood()
{
	for(size_t i = 0 ; i < treeitem.size() ;i++)
	{
		for(size_t j = i+1 ; j < treeitem.size() ;j++)
		{
			if(treeitem[i]->isDuplicate(treeitem[j]) && treeitem[i]->isAlive() && treeitem[j]->isAlive())
			{
				for(size_t k = 0 ; k < treeitem[j]->neighbour.size() ; k++)
				{
					makeNeighbours(treeitem[j]->neighbour[k], treeitem[i]) ;
				}
				
				treeitem[j]->kill(creator) ;
				treeitem[j]->erased = true ;

				treeitem[j]->father->son.erase(
					                    std::find(
					                              treeitem[j]->father->son.begin(), 
					                              treeitem[j]->father->son.end(), 
					                              treeitem[j]
					                              )
					                    ) ;

				treeitem[j]->stepfather->stepson.erase(
					                            std::find(
					                                      treeitem[j]->stepfather->stepson.begin(),
						                              treeitem[j]->stepfather->stepson.end(), 
						                              treeitem[j]
					                                     )
					                            ) ;
				treeitem[j]->father->addSon(treeitem[i]) ;
				treeitem[j]->stepfather->addStepson(treeitem[i]) ;
				treeitem[i]->stepfather = treeitem[j]->stepfather ;
			}
		}
	}	
	
	for(size_t i = 0 ; i < treeitem.size() ;i++)
	{
		for(size_t j = i+1 ; j < treeitem.size() ;j++)
		{
			if(!treeitem[i]->erased && !treeitem[j]->erased )
				makeNeighbours(treeitem[i], treeitem[j]) ;
		}
	}
}


void DelaunayTree_3D::addSharedNodes(size_t nodes_per_side, size_t time_planes, double timestep, const TetrahedralElement * father)
{
	std::vector<DelaunayTetrahedron *> tri = getTetrahedrons() ;
	
	if(nodes_per_side> 1)
		nodes_per_side = 1 ;
	

	for(size_t i = 0 ; i < tri.size() ; i++)
	{
		std::vector<std::pair<Point, Point> > sides ;
		sides.push_back(std::make_pair(*tri[i]->first,*tri[i]->second)) ;
		sides.push_back(std::make_pair(*tri[i]->second,*tri[i]->third)) ;
		sides.push_back(std::make_pair(*tri[i]->third,*tri[i]->fourth)) ;
		sides.push_back(std::make_pair(*tri[i]->fourth,*tri[i]->first)) ;
		sides.push_back(std::make_pair(*tri[i]->fourth,*tri[i]->second)) ;
		sides.push_back(std::make_pair(*tri[i]->first,*tri[i]->third)) ;
		std::vector<size_t> positions ;
		if(nodes_per_side)
		{
			positions.push_back(0) ;
			positions.push_back(1) ;
			positions.push_back(2) ;
			
			positions.push_back(2) ;
			positions.push_back(3) ;
			positions.push_back(4) ;
			
			positions.push_back(4) ;
			positions.push_back(5) ;
			positions.push_back(6) ;
			
			positions.push_back(0) ;
			positions.push_back(7) ;
			positions.push_back(6) ;
			
			positions.push_back(6) ;
			positions.push_back(8) ;
			positions.push_back(2) ;
			
			positions.push_back(0) ;
			positions.push_back(9) ;
			positions.push_back(4) ;
		}
		else
		{
			positions.push_back(0) ;
			positions.push_back(1) ;
			positions.push_back(1) ;
			positions.push_back(2) ;
			positions.push_back(2) ;
			positions.push_back(3) ;
			positions.push_back(3) ;
			positions.push_back(0) ;
			positions.push_back(3) ;
			positions.push_back(1) ;
			positions.push_back(0) ;
			positions.push_back(2) ;
		}
		
		tri[i]->visited = true ;
		
		size_t nodes_per_plane = nodes_per_side*6+4 ;
		
		std::valarray<Point *> newPoints((Point *)NULL,nodes_per_plane*time_planes ) ;
		std::valarray<bool> done(false, nodes_per_plane*time_planes) ;
				
		for(size_t plane = 0 ; plane < time_planes ; plane++)
		{
			size_t current = 0 ;
			for(size_t side = 0 ; side < 6 ; side++)
			{
				Point a(sides[side].first) ;
				Point b(sides[side].second) ;
				if(time_planes> 1)
				{
					a.t = (double)plane*(timestep/(double)(time_planes-1))-timestep/2.;
					b.t = (double)plane*(timestep/(double)(time_planes-1))-timestep/2.;
// 					std::cout << a.t << std::endl ;
				}
				for(size_t node = 0 ; node < nodes_per_side+2 ; node++)
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
								for(size_t k = 0 ; k < tri[i]->neighbourhood[j]->getBoundingPoints().size() ;k++)
								{
									if(tri[i]->neighbourhood[j]->getBoundingPoint(k) ==  proto)
									{
										foundPoint = &tri[i]->neighbourhood[j]->getBoundingPoint(k) ;
										break ;
									}
								}
								
								if(foundPoint)
								{
									break ;
								}
							}
						}
					}
					
					if(!done[positions[current]+plane*nodes_per_plane])
					{
						if(foundPoint)
						{
							newPoints[positions[current]+plane*nodes_per_plane]  = foundPoint ;
						}
						else
						{
							newPoints[positions[current]+plane*nodes_per_plane]  = new Point(proto) ;
							newPoints[positions[current]+plane*nodes_per_plane]->id = global_counter++ ;
						}
						
						done[positions[current]+plane*nodes_per_plane] = true ;
					}
					
					current++ ;
				}
			}
		}

		for(size_t k = 0 ; k< newPoints.size() ; k++)
			if(!newPoints[k])
			{
				std::cout << "ouch !" << std::endl ;
				for(size_t k = 0 ; k< newPoints.size() ; k++)
					if(newPoints[k])
						newPoints[k]->print() ;

				exit(0) ;
			}
		tri[i]->setBoundingPoints(newPoints) ;
	}
	
	
	for(size_t i = 0 ; i < tri.size() ; i++)
	{
		tri[i]->clearVisited() ;
		if(father != NULL)
			refresh(father) ;
	}
}

void DelaunayTree_3D::refresh(const TetrahedralElement *father)
{
	for(size_t i = 0 ; i < this->tree.size() ; i++)
	{
		if(this->tree[i]->isAlive() && this->tree[i]->isTetrahedron)
		{
			static_cast<DelaunayTetrahedron *>(this->tree[i])->refresh(father) ;
			
		}
	}
}

size_t DelaunayTree_3D::numPoints() const
{
	return this->global_counter ;
}

void DelaunayTreeItem_3D::conflicts(std::pair<std::vector<DelaunayTetrahedron *>, std::vector<DelaunayTreeItem_3D *> > & ret, const Geometry *g)
{

	if(visited)
	{
		return  ;
	}
	visited = true ;
	
	ret.second.push_back(this) ;
	
	
	if(!inCircumSphere(g->getCenter()) && (!g->in(*first) && !g->in(*second) && !g->in(*third) && ( !g->in(*fourth) && isTetrahedron ) || ( g->in(*fourth) && isSpace )) )
	{
		return  ;
	}
	
	for (size_t i  = 0 ;  i < stepson.size() ; i++)
	{
		if( !stepson[i]->visited && !(!stepson[i]->inCircumSphere(g->getCenter()) && (!g->in(*stepson[i]->first) && !g->in(*stepson[i]->second) && !g->in(*stepson[i]->third) && ( !g->in(*stepson[i]->fourth) && stepson[i]->isTetrahedron ) || ( g->in(*stepson[i]->fourth) && stepson[i]->isSpace ))))
		{
			std::pair<std::vector<DelaunayTetrahedron *>, std::vector<DelaunayTreeItem_3D *> > temp ;
			stepson[i]->conflicts(temp, g) ;
			ret.first.insert(ret.first.end(), temp.first.begin(), temp.first.end()) ;
			ret.second.insert(ret.second.end(), temp.second.begin(), temp.second.end()) ;
		}
	}
	
	for (size_t i  = 0 ;  i < son.size() ; i++)
	{
		if( !son[i]->visited && !(!son[i]->inCircumSphere(g->getCenter()) && (!g->in(*son[i]->first) && !g->in(*son[i]->second) && !g->in(*son[i]->third) && ( !g->in(*son[i]->fourth) && son[i]->isTetrahedron ) || ( g->in(*son[i]->fourth) && son[i]->isSpace ))))
		{
			std::pair<std::vector<DelaunayTetrahedron *>, std::vector<DelaunayTreeItem_3D *> >  temp  ;
			son[i]->conflicts(temp, g) ;
			ret.first.insert(ret.first.end(), temp.first.begin(), temp.first.end()) ;
			ret.second.insert(ret.second.end(), temp.second.begin(), temp.second.end()) ;
		}
	}
	
	if(isAlive() && isTetrahedron)
	{
		ret.first.push_back(static_cast<DelaunayTetrahedron *>(this)) ;
	}
	
	for (size_t i  = 0 ;  i < neighbour.size() ; i++)
	{
		if( !neighbour[i]->visited && !(!neighbour[i]->inCircumSphere(g->getCenter()) && (!g->in(*neighbour[i]->first) && !g->in(*neighbour[i]->second) && !g->in(*neighbour[i]->third) && ( !g->in(*neighbour[i]->fourth) && neighbour[i]->isTetrahedron ) || ( g->in(*neighbour[i]->fourth) && neighbour[i]->isSpace ))))
		{
			std::pair<std::vector<DelaunayTetrahedron *>, std::vector<DelaunayTreeItem_3D *> > temp  ;
			neighbour[i]->conflicts(temp, g) ;
			ret.first.insert(ret.first.end(), temp.first.begin(), temp.first.end()) ;
			ret.second.insert(ret.second.end(), temp.second.begin(), temp.second.end()) ;
			
		}
	}
	return  ;
}


 void DelaunayTreeItem_3D::conflicts(std::pair<std::vector<DelaunayTreeItem_3D *>, std::vector<DelaunayTreeItem_3D *> > & ret, const Point *p)
{

	if(visited)
	{
		return  ;
	}
	visited = true ;
	
	ret.second.push_back(this) ;
	
	for (size_t i  = 0 ;  i < stepson.size() ; i++)
	{
		bool limit = false ;
		if(!stepson[i]->visited && stepson[i]->isTetrahedron)
		{
			DelaunayTetrahedron * t = static_cast<DelaunayTetrahedron *>(stepson[i]) ;
			limit = std::abs(squareDist3D(t->getCircumCenter(),p)-t->getRadius()*t->getRadius()) 
				< POINT_TOLERANCE*t->getRadius()*t->getRadius()*t->getRadius() ;
		}

		if( (!stepson[i]->visited && stepson[i]->inCircumSphere(*p)) || limit) 
		{
			std::pair<std::vector<DelaunayTreeItem_3D *>, std::vector<DelaunayTreeItem_3D *> > temp  ;
			stepson[i]->conflicts(temp,p) ;
			ret.first.insert(ret.first.end(), temp.first.begin(), temp.first.end()) ;
			ret.second.insert(ret.second.end(), temp.second.begin(), temp.second.end()) ;
		}
	}
	
	for (size_t i  = 0 ;  i < son.size() ; i++)
	{
		
		bool limit = false ;
		if(!son[i]->visited && son[i]->isTetrahedron)
		{
			DelaunayTetrahedron * t = static_cast<DelaunayTetrahedron *>(son[i]) ;
			limit = std::abs(squareDist3D(t->getCircumCenter(),p)-t->getRadius()*t->getRadius()) 
				< POINT_TOLERANCE*t->getRadius()*t->getRadius()*t->getRadius() ;
		}

		if( (!son[i]->visited && son[i]->inCircumSphere(*p)) || limit)
		{
			std::pair<std::vector<DelaunayTreeItem_3D *>, std::vector<DelaunayTreeItem_3D *> > temp  ;
			son[i]->conflicts(temp,p) ;
			ret.first.insert(ret.first.end(), temp.first.begin(), temp.first.end()) ;
			ret.second.insert(ret.second.end(), temp.second.begin(), temp.second.end()) ;
		}
	}
	
	if(!inCircumSphere(*p))
	{
		return  ;
	}
	
	for (size_t i  = 0 ;  i < neighbour.size() ; i++)
	{
		
		bool limit = false ;
		if(!neighbour[i]->visited && neighbour[i]->isTetrahedron)
		{
			DelaunayTetrahedron * t = static_cast<DelaunayTetrahedron *>(neighbour[i]) ;
			limit = std::abs(squareDist3D(t->getCircumCenter(),p)-t->getRadius()*t->getRadius()) 
				< POINT_TOLERANCE*t->getRadius()*t->getRadius()*t->getRadius() ;
		}
		limit = true ;
		
		if( !neighbour[i]->visited && neighbour[i]->inCircumSphere(*p) || limit)
		{
			std::pair<std::vector<DelaunayTreeItem_3D *>, std::vector<DelaunayTreeItem_3D *> > temp ;
			neighbour[i]->conflicts(temp,p) ;
			ret.first.insert(ret.first.end(), temp.first.begin(), temp.first.end()) ;
			ret.second.insert(ret.second.end(), temp.second.begin(), temp.second.end()) ;
		}
	}
	


	if(isAlive())
	{
		ret.first.push_back(this) ;
	}
	

}

void DelaunayTreeItem_3D::removeNeighbour(DelaunayTreeItem_3D * t)
{
	if(std::find(neighbour.begin(), neighbour.end(),t) !=  neighbour.end())
	{
		neighbour.erase(std::find(neighbour.begin(), neighbour.end(), t)) ;
	}

}
	
void DelaunayTreeItem_3D::addNeighbour(DelaunayTreeItem_3D * t)
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

DelaunayTreeItem_3D * DelaunayTreeItem_3D::getNeighbour(size_t i)
{
	return this->neighbour[i] ;
}


void DelaunayTreeItem_3D::erase(const Point * p)
{
	dead = true ;
	m_k  = p ;
}

void DelaunayTreeItem_3D::kill(const Point * p)
{
	dead = true ;
	m_k  = p ;
	
	for(size_t i = 0 ; i < this->neighbour.size() ; i++)
	{
		this->neighbour[i]->removeNeighbour(this) ;
	}
}
	
bool DelaunayTreeItem_3D::isAlive() const
{
	return !this->dead ;
}
	
void DelaunayTreeItem_3D::addStepson(DelaunayTreeItem_3D * s)
{
	if(s == this)
		return ;
	stepson.push_back(s) ;
	s->setStepfather(this) ;
	addNeighbour(s) ;
}

void DelaunayTreeItem_3D::addSon(DelaunayTreeItem_3D * s)
{
	son.push_back(s) ;
}
	
void DelaunayTreeItem_3D::removeSon(DelaunayTreeItem_3D * t)
{
	if(std::find(son.begin(), son.end(),t) !=  son.end())
	{
		son.erase(std::find(son.begin(), son.end(), t)) ;
	}
}

void DelaunayTreeItem_3D::removeStepson(DelaunayTreeItem_3D * t)
{
	if(std::find(stepson.begin(), stepson.end(),t) !=  stepson.end())
	{
		stepson.erase(std::find(stepson.begin(), stepson.end(), t)) ;
	}
}

void DelaunayTreeItem_3D::setStepfather(DelaunayTreeItem_3D * s)
{
	stepfather = s ;
	addNeighbour(s) ;
}
	
void DelaunayTreeItem_3D::clearVisited()
{
	visited = false ;
}

DelaunayTetrahedron::~DelaunayTetrahedron()
{
}

DelaunayTetrahedron::DelaunayTetrahedron(DelaunayTreeItem_3D * father,  Point *p0,  Point *p1, Point *p2, Point *p3,  Point * c) : TetrahedralElement(p0, p1, p2, p3), DelaunayTreeItem_3D(father, c)
{
	
	first = &getBoundingPoint(0) ;
	
	second = &getBoundingPoint(1) ;
	third = &getBoundingPoint(2) ;
	fourth = &getBoundingPoint(3) ;
	
	assert(in(this->getCenter())) ;
	
	neighbour.clear() ;
	isSpace = false ;
	isTetrahedron = true ;
	assert(first->id > -1) ;
	assert(second->id > -1) ;
	assert(third->id > -1) ;
	assert(fourth->id> -1);
}

DelaunayTetrahedron::DelaunayTetrahedron(DelaunayTreeItem_3D * father,  Point *p0,  Point *p1, Point *p2, Point *p3,  Point *p4,  Point *p5, Point *p6, Point *p7,Point * c) : TetrahedralElement(p0, p1, p2, p3, p4, p5, p6, p7), DelaunayTreeItem_3D(father, c)
{
	
	first = &getBoundingPoint(0) ;
	
	second = &getBoundingPoint(1) ;
	third = &getBoundingPoint(2) ;
	fourth = &getBoundingPoint(3) ;
	
	assert(in(this->getCenter())) ;
	
	neighbour.clear() ;
	isSpace = false ;
	isTetrahedron = true ;
	assert(first->id > -1) ;
	assert(second->id > -1) ;
	assert(third->id > -1) ;
	assert(fourth->id> -1);
}

DelaunayTetrahedron::DelaunayTetrahedron() : DelaunayTreeItem_3D(NULL, NULL)
{
	first = &getBoundingPoint(0) ;
	second = &getBoundingPoint(1) ;
	third = &getBoundingPoint(2) ;
	fourth = &getBoundingPoint(3) ;
	assert(in(this->getCenter())) ;
	
	neighbour.clear() ;
	isSpace = false ;
	isTetrahedron = true ;
}

DelaunayDemiSpace::~DelaunayDemiSpace()
{
}
	
bool DelaunayTetrahedron::isVertex(const Point * p) const
{
	return *p == *first || *p == *second || *p == *third || *p == *fourth ;
	return (dist(p,first) < 0.0000001 || dist(p,second) < 0.0000001 || dist(p,third) < 0.0000001||dist(p,fourth) < 0.0000001) ;
}

bool DelaunayDemiSpace::isVertex(const Point * p) const
{
	return *p == *first || *p == *second || *p == *third ;
}

bool DelaunayTetrahedron::isVertexByID(const Point * p) const
{
	return (p->id == first->id || p->id == second->id || p->id == third->id|| p->id == fourth->id) ;
}
	

bool DelaunayTetrahedron::hasVertexByID(const std::valarray<Point *> * p) const
{
	for(size_t i = 0 ; i < p->size() ; i++)
	{
		if((*p)[i]->id == first->id || (*p)[i]->id == second->id || (*p)[i]->id == third->id || (*p)[i]->id == fourth->id)
			return true ;
	}
	return false ;
}

std::pair<Point*, Point*> DelaunayTetrahedron::commonEdge(const DelaunayTreeItem_3D * t) 
{
	if(t == this)
	{
		return std::pair<Point*, Point*>(NULL, NULL) ;
	}
	
	if(t->isTetrahedron)
	{
		if(this->isVertex(t->first) && this->isVertex(t->second))
			return std::pair< Point*,  Point*>(t->first , t->second) ;
		
		if(this->isVertex(t->second) && this->isVertex(t->third))
			return std::pair< Point*,  Point*>(t->second , t->third) ;
		
		if(this->isVertex(t->first) && this->isVertex(t->third))
			return std::pair< Point*,  Point*>(t->first , t->third) ;
		
		if(this->isVertex(t->first) && this->isVertex(t->fourth))
			return std::pair< Point*,  Point*>(t->first , t->fourth) ;
		
		if(this->isVertex(t->second) && this->isVertex(t->fourth))
			return std::pair< Point*,  Point*>(t->second , t->fourth) ;
		
		if(this->isVertex(t->third) && this->isVertex(t->fourth))
			return std::pair< Point*,  Point*>(t->third , t->fourth) ;
	}

	return std::pair< Point*,  Point*>(NULL, NULL) ;
}




inline std::vector<Point*> DelaunayTetrahedron::commonSurface(const DelaunayTreeItem_3D * t) const
{
	std::vector<Point *> ret ;
	if(t == this)
	{
		return ret ;
		ret.push_back(NULL) ;
		ret.push_back(NULL) ;
		ret.push_back(NULL);
	}
	
	if(t->isTetrahedron)
	{
		if(this->isVertexByID(t->first) && this->isVertexByID(t->second)&& this->isVertexByID(t->third))
		{
			ret.push_back(t->first) ;
			ret.push_back(t->second) ;
			ret.push_back(t->third) ;
		}
			
		
		if(this->isVertexByID(t->first) && this->isVertexByID(t->second)&& this->isVertexByID(t->fourth))
		{
			ret.push_back(t->first) ;
			ret.push_back(t->second) ;
			ret.push_back(t->fourth) ;
		}
		
		if(this->isVertexByID(t->first) && this->isVertexByID(t->fourth)&& this->isVertexByID(t->third))
		{
			ret.push_back(t->first) ;
			ret.push_back(t->third) ;
			ret.push_back(t->fourth) ;
		}
		
		if(this->isVertexByID(t->fourth) && this->isVertexByID(t->second)&& this->isVertexByID(t->third))
		{
			ret.push_back(t->third) ;
			ret.push_back(t->second) ;
			ret.push_back(t->fourth) ;
		}
		
		
	}
	else
	{
		Point A(*first-*second) ; 
		Point B(*third-*second) ; 
		Point C(*third-*t->first) ; 
		Point D(*third-*t->second) ; 
		Point E(*third-*t->third) ;
		double fst = std::abs(triProduct(A, B, C)) ;
		fst = std::max(fst, std::abs(triProduct(A, B, D))) ;
		fst = std::max(fst, std::abs(triProduct(A, B, E))) ;
		
		A = (*first-*third) ; 
		B = (*fourth-*third) ; 
		C = (*fourth-*t->first) ; 
		D = (*fourth-*t->second) ; 
		E = (*fourth-*t->third) ;
		
		double ftf = std::abs(triProduct(A, B, C)) ;
		ftf = std::max(ftf, std::abs(triProduct(A, B, D))) ;
		ftf = std::max(ftf, std::abs(triProduct(A, B, E))) ;
		
		A = (*first-*second) ; 
		B = (*fourth-*second) ; 
		C = (*fourth-*t->first) ; 
		D = (*fourth-*t->second) ; 
		E = (*fourth-*t->third) ;
		
		double fsf = std::abs(triProduct(A, B, C)) ;
		fsf = std::max(fsf, std::abs(triProduct(A, B, D))) ;
		fsf = std::max(fsf, std::abs(triProduct(A, B, E))) ;
		
		A = (*third-*second) ; 
		B = (*fourth-*second) ; 
		C = (*fourth-*t->first) ; 
		D = (*fourth-*t->second) ; 
		E = (*fourth-*t->third) ;
		
		double tsf = std::abs(triProduct(A, B, C)) ;
		tsf = std::max(tsf, std::abs(triProduct(A, B, D))) ;
		tsf = std::max(tsf, std::abs(triProduct(A, B, E))) ;
		
		if(   fst < ftf 
		   && fst < fsf 
		   && fst < tsf)
		{
			ret.push_back(first) ;
			ret.push_back(second) ;
			ret.push_back(third) ; 
		}
		
		if(   ftf < fst 
		   && ftf < fsf 
		   && ftf < tsf)
		{
			ret.push_back(first) ;
			ret.push_back(third) ;
			ret.push_back(fourth) ; 
		}
		
		if(   fsf < fst
		   && fsf < ftf 
		   && fsf < tsf)
		{
			ret.push_back(first) ;
			ret.push_back(second) ;
			ret.push_back(fourth) ; 
		}
		
		if(   tsf < fst 
		   && tsf < ftf 
		   && tsf < fsf)
		{
			ret.push_back(third) ;
			ret.push_back(second) ;
			ret.push_back(fourth) ; 
		}
		
// 		if(isCoplanar((*t->first), (*first), (*second),(*third)) && 
// 		   isCoplanar((*t->second), (*first), (*second),(*third)) &&
// 		   isCoplanar((*t->third), (*first), (*second),(*third)))
// 		{
// 			ret.push_back(first) ;
// 			ret.push_back(second) ;
// 			ret.push_back(third) ; 
// 		}
// 		else if(isCoplanar((*t->first), (*first), (*third),(*fourth)) && 
// 		   isCoplanar((*t->second), (*first), (*third),(*fourth)) &&
// 		   isCoplanar((*t->third), (*first), (*third),(*fourth)))
// 		{
// 			ret.push_back(first) ;
// 			ret.push_back(third) ;
// 			ret.push_back(fourth) ; 
// 		}
// 		else if(isCoplanar((*t->first), (*first), (*second),(*fourth)) && 
// 		   isCoplanar((*t->second), (*first), (*second),(*fourth)) &&
// 		   isCoplanar((*t->third), (*first), (*second),(*fourth)))
// 		{
// 			ret.push_back(first) ;
// 			ret.push_back(second) ;
// 			ret.push_back(fourth) ; 
// 		}
// 		else if(isCoplanar((*t->first), (*third), (*second),(*fourth)) && 
// 		   isCoplanar((*t->second), (*third), (*second),(*fourth)) &&
// 		   isCoplanar((*t->third), (*third), (*second),(*fourth)))
// 		{
// 			ret.push_back(third) ;
// 			ret.push_back(second) ;
// 			ret.push_back(fourth) ; 
// 		}
	}
	
	//assert(false) ;
	return ret;
}

inline std::vector<Point*> DelaunayDemiSpace::commonSurface(const DelaunayTreeItem_3D * t) const
{
	std::vector<Point *> ret ;
	if(t == this)
	{
		ret.push_back(NULL) ;
		ret.push_back(NULL) ;
		ret.push_back(NULL);
	}
	
	if(t->isTetrahedron)
	{
		return t->commonSurface(this) ;
// 		if(this->isVertex(t->first) && this->isVertex(t->second)&& this->isVertex(t->third))
// 		{
// 			ret.push_back(t->first) ;
// 			ret.push_back(t->second) ;
// 			ret.push_back(t->third) ;
// 		}
// 		
// 		
// 		if(this->isVertex(t->first) && this->isVertex(t->second)&& this->isVertex(t->fourth))
// 		{
// 			ret.push_back(t->first) ;
// 			ret.push_back(t->second) ;
// 			ret.push_back(t->fourth) ;
// 		}
// 		
// 		if(this->isVertex(t->first) && this->isVertex(t->fourth)&& this->isVertex(t->third))
// 		{
// 			ret.push_back(t->first) ;
// 			ret.push_back(t->third) ;
// 			ret.push_back(t->fourth) ;
// 		}
// 		
// 		if(this->isVertex(t->fourth) && this->isVertex(t->second)&& this->isVertex(t->third))
// 		{
// 			ret.push_back(t->first) ;
// 			ret.push_back(t->second) ;
// 			ret.push_back(t->fourth) ;
// 		}
		
		
	}

	return ret;
}

std::vector<Point*> DelaunayRoot_3D::commonSurface(const DelaunayTreeItem_3D * t) const
{
	std::vector<Point *> ret ;
	
		ret.push_back(NULL) ;
		ret.push_back(NULL) ;
		ret.push_back(NULL);
	
	
	

	return ret;
}



std::pair< Point*,  Point*> DelaunayDemiSpace::commonEdge(const DelaunayTreeItem_3D * t) 
{
	if(t->isTetrahedron)
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
	

void DelaunayDemiSpace::merge(DelaunayDemiSpace * p)
{
	if(isAlive() && p != this && p->isAlive())
	{
		if( (first == p->first && second == p->second && third == p->third) || 
		    (first == p->first && second == p->third && third == p->second) || 
		    (first == p->second && second == p->third && third == p->first)|| 
		    (first == p->second && second == p->first && third == p->third)|| 
		    (first == p->third && second == p->first && third == p->second)||
		    (first == p->third && second == p->second && third == p->first  ))
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
		if(isCoplanar(first, second,third, p->first) && 
		   isCoplanar(first, second,third, p->second) && 
		   isCoplanar(first, second,third, p->third))
		{
// 			std::cout << "merging" << std::endl ;
// 			print() ;
// 			p->print() ;
			for(size_t i = 0 ; i <  p->neighbour.size() ; i++)
			{
				this->addNeighbour(p->neighbour[i]) ;
				p->neighbour[i]->addNeighbour(this) ;
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
				this->addNeighbour(p->stepson[i]) ;
				p->stepson[i]->addNeighbour(this) ;
			}
			p->father->son.push_back(this) ;
			p->son.push_back(this) ;
			p->kill(first) ;
		}
	}
}

std::vector< Point*> DelaunayTetrahedron::nearestSurface(const Point & p)
{
// 	std::cout<<"we came to nearest surfaces"<<std::endl;
	std::vector<Point *> ret ;
	
	std::map<double, Point> cen ;
	
	Point c0(((*first) + (*second)+ (*third))/3.) ;
	cen[squareDist(c0, p)] = c0 ;
	Point c1(((*third) + (*second)+(*fourth))/3.) ;
	cen[squareDist(c1, p)] = c1 ;
	Point c2(((*third) + (*first)+ (*fourth))/3.) ;
	cen[squareDist(c2, p)] = c2 ;
	Point c3(((*second) + (*first)+ (*fourth))/3.) ;
	cen[squareDist(c2, p)] = c2 ;
	
	if(cen.begin()->second == c0)
	{
		ret.push_back(first) ;
		ret.push_back(second) ;
		ret.push_back(third) ;
	}
	if(cen.begin()->second == c1)
	{
		ret.push_back(second) ;
		ret.push_back(third) ;
		ret.push_back(fourth) ;
	}
	if(cen.begin()->second == c2)
	{
		ret.push_back(third) ;
		ret.push_back(fourth) ;
		ret.push_back(first) ;
	}
	
	if(cen.begin()->second == c3)
	{
		ret.push_back(fourth) ;
		ret.push_back(first) ;
		ret.push_back(second) ;
	}
	return ret;
}
std::vector< Point*> DelaunayDemiSpace::nearestSurface(const Point & p)

{
	std::vector<Point *> ret ;
	
	std::map<double, Point> cen ;
	
	Point c0(((*first) + (*second)+ (*third))/3.) ;
	cen[squareDist(c0, p)] = c0 ;
	Point c1(((*third) + (*second)+(*fourth))/3.) ;
	cen[squareDist(c1, p)] = c1 ;
	Point c2(((*third) + (*first)+ (*fourth))/3.) ;
	cen[squareDist(c2, p)] = c2 ;
	Point c3(((*second) + (*first)+ (*fourth))/3.) ;
	cen[squareDist(c2, p)] = c2 ;
	
	if(cen.begin()->second == c0)
	{
		ret.push_back(first) ;
		ret.push_back(second) ;
		ret.push_back(third) ;
	}
	if(cen.begin()->second == c1)
	{
		ret.push_back(second) ;
		ret.push_back(third) ;
		ret.push_back(fourth) ;
	}
	if(cen.begin()->second == c2)
	{
		ret.push_back(third) ;
		ret.push_back(fourth) ;
		ret.push_back(first) ;
	}
	
	if(cen.begin()->second == c3)
	{
		ret.push_back(fourth) ;
		ret.push_back(first) ;
		ret.push_back(second) ;
	}
	return ret;
}
std::vector< Point*> DelaunayRoot_3D::nearestSurface(const Point &p) 
{
	std::vector<Point *> ret ;
	
		ret.push_back(NULL) ;
		ret.push_back(NULL) ;
		ret.push_back(NULL) ;
	
	return ret;
}

std::pair< Point*,  Point*> DelaunayTetrahedron::nearestEdge(const Point & p)
{
	
	
	
	std::map<double, Point> cen ;
	Point c0(((*first) + (*second))/2.) ;
	cen[squareDist(c0, p)] = c0 ;
	Point c1(((*third) + (*second))/2.) ;
	cen[squareDist(c1, p)] = c1 ;
	Point c2(((*third) + (*fourth))/2.) ;
	cen[squareDist(c2, p)] = c2 ;
	Point c3(((*fourth) + (*first))/2.) ;
	cen[squareDist(c3, p)] = c3 ;
	Point c4(((*fourth) + (*second))/2.) ;
	cen[squareDist(c4, p)] = c4 ;
	Point c5(((*third) + (*first))/2.) ;
	cen[squareDist(c5, p)] = c5 ;
	
	
	
	if(cen.begin()->second == c0)
		return std::pair< Point*,  Point*>(first, second) ;
		
	if(cen.begin()->second == c1)
		return std::pair< Point*,  Point*>(second, third) ;
		
	if(cen.begin()->second == c2)
		return std::pair< Point*,  Point*>(third, fourth) ;
	
	if(cen.begin()->second == c3)
		return std::pair< Point*,  Point*>(fourth, first) ;
	
	if(cen.begin()->second == c4)
		return std::pair< Point*,  Point*>(fourth, second) ;
	
	if(cen.begin()->second == c5)
		return std::pair< Point*,  Point*>(first, third) ;
	
	return std::pair<Point *, Point *> (NULL, NULL) ;
}



bool DelaunayTetrahedron::inCircumSphere(const Point &p) const 
{
	return  this->Tetrahedron::inCircumSphere(p) ;
}

bool DelaunayTetrahedron::isNeighbour(const DelaunayTreeItem_3D * t)  const
{
	size_t cv = this->numberOfCommonVertices(t) ;
	if(t->isSpace)
	{
		return (cv == 3 );
	}
	return (cv == 3 || cv == 4);
}

bool DelaunayTetrahedron::isInNeighbourhood(const DelaunayTetrahedron * t) const
{
	return this->numberOfCommonVertices(t) > 0 ;
}

bool DelaunayTreeItem_3D::isDuplicate(const DelaunayTreeItem_3D * t) const
{

	return t!=this 
	&& (
	    (
	      this->isTetrahedron 
	      && t->isTetrahedron 
	      && this->numberOfCommonVertices(t) == 4
	     ) 
	     || 
	     (
	      this->isSpace 
	      && t->isSpace 
// 	      && this->numberOfCommonVertices(t) == 3
	      )
	     )  ;
}


void DelaunayTetrahedron::insert( std::vector<DelaunayTreeItem_3D *> & ret, Point *p,  Star_3D* s)
{
	
// 	if(!isAlive())
// 		return  ;
// 		
// 	if(isVertex(p))
// 		return  ;
// 
// 	
// 	if(son.empty())
// 	{
		for (size_t i = 0 ; i < 4 ; i++)
		{
			bool ins = neighbour[i]->isSpace || ( neighbour[i]->isTetrahedron && !neighbour[i]->inCircumSphere(*p));
			
			if (ins)
			{
				std::vector< Point*> pp = this->commonSurface(neighbour[i]) ;
				if(!isCoplanar(p, pp[0], pp[1], pp[2] ))
				{
					DelaunayTetrahedron *ss = new DelaunayTetrahedron(this, p, pp[0], pp[1],pp[2], p) ;
					son.push_back(ss) ;
					neighbour[i]->addStepson(ss) ;
					ret.push_back(ss) ;
				}
			}
		}
		s->updateNeighbourhood() ;
// 	}
	
}


void DelaunayTetrahedron::print() const
{	
	
	for(size_t i = 0 ; i < this->getBoundingPoints().size() ; i++)
	{
		std::cout << "(" << getBoundingPoint(i).x << ", " << getBoundingPoint(i).y << ", " << getBoundingPoint(i).z << ") " ;
	}
	std::cout <<  ":: "<< isAlive() << std::endl ;
}


// void DelaunayTetrahedron::displace(std::valarray<double> * eps)
// {
// 	(*eps)[first->id*2]+=(*eps)[first->id*2] ;
// 	(*eps)[first->id*2+1]+=(*eps)[first->id*2+1] ;
// }

DelaunayDemiSpace::DelaunayDemiSpace(DelaunayTreeItem_3D * father,  Point  * _one,  Point  * _two, Point  * _three, Point  * p,  Point * c) : DelaunayTreeItem_3D(father, c)
{
	second  = _two ;
	first = _one ;
	third = _three ;
	
	fourth = p ;
	dead = false ;
	vector1 =(*second)- (*first) ;
	vector2 =(*third)- (*first) ;
	pseudonormal = vector1^vector2;
	pseudonormal /= pseudonormal.norm() ;
// 	if(pseudonormal*(*first-*p) < 0)
// 	{
// 		pseudonormal = vector2^vector1 ;
// 	}
	direction = pseudonormal*(*first-*fourth) ;
	isSpace =true ;
	isTetrahedron = false ;
}

std::pair< Point*,  Point*> DelaunayDemiSpace::nearestEdge(const Point & p) 
{
	return std::pair< Point*,  Point*>(first, second) ;
}
	
bool DelaunayDemiSpace::inCircumSphere(const Point & p) const
{

	double planeConst = first->x*pseudonormal.x + first->y*pseudonormal.y + first->z*pseudonormal.z ;
	double signedDistP = p.x*pseudonormal.x + p.y*pseudonormal.y + p.z*pseudonormal.z - planeConst;
	double signedDistF = fourth->x*pseudonormal.x + fourth->y*pseudonormal.y + fourth->z*pseudonormal.z - planeConst;
	return signedDistF*signedDistP < 0   ;
	
// 	Point vector = *fourth-*p ;
// 	Point onPlane = *first-*p ;
// 	
// 	double para = vector*pseudonormal;
// 	if(std::abs(para) < 1e-8)
// 		return false ;
// 	
// 	double lambda = (onPlane*pseudonormal)/para ;
// 	
// 	return lambda > 1e-8 && lambda < 1-1e-8 ;
// 	
// 	Point v = (*first - *p) + (*second - *p) + (*third - *p) ;
// 	Point u = (*second - *fourth) +  (*first - *fourth) + (*third - *fourth);
// 	return (v*pseudonormal)*(u*pseudonormal) < 0 ;

}
	

	

bool  DelaunayDemiSpace::isNeighbour(const DelaunayTreeItem_3D * t) const 
{
	
	if(t->isTetrahedron)
	{
		return t->isNeighbour(this) ;
	}

return false ;
}

void DelaunayDemiSpace::insert( std::vector<DelaunayTreeItem_3D *> & ret, Point *p, Star_3D *s)
{

// 	if(!isAlive())
// 		return  ;
// 		
// 	if(!inCircumSphere(p))
// 		return  ;
// 	
// 	if(son.size() == 0 )
// 	{
		for(size_t i = 0 ; i <neighbour.size() ; i++)
		{
			std::vector< Point*> pp = neighbour[i]->commonSurface(this) ;

			if (!neighbour[i]->inCircumSphere(*p) )
			{
				if(!isCoplanar(p, first, second, third))
				{
					
					DelaunayTetrahedron *ss = new DelaunayTetrahedron(this, p, pp[0], pp[1] ,pp[2], p) ;

					son.push_back(ss) ;
					neighbour[i]->addStepson(ss) ;
		
					ret.push_back(ss) ;
					
					DelaunayDemiSpace *p0 = new DelaunayDemiSpace(this, pp[0], p, pp[1], pp[2],p) ;
					DelaunayDemiSpace *p1 = new DelaunayDemiSpace(this, pp[0], p, pp[2], pp[1],p) ;
					DelaunayDemiSpace *p2 = new DelaunayDemiSpace(this, pp[2], p, pp[1], pp[0],p) ;
					
					
					son.push_back(p0) ;
					son.push_back(p1) ;
					son.push_back(p2) ;
					
					ret.push_back(p0) ;
					ret.push_back(p1) ;
					ret.push_back(p2) ;
				}
			}
		}

		s->updateNeighbourhood() ;
// 	}

	return  ;
}

void DelaunayDemiSpace::print() const 
{
	std::cout << "###############(" << first->x << ", " << first->y << ", "<< first->z << ") (" <<
		second->x << ", " << second->y << ", "<< second->z<<") (" <<
		third->x << ", " << third->y << ", "<< third->z<<")"<< "X (" << fourth->x << ", " << fourth->y << ", "<< fourth->z<< ") :: " << isAlive()  << std::endl ;
}

void makeNeighbours(DelaunayTreeItem_3D *t0, DelaunayTreeItem_3D *t1 ) 
{
	if(t0 == t1 || !t0->isAlive() || !t1->isAlive())
		return ;
	if(t0->isSpace && t1->isSpace)
	   return ;
	if(!t1->isNeighbour( t0)  )
	{
		return ;
	}
	
	t0->addNeighbour(t1) ;
	t1->addNeighbour(t0) ;
} 

void updateNeighbours(std::vector<DelaunayTreeItem_3D *> * t)
{
	for(size_t i = 0 ; i < t->size() ; i++)
	{
		for(size_t j = i+1 ; j < t->size() ; j++)
		{
			makeNeighbours( (*t)[i],(*t)[j] ) ;
		}
	}
} 

DelaunayRoot_3D::DelaunayRoot_3D(Point * p0, Point * p1, Point * p2, Point * p3) : DelaunayTreeItem_3D(NULL, NULL)
{
	isSpace = false ;
	isTetrahedron =false ;
	this->father = NULL ;
	DelaunayTetrahedron *t = new DelaunayTetrahedron(this, p0, p1, p2, p3, NULL) ;
	
	DelaunayDemiSpace * pl0 = new DelaunayDemiSpace(this, p0, p1, p2,p3, NULL);
	DelaunayDemiSpace * pl1 = new DelaunayDemiSpace(this, p0, p1, p3,p2, NULL);
	DelaunayDemiSpace * pl2 = new DelaunayDemiSpace(this, p1, p2, p3,p0, NULL);
	DelaunayDemiSpace * pl3 = new DelaunayDemiSpace(this, p2, p3, p0,p1, NULL);
	
	makeNeighbours(t,pl0 ) ;
	makeNeighbours(t,pl1) ;
	makeNeighbours(t,pl2) ;
	makeNeighbours(t,pl3) ;
	makeNeighbours(pl1,pl0 ) ;
	makeNeighbours(pl2,pl1) ;
	makeNeighbours(pl0,pl2) ;
	makeNeighbours(pl3,pl0) ;
	makeNeighbours(pl3,pl1) ;
	makeNeighbours(pl3,pl2) ;
	
	son.push_back(t) ;
	son.push_back(pl0) ;
	son.push_back(pl1) ;
	son.push_back(pl2) ;
	son.push_back(pl3) ;
	
	kill(p0) ;
}
	
void DelaunayRoot_3D::print() const
{
	std::cout << "I am root !" << std::endl ;
}


DelaunayTreeItem_3D * DelaunayRoot_3D::getSon(size_t i)
{
	return son[i] ;
}
	
bool DelaunayRoot_3D::isVertex(const Point *p) const
{
	return false ;
}
	
bool DelaunayRoot_3D::inCircumSphere(const Point &p) const 
{
	return true ;
}
	
std::pair< Point*,  Point*> DelaunayRoot_3D::nearestEdge(const Point&  p)
{
	assert(false) ;
	return std::pair< Point*,  Point*>(NULL, NULL) ;
}
	
void DelaunayRoot_3D::insert(std::vector<DelaunayTreeItem_3D *>& ret, Point *p, Star_3D *s)
{
	
	for (size_t i  = 0 ;  i < son.size() ; i++)
	{
		
		std::vector<DelaunayTreeItem_3D *> temp ;
		son[i]->insert(temp, p, s) ;
		
		for(size_t j = 0 ; j< temp.size() ; j++)
		{
			ret.push_back(temp[j]) ;
		}
	}
	
	updateNeighbours(&ret) ;
}



void DelaunayRoot_3D::conflicts(std::pair<std::vector<DelaunayTetrahedron *>, std::vector<DelaunayTreeItem_3D *> > & ret, const Geometry *g)
{

	visited = true ;
	
	for (size_t i  = 0 ;  i < son.size() ; i++)
	{
		std::pair<std::vector<DelaunayTetrahedron *>, std::vector<DelaunayTreeItem_3D *> > temp  ;
		son[i]->conflicts(temp,g) ;
		ret.first.insert(ret.first.end(),temp.first.begin(), temp.first.end()) ;
		ret.second.insert(ret.second.end(),temp.second.begin(), temp.second.end()) ;
	}
	
}


void DelaunayRoot_3D::conflicts(std::pair< std::vector<DelaunayTreeItem_3D *>, std::vector<DelaunayTreeItem_3D *> > & ret, const Point *p)
{

	visited = true ;
	
	for (size_t i  = 0 ;  i < son.size() ; i++)
	{
		std::pair<std::vector<DelaunayTreeItem_3D *>, std::vector<DelaunayTreeItem_3D *> > temp  ;
		son[i]->conflicts(temp, p) ;
		ret.first.insert(ret.first.end(),temp.first.begin(), temp.first.end()) ;
		ret.second.insert(ret.second.end(),temp.second.begin(), temp.second.end()) ;
	}
	
	return  ;
}

Star_3D::Star_3D(std::vector<DelaunayTreeItem_3D *> *t, const Point *p) : creator(p)
{
	
	treeitem.insert(treeitem.end(), t->begin(), t->end()) ;
	for(size_t i = 0 ; i < t->size() ; i++)
	{
		if((*t)[i]->isTetrahedron)
		{
			this->edge.push_back((*t)[i]->first) ;
			this->edge.push_back((*t)[i]->second) ;
			this->edge.push_back((*t)[i]->third) ;
			this->edge.push_back((*t)[i]->fourth);
		}
		else if((*t)[i]->isSpace)
		{
			this->edge.push_back((*t)[i]->first) ;
			this->edge.push_back((*t)[i]->second) ;
			this->edge.push_back((*t)[i]->third);
		}
	}
	std::sort(edge.begin(), edge.end()) ;
	edge.erase(unique(edge.begin(), edge.end()), edge.end()) ;
}
	
size_t Star_3D::size()
{
	return edge.size() ;
}
	
const Point * Star_3D::getEdge(size_t i) const
{
	assert(i < edge.size()) ;
	return this->edge[i] ;
}

DelaunayTree_3D::DelaunayTree_3D(Point * p0, Point *p1, Point *p2, Point *p3)
{
	neighbourhood = false ;
	global_counter = 4;
	p0->id = 0 ; p1->id = 1 ; p2->id = 2 ; p3->id = 3;
	DelaunayRoot_3D *root = new DelaunayRoot_3D( p0, p1, p2, p3) ;
	tree.push_back(root) ;
	
	tree.push_back(root->getSon(0)) ;
	tree.push_back(root->getSon(1)) ;
	tree.push_back(root->getSon(2)) ;
	tree.push_back(root->getSon(3)) ;
	tree.push_back(root->getSon(4)) ;
	space.push_back(static_cast<DelaunayDemiSpace *>(root->getSon(1))) ;
	space.push_back(static_cast<DelaunayDemiSpace *>(root->getSon(2))) ;
	space.push_back(static_cast<DelaunayDemiSpace *>(root->getSon(3))) ;
	space.push_back(static_cast<DelaunayDemiSpace *>(root->getSon(4))) ;
}
	
DelaunayTree_3D::~DelaunayTree_3D() 
{ 
	for(size_t i = 0 ;  i < this->tree.size() ; i++)
	{
		delete this->tree[i] ;
	}
	
}


void DelaunayTree_3D::insert(Point *p)
{
	
	std::vector<DelaunayTreeItem_3D *> cons = this->conflicts(p) ;
	
	if(cons.size() == 0)
	{
		for(size_t i =0 ; i < tree.size() ; i++)
		{
			if(tree[i]->isAlive() && tree[i]->inCircumSphere(*p))
				cons.push_back(tree[i]) ;
		}
		
		for(size_t i =0 ; i< cons.size() ; i++)
		{
			std::cout << "when inserting " ;p->print() ; std::cout << " did not get " ;
			cons[i]->print() ;
			std::cout << "fatherline " << std::endl ;
			DelaunayTreeItem_3D * f = cons[i]->father ;
			if(!f->inCircumSphere(*p))
				f = cons[i]->stepfather ;
			
			while(f)
			{
				std::cout << f->inCircumSphere(*p) << "  " ;
				f->print() ;
				if(f->father && f->father->inCircumSphere(*p))
					f = f->father ;
				else if (f->stepfather && f->stepfather->inCircumSphere(*p))
					f = f->stepfather ;
				else
					f = NULL ;
			}
			
			return ;
		}
		
		
	}
	
	neighbourhood = false ;
	
	for(size_t i = 0 ; i < cons.size() ; i++)
	{
		if(cons[i]->isVertex(p)) 
		{
// 			std::cout << "vertex collision" << std::endl ;
			return ;
		}
	}
	
	p->id = this->global_counter ;
	
	this->global_counter++ ;
	
	Star_3D * s = new Star_3D(&cons, p) ;
	
	std::vector<DelaunayTreeItem_3D *> ret ;
	
	for(size_t i = 0 ; i < cons.size() ; i++)
	{
		std::vector<DelaunayTreeItem_3D *> temp ;
		cons[i]->insert(temp,p, s) ;
		ret.insert(ret.end(), temp.begin(), temp.end()) ;
	}

	s->updateNeighbourhood() ;
	
	bool weGotSpaces = false ;
	
	for(size_t j = 0 ; j< ret.size() ; j++)
	{
		if(ret[j]->isAlive() && ret[j]->isSpace )
		{
			weGotSpaces = true ;
			space.push_back(static_cast<DelaunayDemiSpace*>(ret[j])) ;
		}
	}
	
	if(weGotSpaces)
	{
		for(size_t k = 0 ; k< space.size()-1 ; k++)
		{
			for(size_t j = k+1 ; j< space.size() ; j++)
			{
				space[j]->merge(space[k]) ;
			}
		}
	}
	
	for (size_t l=0; l< ret.size(); l++)
	{	
		for (size_t k=0; k< space.size(); k++)
		{
			makeNeighbours(ret[l],space[k]);
		}
	}
	
	for(size_t i = 0 ; i < cons.size() ; i++)
	{
		cons[i]->kill(p) ;
	}

	bool correct = true ;
	
	for(size_t i = 0 ; i < ret.size() ; i++)
	{		

		if(!ret[i]->erased)
		{
			tree.push_back(ret[i]) ;
			if(ret[i]->isTetrahedron && ret[i]->neighbour.size() != 4)
			{
				std::cout << "we have " << ret[i]->neighbour.size() << " neighbours"<< std::endl ;
// 				ret[i]->print() ;
				correct = false ;
				
				Sphere sph(.1, p) ;
				std::vector<DelaunayTetrahedron *> missing = conflicts(&sph) ;
				for(size_t j = 0 ; j < missing.size() ; j++)
				{
					if(!missing[j]->isVertex(p))
					{
						
						if(missing[j]->inCircumSphere(*p))
						{
							std::cout << "found extra tet" << std::endl ;
							missing[j]->insert(ret,p,s) ;
							s->updateNeighbourhood() ;
						}
					}
				}
			}
			
		}
		else
		{
			std::valarray<Point *> nullarray(0) ;
			if(ret[i]->isTetrahedron)
				static_cast<DelaunayTetrahedron *>(ret[i])->setBoundingPoints(nullarray) ;
			delete ret[i] ;
		}
	}


	if(!correct)
		std::cout << "inconsistent state, will crash soon" <<std::endl ;	
	delete s ;
}



std::vector<DelaunayTreeItem_3D *> DelaunayTree_3D::conflicts( const Point *p) const
{
	std::vector<DelaunayTreeItem_3D *> ret  ;
	
// 	for(size_t i = 0 ; i < tree.size() ; i++)
// 	{
// 		if(tree[i]->inCircumSphere(p) && tree[i]->isAlive() && (tree[i]->isTetrahedron || tree[i]->isSpace))
// 			ret.push_back(tree[i]) ;
// 	}
// 	
// 	return ret ;
// 	
	std::pair< std::vector<DelaunayTreeItem_3D *>,std::vector<DelaunayTreeItem_3D *> > cons;
	this->tree[0]->conflicts(cons, p) ;
	
	for(size_t i = 0 ; i < space.size() ; i++)
	{

		std::pair< std::vector<DelaunayTreeItem_3D *>,std::vector<DelaunayTreeItem_3D *> > temp ;
		space[i]->conflicts(temp, p) ;
		
		cons.first.insert(cons.first.end(), temp.first.begin(),temp.first.end()) ;
		cons.second.insert(cons.second.end(), temp.second.begin(),temp.second.end()) ;
		
// 		for(size_t j = 0 ; j < space[i]->son.size() ; j++)
// 		{
// 			std::pair< std::vector<DelaunayTreeItem_3D *>, std::vector<DelaunayTreeItem_3D *> > temp;
// 			space[i]->son[j]->conflicts(temp, p) ;
// 			
// 			cons.first.insert(cons.first.end(), temp.first.begin(),temp.first.end()) ;
// 			cons.second.insert(cons.second.end(), temp.second.begin(),temp.second.end()) ;
// 		}
// 		
// 		for(size_t j = 0 ; j < space[i]->stepson.size() ; j++)
// 		{
// 			std::pair< std::vector<DelaunayTreeItem_3D *>, std::vector<DelaunayTreeItem_3D *> > temp;
// 			space[i]->stepson[j]->conflicts(temp, p) ;
// 			
// 			cons.first.insert(cons.first.end(), temp.first.begin(),temp.first.end()) ;
// 			cons.second.insert(cons.second.end(), temp.second.begin(),temp.second.end()) ;
// 		}
	}
	
	for(size_t i = 0 ; i < cons.second.size() ; i++)
	{
		cons.second[i]->clearVisited() ;
	}
	
	
	ret.insert(ret.end(), cons.first.begin(), cons.first.end()) ;
	
// 	std::cout <<"--------------------------" <<std::endl ;
// 	std::cout << "we have the following conflicts from " << std::flush ;
// 	p->print() ;std::cout << std::endl ;
// 	
// 	for(size_t i = 0 ; i < ret->size() ; i++)
// 		(*ret)[i]->print() ;
// 	std::cout <<"--------------------------" <<std::endl ;
	return ret ;
}


void DelaunayTetrahedron::refresh(const TetrahedralElement * father)
{
	this->TetrahedralElement::refresh(father) ;
// 	this->computeCenter() ;
}

std::vector<DelaunayTetrahedron *> DelaunayTree_3D::conflicts( const Geometry *g) const
{

	std::pair< std::vector<DelaunayTetrahedron *>,std::vector<DelaunayTreeItem_3D *> > cons ;
	this->tree[0]->conflicts(cons, g) ;
	
	for(size_t i = 0 ; i < space.size() ; i++)
	{
		
		if(!space[i]->visited)
		{
			std::pair< std::vector<DelaunayTetrahedron *>,std::vector<DelaunayTreeItem_3D *> > temp ;
			space[i]->conflicts(temp,g) ;
			
			cons.first.insert(cons.first.end(), temp.first.begin(),temp.first.end()) ;
			cons.second.insert(cons.second.end(), temp.second.begin(),temp.second.end()) ;
			
			
			for(size_t j = 0 ; j < space[i]->neighbour.size() ; j++)
			{
				std::pair< std::vector<DelaunayTetrahedron *>, std::vector<DelaunayTreeItem_3D *> > temp ;
				space[i]->neighbour[j]->conflicts(temp,g) ;
				
				cons.first.insert(cons.first.end(), temp.first.begin(),temp.first.end()) ;
				cons.second.insert(cons.second.end(), temp.second.begin(),temp.second.end()) ;
			}
			
		}
		
	}
	
	for(size_t i = 0 ; i < cons.second.size() ; i++)
	{
		cons.second[i]->clearVisited() ;
	}
	
	std::vector<DelaunayTetrahedron *> ret  ;
	ret.insert(ret.end(), cons.first.begin(), cons.first.end()) ;

	return ret ;
}

std::vector<DelaunayDemiSpace *>  DelaunayTree_3D::getConvexHull()
{
	std::vector<DelaunayDemiSpace *> ret ;
	
	for(size_t i = 0 ; i < space.size() ; i++)
	{
		if(space[i]->isAlive())
			ret.push_back(space[i]) ;
	}
	
	return ret ;
}

std::vector<DelaunayTetrahedron *>  DelaunayTree_3D::getTetrahedrons( bool buildNeighbourhood ) 
{
	std::vector<DelaunayTetrahedron *> ret;
	//std::cout<<tree.size();
	for(size_t i = 0 ; i < tree.size() ; i++)
	{
		
		//std::cout<<tree[i]->isAlive()<<std::endl;
// 		tree[i]->print() ;
		if(tree[i]->isAlive() && tree[i]->isTetrahedron)
		{
			ret.push_back(static_cast<DelaunayTetrahedron *>(tree[i])) ;
		}
	}
	
	if(!neighbourhood && buildNeighbourhood)
	{
		std::cout << "\r building neighbourhood... element 0/" << ret.size() << std::flush ;
		for( size_t i = 0 ; i < ret.size() ;i++)
		{
			ret[i]->neighbourhood.clear() ;
			ret[i]->clearVisited() ;
		}
		
		
		
		for( size_t i = 0 ; i < ret.size() ;i++)
		{
			if(i%100 == 0)
				std::cout << "\r building neighbourhood... element " << i <<"/" << ret.size() << std::flush ;
			
			std::vector<DelaunayTetrahedron *> tocheck ;
			std::vector<DelaunayTetrahedron *> toclean ;
			for(size_t j = 0 ; j < ret[i]->neighbour.size() ; j++)
			{
				if(ret[i]->neighbour[j]->isTetrahedron && !ret[i]->neighbour[j]->visited)
				{
					tocheck.push_back(static_cast<DelaunayTetrahedron *>(ret[i]->neighbour[j]));
					ret[i]->neighbour[j]->visited = true ;
					toclean.push_back(*tocheck.rbegin()) ;
					ret[i]->neighbourhood.push_back(*tocheck.rbegin()) ;
				}
			}
			
			while(!tocheck.empty())
			{
				std::vector<DelaunayTetrahedron *> tocheck_temp ;
				for(size_t k = 0 ; k < tocheck.size() ; k++)
				{
					for(size_t j = 0 ; j < tocheck[k]->neighbour.size() ; j++)
					{
						if(
						    tocheck[k]->neighbour[j]->isTetrahedron 
						    && !tocheck[k]->neighbour[j]->visited
						    && tocheck[k]->neighbour[j] != ret[i]
						    && static_cast<DelaunayTetrahedron *>(tocheck[k]->neighbour[j])->isInNeighbourhood(ret[i])
						  )
						{
							tocheck_temp.push_back(static_cast<DelaunayTetrahedron *>(tocheck[k]->neighbour[j]));
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

		}
		
		std::cout << " ...done" << std::endl ;
		neighbourhood = true ;
	}
	
	return ret ;
}

void DelaunayTree_3D::print() const
{
	size_t alive = 0 ;
	std::cout << "we have a total of " << tree.size() << " elements" << std::endl ;
	for(size_t i = 0 ; i < tree.size() ; i++)
	{
		if (tree[i]->isAlive())
		{
			alive++ ;
			
		}
	}
	std::cout << "of which " << alive << "are alive" << std::endl ;
	#ifdef DEBUG
	for(size_t i = 0 ; i < tree.size() ; i++)
	{
		if (tree[i]->isAlive() )
		{
// 		std::cout <<tree[i]->first;
// 		if(tree[i]->isTetrahedron)
// 		{
			tree[i]->print() ;
			
			for(size_t j = 0 ; j< tree[i]->neighbour.size() ; j++)
			{
				
// 				if(tree[i]->neighbour[j]->isTetrahedron)
// 				{
					std::cout << "\t --> " ;
					tree[i]->neighbour[j]->print() ;
// 				}
			}
		}
		//}
	}
	#endif
}

std::vector<std::vector<Matrix> > DelaunayTetrahedron::getElementaryMatrix() const 
{
	
	std::vector<std::vector<Matrix > > mother ;
	GaussPointArray gp ; 
	std::valarray<Matrix> Jinv ;
	std::vector<size_t> dofs = getDofIds() ;
	if(getEnrichmentFunctions().size() > 0 /*&& getEnrichmentFunction(0).second.getIntegrationHint().size() > 0*/ )
	{
		if(!(getEnrichmentFunctions().size() == 1 ||
		     getEnrichmentFunctions().size() == 2 ||
		     getEnrichmentFunctions().size() == 3 ||
		     getEnrichmentFunctions().size() == 4 ||
		     getEnrichmentFunctions().size() == 5 ||
		     getEnrichmentFunctions().size() == 6 ||
		     getEnrichmentFunctions().size() == 8 ||
		     getEnrichmentFunctions().size() == 9 ||
		     getEnrichmentFunctions().size() == 12 ))
			std::cout << "danger, will robinson, danger !" << std::endl ;
		std::vector<std::pair<Point, double> > gp_alternative ;
		std::vector<Point> to_add ;
		std::vector<Point> to_add_extra ;
		to_add.push_back(Point(0,0,0)) ;
		to_add.push_back(Point(0,0,1)) ;
		to_add.push_back(Point(1,0,0)) ;
		to_add.push_back(Point(0,1,0)) ;
		for(size_t i = 0 ; i <  getEnrichmentFunctions().size() ; i++)
		{
			for(size_t j = 0 ; j < getEnrichmentFunction(i).getIntegrationHint().size() ; j++)
			{
				if(getEnrichmentFunction(i).getIntegrationHint(j) != to_add[0] && 
				   getEnrichmentFunction(i).getIntegrationHint(j) != to_add[1] && 
				   getEnrichmentFunction(i).getIntegrationHint(j) != to_add[2]  &&   getEnrichmentFunction(i).getIntegrationHint(j) != to_add[3])
				{
					bool ok = true ;
					for(size_t k = 0 ; k < to_add_extra.size() ; k++)
					{
						if(getEnrichmentFunction(i).getIntegrationHint(j) == to_add_extra[k])
						{
							ok = false ;
							break ;
						}
					}
					if(ok)
						to_add_extra.push_back(getEnrichmentFunction(i).getIntegrationHint(j)) ;
				}
			}
		}
		
		to_add.insert(to_add.end(), to_add_extra.begin(),to_add_extra.end() ) ;
		DelaunayTree_3D dt(&to_add[0], &to_add[1], &to_add[2],&to_add[3]) ;
		for(size_t i = 4 ; i < to_add.size() ; i++)
		{
			dt.insert(&to_add[i]) ;
		}
		
// 		TetrahedralElement father(QUADRATIC) ;
// 		dt.addSharedNodes(1) ;
// 		dt.refresh( &father) ;
		
		std::vector<DelaunayTetrahedron *> tetra = dt.getTetrahedrons() ;
		
		Jinv.resize(tetra[0]->getGaussPoints().gaussPoints.size()*tetra.size()) ;
		
		if(moved)
		{
			for(size_t i = 0 ; i < tetra.size() ; i++)
			{
				
				Function x = tetra[i]->getXTransform() ;
				Function y = tetra[i]->getYTransform() ;
				Function z = tetra[i]->getZTransform() ;
				tetra[i]->order = QUADRATIC ;
				GaussPointArray gp_temp = tetra[i]->getGaussPoints() ;
				VirtualMachine vm ;
				
				for(size_t j = 0 ; j < gp_temp.gaussPoints.size() ; j++)
				{
					
					gp_temp.gaussPoints[j].second *= this->jacobianAtPoint(gp_temp.gaussPoints[j].first) ;
					gp_temp.gaussPoints[j].first.set(vm.eval(x, gp_temp.gaussPoints[j].first), vm.eval(y, gp_temp.gaussPoints[j].first), vm.eval(z, gp_temp.gaussPoints[j].first)) ;
					Jinv[i*gp_temp.gaussPoints.size()+j] = this->getInverseJacobianMatrix( gp_temp.gaussPoints[j].first ) ;
					
					gp_alternative.push_back(gp_temp.gaussPoints[j]) ;
				}
				
			}
		}
		else
		{
			Matrix J = this->getInverseJacobianMatrix( Point(1./4.,1./4.,1./4.)) ;
			double ja = this->jacobianAtPoint(Point(1./4.,1./4.,1./4.)) ;
			for(size_t i = 0 ; i < tetra.size() ; i++)
			{
				
				Function x = tetra[i]->getXTransform() ;
				Function y = tetra[i]->getYTransform() ;
				Function z = tetra[i]->getZTransform() ;
				GaussPointArray gp_temp = tetra[i]->getGaussPoints() ;
				VirtualMachine vm ;
				
				
				for(size_t j = 0 ; j < gp_temp.gaussPoints.size() ; j++)
				{
					gp_temp.gaussPoints[j].second *= ja ;
					gp_temp.gaussPoints[j].first.set(vm.eval(x, gp_temp.gaussPoints[j].first), vm.eval(y, gp_temp.gaussPoints[j].first), vm.eval(z, gp_temp.gaussPoints[j].first)) ;
					Jinv[i*gp_temp.gaussPoints.size()+j] = J ;
					gp_alternative.push_back(gp_temp.gaussPoints[j]) ;
				}
				
			}
		}
		
		gp.gaussPoints.resize(gp_alternative.size()) ;
		std::copy(gp_alternative.begin(), gp_alternative.end(), &gp.gaussPoints[0]);
	}
	else
	{
		if(moved)
		{
			GaussPointArray gp_alt(this->getGaussPoints()) ;
			gp.gaussPoints.resize(gp_alt.gaussPoints.size()) ;
			std::copy(&gp_alt.gaussPoints[0], &gp_alt.gaussPoints[gp_alt.gaussPoints.size()], &gp.gaussPoints[0]);
			Jinv.resize(gp_alt.gaussPoints.size()) ;
			for(size_t i = 0 ; i < gp.gaussPoints.size() ;  i++)
			{
				Jinv[i] = getInverseJacobianMatrix( gp.gaussPoints[i].first ) ;
			}
		}
		else
		{
			Matrix J = this->getInverseJacobianMatrix( Point(.25,.25,.25)) ;
			GaussPointArray gp_alt(this->getGaussPoints()) ;
			gp.gaussPoints.resize(gp_alt.gaussPoints.size()) ;
			Jinv.resize(gp_alt.gaussPoints.size()) ;
			std::copy(&gp_alt.gaussPoints[0], &gp_alt.gaussPoints[gp_alt.gaussPoints.size()], &gp.gaussPoints[0]);
			
			for(size_t i = 0 ; i < gp.gaussPoints.size() ;  i++)
			{
				Jinv[i] = J ;
			}
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


 std::vector<std::vector<Matrix> > DelaunayTetrahedron::getNonLinearElementaryMatrix() 
{
	std::vector<size_t> dofs = getDofIds() ;
	std::vector<std::vector<Matrix> > mother ;
	
	Vector dsp = this->getState().getDisplacements() ;
	this->getState().step(0, &dsp) ;
	
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
	GaussPointArray gp ; 
	std::vector<std::pair<Point, double> > gp_alternative ;
	
	
	VirtualMachine vm ;
	
	std::vector<Point> to_add ;
	std::vector<Point> to_add_extra ;
	
	to_add.push_back(Point(0,1,0)) ;
	to_add.push_back(Point(0,0,0)) ;
	to_add.push_back(Point(1,0,0)) ;
	to_add.push_back(Point(0,0,1)) ;
	for(size_t i = 0 ; i <  nonlinbehaviour->getIntegrationHints().size() ; i++)
	{
		
		if( nonlinbehaviour->getIntegrationHint(i) != to_add[0] && 
		    nonlinbehaviour->getIntegrationHint(i) != to_add[1] && 
		    nonlinbehaviour->getIntegrationHint(i) != to_add[2])
		{
			bool ok = true ;
			for(size_t k = 0 ; k < to_add_extra.size() ; k++)
			{
				if( nonlinbehaviour->getIntegrationHint(i) == to_add_extra[k])
				{
					ok = false ;
					break ;
				}
			}
			if(ok)
				to_add_extra.push_back(nonlinbehaviour->getIntegrationHint(i)) ;
		}
	}
	
	to_add.insert(to_add.end(), to_add_extra.begin(),to_add_extra.end() ) ;
	DelaunayTree_3D dt(&to_add[0], &to_add[1], &to_add[2], &to_add[4]) ;
	for(size_t i = 4 ; i < to_add.size() ; i++)
	{
		dt.insert(&to_add[i]) ;
	}
	
// 	TetrahedralElement father(QUADRATIC) ;
// 	dt.addSharedNodes(1) ;
// 	dt.refresh( &father) ;
	
	std::vector<DelaunayTetrahedron *> tetra = dt.getTetrahedrons() ;
	
	Jinv.resize(tetra.size()*tetra[0]->getGaussPoints().gaussPoints.size()) ;
	if(moved)
	{
		for(size_t i = 0 ; i < tetra.size() ; i++)
		{
			
			Function x = tetra[i]->getXTransform() ;
			Function y = tetra[i]->getYTransform() ;
			Function z = tetra[i]->getZTransform() ;
			tetra[i]->order = QUADRATIC ;
			GaussPointArray gp_temp = tetra[i]->getGaussPoints() ;
			
			if(moved)
			{
				for(size_t j = 0 ; j < gp_temp.gaussPoints.size() ; j++)
				{
					
					gp_temp.gaussPoints[j].second *= this->jacobianAtPoint(gp_temp.gaussPoints[j].first) ;
					gp_temp.gaussPoints[j].first.set(vm.eval(x, gp_temp.gaussPoints[j].first), vm.eval(y, gp_temp.gaussPoints[j].first), vm.eval(z, gp_temp.gaussPoints[j].first)) ;
					Jinv[i*gp_temp.gaussPoints.size()+j] = this->getInverseJacobianMatrix( gp_temp.gaussPoints[j].first ) ;
					
					gp_alternative.push_back(gp_temp.gaussPoints[j]) ;
				}
			}
			
		}
	}
	else
	{
		Matrix J = this->getInverseJacobianMatrix( Point(1./4., 1./4.,1./4.) ) ;
		double ja = this->jacobianAtPoint(Point(1./4., 1./4.,1./4.)) ;
		for(size_t i = 0 ; i < tetra.size() ; i++)
		{
			
			Function x = tetra[i]->getXTransform() ;
			Function y = tetra[i]->getYTransform() ;
			Function z = tetra[i]->getZTransform() ;
			
			GaussPointArray gp_temp = tetra[i]->getGaussPoints() ;
			VirtualMachine vm ;
			
			
			for(size_t j = 0 ; j < gp_temp.gaussPoints.size() ; j++)
			{
				gp_temp.gaussPoints[j].second *= ja ;
				gp_temp.gaussPoints[j].first.set(vm.eval(x, gp_temp.gaussPoints[j].first), vm.eval(y, gp_temp.gaussPoints[j].first), vm.eval(z, gp_temp.gaussPoints[j].first)) ;
				Jinv[i*gp_temp.gaussPoints.size()+j] = J ;
				
				gp_alternative.push_back(gp_temp.gaussPoints[j]) ;
			}
			
		}
	}
	
	gp.gaussPoints.resize(gp_alternative.size()) ;
	std::copy(gp_alternative.begin(), gp_alternative.end(), &gp.gaussPoints[0]);
	std::valarray<Point> gp_points(gp.gaussPoints.size()) ;
	
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
		gp_points[i] = gp.gaussPoints[i].first ;
	
	Vector displacement_state = this->getState().getDisplacements(gp_points) ;

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
		mother[i][i] = nonlinbehaviour->apply(getShapeFunction(i), getShapeFunction(i),gp, Jinv) ;
		
		for(size_t j = i+1 ; j < getShapeFunctions().size() ; j++)
		{
			mother[i][j] = nonlinbehaviour->apply(getShapeFunction(i), getShapeFunction(j),gp, Jinv) ;
			mother[j][i] = nonlinbehaviour->apply(getShapeFunction(j), getShapeFunction(i),gp, Jinv) ;
		}
		for(size_t j = 0 ; j < getEnrichmentFunctions().size() ; j++)
		{
			mother[i][j+getShapeFunctions().size()] = nonlinbehaviour->apply(getShapeFunction(i), getEnrichmentFunction(j),gp, Jinv) ;
			mother[j+getShapeFunctions().size()][i] = nonlinbehaviour->apply(getEnrichmentFunction(j), getShapeFunction(i),gp, Jinv) ;
		}
	}
	for(size_t i = 0 ; i < getEnrichmentFunctions().size() ; i++)
	{
		mother[i+getShapeFunctions().size()][i+getShapeFunctions().size()] = nonlinbehaviour->apply(getEnrichmentFunction(i), getEnrichmentFunction(i),gp, Jinv) ;
		
		for(size_t j = 0 ; j < getEnrichmentFunctions().size() ; j++)
		{
			mother[i+getShapeFunctions().size()][j+getShapeFunctions().size()] = nonlinbehaviour->apply(getEnrichmentFunction(i), getEnrichmentFunction(j),gp, Jinv) ;
			mother[j+getShapeFunctions().size()][i+getShapeFunctions().size()] = nonlinbehaviour->apply(getEnrichmentFunction(j), getEnrichmentFunction(i),gp, Jinv) ;
		}
	}
	
	return mother ;
}

Vector DelaunayTetrahedron::getNonLinearForces() const
{
	std::vector<size_t> dofs = getDofIds() ;
	Vector forces(dofs.size()*3) ;
	
	if(!this->getNonLinearBehaviour()->isActive())
	{
		forces = 0 ;
		return forces ;
	}
	
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
		Matrix J = this->getInverseJacobianMatrix(Point( 1./4.,1./4., 1./4.) ) ;
		Jinv.resize(gp.gaussPoints.size()) ;
		for(size_t i = 0 ; i < gp.gaussPoints.size() ;  i++)
		{
			Jinv[i] = J ;
		}
	}
	
	
		for(size_t i = 0 ; i < getShapeFunctions().size() ; i++)
		{

				Vector f = behaviour->getForces(this->getState(), getShapeFunction(i),gp, Jinv) ;
				
				forces[i*3]+=f[0];
				forces[i*3+1]+=f[1];
				forces[i*3+2]+=f[2];
		}

		for(size_t i = 0 ; i < getEnrichmentFunctions().size() ; i++)
		{
			Vector f = behaviour->getForces(this->getState(), getEnrichmentFunction(i),gp, Jinv) ;
			
			forces[(i+getShapeFunctions().size())*3]+=f[0];
			forces[(i+getShapeFunctions().size())*3+1]+=f[1];
			forces[(i+getShapeFunctions().size())*3+2]+=f[2];
		}

	return forces ;
}

Vector DelaunayTetrahedron::getForces() const
{
	std::vector<size_t> dofs = getDofIds() ;
	Vector forces(dofs.size()*3) ;
	
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
		Matrix J = this->getInverseJacobianMatrix(Point( 1./4.,1./4., 1./4.) ) ;
		Jinv.resize(gp.gaussPoints.size()) ;
		for(size_t i = 0 ; i < gp.gaussPoints.size() ;  i++)
		{
			Jinv[i] = J ;
		}
	}
	
	
	for(size_t i = 0 ; i < getShapeFunctions().size() ; i++)
	{

		Vector f = behaviour->getForces(this->getState(), getShapeFunction(i) ,gp, Jinv) ;
		
		forces[i*3]+=f[0];
		forces[i*3+1]+=f[1];
		forces[i*3+2]+=f[2];
	}

	for(size_t i = 0 ; i < getEnrichmentFunctions().size() ; i++)
	{
		Vector f = behaviour->getForces(this->getState(), getEnrichmentFunction(i) ,gp, Jinv) ;
		
		forces[(i+getShapeFunctions().size())*3]+=f[0];
		forces[(i+getShapeFunctions().size())*3+1]+=f[1];
		forces[(i+getShapeFunctions().size())*3+2]+=f[2];
	}


	return forces ;
}


GaussPointArray DelaunayTetrahedron::getSubTriangulatedGaussPoints() const
{
	GaussPointArray gp = getGaussPoints() ; 
	
	VirtualMachine vm ;
	
	if(getEnrichmentFunctions().size() > 0 )
	{
		std::vector<std::pair<Point, double> > gp_alternative ;
		std::vector<Point> to_add ;
		std::vector<Point> to_add_extra ;
		to_add.push_back(Point(0,1)) ;
		to_add.push_back(Point(0,0)) ;
		to_add.push_back(Point(1,0)) ;
		to_add.push_back(Point(0,0,1)) ;
		for(size_t i = 0 ; i <  getEnrichmentFunctions().size() ; i++)
		{
			for(size_t j = 0 ; j < getEnrichmentFunction(i).getIntegrationHint().size() ; j++)
			{
				if(getEnrichmentFunction(i).getIntegrationHint(j) != to_add[0] && 
				   getEnrichmentFunction(i).getIntegrationHint(j) != to_add[1] && 
				   getEnrichmentFunction(i).getIntegrationHint(j) != to_add[2])
				{
					bool ok = true ;
					for(size_t k = 0 ; k < to_add_extra.size() ; k++)
					{
						if(getEnrichmentFunction(i).getIntegrationHint(j) == to_add_extra[k])
						{
							ok = false ;
							break ;
						}
					}
					if(ok)
						to_add_extra.push_back(getEnrichmentFunction(i).getIntegrationHint(j)) ;
				}
			}
		}

		to_add.insert(to_add.end(), to_add_extra.begin(),to_add_extra.end() ) ;
		DelaunayTree_3D dt(&to_add[0], &to_add[1], &to_add[2],&to_add[3] ) ;
		for(size_t i = 4 ; i < to_add.size() ; i++)
		{
			dt.insert(&to_add[i]) ;
		}
		
// 		TetrahedralElement father(QUADRATIC) ;
// 		dt.addSharedNodes(1) ;
// 		dt.refresh( &father) ;
		
		std::vector<DelaunayTetrahedron *> tri = dt.getTetrahedrons() ;
				
		if(moved)
		{
			for(size_t i = 0 ; i < tri.size() ; i++)
			{

				Function x = tri[i]->getXTransform() ;
				Function y = tri[i]->getYTransform() ;
				Function z = tri[i]->getZTransform() ;
				tri[i]->order = QUADRATIC ;
				GaussPointArray gp_temp = tri[i]->getGaussPoints() ;

				for(size_t j = 0 ; j < gp_temp.gaussPoints.size() ; j++)
				{
	
					gp_temp.gaussPoints[j].first.set(vm.eval(x, gp_temp.gaussPoints[j].first), vm.eval(y, gp_temp.gaussPoints[j].first),  vm.eval(z, gp_temp.gaussPoints[j].first)) ;
					gp_temp.gaussPoints[j].second *= this->jacobianAtPoint(gp_temp.gaussPoints[j].first) ;
					
					
					gp_alternative.push_back(gp_temp.gaussPoints[j]) ;
				}

			}
		}
		else
		{
			double ja = this->jacobianAtPoint(Point(1./4.,1./4., 1./4.)) ;
			for(size_t i = 0 ; i < tri.size() ; i++)
			{

				Function x = tri[i]->getXTransform() ;
				Function y = tri[i]->getYTransform() ;
				Function z = tri[i]->getZTransform() ;
				
				GaussPointArray gp_temp = tri[i]->getGaussPoints() ;
				
				for(size_t j = 0 ; j < gp_temp.gaussPoints.size() ; j++)
				{
					gp_temp.gaussPoints[j].second *= ja ;
					gp_temp.gaussPoints[j].first.set(vm.eval(x, gp_temp.gaussPoints[j].first), vm.eval(y, gp_temp.gaussPoints[j].first), vm.eval(z, gp_temp.gaussPoints[j].first)) ;
					gp_alternative.push_back(gp_temp.gaussPoints[j]) ;
				}

			}
		}
		
		gp.gaussPoints.resize(gp_alternative.size()) ;
		std::copy(gp_alternative.begin(), gp_alternative.end(), &gp.gaussPoints[0]);
	}
	
	return gp ;
}




