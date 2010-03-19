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
#include <limits>

#define DEBUG
// #undef DEBUG

using namespace Mu ;

DelaunayTreeItem3D::DelaunayTreeItem3D( Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> *t, DelaunayTreeItem3D * father,  const Point * c) : stepson(0), neighbour(0), son(0)
{
	tree = t ;
	this->stepfather = -1 ;
	if(father)
		this->father = father->index;
	else
		this->father = -1 ;

	this->creator  = c ;
	this->dead() = false ;
	this->erased() = false ;
	this->visited() = false ;
	
	this->isSpace() = false ;
	this->isTetrahedron() =false ;
	this->isDeadTetrahedron() = false ;
	this->first = NULL ;
	this->second = NULL ;
	this->third = NULL ;
	this->fourth = NULL ;
	index = 0 ;
	if(tree)
	{
		index = tree->tree.size() ;
		tree->tree.push_back(this) ;
	}
}

bool & DelaunayTreeItem3D::dead()
{
	return m_dead ;
}

bool & DelaunayTreeItem3D::isSpace()
{
	return m_isSpace ;
}

bool & DelaunayTreeItem3D::isTetrahedron()
{
	return m_isTetrahedron ;
}

bool & DelaunayTreeItem3D::isDeadTetrahedron()
{
	return m_isDeadTetrahedron ;
}

bool & DelaunayTreeItem3D::visited()
{
	return m_visited ;
}

bool & DelaunayTreeItem3D::erased()
{
	return m_erased ;
}

bool DelaunayTreeItem3D::dead() const
{
	return m_dead ;
}

bool DelaunayTreeItem3D::isSpace() const
{
	return m_isSpace ;
}

bool DelaunayTreeItem3D::isTetrahedron() const
{
	return m_isTetrahedron ;
}

bool DelaunayTreeItem3D::isDeadTetrahedron() const
{
	return m_isDeadTetrahedron ;
}

bool DelaunayTreeItem3D::visited() const
{
	return m_visited ;
}

bool DelaunayTreeItem3D::erased() const
{
	return m_erased ;
}

bool DelaunayTreeItem3D::isAlive() const
{
	return !m_dead ;
}


DelaunayTreeItem3D::~DelaunayTreeItem3D()
{
// 	tree->tree[index] = NULL ;
}
	

size_t DelaunayTreeItem3D::numberOfCommonVertices(const DelaunayTreeItem3D * s) const
{
	if(this->isTetrahedron() && s->isTetrahedron())
	{
		Point * inter[4] ;
		Point ** e = std::set_intersection(&first, &first+4, &s->first, &s->first+4, &inter[0])  ;
		return e-&inter[0] ;
	}
	if(this->isTetrahedron() && s->isSpace())
	{
		size_t ret = 0 ;
		
		if(isCoplanar(first,s->first, s->second, s->third))
			ret++ ;
		if(isCoplanar(second,s->first, s->second, s->third ))
			ret++ ;
		if(isCoplanar(third,s->first, s->second, s->third ))
			ret++ ;
		if(ret != 3 && isCoplanar(fourth,s->first, s->second, s->third ))
			ret++ ;
		
		assert(ret < 4) ;
		
		return ret ;
	}
	if(this->isSpace() && s->isSpace())
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
	
	if(this->isSpace() && s->isTetrahedron())
	{
		size_t ret = 0 ;
		
		if(isCoplanar(s->first,this->first, this->second, this->third  ))
			ret++ ;
		if(isCoplanar(s->second,this->first, this->second, this->third ))
			ret++ ;
		if(isCoplanar(s->third,this->first, this->second, this->third ))
			ret++ ;
		if(ret != 3 && isCoplanar(s->fourth,this->first, this->second, this->third ))
			ret++ ;
		
		assert(ret < 4) ;
		
		return ret ;
	}
	
	return 0 ;
}


void Star3D::updateNeighbourhood()
{
	int count = 0 ;
	int soncount = 0 ;
	for(std::vector<DelaunayTreeItem3D *>::const_iterator i = treeitem.begin() ; i != treeitem.end() ;++i)
	{
		count += (*i)->son.size() + (*i)->stepson.size() + (*i)->neighbour.size() ;
		soncount += (*i)->son.size() ;
	}
	std::valarray<DelaunayTreeItem3D *> items(count) ;
	count = 0 ;
	for(std::vector<DelaunayTreeItem3D *>::const_iterator i = treeitem.begin() ; i != treeitem.end() ;++i)
	{	
		for(size_t j = 0 ; j < (*i)->son.size() ; j++)
		{
			items[count++] = (*i)->getSon(j) ;
		}
		
	}
	for(std::vector<DelaunayTreeItem3D *>::const_iterator i = treeitem.begin() ; i != treeitem.end() ;++i)
	{	
		for(size_t j = 0 ; j < (*i)->stepson.size() ; j++)
		{
			if((*i)->getStepson(j)->isAlive()&& ((*i)->getStepson(j)->onCircumSphere(*creator) || (*i)->getStepson(j)->inCircumSphere(*creator)))
				items[count++] = (*i)->getStepson(j) ;
		}
		for(size_t j = 0 ; j < (*i)->neighbour.size() ; j++)
		{
			if((*i)->getNeighbour(j)->isAlive() && ((*i)->getNeighbour(j)->onCircumSphere(*creator) || (*i)->getNeighbour(j)->inCircumSphere(*creator))) ;
				items[count++] = (*i)->getNeighbour(j) ;
		}
		
	}
	
	std::sort(&items[soncount], &items[count]) ;
	DelaunayTreeItem3D ** e = std::unique(&items[0], &items[count]) ;
	size_t end =  e-&items[0] ;
	
	if(!items.size())
		return ;
// 	std::cout << end << std::endl ;
// 	std::vector<DelaunayTreeItem3D *> & tree = items[0]->tree->tree ;
// 	int dcount = 0 ;
// 	for(size_t i = 0 ; i < end; ++i)
// 	{
// 		if(items[i]->isAlive())
// 		{
// 			for(size_t j = i+1 ; j < end ; ++j)
// 			{
// 				DelaunayTreeItem3D * ii = items[i] ;
// 				DelaunayTreeItem3D * jj = items[j] ;
// 				if(jj->isAlive() && ii->isDuplicate(jj) )
// 				{
// 					dcount++ ;
// 					for(size_t k = 0 ; k < jj->neighbour.size() ; k++)
// 					{
// 						if(jj->neighbour.size() != 4 && ii->neighbour.size() != 4)
// 							makeNeighbours(jj->getNeighbour(k), ii) ;
// 					}
// 					
// 					jj->kill(creator) ;
// 					jj->erased() = true ;
// 					std::valarray<unsigned int> newSons(tree[jj->father]->son.size()-1);
// 					bool found = 0 ;
// 					for(size_t k = 0 ; k < tree[jj->father]->son.size() ; k++)
// 					{
// 						if(tree[jj->father]->son[k] == jj->index)
// 							found = true ;
// 						else
// 							newSons[k-found] = tree[jj->father]->son[k] ;
// 						
// 					}
// 					tree[jj->father]->son.resize(newSons.size()-1);
// 					tree[jj->father]->son = newSons;
// 					if(jj->stepfather != -1)
// 					{
// 						std::valarray<unsigned int> newStepsons(tree[jj->stepfather]->stepson.size()-1);
// 						bool found = 0 ;
// 						for(size_t k = 0 ; k < tree[jj->stepfather]->stepson.size() ; k++)
// 						{
// 							if(tree[jj->stepfather]->stepson[k] == jj->index)
// 								found = true ;
// 							else
// 								newStepsons[k-found] = tree[jj->stepfather]->stepson[k] ;
// 							
// 						}
// 						tree[jj->stepfather]->stepson.resize(newStepsons.size());
// 						tree[jj->stepfather]->stepson =  newStepsons;
// 					}
// 				}
// 			}
// 		}
// 	}
// 	std::cout << "\n" << dcount << std::endl ;
	for(DelaunayTreeItem3D ** i = &items[0] ; i != e+1 && i != &items[count] ; ++i)
	{
		if(!(*i)->isSpace())
		{
			DelaunayTreeItem3D * ii = *i ;
			size_t ins = ii->neighbour.size() ;
			if(ins != 4 )
			{
				for(DelaunayTreeItem3D ** j = i+1 ; j != e+1 && j != &items[count] ; ++j)
				{
					
					DelaunayTreeItem3D * jj = *j ;

					size_t jns = jj->neighbour.size() ;
					if(!jj->isSpace() && ii->numberOfCommonVertices(jj) == 3 && jns != 4)
					{
						
						if(std::find(&ii->neighbour[0], &ii->neighbour[ins],jj->index) ==  &ii->neighbour[ins])
						{
							std::valarray<unsigned int>  newneighbours(ins+1) ;
							std::copy(&ii->neighbour[0], &ii->neighbour[ins], &newneighbours[0]) ;
							newneighbours[ins] = jj->index ;
							ii->neighbour.resize(ins+1) ;
							ii->neighbour = newneighbours ;
						}
						
						if(std::find(&jj->neighbour[0], &jj->neighbour[jns],ii->index) ==  &jj->neighbour[jns])
						{
							std::valarray<unsigned int>  newneighbours(jns+1) ;
							std::copy(&jj->neighbour[0], &jj->neighbour[jns], &newneighbours[0]) ;
							newneighbours[jns] = ii->index ;
							jj->neighbour.resize(jns+1) ;
							jj->neighbour = newneighbours ;
						}
						ins = ii->neighbour.size() ;
						if(ins == 4)
							break ;
					}
				}
			}
		}
	}

}


void DelaunayTree3D::addSharedNodes(size_t nodes_per_side, size_t time_planes, double timestep, const TetrahedralElement * father)
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
		
		tri[i]->visited() = true ;
		
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
							if(tri[i]->getNeighbourhood(j)->visited())
							{
								for(size_t k = 0 ; 
								    k < tri[i]->getNeighbourhood(j)->getBoundingPoints().size() ;
								    k++)
								{
									if(tri[i]->getNeighbourhood(j)->getBoundingPoint(k) == proto)
									{
										foundPoint = &tri[i]->getNeighbourhood(j)->getBoundingPoint(k) ;
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
							additionalPoints.push_back(new Point(proto)) ;
							newPoints[positions[current]+plane*nodes_per_plane]  = additionalPoints.back() ;
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

void DelaunayTree3D::refresh(const TetrahedralElement *father)
{
/*	std::cout << "debug - assemble - setelementbehaviour - refresh" << std::endl ;
	std::cout << this->tree.size() << std::endl ;
	std::vector<int> pb ;
	for(size_t i = 0 ; i < this->tree.size() ; i++)
	{
		if(this->tree[i] == NULL)
			pb.push_back(i) ;
	}
	std::cout << pb.size() << " NULL elements in tree" << std::endl ;
	for(size_t i = pb.size() ; i > 0 ; i--)
	{
		this->tree.erase(this->tree.begin()+pb[i-1]) ;
	}
	std::cout << this->tree.size() << std::endl ;
*/
	for(size_t i = 0 ; i < this->tree.size() ; i++)
	{
//		if(this->tree[i] == NULL)
//			std::cout << "NULL element in tree" << std::endl ;
		if(this->tree[i]->isAlive() && this->tree[i]->isTetrahedron())
		{
			static_cast<DelaunayTetrahedron *>(this->tree[i])->refresh(father) ;
			
		}
	}
}

size_t DelaunayTree3D::numPoints() const
{
	return this->global_counter ;
}

void DelaunayTreeItem3D::conflicts(std::pair<std::vector<DelaunayTetrahedron *>, std::vector<DelaunayTreeItem3D *> > & ret, const Geometry *g)
{

	if(visited())
	{
		return  ;
	}
	visited() = true ;
	
	ret.second.push_back(this) ;
	
	
	if(
	    !inCircumSphere(g->getCenter()) 
	    && 
	    (
	      (!g->in(*first) 
	      && !g->in(*second) 
	      && !g->in(*third) 
	      && !g->in(*fourth) 
	       && isTetrahedron() )
	      || 
	      ( 
	        g->in(*fourth) 
	        && isSpace() 
	      )
	    ) 
	  )
	{
		return  ;
	}
	
	for (size_t i  = 0 ;  i < stepson.size() ; i++)
	{
		if( !getStepson(i)->visited() 
		    && 
		    !(
		          !getStepson(i)->inCircumSphere(g->getCenter()) 
		          && 
		          (
		            (!g->in(*getStepson(i)->first) 
		            && !g->in(*getStepson(i)->second) 
		            && !g->in(*getStepson(i)->third) 
		            && !g->in(*getStepson(i)->fourth) 
		            && getStepson(i)->isTetrahedron() 
		            ) 
		            || 
		            ( 
		              g->in(*getStepson(i)->fourth) 
		              && getStepson(i)->isSpace() 
		            )
		          )
		        )
		  )
		{
			getStepson(i)->conflicts(ret, g) ;
		}
	}
	
	for (size_t i  = 0 ;  i < son.size() ; i++)
	{
		if( 
		    !getSon(i)->visited() 
		    && 
		    !(
		       !getSon(i)->inCircumSphere(g->getCenter()) 
		       && 
		       (
		         (!g->in(*getSon(i)->first) 
		          && !g->in(*getSon(i)->second) 
		          && !g->in(*getSon(i)->third) 
		          && !g->in(*getSon(i)->fourth) 
		          && getSon(i)->isTetrahedron()
		         ) 
		         || 
		         ( 
		           g->in(*getSon(i)->fourth) 
		           && getSon(i)->isSpace() 
		         )
		       )
		     )
		  )
		{
			getSon(i)->conflicts(ret, g) ;
		}
	}
	
	if(isAlive() && isTetrahedron())
	{
		ret.first.push_back(static_cast<DelaunayTetrahedron *>(this)) ;
	}
	
	for (size_t i  = 0 ;  i < neighbour.size() ; i++)
	{
		if( 
		    !getNeighbour(i)->visited() 
		    && 
		    !(
		       !getNeighbour(i)->inCircumSphere(g->getCenter()) 
		       && 
		       (
		         (
		           !g->in(*getNeighbour(i)->first) 
		           && !g->in(*getNeighbour(i)->second) 
		           && !g->in(*getNeighbour(i)->third) 
		           && !g->in(*getNeighbour(i)->fourth) 
		           && getNeighbour(i)->isTetrahedron()  
		         )
		         || 
		         ( 
		           g->in(*getNeighbour(i)->fourth) 
		           && getNeighbour(i)->isSpace() 
		         )
		       )
		     )
		  )
		{
			getNeighbour(i)->conflicts(ret, g) ;
		}
	}
}


void DelaunayTreeItem3D::flatConflicts(std::vector<DelaunayTreeItem3D *> & toTest,  std::pair<std::vector<DelaunayTetrahedron *>, std::vector<DelaunayTreeItem3D *> > & ret, const Geometry *g)
{
	if(visited())
	{
		return ;
	}
	visited() = true ;
	
	ret.second.push_back(this) ;
	
	
	if(
	    !inCircumSphere(g->getCenter()) 
	    && 
	    (
	      (!g->in(*first) 
	      && !g->in(*second) 
	      && !g->in(*third) 
	      && !g->in(*fourth) 
	       && isTetrahedron() )
	      || 
	      ( 
	        g->in(*fourth) 
	        && isSpace() 
	      )
	    ) 
	  )
	{
		return ;
	}
	
	for (size_t i  = 0 ;  i < stepson.size() ; i++)
	{
		if( !getStepson(i)->visited() 
		    && 
		    !(
		          !getStepson(i)->inCircumSphere(g->getCenter()) 
		          && 
		          (
		            (!g->in(*getStepson(i)->first) 
		            && !g->in(*getStepson(i)->second) 
		            && !g->in(*getStepson(i)->third) 
		            && !g->in(*getStepson(i)->fourth) 
		            && getStepson(i)->isTetrahedron() 
		            ) 
		            || 
		            ( 
		              g->in(*getStepson(i)->fourth) 
		              && getStepson(i)->isSpace() 
		            )
		          )
		        )
		  )
		{
			toTest.push_back(getStepson(i)) ;
		}
	}
	
	for (size_t i  = 0 ;  i < son.size() ; i++)
	{
		if( 
		    !getSon(i)->visited() 
		    && 
		    !(
		       !getSon(i)->inCircumSphere(g->getCenter()) 
		       && 
		       (
		         (!g->in(*getSon(i)->first) 
		          && !g->in(*getSon(i)->second) 
		          && !g->in(*getSon(i)->third) 
		          && !g->in(*getSon(i)->fourth) 
		          && getSon(i)->isTetrahedron()
		         ) 
		         || 
		         ( 
		           g->in(*getSon(i)->fourth) 
		           && getSon(i)->isSpace() 
		         )
		       )
		     )
		  )
		{
			toTest.push_back(getSon(i)) ;
		}
	}
	
	if(isAlive() && isTetrahedron())
	{
		ret.first.push_back(static_cast<DelaunayTetrahedron *>(this)) ;
	}
	
	for (size_t i  = 0 ;  i < neighbour.size() ; i++)
	{
		if( 
		    !getNeighbour(i)->visited() 
		    && 
		    !(
		       !getNeighbour(i)->inCircumSphere(g->getCenter()) 
		       && 
		       (
		         (
		           !g->in(*getNeighbour(i)->first) 
		           && !g->in(*getNeighbour(i)->second) 
		           && !g->in(*getNeighbour(i)->third) 
		           && !g->in(*getNeighbour(i)->fourth) 
		           && getNeighbour(i)->isTetrahedron()  
		         )
		         || 
		         ( 
		           g->in(*getNeighbour(i)->fourth) 
		           && getNeighbour(i)->isSpace() 
		         )
		       )
		     )
		  )
		{
			toTest.push_back(getNeighbour(i)) ;
		}
	}
}


 
 
 void DelaunayTreeItem3D::conflicts(std::pair<std::vector<DelaunayTreeItem3D *>, std::vector<DelaunayTreeItem3D *> > & ret, const Point *p)
{

	if(visited())
	{
		return  ;
	}
	visited() = true ;
	
	ret.second.push_back(this) ;
	for (size_t i  = 0 ;  i < stepson.size() ; i++)
	{
		if( (!getStepson(i)->visited() && getStepson(i)->inCircumSphere(*p)) ) 
		{
			getStepson(i)->conflicts(ret,p) ;
		}
	}
	
	for (size_t i  = 0 ;  i < son.size() ; i++)
	{

		if( (!getSon(i)->visited() && getSon(i)->inCircumSphere(*p)))
		{
			getSon(i)->conflicts(ret,p) ;
		}
	}
	
	if(!inCircumSphere(*p))
	{
		return  ;
	}
	
	for (size_t i  = 0 ;  i < neighbour.size() ; i++)
	{
		
		if( (!getNeighbour(i)->visited() && getNeighbour(i)->inCircumSphere(*p)))
		{
			getNeighbour(i)->conflicts(ret,p) ;
		}
	}
	


	if(isAlive())
	{
		ret.first.push_back(this) ;
	}
	

}

 void DelaunayTreeItem3D::flatConflicts( std::vector<DelaunayTreeItem3D *> & toTest, std::pair<std::vector<DelaunayTreeItem3D *>, std::vector<DelaunayTreeItem3D *> > & ret, const Point *p)
{
	if(visited())
	{
		return ;
	}
	visited() = true ;
	
	ret.second.push_back(this) ;
	for (size_t i  = 0 ;  i < stepson.size() ; i++)
	{
		if( (!getStepson(i)->visited() && getStepson(i)->inCircumSphere(*p)) ) 
		{
			toTest.push_back(getStepson(i)) ;
		}
	}
	
	for (size_t i  = 0 ;  i < son.size() ; i++)
	{

		if( (!getSon(i)->visited() && getSon(i)->inCircumSphere(*p)))
		{
			toTest.push_back(getSon(i)) ;
		}
	}
	
	if(!inCircumSphere(*p))
	{
		return ;
	}
	
	for (size_t i  = 0 ;  i < neighbour.size() ; i++)
	{
		
		if( (!getNeighbour(i)->visited() && getNeighbour(i)->inCircumSphere(*p)))
		{
			toTest.push_back(getNeighbour(i)) ;
		}
	}
	


	if(isAlive())
	{
		ret.first.push_back(this) ;
	}
}


void DelaunayTreeItem3D::removeNeighbour(DelaunayTreeItem3D * t)
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

void DelaunayTetrahedron::removeNeighbourhood(DelaunayTetrahedron * t)
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
	
void DelaunayTreeItem3D::addNeighbour(DelaunayTreeItem3D * t)
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

void DelaunayTetrahedron::addNeighbourhood(DelaunayTetrahedron * t)
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

DelaunayTetrahedron * DelaunayTetrahedron::getNeighbourhood(size_t i) const
{
	return static_cast<DelaunayTetrahedron *>(tree->tree[neighbourhood[i]]) ;
}

DelaunayTreeItem3D * DelaunayTreeItem3D::getNeighbour(size_t i) const
{
	return tree->tree[neighbour[i]] ;
}

DelaunayTreeItem3D * DelaunayTreeItem3D::getSon(size_t i) const
{
	return tree->tree[son[i]] ;
}

DelaunayTreeItem3D * DelaunayTreeItem3D::getStepson(size_t i) const
{
	return tree->tree[stepson[i]] ;
}


void DelaunayTreeItem3D::erase(const Point * p)
{
	dead() = true ;
	killer  = p ;
}

void DelaunayTreeItem3D::kill(const Point * p)
{
	dead() = true ;
	killer  = p ;
	
	for(size_t i = 0 ; i < this->neighbour.size() ; i++)
	{
		this->getNeighbour(i)->removeNeighbour(this) ;
	}
}

void DelaunayTreeItem3D::addStepson(DelaunayTreeItem3D * s)
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

void DelaunayTreeItem3D::addSon(DelaunayTreeItem3D * s)
{
	std::valarray<unsigned int>  newson(son.size()+1) ;
	std::copy(&son[0], &son[son.size()], &newson[0]) ;
	newson[son.size()] = s->index ;
	son.resize(son.size()+1) ;
	son = newson ;
}
	
void DelaunayTreeItem3D::removeSon(DelaunayTreeItem3D * t)
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

void DelaunayTreeItem3D::removeStepson(DelaunayTreeItem3D * t)
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
		if(newstepson.size() == 0)
			t->stepfather = -1 ;

	}
}

void DelaunayTreeItem3D::setStepfather(DelaunayTreeItem3D * s)
{
	if(stepfather != -1)
		tree->tree[stepfather]->removeStepson(this) ;
	stepfather = s->index ;
	addNeighbour(s) ;
}

void DelaunayTreeItem3D::setFather(DelaunayTreeItem3D * s)
{
	if(father != -1)
		tree->tree[father]->removeSon(this) ;
	father = s->index ;
}
	
void DelaunayTreeItem3D::clearVisited()
{
	visited() = false ;
}

DelaunayTetrahedron::~DelaunayTetrahedron()
{
}

void DelaunayTetrahedron::kill(const Point * p)
{
	this->DelaunayTreeItem3D::kill(p) ;
	getBoundingPoints().resize(0) ;
	getInPoints().resize(0) ;
}

DelaunayTetrahedron::DelaunayTetrahedron(Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> *t, DelaunayTreeItem3D * father,  Point *p0,  Point *p1, Point *p2, Point *p3,  Point * c) : TetrahedralElement(p0, p1, p2, p3), DelaunayTreeItem3D(t, father, c)
{
	first = &getBoundingPoint(0) ;
	second = &getBoundingPoint(1) ;
	third = &getBoundingPoint(2) ;
	fourth = &getBoundingPoint(3) ;
	
	std::sort(&first, &first+4) ;
	
	assert(in(this->getCenter())) ;
	
	isSpace() = false ;
	isTetrahedron() = true ;
	isDeadTetrahedron() = false ;
	visited() =false ;
	assert(first->id > -1) ;
	assert(second->id > -1) ;
	assert(third->id > -1) ;
	assert(fourth->id> -1);
}

DelaunayTetrahedron::DelaunayTetrahedron(Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t, DelaunayTreeItem3D * father,  Point *p0,  Point *p1, Point *p2, Point *p3,  Point *p4,  Point *p5, Point *p6, Point *p7,Point * c) : TetrahedralElement(p0, p1, p2, p3, p4, p5, p6, p7), DelaunayTreeItem3D(t, father, c)
{
	
	first = &getBoundingPoint(0) ;
	second = &getBoundingPoint(1) ;
	third = &getBoundingPoint(2) ;
	fourth = &getBoundingPoint(3) ;
	std::sort(&first, &first+4) ;
	assert(in(this->getCenter())) ;
	
	isSpace() = false ;
	isTetrahedron() = true ;
	isDeadTetrahedron() = false ;
	visited() =false ;
	
	assert(first->id > -1) ;
	assert(second->id > -1) ;
	assert(third->id > -1) ;
	assert(fourth->id> -1);
	
}

DelaunayTetrahedron::DelaunayTetrahedron() : DelaunayTreeItem3D(NULL, NULL, NULL)
{
	first = &getBoundingPoint(0) ;
	second = &getBoundingPoint(1) ;
	third = &getBoundingPoint(2) ;
	fourth = &getBoundingPoint(3) ;
	assert(in(this->getCenter())) ;
	
	isSpace() = false ;
	isTetrahedron() = true ;
	isDeadTetrahedron() = false ;
	visited() =false ;
	index = 0 ;
}

DelaunayDemiSpace::~DelaunayDemiSpace()
{
}
	
bool DelaunayTetrahedron::isVertex(const Point * p) const
{
	return squareDist3D(*p, *first) < 128.*POINT_TOLERANCE*POINT_TOLERANCE
			|| squareDist3D(*p, *second) < 128.*POINT_TOLERANCE*POINT_TOLERANCE 
			|| squareDist3D(*p, *third) < 128.*POINT_TOLERANCE*POINT_TOLERANCE 
			|| squareDist3D(*p, *fourth) < 128.*POINT_TOLERANCE*POINT_TOLERANCE ;
	return (dist(p,first) < 0.0000001 || dist(p,second) < 0.0000001 || dist(p,third) < 0.0000001||dist(p,fourth) < 0.0000001) ;
}

bool DelaunayDemiSpace::isVertex(const Point * p) const
{
	return squareDist3D(*p, *first) < 128.*POINT_TOLERANCE*POINT_TOLERANCE
			|| squareDist3D(*p, *second) < 128.*POINT_TOLERANCE*POINT_TOLERANCE 
			|| squareDist3D(*p, *third) < 128.*POINT_TOLERANCE*POINT_TOLERANCE ;
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

DelaunayDeadTetrahedron::DelaunayDeadTetrahedron( DelaunayTetrahedron * parent) : DelaunayTreeItem3D(*parent),
x(parent->getCircumCenter()->x), y(parent->getCircumCenter()->y), z(parent->getCircumCenter()->z),
radius(parent->getRadius())
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

	isTetrahedron() = true ;
	isSpace() = false ;
	isDeadTetrahedron() = true ;
	visited() =false ;
	
	first = parent->first ;
	second = parent->second ;
	third = parent->third ;
	fourth = parent->fourth ;
	dead() = true ;
}
	
DelaunayDeadTetrahedron::~DelaunayDeadTetrahedron() { } ;

	
std::vector< Point*> DelaunayDeadTetrahedron::commonSurface( DelaunayTreeItem3D * t)
{
	std::vector<Point *> ret ;
	if(t == this)
	{
		return ret ;
		ret.push_back(NULL) ;
		ret.push_back(NULL) ;
		ret.push_back(NULL);
	}
	
	if(t->isTetrahedron())
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
		
		if(   fst <= ftf 
		   && fst <= fsf 
		   && fst <= tsf)
		{
			ret.push_back(first) ;
			ret.push_back(second) ;
			ret.push_back(third) ; 
		}
		
		else if(   ftf <= fst 
		   && ftf <= fsf 
		   && ftf <= tsf)
		{
			ret.push_back(first) ;
			ret.push_back(third) ;
			ret.push_back(fourth) ; 
		}
		
		else if(   fsf <= fst
		   && fsf <= ftf 
		   && fsf <= tsf)
		{
			ret.push_back(first) ;
			ret.push_back(second) ;
			ret.push_back(fourth) ; 
		}
		
		else
		{
			ret.push_back(third) ;
			ret.push_back(second) ;
			ret.push_back(fourth) ; 
		}

	}
	
	return ret;
}


bool DelaunayDeadTetrahedron::inCircumSphere(const Point & p) const
{
	if(p.x > x+1.0001*radius)
		return false ;
	if(p.x < x-1.0001*radius)
		return false ;
	if(p.y > y+1.0001*radius)
		return false ;
	if(p.y < y-1.0001*radius)
		return false ;
	if(p.z > z+1.0001*radius)
		return false ;
	if(p.z < z-1.0001*radius)
		return false ;
	
	double d = squareDist3D(Point(x, y, z), p) ;
	return  d/(radius*radius)-1 < POINT_TOLERANCE/radius ;
}

bool DelaunayDeadTetrahedron::onCircumSphere(const Point & p) const
{
	if(p.x > x+1.0001*radius)
		return false ;
	if(p.x < x-1.0001*radius)
		return false ;
	if(p.y > y+1.0001*radius)
		return false ;
	if(p.y < y-1.0001*radius)
		return false ;
	if(p.z > z+1.0001*radius)
		return false ;
	if(p.z < z-1.0001*radius)
		return false ;
	
	double d = squareDist3D(Point(x, y, z), p) ;
	return  std::abs(d/(radius*radius)-1) < POINT_TOLERANCE/radius ;
}

bool DelaunayDeadTetrahedron::isNeighbour( const DelaunayTreeItem3D * t) const
{
	size_t cv = this->numberOfCommonVertices(t) ;

	return (cv == 3 );

}

bool DelaunayDeadTetrahedron::isVertex(const Point *p) const
{
		return squareDist3D(*p, *first) < 128.*POINT_TOLERANCE*POINT_TOLERANCE
			|| squareDist3D(*p, *second) < 128.*POINT_TOLERANCE*POINT_TOLERANCE 
			|| squareDist3D(*p, *third) < 128.*POINT_TOLERANCE*POINT_TOLERANCE 
			|| squareDist3D(*p, *fourth) < 128.*POINT_TOLERANCE*POINT_TOLERANCE ;
}

bool DelaunayDeadTetrahedron::isVertexByID(const Point *p)
{
	return p == first || p == second || p == third || p == fourth ;
}

bool DelaunayDeadTetrahedron::in( const Point & p) const
{
	return Tetrahedron(third, first, second, fourth).in(p) ;
}

void DelaunayDeadTetrahedron::print() const
{

	std::cout << "(" << first->x << ", " << first->y << ", " << first->z << ") " ;
	std::cout << "(" << second->x << ", " << second->y << ", " << second->z << ") " ;
	std::cout << "(" << third->x << ", " << third->y << ", " << third->z << ") " ;
	std::cout << "(" << fourth->x << ", " << fourth->y << ", " << fourth->z << ") " ;
	std::cout <<  ":: "<< isAlive() << std::endl ;
}
	
std::vector<Point*> DelaunayTetrahedron::commonSurface( DelaunayTreeItem3D * t)
{
	std::vector<Point *> ret ;
	if(t == this)
	{
		return ret ;
		ret.push_back(NULL) ;
		ret.push_back(NULL) ;
		ret.push_back(NULL);
	}
	
	if(t->isTetrahedron())
	{
		std::vector<Point*> inter(4) ;
		std::vector<Point*>::iterator e = std::set_intersection(&first, &first+4, &t->first, &t->first+4, inter.begin())  ;
		inter.erase(e, inter.end()) ;
		return inter ;
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
		
		if(   fst <= ftf 
		   && fst <= fsf 
		   && fst <= tsf)
		{
			ret.push_back(first) ;
			ret.push_back(second) ;
			ret.push_back(third) ; 
		}
		
		else if(   ftf <= fst 
		   && ftf <= fsf 
		   && ftf <= tsf)
		{
			ret.push_back(first) ;
			ret.push_back(third) ;
			ret.push_back(fourth) ; 
		}
		
		else if(   fsf <= fst
		   && fsf <= ftf 
		   && fsf <= tsf)
		{
			ret.push_back(first) ;
			ret.push_back(second) ;
			ret.push_back(fourth) ; 
		}
		
		else 
		{
			ret.push_back(third) ;
			ret.push_back(second) ;
			ret.push_back(fourth) ; 
		}

	}
	
	return ret;
}

inline std::vector<Point*> DelaunayDemiSpace::commonSurface( DelaunayTreeItem3D * t)
{
	std::vector<Point *> ret ;
	if(t == this)
	{
		ret.push_back(NULL) ;
		ret.push_back(NULL) ;
		ret.push_back(NULL);
	}
	
	if(t->isTetrahedron())
	{
		return t->commonSurface(this) ;
	}

	return ret;
}

std::vector<Point*> DelaunayRoot3D::commonSurface( DelaunayTreeItem3D * t) 
{
	std::vector<Point *> ret ;
	
	ret.push_back(NULL) ;
	ret.push_back(NULL) ;
	ret.push_back(NULL);

	return ret;
}


void DelaunayDemiSpace::merge(DelaunayDemiSpace * p)
{
	if(isAlive() && p != this && p->isAlive())
	{
		if( std::binary_search(&first, &first+3, p->first) && 
		    std::binary_search(&first, &first+3, p->second) && 
		    std::binary_search(&first, &first+3, p->third) )
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
		if(isCoplanar(p->first,first, second,third ) && 
		   isCoplanar(p->second,first, second,third ) && 
		   isCoplanar( p->third,first, second,third))
		{
			for(size_t i = 0 ; i <  p->neighbour.size() ; i++)
			{
				this->addNeighbour(p->getNeighbour(i)) ;
				p->getNeighbour(i)->addNeighbour(this) ;
			}
			for(size_t i = 0 ; i <  p->son.size() ; i++)
			{
				p->getSon(i)->setFather(this) ;
				addSon(p->getSon(i)) ;
			}
			for(size_t i = 0 ; i <  p->stepson.size() ; i++)
			{
				std::valarray<unsigned int>  newstepson(stepson.size()+1) ;
				if(stepson.size() != 0)
					std::copy(&stepson[0], &stepson[stepson.size()], &newstepson[0]) ;
				newstepson[stepson.size()] = p->stepson[i] ;
				stepson.resize(newstepson.size()) ;
				stepson = newstepson ;
				p->getStepson(i)->stepfather = index ;

				this->addNeighbour(p->getStepson(i)) ;
				p->getStepson(i)->addNeighbour(this) ;
			}
			tree->tree[p->father]->addSon(this) ;
			p->addSon(this) ;
			p->kill(first) ;
		}
	}
}

bool DelaunayTetrahedron::inCircumSphere(const Point &p) const 
{
	if(p.x > circumCenter.x+1.0001*radius)
		return false ;
	if(p.x < circumCenter.x-1.0001*radius)
		return false ;
	if(p.y > circumCenter.y+1.0001*radius)
		return false ;
	if(p.y < circumCenter.y-1.0001*radius)
		return false ;
	if(p.z > circumCenter.z+1.0001*radius)
		return false ;
	if(p.z < circumCenter.z-1.0001*radius)
		return false ;
	
	double d = dist(circumCenter, p) ;
	return  d/radius-1 < POINT_TOLERANCE ;
}

bool DelaunayTetrahedron::onCircumSphere(const Point &p) const 
{
	if(p.x > circumCenter.x+1.0001*radius)
		return false ;
	if(p.x < circumCenter.x-1.0001*radius)
		return false ;
	if(p.y > circumCenter.y+1.0001*radius)
		return false ;
	if(p.y < circumCenter.y-1.0001*radius)
		return false ;
	if(p.z > circumCenter.z+1.0001*radius)
		return false ;
	if(p.z < circumCenter.z-1.0001*radius)
		return false ;
	
	double d = dist(*getCircumCenter(), p) ;
	return  std::abs(d/radius-1) < POINT_TOLERANCE ;
}


bool DelaunayTetrahedron::isNeighbour( const DelaunayTreeItem3D * t) const
{

// 	for(size_t i = 0 ; i < neighbour.size() ; i++)
// 		if(neighbour[i] == t->index) 
// 			return true ;
// 	for(size_t i = 0 ; i < t->neighbour.size() ; i++)
// 		if(t->neighbour[i] == index)
// 			return true ;

	size_t cv = numberOfCommonVertices(t) ;

	return (cv == 3 );

}

bool DelaunayTreeItem3D::isDuplicate( const DelaunayTreeItem3D * t) const
{

	if(t == this)
		return false;
	
	if(!isTetrahedron() || isDeadTetrahedron())
		return false ;
	
	if(!t->isTetrahedron() || t->isDeadTetrahedron())
		return false ;

	if(squareDist3D( static_cast<const DelaunayTetrahedron *>(t)->getCenter(), static_cast<const DelaunayTetrahedron *>(this)->getCenter() ) > POINT_TOLERANCE)
		return false ;
	
	if(squareDist3D( static_cast<const DelaunayTetrahedron *>(t)->getCircumCenter() ,static_cast<const DelaunayTetrahedron *>(this)->getCircumCenter() )  > POINT_TOLERANCE)
		return false ;
	
	if(!std::binary_search(&first, &first+4, t->first))
	{
		return false ;
	}
	if(!std::binary_search(&first, &first+4, t->second))
	{
		return false ;
	}
	if(!std::binary_search(&first, &first+4, t->third))
	{
		return false ;
	}
	if(!std::binary_search(&first, &first+4, t->fourth))
	{
		return false ;
	}

	return true ;
}


void DelaunayTetrahedron::insert( std::vector<DelaunayTreeItem3D *> & ret, Point *p,  Star3D* s)
{
	if (visited())
		return ;

	visited() = true ;

	for (size_t i = 0 ; i < 4 ; i++)
	{
		bool ins =  !getNeighbour(i)->visited() && (!getNeighbour(i)->inCircumSphere(*p) || getNeighbour(i)->onCircumSphere(*p)) ;
		
		if (ins)
		{
			
			std::vector< Point*> pp = this->commonSurface(getNeighbour(i)) ;

			if(Tetrahedron(*p, *pp[0], *pp[1], *pp[2]).volume() > POINT_TOLERANCE)
			{
				DelaunayTetrahedron *ss = new DelaunayTetrahedron(this->tree, this, p, pp[0], pp[1],pp[2], p) ;

				addSon(ss) ;
				getNeighbour(i)->addStepson(ss) ;
				ret.push_back(ss) ;
				
			}
		}
	}
// 	s->updateNeighbourhood() ;

	
}

void DelaunayTetrahedron::print() const
{	
	
	std::cout << "(" << first->x << ", " << first->y << ", " << first->z << ", " << first->id << ") " ;
	std::cout << "(" << second->x << ", " << second->y << ", " << second->z << ", " << second->id << ") " ;
	std::cout << "(" << third->x << ", " << third->y << ", " << third->z << ", " << third->id<< ") " ;
	std::cout << "(" << fourth->x << ", " << fourth->y << ", " << fourth->z << ", " << fourth->id <<") " ;
	std::cout <<  ":: "<< isAlive() << std::endl ;
}

// void DelaunayTetrahedron::displace(std::valarray<double> * eps)
// {
// 	(*eps)[first->id*2]+=(*eps)[first->id*2] ;
// 	(*eps)[first->id*2+1]+=(*eps)[first->id*2+1] ;
// }

DelaunayDemiSpace::DelaunayDemiSpace(Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t, DelaunayTreeItem3D * father,  Point  * _one,  Point  * _two, Point  * _three, Point  * p,  Point * c) : DelaunayTreeItem3D(t, father, c)
{
	second  = _two ;
	first = _one ;
	third = _three ;
	std::sort(&first, &first+3) ;
	fourth = p ;
	dead() = false ;
	pseudonormal = ((*second)- (*first))^((*third)- (*first));
	pseudonormal /= pseudonormal.norm() ;
// 	if(pseudonormal*(*first-*p) < 0)
// 	{
// 		pseudonormal = vector2^vector1 ;
// 	}
	isSpace() =true ;
	isTetrahedron() = false ;
	isDeadTetrahedron() = false;
	visited() =false ;
	

}

	
bool DelaunayDemiSpace::inCircumSphere(const Point & p) const
{
	double planeConst = first->x*pseudonormal.x + first->y*pseudonormal.y + first->z*pseudonormal.z ;
	double signedDistP = p.x*pseudonormal.x + p.y*pseudonormal.y + p.z*pseudonormal.z - planeConst;
	double signedDistF = fourth->x*pseudonormal.x + fourth->y*pseudonormal.y + fourth->z*pseudonormal.z - planeConst;
	return signedDistF*signedDistP < POINT_TOLERANCE*POINT_TOLERANCE || isCoplanar(&p, first, second, third);

}

bool DelaunayDemiSpace::onCircumSphere(const Point & p) const
{
	double planeConst = first->x*pseudonormal.x + first->y*pseudonormal.y + first->z*pseudonormal.z ;
	double signedDistP = p.x*pseudonormal.x + p.y*pseudonormal.y + p.z*pseudonormal.z - planeConst;
	double signedDistF = fourth->x*pseudonormal.x + fourth->y*pseudonormal.y + fourth->z*pseudonormal.z - planeConst;
	return std::abs(signedDistP) < POINT_TOLERANCE ;

}

	

bool  DelaunayDemiSpace::isNeighbour(const DelaunayTreeItem3D * t)  const 
{
// 	for(size_t i = 0 ; i < neighbour.size() ; i++)
// 		if(neighbour[i] == t->index) 
// 			return true ;
// 	for(size_t i = 0 ; i < t->neighbour.size() ; i++)
// 		if(t->neighbour[i] == index)
// 			return true ;
// 	if(std::find(&neighbour[0], &neighbour[neighbour.size()],t->index) !=  &neighbour[neighbour.size()])
// 		return true ;
// 	if(std::find(&t->neighbour[0], &t->neighbour[t->neighbour.size()],index) !=  &t->neighbour[t->neighbour.size()])
// 		return true ;
	
	if(t->isTetrahedron())
	{
		return numberOfCommonVertices(t) == 3;
	}

	return numberOfCommonVertices(t) >=2 ;
}

void DelaunayDemiSpace::insert( std::vector<DelaunayTreeItem3D *> & ret, Point *p, Star3D *s)
{

	if (visited())
		return ;
	
	visited() = true ;

	for(size_t i = 0 ; i <neighbour.size() ; i++)
	{
		std::vector< Point*> pp = getNeighbour(i)->commonSurface(this) ;

		if (!getNeighbour(i)->visited() && (!getNeighbour(i)->inCircumSphere(*p) || getNeighbour(i)->onCircumSphere(*p)))
		{
			if(!isCoplanar(p, pp[0], pp[1] ,pp[2]))
			{

				DelaunayTetrahedron *ss = new DelaunayTetrahedron(tree, this, p, pp[0], pp[1] ,pp[2], p) ;

				addSon(ss) ;
				getNeighbour(i)->addStepson(ss) ;
	
				ret.push_back(ss) ;
				
				DelaunayDemiSpace *p0 = new DelaunayDemiSpace(tree, this, pp[0], p, pp[1], pp[2],p) ;
				DelaunayDemiSpace *p1 = new DelaunayDemiSpace(tree, this, pp[0], p, pp[2], pp[1],p) ;
				DelaunayDemiSpace *p2 = new DelaunayDemiSpace(tree, this, pp[2], p, pp[1], pp[0],p) ;

				addSon(p0) ;
				addSon(p1) ;
				addSon(p2) ;
				
				ret.push_back(p0) ;
				ret.push_back(p1) ;
				ret.push_back(p2) ;
			}
		}
	}

// 	s->updateNeighbourhood() ;


	return  ;
}

void DelaunayDemiSpace::print() const 
{
	std::cout << "###############(" << first->x << ", " << first->y << ", "<< first->z << ") (" <<
		second->x << ", " << second->y << ", "<< second->z<<") (" <<
		third->x << ", " << third->y << ", "<< third->z<<")"<< "X (" << fourth->x << ", " << fourth->y << ", "<< fourth->z<< ") :: " << isAlive()  << std::endl ;
}

void makeNeighbours(DelaunayTreeItem3D *t0, DelaunayTreeItem3D *t1 ) 
{
	if(t0 == t1 || !t0->isAlive() || !t1->isAlive())
		return ;
	if(t0->isSpace() && t1->isSpace())
	   return ;
	if(!t1->isNeighbour( t0)  )
	{
		return ;
	}
	if(std::find(&t0->neighbour[0], &t0->neighbour[t0->neighbour.size()],t1->index) ==  &t0->neighbour[t0->neighbour.size()])
	{
		std::valarray<unsigned int>  newneighbours(t0->neighbour.size()+1) ;
		std::copy(&t0->neighbour[0], &t0->neighbour[t0->neighbour.size()], &newneighbours[0]) ;
		newneighbours[t0->neighbour.size()] = t1->index ;
		t0->neighbour.resize(t0->neighbour.size()+1) ;
		t0->neighbour = newneighbours ;
	}
	
	if(std::find(&t1->neighbour[0], &t1->neighbour[t1->neighbour.size()],t0->index) ==  &t1->neighbour[t1->neighbour.size()])
	{
		std::valarray<unsigned int>  newneighbours(t1->neighbour.size()+1) ;
		std::copy(&t1->neighbour[0], &t1->neighbour[t1->neighbour.size()], &newneighbours[0]) ;
		newneighbours[t1->neighbour.size()] = t0->index ;
		t1->neighbour.resize(t1->neighbour.size()+1) ;
		t1->neighbour = newneighbours ;
	}

} 

void updateNeighbours(std::vector<DelaunayTreeItem3D *> * t)
{
	for(size_t i = 0 ; i < t->size() ; i++)
	{
		for(size_t j = i+1 ; j < t->size() ; j++)
		{
			makeNeighbours( (*t)[i],(*t)[j] ) ;
		}
	}
} 

DelaunayRoot3D::DelaunayRoot3D(Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> *t, Point * p0, Point * p1, Point * p2, Point * p3) : DelaunayTreeItem3D(t, NULL, NULL)
{
	isSpace() = false ;
	isTetrahedron() = false ;
	isDeadTetrahedron() =false ;
	this->father = -1 ;
	DelaunayTetrahedron *tet = new DelaunayTetrahedron(t, this, p0, p1, p2, p3, NULL) ;
	DelaunayDemiSpace * pl0 = new DelaunayDemiSpace(t, this, p0, p1, p2,p3, NULL);
	DelaunayDemiSpace * pl1 = new DelaunayDemiSpace(t, this, p0, p1, p3,p2, NULL);
	DelaunayDemiSpace * pl2 = new DelaunayDemiSpace(t, this, p1, p2, p3,p0, NULL);
	DelaunayDemiSpace * pl3 = new DelaunayDemiSpace(t, this, p2, p3, p0,p1, NULL);
	
	makeNeighbours(tet,pl0 ) ;
	makeNeighbours(tet,pl1) ;
	makeNeighbours(tet,pl2) ;
	makeNeighbours(tet,pl3) ;
	makeNeighbours(pl1,pl0 ) ;
	makeNeighbours(pl2,pl1) ;
	makeNeighbours(pl0,pl2) ;
	makeNeighbours(pl3,pl0) ;
	makeNeighbours(pl3,pl1) ;
	makeNeighbours(pl3,pl2) ;
	
	addSon(tet) ;
	addSon(pl0) ;
	addSon(pl1) ;
	addSon(pl2) ;
	addSon(pl3) ;
	
	kill(p0) ;
}
	
void DelaunayRoot3D::print() const
{
	std::cout << "I am root !" << std::endl ;
}
	
bool DelaunayRoot3D::isVertex(const Point *p) const
{
	return false ;
}
	
bool DelaunayRoot3D::inCircumSphere(const Point &p) const 
{
	return true ;
}

bool DelaunayRoot3D::onCircumSphere(const Point &p) const 
{
	return true ;
}
	
	
void DelaunayRoot3D::insert(std::vector<DelaunayTreeItem3D *>& ret, Point *p, Star3D *s)
{
	
	for (size_t i  = 0 ;  i < son.size() ; i++)
	{
		
		std::vector<DelaunayTreeItem3D *> temp ;
		getSon(i)->insert(temp, p, s) ;
		
		for(size_t j = 0 ; j< temp.size() ; j++)
		{
			ret.push_back(temp[j]) ;
		}
	}
	
	updateNeighbours(&ret) ;
}



void DelaunayRoot3D::conflicts(std::pair<std::vector<DelaunayTetrahedron *>, std::vector<DelaunayTreeItem3D *> > & ret, const Geometry *g)
{
	visited() = true ;
	
	for (size_t i  = 0 ;  i < son.size() ; i++)
	{
		std::pair<std::vector<DelaunayTetrahedron *>, std::vector<DelaunayTreeItem3D *> > temp  ;
		std::vector<DelaunayTreeItem3D *> toTest ;
		getSon(i)->flatConflicts(toTest, temp,g) ;
		while(!toTest.empty())
		{
			std::vector<DelaunayTreeItem3D *> tempToTest ;
			for(size_t j  = 0 ;  j < toTest.size() ; j++)
			{
				toTest[j]->flatConflicts(tempToTest, temp,g) ;
			}
			
			toTest = tempToTest ;
		}
		
		ret.first.insert(ret.first.end(),temp.first.begin(), temp.first.end()) ;
		ret.second.insert(ret.second.end(),temp.second.begin(), temp.second.end()) ;
	}
	
// 	visited() = true ;
// 	
// 	for (size_t i  = 0 ;  i < son.size() ; i++)
// 	{
// 		std::pair<std::vector<DelaunayTetrahedron *>, std::vector<DelaunayTreeItem3D *> > temp  ;
// 		getSon(i)->conflicts(temp,g) ;
// 		ret.first.insert(ret.first.end(),temp.first.begin(), temp.first.end()) ;
// 		ret.second.insert(ret.second.end(),temp.second.begin(), temp.second.end()) ;
// 	}
	
}


void DelaunayRoot3D::conflicts(std::pair< std::vector<DelaunayTreeItem3D *>, std::vector<DelaunayTreeItem3D *> > & ret, const Point *p)
{

// 	visited() = true ;
// 	
// 	for (size_t i  = 0 ;  i < son.size() ; i++)
// 	{
// 		std::pair<std::vector<DelaunayTreeItem3D *>, std::vector<DelaunayTreeItem3D *> > temp  ;
// 		std::vector<DelaunayTreeItem3D *> toTest ;
// 		getSon(i)->flatConflicts(toTest, temp,p) ;
// 		while(!toTest.empty())
// 		{
// 			std::vector<DelaunayTreeItem3D *> tempToTest ;
// 			for(size_t j  = 0 ;  j < toTest.size() ; j++)
// 			{
// 				toTest[j]->flatConflicts(tempToTest, temp,p) ;
// 			}
// 			
// 			toTest = tempToTest ;
// 		}
// 		
// 		ret.first.insert(ret.first.end(),temp.first.begin(), temp.first.end()) ;
// 		ret.second.insert(ret.second.end(),temp.second.begin(), temp.second.end()) ;
// 	}
		
	visited() = true ;
	
	for (size_t i  = 0 ;  i < son.size() ; i++)
	{
		std::pair<std::vector<DelaunayTreeItem3D *>, std::vector<DelaunayTreeItem3D *> > temp  ;
		getSon(i)->conflicts(temp, p) ;
		ret.first.insert(ret.first.end(),temp.first.begin(), temp.first.end()) ;
		ret.second.insert(ret.second.end(),temp.second.begin(), temp.second.end()) ;
	}
	
	return  ;
}

Star3D::Star3D(std::vector<DelaunayTreeItem3D *> *t, const Point *p) :  treeitem(*t), creator(p)
{
	for(size_t i = 0 ; i < t->size() ; i++)
	{
		if((*t)[i]->isTetrahedron())
		{
			this->edge.push_back((*t)[i]->first) ;
			this->edge.push_back((*t)[i]->second) ;
			this->edge.push_back((*t)[i]->third) ;
			this->edge.push_back((*t)[i]->fourth);
		}
		else if((*t)[i]->isSpace())
		{
			this->edge.push_back((*t)[i]->first) ;
			this->edge.push_back((*t)[i]->second) ;
			this->edge.push_back((*t)[i]->third);
		}
	}
	std::sort(edge.begin(), edge.end()) ;
	edge.erase(std::unique(edge.begin(), edge.end()), edge.end()) ;
}
	
size_t Star3D::size()
{
	return edge.size() ;
}
	
const Point * Star3D::getEdge(size_t i) const
{
	assert(i < edge.size()) ;
	return this->edge[i] ;
}

DelaunayTree3D::DelaunayTree3D(Point * p0, Point *p1, Point *p2, Point *p3)
{
	neighbourhood = false ;
	global_counter = 4;
	p0->id = 0 ; p1->id = 1 ; p2->id = 2 ; p3->id = 3;
	DelaunayRoot3D *root = new DelaunayRoot3D(this, p0, p1, p2, p3) ;
	
	space.push_back(static_cast<DelaunayDemiSpace *>(root->getSon(1))) ;
	space.push_back(static_cast<DelaunayDemiSpace *>(root->getSon(2))) ;
	space.push_back(static_cast<DelaunayDemiSpace *>(root->getSon(3))) ;
	space.push_back(static_cast<DelaunayDemiSpace *>(root->getSon(4))) ;
	
}
	
DelaunayTree3D::~DelaunayTree3D() 
{ 

	std::valarray<Point *> nularray(0) ;
	for(size_t i = 0 ;  i < this->tree.size() ; i++)
	{

		if(this->tree[i] && this->tree[i]->isTetrahedron() && !this->tree[i]->isDeadTetrahedron() )
		{
			DelaunayTetrahedron * t = dynamic_cast<DelaunayTetrahedron *>(tree[i]) ;
			
			t->setBoundingPoints(nularray) ;
		}

	}
	
	for(size_t i = 0 ;  i < this->tree.size() ; i++)
	{
		delete this->tree[i] ;
	}
	
	for(size_t i = 0 ; i < additionalPoints.size() ; i++)
		delete additionalPoints[i] ;
	
}


void DelaunayTree3D::insert(Point *p)
{
	std::vector<DelaunayTreeItem3D *> cons = this->conflicts(p) ;
	neighbourhood = false ;

	for(size_t i = 0 ; i < cons.size() ; i++)
	{
		if(cons[i]->isVertex(p)) 
		{
			std::cout << "vertex collision" << std::endl ;
			return ;
		}
	}
	
	p->id = this->global_counter ;
	
	this->global_counter++ ;
	
	Star3D * s = new Star3D(&cons, p) ;
	
	std::vector<DelaunayTreeItem3D *> ret ;
	
	for(size_t i = 0 ; i < cons.size() ; i++)
	{
		if(!cons[i]->onCircumSphere(*p))
		{
			std::vector<DelaunayTreeItem3D *> temp ;
			cons[i]->insert(temp,p, s) ;
			ret.insert(ret.end(), temp.begin(), temp.end()) ;
		}
	}
	s->updateNeighbourhood() ;
	bool weGotSpaces = false ;
	
	for(size_t i = 0 ; i < ret.size() ; i++)
	{
		ret[i]->clearVisited() ;
		
	}
	
	for(size_t j = 0 ; j< ret.size() ; j++)
	{
		if(ret[j]->isAlive() && ret[j]->isSpace() )
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
		if(!cons[i]->onCircumSphere(*p))
			cons[i]->kill(p) ;
	}

	bool correct = true ;
	
	for(size_t i = 0 ; i < ret.size() ; i++)
	{		

		if(!ret[i]->erased() && ((ret[i]->isAlive() && ret[i]->isTetrahedron()) || ret[i]->isSpace()) )
		{
// #ifndef NDEBUG
			if(ret[i]->isTetrahedron() && ret[i]->neighbour.size() != 4)
			{
				
				std::cout << "we have " << ret[i]->neighbour.size() << " neighbours"<< std::endl ;
				p->print() ;
				tree[ret[i]->father]->print() ;
				ret[i]->print() ;
				for(size_t k = 0 ; k < ret[i]->neighbour.size() ;  k++)
				{
					std::cout << "   --> " << std::flush ;ret[i]->getNeighbour(k)->print() ;
				}
				correct = false ;
			}
// #endif
		}
		else
		{

			std::valarray<Point *> nullarray(0) ;
			if(ret[i]->isTetrahedron())
				static_cast<DelaunayTetrahedron *>(ret[i])->setBoundingPoints(nullarray) ;
			tree[ret[i]->index] = NULL ;
			delete ret[i] ;
		}
	}
	
	for(std::vector<DelaunayDemiSpace *>::iterator i = space.begin() ; i != space.end() ; i++)
	{
		if(!(*i)->isAlive())
		{
			space.erase(i) ;
			i-- ;
		}
	}

	if(!correct)
		std::cout << "inconsistent state, will crash soon" <<std::endl ;	
	
	for(size_t i = 0 ; i < cons.size() ; i++)
	{
	
			if(!cons[i]->isAlive() &&cons[i]->isTetrahedron() && !cons[i]->isDeadTetrahedron())
			{
				DelaunayDeadTetrahedron* dt = new DelaunayDeadTetrahedron(static_cast<DelaunayTetrahedron *>(cons[i])) ;
				dt->clearVisited() ;
				tree[cons[i]->index] = dt ;
				delete cons[i] ;
			}
	}
// // 	
	delete s ;
// 	print() ;
}

Star3D::~Star3D()
{
	for( size_t i = 0 ; i < cleanup.size() ; i++)
		cleanup[i]->visited() = false ;
}

std::vector<DelaunayTreeItem3D *> DelaunayTree3D::conflicts( const Point *p) const
{
	std::vector<DelaunayTreeItem3D *> ret  ;
	
	std::pair< std::vector<DelaunayTreeItem3D *>,std::vector<DelaunayTreeItem3D *> > cons;
	this->tree[0]->conflicts(cons, p) ;
	
	if(cons.first.empty() || cons.second.empty())
	{
		for(size_t i = 0 ; i < space.size() ; i++)
		{
	
			std::pair< std::vector<DelaunayTreeItem3D *>,std::vector<DelaunayTreeItem3D *> > temp ;
			space[i]->conflicts(temp, p) ;
			
			cons.first.insert(cons.first.end(), temp.first.begin(),temp.first.end()) ;
			cons.second.insert(cons.second.end(), temp.second.begin(),temp.second.end()) ;
		}
	}
	
	for(size_t i = 0 ; i < cons.second.size() ; i++)
	{
		cons.second[i]->clearVisited() ;
	}
	
	ret.insert(ret.end(), cons.first.begin(), cons.first.end()) ;

	return ret ;
}

void DelaunayTree3D::purge()
{
	int originalSize = tree.size() ;
	std::cerr << "purging..." << std::flush ;
	std::vector<DelaunayTetrahedron *> tets = getTetrahedrons() ;
	std::valarray<int> indexCorrespondance(tree.size());
	std::vector<DelaunayTreeItem3D *> newTree ;
	std::vector<DelaunayDemiSpace *> newSpace ;
	newTree.push_back(tree[0]) ;
	indexCorrespondance[tree[0]->index] = 0 ;
	int counter = 1 ;
	for(size_t i = 0 ; i < tets.size() ; i++)
	{
		if(i%1000 == 0)
			std::cerr << "\r getting elements element " << i << "/" << tets.size()+space.size() << std::flush ;
		
		newTree.push_back(tets[i]) ;
		indexCorrespondance[tets[i]->index] = counter++ ;
	}
	for(size_t i = 0 ; i < space.size() ; i++)
	{
		if(i%1000 == 0)
			std::cerr << "\r getting elements element " << i+tets.size() << "/" << tets.size()+space.size() << std::flush ;
		if(space[i]->isAlive())
		{
			newTree.push_back(space[i]) ;
			indexCorrespondance[space[i]->index] = counter++ ;
			newSpace.push_back(space[i]) ;
		}
	}
	std::cerr << std::endl ;
	
	tree[0]->neighbour.resize(0) ;
	tree[0]->son.resize(0) ;
	tree[0]->stepson.resize(0) ;
	std::valarray<unsigned int>  newson(newTree.size()-1) ;
	for(size_t i = 1 ; i < newTree.size() ; i++)
	{
		if(i%1000 == 0)
			std::cerr << "\r reconstructing element " << i << "/" << newTree.size() << std::flush ;
		for(size_t j = 0 ; j < newTree[i]->neighbour.size() ; j++)
		{
			newTree[i]->neighbour[j] = indexCorrespondance[newTree[i]->neighbour[j]] ;
		}
		
		for(size_t j = 0 ; j < newTree[i]->son.size() ; j++)
		{
			newTree[i]->son[j] = indexCorrespondance[newTree[i]->son[j]] ;
		}
		
		for(size_t j = 0 ; j < newTree[i]->stepson.size() ; j++)
		{
			newTree[i]->stepson[j] = indexCorrespondance[newTree[i]->stepson[j]] ;
		}
		
		newTree[i]->father = 0 ;
		newTree[i]->stepfather = -1 ;
		newTree[i]->index = i ;
		newson[i-1] = i ;
	}
	
	tree[0]->son.resize(newson.size()) ;
	std::copy(&tree[0]->son[0], &tree[0]->son[tree[0]->son.size()], &newson[0]) ;
	
	std::cerr << std::endl ;
	for(size_t i = 1 ; i < tree.size() ; i++)
	{
		if(!tree[i]->isAlive())
			delete tree[i] ;
	}
	tree = newTree ;
	space = newSpace ;
	std::cerr << " ...done. " << tree.size() << " elements kept from "<< originalSize << "."<< std::endl ;
}


void DelaunayTetrahedron::refresh(const TetrahedralElement * father)
{
	this->TetrahedralElement::refresh(father) ;
// 	this->computeCenter() ;
}

std::vector<DelaunayTetrahedron *> DelaunayTree3D::conflicts( const Geometry *g) const
{
	std::pair< std::vector<DelaunayTetrahedron *>,std::vector<DelaunayTreeItem3D *> > cons ;
	this->tree[0]->conflicts(cons, g) ;
	
	for(size_t i = 0 ; i < space.size() ; i++)
	{
		
		if(!space[i]->visited())
		{
			std::pair< std::vector<DelaunayTetrahedron *>,std::vector<DelaunayTreeItem3D *> > temp ;
			space[i]->conflicts(temp,g) ;
			
			cons.first.insert(cons.first.end(), temp.first.begin(),temp.first.end()) ;
			cons.second.insert(cons.second.end(), temp.second.begin(),temp.second.end()) ;
			
			
			for(size_t j = 0 ; j < space[i]->neighbour.size() ; j++)
			{
				std::pair< std::vector<DelaunayTetrahedron *>, std::vector<DelaunayTreeItem3D *> > temp ;
				space[i]->getNeighbour(j)->conflicts(temp,g) ;
				
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

std::vector<DelaunayDemiSpace *>  DelaunayTree3D::getConvexHull()
{
	std::vector<DelaunayDemiSpace *> ret ;
	
	for(size_t i = 0 ; i < space.size() ; i++)
	{
		if(space[i]->isAlive())
			ret.push_back(space[i]) ;
	}
	
	return ret ;
}

std::vector<DelaunayTetrahedron *>  DelaunayTree3D::getTetrahedrons( bool buildNeighbourhood ) 
{
	std::vector<DelaunayTetrahedron *> ret;
	//std::cout<<tree.size();
	for(size_t i = 0 ; i < tree.size() ; i++)
	{

		if(tree[i] && tree[i]->isAlive() && tree[i]->isTetrahedron())
		{
			ret.push_back(static_cast<DelaunayTetrahedron *>(tree[i])) ;
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
			if(i%10000 == 0)
				std::cerr << "\r building neighbourhood... element " << i <<"/" << ret.size() << std::flush ;
			
			std::vector<DelaunayTetrahedron *> tocheck ;
			std::vector<DelaunayTetrahedron *> toclean ;
for(size_t j = 0 ; j < ret[i]->neighbour.size() ; j++)
			{
				if(ret[i]->getNeighbour(j)->isTetrahedron() && !ret[i]->getNeighbour(j)->visited())
				{
					tocheck.push_back(static_cast<DelaunayTetrahedron *>(ret[i]->getNeighbour(j)));
					ret[i]->getNeighbour(j)->visited() = true ;
					toclean.push_back(tocheck.back()) ;
					ret[i]->addNeighbourhood(tocheck.back()) ;
				}
			}
			
			for(size_t k = 0 ; k < tocheck.size() ; k++)
			{
				if(tocheck[k]->numberOfCommonVertices(ret[i]) > 0)
					ret[i]->addNeighbourhood(tocheck[k]) ;
			}
			
			while(!tocheck.empty())
			{
				std::vector<DelaunayTetrahedron *> tocheck_temp ;
				for(size_t k = 0 ; k < tocheck.size() ; k++)
				{
					for(size_t j = 0 ; j < tocheck[k]->neighbour.size() ; j++)
					{
						if(
						    tocheck[k]->getNeighbour(j)->isTetrahedron() 
						    && !tocheck[k]->getNeighbour(j)->visited()
						    && tocheck[k]->getNeighbour(j) != ret[i]
						    && static_cast<DelaunayTetrahedron *>(tocheck[k]->getNeighbour(j))->numberOfCommonVertices(ret[i]) > 0
						  )
						{
							tocheck_temp.push_back(static_cast<DelaunayTetrahedron *>(tocheck[k]->getNeighbour(j)));
							tocheck[k]->getNeighbour(j)->visited() = true ;
							toclean.push_back(tocheck_temp.back()) ;
							ret[i]->addNeighbourhood(tocheck_temp.back()) ;
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
		
		std::cerr << " ...done" << std::endl ;
		neighbourhood = true ;
	}
	
	return ret ;
}

void DelaunayTree3D::print() const
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
			tree[i]->print() ;
			
			for(size_t j = 0 ; j< tree[i]->neighbour.size() ; j++)
			{
				std::cout << "\t --> " ;
				tree[i]->getNeighbour(j)->print() ;
			}
		}
	}
	#endif
}

void DelaunayTetrahedron::clearElementaryMatrix()
{
	for(size_t i = 0 ; i <  cachedElementaryMatrix.size() ; i++)
	{
		cachedElementaryMatrix[i].clear() ;
	}
	cachedElementaryMatrix.clear() ;
}

std::vector<std::vector<Matrix> > & DelaunayTetrahedron::getElementaryMatrix() 
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
		getInverseJacobianMatrix(Point( .25, .25, .25), J ) ;
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
	DelaunayTree3D dt(&to_add[0], &to_add[1], &to_add[2], &to_add[4]) ;
	for(size_t i = 4 ; i < to_add.size() ; i++)
	{
		dt.insert(&to_add[i]) ;
	}
	
// 	TetrahedralElement father(QUADRATIC) ;
// 	dt.addSharedNodes(1) ;
// 	dt.refresh( &father) ;
	
	std::vector<DelaunayTetrahedron *> tetra = dt.getTetrahedrons(false) ;
	
	Jinv.resize(tetra.size()*tetra[0]->getGaussPoints().gaussPoints.size(), Matrix()) ;
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
					getInverseJacobianMatrix( gp_temp.gaussPoints[j].first, Jinv[i*gp_temp.gaussPoints.size()+j]) ;
					
					gp_alternative.push_back(gp_temp.gaussPoints[j]) ;
				}
			}
			
		}
	}
	else
	{
		Matrix J ;
		getInverseJacobianMatrix( Point(1./4., 1./4.,1./4.), J ) ;
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
		Jinv.resize(gp.gaussPoints.size(), Matrix()) ;
		for(size_t i = 0 ; i < gp.gaussPoints.size() ;  i++)
		{
			 getInverseJacobianMatrix( gp.gaussPoints[i].first, Jinv[i]) ;
		}
	}
	else
	{
		Matrix J ;
		getInverseJacobianMatrix(Point( 1./4.,1./4., 1./4.) , J);
		Jinv.resize(gp.gaussPoints.size(), J) ;
	}
	
	Vector f(0., 3) ;
		for(size_t i = 0 ; i < getShapeFunctions().size() ; i++)
		{

				behaviour->getForces(this->getState(), getShapeFunction(i),gp, Jinv, f) ;
				
				forces[i*3]+=f[0];
				forces[i*3+1]+=f[1];
				forces[i*3+2]+=f[2];
		}

		for(size_t i = 0 ; i < getEnrichmentFunctions().size() ; i++)
		{
			behaviour->getForces(this->getState(), getEnrichmentFunction(i),gp, Jinv, f) ;
			
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
		Jinv.resize(gp.gaussPoints.size(), Matrix()) ;
		for(size_t i = 0 ; i < gp.gaussPoints.size() ;  i++)
		{
			getInverseJacobianMatrix( gp.gaussPoints[i].first, Jinv[i] ) ;
		}
	}
	else
	{
		Matrix J ;
		getInverseJacobianMatrix(Point( 1./4.,1./4., 1./4.), J );
		Jinv.resize(gp.gaussPoints.size(), J) ;
	}
	
	Vector f(0., 3) ;
	for(size_t i = 0 ; i < getShapeFunctions().size() ; i++)
	{
		behaviour->getForces(this->getState(), getShapeFunction(i) ,gp, Jinv, f) ;
		
		forces[i*3]+=f[0];
		forces[i*3+1]+=f[1];
		forces[i*3+2]+=f[2];
	}

	for(size_t i = 0 ; i < getEnrichmentFunctions().size() ; i++)
	{
		behaviour->getForces(this->getState(), getEnrichmentFunction(i) ,gp, Jinv, f) ;
		
		forces[(i+getShapeFunctions().size())*3]+=f[0];
		forces[(i+getShapeFunctions().size())*3+1]+=f[1];
		forces[(i+getShapeFunctions().size())*3+2]+=f[2];
	}


	return forces ;
}


std::vector<Point *> DelaunayTetrahedron::getIntegrationHints() const
{
	std::vector<Point *> to_add ;
	to_add.push_back(new Point(0,1,0)) ;
	to_add.push_back(new Point(0,0,0)) ;
	to_add.push_back(new Point(1,0.0)) ;
	to_add.push_back(new Point(0,0,1)) ;
	
/*
	std::vector<std::vector<Segment> > segments ;
	
	for(int i = 0 ; i <  (int)getEnrichmentFunctions().size() ; i++)
	{
		std::vector<Segment> temp ;
		for(int j = 0 ; j < (int)getEnrichmentFunction(i).getIntegrationHint().size()-1 ; j++)
		{
			temp.push_back(Segment(getEnrichmentFunction(i).getIntegrationHint(j), getEnrichmentFunction(i).getIntegrationHint(j+1))) ;
		}
		
		segments.push_back(temp) ;
	}
	for(size_t i = 0 ; i < segments.size() ; i++) // per shape function
	{
		for(size_t j = i+1 ; j < segments.size() ; j++) // per other shape function
		{
			for(size_t k = 0 ; k < segments[i].size() ; k++) //per segment in first SF
			{
				for(size_t l = 0 ; l < segments[j].size() ; l++) //per segment in second SF
				{
					if(segments[i][k].intersects(segments[j][l]))
					{
						bool go = true ;
						Point test = segments[i][k].intersection(segments[j][l]) ;
						for(int k = 0 ; k < to_add.size()  ; k++ )
						{
							if(squareDist3D(&test, to_add[k]) 
							   < 4.*POINT_TOLERANCE*POINT_TOLERANCE)
							{
								go = false ;
								break ;
							}
						}
						if(go)
						{
							to_add.push_back(new Point(test)) ;
							if(to_add.back()->x< 0)
								to_add.back()->x = 0 ;
							if(to_add.back()->y < 0)
								to_add.back()->y = 0 ;
							if(to_add.back()->z < 0)
								to_add.back()->z = 0 ;
						}
					}
				}
			}
		}
	}*/
	
	for(size_t i = 0 ; i <  getEnrichmentFunctions().size() ; i++)
	{
		
		for(size_t j = 0 ; j < getEnrichmentFunction(i).getIntegrationHint().size() ; j++)
		{
			bool go = true ;
			for(int k = 0 ; k < to_add.size()  ; k++ )
			{
				if(dist(getEnrichmentFunction(i).getIntegrationHint(j), *to_add[k]) 
				   < 1e-8)
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
				if(to_add.back()->z < 0)
					to_add.back()->z = 0 ;
			}
		}
	}
	
	return to_add ;
}


GaussPointArray DelaunayTetrahedron::getSubTriangulatedGaussPoints() const
{
	if(!enrichmentUpdated)
		return cachedGaussPoints ;

	GaussPointArray gp = getGaussPoints() ; 
	size_t numberOfRefinements = 1;
	
	double tol = std::numeric_limits<double>::epsilon() ;
	double position_tol = 4.*POINT_TOLERANCE ;
	double infinity = .15 ;
	VirtualMachine vm ;
	if(getEnrichmentFunctions().size() > 0)
	{
		std::vector<std::pair<Point, double> > gp_alternative ;
		VirtualMachine vm ;
		std::vector<Point *> to_add = getIntegrationHints();
		std::vector<Point *> pointsToCleanup = to_add;
		std::vector<DelaunayTetrahedron *> triangleToCleanup;
		std::vector<DelaunayTetrahedron *> tri ;
		std::vector<bool> pass ;
		int passNum = 0;
		double J = jacobianAtPoint(Point(1./4., 1./4., 1/4.)) ;
		double lastError = 5 ;
		size_t maxGradientIndex = 0 ;
		std::vector<double> grads(getEnrichmentFunctions().size(), 0.) ;
		TetrahedralElement f(LINEAR) ;
		
		double ndivs = 20 ;
		for(double k = 0  ; k < ndivs ; k++)
		{
			for(double l = 0  ; l < ndivs ; l++)
			{
				for(double m = 0  ; m < ndivs ; m++)
				{
					if( k+l+m < ndivs )
						gp_alternative.push_back(std::make_pair(Point(k/ndivs, l/ndivs, m/ndivs), 1.)) ;
				}
			}
		}
		for(size_t i =0 ; i < gp_alternative.size() ; i++)
		{
			gp_alternative[i].second *= gp.gaussPoints[0].second/gp_alternative.size() ;
		}
		if(gp.gaussPoints.size() < gp_alternative.size())
		{
			
			gp.gaussPoints.resize(gp_alternative.size()) ;
			std::copy(gp_alternative.begin(), gp_alternative.end(), &gp.gaussPoints[0]);
			gp.id = -1 ;
		}
		std::cout << "." << std::flush ;
		return gp ;

		DelaunayTree3D * dt = new DelaunayTree3D(to_add[0], to_add[1], to_add[2], to_add[3]) ;
		if(to_add.size() > 5)
			std::random_shuffle(to_add.begin()+4, to_add.end()) ;

// 		for(size_t i = 4 ; i < 5to_add.size() ; i++)
// 		{
// 			dt->insert(to_add[i]) ;
// 		}

		tri = dt->getTetrahedrons(false) ;
		// original set of tetrahedrons
		
// 		
		pass.resize(tri.size(), false) ;
		
		for(size_t j = 0 ; j < tri.size() ; j++)
		{
			for(size_t k = 0 ; k <  getEnrichmentFunctions().size() ; k++)
			{
				
				double dx = vm.deval(getEnrichmentFunction(k),XI, tri[j]->getCenter()) ;
				double dy = vm.deval(getEnrichmentFunction(k),ETA, tri[j]->getCenter()) ;
				double dz = vm.deval(getEnrichmentFunction(k),ZETA, tri[j]->getCenter()) ;
				grads[k] = std::max((dx*dx+dy*dy+dz*dz), grads[k]) ;
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
			double dz = vm.deval(getEnrichmentFunction(maxGradientIndex),ZETA, tri[j]->getCenter()) ;
			double t_integral = (dx*dx+dy*dy+dz*dz)*tri[j]->volume()*J/.15 ;

			triIntegral.push_back(t_integral) ;
			integral += t_integral ;
			
		}

// 		vm.print(getEnrichmentFunction(maxGradientIndex)) ;// 		if(to_add.size() < 5)
		numberOfRefinements = 3 ;
		for(size_t i = 0 ; i < numberOfRefinements ; i++)
		{

			double newIntegral = 0;

			for(size_t j = 0 ; j < tri.size() ; j++)
			{
				if(!pass[j] )
				{
					Function x = XTransform(tri[j]->getBoundingPoints(), f.getShapeFunctions()) ;
					Function y = YTransform(tri[j]->getBoundingPoints(), f.getShapeFunctions()) ;
					Function z = ZTransform(tri[j]->getBoundingPoints(), f.getShapeFunctions()) ;
					
					double error = 0 ;

					tri[j]->setOrder(QUADRATIC) ;
					GaussPointArray gpquad = tri[j]->getGaussPoints() ;
					tri[j]->setOrder(LINEAR) ;
// 					GaussPointArray gplin = tri[j]->getGaussPoints() ;
					
					
					for(size_t m = 0 ; m < gpquad.gaussPoints.size() ; m++)
					{
						gpquad.gaussPoints[m].first.set(vm.eval(x, gpquad.gaussPoints[m].first), vm.eval(y, gpquad.gaussPoints[m].first), vm.eval(z, gpquad.gaussPoints[m].first)) ;
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
						double dz = vm.deval(getEnrichmentFunction(maxGradientIndex),ZETA, gpquad.gaussPoints[m].first) ;
						newTriIntegral += (dx*dx+dy*dy+dz*dz)*gpquad.gaussPoints[m].second*J ;
					}
					
					newIntegral += newTriIntegral ;
					error = std::abs((triIntegral[j]-newTriIntegral)) ;
					if(error < tol )
					{
						pass[j] = true ;
					}
				}
			}

			
			std::vector<DelaunayTetrahedron *> newTris ;
			std::vector<bool> newPass ;

			for(size_t j = 0 ; j < tri.size() ; j++)
			{
			std::cout << "." << std::flush ;
				if(!pass[j] && tri[j]->area() > sqrt(std::numeric_limits<double>::epsilon()))
				{

					std::pair<std::vector<DelaunayTetrahedron *>, std::vector<Point *> > q =
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
		}
// 		
		for(size_t i = 0 ; i < tri.size() ; i++)
		{

			Function x = XTransform(tri[i]->getBoundingPoints(), f.getShapeFunctions()) ;
			Function y = YTransform(tri[i]->getBoundingPoints(), f.getShapeFunctions()) ;
			Function z = ZTransform(tri[i]->getBoundingPoints(), f.getShapeFunctions()) ;
			tri[i]->setOrder(QUADRATIC) ;

			GaussPointArray gp_temp = tri[i]->getGaussPoints() ;
			
			for(size_t j = 0 ; j < gp_temp.gaussPoints.size() ; j++)
			{

				gp_temp.gaussPoints[j].first.set(vm.eval(x, gp_temp.gaussPoints[j].first), vm.eval(y, gp_temp.gaussPoints[j].first), vm.eval(z, gp_temp.gaussPoints[j].first)) ;
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
// 	std::cout << gp.gaussPoints.size() << std::endl ;
	return gp ;
}

void DelaunayTree3D::setElementOrder(Order elemOrder)
{
	switch(elemOrder)
	{
	case CONSTANT:
		{
			break ;
		}
	case LINEAR:
		{
			break ;
		}
	case QUADRATIC:
		{
			addSharedNodes(1,1,0) ;
			break ;
		}
	case CUBIC:
		{
			addSharedNodes(2,1,0) ;
			break ;
		}
	case QUADRIC:
		{
			addSharedNodes(3,1,0) ;
			break ;
		}
	case QUINTIC:
		{
			addSharedNodes(3,1,0) ;
			break ;
		}
	case CONSTANT_TIME_LINEAR:
		{
			addSharedNodes(0,2,2) ;
			break ;
		}
	case CONSTANT_TIME_QUADRATIC:
		{
			addSharedNodes(0,3,2) ;
			break ;
		}
	case LINEAR_TIME_LINEAR:
		{
			addSharedNodes(0,2,2) ;
			break ;
		}
	case LINEAR_TIME_QUADRATIC:
		{
			addSharedNodes(0,3,2) ;
			break ;
		}
	case QUADRATIC_TIME_LINEAR:
		{
			addSharedNodes(1,2,2) ;
			break ;
		}
	case QUADRATIC_TIME_QUADRATIC:
		{
			addSharedNodes(1,3,2) ;
			break ;
		}
	case CUBIC_TIME_LINEAR:
		{
			addSharedNodes(2,2,2) ;
			break ;
		}
	case CUBIC_TIME_QUADRATIC:
		{
			addSharedNodes(2,3,2) ;
			break ;
		}
	case QUADRIC_TIME_LINEAR:
		{
			addSharedNodes(3,2,2) ;
			break ;
		}
	case QUADRIC_TIME_QUADRATIC:
		{
			addSharedNodes(3,3,2) ;
			break ;
		}
	case QUINTIC_TIME_LINEAR:
		{
			addSharedNodes(3,2,2) ;
			break ;
		}
	case QUINTIC_TIME_QUADRATIC:
		{
			addSharedNodes(3,3,2) ;
			break ;
		}
	default:
		break ;
		
	}
}



// GaussPointArray DelaunayTetrahedron::getSubTriangulatedGaussPoints() const
// {
// 	GaussPointArray gp = getGaussPoints() ; 
// 	
// 	VirtualMachine vm ;
// 	
// 	if(getEnrichmentFunctions().size() > 0 )
// 	{
// 		std::vector<std::pair<Point, double> > gp_alternative ;
// 		std::vector<Point> to_add ;
// 		std::vector<Point> to_add_extra ;
// 		to_add.push_back(Point(0,1,0)) ;
// 		to_add.push_back(Point(0,0,0)) ;
// 		to_add.push_back(Point(1,0,0)) ;
// 		to_add.push_back(Point(0,0,1)) ;
// 		for(size_t i = 0 ; i <  getEnrichmentFunctions().size() ; i++)
// 		{
// 			for(size_t j = 0 ; j < getEnrichmentFunction(i).getIntegrationHint().size() ; j++)
// 			{
// 				if(getEnrichmentFunction(i).getIntegrationHint(j) != to_add[0] && 
// 				   getEnrichmentFunction(i).getIntegrationHint(j) != to_add[1] && 
// 				   getEnrichmentFunction(i).getIntegrationHint(j) != to_add[2] && 
// 				   getEnrichmentFunction(i).getIntegrationHint(j) != to_add[3])
// 				{
// 					bool ok = true ;
// 					for(size_t k = 0 ; k < to_add_extra.size() ; k++)
// 					{
// 						if(getEnrichmentFunction(i).getIntegrationHint(j) == to_add_extra[k])
// 						{
// 							ok = false ;
// 							break ;
// 						}
// 					}
// 					if(ok)
// 						to_add_extra.push_back(getEnrichmentFunction(i).getIntegrationHint(j)) ;
// 				}
// 			}
// 		}
// 
// 		to_add.insert(to_add.end(), to_add_extra.begin(),to_add_extra.end() ) ;
// 		DelaunayTree3D dt(&to_add[0], &to_add[1], &to_add[2],&to_add[3] ) ;
// 		for(size_t i = 4 ; i < to_add.size() ; i++)
// 		{
// 			dt.insert(&to_add[i]) ;
// 		}
// 		
// // 		TetrahedralElement father(QUADRATIC) ;
// // 		dt.addSharedNodes(1) ;
// // 		dt.refresh( &father) ;
// 		
// 		std::vector<DelaunayTetrahedron *> tri = dt.getTetrahedrons(false) ;
// 				
// 		if(moved)
// 		{
// 			for(size_t i = 0 ; i < tri.size() ; i++)
// 			{
// 
// 				Function x = tri[i]->getXTransform() ;
// 				Function y = tri[i]->getYTransform() ;
// 				Function z = tri[i]->getZTransform() ;
// 				tri[i]->order = QUADRATIC ;
// 				GaussPointArray gp_temp = tri[i]->getGaussPoints() ;
// 
// 				for(size_t j = 0 ; j < gp_temp.gaussPoints.size() ; j++)
// 				{
// 	
// 					gp_temp.gaussPoints[j].first.set(vm.eval(x, gp_temp.gaussPoints[j].first), vm.eval(y, gp_temp.gaussPoints[j].first),  vm.eval(z, gp_temp.gaussPoints[j].first)) ;
// 					gp_temp.gaussPoints[j].second *= this->jacobianAtPoint(gp_temp.gaussPoints[j].first) ;
// 					
// 					
// 					gp_alternative.push_back(gp_temp.gaussPoints[j]) ;
// 				}
// 
// 			}
// 		}
// 		else
// 		{
// 			double ja = this->jacobianAtPoint(Point(1./4.,1./4., 1./4.)) ;
// 			for(size_t i = 0 ; i < tri.size() ; i++)
// 			{
// 
// 				Function x = tri[i]->getXTransform() ;
// 				Function y = tri[i]->getYTransform() ;
// 				Function z = tri[i]->getZTransform() ;
// 				
// 				GaussPointArray gp_temp = tri[i]->getGaussPoints() ;
// 				
// 				for(size_t j = 0 ; j < gp_temp.gaussPoints.size() ; j++)
// 				{
// 					gp_temp.gaussPoints[j].second *= ja ;
// 					gp_temp.gaussPoints[j].first.set(vm.eval(x, gp_temp.gaussPoints[j].first), vm.eval(y, gp_temp.gaussPoints[j].first), vm.eval(z, gp_temp.gaussPoints[j].first)) ;
// 					gp_alternative.push_back(gp_temp.gaussPoints[j]) ;
// 				}
// 
// 			}
// 		}
// 		
// 		gp.gaussPoints.resize(gp_alternative.size()) ;
// 		std::copy(gp_alternative.begin(), gp_alternative.end(), &gp.gaussPoints[0]);
// 	}
// 	
// 	return gp ;
// }

void Mu::DelaunayTreeItem3D::print() const
{
}


std::pair<std::vector<DelaunayTetrahedron *>, std::vector<Point *> > Mu::quad(const DelaunayTetrahedron * t)
{
	std::vector<DelaunayTetrahedron* > tris ;
	std::vector<Point *> points ;

	points.push_back(new Point(*t->first + (*t->second - *t->first)*.5) ) ;
	points.push_back(new Point(*t->first + (*t->third - *t->first)*.5 )) ;
	points.push_back(new Point(*t->first + (*t->fourth - *t->first)*.5 )) ;
	points.push_back(new Point(*t->second + (*t->third - *t->second)*.5 )) ;
	points.push_back(new Point(*t->second + (*t->fourth - *t->second)*.5 )) ;
	points.push_back(new Point(*t->third + (*t->fourth - *t->third)*.5 )) ;
	
	tris.push_back(new DelaunayTetrahedron(NULL,NULL, points[0], points[1], points[2], t->first, NULL)) ;
	tris.push_back(new DelaunayTetrahedron(NULL,NULL,points[0], points[3], points[4], t->second, NULL)) ;
	tris.push_back(new DelaunayTetrahedron(NULL,NULL,points[1], points[3], points[5], t->third, NULL)) ;
	tris.push_back(new DelaunayTetrahedron(NULL,NULL,points[2], points[4], points[5], t->fourth, NULL)) ;
	tris.push_back(new DelaunayTetrahedron(NULL,NULL,points[0], points[1], points[2], points[3], NULL)) ;
	tris.push_back(new DelaunayTetrahedron(NULL,NULL,points[0], points[3], points[4], points[1], NULL)) ;
	tris.push_back(new DelaunayTetrahedron(NULL,NULL,points[1], points[3], points[5], points[0], NULL)) ;
	tris.push_back(new DelaunayTetrahedron(NULL,NULL,points[2], points[4], points[5], points[1], NULL)) ;
	return std::make_pair(tris, points) ;
}
