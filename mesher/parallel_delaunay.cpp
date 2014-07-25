//
// C++ Implementation: delaunay
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "parallel_delaunay.h"
#include <omp.h>
#include <limits>
#include <algorithm>

// #define DEBUG
// #undef DEBUG

using namespace Mu ;


ParallelDelaunayTree::ParallelDelaunayTree(Point * p0,  Point *p1,  Point *p2, const std::vector<Geometry *> & domains) : domains(domains)
{
	global_counter = 0 ;
	for(size_t i = 0 ; i < domains.size() ; i++)
	{
		meshes.push_back(new DelaunayTree(p0, p1, p2));
		elementMap.push_back(std::vector<int>());
	}
}

std::vector<DelaunayTriangle *> ParallelDelaunayTree::getConflictingElements(const Point  * p) 
{
	std::valarray<std::vector<DelaunayTriangle *> > conflicts(meshes.size()) ;
	std::valarray<std::vector<DelaunayTriangle *> > boundaries(meshes.size()) ;
#pragma omp parallel for
	for(size_t i = 0 ; i < meshes.size() ; i++)
	{
		std::vector<DelaunayTriangle *> tmpConflicts = meshes[i]->getConflictingElements(p) ;
		for(size_t j = 0 ; j < tmpConflicts.size() ;j++)
		{
			if(elementMap[i][tmpConflicts[j]->index] >= 0)
				conflicts[i].push_back(tmpConflicts[j]) ;
		}
	}
	std::vector<DelaunayTriangle *> ret = conflicts[0] ;
	for(size_t i = 1 ; i < conflicts.size() ; i++)
		ret.insert(ret.end(), conflicts[i].begin(), conflicts[i].end());
	
	return ret ;
}

std::vector<DelaunayTriangle *> ParallelDelaunayTree::getConflictingElements(const Geometry  * p) 
{
	std::valarray<std::vector<DelaunayTriangle *> > conflicts(meshes.size()) ;
	std::valarray<std::vector<DelaunayTriangle *> > boundaries(meshes.size()) ;
	#pragma omp parallel for
	for(size_t i = 0 ; i < meshes.size() ; i++)
	{
		std::vector<DelaunayTriangle *> tmpConflicts = meshes[i]->getConflictingElements(p) ;
		for(size_t j = 0 ; j < tmpConflicts.size() ;j++)
		{
			if(elementMap[i][tmpConflicts[j]->index] >= 0)
				conflicts[i].push_back(tmpConflicts[j]) ;
		}
	}

	std::vector<DelaunayTriangle *> ret = conflicts[0] ;
	for(size_t i = 1 ; i < conflicts.size() ; i++)
		ret.insert(ret.end(), conflicts[i].begin(), conflicts[i].end());
	
	return ret ;
}

void ParallelDelaunayTree::insert(Point * p)
{
	bool insertion = false ;
	std::valarray<std::vector<DelaunayTreeItem *>> newElems(meshes.size()) ;
#pragma omp parallel for
	for(size_t i = 0 ; i < meshes.size() ; i++)
	{
		std::vector<DelaunayTreeItem *> cons = meshes[i]->conflicts(p) ;
		if(cons.empty())
		{
			std::cout << "Failed insertion : in nothing !" << std::endl ;
			exit(0) ;
			continue ;
		}
		
		meshes[i]->neighbourhood = false ;
		
		if(!domains[i]->in(*p))
		{
			bool allout = true ;
			for(size_t j = 0 ; j < cons.size() ; j++)
			{
				if(domains[i]->in(*cons[j]->first) || domains[i]->in(*cons[j]->second) || domains[i]->in(*cons[j]->third))
					allout = false ;
				
				for(size_t k = 0 ; k < cons[j]->neighbour.size() ; k++)
				{
					if(domains[i]->in(*cons[j]->getNeighbour(k)->first) || domains[i]->in(*cons[j]->getNeighbour(k)->second) || domains[i]->in(*cons[j]->getNeighbour(k)->third) )
						allout = false ;
				}
				
				if(cons[j]->father)
				{
					for(size_t k = 0 ; k < cons[j]->getFather()->son.size() ; k++)
					{
						if(domains[i]->in(*cons[j]->getFather()->getSon(k)->first) || 
							domains[i]->in(*cons[j]->getFather()->getSon(k)->second) || 
							domains[i]->in(*cons[j]->getFather()->getSon(k)->third) )
							allout = false ;
					}
					for(size_t k = 0 ; k < cons[j]->getFather()->stepson.size() ; k++)
					{
						if(domains[i]->in(*cons[j]->getFather()->getStepson(k)->first) || 
							domains[i]->in(*cons[j]->getFather()->getStepson(k)->second) || 
							domains[i]->in(*cons[j]->getFather()->getStepson(k)->third) )
							allout = false ;
					}
				}
				
				if(cons[j]->stepfather)
				{
					for(size_t k = 0 ; k < cons[j]->getStepfather()->son.size() ; k++)
					{
						if(cons[j]->getStepfather()->getSon(k))
						{
							if(domains[i]->in(*cons[j]->getStepfather()->getSon(k)->first) || 
								domains[i]->in(*cons[j]->getStepfather()->getSon(k)->second) || 
								domains[i]->in(*cons[j]->getStepfather()->getSon(k)->third) )
								allout = false ;
						}
					}
					for(size_t k = 0 ; k < cons[j]->getFather()->stepson.size() ; k++)
					{
						if(cons[j]->getStepfather()->getStepson(k))
						{
							if(domains[i]->in(*cons[j]->getStepfather()->getStepson(k)->first) || 
								domains[i]->in(*cons[j]->getStepfather()->getStepson(k)->second) || 
								domains[i]->in(*cons[j]->getStepfather()->getStepson(k)->third) )
								allout = false ;
						}
					}
				}
				
				if(allout)
				{
					continue ;
				}
			}
		}
		
		Star * s = new Star(&cons, p) ;
		
		std::vector<DelaunayTreeItem *> ret ;
		
		for(size_t j = 0 ; j < cons.size() ; j++)
		{
			
			if(!cons[j]->onCircumCircle(*p))
			{
				std::vector<DelaunayTreeItem *> temp ;
				cons[j]->insert(temp,p, s) ;
				ret.insert(ret.end(), temp.begin(), temp.end()) ;
			}
		}
		
		s->updateNeighbourhood() ;
		
		bool weGotPlanes = false ;
		
		for(size_t j = 0 ; j< ret.size() ; j++)
		{
			if(ret[j]->isAlive() && ret[j]->isPlane )
			{
				weGotPlanes = true ;
				meshes[i]->plane.push_back((DelaunayDemiPlane*)(ret[j])) ;
			}
		}
		
		if(weGotPlanes)
		{
			for(int k = 0 ; k< (int)meshes[i]->plane.size()-1 ; k++)
			{
				for(int j = k ; j< (int)meshes[i]->plane.size() ; j++)
				{
					meshes[i]->plane[j]->merge(meshes[i]->plane[k]) ;
				}
			}
		}

		for(size_t j = 0 ; j < ret.size() ; j++)
		{
			ret[j]->clearVisited() ;

		}
		for(size_t j = 0 ; j < cons.size() ; j++)
		{
			cons[j]->clearVisited() ;
		}
		
		for(size_t j = 0 ; j < cons.size() ; j++)
		{
			if(!cons[j]->onCircumCircle(*p))
				cons[j]->kill(p) ;
		}
		std::vector<DelaunayDemiPlane *> * hull = meshes[i]->getConvexHull() ;
		meshes[i]->plane.clear() ;
		meshes[i]->plane.insert(meshes[i]->plane.end(), hull->begin(), hull->end()) ;

		for(size_t j = 0 ; j < cons.size() ; j++)
		{
			if(!cons[j]->isAlive() && cons[j]->isTriangle && !cons[j]->isDeadTriangle)
			{
				DelaunayDeadTriangle* dt = new DelaunayDeadTriangle(static_cast<DelaunayTriangle *>(cons[j])) ;
				dt->clearVisited() ;
				std::valarray<Point *> nularray(0) ;
				static_cast<DelaunayTriangle *>(cons[j])->setBoundingPoints(nularray) ;
				meshes[i]->tree[cons[j]->index] = dt ;
				delete cons[j] ;
			}
		}
		delete hull ;
		delete s ;
		newElems[i] = ret ;
		#pragma omp critical
		insertion = true ;
	}
	
	for(size_t i = 0 ; i < newElems.size() ; i++)
	{
		for(size_t j = 0 ; j < newElems[i].size() ; j++)
			elementMap[i].push_back(0) ;
		for(size_t j = 0 ; j < newElems[i].size() ; j++)
		{
			elementMap[i][newElems[i][j]->index] = newElems[i][j]->index ;
			if(!domains[i]->in(*newElems[i][j]->first) || !domains[i]->in(*newElems[i][j]->second) || !domains[i]->in(*newElems[i][j]->third))
			{
				elementMap[i][newElems[i][j]->index] *= -1 ;
			}
			else
			{
				for(size_t k = i+1 ; k < newElems.size() ; k++)
				{
					for(size_t l = 0 ; l < newElems[k].size() ; l++)
					{
						if(*newElems[i][j]->first == *newElems[k][l]->first && 
							*newElems[i][j]->second == *newElems[k][l]->second && 
							*newElems[i][j]->third == *newElems[k][l]->third&& 
							elementMap[i][newElems[i][j]->index] >= 0)
							elementMap[i][newElems[i][j]->index] *= -1 ;
					}
				}
			}
		}
	}
	
	p->id = global_counter++ ;
	
}

std::vector<DelaunayTriangle *> ParallelDelaunayTree::getElements() 
{
	std::vector<std::vector<DelaunayTriangle *>> tris(meshes.size()) ;
	
	#pragma omp parallel for
	for(size_t i = 0 ; i < meshes.size() ;  i++)
	{
		std::vector<DelaunayTriangle *> tmp = meshes[i]->getElements() ;
		for(size_t j = 0 ; j < tmp.size() ;  j++)
		{
			if(elementMap[i][tmp[j]->index] >= 0)
				tris[i].push_back(tmp[j]);
		}
	}
	
	std::vector<DelaunayTriangle *> ret ;
	for(size_t i = 0 ; i < tris.size() ;  i++)
		ret.insert(ret.end(), tris[i].begin(), tris[i].end());
	
	return ret ;
}

void ParallelDelaunayTree::setElementOrder(Order o, double dt )
{
	#pragma omp parallel for
	for(size_t i = 0 ; i < meshes.size() ;  i++)
	{
		meshes[i]->global_counter = global_counter ;
		meshes[i]->setElementOrder(o, dt) ;
	}
	int initial_global_counter = global_counter ;
	for(size_t i = 0 ; i < meshes.size() ;  i++)
	{
		std::vector<DelaunayTriangle *> tmp = meshes[i]->getElements() ;
		for(size_t j = 0 ; j < tmp.size() ;  j++)
		{
			if(elementMap[i][tmp[j]->index] >= 0)
			{
				for(size_t k = 0 ; k < tmp[j]->getBoundingPoints().size() ; k++)
				{
					if(tmp[j]->getBoundingPoint(k).id >=initial_global_counter)
						tmp[j]->getBoundingPoint(k).id = global_counter++ ;
				}
			}
		}
	}
}

void ParallelDelaunayTree::extrude(double dt)
{
	#pragma omp parallel for
	for(size_t i = 0 ; i < meshes.size() ;  i++)
	{
		meshes[i]->global_counter = global_counter ;
		meshes[i]->extrude(dt) ;
	}
	
	int initial_global_counter = global_counter ;
	for(size_t i = 0 ; i < meshes.size() ;  i++)
	{
		std::vector<DelaunayTriangle *> tmp = meshes[i]->getElements() ;
		for(size_t j = 0 ; j < tmp.size() ;  j++)
		{
			if(elementMap[i][tmp[j]->index] >= 0)
			{
				for(size_t k = 0 ; k < tmp[j]->getBoundingPoints().size() ; k++)
				{
					if(tmp[j]->getBoundingPoint(k).id >=initial_global_counter)
						tmp[j]->getBoundingPoint(k).id = global_counter++ ;
				}
			}
		}
	}
}

void ParallelDelaunayTree::extrude(const Vector & dt)
{
	#pragma omp parallel for
	for(size_t i = 0 ; i < meshes.size() ;  i++)
	{
		meshes[i]->global_counter = global_counter ;
		meshes[i]->extrude(dt) ;
	}
	
	int initial_global_counter = global_counter ;
	for(size_t i = 0 ; i < meshes.size() ;  i++)
	{
		std::vector<DelaunayTriangle *> tmp = meshes[i]->getElements() ;
		for(size_t j = 0 ; j < tmp.size() ;  j++)
		{
			if(elementMap[i][tmp[j]->index] >= 0)
			{
				for(size_t k = 0 ; k < tmp[j]->getBoundingPoints().size() ; k++)
				{
					if(tmp[j]->getBoundingPoint(k).id >=initial_global_counter)
						tmp[j]->getBoundingPoint(k).id = global_counter++ ;
				}
			}
		}
	}
}

std::vector<Point * > & ParallelDelaunayTree::getAdditionalPoints() 
{ 
	return additionalPoints ; 
}

const std::vector<Point * > & ParallelDelaunayTree::getAdditionalPoints() const 
{
	return additionalPoints ; 
}

std::vector<DelaunayTreeItem *> & ParallelDelaunayTree::getTree()
{
	return tree ;
}

const std::vector<DelaunayTreeItem *> & ParallelDelaunayTree::getTree() const
{
	return tree ;
}

