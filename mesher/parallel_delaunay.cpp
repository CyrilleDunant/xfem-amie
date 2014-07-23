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
		meshes.push_back(new DelaunayTree(p0, p1, p2, p3));
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

	return conflicts ;
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

	return conflicts ;
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
			return ;
		}
		
		meshes[i]->neighbourhood = false ;
		bool allout= true ;
		
		for(size_t j = 0 ; j < cons.size() ; j++)
		{
			if(cons[j]->isVertex(p)) 
			{
				return ;
			}
			
			if(!allout)
			{
				allout = allout && !domains[i]->in(cons[j]->getCenter()) && !domains[i]->in(cons[j]->first) && !domains[i]->in(cons[j]->second) && !domains[i]->in(cons[j]->third) ;
				for(size_t k = 0 ; k < cons[j]->neighbour.size() ; k++)
				{
					if(domains[i]->in(cons[j]->getNeighbour(k)->first) || domains[i]->in(cons[j]->getNeighbour(k)->second) || domains[i]->in(cons[j]->getNeighbour(k)->third) || domains[i]->in(cons[j]->getNeighbour(k)->getCenter()))
						allout = false ;
				}
			}
		}
		
		if(allout)
			return ;
		
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
			for(int k = 0 ; k< (int)plane.size()-1 ; k++)
			{
				for(int j = k ; j< (int)plane.size() ; j++)
				{
					meshes[i]->plane[j]->merge(plane[k]) ;
				}
			}
		}

		for(size_t j = 0 ; j < ret.size() ; j++)
		{
			ret[i]->clearVisited() ;

		}
		for(size_t j = 0 ; j < cons.size() ; j++)
		{
			cons[i]->clearVisited() ;
		}
		
		for(size_t j = 0 ; i < cons.size() ; j++)
		{
			if(!cons[i]->onCircumCircle(*p))
				cons[i]->kill(p) ;
		}
		std::vector<DelaunayDemiPlane *> * hull = meshes[i]->getConvexHull() ;
		meshes[i]->plane.clear() ;
		meshes[i]->plane.insert(plane.end(), hull->begin(), hull->end()) ;

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
			for(size_t k = i+1 ; k < newElems.size() ; k++)
			{
				for(size_t l = 0 ; l < newElems[k].size() ; l++)
				{
					if(*newElems[i][j]->first == *newElems[k][l]->first && *newElems[i][j]->second == *newElems[k][l]->second && *newElems[i][j]->third == *newElems[k][l]->third)
						elementMap[i][newElems[i][j]->index] *= -1 ;
				}
			}
		}
	}
	
	p->id = global_counter++ ;
	
}

std::vector<DelaunayTriangle *> ParallelDelaunayTree::getElements() 
{
	std::vector<DelaunayTriangle *> tris ;
	for(size_t i = 0 ; i < meshes.size() ;  i++)
	{
		std::vector<DelaunayTriangle *> tmp = meshes[i]->getElements() ;
		for(size_t j = 0 ; j < tmp.size() ;  j++)
		{
			if(elementMap[i][tmp[j]->index] >= 0)
				tris.push_back(tmp[j]);
		}
		
	}
	return tris ;
}

void ParallelDelaunayTree::setElementOrder(Order o, double dt )
{
	#pragma omp parallel for
	for(size_t i = 0 ; i < meshes.size() ;  i++)
		meshes[i]->setElementOrder(o, dt) ;
	
	additionalPoints = meshes[0]->getAdditionalPoints() ;
	for(size_t i = 0 ; i < additionalPoints.size() ;  i++)
	{
		additionalPoints[i]->id = global_counter++ ;
	}
	meshes[0]->getAdditionalPoints().clear() ;
	for(size_t i = 1 ; i < meshes.size() ;  i++)
	{
		for(size_t j = 0 ; j < additionalPoints.size() ;  j++)
		{
			if(meshes[i]->additionalPoints)
		}
	}
	
// 	std::valarray<std::vector<DelaunayTriangle *> > boundaries(meshes.size()) ;
// 	#pragma omp parallel for
// 	for(size_t i = 0 ; i < meshes.size() ; i++)
// 	{
// 		std::vector<DelaunayTriangle *> tmpConflicts = meshes[i]->getConflictingElements(p) ;
// 		for(size_t j = 0 ; j < tmpConflicts.size() ;j++)
// 		{
// 			if(!(domains[i]->in(tmpConflicts[j]->first) && domains[i]->in(tmpConflicts[j]->second)&& domains[i]->in(tmpConflicts[j]->third)))
// 				boundaries[i].push_back(tmpConflicts[j]) ;
// 		}
// 	}
// 	std::vector<DelaunayTriangle *> trueBounds ;
// 	for(size_t i = 0 ; i < boundaries.size() ; i++)
// 	{
// 		for(size_t j = i+1 ; j < boundaries.size() ; j++)
// 		{
// 			if(boundaries[i]->first->id == boundaries[j]->first->id &&
// 				 boundaries[i]->second->id == boundaries[j]->second->id &&
// 				 boundaries[i]->third->id == boundaries[j]->third->id
// 			)
// 			{
// 				bool already = false ;
// 				for(size_t k = 0 ; k < trueBounds.size() ; k++)
// 				{
// 					if(boundaries[i]->first->id == trueBounds[k]->first->id &&
// 				     boundaries[i]->second->id == trueBounds[k]->second->id &&
// 				     boundaries[i]->third->id == trueBounds[k]->third->id)
// 					{
// 						already = true ;
// 						break ;
// 					}
// 				}
// 				if(!already)
// 				{
// 					trueBounds.push_back(boundaries[i]);
// 					break ;
// 				}
// 			}
// 		}
// 	}
// 	
// 	std::vector<Point * > extraPoints ;
// 	for()
}

void ParallelDelaunayTree::extrude(double dt)
{
	#pragma omp parallel for
	for(size_t i = 0 ; i < meshes.size() ;  i++)
		meshes[i]->extrude(dt) ;
}

void ParallelDelaunayTree::extrude(const Vector & dt)
{
	#pragma omp parallel for
	for(size_t i = 0 ; i < meshes.size() ;  i++)
		meshes[i]->extrude(dt) ;
}

