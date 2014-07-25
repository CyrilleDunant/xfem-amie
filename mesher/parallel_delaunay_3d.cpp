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

#include "parallel_delaunay_3d.h"
#include <omp.h>
#include <limits>
#include <algorithm>

// #define DEBUG
// #undef DEBUG

using namespace Mu ;


ParallelDelaunayTree3D::ParallelDelaunayTree3D(Point * p0,  Point *p1,  Point *p2,  Point *p3, const std::vector<Geometry *> & domains) : domains(domains)
{
	global_counter = 0 ;
	for(size_t i = 0 ; i < domains.size() ; i++)
	{
		meshes.push_back(new DelaunayTree3D(p0, p1, p2, p3));
		elementMap.push_back(std::vector<int>());
	}
}

std::vector<DelaunayTetrahedron *> ParallelDelaunayTree3D::getConflictingElements(const Point  * p) 
{
	std::valarray<std::vector<DelaunayTetrahedron *> > conflicts(meshes.size()) ;
	std::valarray<std::vector<DelaunayTetrahedron *> > boundaries(meshes.size()) ;
#pragma omp parallel for
	for(size_t i = 0 ; i < meshes.size() ; i++)
	{
		std::vector<DelaunayTetrahedron *> tmpConflicts = meshes[i]->getConflictingElements(p) ;
		for(size_t j = 0 ; j < tmpConflicts.size() ;j++)
		{
			if(elementMap[i][tmpConflicts[j]->index] >= 0)
				conflicts[i].push_back(tmpConflicts[j]) ;
		}
	}
	std::vector<DelaunayTetrahedron *> ret = conflicts[0] ;
	for(size_t i = 1 ; i < conflicts.size() ; i++)
		ret.insert(ret.end(), conflicts[i].begin(), conflicts[i].end());
	
	return ret ;
}

std::vector<DelaunayTetrahedron *> ParallelDelaunayTree3D::getConflictingElements(const Geometry  * p) 
{
	std::valarray<std::vector<DelaunayTetrahedron *> > conflicts(meshes.size()) ;
	std::valarray<std::vector<DelaunayTetrahedron *> > boundaries(meshes.size()) ;
	#pragma omp parallel for
	for(size_t i = 0 ; i < meshes.size() ; i++)
	{
		std::vector<DelaunayTetrahedron *> tmpConflicts = meshes[i]->getConflictingElements(p) ;
		for(size_t j = 0 ; j < tmpConflicts.size() ;j++)
		{
			if(elementMap[i][tmpConflicts[j]->index] >= 0)
				conflicts[i].push_back(tmpConflicts[j]) ;
		}
	}

	std::vector<DelaunayTetrahedron *> ret = conflicts[0] ;
	for(size_t i = 1 ; i < conflicts.size() ; i++)
		ret.insert(ret.end(), conflicts[i].begin(), conflicts[i].end());
	
	return ret ;
}

void ParallelDelaunayTree3D::insert(Point * p)
{
	std::valarray<std::vector<DelaunayTreeItem3D *>> newElems(meshes.size()) ;
	#pragma omp parallel for
	for(size_t i = 0 ; i < meshes.size() ; i++)
	{
		std::vector<DelaunayTreeItem3D *> cons = meshes[i]->conflicts(p) ;
		if(cons.empty())
		{
			std::cout << "Failed insertion : in nothing !" << std::endl ;
			exit(0) ;
			continue ;
		}
		
		meshes[i]->neighbourhood = false ;
		bool isVertex = false ;
		

		
		for(size_t j = 0 ; j < cons.size() ; j++)
		{
			if(cons[j]->isVertex(p)) 
				isVertex = true ;
		}
		
		if(isVertex)
			continue ;
		
		if(!domains[i]->in(*p))
		{
			bool allout = true ;
			for(size_t j = 0 ; j < cons.size() ; j++)
			{
				if(domains[i]->in(*cons[j]->first) || domains[i]->in(*cons[j]->second) || domains[i]->in(*cons[j]->third) || domains[i]->in(*cons[j]->fourth))
					allout = false ;
				
				for(size_t k = 0 ; k < cons[j]->neighbour.size() ; k++)
				{
					if(domains[i]->in(*cons[j]->getNeighbour(k)->first) || domains[i]->in(*cons[j]->getNeighbour(k)->second) || domains[i]->in(*cons[j]->getNeighbour(k)->third) || domains[i]->in(*cons[j]->getNeighbour(k)->fourth))
						allout = false ;
				}
				
				if(cons[j]->father != -1)
				{
					for(size_t k = 0 ; k < cons[j]->getFather()->son.size() ; k++)
					{
						if(domains[i]->in(*cons[j]->getFather()->getSon(k)->first) || 
							domains[i]->in(*cons[j]->getFather()->getSon(k)->second) || 
							domains[i]->in(*cons[j]->getFather()->getSon(k)->third) || 
							domains[i]->in(*cons[j]->getFather()->getSon(k)->fourth))
							allout = false ;
					}
					for(size_t k = 0 ; k < cons[j]->getFather()->stepson.size() ; k++)
					{
						if(domains[i]->in(*cons[j]->getFather()->getStepson(k)->first) || 
							domains[i]->in(*cons[j]->getFather()->getStepson(k)->second) || 
							domains[i]->in(*cons[j]->getFather()->getStepson(k)->third) || 
							domains[i]->in(*cons[j]->getFather()->getStepson(k)->fourth))
							allout = false ;
					}
				}
				
				if(cons[j]->stepfather != -1)
				{
					for(size_t k = 0 ; k < cons[j]->getStepfather()->son.size() ; k++)
					{
						if(cons[j]->getStepfather()->getSon(k))
						{
							if(domains[i]->in(*cons[j]->getStepfather()->getSon(k)->first) || 
								domains[i]->in(*cons[j]->getStepfather()->getSon(k)->second) || 
								domains[i]->in(*cons[j]->getStepfather()->getSon(k)->third) || 
								domains[i]->in(*cons[j]->getStepfather()->getSon(k)->fourth))
								allout = false ;
						}
					}
					for(size_t k = 0 ; k < cons[j]->getFather()->stepson.size() ; k++)
					{
						if(cons[j]->getStepfather()->getStepson(k))
						{
							if(domains[i]->in(*cons[j]->getStepfather()->getStepson(k)->first) || 
								domains[i]->in(*cons[j]->getStepfather()->getStepson(k)->second) || 
								domains[i]->in(*cons[j]->getStepfather()->getStepson(k)->third) || 
								domains[i]->in(*cons[j]->getStepfather()->getStepson(k)->fourth))
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

		Star3D *s = new Star3D( &cons, p ) ;

		std::vector<DelaunayTreeItem3D *> ret ;

		for( size_t i = 0 ; i < cons.size() ; i++ )
		{
			if( !cons[i]->onCircumSphere( *p ) )
			{
				std::vector<DelaunayTreeItem3D *> temp ;
				cons[i]->insert( temp, p, s ) ;
				ret.insert( ret.end(), temp.begin(), temp.end() ) ;
			}
		}

		s->updateNeighbourhood() ;
		bool weGotSpaces = false ;

		for( size_t j = 0 ; j < ret.size() ; j++ )
		{
			ret[j]->clearVisited() ;

		}

		for( size_t j = 0 ; j < ret.size() ; j++ )
		{
			if( ret[j]->isAlive() && ret[j]->isSpace() )
			{
				weGotSpaces = true ;
				meshes[i]->space.push_back( static_cast<DelaunayDemiSpace *>( ret[j] ) ) ;
			}
		}

		if( weGotSpaces )
		{
			for( size_t k = 0 ; k < meshes[i]->space.size() - 1 ; k++ )
			{
				for( size_t j = k + 1 ; j < meshes[i]->space.size() ; j++ )
				{
					meshes[i]->space[j]->merge( meshes[i]->space[k] ) ;
				}
			}
		}

		for( size_t l = 0; l < ret.size(); l++ )
		{
			for( size_t k = 0; k < meshes[i]->space.size(); k++ )
			{
				makeNeighbours( ret[l], meshes[i]->space[k] );
			}
		}

		for( size_t j = 0 ; j < cons.size() ; j++ )
		{
			if( !cons[j]->onCircumSphere( *p ) )
				cons[j]->kill( p ) ;
		}


		for( size_t j = 0 ; j < ret.size() ; j++ )
		{
			if (! (!ret[j]->erased() && ( ( ret[j]->isAlive() && ret[j]->isTetrahedron() ) || ret[j]->isSpace() )) )
			{

				std::valarray<Point *> nullarray( 0 ) ;

				if( ret[j]->isTetrahedron() )
					static_cast<DelaunayTetrahedron *>( ret[j] )->setBoundingPoints( nullarray ) ;

				tree[ret[j]->index] = nullptr ;
				delete ret[j] ;
			}
		}

		for( auto j = meshes[i]->space.begin() ; j != meshes[i]->space.end() ; j++ )
		{
			if( !( *j )->isAlive() )
			{
				meshes[i]->space.erase( j ) ;
				j-- ;
			}
		}

		for( size_t j = 0 ; j < cons.size() ; j++ )
		{

			if( !cons[j]->isAlive() && cons[j]->isTetrahedron() && !cons[j]->isDeadTetrahedron() )
			{
				DelaunayDeadTetrahedron *dt = new DelaunayDeadTetrahedron( static_cast<DelaunayTetrahedron *>( cons[j] ) ) ;
				dt->clearVisited() ;
				meshes[i]->tree[cons[j]->index] = dt ;
				delete cons[j] ;
			}
		}

		delete s ;
		newElems[i] = ret ;
	}
	
	p->id = global_counter++ ;
	for(size_t i = 0 ; i < newElems.size() ; i++)
	{
		int maxIdx = 0 ;
		for(size_t j = 0 ; j < newElems[i].size() ; j++)
			maxIdx = std::max(maxIdx, newElems[i][j]->index) ;
		
		for(int j = elementMap[i].size()-1 ; j < maxIdx ; j++)
			elementMap[i].push_back(0) ;
		
		for(size_t j = 0 ; j < newElems[i].size() ; j++)
		{
			elementMap[i][newElems[i][j]->index] = newElems[i][j]->index ;

			for(size_t k = i+1 ; k < newElems.size() ; k++)
			{
				for(size_t l = 0 ; l < newElems[k].size() ; l++)
				{
					if(newElems[i][j]->isVertex( newElems[k][l]->first) && 
						newElems[i][j]->isVertex( newElems[k][l]->second) && 
						newElems[i][j]->isVertex( newElems[k][l]->third) && 
						newElems[i][j]->isVertex( newElems[k][l]->fourth) && 
						elementMap[i][newElems[i][j]->index] >= 0)
						elementMap[i][newElems[i][j]->index] *= -1 ;
				}
			}
		}
	}
	
	
	
}

std::vector<DelaunayTetrahedron *> ParallelDelaunayTree3D::getElements() 
{
	std::vector<std::vector<DelaunayTetrahedron *>> tris(meshes.size()) ;
	
	#pragma omp parallel for
	for(size_t i = 0 ; i < meshes.size() ;  i++)
	{
		std::vector<DelaunayTetrahedron *> tmp = meshes[i]->getElements() ;
		for(size_t j = 0 ; j < tmp.size() ;  j++)
		{
			if(elementMap[i][tmp[j]->index] >= 0)
				tris[i].push_back(tmp[j]);
		}
	}
	
	std::vector<DelaunayTetrahedron *> ret ;
	for(size_t i = 0 ; i < tris.size() ;  i++)
		ret.insert(ret.end(), tris[i].begin(), tris[i].end());
	
	return ret ;
}

void ParallelDelaunayTree3D::setElementOrder(Order o, double dt )
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
		std::vector<DelaunayTetrahedron *> tmp = meshes[i]->getElements() ;
		for(size_t j = 0 ; j < tmp.size() ;  j++)
		{
			if(elementMap[i][tmp[j]->index] >= 0)
			{
				for(size_t k = 0 ; k < tmp[j]->getBoundingPoints().size() ; k++)
				{
					if(tmp[j]->getBoundingPoint(k).id >= initial_global_counter)
						tmp[j]->getBoundingPoint(k).id = global_counter++ ;
				}
			}
		}
	}
}

void ParallelDelaunayTree3D::extrude(double dt)
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
		std::vector<DelaunayTetrahedron *> tmp = meshes[i]->getElements() ;
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

void ParallelDelaunayTree3D::extrude(const Vector & dt)
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
		std::vector<DelaunayTetrahedron *> tmp = meshes[i]->getElements() ;
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

std::vector<Point * > & ParallelDelaunayTree3D::getAdditionalPoints() 
{ 
	return additionalPoints ; 
}

const std::vector<Point * > & ParallelDelaunayTree3D::getAdditionalPoints() const 
{
	return additionalPoints ; 
}

std::vector<DelaunayTreeItem3D *> & ParallelDelaunayTree3D::getTree()
{
	return tree ;
}

const std::vector<DelaunayTreeItem3D *> & ParallelDelaunayTree3D::getTree() const
{
	return tree ;
}

