//
// C++ Implementation: fracturecriterion
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "fracturecriterion.h"
#include "../../mesher/delaunay.h"
#include "../../mesher/delaunay_3d.h"
namespace Mu {

FractureCriterion::FractureCriterion() : eps(.0005)
{
}

void FractureCriterion::initialiseCache(const ElementState & s)
{
	DelaunayTriangle * testedTri = dynamic_cast<DelaunayTriangle *>(s.getParent()) ;
	DelaunayTetrahedron * testedTet = dynamic_cast<DelaunayTetrahedron *>(s.getParent()) ;
	if(testedTri)
	{
		if(!cache.empty())
			return ;
		Circle epsilon(eps,testedTri->getCenter()) ;
		cache = testedTri->tree->getConflictingElements(&epsilon);
	}
	else if(testedTet)
	{
		if(!cache3d.empty())
			return ;
		Sphere epsilon(eps,testedTet->getCenter()) ;
		cache3d = testedTet->tree->getConflictingElements(&epsilon);
	}
}

FractureCriterion::~FractureCriterion()
{
}

void FractureCriterion::setNeighbourhoodRadius(double r)
{
	eps = r ;
	cache.clear() ;
}

bool FractureCriterion::met(const ElementState &s)
{
	if( s.getParent()->getBehaviour()->fractured())
		return false ;
	double tol = .9999999 ;
	DelaunayTriangle * testedTri = dynamic_cast<DelaunayTriangle *>(s.getParent()) ;
	DelaunayTetrahedron * testedTet = dynamic_cast<DelaunayTetrahedron *>(s.getParent()) ;
	HexahedralElement * testedHex = dynamic_cast<HexahedralElement *>(s.getParent()) ;
	if(testedTri)
	{
		if(cache.empty())
		{
			Circle epsilon(eps,testedTri->getCenter()) ;
			std::vector<DelaunayTriangle *> tris = testedTri->tree->getElements() ;
			std::vector<bool> vis(tris.size()) ;
			for(size_t i = 0 ; i != tris.size() ; ++i)
			{
				vis[i] = tris[i]->visited ;
				tris[i]->clearVisited() ;
			}
			cache = testedTri->tree->getConflictingElements(&epsilon);
			for(size_t i = 0 ; i != tris.size() ; ++i)
			{
				tris[i]->visited = vis[i];
			}
		}
// 		cache = testedTri->tree->getTriangles() ;
		if(testedTri->visited)
			return false ;
		
		double score = grade(s) ;
		
		if (score == 0)
		{
			testedTri->visited = true ;
			return false ;
		}

		double maxNeighbourhoodScore = 0 ;
		
		std::vector<double> scores(cache.size(), double(0)) ;
		if(!cache.empty())
		{
			for(size_t i = 0 ; i< cache.size() ; i++)
			{

				if( cache[i]->getBehaviour()->getFractureCriterion() 
					&& !cache[i]->getBehaviour()->fractured())
				{
					scores[i] = cache[i]->getBehaviour()->getFractureCriterion()->grade(cache[i]->getState()) ;
					maxNeighbourhoodScore = std::max(maxNeighbourhoodScore,scores[i]) ;
				}

				
				if ((maxNeighbourhoodScore*tol) > score)
					return false ;
					
			}
		}
		
		if(score >= maxNeighbourhoodScore*tol)
		{
			for(size_t i = 0 ; i< cache.size() ; i++)
			{
				if(scores[i] < maxNeighbourhoodScore*tol)
					cache[i]->visited = true ;
			}
			testedTri->visited = true ;
			return true ;
		}
		return false ;
	}
	if(testedTet)
	{
		if(testedTet->visited())
			return false ;
		
		double score = grade(s) ;
		
		if (score == 0)
		{
			testedTet->visited() = true ;
			return false ;
		}
		
		if(cache3d.empty())
		{
			Sphere epsilon(0.0002,testedTet->getCenter()) ;
			cache3d = testedTet->tree->getConflictingElements(&epsilon);
		}

		double maxNeighbourhoodScore = 0 ;

		std::vector<double> scores(cache3d.size(), double(0)) ;
		if(!cache3d.empty())
		{
			for(size_t i = 0 ; i< cache3d.size() ; i++)
			{
				
				if( cache3d[i]->getBehaviour()->getFractureCriterion() 
					&& !cache3d[i]->getBehaviour()->fractured())
				{
					scores[i] = cache3d[i]->getBehaviour()->getFractureCriterion()->grade(cache3d[i]->getState()) ;
					maxNeighbourhoodScore = std::max(maxNeighbourhoodScore,scores[i]) ;
				}

				if ((maxNeighbourhoodScore*tol) > score)
					return false ;
				
			}
		}
		
		
		if(score >= maxNeighbourhoodScore*tol)
		{
			for(size_t i = 0 ; i< cache3d.size() ; i++)
			{
				if(scores[i] < maxNeighbourhoodScore*tol)
					cache3d[i]->visited() = true;
			}
			testedTet->visited() = true ;
			return true ;
		}
		return false ;
	}
	else if(testedHex)
	{
		std::set<HexahedralElement *> neighbourhood ;
		std::vector<HexahedralElement *> neighbours = testedHex->neighbourhood ;
		for(size_t i = 0 ; i < neighbours.size() ; i++)
		{
			for(size_t j = 0 ; j <  neighbours[i]->neighbourhood.size() ; j++)
			{
				if(neighbours[i]->neighbourhood[j] != testedHex 
				   && !neighbours[i]->neighbourhood[j]->getBehaviour()->fractured())
					neighbourhood.insert(neighbours[i]->neighbourhood[j]) ;
			}
		}
		double score = grade(s) ;
		double maxNeighbourhoodScore = 0 ;
		if(!neighbourhood.empty())
		{
			for(std::set<HexahedralElement *>::iterator i= neighbourhood.begin() ; i != neighbourhood.end() ; ++i)
			{
				if((*i)->getBehaviour()->getFractureCriterion() 
					&& !(*i)->getBehaviour()->fractured())
					maxNeighbourhoodScore = std::max(maxNeighbourhoodScore,
					                                 (*i)->getBehaviour()->getFractureCriterion()->grade((*i)->getState())) ;
				if((*i)->getBehaviour()->changed())
				{
					maxNeighbourhoodScore = 10.*score ;
					break ;
				}
				
				if (maxNeighbourhoodScore > score)
					break ;
				
			}
		}
		
		if( score > 0 )
		{
			if(score > maxNeighbourhoodScore)
			{
				return true ;
			}
		}

		return false ;
	}
	else
	{
		std::cout << " criterion not implemented for this kind of element" << std::endl ;
		return false ;
	}
	
	//shut up the compiler
	return false ;
}
}

