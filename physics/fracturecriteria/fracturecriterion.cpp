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

FractureCriterion::FractureCriterion() : neighbourhoodradius(.0005), physicalCharacteristicRadius(.008), scoreAtState(0)
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
		Circle epsilon(neighbourhoodradius,testedTri->getCenter()) ;
		std::vector<DelaunayTriangle *> tempcache = testedTri->tree->getConflictingElements(&epsilon);
		for(size_t i = 0 ; i < tempcache.size() ; i++)
		{
			if(tempcache[i]->getBehaviour()->getFractureCriterion())
			{
				cache.push_back(tempcache[i]);
				area.push_back(cache.back()->area());
			}
		}
	}
	else if(testedTet)
	{
		if(!cache3d.empty())
			return ;
		Sphere epsilon(neighbourhoodradius,testedTet->getCenter()) ;
		std::vector<DelaunayTetrahedron *> tempcache3d = testedTet->tree->getConflictingElements(&epsilon);
		for(size_t i = 0 ; i < tempcache3d.size() ; i++)
		{
			if(tempcache3d[i]->getBehaviour()->getFractureCriterion())
			{
				cache3d.push_back(tempcache3d[i]);
				area.push_back(cache3d.back()->volume());
			}
		}
	}
}

void FractureCriterion::step(const ElementState &s)
{
	scoreAtState = grade(s) ;
}

FractureCriterion::~FractureCriterion()
{
}

void FractureCriterion::setNeighbourhoodRadius(double r)
{
	neighbourhoodradius = r ;
	cache.clear() ;
}

void FractureCriterion::setMaterialCharacteristicRadius(double r)
{
	physicalCharacteristicRadius = r ;
}

bool FractureCriterion::met(const ElementState &s)
{
	if( s.getParent()->getBehaviour()->fractured())
		return false ;
	double tol = 0 ;
	DelaunayTriangle * testedTri = dynamic_cast<DelaunayTriangle *>(s.getParent()) ;
	DelaunayTetrahedron * testedTet = dynamic_cast<DelaunayTetrahedron *>(s.getParent()) ;
	HexahedralElement * testedHex = dynamic_cast<HexahedralElement *>(s.getParent()) ;
	if(testedTri)
	{
		if(cache.empty())
			initialiseCache(s) ;

		if(testedTri->visited)
			return false ;
				
		if (scoreAtState <= 0)
			return false ;

		double maxNeighbourhoodScore = 0 ;
		double matchedArea = 0 ;
		std::map<double, DelaunayTriangle *> scores ;
		std::vector<double> unsortedScores ;
		std::map<DelaunayTriangle *, double> areatemp ;
		DelaunayTriangle * maxLocus = NULL;
		
		if(!cache.empty())
		{
			for(size_t i = 0 ; i< cache.size() ; i++)
			{

				if( cache[i]->getBehaviour()->getFractureCriterion() 
					&& !cache[i]->getBehaviour()->fractured())
				{
					double s = cache[i]->getBehaviour()->getFractureCriterion()->getSteppedScore() ;
					scores[-s] =  cache[i];
					unsortedScores.push_back(s);
					if(s > maxNeighbourhoodScore)
					{
						maxNeighbourhoodScore = s ;
						maxLocus = cache[i] ;
					}
				}
				else if(cache[i]->getBehaviour()->getFractureCriterion() 
					&& cache[i]->getBehaviour()->fractured())
				{
					double s = 1 ;
					scores[-s] =  cache[i];
					unsortedScores.push_back(s);
					if(s > maxNeighbourhoodScore)
					{
						maxNeighbourhoodScore = s ;
						maxLocus = cache[i] ;
					}
					
				}
				areatemp[cache[i]] = area[i] ;
				
// 				if ((maxNeighbourhoodScore*tol) > score)
// 					return false ;
					
			}
		}
		
		if(!maxLocus)
			return false ;
		
		bool foundcutoff = false ;
		double thresholdscore = maxNeighbourhoodScore ;
		
		for(std::map<double, DelaunayTriangle *>::iterator i = scores.begin() ; i != scores.end() ; ++i)
		{
			
			if(!foundcutoff)
			{
				if(-i->first > 0 )
				{
					matchedArea += areatemp[i->second] ;
				}
				if (matchedArea > physicalCharacteristicRadius*physicalCharacteristicRadius*M_PI)
				{
					thresholdscore = -i->first ;
					foundcutoff  = true ;
					break ;
				}
			}
// 			else
// 				i->second->visited = true ;
		}
		if (!foundcutoff )
			return false ;

		if (squareDist2D(maxLocus->getCenter(), s.getParent()->getCenter()) > physicalCharacteristicRadius*physicalCharacteristicRadius)
			return false ;
		
		return true ;

	}
	if(testedTet)
	{
		if(testedTet->visited())
			return false ;
		
		double score = grade(s) ;
		
		if (score <= 0)
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

Material FractureCriterion::toMaterial()
{
	Material mat ;
	return mat ;
}

}

