//
// C++ Implementation: crackinitiation
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "crackinitiation.h"

namespace Mu {

CrackInitiation::CrackInitiation()
 : EnrichmentBehaviour()
{
}


CrackInitiation::~CrackInitiation()
{
}

std::vector<EnrichmentFeature *> CrackInitiation::step(double delta, DelaunayTree * dt) const
{
	std::vector<EnrichmentFeature *> ret ;
	
	std::vector<DelaunayTriangle *> triangles = dt->getTriangles() ;
	std::vector<DelaunayTriangle *> brokenTriangles ;
	for(size_t i = 0 ; i < triangles.size() ;i++)
	{
		if(triangles[i]->getBehaviour()->fractured())
			brokenTriangles.push_back(triangles[i]) ;
	}
	std::vector<std::vector<DelaunayTriangle *> > clusters ;
	for(size_t i = 0 ; i < brokenTriangles.size() ; i++)
	{
		std::vector<DelaunayTriangle *> cluster;

		while(i < brokenTriangles.size() && brokenTriangles[i]->visited)
			i++ ;

		if(i >= brokenTriangles.size())
			break ;
		
		cluster.push_back(brokenTriangles[i]) ;
		brokenTriangles[i]->visited = true ;
		std::vector<DelaunayTriangle *> toTest ;
		for(size_t j = 0 ; j < brokenTriangles[i]->neighbourhood.size() ; j++)
			toTest.push_back(brokenTriangles[i]->getNeighbourhood(j));
		
		while(!toTest.empty())
		{
			std::vector<DelaunayTriangle *> newToTest ;

			for(size_t j = 0 ; j < toTest.size() ; j++)
			{
				if(toTest[j]->getBehaviour()->fractured() && !toTest[j]->visited)
				{
					toTest[j]->visited = true ;
					cluster.push_back(toTest[j]) ;
					for(size_t k = 0 ; k < toTest[j]->neighbourhood.size() ; k++)
						newToTest.push_back(toTest[j]->getNeighbourhood(k));
				}
			}
			auto e = std::unique(newToTest.begin(),newToTest.end() ) ;
			newToTest.erase(e, newToTest.end()) ;
			toTest = newToTest ;
			
		}

		clusters.push_back(cluster) ;
	}

	for(size_t i = 0 ; i < brokenTriangles.size() ; i++)
		brokenTriangles[i]->clearVisited() ;

	for(size_t i = 0 ; i < clusters.size() ; i++)
	{
		if(clusters[i].size() > 3)
		{
			double cx_bottom = clusters[i][0]->getCenter().getX() ;
			double cy_bottom = clusters[i][0]->getCenter().getY() ;
			double cx_left = clusters[i][0]->getCenter().getX() ;
			double cy_left = clusters[i][0]->getCenter().getY() ;
			double count = 1 ;
			for(size_t j = 1 ; j < clusters[i].size() ; j++)
			{
				if(clusters[i][j]->getCenter().getY() < cx_bottom)
				{
					cx_bottom = clusters[i][j]->getCenter().getX() ;
					cy_bottom = clusters[i][j]->getCenter().getY() ;
				}
				if(clusters[i][j]->getCenter().getX() < cx_left)
				{
					cx_left = clusters[i][j]->getCenter().getX() ;
					cy_left = clusters[i][j]->getCenter().getY() ;
				}

				count++ ;
			}

			cx_left /= count ;
			cy_left /= count ;
			cx_bottom /= count ;
			cy_bottom /= count ;
			
			double d_l = 0 ;
			double d_b = 0 ;
			double dx_l = 0 ;
			double dy_l = 0 ;
			double dx_b = 0 ;
			double dy_b = 0 ;

			double cx = 0 ;
			double cy = 0 ;
			
			for(size_t j = 0 ; j < clusters[i].size() ; j++)
			{
				dx_l += clusters[i][j]->getCenter().getX() - cx_left;
				dy_l += clusters[i][j]->getCenter().getY() - cy_left;
				dx_b += clusters[i][j]->getCenter().getX() - cx_bottom;
				dy_b += clusters[i][j]->getCenter().getY() - cy_bottom;
				cx += clusters[i][j]->getCenter().getX() ;
				cy += clusters[i][j]->getCenter().getY() ;
				d_l += sqrt(dx_l*dx_l+dy_l*dy_l) ;
				d_b += sqrt(dx_l*dx_l+dy_l*dy_l) ;
			}

			d_l /= count ;
			dx_l /= count ;
			dy_l /= count ;
			d_b /= count ;
			dx_b /= count ;
			dy_b /= count ;

			cx /= count ;
			cy /= count ;
			std::cout << "vector is..." << (dx_l+dx_b) << ", " << (dy_l+dy_b) << ", center " << cx << ", " << cy << std::endl ;
		}


	}

	return ret ;

}


}
