//
// C++ Implementation: voxelfilter
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "voxelporefilter.h"
#include <fstream>
#include <string>

#include "../polynomial/vm_base.h"

using namespace Mu ;

VoxelPoreFilter::VoxelPoreFilter()
{
	poreIndex = 0 ;
}


VoxelPoreFilter::~VoxelPoreFilter()
{
}


void VoxelPoreFilter::read(const char * filename)
{
	points.clear() ;
	elems.clear() ;
	if(behaviourMap.empty())
	{
		std::cerr << "no behaviours !" << std::endl ;
		return ;
	}
	
	std::ifstream file(filename) ;
	int r ;
	int c ;
	int s ;
	
	
	std::string dummy ;
	
	if(file.is_open())
	{
		file >> dummy  ;
		file >> dummy  ;
		file >> r  ;
		file >> c  ;
		file >> s  ;
	}
	else
	{
		std::cerr << "could not open file !" << std::endl ;
		return ;
	}
	
	std::cout << "volume is " << r << " x " << c << " x " << s << std::endl ;
	int index = 0 ;
	for( int i = 0 ; i < r+1 ; i++)
	{
		for( int j = 0 ; j < c+1 ; j++)
		{
			for( int k = 0 ; k < s+1 ; k++)
			{
				double xr = ((double)random()/RAND_MAX*2.-1.)*(.05*10./r) ;
				double yr = ((double)random()/RAND_MAX*2.-1.)*(.05*10./c) ;
				double zr = ((double)random()/RAND_MAX*2.-1.)*(.05*10./s) ;
// 				if(i == 0 || i == r || j == 0 || j== c || k == 0 || k == s)
// 				{
					xr = 0 ;
					yr = 0 ;
					zr = 0 ;
// 				}
				
				points.push_back(new Point(7.5*((double)i/r) + xr, 7.5*(double)j/c +yr ,7.5*(double)k/s + zr, -1.)) ;
				(*points.rbegin())->id = index++ ;
			}
		}
	}
	
	for( int i = 0 ; i < r+1 ; i++)
	{
		for( int j = 0 ; j < c+1 ; j++)
		{
			for( int k = 0 ; k < s+1 ; k++)
			{
				double xr = ((double)random()/RAND_MAX*2.-1.)*(.05*10./r) ;
				double yr = ((double)random()/RAND_MAX*2.-1.)*(.05*10./c) ;
				double zr = ((double)random()/RAND_MAX*2.-1.)*(.05*10./s) ;
// 				if(i == 0 || i == r || j == 0 || j== c || k == 0 || k == s)
// 				{
					xr = 0 ;
					yr = 0 ;
					zr = 0 ;
// 				}
				
				points.push_back(new Point(7.5*((double)i/r) + xr, 7.5*(double)j/c +yr ,7.5*(double)k/s + zr, 1.)) ;
				(*points.rbegin())->id = index++ ;
			}
		}
	}
	
	std::cout << "generated points" << std::endl ;
	
	index = 0 ;

	TetrahedralElement * father = new TetrahedralElement(LINEAR_TIME_LINEAR) ;
	for( int i = 0 ; i < r ; i++)
	{
		for( int j = 0 ; j < c ; j++)
		{
			for( int k = 0 ; k < s ; k++)
			{
				if(!file.eof())
				{
					int behaviourKey ;
					file >> behaviourKey ;
					
					std::vector<Point *> corner ;
					corner.push_back(points[i*(r+1)*(c+1)+j*(s+1)+k]) ;
					corner.push_back(points[i*(r+1)*(c+1)+j*(s+1)+k+1]) ;
					corner.push_back(points[i*(r+1)*(c+1)+(j+1)*(s+1)+k]) ;
					corner.push_back(points[i*(r+1)*(c+1)+(j+1)*(s+1)+k+1]) ;
					corner.push_back(points[(i+1)*(r+1)*(c+1)+j*(s+1)+k]) ;
					corner.push_back(points[(i+1)*(r+1)*(c+1)+j*(s+1)+k+1]) ;
					corner.push_back(points[(i+1)*(r+1)*(c+1)+(j+1)*(s+1)+k]) ;
					corner.push_back(points[(i+1)*(r+1)*(c+1)+(j+1)*(s+1)+k+1]) ;
					corner.push_back(points[points.size()/2 + i*(r+1)*(c+1)+j*(s+1)+k]) ;
					corner.push_back(points[points.size()/2 +i*(r+1)*(c+1)+j*(s+1)+k+1]) ;
					corner.push_back(points[points.size()/2 +i*(r+1)*(c+1)+(j+1)*(s+1)+k]) ;
					corner.push_back(points[points.size()/2 +i*(r+1)*(c+1)+(j+1)*(s+1)+k+1]) ;
					corner.push_back(points[points.size()/2 +(i+1)*(r+1)*(c+1)+j*(s+1)+k]) ;
					corner.push_back(points[points.size()/2 +(i+1)*(r+1)*(c+1)+j*(s+1)+k+1]) ;
					corner.push_back(points[points.size()/2 +(i+1)*(r+1)*(c+1)+(j+1)*(s+1)+k]) ;
					corner.push_back(points[points.size()/2 +(i+1)*(r+1)*(c+1)+(j+1)*(s+1)+k+1]) ;
					// 0 1 3 2 4 5 7 6
// 					Hexahedron * hex = new Hexahedron(corner[0], corner[1], corner[2], corner[3], corner[4], corner[5], corner[6], corner[7]) ;

					if(behavio urKey == poreIndex)
					{
						DelaunayTetrahedron * tet = new DelaunayTetrahedron( NULL, 
							corner[1], corner[5], corner[4], corner[7],corner[1+8], 
							corner[5+8], corner[4+8], corner[7+8], NULL) ;
						tet->refresh(father) ;
						elems.push_back(tet) ;
						(*elems.rbegin())->setBehaviour(behaviour) ;
						tet = new DelaunayTetrahedron( NULL, corner[0], corner[2], corner[3], corner[6], 
						                               corner[0+8], corner[2+8], corner[3+8], corner[6+8], NULL) ;
						tet->refresh(father) ;
						elems.push_back(tet) ;
						(*elems.rbegin())->setBehaviour(behaviour) ;
						tet = new DelaunayTetrahedron( NULL, corner[0], corner[4], corner[6], corner[7], 
						                               corner[0+8], corner[4+8], corner[6+8], corner[7+8], NULL) ;
						tet->refresh(father) ;
						elems.push_back(tet) ;
						(*elems.rbegin())->setBehaviour(behaviour) ;
						tet = new DelaunayTetrahedron( NULL, corner[0], corner[1], corner[4], corner[7],
						                               corner[0+8], corner[1+8], corner[4+8], corner[7+8], NULL) ;
						tet->refresh(father) ;
						elems.push_back(tet) ;
						(*elems.rbegin())->setBehaviour(behaviour) ;
						tet = new DelaunayTetrahedron( NULL, corner[0], corner[1], corner[3], corner[7],
						                               corner[0+8], corner[1+8], corner[3+8], corner[7+8], NULL) ;
						tet->refresh(father) ;
						elems.push_back(tet) ;
						(*elems.rbegin())->setBehaviour(behaviour) ;
					}

				}
				
			}
			
		}
		
	}
	
	
}

std::vector<Point *> & VoxelPoreFilter::getPoints()
{
	return points ;
}

std::vector<DelaunayTetrahedron *> & VoxelPoreFilter::getElements()
{
	return elems ;
}

const std::vector<Point *> & VoxelPoreFilter::getPoints() const
{
	return points ;
}

const std::vector<DelaunayTetrahedron *> & VoxelPoreFilter::getElements() const
{
	return elems ;
}

