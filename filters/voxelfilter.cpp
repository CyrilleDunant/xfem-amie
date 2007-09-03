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
#include "voxelfilter.h"
#include <fstream>
#include <string>

namespace Mu {

VoxelFilter::VoxelFilter()
{
}


VoxelFilter::~VoxelFilter()
{
}


void VoxelFilter::read(const char * filename)
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
				points.push_back(new Point(10.*((double)i/r), 10.*(double)j/c,10.*(double)k/s)) ;
				(*points.rbegin())->id = index++ ;
			}
		}
	}
	std::cout << "generated points" << std::endl ;
	
	index = 0 ;


	HexahedralElement * father = new HexahedralElement(LINEAR) ;
	for( int i = 0 ; i < r/4 ; i++)
	{
		for( int j = 0 ; j < c/4 ; j++)
		{
			for( int k = 0 ; k < s/4 ; k++)
			{
				if(!file.eof())
				{
					int behaviourKey ;
					file >> behaviourKey ;
					
					std::vector<Point *> corner ;
					corner.push_back(points[k*(r+1)*(c+1)+j*(c+1)+i]) ;
					corner.push_back(points[k*(r+1)*(c+1)+j*(c+1)+i+1]) ;
					corner.push_back(points[k*(r+1)*(c+1)+(j+1)*(c+1)+i]) ;
					corner.push_back(points[k*(r+1)*(c+1)+(j+1)*(c+1)+i+1]) ;
					corner.push_back(points[(k+1)*(r+1)*(c+1)+j*(c+1)+i]) ;
					corner.push_back(points[(k+1)*(r+1)*(c+1)+j*(c+1)+i+1]) ;
					corner.push_back(points[(k+1)*(r+1)*(c+1)+(j+1)*(c+1)+i]) ;
					corner.push_back(points[(k+1)*(r+1)*(c+1)+(j+1)*(c+1)+i+1]) ;
					Hexahedron * hex = new Hexahedron(corner[0], corner[1], corner[3], corner[2], corner[4], corner[5], corner[7], corner[6]) ;
					elems.push_back(new HexahedralElement(father, hex)) ;
					LinearForm * behaviour = behaviourMap[behaviourKey] ;
					(*elems.rbegin())->setBehaviour(behaviour) ;

				}
			}
		}
	}
	
	
}

std::vector<Point *> & VoxelFilter::getPoints()
{
	return points ;
}

std::vector<HexahedralElement *> & VoxelFilter::getElements()
{
	return elems ;
}

const std::vector<Point *> & VoxelFilter::getPoints() const
{
	return points ;
}

const std::vector<HexahedralElement *> & VoxelFilter::getElements() const
{
	return elems ;
}

}
