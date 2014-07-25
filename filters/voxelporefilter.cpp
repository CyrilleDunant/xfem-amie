//
// C++ Implementation: voxelfilter
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
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

bool VoxelPoreFilter::existsPath(std::vector<std::vector<std::vector<int> > > & phase,
	                int isource, int jsource, int ksource,
	                int itarget, int jtarget, int ktarget,
			int istart, int jstart, int kstart,
	                int iend, int jend, int kend
	               ) const
{
	
	std::vector<ConnectedNode *> local ;
			
	ConnectedNode * start = nullptr;
	ConnectedNode * end = nullptr;
	for(int i =istart ; i < iend+1 ; i++)
	{
		for(int j = jstart ; j < jend+1 ; j++)
		{
			for(int k = kstart ; k < kend+1 ; k++)
			{
				if(i >-1 && j>-1 && k >-1&& 
				   i < (int)phase.size()&& 
				   j < (int)phase[i].size() && 
				   k < (int)phase[i][j].size() && 
					phase[i][j][k] == poreIndex)
				{
					local.push_back(new ConnectedNode(i, j, k)) ;
				
					if(i == isource && j == jsource && k== ksource)
						start = local[local.size()-1] ;
					
					if(i == itarget && j == jtarget && k== ktarget)
						end = local[local.size()-1] ;
				}
			}
		}
	}
	
	if(!end)
		return true ;
	
	for(size_t i = 0 ; i < local.size() ; i++)
	{
		for(size_t j = i+1 ; j < local.size() ; j++)
		{
			if(local[i]->isNeighbour(local[j]))
			{
				local[i]->neighbour.push_back(local[j]) ;
				local[j]->neighbour.push_back(local[i]) ;
			}
		}
	}
	
	std::vector<ConnectedNode *> toCheck = start->neighbour ;
	start->visited = true ;
	
	while(!toCheck.empty())
	{
		std::vector<ConnectedNode *> temp ;
		
		for(size_t i = 0 ; i < toCheck.size() ;i++)
		{
			if(!toCheck[i]->visited)
			{
				toCheck[i]->visited = true ;
				if(toCheck[i] == end)
					return true ;
				
				for(size_t j = 0 ; j < toCheck[i]->neighbour.size() ; j++)
				{
					temp.push_back(toCheck[i]->neighbour[j]) ;
				}
			}
		}
		
		toCheck = temp ;
	}

	for(size_t i = 0 ; i < local.size() ; i++)
		delete local[i] ;
	
	return false ;
	
}


void VoxelPoreFilter::read(const char * filename)
{
	points.clear() ;
	elems.clear() ;
	
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
				double xr = ((double)rand()/RAND_MAX*2.-1.)*(.05*10./r) ;
				double yr = ((double)rand()/RAND_MAX*2.-1.)*(.05*10./c) ;
				double zr = ((double)rand()/RAND_MAX*2.-1.)*(.05*10./s) ;
// 				if(i == 0 || i == r || j == 0 || j== c || k == 0 || k == s)
// 				{
					xr = 0 ;
					yr = 0 ;
					zr = 0 ;
// 				}
				
				points.push_back(new Point(7.5*((double)i/r) + xr, 7.5*(double)j/c +yr ,7.5*(double)k/s + zr, -1)) ;
				(*points.rbegin())->setId( index++) ;
			}
		}
	}
	
	for( int i = 0 ; i < r+1 ; i++)
	{
		for( int j = 0 ; j < c+1 ; j++)
		{
			for( int k = 0 ; k < s+1 ; k++)
			{
				double xr = ((double)rand()/RAND_MAX*2.-1.)*(.05*10./r) ;
				double yr = ((double)rand()/RAND_MAX*2.-1.)*(.05*10./c) ;
				double zr = ((double)rand()/RAND_MAX*2.-1.)*(.05*10./s) ;
// 				if(i == 0 || i == r || j == 0 || j== c || k == 0 || k == s)
// 				{
					xr = 0 ;
					yr = 0 ;
					zr = 0 ;
// 				}
				
				points.push_back(new Point(7.5*((double)i/r) + xr, 7.5*(double)j/c +yr ,7.5*(double)k/s + zr, 1)) ;
				(*points.rbegin())->setId( index++) ;
			}
		}
	}
	
	std::cout << "generated points" << std::endl ;
	
	index = 0 ;

	TetrahedralElement * father = new TetrahedralElement(LINEAR_TIME_LINEAR) ;
	father->compileAndPrecalculate() ;
	
	std::vector<std::vector<std::vector<int> > > phase ;
	
	for( int i = 0 ; i < r ; i++)
	{
		phase.push_back(std::vector<std::vector<int> >(0)) ;
		for( int j = 0 ; j < c ; j++)
		{
			phase[i].push_back(std::vector<int>(0)) ;
			for( int k = 0 ; k < s ; k++)
			{
				if(!file.eof())
				{
					int behaviourKey ;
					file >> behaviourKey ;
					phase[i][j].push_back(behaviourKey) ;
				}
			}
		}
	}
	std::cout << "generated phases" << std::endl ;
	
	for( int i = 0 ; i < r ; i++)
	{
		for( int j = 0 ; j < c ; j++)
		{
			for( int k = 0 ; k < s ; k++)
			{				
				
				if(phase[i][j][k] == poreIndex)
				{
					bool aaa = true ;
					
					if(!(
					      existsPath(phase, i, j, k, i-1, j-1, k-1, i-1, j-1, k-1, i, j, k) &&
					      existsPath(phase, i, j, k, i-1, j-1, k, i-1, j-1, k-1, i, j, k) &&
					      existsPath(phase, i, j, k, i-1, j, k-1, i-1, j-1, k-1, i, j, k) &&
					      existsPath(phase, i, j, k, i, j-1, k-1, i-1, j-1, k-1, i, j, k)
					    )
					  )
						aaa = false ;
					
					bool aab = true ;
					
					if(!(
					      existsPath(phase, i, j, k, i-1, j-1, k+1,i-1,j-1,k,i,j,k+1) &&
					      existsPath(phase, i, j, k, i-1, j-1, k,i-1,j-1,k,i,j,k+1) &&
					      existsPath(phase, i, j, k, i-1, j, k+1,i-1,j-1,k,i,j,k+1) &&
					      existsPath(phase, i, j, k, i, j-1, k+1,i-1,j-1,k,i,j,k+1)
					    )
					  )
						aab = false ;
					
					bool aba = true ;

					if(!(
					      existsPath(phase, i, j, k, i-1, j+1, k-1, i-1, j, k-1,i,j+1, k) &&
					      existsPath(phase, i, j, k, i-1, j+1, k,i-1, j, k-1,i,j+1, k) &&
					      existsPath(phase, i, j, k, i-1, j, k-1,i-1, j, k-1,i,j+1, k) &&
					      existsPath(phase, i, j, k, i, j+1, k-1,i-1, j, k-1,i,j+1, k)
					    )
					  )
						aba = false ;
					

					bool abb = true ;
					if(!(
					      existsPath(phase, i, j, k, i-1, j+1, k+1,i-1,j,k,i,j+1,k+1) &&
					      existsPath(phase, i, j, k, i-1, j+1, k,i-1,j,k,i,j+1,k+1) &&
					      existsPath(phase, i, j, k, i-1, j, k+1,i-1,j,k,i,j+1,k+1) &&
					      existsPath(phase, i, j, k, i, j+1, k+1,i-1,j,k,i,j+1,k+1)
					    )
					  )
						abb = false ;
					
					bool baa = true ;
					if(!(
					      existsPath(phase, i, j, k, i+1, j-1, k-1,i,j-1, k-1,i+1,j,k) &&
					      existsPath(phase, i, j, k, i+1, j-1, k,i,j-1, k-1,i+1,j,k) &&
					      existsPath(phase, i, j, k, i+1, j, k-1,i,j-1, k-1,i+1,j,k) &&
					      existsPath(phase, i, j, k, i, j-1, k-1,i,j-1, k-1,i+1,j,k)
					    )
					  )
						baa = false ;
					
					bool bab = true ;
					if(!(
					      existsPath(phase, i, j, k, i+1, j-1, k+1, i,j-1, k, i+1, j,k+1) &&
					      existsPath(phase, i, j, k, i+1, j-1, k, i,j-1, k, i+1, j,k+1) &&
					      existsPath(phase, i, j, k, i+1, j, k+1, i,j-1, k, i+1, j,k+1) &&
					      existsPath(phase, i, j, k, i, j-1, k+1, i,j-1, k, i+1, j,k+1)
					    )
					  )
						bab = false ;
					
					bool bba = true ;
					if(!(
					      existsPath(phase, i, j, k, i+1, j+1, k-1,i,j,k-1,i+1, j+1, k) &&
					      existsPath(phase, i, j, k, i+1, j+1, k,i,j,k-1,i+1, j+1, k) &&
					      existsPath(phase, i, j, k, i+1, j, k-1,i,j,k-1,i+1, j+1, k) &&
					      existsPath(phase, i, j, k, i, j+1, k-1,i,j,k-1,i+1, j+1, k)
					    )
					  )
						bba = false ;
					
					bool bbb = true ;
					if(!(
					      existsPath(phase, i, j, k, i+1, j+1, k+1,i,j,k, i+1, j+1, k+1) &&
					      existsPath(phase, i, j, k, i+1, j+1, k,i,j,k, i+1, j+1, k+1) &&
					      existsPath(phase, i, j, k, i+1, j, k+1,i,j,k, i+1, j+1, k+1) &&
					      existsPath(phase, i, j, k, i, j+1, k+1,i,j,k, i+1, j+1, k+1)
					    )
					  )
						bbb = false ;

// 				std::cout << aaa << aab << aba << abb << baa << bab << bba << bbb << std::endl ;

				std::vector<Point *> corner ;
				if(aaa)
					corner.push_back(points[i*(r+1)*(c+1)+j*(s+1)+k]) ;
				else
				{
					corner.push_back(new Point(*points[i*(r+1)*(c+1)+j*(s+1)+k])) ;
					(*corner.rbegin())->setId( index++) ;
				}
				if(aab)
					corner.push_back(points[i*(r+1)*(c+1)+j*(s+1)+k+1]) ;
				else
				{
					corner.push_back(new Point(*points[i*(r+1)*(c+1)+j*(s+1)+k+1])) ;
					(*corner.rbegin())->setId( index++) ;
				}
				if(aba)
					corner.push_back(points[i*(r+1)*(c+1)+(j+1)*(s+1)+k]) ;
				else
				{
					corner.push_back(new Point(*points[i*(r+1)*(c+1)+(j+1)*(s+1)+k])) ;
					(*corner.rbegin())->setId( index++) ;
				}
				if(abb)
					corner.push_back(points[i*(r+1)*(c+1)+(j+1)*(s+1)+k+1]) ;
				else
				{
					corner.push_back(new Point(*points[i*(r+1)*(c+1)+(j+1)*(s+1)+k+1])) ;
					(*corner.rbegin())->setId( index++) ;
				}
				if(baa)
					corner.push_back(points[(i+1)*(r+1)*(c+1)+j*(s+1)+k]) ;
				else
				{
					corner.push_back(new Point(*points[(i+1)*(r+1)*(c+1)+j*(s+1)+k])) ;
					(*corner.rbegin())->setId( index++) ;
				}
			
				if(bab)
					corner.push_back(points[(i+1)*(r+1)*(c+1)+j*(s+1)+k+1]) ;
				else
				{
					corner.push_back(new Point(*points[(i+1)*(r+1)*(c+1)+j*(s+1)+k+1])) ;
					(*corner.rbegin())->setId( index++) ;
				}
				if(bba)
					corner.push_back(points[(i+1)*(r+1)*(c+1)+(j+1)*(s+1)+k]) ;
				else
				{
					corner.push_back(new Point(*points[(i+1)*(r+1)*(c+1)+(j+1)*(s+1)+k])) ;
					(*corner.rbegin())->setId( index++) ;
				}
				if(bbb)
					corner.push_back(points[(i+1)*(r+1)*(c+1)+(j+1)*(s+1)+k+1]) ;
				else
				{
					corner.push_back(new Point(*points[(i+1)*(r+1)*(c+1)+(j+1)*(s+1)+k+1])) ;
					(*corner.rbegin())->setId( index++) ;
				}
					
				if(aaa)
					corner.push_back(points[points.size()/2 + i*(r+1)*(c+1)+j*(s+1)+k]) ;
				else
				{
					corner.push_back(new Point(*points[points.size()/2 + i*(r+1)*(c+1)+j*(s+1)+k])) ;
					(*corner.rbegin())->setId( index++) ;
				}
				if(aab)
					corner.push_back(points[points.size()/2 +i*(r+1)*(c+1)+j*(s+1)+k+1]) ;
				else
				{
					corner.push_back(new Point(*points[points.size()/2 +i*(r+1)*(c+1)+j*(s+1)+k+1])) ;
					(*corner.rbegin())->setId( index++) ;
				}
				if(aba)
					corner.push_back(points[points.size()/2 +i*(r+1)*(c+1)+(j+1)*(s+1)+k]) ;
				else
				{
					corner.push_back(new Point(*points[points.size()/2 +i*(r+1)*(c+1)+(j+1)*(s+1)+k])) ;
					(*corner.rbegin())->setId( index++) ;
				}
				if(abb)
					corner.push_back(points[points.size()/2 +i*(r+1)*(c+1)+(j+1)*(s+1)+k+1]) ;
				else
				{
					corner.push_back(new Point(*points[points.size()/2 +i*(r+1)*(c+1)+(j+1)*(s+1)+k+1])) ;
					(*corner.rbegin())->setId( index++) ;
				}
				if(baa)
					corner.push_back(points[points.size()/2 +(i+1)*(r+1)*(c+1)+j*(s+1)+k]) ;
				else
				{
					corner.push_back(new Point(*points[points.size()/2 +(i+1)*(r+1)*(c+1)+j*(s+1)+k])) ;
					(*corner.rbegin())->setId( index++) ;
				}
				if(bab)
					corner.push_back(points[points.size()/2 +(i+1)*(r+1)*(c+1)+j*(s+1)+k+1]) ;
				else
				{
					corner.push_back(new Point(*points[points.size()/2 +(i+1)*(r+1)*(c+1)+j*(s+1)+k+1])) ;
					(*corner.rbegin())->setId( index++) ;
				}
				if(bba)
					corner.push_back(points[points.size()/2 +(i+1)*(r+1)*(c+1)+(j+1)*(s+1)+k]) ;
				else
				{
					corner.push_back(new Point(*points[points.size()/2 +(i+1)*(r+1)*(c+1)+(j+1)*(s+1)+k])) ;
					(*corner.rbegin())->setId( index++) ;
				}
				if(bbb)
					corner.push_back(points[points.size()/2 +(i+1)*(r+1)*(c+1)+(j+1)*(s+1)+k+1]) ;
				else
				{
					corner.push_back(new Point(*points[points.size()/2 +(i+1)*(r+1)*(c+1)+(j+1)*(s+1)+k+1])) ;
					(*corner.rbegin())->setId( index++) ;
				}
					DelaunayTetrahedron * tet = new DelaunayTetrahedron(nullptr, nullptr, 
						corner[1], corner[5], corner[4], corner[7],
						corner[1+8], corner[5+8], corner[4+8], corner[7+8], nullptr) ;
					tet->refresh(father) ;
					tet->setOrder(QUADRATIC_TIME_QUADRATIC) ;
					elems.push_back(tet) ;
					(*elems.rbegin())->setBehaviour(behaviour) ;
					tet = new DelaunayTetrahedron(nullptr, nullptr, corner[0], corner[2], corner[3], corner[6], 
												corner[0+8], corner[2+8], corner[3+8], corner[6+8], nullptr) ;
					tet->refresh(father) ;
					tet->setOrder(QUADRATIC_TIME_QUADRATIC) ;
					elems.push_back(tet) ;
					(*elems.rbegin())->setBehaviour(behaviour) ;
					tet = new DelaunayTetrahedron( nullptr,nullptr, corner[0], corner[4], corner[6], corner[7], 
												corner[0+8], corner[4+8], corner[6+8], corner[7+8], nullptr) ;
					tet->refresh(father) ;
					tet->setOrder(QUADRATIC_TIME_QUADRATIC) ;
					elems.push_back(tet) ;
					(*elems.rbegin())->setBehaviour(behaviour) ;
					tet = new DelaunayTetrahedron( nullptr,nullptr, corner[0], corner[1], corner[4], corner[7],
												corner[0+8], corner[1+8], corner[4+8], corner[7+8], nullptr) ;
					tet->refresh(father) ;
					tet->setOrder(QUADRATIC_TIME_QUADRATIC) ;
					elems.push_back(tet) ;
					(*elems.rbegin())->setBehaviour(behaviour) ;
					tet = new DelaunayTetrahedron( nullptr,nullptr, corner[0], corner[1], corner[3], corner[7],
												corner[0+8], corner[1+8], corner[3+8], corner[7+8], nullptr) ;
					tet->refresh(father) ;
					tet->setOrder(QUADRATIC_TIME_QUADRATIC) ;
					elems.push_back(tet) ;
					(*elems.rbegin())->setBehaviour(behaviour) ;
				}
			}
			
// 			for( int k = 16 ; k < s ; k++)
// 			{
// 				if(!file.eof())
// 				{
// 					int behaviourKey ;
// 					file >> behaviourKey ;
// 				}
// 			}
			
		}
		
// 		for( int j = 16 ; j < c ; j++)
// 		{
// 			
// 			for( int k = 0 ; k < s ; k++)
// 			{
// 				if(!file.eof())
// 				{
// 					int behaviourKey ;
// 					file >> behaviourKey ;
// 				}
// 			}
// 		}
		
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

