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
#include "voxelfilter.h"
#include <fstream>
#include <string>

#include "../polynomial/vm_base.h"

using namespace Amie ;

VoxelFilter::VoxelFilter()
{
}


VoxelFilter::~VoxelFilter()
{
}

bool VoxelFilter::existsPath(std::vector<std::vector<std::vector<unsigned char> > > & phase,
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
					phase[i][j][k])
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



void VoxelFilter::read(const char * filename)
{
	points.clear() ;
	elems.clear() ;
	if(behaviourMap.empty())
	{
		std::cerr << "no behaviours !" << std::endl ;
		return ;
	}
	
	std::fstream file(filename) ;
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

	long position = file.tellp() ;
	file.close() ;
	file.open(filename, std::ios::binary|std::ios::in) ;
	file.seekp(position)  ;

	std::cout << "volume is " << r << " x " << c << " x " << s << std::endl ;
	int index = 0 ;
	for( int i = 0 ; i < r+1 ; i++)
	{
		for( int j = 0 ; j < c+1 ; j++)
		{
			for( int k = 0 ; k < s+1 ; k++)
			{

				points.push_back(new Point(100.*((double)i/r), 100.*(double)j/c  ,100.*(double)k/s)) ;
				(*points.rbegin())->setId(index++) ;
			}
		}
	}
	std::cerr << "generated points" << std::endl ;
	
	index = 0 ;

	TetrahedralElement * father = new TetrahedralElement(LINEAR) ;

	father->compileAndPrecalculate() ;
	
	std::vector<std::vector<std::vector<unsigned char> > > phase ;
	std::vector<std::vector<std::vector<std::pair<bool, bool> > > > connected_visited ;
	std::vector< std::valarray<size_t> > connected_to_check ;
	for( int i = 0 ; i < r ; i++)
	{
		phase.push_back(std::vector<std::vector<unsigned char> >(0)) ;
		connected_visited.push_back(std::vector<std::vector<std::pair<bool, bool> > >(0)) ;
		for( int j = 0 ; j < c ; j++)
		{
			phase[i].push_back(std::vector<unsigned char>(0)) ;
			connected_visited[i].push_back(std::vector<std::pair<bool, bool> >(0)) ;
			for( int k = 0 ; k < s ; k++)
			{
				if(!file.eof())
				{
					unsigned char behaviourKey ;
					file >> behaviourKey ;
					phase[i][j].push_back(behaviourKey) ;
					connected_visited[i][j].push_back(std::make_pair((i == 0 
						|| i == r-1 
						|| j == 0 
						|| j == c-1 
						|| k == 0 
						|| k == s-1) && behaviourMap[behaviourKey]->type != VOID_BEHAVIOUR, (i == 0 
						               || i == r-1 
						               || j == 0 
						               || j == c-1 
						               || k == 0 
						               || k == s-1))) ;
					if(connected_visited[i][j][k].first)
					{
						std::valarray<size_t> coord(3) ;
						coord[0] = i ;
						coord[1] = j ;
						coord[2] = k ;
						
						connected_to_check.push_back(coord) ;
					}
				}
			}
		}
	}
	
	while(!connected_to_check.empty())
	{
		std::vector< std::valarray<size_t> > connected_to_check_temp ;
		
		for(size_t i = 0 ; i < connected_to_check.size() ; i++)
		{
			if(connected_to_check[i][0] > 0 
			   && !connected_visited
			   [connected_to_check[i][0]-1]
			   [connected_to_check[i][1]]
			   [connected_to_check[i][2]].second 
			   && !connected_visited
			   [connected_to_check[i][0]-1]
			   [connected_to_check[i][1]]
			   [connected_to_check[i][2]].first 
			   && !behaviourMap[phase[connected_to_check[i][0]-1]
			   [connected_to_check[i][1]]
			   [connected_to_check[i][2]]]->type != VOID_BEHAVIOUR
			  )
			{
				connected_visited
				[connected_to_check[i][0]-1]
				[connected_to_check[i][1]]
				[connected_to_check[i][2]].first = true ;
				connected_visited
				[connected_to_check[i][0]-1]
				[connected_to_check[i][1]]
				[connected_to_check[i][2]].second = true ;
				
				std::valarray<size_t> coord(3) ;
				coord[0] = connected_to_check[i][0]-1 ;
				coord[1] = connected_to_check[i][1] ;
				coord[2] = connected_to_check[i][2] ;
				
				connected_to_check_temp.push_back(coord) ;
			}
			if(connected_to_check[i][1] > 0 
			   && !connected_visited
			   [connected_to_check[i][0]]
			   [connected_to_check[i][1]-1]
			   [connected_to_check[i][2]].second 
			   && !connected_visited
			   [connected_to_check[i][0]]
			   [connected_to_check[i][1]-1]
			   [connected_to_check[i][2]].first 
			   && !behaviourMap[phase[connected_to_check[i][0]]
			                    [connected_to_check[i][1]-1]
			                    [connected_to_check[i][2]]]->type != VOID_BEHAVIOUR
			  )
			{
				connected_visited
					[connected_to_check[i][0]]
					[connected_to_check[i][1]-1]
					[connected_to_check[i][2]].first = true ;
				connected_visited
					[connected_to_check[i][0]]
					[connected_to_check[i][1]-1]
					[connected_to_check[i][2]].second = true ;
				
				std::valarray<size_t> coord(3) ;
				coord[0] = connected_to_check[i][0] ;
				coord[1] = connected_to_check[i][1]-1 ;
				coord[2] = connected_to_check[i][2] ;
				
				connected_to_check_temp.push_back(coord) ;
			}
			if(connected_to_check[i][2] > 0 
			   && !connected_visited
			   [connected_to_check[i][0]]
			   [connected_to_check[i][1]]
			   [connected_to_check[i][2]-1].second 
			   && !connected_visited
			   [connected_to_check[i][0]]
			   [connected_to_check[i][1]]
			   [connected_to_check[i][2]-1].first 
			   && !behaviourMap[phase[connected_to_check[i][0]]
			                    [connected_to_check[i][1]]
			                    [connected_to_check[i][2]]-1]->type != VOID_BEHAVIOUR
			  )
			{
				connected_visited
					[connected_to_check[i][0]]
					[connected_to_check[i][1]]
					[connected_to_check[i][2]-1].first = true ;
				connected_visited
					[connected_to_check[i][0]]
					[connected_to_check[i][1]]
					[connected_to_check[i][2]-1].second = true ;
				
				std::valarray<size_t> coord(3) ;
				coord[0] = connected_to_check[i][0] ;
				coord[1] = connected_to_check[i][1] ;
				coord[2] = connected_to_check[i][2]-1 ;
				
				connected_to_check_temp.push_back(coord) ;
			}
			if(connected_to_check[i][0] < r-1 
			   && !connected_visited
			   [connected_to_check[i][0]+1]
			   [connected_to_check[i][1]]
			   [connected_to_check[i][2]].second 
			   && !connected_visited
			   [connected_to_check[i][0]+1]
			   [connected_to_check[i][1]]
			   [connected_to_check[i][2]].first 
			   && !behaviourMap[phase[connected_to_check[i][0]+1]
			                    [connected_to_check[i][1]]
			                    [connected_to_check[i][2]]]->type != VOID_BEHAVIOUR
			  )
			{
				connected_visited
					[connected_to_check[i][0]+1]
					[connected_to_check[i][1]]
					[connected_to_check[i][2]].first = true ;
				connected_visited
					[connected_to_check[i][0]+1]
					[connected_to_check[i][1]]
					[connected_to_check[i][2]].second = true ;
				
				std::valarray<size_t> coord(3) ;
				coord[0] = connected_to_check[i][0]+1 ;
				coord[1] = connected_to_check[i][1] ;
				coord[2] = connected_to_check[i][2] ;
				
				connected_to_check_temp.push_back(coord) ;
			}
			if(connected_to_check[i][1] < c-1 
			   && !connected_visited
			   [connected_to_check[i][0]]
			   [connected_to_check[i][1]+1]
			   [connected_to_check[i][2]].second 
			   && !connected_visited
			   [connected_to_check[i][0]]
			   [connected_to_check[i][1]+1]
			   [connected_to_check[i][2]].first 
			   && !behaviourMap[phase[connected_to_check[i][0]]
			                    [connected_to_check[i][1]+1]
			                    [connected_to_check[i][2]]]->type != VOID_BEHAVIOUR
			  )
			{
				connected_visited
					[connected_to_check[i][0]]
					[connected_to_check[i][1]+1]
					[connected_to_check[i][2]].first = true ;
				connected_visited
					[connected_to_check[i][0]]
					[connected_to_check[i][1]+1]
					[connected_to_check[i][2]].second = true ;
				
				std::valarray<size_t> coord(3) ;
				coord[0] = connected_to_check[i][0] ;
				coord[1] = connected_to_check[i][1]+1 ;
				coord[2] = connected_to_check[i][2] ;
				
				connected_to_check_temp.push_back(coord) ;
			}
			if(connected_to_check[i][2] < s-1 
			   && !connected_visited
			   [connected_to_check[i][0]]
			   [connected_to_check[i][1]]
			   [connected_to_check[i][2]+1].second 
			   && !connected_visited
			   [connected_to_check[i][0]]
			   [connected_to_check[i][1]]
			   [connected_to_check[i][2]+1].first 
			   && !behaviourMap[phase[connected_to_check[i][0]]
			                    [connected_to_check[i][1]]
			                    [connected_to_check[i][2]]+1]->type != VOID_BEHAVIOUR
			  )
			{
				connected_visited
					[connected_to_check[i][0]]
					[connected_to_check[i][1]]
					[connected_to_check[i][2]+1].first = true ;
				connected_visited
					[connected_to_check[i][0]]
					[connected_to_check[i][1]]
					[connected_to_check[i][2]+1].second = true ;
				
				std::valarray<size_t> coord(3) ;
				coord[0] = connected_to_check[i][0] ;
				coord[1] = connected_to_check[i][1] ;
				coord[2] = connected_to_check[i][2]+1 ;
				
				connected_to_check_temp.push_back(coord) ;
			}
			
			
		}
		
		connected_to_check = connected_to_check_temp ;
	}
	std::cout << "generated phases" << std::endl ;
	int eindex = 0 ;
	for( int i = 0 ; i < r ; i++)
	{
		for( int j = 0 ; j < c ; j++)
		{
			for( int k = 0 ; k < s ; k++)
			{
				if(behaviourMap[phase[i][j][k]]->type != VOID_BEHAVIOUR && connected_visited[i][j][k].first)
				{
					std::vector<Point *> corner ;

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
					
					
					DelaunayTetrahedron * tet = new DelaunayTetrahedron(nullptr, nullptr, corner[1], corner[5], corner[4], corner[7], nullptr) ;
					tet->refresh(father) ;
					elems.push_back(tet) ;
                    elems.back()->index = eindex++ ;
					LinearForm * behaviour = behaviourMap[phase[i][j][k]] ;
					elems.back()->setBehaviour(nullptr,behaviour) ;
					tet = new DelaunayTetrahedron( nullptr,nullptr, corner[0], corner[2], corner[3], corner[6], nullptr) ;
					tet->refresh(father) ;
					elems.push_back(tet) ;
                    elems.back()->index = eindex++ ;
					elems.back()->setBehaviour(nullptr,behaviour) ;
					tet = new DelaunayTetrahedron( nullptr,nullptr, corner[0], corner[4], corner[6], corner[7], nullptr) ;
					tet->refresh(father) ;
					elems.push_back(tet) ;
                    elems.back()->index = eindex++ ;
					elems.back()->setBehaviour(nullptr,behaviour) ;
					tet = new DelaunayTetrahedron( nullptr,nullptr, corner[0], corner[1], corner[4], corner[7], nullptr) ;
					tet->refresh(father) ;
					elems.push_back(tet) ;
                    elems.back()->index = eindex++ ;
					elems.back()->setBehaviour(nullptr,behaviour) ;
					tet = new DelaunayTetrahedron( nullptr,nullptr, corner[0], corner[1], corner[3], corner[7], nullptr) ;
					tet->refresh(father) ;
					elems.push_back(tet) ;
                    elems.back()->index = eindex++ ;
					elems.back()->setBehaviour(nullptr,behaviour) ;
				}
            }
		}
	}
	
	for( int i = 0 ; i < r ; i++)
    {
        for( int j = 0 ; j < c ; j++)
        {
            for( int k = 0 ; k < s ; k++)
            {
                int elemIndex = (j*r+k*r*c+i)*5 ;
                for(int delta = 0 ; delta < 5 ; delta++ )
                {
                    for(int l = std::max(i-1, 0) ; l <= std::min(i+1, r-1) ; l++)
                    {
                        for(int m = std::max(j-1, 0) ; m <= std::min(j+1, c-1) ; m++)
                        {
                            for(int n = std::max(k-1, 0) ; n <= std::min(k+1, s-1) ; n++)
                            {
                                int auxElemIndex((n*r+m*r*c+l)*5) ;
                                
                                for(int auxdelta = 0 ; auxdelta < 5 ; auxdelta++ )
                                {
                                    if(elems[elemIndex+delta]->isNeighbour(elems[auxElemIndex+auxdelta]))
                                        elems[elemIndex+delta]->addNeighbour(elems[auxElemIndex+auxdelta]) ;
                                    if(elems[elemIndex+delta]->numberOfCommonVertices(elems[auxElemIndex+auxdelta]))
                                        elems[elemIndex+delta]->addNeighbourhood(elems[auxElemIndex+auxdelta]) ;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

std::vector<Point *> & VoxelFilter::getPoints()
{
	return points ;
}

std::vector<DelaunayTetrahedron *> & VoxelFilter::getElements()
{
	return elems ;
}

const std::vector<Point *> & VoxelFilter::getPoints() const
{
	return points ;
}

const std::vector<DelaunayTetrahedron *> & VoxelFilter::getElements() const
{
	return elems ;
}

