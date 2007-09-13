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
	
	for(size_t j = 0 ; j < father->getShapeFunctions().size()  ; j++)
	{
		father->getShapeFunction(j).compile() ;
	}
		
	
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
	int idx = 0 ;
	
	for( int i = 0 ; i < 2 ; i++)
	{
		for( int j = 0 ; j < c ; j++)
		{
			for( int k = 0 ; k < s ; k++)
			{				
				
				if(phase[i][j][k] == poreIndex)
				{
					bool aaa = true ;
					bool aab = true ;
					bool aba = true ;
					bool abb = true ;
					bool baa = true ;
					bool bab = true ;
					bool bba = true ;
					bool bbb = true ;
				
				

					if(i && j && k && phase[i-1][j-1][k-1] == poreIndex)
					{
						if(!(
						      phase[i-1][j-1][k] == poreIndex && phase[i-1][j][k]== poreIndex 
						      || phase[i-1][j-1][k] == poreIndex && phase[i][j-1][k]== poreIndex 
						      || phase[i][j-1][k-1] == poreIndex && phase[i][j-1][k]== poreIndex 
						      || phase[i][j-1][k-1] == poreIndex && phase[i][j][k-1]== poreIndex 
						      || phase[i-1][j][k-1] == poreIndex && phase[i-1][j][k]== poreIndex 
						      || phase[i-1][j][k-1] == poreIndex && phase[i][j][k-1]== poreIndex 
						    )
						  )
							aaa = false ;
					}
					
					if(i && j && k != s-1 && phase[i-1][j-1][k+1] == poreIndex)
					{
						if(!(
						      phase[i-1][j-1][k] == poreIndex && phase[i-1][j][k]== poreIndex 
						      || phase[i-1][j-1][k] == poreIndex && phase[i][j-1][k]== poreIndex 
						      || phase[i][j-1][k+1] == poreIndex && phase[i][j-1][k]== poreIndex 
						      || phase[i][j-1][k+1] == poreIndex && phase[i][j][k+1]== poreIndex 
						      || phase[i-1][j][k+1] == poreIndex && phase[i-1][j][k]== poreIndex 
						      || phase[i-1][j][k+1] == poreIndex && phase[i][j][k+1]== poreIndex 
						    )
						  )
							aab = false ;
					}
					
					if(i && j != c-1 && k &&phase[i-1][j+1][k-1] != poreIndex)
					{
						if(!(
						      phase[i-1][j+1][k] == poreIndex && phase[i-1][j][k]== poreIndex 
						      || phase[i-1][j+1][k] == poreIndex && phase[i][j+1][k]== poreIndex 
						      || phase[i][j+1][k-1] == poreIndex && phase[i][j+1][k]== poreIndex 
						      || phase[i][j+1][k-1] == poreIndex && phase[i][j][k-1]== poreIndex 
						      || phase[i-1][j][k-1] == poreIndex && phase[i-1][j][k]== poreIndex 
						      || phase[i-1][j][k-1] == poreIndex && phase[i][j][k-1]== poreIndex 
						    )
						  )
						aba = false ;
					}
					if(i && j != c-1 && k != s-1 && phase[i-1][j+1][k+1] == poreIndex)
					{
						if(!(
						      phase[i-1][j+1][k] == poreIndex && phase[i-1][j][k]== poreIndex 
						      || phase[i-1][j+1][k] == poreIndex && phase[i][j+1][k]== poreIndex 
						      || phase[i][j+1][k+1] == poreIndex && phase[i][j+1][k]== poreIndex 
						      || phase[i][j+1][k+1] == poreIndex && phase[i][j][k+1]== poreIndex 
						      || phase[i-1][j][k+1] == poreIndex && phase[i-1][j][k]== poreIndex 
						      || phase[i-1][j][k+1] == poreIndex && phase[i][j][k+1]== poreIndex 
						    )
						  )
						abb = false ;
					}
					
					if(i != r-1 && j  && k  && phase[i+1][j-1][k-1] == poreIndex)
					{
						if(!(
						      phase[i+1][j-1][k] == poreIndex && phase[i+1][j][k]== poreIndex 
						      || phase[i+1][j-1][k] == poreIndex && phase[i][j-1][k]== poreIndex 
						      || phase[i][j-1][k-1] == poreIndex && phase[i][j-1][k]== poreIndex 
						      || phase[i][j-1][k-1] == poreIndex && phase[i][j][k-1]== poreIndex 
						      || phase[i+1][j][k-1] == poreIndex && phase[i+1][j][k]== poreIndex 
						      || phase[i+1][j][k-1] == poreIndex && phase[i][j][k-1]== poreIndex 
						    )
						  )
							baa = false ;
					}
					
					if(i != r-1 && j  && k != s-1 && phase[i+1][j-1][k+1] == poreIndex)
					{
						if(!(
						      phase[i+1][j-1][k] == poreIndex && phase[i+1][j][k]== poreIndex 
						      || phase[i+1][j-1][k] == poreIndex && phase[i][j-1][k]== poreIndex 
						      || phase[i][j-1][k+1] == poreIndex && phase[i][j-1][k]== poreIndex 
						      || phase[i][j-1][k+1] == poreIndex && phase[i][j][k+1]== poreIndex 
						      || phase[i+1][j][k+1] == poreIndex && phase[i+1][j][k]== poreIndex 
						      || phase[i+1][j][k+1] == poreIndex && phase[i][j][k+1]== poreIndex 
						    )
						  )
							bab = false ;
					}
					
					if(i != r-1 && j != c-1 && k && phase[i+1][j+1][k-1] == poreIndex)
					{
						if(!(
						      phase[i+1][j+1][k] == poreIndex && phase[i+1][j][k]== poreIndex 
						      || phase[i+1][j+1][k] == poreIndex && phase[i][j+1][k]== poreIndex 
						      || phase[i][j+1][k-1] == poreIndex && phase[i][j+1][k]== poreIndex 
						      || phase[i][j+1][k-1] == poreIndex && phase[i][j][k-1]== poreIndex 
						      || phase[i+1][j][k-1] == poreIndex && phase[i+1][j][k]== poreIndex 
						      || phase[i+1][j][k-1] == poreIndex && phase[i][j][k-1]== poreIndex 
						    )
						  )
							bba = false ;
					}
					if(i != r-1 && j != c-1 && k != s-1&& phase[i+1][j+1][k+1] == poreIndex)
					{
						if(!(
						      phase[i+1][j+1][k] == poreIndex && phase[i+1][j][k]== poreIndex 
						      || phase[i+1][j+1][k] == poreIndex && phase[i][j+1][k]== poreIndex 
						      || phase[i][j+1][k+1] == poreIndex && phase[i][j+1][k]== poreIndex 
						      || phase[i][j+1][k+1] == poreIndex && phase[i][j][k+1]== poreIndex 
						      || phase[i+1][j][k+1] == poreIndex && phase[i+1][j][k]== poreIndex 
						      || phase[i+1][j][k+1] == poreIndex && phase[i][j][k+1]== poreIndex 
						    )
						  )
							bbb = false ;
					}
					
					if(j != c-1 && k != s-1&& phase[i][j+1][k+1] == poreIndex)
					{
						if(!(
						      phase[i][j][k+1] == poreIndex || phase[i][j+1][k]== poreIndex 
						    )
						  )
						{
							std::cout << "[j+1][k+1]" << std::endl ;
							abb = false ;
							bbb = false ;
						}
					}
					if(i != r-1 && k != s-1&& phase[i+1][j][k+1] == poreIndex)
					{
						if(!(
						      phase[i][j][k+1] == poreIndex || phase[i+1][j][k]== poreIndex 
						    )
						  )
						{
							bab = false ;
							bbb = false ;
						}
					}
					if(i != r-1 && j != c-1&& phase[i+1][j+1][k] == poreIndex)
					{
						if(!(
						      phase[i][j+1][k] == poreIndex || phase[i+1][j][k]== poreIndex 
						    )
						  )
						{
							bba = false ;
							bbb = false ;
						}
					}
					
					if(j && k&&  phase[i][j-1][k-1] == poreIndex)
					{
						if(!(
						      phase[i][j][k-1] == poreIndex || phase[i][j-1][k]== poreIndex 
						    )
						  )
						{
							std::cout << "[j-1][k-1]" << std::endl ;
							aaa = false ;
							baa = false ;
						}
					}
					if(i&& k&& phase[i-1][j][k-1] == poreIndex)
					{
						if(!(
						      phase[i][j][k-1] == poreIndex || phase[i-1][j][k]== poreIndex 
						    )
						  )
						{
							aaa = false ;
							aba = false ;
						}
					}
					if(i && j && phase[i-1][j-1][k] == poreIndex)
					{
						if(!(
						      phase[i][j-1][k] == poreIndex || phase[i-1][j][k]== poreIndex 
						    )
						  )
						{
							aaa = false ;
							aab = false ;
						}
					}
					
					if(k != s-1 && j &&  phase[i][j-1][k+1] == poreIndex)
					{
						if(!(
						      phase[i][j][k+1] == poreIndex || phase[i][j-1][k]== poreIndex 
						    )
						  )
						{
							std::cout << "[j-1][k+1]" << std::endl ;
							aab = false ;
							bab = false ;
						}
					}
					if( i && k != s-1 && phase[i-1][j][k+1] == poreIndex)
					{
						if(!(
						      phase[i][j][k+1] == poreIndex || phase[i-1][j][k]== poreIndex 
						    )
						  )
						{
							aab = false ;
							abb = false ;
						}
					}
					if(i && j != c-1 && phase[i-1][j+1][k] == poreIndex)
					{
						if(!(
						      phase[i][j+1][k] == poreIndex || phase[i-1][j][k]== poreIndex 
						    )
						  )
						{
							aba = false ;
							abb = false ;
						}
					}
					
					if(k && j != c-1 && phase[i][j+1][k-1] == poreIndex)
					{
						if(!(
						      phase[i][j][k-1] == poreIndex || phase[i][j+1][k] == poreIndex 
						    )
						  )
						{
							std::cout << "[j+1][k-1]" << std::endl ;
							aba = false ;
							bba = false ;
						}
					}
					if(k && i != r-1 && phase[i+1][j][k-1] == poreIndex)
					{
						if(!(
						      phase[i][j][k-1] == poreIndex || phase[i+1][j][k]== poreIndex 
						    )
						  )
						{
							baa = false ;
							bba = false ;
						}
					}
					if( i != r -1 && j && phase[i+1][j-1][k] == poreIndex)
					{
						if(!(
						      phase[i][j-1][k] == poreIndex || phase[i+1][j][k]== poreIndex 
						    )
						  )
						{
							baa = false ;
							bab = false ;
						}
					}
				
				
// 					 aaa = false ;
// 					 aab = false ;
// 					 aba = false ;
// 					 abb = false ;
// 					 baa = false ;
// 					 bab = false ;
// 					 bba = false ;
// 					 bbb = false ;
				
					
				std::vector<Point *> corner ;
				if(aaa)
					corner.push_back(points[i*(r+1)*(c+1)+j*(s+1)+k]) ;
				else
				{
					corner.push_back(new Point(*points[i*(r+1)*(c+1)+j*(s+1)+k])) ;
					(*corner.rbegin())->id = index++ ;
				}
				if(aab)
					corner.push_back(points[i*(r+1)*(c+1)+j*(s+1)+k+1]) ;
				else
				{
					corner.push_back(new Point(*points[i*(r+1)*(c+1)+j*(s+1)+k+1])) ;
					(*corner.rbegin())->id = index++ ;
				}
				if(aba)
					corner.push_back(points[i*(r+1)*(c+1)+(j+1)*(s+1)+k]) ;
				else
				{
					corner.push_back(new Point(*points[i*(r+1)*(c+1)+(j+1)*(s+1)+k])) ;
					(*corner.rbegin())->id = index++ ;
				}
				if(abb)
					corner.push_back(points[i*(r+1)*(c+1)+(j+1)*(s+1)+k+1]) ;
				else
				{
					corner.push_back(new Point(*points[i*(r+1)*(c+1)+(j+1)*(s+1)+k+1])) ;
					(*corner.rbegin())->id = index++ ;
				}
				if(baa)
					corner.push_back(points[(i+1)*(r+1)*(c+1)+j*(s+1)+k]) ;
				else
				{
					corner.push_back(new Point(*points[(i+1)*(r+1)*(c+1)+j*(s+1)+k])) ;
					(*corner.rbegin())->id = index++ ;
				}
			
				if(bab)
					corner.push_back(points[(i+1)*(r+1)*(c+1)+j*(s+1)+k+1]) ;
				else
				{
					corner.push_back(new Point(*points[(i+1)*(r+1)*(c+1)+j*(s+1)+k+1])) ;
					(*corner.rbegin())->id = index++ ;
				}
				if(bba)
					corner.push_back(points[(i+1)*(r+1)*(c+1)+(j+1)*(s+1)+k]) ;
				else
				{
					corner.push_back(new Point(*points[(i+1)*(r+1)*(c+1)+(j+1)*(s+1)+k])) ;
					(*corner.rbegin())->id = index++ ;
				}
				if(bbb)
					corner.push_back(points[(i+1)*(r+1)*(c+1)+(j+1)*(s+1)+k+1]) ;
				else
				{
					corner.push_back(new Point(*points[(i+1)*(r+1)*(c+1)+(j+1)*(s+1)+k+1])) ;
					(*corner.rbegin())->id = index++ ;
				}
					
				if(aaa)
					corner.push_back(points[points.size()/2 + i*(r+1)*(c+1)+j*(s+1)+k]) ;
				else
				{
					corner.push_back(new Point(*points[points.size()/2 + i*(r+1)*(c+1)+j*(s+1)+k])) ;
					(*corner.rbegin())->id = index++ ;
				}
				if(aab)
					corner.push_back(points[points.size()/2 +i*(r+1)*(c+1)+j*(s+1)+k+1]) ;
				else
				{
					corner.push_back(new Point(*points[points.size()/2 +i*(r+1)*(c+1)+j*(s+1)+k+1])) ;
					(*corner.rbegin())->id = index++ ;
				}
				if(aba)
					corner.push_back(points[points.size()/2 +i*(r+1)*(c+1)+(j+1)*(s+1)+k]) ;
				else
				{
					corner.push_back(new Point(*points[points.size()/2 +i*(r+1)*(c+1)+(j+1)*(s+1)+k])) ;
					(*corner.rbegin())->id = index++ ;
				}
				if(abb)
					corner.push_back(points[points.size()/2 +i*(r+1)*(c+1)+(j+1)*(s+1)+k+1]) ;
				else
				{
					corner.push_back(new Point(*points[points.size()/2 +i*(r+1)*(c+1)+(j+1)*(s+1)+k+1])) ;
					(*corner.rbegin())->id = index++ ;
				}
				if(baa)
					corner.push_back(points[points.size()/2 +(i+1)*(r+1)*(c+1)+j*(s+1)+k]) ;
				else
				{
					corner.push_back(new Point(*points[points.size()/2 +(i+1)*(r+1)*(c+1)+j*(s+1)+k])) ;
					(*corner.rbegin())->id = index++ ;
				}
				if(bab)
					corner.push_back(points[points.size()/2 +(i+1)*(r+1)*(c+1)+j*(s+1)+k+1]) ;
				else
				{
					corner.push_back(new Point(*points[points.size()/2 +(i+1)*(r+1)*(c+1)+j*(s+1)+k+1])) ;
					(*corner.rbegin())->id = index++ ;
				}
				if(bba)
					corner.push_back(points[points.size()/2 +(i+1)*(r+1)*(c+1)+(j+1)*(s+1)+k]) ;
				else
				{
					corner.push_back(new Point(*points[points.size()/2 +(i+1)*(r+1)*(c+1)+(j+1)*(s+1)+k])) ;
					(*corner.rbegin())->id = index++ ;
				}
				if(bbb)
					corner.push_back(points[points.size()/2 +(i+1)*(r+1)*(c+1)+(j+1)*(s+1)+k+1]) ;
				else
				{
					corner.push_back(new Point(*points[points.size()/2 +(i+1)*(r+1)*(c+1)+(j+1)*(s+1)+k+1])) ;
					(*corner.rbegin())->id = index++ ;
				}
					// 0 1 3 2 4 5 7 6
// 					Hexahedron * hex = new Hexahedron(corner[0], corner[1], corner[2], corner[3], corner[4], corner[5], corner[6], corner[7]) ;


					DelaunayTetrahedron * tet = new DelaunayTetrahedron( NULL, 
						corner[1], corner[5], corner[4], corner[7],
						corner[1+8], corner[5+8], corner[4+8], corner[7+8], NULL) ;
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
				
				idx++ ;
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

