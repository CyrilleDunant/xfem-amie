//
// C++ Interface: voxel writer
//
// Description: 
//
//
// Author: Alain Giorla
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "voxel_writer.h"
#include "../physics/stiffness.h"
#include "../geometry/geometry_3D.h"
#include <iostream>
#include <fstream>


namespace Mu
{

VoxelWriter::VoxelWriter(std::string f, int n, int nf)
{
	filename = f ;
	nVoxelX = n ;
	nVoxelY = n ;
	nVoxelZ = n ;
	nFields = nf ;
	cFields = 0 ;
	fullSample = true ;
	writeHeader() ;
}

VoxelWriter::VoxelWriter(std::string f, int nx, int ny, int nz, int nf)
{
	filename = f ;
	nVoxelX = nx ;
	nVoxelY = ny ;
	nVoxelZ = nz ;
	nFields = nf ;
	cFields = 0 ;
	fullSample = true ;
	writeHeader() ;
}


VoxelWriter::VoxelWriter(std::string f, Point bl, Point tr, int n, int nf)
{
	filename = f ;
	nVoxelX = n ;
	nVoxelY = n ;
	nVoxelZ = n ;
	nFields = nf ;
	cFields = 0 ;
	fullSample = false ;
	bottom_left = bl ;
	top_right = tr ;
	writeHeader() ;
}

VoxelWriter::VoxelWriter(std::string f, Point bl, Point tr, int nx, int ny, int nz, int nf)
{
	filename = f ;
	nVoxelX = nx ;
	nVoxelY = ny ;
	nVoxelZ = nz ;
	nFields = nf ;
	cFields = 0 ;
	fullSample = false ;
	bottom_left = bl ;
	top_right = tr ;
	writeHeader() ;
}

void VoxelWriter::write(FeatureTree * F, VWFieldType field)
{
	std::fstream outbin ;
	outbin.open(filename.c_str(), std::ios::out|std::ios::binary|std::ios::app) ;
	if(cFields < nFields)
	{
		std::valarray<double> val_all = getDoubleValues(F, field) ;
		std::cout << val_all.size() << std::endl ;
		int i = 0 ;
		int nf = numberOfFields(field) ;
		while(cFields < nFields && i < nf)
		{
			std::valarray<double> val_double(val_all.size()/nf) ;
//			std::cout << val_double.size() << std::endl ;
			for(int j = 0 ; j < val_double.size() ; j++)
				val_double[j] = val_all[j+i*val_double.size()] ;
			std::valarray<unsigned short int> val_int = normalizeArray(val_double) ;
//			std::cout << val_int.size() << std::endl ;
			for(size_t k = 0 ; k < val_int.size() ; k++)
			{
				char val = val_int[k] ;
				outbin.put(val) ;
			}
			i++ ;
			cFields++ ;
		}
	}
	outbin.close() ;
}

std::valarray<double> VoxelWriter::getDoubleValues(FeatureTree * F, VWFieldType field)
{
	int max = nVoxelX*nVoxelY*nVoxelZ ;
	int count = 0 ;
	std::valarray<double> ret(max*numberOfFields(field)) ;
	if(fullSample)
	{
		Hexahedron * box = dynamic_cast<Hexahedron *>(F->getFeature(0)) ;
		Point c = box->getCenter() ;
		double sx = box->getXSize() ;
		double sy = box->getYSize() ;
		double sz = box->getZSize() ;
		Point vec(sx*(nVoxelX-1)/(nVoxelX),sy*(nVoxelY-1)/(nVoxelY),sz*(nVoxelZ-1)/(nVoxelZ)) ;
		vec *= 0.5 ;
		bottom_left = c-vec ;
		top_right = c+vec ;
	}
	for(int i = 0 ; i < nVoxelX ; i++)
	{
		for(int j = 0 ; j < nVoxelY ; j++)
		{
			for(int k = 0 ; k < nVoxelZ ; k++)
			{
				Point p(bottom_left) ;
				p.x += ((top_right.x)-(bottom_left.x))*((double)(i))/(double(nVoxelX-1)) ;
				p.y += ((top_right.y)-(bottom_left.y))*((double)(j))/(double(nVoxelY-1)) ;
				p.z += ((top_right.z)-(bottom_left.z))*((double)(k))/(double(nVoxelZ-1)) ;
				std::vector<DelaunayTetrahedron *> tris = F->get3DMesh()->getConflictingElements(&p) ;
				bool done = false ;
				if(!tris.empty())
				{
					for(size_t l = 0 ; l < tris.size() ; l++)
					{
						if(tris[l]->in(p))
						{
							std::pair<bool, std::vector<double> > val = getDoubleValue(tris[l],p,field) ;
							if(val.first)
							{	
								for(int m = 0 ; m < numberOfFields(field) ; m++)
								{
									ret[count*numberOfFields(field)+m] = val.second[m] ;
								}
								count++ ;
								done = true ;								
								break ;
							}
						}
					}
				}
				if(!done)
				{
					std::cout << "not done" << std::endl ;
					for(int m = 0 ; m < numberOfFields(field) ; m++)
					{
						ret[count*numberOfFields(field)+m] = 0 ;
					}
					count++ ;
				}
			}
		}
	}
	return ret ;
}

std::pair<bool,std::vector<double> > VoxelWriter::getDoubleValue(DelaunayTetrahedron * tet, const Point & p, VWFieldType field)
{
	std::vector<double> ret(numberOfFields(field)) ;
	Vector tmp ;
	bool found = false ;
	switch(field)
	{
		case F_STIFFNESS:
		{
			Stiffness * b = dynamic_cast<Stiffness *>(tet->getBehaviour()) ;
			if(b)
			{
				ret[0]=b->getTensor(Point(0.3,0.3,0.3))[0][0] ;
				found = true ;
			}
			break ;
		}
			
		case F_STRAIN:
		{
			tmp = tet->getState().getStrain(p,false) ;
			for(int i = 0 ; i < numberOfFields(field) ; i++)
				ret[i] = tmp[i] ;
			found = true ;
			break ;
		}
			
		case F_STRESS:
		{
			tmp = tet->getState().getStress(p,false) ;
			for(int i = 0 ; i < numberOfFields(field) ; i++)
				ret[i] = tmp[i] ;
			found = true ;
			break ;
		}
			
	}
	return std::make_pair(found, ret) ;
}

void VoxelWriter::writeHeader()
{
	std::fstream outstream ;
	outstream.open(filename.c_str(), std::ios::out) ;
	outstream << "VOXELS" << std::endl ;
	outstream << nFields << std::endl ;
	outstream << nVoxelX << std::endl ;
	outstream << nVoxelY << std::endl ;
	outstream << nVoxelZ << std::endl ;
	outstream.close() ;
}

std::valarray<unsigned short int> normalizeArray(std::valarray<double> val, unsigned short int min, unsigned short int max)
{
	double vmax = val.max() ;
	double vmin = val.min() ;
	std::valarray<unsigned short int> norm(val.size()) ;
	for(size_t i = 0 ; i < val.size() ; i++)
	{
		norm[i] = (unsigned short int) std::floor((double) min + (max-min)*((val[i]-vmin)/(vmax-vmin))) ;
	}	
	return norm ;
}

int numberOfFields(VWFieldType field)
{
	switch(field)
	{
		case F_STIFFNESS:
			return 1 ;
		case F_STRESS:
			return 6 ;
		case F_STRAIN:
			return 6 ;
	}
	return 1 ;
}


}




























