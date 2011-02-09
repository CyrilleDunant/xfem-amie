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
#include "../../physics/stiffness.h"
#include "../../geometry/geometry_3D.h"
#include <iostream>
#include <fstream>


namespace Mu
{

VoxelWriter::VoxelWriter(std::string f, int n)  : filename(f), nVoxelX(n), nVoxelY(n), nVoxelZ(n), fullSample(true) { }

VoxelWriter::VoxelWriter(std::string f, int nx, int ny, int nz): filename(f), nVoxelX(nx), nVoxelY(ny), nVoxelZ(nz), fullSample(true) { }


VoxelWriter::VoxelWriter(std::string f, Point bl, Point tr, int n)
{
	filename = f ;
	nVoxelX = n ;
	nVoxelY = n ;
	nVoxelZ = n ;
	fullSample = false ;
	bottom_left = bl ;
	top_right = tr ;
}

VoxelWriter::VoxelWriter(std::string f, Point bl, Point tr, int nx, int ny, int nz)
{
	filename = f ;
	nVoxelX = nx ;
	nVoxelY = ny ;
	nVoxelZ = nz ;
	fullSample = false ;
	bottom_left = bl ;
	top_right = tr ;
}

void VoxelWriter::write()
{
	writeHeader() ;
	std::fstream outbin ;
	outbin.open(filename.c_str(), std::ios::out|std::ios::binary|std::ios::app) ;
	for(int i = 0 ; i < nPoints() ; i++)
	{
		for(size_t j = 0 ; j < values.size() ; j++)
		{
			unsigned char val = values[j][i] ;
			outbin.put(val) ;
		}
	}
	outbin.close() ;
}

void VoxelWriter::getField(FeatureTree * F, VWFieldType field)
{
	values.clear() ;
	std::vector<std::valarray<double> > val = getDoubleValues(F, field) ;
	for(int i = 0 ; i < numberOfFields(field) ; i++)
	{
		values.push_back(normalizeArray(val.back())) ;
		val.pop_back() ;
	}
}


std::vector<std::valarray<double> > VoxelWriter::getDoubleValues(FeatureTree * F, VWFieldType field)
{
	
	int max = nVoxelX*nVoxelY*nVoxelZ ;
	int count = 0 ;
	std::vector<std::valarray<double> > ret ;
	for(int i = 0 ; i < numberOfFields(field) ; i++)
	{
		std::valarray<double> reti(0., max) ;
		ret.push_back(reti) ;
	}
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
	std::cerr << "generating values... "<<count<<"/" << max << std::flush ;
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
				std::vector<DelaunayTetrahedron *> tris = F->getElements3D(&p) ;
				bool done = false ;
				for(size_t l = 0 ; l < tris.size() ; l++)
				{
					if(tris[l]->in(p))
					{
						std::pair<bool, std::vector<double> > val = getDoubleValue(tris[l],p,field) ;
						if(val.first)
						{
							for(int m = 0 ; m < numberOfFields(field) ; m++)
							{
								ret[m][count] = val.second[m] ;
							}
							count++ ;
							done = true ;
							break ;
						}
					}
				}

				
				if(!done)
				{
					for(int m = 0 ; m < numberOfFields(field) ; m++)
					{
						ret[m][count] = 0 ;
					}
					count++ ;
				}
				
				if(count %10000 == 0)
					std::cerr << "\rgenerating values... "<<count<<"/" << max << std::flush ;
			}
		}
	}
	std::cerr << "generating values... "<<count<<"/" << max << " ...done." << std::endl ;
	return ret ;
}

std::pair<bool,std::vector<double> > VoxelWriter::getDoubleValue(DelaunayTetrahedron * tet, const Point & p, VWFieldType field)
{
	std::vector<double> ret(numberOfFields(field)) ;
	
	bool found = false ;
	switch(field)
	{
		case VWFT_PRINCIPAL_ANGLE:
		{
			ret[0]=tet->getState().getPrincipalAngle(tet->getCenter()) ;			
			found = true ;	
			break ;
		}

		case VWFT_STIFFNESS:
		{
			Stiffness * b = dynamic_cast<Stiffness *>(tet->getBehaviour()) ;
			if(b)
			{
				ret[0]=b->getTensor(Point(0.3,0.3,0.3))[0][0] ;
				found = true ;
			}
			break ;
		}
			
		case VWFT_STRAIN:
		{
			Vector tmp = tet->getState().getStrain(p,false) ;
			for(int i = 0 ; i < 6 ; i++)
				ret[i] = tmp[5-i] ;
			found = true ;
			break ;
		}
			
		case VWFT_STRAIN_AND_STRESS:
		{
			Vector tmp = tet->getState().getStrain(p,false) ;
			for(int i = 0 ; i < 6 ; i++)
				ret[i] = tmp[5-i] ;
			tmp = tet->getState().getStress(p,false) ;
			for(int i = 0 ; i < 6 ; i++)
				ret[6+i] = tmp[5-i] ;
			found = true ;
			break ;
		}
			
		case VWFT_STRESS:
		{
			Vector tmp = tet->getState().getStress(p,false) ;
			for(int i = 0 ; i < 6 ; i++)
				ret[i] = tmp[5-i] ;
			found = true ;
			break ;
		}
			
		case VWFT_GRADIENT:
		{
			Vector tmp = tet->getState().getGradient(p,false) ;
			for(int i = 0 ; i < 3 ; i++)
				ret[i] = tmp[2-i] ;
			found = true ;
			break ;
		}
			
		case VWFT_GRADIENT_AND_FLUX:
		{
			Vector tmp = tet->getState().getGradient(p,false) ;
			for(int i = 0 ; i < 3 ; i++)
				ret[i] = tmp[2-i] ;
			tmp = tet->getState().getFlux(p,false) ;
			for(int i = 0 ; i < 3 ; i++)
				ret[3+i] = tmp[2-i] ;
			found = true ;
			break ;
		}
			
		case VWFT_FLUX:
		{
			Vector tmp = tet->getState().getFlux(p,false) ;
			for(int i = 0 ; i < 3 ; i++)
				ret[i] = tmp[2-i] ;
			found = true ;
			break ;
		}
			
		case VWFT_VON_MISES:
		{
			ret[0]=tet->getState().getMaximumVonMisesStress() ;			
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
	outstream << (int) values.size() << std::endl ;
	outstream << nVoxelX << std::endl ;
	outstream << nVoxelY << std::endl ;
	outstream << nVoxelZ << std::endl ;
	outstream.close() ;
}

void VoxelWriter::writeMap(std::string filename, FeatureTree * F, Variable axis, double pos, int n, VWFieldType field, int k, int min, int max)
{
	std::valarray<double> values((n+1)*(n+1)) ;
	
	Hexahedron * box = dynamic_cast<Hexahedron *>(F->getFeature(0)) ;
	Point c = box->getCenter() ;
	double sx = box->getXSize() ;
	double sy = box->getYSize() ;
	double sz = box->getZSize() ;

	Point origin = c - Point(sx*0.5,sy*0.5,sz*0.5) ;
	Point xlocal ;
	Point ylocal ;
	
	switch(axis)
	{
		case XI:
			origin.x = pos ;
			xlocal = Point(0.,sy/n,0.) ;
			ylocal = Point(0.,0.,sz/n) ;
			break ;

		case ETA:
			origin.y = pos ;
			xlocal = Point(0.,0.,sz/n) ;
			ylocal = Point(sx/n,0.,0.) ;
			break ;

		case ZETA:
			origin.z = pos ;
			xlocal = Point(sx/n,0.,0.) ;
			ylocal = Point(0.,sy/n,0.) ;
			break ;
	}
	
	int count = 0 ;
	VoxelWriter * dummy = new VoxelWriter("", 0) ;
	for(int i = 0 ; i < n+1 ; i++)
	{
		for(int j = 0 ; j < n+1 ; j++)
		{
			Point p(origin) ;
			p += (xlocal*(double) i) ;
			p += (ylocal*(double) j) ;
			
			std::vector<DelaunayTetrahedron *> tris = F->getElements3D(&p) ;
			bool done = false ;
			if(!tris.empty())
			{
				for(size_t l = 0 ; l < tris.size() ; l++)
				{
					if(tris[l]->in(p))
					{
						std::pair<bool, std::vector<double> > val = dummy->getDoubleValue(tris[l],p,field) ;
						if(val.first)
						{	
							values[count] = val.second[k] ;
							count++ ;
							done = true ;								
							break ;
						}
					}
				}
			}
			if(!done)
			{
				values[count] = 0 ;
				count++ ;
			}
		}
	}
	delete dummy ;
	
	std::valarray<unsigned short int> val_int = normalizeArray(values,min,max) ;
	values.resize(0) ;
	
	std::fstream outfile ;
	outfile.open(filename.c_str(), std::ios::out) ;
	count = 0 ;
	for(int i = 0 ; i < n+1 ; i++)
	{
		for(int j = 0 ; j < n+1 ; j++)
		{
			outfile << val_int[count] << " " ;
			count++ ;
		}
		outfile << std::endl ;
	}
}

std::valarray<unsigned short int> normalizeArray(const std::valarray<double> & val, unsigned short int min, unsigned short int max)
{
	double vmax = val.max() ;
	double vmin = val.min() ;
	std::valarray<unsigned short int> norm(val.size()) ;
	for(size_t i = 0 ; i < val.size() ; i++)
	{
		norm[i] = (unsigned short int) std::floor(round((double) min + (max-min)*((val[i]-vmin)/(vmax-vmin)))) ;
	}	
	return norm ;
}

int numberOfFields(VWFieldType field)
{
	switch(field)
	{
		case VWFT_PRINCIPAL_ANGLE:
			return 1 ;
		case VWFT_STIFFNESS:
			return 1 ;
		case VWFT_STRAIN:
			return 6 ;
		case VWFT_STRAIN_AND_STRESS:
			return 12 ;
		case VWFT_STRESS:
			return 6 ;
		case VWFT_GRADIENT:
			return 3 ;
		case VWFT_GRADIENT_AND_FLUX:
			return 6 ;
		case VWFT_FLUX:
			return 3 ;
		case VWFT_VON_MISES:
			return 1 ;
	}
	return 1 ;
}


} ;




























