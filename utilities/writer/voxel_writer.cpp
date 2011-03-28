//
// C++ Interface: voxel writer
//
// Description: 
//
//
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "voxel_writer.h"
#include "../../physics/stiffness.h"
#include "../../geometry/geometry_3D.h"
#include <iostream>
#include <fstream>
#include <sstream>


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
	std::string originalName = filename ;
	for(size_t j = 0 ; j < values.size() ; j++)
	{
		std::ostringstream newname ;
		newname << originalName << "_"<< nVoxelX << "_" << j << ".vox";
		filename = newname.str() ;
		writeHeader() ;
		std::fstream outbin ;
		outbin.open(filename.c_str(), std::ios::out|std::ios::binary|std::ios::app) ;
	// 	for(size_t j = 0 ; j < values.size() ; j++)
	// 	{
			for(int i = 0 ; i < nPoints() ; i++)
			{
				unsigned char val = values[j][i] ;
				outbin.put(val) ;
			}
	// 	}
		outbin.close() ;
	}
	filename = originalName;
}

void VoxelWriter::getField(FeatureTree * F, VWFieldType field)
{
	values.clear() ;
	std::vector<std::valarray<double> > val = getDoubleValues(F, field) ;
	for(int i = 0 ; i < numberOfFields(field) ; i++)
	{
		values.push_back(normalizeArray(val.back(), voids)) ;
		val.pop_back() ;
	}
}


std::vector<std::valarray<double> > VoxelWriter::getDoubleValues(FeatureTree * F, VWFieldType field)
{
	
	int max = nVoxelX*nVoxelY*nVoxelZ ;
	int count = 0 ;
	voids.resize(max, false);
	std::vector<std::valarray<double> > ret ;
	for(int i = 0 ; i < numberOfFields(field) ; i++)
	{
		Vector reti(-1e9, max) ;
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
	
	
	std::vector<DelaunayTetrahedron *> tris = F->getElements3D() ;
	std::cerr << "generating values ( "<<filename << " )... " << count << "/" << tris.size() << std::flush ;
	
	for(size_t t = 0 ; t < tris.size() ; t++)
	{
		double minx = tris[t]->getCircumCenter().x - tris[t]->getRadius()*1.01 ;
		double miny = tris[t]->getCircumCenter().y - tris[t]->getRadius()*1.01 ;
		double minz = tris[t]->getCircumCenter().z - tris[t]->getRadius()*1.01 ;
		
		double maxx = tris[t]->getCircumCenter().x + tris[t]->getRadius()*1.01 ;
		double maxy = tris[t]->getCircumCenter().y + tris[t]->getRadius()*1.01 ;
		double maxz = tris[t]->getCircumCenter().z + tris[t]->getRadius()*1.01 ;
		
		for(int i = nVoxelX*(minx-bottom_left.x)/((top_right.x)-(bottom_left.x)) ; i < nVoxelX*(maxx-bottom_left.x)/((top_right.x)-(bottom_left.x)) ; i++)
		{
			if(i >= 0 && i < nVoxelX)
			{
				for(int j = nVoxelY*(miny-bottom_left.y)/((top_right.y)-(bottom_left.y)) ; j < nVoxelY*(maxy-bottom_left.y)/((top_right.y)-(bottom_left.y)) ; j++)
				{
					if(j >= 0 && j < nVoxelY)
					{
						for(int k = nVoxelZ*(minz-bottom_left.z)/((top_right.z)-(bottom_left.z)) ; k < nVoxelZ*(maxz-bottom_left.z)/((top_right.z)-(bottom_left.z)) ; k++)
						{
							if(k >= 0 && k < nVoxelZ)
							{
								Point p(bottom_left) ;
								p.x += ((top_right.x)-(bottom_left.x))*((double)(i))/(double(nVoxelX-1)) ;
								p.y += ((top_right.y)-(bottom_left.y))*((double)(j))/(double(nVoxelY-1)) ;
								p.z += ((top_right.z)-(bottom_left.z))*((double)(k))/(double(nVoxelZ-1)) ;
								
								if(tris[t]->in(p) && tris[t]->getBehaviour()->type != VOID_BEHAVIOUR)
								{
									std::pair<bool, std::vector<double> > val = getDoubleValue(tris[t],p,field) ;
									std::pair<bool, std::vector<double> > valAlternate ;
									
									if(val.first)
									{
										bool hasAlternate = false ;
										for(size_t l = 0 ; l < tris[t]->neighbourhood.size() ; l++)
										{
											if( tris[t]->getNeighbourhood(l)->in(p) && tris[t]->getNeighbourhood(l)->getBehaviour()->type != VOID_BEHAVIOUR)
											{
												valAlternate = getDoubleValue(tris[t]->getNeighbourhood(l),p,field) ; 
												hasAlternate = true ;
												break ;
											}
										}
										
										for(int m = 0 ; m < numberOfFields(field) ; m++)
										{
											ret[m][k+nVoxelY*j+nVoxelY*nVoxelX*i] = val.second[m] ;
// 											if(hasAlternate)
// 											{
// 												ret[m][k+nVoxelY*j+nVoxelY*nVoxelX*i] = (ret[m][k+nVoxelY*j+nVoxelY*nVoxelX*i]+valAlternate.second[m])*.5 ;
// 											}
										}
									}
								}
								else if(tris[t]->in(p) && tris[t]->getBehaviour()->type == VOID_BEHAVIOUR)
								{
									for(int m = 0 ; m < numberOfFields(field) ; m++)
									{
										ret[m][k+nVoxelY*j+nVoxelY*nVoxelX*i] = 0 ;
										voids[k+nVoxelY*j+nVoxelY*nVoxelX*i] = true ;
									}
								}
							}
						}
					}
				}
			}
		}
		if(t %100 == 0)
					std::cerr << "\rgenerating values ( "<<filename << " )... " << t <<"/" << tris.size() << std::flush ;
	}

	std::cerr << "\rgenerating values ( "<<filename << " )... " << tris.size()<< "/" << tris.size() << " ...done." << std::endl ;
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
			ret[0]=tet->getState().getPrincipalAngle(tet->inLocalCoordinates(p))[0] ;
			found = true ;	
			break ;
		}

		case VWFT_STIFFNESS:
		{
			Matrix m = tet->getBehaviour()->getTensor(tet->inLocalCoordinates(p)) ;
			if(!m.isNull())
				ret[0] = m[0][0] ;
			else
				ret[0] = 0 ;
			
			found = true ;
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
			
		case VWFT_CONCENTRATION:
		{
			Vector tmp = tet->getState().getConcentrations(p,false) ;
			ret[0] = tmp[0] ;
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
		case VWFT_ENRICHEMENT:
		{
			ret[0]=tet->getEnrichmentFunctions().size() ;
			found = true ;	
			break ;
		}
		case VWFT_DAMAGE:
		{
			if(tet->getBehaviour()->getDamageModel())
			{
				Vector s = tet->getBehaviour()->getDamageModel()->getState() ;
				double v = std::inner_product(&s[0], &s[s.size()], &s[0], double(0)) ;
				ret[0]= v ;
				found = true ;
				
			}
			else
			{
				ret[0]= 0 ;
				found = true ;
			}
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
	outstream << 1 << std::endl ;
	outstream << nVoxelX << std::endl ;
	outstream << nVoxelY << std::endl ;
	outstream << nVoxelZ << std::endl ;
	outstream.close() ;
}

void VoxelWriter::writeMap(std::string filename, FeatureTree * F, Variable axis, double pos, int n, VWFieldType field, int k, int min, int max)
{
	std::valarray<double> vals((n+1)*(n+1)) ;
	voids.resize(false, (n+1)*(n+1));
	
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
					if(tris[l]->in(p) && tris[l]->getBehaviour()->type == VOID_BEHAVIOUR)
					{
						std::pair<bool, std::vector<double> > val = dummy->getDoubleValue(tris[l],p,field) ;
						if(val.first)
						{	
							vals[count] = val.second[k] ;
							count++ ;
							done = true ;
							break ;
						}
					}
					else if(tris[l]->in(p) && tris[l]->getBehaviour()->type == VOID_BEHAVIOUR)
					{
							vals[count] = 0 ;
							voids[count] = 0 ;
							count++ ;
							
							done = true ;
							break ;
					}
				}
			}
			if(!done)
			{
				vals[count] = 0 ;
				count++ ;
			}
		}
	}
	delete dummy ;
	
	std::valarray<unsigned char> val_int = normalizeArray(vals, voids, min,max) ;
	vals.resize(0) ;
	
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

std::valarray<unsigned char> normalizeArray(const std::valarray<double> & val, const std::valarray<bool> & voids, unsigned short int min, unsigned short int max)
{
	Vector sortedArray = val ;
	std::sort(&sortedArray[0], &sortedArray[sortedArray.size()]) ;
	double vmax = sortedArray[std::min(sortedArray.size()*1-1, (size_t)(sortedArray.size()*.99))] ;
	double vmin = sortedArray[sortedArray.size()*.01] ;
	std::valarray<unsigned char> norm(val.size()) ;
	for(size_t i = 0 ; i < val.size() ; i++)
	{
		if(!voids[i])
			norm[i] = (unsigned char) std::min(std::max(round((double) min + (double)(max-min)*((val[i]-vmin)/(vmax-vmin))), (double)min), (double)max) ;
		else
			norm[i] = 0 ;
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
		case VWFT_CONCENTRATION:
			return 1 ;
		case VWFT_GRADIENT:
			return 3 ;
		case VWFT_GRADIENT_AND_FLUX:
			return 6 ;
		case VWFT_FLUX:
			return 3 ;
		case VWFT_VON_MISES:
			return 1 ;
		default:
			return 1 ;
	}
	return 1 ;
}


} ;




























