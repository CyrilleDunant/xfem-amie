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

#include "triangle_writer.h"
#include "voxel_writer.h"
#include "../../physics/stiffness.h"
#include <iostream>
#include <fstream>


namespace Mu
{

TriangleWriter::TriangleWriter(std::string f, FeatureTree * F, int t)
{
	filename = f ;
	source = F ;
	
	if(source != NULL)
	{
		layers = source->listLayers() ;
		for(size_t j = 0 ; j < layers.size() ; j++)
		{
	    std::vector<DelaunayTriangle *> tri =  source->getElements2DInLayer(layers[j]) ;
	    nTriangles.push_back(tri.size());
	    int count = 0 ;
	    for(int i = 0 ; i < nTriangles.back() ; i++)
		    if(tri[i]->getBehaviour() && tri[i]->getBehaviour()->type != VOID_BEHAVIOUR)
			    count++ ;
	    nTriangles.back() = count ;
	    timePlane.push_back(t) ;
	    if(timePlane.back() < 0)
		    timePlane.back() = 0 ;
	    if(timePlane.back() >= tri[0]->timePlanes())
		    timePlane.back() = tri[0]->timePlanes()-1 ;
			layerTranslator[layers[j]] = j ;
			values.push_back(std::vector<std::valarray<double> >(0));

		}
	  getField(TWFT_COORDINATE, false) ;
	  getField(TWFT_DISPLACEMENTS, false) ;
	}
}

BinaryTriangleWriter::BinaryTriangleWriter(std::string f, FeatureTree * F, int t) : TriangleWriter(f, F, t) { }

void TriangleWriter::reset(FeatureTree * F, int t)
{
	values.clear() ;
	nTriangles.clear();
	layers.clear();
	timePlane.clear() ;
	extraFields.clear() ;
	layerTranslator.clear() ;
	source = F ;
	if(source != NULL)
	{
		layers = source->listLayers() ;
		for(size_t j = 0 ; j < layers.size() ; j++)
		{
	    std::vector<DelaunayTriangle *> tri =  source->getElements2DInLayer(layers[j]) ;
	    nTriangles.push_back(tri.size());
	    int count = 0 ;
	    for(int i = 0 ; i < nTriangles.back() ; i++)
		    if(tri[i]->getBehaviour() && tri[i]->getBehaviour()->type != VOID_BEHAVIOUR)
			    count++ ;
	    nTriangles.back() = count ;
	    timePlane.push_back(t) ;
	    if(timePlane.back() < 0)
		    timePlane.back() = 0 ;
	    if(timePlane.back() >= tri[0]->timePlanes())
		    timePlane.back() = tri[0]->timePlanes()-1 ;
			layerTranslator[layers[j]] = j ;
			values.push_back(std::vector<std::valarray<double> >(0));

		}
	    getField(TWFT_COORDINATE, false) ;
	    getField(TWFT_DISPLACEMENTS, false) ;
	}
}

void TriangleWriter::write()
{
	
	writeHeader(layers[0],false) ;
	std::fstream outfile  ;
	outfile.open(filename.c_str(), std::ios::out|std::ios::app) ;
	for(int i = 0 ; i < nTriangles[0] ; i++)
	{
		for(size_t j = 0 ; j < values[0].size() ; j++)
		{
			outfile << values[0][j][i] << " " ;
		}
		outfile << std::endl ;
	}
	outfile.close();
	
	std::string filename_orig = filename;
	for(size_t k = 1 ; k < layers.size() ; k++)
	{
		writeHeader(layers[k], true) ;
		outfile.open(filename_orig.c_str(), std::ios::out|std::ios::app) ;
		for(int i = 0 ; i < nTriangles[k] ; i++)
		{
			for(size_t j = 0 ; j < values[k].size() ; j++)
			{
				outfile << values[k][j][i] << " " ;
			}
			outfile << std::endl ;
		}
		outfile.close();
	}
}

void BinaryTriangleWriter::write()
{
	writeHeader(false) ;

	std::vector<std::valarray<unsigned char> > norm ;
	std::valarray<bool> voids(false, values[0][0].size()) ;
	for(size_t i = 6 ; i < values[0].size() ; i++)
	{
		std::valarray<unsigned char> n = normalizeArray(values[0][i], voids) ;
		norm.push_back(n) ;
	}

	std::fstream outbin ;
	outbin.open(filename.c_str(), std::ios::out|std::ios::app) ;
	for(int i = 0 ; i < nTriangles[0] ; i++)
	{
		for(size_t j = 0 ; j < 6 ; j++)
			outbin << values[0][j][i] << " " ;

		for(size_t j = 6 ; j < values[0].size() ; j++)
		{
			unsigned char val = norm[j-6][i] ;
			outbin << (unsigned short int) val << " " ;
		}
		outbin << std::endl ;
	}
	outbin.close();
}

void TriangleWriter::append()
{

	for(size_t k = 0 ; k < nTriangles.size() ; k++)
	{	
		writeHeader(true) ;
		std::fstream outfile  ;
		outfile.open(filename.c_str(), std::ios::out|std::ios::app) ;
		for(int i = 0 ; i < nTriangles[k] ; i++)
		{
			for(size_t j = 0 ; j < values[k].size() ; j++)
			{
				outfile << values[k][j][i] << " " ;
			}
			outfile << std::endl ;
		}
		outfile.close();
	}
	
}

void BinaryTriangleWriter::append()
{
	writeHeader(true) ;

	std::vector<std::valarray<unsigned char> > norm ;
	std::valarray<bool> voids(false, values[0].size()) ;
	for(size_t i = 6 ; i < values.size() ; i++)
	{
		std::valarray<unsigned char> n = normalizeArray(values[0][i], voids) ;
		norm.push_back(n) ;
	}

	std::fstream outbin ;
	outbin.open(filename.c_str(), std::ios::out|std::ios::app) ;
	for(int i = 0 ; i < nTriangles[0] ; i++)
	{
		for(size_t j = 0 ; j < 6 ; j++)
			outbin << values[0][j][i] << " " ;

		for(size_t j = 6 ; j < values[0].size() ; j++)
		{
			unsigned char val = norm[j-6][i] ;
			outbin << (unsigned short int) val << " " ;
		}
		outbin << std::endl ;
	}
	outbin.close();
}


void TriangleWriter::getField(TWFieldType field, bool extra)
{
	for(size_t j = 0 ; j < layers.size() ; j++)
	{
		std::vector<std::valarray<double> > val = getDoubleValues(field, layers[j]) ;
		std::reverse(val.begin(),val.end());
		values[layerTranslator[layers[j]]].insert( values[layerTranslator[layers[j]]].end(), val.begin(), val.end()) ;
	}

}

std::vector<std::valarray<double> > TriangleWriter::getDoubleValues(TWFieldType field, int layer)
{
	std::vector<std::valarray<double> > ret ;
	int iterator = 0 ;
	for(int i = 0 ; i < numberOfFields(field) ; i++)
	{
		std::valarray<double> reti(nTriangles[layerTranslator[layer]]) ;
		ret.push_back(reti) ;
	}
	if(field == TWFT_STRAIN || field == TWFT_STRAIN_AND_STRESS || field == TWFT_STRESS)
	{

		std::pair<Vector, Vector> stress_strain = source->getStressAndStrainInLayer(layer) ;
		std::vector<DelaunayTriangle *> triangles = source->getElements2DInLayer(layer) ;
		int pointsPerTri = triangles[0]->getBoundingPoints().size() ;
		int pointsPerPlane = pointsPerTri / triangles[0]->timePlanes() ;
		int factor = 1 ;
		if(triangles[0]->getBoundingPoints().size() % 6 == 0)
			factor = 2 ;

		int time_offset = timePlane[layerTranslator[layer]] * pointsPerTri / triangles[0]->timePlanes() ;

		switch(field)
		{
			case TWFT_STRAIN:
				stress_strain.first.resize(0) ;
				for(int i = 0 ; i < triangles.size() ; i++)
				{
					if(triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR && !triangles[i]->getBehaviour()->fractured())
					{
						// epsilon11
						ret[8][iterator]   = stress_strain.second[i*3*pointsPerTri+pointsPerPlane*0*factor+0+3*time_offset] ;
						ret[7][iterator]   = stress_strain.second[i*3*pointsPerTri+pointsPerPlane*1*factor+0+3*time_offset] ;
						ret[6][iterator]   = stress_strain.second[i*3*pointsPerTri+pointsPerPlane*2*factor+0+3*time_offset] ;

						// epsilon12
						ret[5][iterator]   = stress_strain.second[i*3*pointsPerTri+pointsPerPlane*0*factor+1+3*time_offset] ;
						ret[4][iterator]   = stress_strain.second[i*3*pointsPerTri+pointsPerPlane*1*factor+1+3*time_offset] ;
						ret[3][iterator]   = stress_strain.second[i*3*pointsPerTri+pointsPerPlane*2*factor+1+3*time_offset] ;

						// epsilon22
						ret[2][iterator]   = stress_strain.second[i*3*pointsPerTri+pointsPerPlane*0*factor+2+3*time_offset] ;
						ret[1][iterator]   = stress_strain.second[i*3*pointsPerTri+pointsPerPlane*1*factor+2+3*time_offset] ;
						ret[0][iterator++] = stress_strain.second[i*3*pointsPerTri+pointsPerPlane*2*factor+2+3*time_offset] ;
					}
					else if (triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						// epsilon11
						ret[8][iterator] = 0 ;
						ret[7][iterator] = 0 ;
						ret[6][iterator] = 0 ;

						// epsilon12
						ret[5][iterator] = 0 ;
						ret[4][iterator] = 0 ;
						ret[3][iterator] = 0 ;

						// epsilon22
						ret[2][iterator] = 0 ;
						ret[1][iterator] = 0 ;
						ret[0][iterator++] = 0 ;
					}
				}
				break ;

			case TWFT_STRAIN_AND_STRESS:
				for(int i = 0 ; i < triangles.size() ; i++)
				{
					if(triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR && !triangles[i]->getBehaviour()->fractured())
					{
						// epsilon11
						ret[17][iterator]   = stress_strain.second[i*3*pointsPerTri+pointsPerPlane*0*factor+0+3*time_offset] ;
						ret[16][iterator]   = stress_strain.second[i*3*pointsPerTri+pointsPerPlane*1*factor+0+3*time_offset] ;
						ret[15][iterator]   = stress_strain.second[i*3*pointsPerTri+pointsPerPlane*2*factor+0+3*time_offset] ;

						// epsilon12
						ret[14][iterator]   = stress_strain.second[i*3*pointsPerTri+pointsPerPlane*0*factor+1+3*time_offset] ;
						ret[13][iterator]   = stress_strain.second[i*3*pointsPerTri+pointsPerPlane*1*factor+1+3*time_offset] ;
						ret[12][iterator]   = stress_strain.second[i*3*pointsPerTri+pointsPerPlane*2*factor+1+3*time_offset] ;

						// epsilon22
						ret[11][iterator]   = stress_strain.second[i*3*pointsPerTri+pointsPerPlane*0*factor+2+3*time_offset] ;
						ret[10][iterator]   = stress_strain.second[i*3*pointsPerTri+pointsPerPlane*1*factor+2+3*time_offset] ;
						ret[9][iterator]    = stress_strain.second[i*3*pointsPerTri+pointsPerPlane*2*factor+2+3*time_offset] ;

						// sigma11
						ret[8][iterator]   = stress_strain.first[i*3*pointsPerTri+pointsPerPlane*0*factor+0+3*time_offset] ;
						ret[7][iterator]   = stress_strain.first[i*3*pointsPerTri+pointsPerPlane*1*factor+0+3*time_offset] ;
						ret[6][iterator]   = stress_strain.first[i*3*pointsPerTri+pointsPerPlane*2*factor+0+3*time_offset] ;

						// sigma12
						ret[5][iterator]   = stress_strain.first[i*3*pointsPerTri+pointsPerPlane*0*factor+1+3*time_offset] ;
						ret[4][iterator]   = stress_strain.first[i*3*pointsPerTri+pointsPerPlane*1*factor+1+3*time_offset] ;
						ret[3][iterator]   = stress_strain.first[i*3*pointsPerTri+pointsPerPlane*2*factor+1+3*time_offset] ;

						// sigma22
						ret[2][iterator]   = stress_strain.first[i*3*pointsPerTri+pointsPerPlane*0*factor+2+3*time_offset] ;
						ret[1][iterator]   = stress_strain.first[i*3*pointsPerTri+pointsPerPlane*1*factor+2+3*time_offset] ;
						ret[0][iterator++] = stress_strain.first[i*3*pointsPerTri+pointsPerPlane*2*factor+2+3*time_offset] ;
					}
					else if(triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						// epsilon11
						ret[17][iterator] = 0 ;
						ret[16][iterator] = 0 ;
						ret[15][iterator] = 0 ;

						// epsilon12
						ret[14][iterator] = 0 ;
						ret[13][iterator] = 0 ;
						ret[12][iterator] = 0 ;

						// epsilon22
						ret[11][iterator] = 0 ;
						ret[10][iterator] = 0 ;
						ret[9][iterator] = 0 ;

						// sigma11
						ret[8][iterator] = 0 ;
						ret[7][iterator] = 0 ;
						ret[6][iterator] = 0 ;

						// sigma12
						ret[5][iterator] = 0 ;
						ret[4][iterator] = 0 ;
						ret[3][iterator] = 0 ;

						// sigma22
						ret[2][iterator] = 0 ;
						ret[1][iterator] = 0 ;
						ret[0][iterator++] = 0 ;
					}
				}
				break ;

			case TWFT_STRESS:
				stress_strain.second.resize(0) ;
				for(int i = 0 ; i < triangles.size() ; i++)
				{
					if(triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR && !triangles[i]->getBehaviour()->fractured())
					{
						// sigma11
						ret[8][iterator]   = stress_strain.first[i*3*pointsPerTri+pointsPerPlane*0*factor+0+3*time_offset] ;
						ret[7][iterator]   = stress_strain.first[i*3*pointsPerTri+pointsPerPlane*1*factor+0+3*time_offset] ;
						ret[6][iterator]   = stress_strain.first[i*3*pointsPerTri+pointsPerPlane*2*factor+0+3*time_offset] ;

						// sigma12
						ret[5][iterator]   = stress_strain.first[i*3*pointsPerTri+pointsPerPlane*0*factor+1+3*time_offset] ;
						ret[4][iterator]   = stress_strain.first[i*3*pointsPerTri+pointsPerPlane*1*factor+1+3*time_offset] ;
						ret[3][iterator]   = stress_strain.first[i*3*pointsPerTri+pointsPerPlane*2*factor+1+3*time_offset] ;

						// sigma22
						ret[2][iterator]   = stress_strain.first[i*3*pointsPerTri+pointsPerPlane*0*factor+2+3*time_offset] ;
						ret[1][iterator]   = stress_strain.first[i*3*pointsPerTri+pointsPerPlane*1*factor+2+3*time_offset] ;
						ret[0][iterator++] = stress_strain.first[i*3*pointsPerTri+pointsPerPlane*2*factor+2+3*time_offset] ;
					}
					else if(triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						// sigma11
						ret[8][iterator] = 0 ;
						ret[7][iterator] = 0 ;
						ret[6][iterator] = 0 ;

						// sigma12
						ret[5][iterator] = 0 ;
						ret[4][iterator] = 0 ;
						ret[3][iterator] = 0 ;

						// sigma22
						ret[2][iterator] = 0 ;
						ret[1][iterator] = 0 ;
						ret[0][iterator++] = 0 ;
					}
				}
				break ;
		}
	}
	else
	{
		if(field == TWFT_GRADIENT || field == TWFT_GRADIENT_AND_FLUX || field == TWFT_FLUX)
		{
			std::pair<Vector, Vector> gradient_flux = source->getGradientAndFluxInLayer(layer) ;
			std::vector<DelaunayTriangle *> triangles = source->getElements2DInLayer(layer) ;
			switch(field)
			{

				case TWFT_FLUX:
					gradient_flux.first.resize(0) ;
					for(int i = 0 ; i < triangles.size() ; i++)
					{
						if(triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR && !triangles[i]->getBehaviour()->fractured())
						{
							// j11
							ret[5][iterator] = gradient_flux.second[i*3*2+0] ;
							ret[4][iterator] = gradient_flux.second[i*3*2+2] ;
							ret[3][iterator] = gradient_flux.second[i*3*2+4] ;

							// j22
							ret[2][iterator] = gradient_flux.second[i*3*2+1] ;
							ret[1][iterator] = gradient_flux.second[i*3*2+3] ;
							ret[0][iterator++] = gradient_flux.second[i*3*2+5] ;
						}
						else if(triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR)
						{
							// j11
							ret[5][iterator] = 0 ;
							ret[4][iterator] = 0 ;
							ret[3][iterator++] = 0 ;

							// j22
							ret[2][iterator] = 0 ;
							ret[1][iterator] = 0 ;
							ret[0][iterator++] = 0 ;
						}
					}
					break ;

				case TWFT_GRADIENT_AND_FLUX:
					for(int i = 0 ; i < triangles.size() ; i++)
					{
						if(triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR && !triangles[i]->getBehaviour()->fractured())
						{
							// d11
							ret[11][iterator] = gradient_flux.first[i*3*2+0] ;
							ret[10][iterator] = gradient_flux.first[i*3*2+2] ;
							ret[9][iterator] = gradient_flux.first[i*3*2+4] ;

							// d22
							ret[8][iterator] = gradient_flux.first[i*3*2+1] ;
							ret[7][iterator] = gradient_flux.first[i*3*2+3] ;
							ret[6][iterator] = gradient_flux.first[i*3*2+5] ;

							// j11
							ret[5][iterator] = gradient_flux.second[i*3*2+0] ;
							ret[4][iterator] = gradient_flux.second[i*3*2+2] ;
							ret[3][iterator] = gradient_flux.second[i*3*2+4] ;

							// j22
							ret[2][iterator] = gradient_flux.second[i*3*2+1] ;
							ret[1][iterator] = gradient_flux.second[i*3*2+3] ;
							ret[0][iterator++] = gradient_flux.second[i*3*2+5] ;
						}
						else if(triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR)
						{
							// d11
							ret[11][iterator] = 0 ;
							ret[10][iterator] = 0 ;
							ret[9][iterator] = 0 ;

							// d22
							ret[8][iterator] = 0 ;
							ret[7][iterator] = 0 ;
							ret[6][iterator] = 0 ;

							// j11
							ret[5][iterator] = 0 ;
							ret[4][iterator] = 0 ;
							ret[3][iterator] = 0 ;

							// j22
							ret[2][iterator] = 0 ;
							ret[1][iterator] = 0 ;
							ret[0][iterator++] = 0 ;
						}
					}
					break ;

				case TWFT_GRADIENT:
					gradient_flux.second.resize(0) ;
					for(int i = 0 ; i < triangles.size() ; i++)
					{
						if(triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR && !triangles[i]->getBehaviour()->fractured())
						{
							// d11
							ret[5][iterator] = gradient_flux.first[i*3*2+0] ;
							ret[4][iterator] = gradient_flux.first[i*3*2+2] ;
							ret[3][iterator] = gradient_flux.first[i*3*2+4] ;

							// d22
							ret[2][iterator] = gradient_flux.first[i*3*2+1] ;
							ret[1][iterator] = gradient_flux.first[i*3*2+3] ;
							ret[0][iterator++] = gradient_flux.first[i*3*2+5] ;
						}
						else if(triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR)
						{
							// d11
							ret[5][iterator] = 0 ;
							ret[4][iterator] = 0 ;
							ret[3][iterator] = 0 ;

							// d22
							ret[2][iterator] = 0 ;
							ret[1][iterator] = 0 ;
							ret[0][iterator++] = 0 ;
						}
					}
					break ;
			}

		}
		else
		{
			if(field == TWFT_DISPLACEMENTS)
			{
				Vector x = source->getDisplacements() ;
				std::vector<DelaunayTriangle *> triangles = source->getElements2DInLayer(layer) ;

				int pointsPerTri = triangles[0]->getBoundingPoints().size() ;
				int factor = 1 ;
				if(triangles[0]->getBoundingPoints().size() == 6)
					factor = 2 ;

				if(timePlane[layerTranslator[layer]] >= triangles[0]->timePlanes())
					timePlane[layerTranslator[layer]] = triangles[0]->timePlanes()-1 ;

				int time_offset = timePlane[layerTranslator[layer]] * pointsPerTri / triangles[0]->timePlanes() ;
				std::cerr << time_offset << std::endl ;

				for(int i = 0 ; i < triangles.size() ; i++)
				{
					if(triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						size_t id1 = triangles[i]->getBoundingPoint(factor*0 + time_offset).id ;
						size_t id2 = triangles[i]->getBoundingPoint(factor*1 + time_offset).id ;
						size_t id3 = triangles[i]->getBoundingPoint(factor*2 + time_offset).id ;

						ret[5][iterator] = x[id1*2] ;
						ret[4][iterator] = x[id2*2] ;
						ret[3][iterator] = x[id3*2] ;
						ret[2][iterator] = x[id1*2+1] ;
						ret[1][iterator] = x[id2*2+1] ;
						ret[0][iterator++] = x[id3*2+1] ;
					}
					else if(triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						ret[5][iterator] = 0 ;
						ret[4][iterator] = 0 ;
						ret[3][iterator] = 0 ;
						ret[2][iterator] = 0 ;
						ret[1][iterator] = 0 ;
						ret[0][iterator++] = 0 ;
					}
				}
			}
			else if(field == TWFT_DAMAGE)
			{
				std::vector<DelaunayTriangle *> triangles = source->getElements2DInLayer(layer) ;
				Vector x(triangles.size()*3) ;
				for(int i = 0 ; i < triangles.size() ; i++)
				{
					if(triangles[i]->getBehaviour() && triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR &&  triangles[i]->getBehaviour()->getDamageModel())
					{
						double d = triangles[i]->getBehaviour()->getDamageModel()->getState().max();
						if( triangles[i]->getBehaviour()->getDamageModel()->fractured())
							d = 1 ;
						ret[0][iterator] = d;
						ret[1][iterator] = d ;
						ret[2][iterator++] = d ;

					}
					else if(triangles[i]->getBehaviour() && triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						ret[2][iterator] = 0 ;
						ret[1][iterator] = 0 ;
						ret[0][iterator++] = 0 ;
					}
				}
			}
			else
			{
				std::vector<DelaunayTriangle *> tri = source->getElements2DInLayer(layer) ;
				for(size_t i = 0 ; i < tri.size() ; i++)
				{
					std::pair<bool, std::vector<double> > val = getDoubleValue(tri[i], field) ;
					if(tri[i]->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						if(val.first)
						{
							for(size_t j = 0 ; j < numberOfFields(field) ; j++)
								ret[j][iterator] = val.second[j] ;
							iterator++ ;
						}
						else
						{
							for(size_t j = 0 ; j < numberOfFields(field) ; j++)
								ret[j][iterator] = 0 ;
							iterator++ ;
						}
					}
				}
			}
		}
	}
	return ret ;
}

std::pair<bool, std::vector<double> > TriangleWriter::getDoubleValue(DelaunayTriangle * tri, TWFieldType field)
{
	bool found = false ;
	std::vector<double> ret(numberOfFields(field)) ;
	if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
	{
		switch(field)
		{
			case TWFT_COORDINATE:
				ret[5]=tri->first->x ;
				ret[4]=tri->first->y ;
				ret[3]=tri->second->x ;
				ret[2]=tri->second->y ;
				ret[1]=tri->third->x ;
				ret[0]=tri->third->y ;
				found = true ;
				break ;

			case TWFT_PRINCIPAL_ANGLE:
			{
				ret[2] = tri->getState().getPrincipalAngle(*tri->first)[0] ;
				ret[1] = tri->getState().getPrincipalAngle(*tri->second)[0] ;
				ret[0] = tri->getState().getPrincipalAngle(*tri->third)[0] ;

				found = true ;
				break ;
			}
			
			case TWFT_PRINCIPAL_STRESS:
			{
				ret[5] = tri->getState().getPrincipalStresses(*tri->first, false)[1] ;
				ret[4] = tri->getState().getPrincipalStresses(*tri->second, false)[1] ;
				ret[3] = tri->getState().getPrincipalStresses(*tri->third, false)[1] ;
				ret[2] = tri->getState().getPrincipalStresses(*tri->first, false)[0] ;
				ret[1] = tri->getState().getPrincipalStresses(*tri->second, false)[0] ;
				ret[0] = tri->getState().getPrincipalStresses(*tri->third, false)[0] ;
				found = true ;
				break ;
			}
			case TWFT_PRINCIPAL_STRAIN:
			{
				ret[5] = tri->getState().getPrincipalStrains(*tri->first, false)[1] ;
				ret[4] = tri->getState().getPrincipalStrains(*tri->second, false)[1] ;
				ret[3] = tri->getState().getPrincipalStrains(*tri->third, false)[1] ;
				ret[2] = tri->getState().getPrincipalStrains(*tri->first, false)[0] ;
				ret[1] = tri->getState().getPrincipalStrains(*tri->second, false)[0] ;
				ret[0] = tri->getState().getPrincipalStrains(*tri->third, false)[0] ;
				found = true ;
				break ;
			}
			case TWFT_STIFFNESS:
			{
				LinearForm * b = dynamic_cast<LinearForm *>(tri->getBehaviour()) ;
				double t = 0 ;
				if(tri->timePlanes() > 1)
					t = -1 + timePlane[0] * 2 / (tri->timePlanes()-1) ;

				if(b)
				{
					ret[2]=b->getTensor(Point(1,0,0,t))[0][0] ;
					ret[1]=b->getTensor(Point(0,0,0,t))[0][0] ;
					ret[0]=b->getTensor(Point(0,1,0,t))[0][0] ;
					found = true ;
				}
				break ;
			}

			case TWFT_VON_MISES:
			{
				ret[2]=tri->getState().getMaximumVonMisesStress() ;
				ret[1]=tri->getState().getMaximumVonMisesStress() ;
				ret[0]=tri->getState().getMaximumVonMisesStress() ;
				found = true ;
				break ;
			}
			
		}
	}
	return std::make_pair(found,ret) ;

}

void TriangleWriter::writeHeader(int layer, bool append)
{
	std::fstream outstream ;
	if(append)
		outstream.open(filename.c_str(), std::ios::out|std::ios::app) ;
	else
		outstream.open(filename.c_str(), std::ios::out) ;

	outstream << "TRIANGLES" << std::endl ;
	outstream << (int) nTriangles[layerTranslator[layer]] << std::endl ;
	outstream << 3 << std::endl ;
	outstream << ((int) values[layerTranslator[layer]].size()-6)/3 << std::endl ;
	outstream.close() ;
}

void BinaryTriangleWriter::writeHeader(bool append)
{
	std::fstream outstream ;
	if(append)
		outstream.open(filename.c_str(), std::ios::out|std::ios::app) ;
	else
		outstream.open(filename.c_str(), std::ios::out) ;
	outstream << "BIN_TRIANGLES" << std::endl ;
	outstream << (int) values[0].size() << std::endl ;
	outstream << 3 << std::endl ;
	outstream << ((int) values.size()-6)/3 << std::endl ;
	outstream.close() ;
}

int numberOfFields(TWFieldType field)
{
	switch(field)
	{
		case TWFT_COORDINATE:
			return 6 ;
		case TWFT_DISPLACEMENTS:
			return 6 ;
		case TWFT_PRINCIPAL_ANGLE:
			return 3 ;
		case TWFT_STIFFNESS:
			return 3 ;
		case TWFT_STRAIN:
			return 9 ;
		case TWFT_PRINCIPAL_STRAIN:
			return 6 ;
		case TWFT_PRINCIPAL_STRESS:
			return 6 ;
		case TWFT_STRAIN_AND_STRESS:
			return 18 ;
		case TWFT_STRESS:
			return 9 ;
		case TWFT_GRADIENT:
			return 6 ;
		case TWFT_GRADIENT_AND_FLUX:
			return 12 ;
		case TWFT_FLUX:
			return 6 ;
		case TWFT_VON_MISES:
			return 3 ;
		case TWFT_DAMAGE:
			return 3 ;
	}
	return 3 ;
}


}



























