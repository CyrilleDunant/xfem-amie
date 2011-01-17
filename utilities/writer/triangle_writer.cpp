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

#include "triangle_writer.h"
#include "../../physics/stiffness.h"
#include <iostream>
#include <fstream>


namespace Mu
{

TriangleWriter::TriangleWriter(std::string f, FeatureTree * F)
{
	filename = f ;
	source = F ;
	nTriangles = source->getTriangles().size() ;
	getField(TWFT_COORDINATE) ;
}

void TriangleWriter::write()
{
	writeHeader() ;
	std::fstream outfile  ;
	outfile.open(filename.c_str(), std::ios::out) ;
	for(int i = 0 ; i < nTriangles ; i++)
	{
		for(size_t j = 0 ; j < values.size() ; j++)
		{
			outfile << values[j][i] << " " ;
		}
		outfile << std::endl ;
	}
}

void TriangleWriter::getField(TWFieldType field)
{
	std::vector<std::valarray<double> > val = getDoubleValues(field) ;
	while((int) val.size() > 0)
	{
		values.push_back(val[0]) ;
		val.pop_back() ;
	}
}

std::vector<std::valarray<double> > TriangleWriter::getDoubleValues(TWFieldType field)
{
	std::vector<std::valarray<double> > ret ;
	for(int i = 0 ; i < numberOfFields(field) ; i++)
	{
		std::valarray<double> reti(nTriangles) ;
		ret.push_back(reti) ;
	}
	if(field == TWFT_STRAIN || field == TWFT_STRAIN_AND_STRESS || field == TWFT_STRESS)
	{
		std::pair<Vector, Vector> stress_strain = source->getStressAndStrain() ;
		switch(field)
		{
			case TWFT_STRAIN:
				stress_strain.first.resize(0) ;
				for(int i = 0 ; i < nTriangles ; i++)
				{
					if(!source->getTriangles()[i]->getBehaviour()->fractured())
					{
						// epsilon11
						ret[8][i] = stress_strain.second[i*3*3+0] ;
						ret[7][i] = stress_strain.second[i*3*3+3] ;
						ret[6][i] = stress_strain.second[i*3*3+6] ;
					
						// epsilon12
						ret[5][i] = stress_strain.second[i*3*3+1] ;
						ret[4][i] = stress_strain.second[i*3*3+4] ;
						ret[3][i] = stress_strain.second[i*3*3+7] ;
					
						// epsilon22
						ret[2][i] = stress_strain.second[i*3*3+2] ;
						ret[1][i] = stress_strain.second[i*3*3+5] ;
						ret[0][i] = stress_strain.second[i*3*3+8] ;
					}
					else
					{
						// epsilon11
						ret[8][i] = 0 ;
						ret[7][i] = 0 ;
						ret[6][i] = 0 ;
					
						// epsilon12
						ret[5][i] = 0 ;
						ret[4][i] = 0 ;
						ret[3][i] = 0 ;
					
						// epsilon22
						ret[2][i] = 0 ;
						ret[1][i] = 0 ;
						ret[0][i] = 0 ;					
					}
				}
				break ;
				
			case TWFT_STRAIN_AND_STRESS:
				for(int i = 0 ; i < nTriangles ; i++)
				{
					if(!source->getTriangles()[i]->getBehaviour()->fractured())
					{
						// epsilon11
						ret[17][i] = stress_strain.second[i*3*3+0] ;
						ret[16][i] = stress_strain.second[i*3*3+3] ;
						ret[15][i] = stress_strain.second[i*3*3+6] ;
					
						// epsilon12
						ret[14][i] = stress_strain.second[i*3*3+1] ;
						ret[13][i] = stress_strain.second[i*3*3+4] ;
						ret[12][i] = stress_strain.second[i*3*3+7] ;
					
						// epsilon22
						ret[11][i] = stress_strain.second[i*3*3+2] ;
						ret[10][i] = stress_strain.second[i*3*3+5] ;
						ret[9][i] = stress_strain.second[i*3*3+8] ;

						// sigma11
						ret[8][i] = stress_strain.first[i*3*3+0] ;
						ret[7][i] = stress_strain.first[i*3*3+3] ;
						ret[6][i] = stress_strain.first[i*3*3+6] ;
					
						// sigma12
						ret[5][i] = stress_strain.first[i*3*3+1] ;
						ret[4][i] = stress_strain.first[i*3*3+4] ;
						ret[3][i] = stress_strain.first[i*3*3+7] ;
					
						// sigma22
						ret[2][i] = stress_strain.first[i*3*3+2] ;
						ret[1][i] = stress_strain.first[i*3*3+5] ;
						ret[0][i] = stress_strain.first[i*3*3+8] ;
					}	
					else
					{
						// epsilon11
						ret[17][i] = 0 ;
						ret[16][i] = 0 ;
						ret[15][i] = 0 ;
					
						// epsilon12
						ret[14][i] = 0 ;
						ret[13][i] = 0 ;
						ret[12][i] = 0 ;
					
						// epsilon22
						ret[11][i] = 0 ;
						ret[10][i] = 0 ;
						ret[9][i] = 0 ;

						// sigma11
						ret[8][i] = 0 ;
						ret[7][i] = 0 ;
						ret[6][i] = 0 ;
					
						// sigma12
						ret[5][i] = 0 ;
						ret[4][i] = 0 ;
						ret[3][i] = 0 ;
					
						// sigma22
						ret[2][i] = 0 ;
						ret[1][i] = 0 ;
						ret[0][i] = 0 ;					
					}
				}
				break ;

			case TWFT_STRESS:
				stress_strain.second.resize(0) ;
				for(int i = 0 ; i < nTriangles ; i++)
				{
					if(!source->getTriangles()[i]->getBehaviour()->fractured())
					{
						// sigma11
						ret[8][i] = stress_strain.first[i*3*3+0] ;
						ret[7][i] = stress_strain.first[i*3*3+3] ;
						ret[6][i] = stress_strain.first[i*3*3+6] ;
					
						// sigma12
						ret[5][i] = stress_strain.first[i*3*3+1] ;
						ret[4][i] = stress_strain.first[i*3*3+4] ;
						ret[3][i] = stress_strain.first[i*3*3+7] ;
					
						// sigma22
						ret[2][i] = stress_strain.first[i*3*3+2] ;
						ret[1][i] = stress_strain.first[i*3*3+5] ;
						ret[0][i] = stress_strain.first[i*3*3+8] ;
					}
					else
					{
						// sigma11
						ret[8][i] = 0 ;
						ret[7][i] = 0 ;
						ret[6][i] = 0 ;
					
						// sigma12
						ret[5][i] = 0 ;
						ret[4][i] = 0 ;
						ret[3][i] = 0 ;
					
						// sigma22
						ret[2][i] = 0 ;
						ret[1][i] = 0 ;
						ret[0][i] = 0 ;										
					}
				}
				break ;
		}
		
	}
	else
	{
		if(field == TWFT_GRADIENT || field == TWFT_GRADIENT_AND_FLUX || field == TWFT_FLUX)
		{
			std::pair<Vector, Vector> gradient_flux = source->getGradientAndFlux() ;
			switch(field)
			{
				case TWFT_FLUX:
					gradient_flux.first.resize(0) ;
					for(int i = 0 ; i < nTriangles ; i++)
					{
						if(!source->getTriangles()[i]->getBehaviour()->fractured())
						{
							// j11
							ret[5][i] = gradient_flux.second[i*3*2+0] ;
							ret[4][i] = gradient_flux.second[i*3*2+2] ;
							ret[3][i] = gradient_flux.second[i*3*2+4] ;
					
							// j22
							ret[2][i] = gradient_flux.second[i*3*2+1] ;
							ret[1][i] = gradient_flux.second[i*3*2+3] ;
							ret[0][i] = gradient_flux.second[i*3*2+5] ;
						}
						else
						{
							// j11
							ret[5][i] = 0 ;
							ret[4][i] = 0 ;
							ret[3][i] = 0 ;
					
							// j22
							ret[2][i] = 0 ;
							ret[1][i] = 0 ;
							ret[0][i] = 0 ;
						}
					}
					break ;
				
				case TWFT_GRADIENT_AND_FLUX:
					for(int i = 0 ; i < nTriangles ; i++)
					{
						if(!source->getTriangles()[i]->getBehaviour()->fractured())
						{
							// d11
							ret[11][i] = gradient_flux.first[i*3*2+0] ;
							ret[10][i] = gradient_flux.first[i*3*2+2] ;
							ret[9][i] = gradient_flux.first[i*3*2+4] ;
					
							// d22
							ret[8][i] = gradient_flux.first[i*3*2+1] ;
							ret[7][i] = gradient_flux.first[i*3*2+3] ;
							ret[6][i] = gradient_flux.first[i*3*2+5] ;

							// j11
							ret[5][i] = gradient_flux.second[i*3*2+0] ;
							ret[4][i] = gradient_flux.second[i*3*2+2] ;
							ret[3][i] = gradient_flux.second[i*3*2+4] ;
					
							// j22
							ret[2][i] = gradient_flux.second[i*3*2+1] ;
							ret[1][i] = gradient_flux.second[i*3*2+3] ;
							ret[0][i] = gradient_flux.second[i*3*2+5] ;
						}
						else
						{
							// d11
							ret[11][i] = 0 ;
							ret[10][i] = 0 ;
							ret[9][i] = 0 ;
					
							// d22
							ret[8][i] = 0 ;
							ret[7][i] = 0 ;
							ret[6][i] = 0 ;
						
							// j11
							ret[5][i] = 0 ;
							ret[4][i] = 0 ;
							ret[3][i] = 0 ;
					
							// j22
							ret[2][i] = 0 ;
							ret[1][i] = 0 ;
							ret[0][i] = 0 ;
						}
					}
					break ;

				case TWFT_GRADIENT:
					gradient_flux.second.resize(0) ;
					for(int i = 0 ; i < nTriangles ; i++)
					{
						if(!source->getTriangles()[i]->getBehaviour()->fractured())
						{
							// d11
							ret[5][i] = gradient_flux.first[i*3*2+0] ;
							ret[4][i] = gradient_flux.first[i*3*2+2] ;
							ret[3][i] = gradient_flux.first[i*3*2+4] ;
					
							// d22
							ret[2][i] = gradient_flux.first[i*3*2+1] ;
							ret[1][i] = gradient_flux.first[i*3*2+3] ;
							ret[0][i] = gradient_flux.first[i*3*2+5] ;
						}
						else
						{
							// d11
							ret[5][i] = 0 ;
							ret[4][i] = 0 ;
							ret[3][i] = 0 ;
					
							// d22
							ret[2][i] = 0 ;
							ret[1][i] = 0 ;
							ret[0][i] = 0 ;
						}
					}
					break ;
			}
		
		}
		else
		{
			std::vector<DelaunayTriangle *> tri = source->getTriangles() ;
			for(size_t i = 0 ; i < tri.size() ; i++)
			{
				std::pair<bool, std::vector<double> > val = getDoubleValue(tri[i], field) ;
				if(val.first)
				{
					for(size_t j = 0 ; j < numberOfFields(field) ; j++)
						ret[j][i] = val.second[j] ;
				}
				else
				{
					for(size_t j = 0 ; j < numberOfFields(field) ; j++)
						ret[j][i] = 0 ;
				}
			}
		}
	}
}

std::pair<bool, std::vector<double> > TriangleWriter::getDoubleValue(DelaunayTriangle * tri, TWFieldType field)
{
	bool found = false ;
	std::vector<double> ret(numberOfFields(field)) ;
	switch(field)
	{
		case TWFT_COORDINATE:
			ret[5]=tri->getBoundingPoint(0).x ;
			ret[4]=tri->getBoundingPoint(0).y ;
			ret[3]=tri->getBoundingPoint(1).x ;
			ret[2]=tri->getBoundingPoint(1).y ;
			ret[1]=tri->getBoundingPoint(2).x ;
			ret[0]=tri->getBoundingPoint(2).y ;
			found = true ;
			break ;
			
		case TWFT_PRINCIPAL_ANGLE:
			ret[0] = tri->getState().getPrincipalAngle(tri->getCenter()) ;
			found = true ;
			break ;

		case TWFT_STIFFNESS:
		{
			Stiffness * b = dynamic_cast<Stiffness *>(tri->getBehaviour()) ;
			if(b)
			{
				ret[0]=b->getTensor(Point(0.3,0.3))[0][0] ;
				found = true ;
			}
			break ;
		}
			
		case TWFT_VON_MISES:
		{
			ret[0]=tri->getState().getMaximumVonMisesStress() ;		
			found = true ;	
			break ;
		}

	}
	return std::make_pair(found,ret) ;

}

void TriangleWriter::writeHeader()
{
	std::fstream outstream ;
	outstream.open(filename.c_str(), std::ios::out) ;
	outstream << "TRIANGLES" << std::endl ;
	outstream << (int) values[0].size() << std::endl ;
	outstream << 3 << std::endl ;
	outstream << (int) values.size()-3 << std::endl ;
	outstream.close() ;
}

int numberOfFields(TWFieldType field)
{
	switch(field)
	{
		case TWFT_COORDINATE:
			return 6 ;
		case TWFT_PRINCIPAL_ANGLE:
			return 1 ;
		case TWFT_STIFFNESS:
			return 1 ;
		case TWFT_STRAIN:
			return 9 ;
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
			return 1 ;
	}
	return 1 ;
}


}




























