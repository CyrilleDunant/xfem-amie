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
	std::vector<DelaunayTriangle *> tri =  source->getElements2D() ;
	nTriangles = tri.size();
	int count = 0 ;
	for(int i = 0 ; i < nTriangles ; i++)
		if(tri[i]->getBehaviour() && tri[i]->getBehaviour()->type != VOID_BEHAVIOUR)
			count++ ;
	nTriangles = count ;
	getField(TWFT_COORDINATE) ;
	getField(TWFT_DISPLACEMENTS) ;
}

void TriangleWriter::write()
{
	
	writeHeader() ;
	std::fstream outfile  ;
	outfile.open(filename.c_str(), std::ios::out|std::ios::app) ;
	for(int i = 0 ; i < nTriangles ; i++)
	{
		for(size_t j = 0 ; j < values.size() ; j++)
		{
			outfile << values[j][i] << " " ;
		}
		outfile << std::endl ;
	}
	outfile.close();
}

void TriangleWriter::getField(TWFieldType field)
{
	std::vector<std::valarray<double> > val = getDoubleValues(field) ;
	std::reverse(val.begin(),val.end());
	values.insert( values.end(), val.begin(), val.end()) ;

}

std::vector<std::valarray<double> > TriangleWriter::getDoubleValues(TWFieldType field)
{
	std::vector<std::valarray<double> > ret ;
	int iterator = 0 ;
	for(int i = 0 ; i < numberOfFields(field) ; i++)
	{
		std::valarray<double> reti(nTriangles) ;
		ret.push_back(reti) ;
	}
	if(field == TWFT_STRAIN || field == TWFT_STRAIN_AND_STRESS || field == TWFT_STRESS)
	{
		std::pair<Vector, Vector> stress_strain = source->getStressAndStrain() ;
		std::vector<DelaunayTriangle *> triangles = source->getElements2D() ;
		switch(field)
		{
			case TWFT_STRAIN:
				stress_strain.first.resize(0) ;
				for(int i = 0 ; i < triangles.size() ; i++)
				{
					if(triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR && !triangles[i]->getBehaviour()->fractured())
					{
						// epsilon11
						ret[8][iterator] = stress_strain.second[i*3*3+0] ;
						ret[7][iterator] = stress_strain.second[i*3*3+3] ;
						ret[6][iterator] = stress_strain.second[i*3*3+6] ;
					
						// epsilon12
						ret[5][iterator] = stress_strain.second[i*3*3+1] ;
						ret[4][iterator] = stress_strain.second[i*3*3+4] ;
						ret[3][iterator] = stress_strain.second[i*3*3+7] ;
					
						// epsilon22
						ret[2][iterator] = stress_strain.second[i*3*3+2] ;
						ret[1][iterator] = stress_strain.second[i*3*3+5] ;
						ret[0][iterator++] = stress_strain.second[i*3*3+8] ;
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
						ret[17][iterator] = stress_strain.second[i*3*3+0] ;
						ret[16][iterator] = stress_strain.second[i*3*3+3] ;
						ret[15][iterator] = stress_strain.second[i*3*3+6] ;
					
						// epsilon12
						ret[14][iterator] = stress_strain.second[i*3*3+1] ;
						ret[13][iterator] = stress_strain.second[i*3*3+4] ;
						ret[12][iterator] = stress_strain.second[i*3*3+7] ;
					
						// epsilon22
						ret[11][iterator] = stress_strain.second[i*3*3+2] ;
						ret[10][iterator] = stress_strain.second[i*3*3+5] ;
						ret[9][iterator] = stress_strain.second[i*3*3+8] ;

						// sigma11
						ret[8][iterator] = stress_strain.first[i*3*3+0] ;
						ret[7][iterator] = stress_strain.first[i*3*3+3] ;
						ret[6][iterator] = stress_strain.first[i*3*3+6] ;
					
						// sigma12
						ret[5][iterator] = stress_strain.first[i*3*3+1] ;
						ret[4][iterator] = stress_strain.first[i*3*3+4] ;
						ret[3][iterator] = stress_strain.first[i*3*3+7] ;
					
						// sigma22
						ret[2][iterator] = stress_strain.first[i*3*3+2] ;
						ret[1][iterator] = stress_strain.first[i*3*3+5] ;
						ret[0][iterator++] = stress_strain.first[i*3*3+8] ;
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
						ret[8][iterator] = stress_strain.first[i*3*3+0] ;
						ret[7][iterator] = stress_strain.first[i*3*3+3] ;
						ret[6][iterator] = stress_strain.first[i*3*3+6] ;
					
						// sigma12
						ret[5][iterator] = stress_strain.first[i*3*3+1] ;
						ret[4][iterator] = stress_strain.first[i*3*3+4] ;
						ret[3][iterator] = stress_strain.first[i*3*3+7] ;
					
						// sigma22
						ret[2][iterator] = stress_strain.first[i*3*3+2] ;
						ret[1][iterator] = stress_strain.first[i*3*3+5] ;
						ret[0][iterator++] = stress_strain.first[i*3*3+8] ;
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
			std::pair<Vector, Vector> gradient_flux = source->getGradientAndFlux() ;
			std::vector<DelaunayTriangle *> triangles = source->getElements2D() ;
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
				std::vector<DelaunayTriangle *> triangles = source->getElements2D() ;
				
				for(int i = 0 ; i < triangles.size() ; i++)
				{
					if(triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						ret[5][iterator] = x[triangles[i]->getBoundingPoint(0).id*2] ;
						ret[4][iterator] = x[triangles[i]->getBoundingPoint(1).id*2] ;
						ret[3][iterator] = x[triangles[i]->getBoundingPoint(2).id*2] ;
						ret[2][iterator] = x[triangles[i]->getBoundingPoint(0).id*2+1] ;
						ret[1][iterator] = x[triangles[i]->getBoundingPoint(1).id*2+1] ;
						ret[0][iterator++] = x[triangles[i]->getBoundingPoint(2).id*2+1] ;
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
				std::vector<DelaunayTriangle *> triangles = source->getElements2D() ;
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
				std::vector<DelaunayTriangle *> tri = source->getElements2D() ;
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
				ret[5]=tri->getBoundingPoint(0).x ;
				ret[4]=tri->getBoundingPoint(0).y ;
				ret[3]=tri->getBoundingPoint(1).x ;
				ret[2]=tri->getBoundingPoint(1).y ;
				ret[1]=tri->getBoundingPoint(2).x ;
				ret[0]=tri->getBoundingPoint(2).y ;
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
			case TWFT_STIFFNESS:
			{
				LinearForm * b = dynamic_cast<LinearForm *>(tri->getBehaviour()) ;
				if(b)
				{
					ret[2]=b->getTensor(Point(0.3,0.3))[0][0] ;
					ret[1]=b->getTensor(Point(0.3,0.3))[0][0] ;
					ret[0]=b->getTensor(Point(0.3,0.3))[0][0] ;
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

void TriangleWriter::writeHeader()
{
	std::fstream outstream ;
	outstream.open(filename.c_str(), std::ios::out) ;
	outstream << "TRIANGLES" << std::endl ;
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




























