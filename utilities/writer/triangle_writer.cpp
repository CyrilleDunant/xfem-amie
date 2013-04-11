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

TriangleWriter::TriangleWriter( std::string f, FeatureTree *F, int t )
{
	filename = f ;
	source = F ;

	if( source != nullptr )
	{
		layers = source->listLayers() ;

		for( size_t j = 0 ; j < layers.size() ; j++ )
		{
			std::vector<DelaunayTriangle *> tri =  source->getElements2DInLayer( layers[j] ) ;
			nTriangles.push_back( tri.size() );
/*			int count = 0 ;

			for( int i = 0 ; i < nTriangles.back() ; i++ )
				if( tri[i]->getBehaviour() && tri[i]->getBehaviour()->type != VOID_BEHAVIOUR )
					count++ ;

			nTriangles.back() = count ;*/
			timePlane.push_back( t ) ;

			if( timePlane.back() < 0 )
				timePlane.back() = 0 ;

			if( timePlane.back() >= tri[0]->timePlanes() )
				timePlane.back() = tri[0]->timePlanes() - 1 ;

			layerTranslator[layers[j]] = j ;
			values.push_back( std::vector<std::valarray<double> >( 0 ) );

		}

		getField( TWFT_COORDINATE, false ) ;
		getField( TWFT_DISPLACEMENTS, false ) ;
	}
}

BinaryTriangleWriter::BinaryTriangleWriter( std::string f, FeatureTree *F, int t ) : TriangleWriter( f, F, t ) { }

MultiTriangleWriter::MultiTriangleWriter(std::string h, std::string b, FeatureTree *F, int t) : TriangleWriter(b, F, t)
{
    counter = 0 ;
    head = h ;
    base = b ;

    std::fstream headstream ;

    headstream.open( head.c_str(), std::ios::out ) ;
    headstream << "MULTI" << std::endl ;
    headstream.close();

}

void TriangleWriter::reset( FeatureTree *F, int t )
{
	values.clear() ;
	nTriangles.clear();
	layers.clear();
	timePlane.clear() ;
	extraFields.clear() ;
	layerTranslator.clear() ;
	source = F ;
	if( source  )
	{
		layers = source->listLayers() ;
		for( size_t j = 0 ; j < layers.size() ; j++ )
		{
			std::vector<DelaunayTriangle *> tri =  source->getElements2DInLayer( layers[j] ) ;
			if(tri.empty())
				continue ;
			nTriangles.push_back( tri.size() );
			int count = 0 ;

			for( int i = 0 ; i < nTriangles.back() ; i++ )
				if( tri[i]->getBehaviour() && tri[i]->getBehaviour()->type != VOID_BEHAVIOUR )
					count++ ;

			nTriangles.back() = count ;
			timePlane.push_back( t ) ;

			if( timePlane.back() < 0 )
				timePlane.back() = 0 ;

			if( timePlane.back() >= tri[0]->timePlanes() )
				timePlane.back() = tri[0]->timePlanes() - 1 ;

			layerTranslator[layers[j]] = j ;
			values.push_back( std::vector<std::valarray<double> >( 0 ) );

		}
		getField( TWFT_COORDINATE, false ) ;
		getField( TWFT_DISPLACEMENTS, false ) ;
	}
}

void TriangleWriter::write()
{

	writeHeader( layers[0], false ) ;
	std::fstream outfile  ;
	outfile.open( filename.c_str(), std::ios::out | std::ios::app ) ;

	for( int i = 0 ; i < nTriangles[0] ; i++ )
	{
		for( size_t j = 0 ; j < values[0].size() ; j++ )
		{
			outfile << values[0][j][i] << " " ;
		}

		outfile << std::endl ;
	}

	outfile.close();

	std::string filename_orig = filename;

	for( size_t k = 1 ; k < layers.size() ; k++ )
	{
		writeHeader( layers[k], true ) ;
		outfile.open( filename_orig.c_str(), std::ios::out | std::ios::app ) ;

		for( int i = 0 ; i < nTriangles[k] ; i++ )
		{
			for( size_t j = 0 ; j < values[k].size() ; j++ )
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
	writeHeader( false ) ;

	std::vector<std::valarray<unsigned char> > norm ;
	std::valarray<bool> voids( false, values[0][0].size() ) ;

	for( size_t i = 6 ; i < values[0].size() ; i++ )
	{
		std::valarray<unsigned char> n = normalizeArray( values[0][i], voids ) ;
		norm.push_back( n ) ;
	}

	std::fstream outbin ;
	outbin.open( filename.c_str(), std::ios::out | std::ios::app ) ;

	for( int i = 0 ; i < nTriangles[0] ; i++ )
	{
		for( size_t j = 0 ; j < 6 ; j++ )
			outbin << values[0][j][i] << " " ;

		for( size_t j = 6 ; j < values[0].size() ; j++ )
		{
			unsigned char val = norm[j - 6][i] ;
			outbin << ( unsigned short int ) val << " " ;
		}

		outbin << std::endl ;
	}

	outbin.close();
}


void TriangleWriter::append()
{

	for( size_t k = 0 ; k < nTriangles.size() ; k++ )
	{
		writeHeader( true ) ;
		std::fstream outfile  ;
		outfile.open( filename.c_str(), std::ios::out | std::ios::app ) ;

		for( int i = 0 ; i < nTriangles[k] ; i++ )
		{
			for( size_t j = 0 ; j < values[k].size() ; j++ )
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
	writeHeader( true ) ;

	std::vector<std::valarray<unsigned char> > norm ;
	std::valarray<bool> voids( false, values[0].size() ) ;

	for( size_t i = 6 ; i < values.size() ; i++ )
	{
		std::valarray<unsigned char> n = normalizeArray( values[0][i], voids ) ;
		norm.push_back( n ) ;
	}

	std::fstream outbin ;
	outbin.open( filename.c_str(), std::ios::out | std::ios::app ) ;

	for( int i = 0 ; i < nTriangles[0] ; i++ )
	{
		for( size_t j = 0 ; j < 6 ; j++ )
			outbin << values[0][j][i] << " " ;

		for( size_t j = 6 ; j < values[0].size() ; j++ )
		{
			unsigned char val = norm[j - 6][i] ;
			outbin << ( unsigned short int ) val << " " ;
		}

		outbin << std::endl ;
	}

	outbin.close();
}

void MultiTriangleWriter::append()
{

	writeHeader( layers[0], true ,true ) ;
	std::fstream outfile  ;
	outfile.open( filename.c_str(), std::ios::out | std::ios::app ) ;

	for( int i = 0 ; i < nTriangles[0] ; i++ )
	{
		for( size_t j = 0 ; j < values[0].size() ; j++ )
		{
			outfile << values[0][j][i] << " " ;
		}

		outfile << std::endl ;
	}

	outfile.close();

	std::string filename_orig = filename;

	for( size_t k = 1 ; k < layers.size() ; k++ )
	{
		writeHeader( layers[k], false, true ) ;
		outfile.open( filename_orig.c_str(), std::ios::out | std::ios::app ) ;

		for( int i = 0 ; i < nTriangles[k] ; i++ )
		{
			for( size_t j = 0 ; j < values[k].size() ; j++ )
			{
				outfile << values[k][j][i] << " " ;
			}

			outfile << std::endl ;
		}

		outfile.close();
	}
}


void TriangleWriter::getField( TWFieldType field, bool extra )
{
	for( size_t j = 0 ; j < layers.size() ; j++ )
	{
		std::vector<std::valarray<double> > val = getDoubleValues( field, layers[j] ) ;
		std::reverse( val.begin(), val.end() );
		values[layerTranslator[layers[j]]].insert( values[layerTranslator[layers[j]]].end(), val.begin(), val.end() ) ;
	}

}

void TriangleWriter::getField( FieldType field, bool extra )
{
	for( size_t j = 0 ; j < layers.size() ; j++ )
	{
		std::vector<std::valarray<double> > val = getDoubleValues( field, layers[j] ) ;
		std::reverse( val.begin(), val.end() );
		values[layerTranslator[layers[j]]].insert( values[layerTranslator[layers[j]]].end(), val.begin(), val.end() ) ;
	}

}

std::vector<std::valarray<double> > TriangleWriter::getDoubleValues( TWFieldType field, int layer )
{
	std::vector<std::valarray<double> > ret ;
	int iterator = 0 ;

	for( int i = 0 ; i < numberOfFields( field ) ; i++ )
	{
		std::valarray<double> reti( nTriangles[layerTranslator[layer]] ) ;
		ret.push_back( reti ) ;
	}
	
	if( field == TWFT_STRAIN || field == TWFT_STRAIN_AND_STRESS || field == TWFT_STRESS )
	{

		std::pair<Vector, Vector> stress_strain ;//= source->getStressAndStrainInLayer( layer, false ) ;
		std::vector<DelaunayTriangle *> triangles = source->getElements2DInLayer( layer ) ;
		int pointsPerTri = triangles[0]->getBoundingPoints().size() ;
		int pointsPerPlane = pointsPerTri / triangles[0]->timePlanes() ;
		int factor = 1 ;

		if( pointsPerPlane % 6 == 0 )
			factor = 2 ;

		int time_offset = timePlane[layerTranslator[layer]] * pointsPerTri / triangles[0]->timePlanes() ;

		switch( field )
		{
			case TWFT_STRAIN:
				stress_strain.first.resize( 0 ) ;

				for( int i = 0 ; i < triangles.size() ; i++ )
				{
					if( triangles[i]->getBehaviour() && triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR && !triangles[i]->getBehaviour()->fractured() )
					{
						Vector strain0(3) ;
						Vector strain1(3) ;
						Vector strain2(3) ;
						
						triangles[i]->getState().getField(STRAIN_FIELD,  triangles[i]->getBoundingPoint(time_offset + 0), strain0, false);
						triangles[i]->getState().getField(STRAIN_FIELD,  triangles[i]->getBoundingPoint(time_offset + factor), strain1, false);
						triangles[i]->getState().getField(STRAIN_FIELD,  triangles[i]->getBoundingPoint(time_offset + factor*2), strain2, false);
						
						// epsilon11
						ret[8][iterator]   = strain0[0];
						ret[7][iterator]   = strain1[0];
						ret[6][iterator]   = strain2[0];

						// epsilon12
						ret[5][iterator]   = strain0[1];
						ret[4][iterator]   = strain1[1];
						ret[3][iterator]   = strain2[1];

						// epsilon22
						ret[2][iterator]   = strain0[2];
						ret[1][iterator]   = strain1[2];
						ret[0][iterator++]   = strain2[2];
					  
					}

				}

				break ;

			case TWFT_STRAIN_AND_STRESS:

				for( int i = 0 ; i < triangles.size() ; i++ )
				{
					if( triangles[i]->getBehaviour() && triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR && !triangles[i]->getBehaviour()->fractured() )
					{
						Vector strain0(3) ;
						Vector strain1(3) ;
						Vector strain2(3) ;
						triangles[i]->getState().getField(STRAIN_FIELD, triangles[i]->getBoundingPoint(time_offset + 0), strain0, false);
						triangles[i]->getState().getField(STRAIN_FIELD, triangles[i]->getBoundingPoint(time_offset + factor), strain1, false);
						triangles[i]->getState().getField(STRAIN_FIELD, triangles[i]->getBoundingPoint(time_offset + 2*factor), strain2, false);
						
						// epsilon11
						ret[17][iterator]   = strain0[0];
						ret[16][iterator]   = strain1[0];
						ret[15][iterator]   = strain2[0];

						// epsilon12
						ret[14][iterator]   = strain0[1];
						ret[13][iterator]   = strain1[1];
						ret[12][iterator]   = strain2[1];

						// epsilon22
						ret[11][iterator]   = strain0[2];
						ret[10][iterator]   = strain1[2];
						ret[9][iterator]    = strain2[2];
						// epsilon11

						triangles[i]->getState().getField(REAL_STRESS_FIELD, triangles[i]->getBoundingPoint(time_offset + 0), strain0, false);
						triangles[i]->getState().getField(REAL_STRESS_FIELD, triangles[i]->getBoundingPoint(time_offset + factor), strain1, false);
						triangles[i]->getState().getField(REAL_STRESS_FIELD, triangles[i]->getBoundingPoint(time_offset + 2*factor), strain2, false);
						// sigma11
						ret[8][iterator]   = strain0[0];
						ret[7][iterator]   = strain1[0];
						ret[6][iterator]   = strain2[0];

						// sigma12
						ret[5][iterator]   = strain0[1];
						ret[4][iterator]   = strain1[1];
						ret[3][iterator]   = strain2[1];

						// sigma22
						ret[2][iterator]   = strain0[2];
						ret[1][iterator]   = strain1[2];
						ret[0][iterator++] = strain2[2];
					}
				}
				break ;

			case TWFT_STRESS:
				stress_strain.second.resize( 0 ) ;

				for( int i = 0 ; i < triangles.size() ; i++ )
				{
					if( triangles[i]->getBehaviour()&& triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR && !triangles[i]->getBehaviour()->fractured() )
					{
						Vector strain0(3) ;
						Vector strain1(3) ;
						Vector strain2(3) ;
						triangles[i]->getState().getField(REAL_STRESS_FIELD, triangles[i]->getBoundingPoint(time_offset + 0), strain0, false);
						triangles[i]->getState().getField(REAL_STRESS_FIELD, triangles[i]->getBoundingPoint(time_offset + factor), strain1, false);
						triangles[i]->getState().getField(REAL_STRESS_FIELD, triangles[i]->getBoundingPoint(time_offset + 2*factor), strain2, false);

						// epsilon11
						ret[8][iterator]   = strain0[0];
						ret[7][iterator]   = strain1[0];
						ret[6][iterator]   = strain2[0];

						// epsilon12
						ret[5][iterator]   = strain0[1];
						ret[4][iterator]   = strain1[1];
						ret[3][iterator]   = strain2[1];

						// epsilon22
						ret[2][iterator]   = strain0[2];
						ret[1][iterator]   = strain1[2];
						ret[0][iterator++] = strain2[2];
					}
				}

				break ;
		}
	}
	else
	{
		if( field == TWFT_GRADIENT || field == TWFT_GRADIENT_AND_FLUX || field == TWFT_FLUX )
		{
			std::pair<Vector, Vector> gradient_flux = source->getGradientAndFluxInLayer( layer, false ) ;
			std::vector<DelaunayTriangle *> triangles = source->getElements2DInLayer( layer ) ;

			switch( field )
			{

				case TWFT_FLUX:
					gradient_flux.first.resize( 0 ) ;

					for( int i = 0 ; i < triangles.size() ; i++ )
					{
						if(  triangles[i]->getBehaviour() && triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR && !triangles[i]->getBehaviour()->fractured() )
						{
							// j11
							ret[5][iterator] = gradient_flux.second[i * 3 * 2 + 0] ;
							ret[4][iterator] = gradient_flux.second[i * 3 * 2 + 2] ;
							ret[3][iterator] = gradient_flux.second[i * 3 * 2 + 4] ;

							// j22
							ret[2][iterator] = gradient_flux.second[i * 3 * 2 + 1] ;
							ret[1][iterator] = gradient_flux.second[i * 3 * 2 + 3] ;
							ret[0][iterator++] = gradient_flux.second[i * 3 * 2 + 5] ;
						}
					}

					break ;

				case TWFT_GRADIENT_AND_FLUX:

					for( int i = 0 ; i < triangles.size() ; i++ )
					{
						if(  triangles[i]->getBehaviour() && triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR && !triangles[i]->getBehaviour()->fractured() )
						{
							// d11
							ret[11][iterator] = gradient_flux.first[i * 3 * 2 + 0] ;
							ret[10][iterator] = gradient_flux.first[i * 3 * 2 + 2] ;
							ret[9][iterator] = gradient_flux.first[i * 3 * 2 + 4] ;

							// d22
							ret[8][iterator] = gradient_flux.first[i * 3 * 2 + 1] ;
							ret[7][iterator] = gradient_flux.first[i * 3 * 2 + 3] ;
							ret[6][iterator] = gradient_flux.first[i * 3 * 2 + 5] ;

							// j11
							ret[5][iterator] = gradient_flux.second[i * 3 * 2 + 0] ;
							ret[4][iterator] = gradient_flux.second[i * 3 * 2 + 2] ;
							ret[3][iterator] = gradient_flux.second[i * 3 * 2 + 4] ;

							// j22
							ret[2][iterator] = gradient_flux.second[i * 3 * 2 + 1] ;
							ret[1][iterator] = gradient_flux.second[i * 3 * 2 + 3] ;
							ret[0][iterator++] = gradient_flux.second[i * 3 * 2 + 5] ;
						}
					}

					break ;

				case TWFT_GRADIENT:
					gradient_flux.second.resize( 0 ) ;

					for( int i = 0 ; i < triangles.size() ; i++ )
					{
						if(  triangles[i]->getBehaviour() &&triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR && !triangles[i]->getBehaviour()->fractured() )
						{
							// d11
							ret[5][iterator] = gradient_flux.first[i * 3 * 2 + 0] ;
							ret[4][iterator] = gradient_flux.first[i * 3 * 2 + 2] ;
							ret[3][iterator] = gradient_flux.first[i * 3 * 2 + 4] ;

							// d22
							ret[2][iterator] = gradient_flux.first[i * 3 * 2 + 1] ;
							ret[1][iterator] = gradient_flux.first[i * 3 * 2 + 3] ;
							ret[0][iterator++] = gradient_flux.first[i * 3 * 2 + 5] ;
						}
					}

					break ;
			}

		}
		else
		{
			if( field == TWFT_DISPLACEMENTS )
			{
				Vector x = source->getDisplacements(-1, false) ;
				std::vector<DelaunayTriangle *> triangles = source->getElements2DInLayer( layer ) ;
				int pointsPerTri = triangles[0]->getBoundingPoints().size() ;
				int pointsPerTimePlanes = pointsPerTri / triangles[0]->timePlanes() ;
				int factor = pointsPerTimePlanes / 3 ;
				if( timePlane[layerTranslator[layer]] >= triangles[0]->timePlanes() )
					timePlane[layerTranslator[layer]] = triangles[0]->timePlanes() - 1 ;

				int time_offset = timePlane[layerTranslator[layer]] * pointsPerTri / triangles[0]->timePlanes() ;

				for( int i = 0 ; i < triangles.size() ; i++ )
				{
					if(  triangles[i]->getBehaviour() && triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR )
					{
						size_t dof = triangles[i]->getBehaviour()->getNumberOfDegreesOfFreedom() ;
					  
						size_t id1 = triangles[i]->getBoundingPoint( factor * 0 + time_offset ).id ;
						size_t id2 = triangles[i]->getBoundingPoint( factor * 1 + time_offset ).id ;
						size_t id3 = triangles[i]->getBoundingPoint( factor * 2 + time_offset ).id ;

						ret[5][iterator] = x[id1 * dof] ;
						ret[4][iterator] = x[id2 * dof] ;
						ret[3][iterator] = x[id3 * dof] ;
						ret[2][iterator] = x[id1 * dof + 1] ;
						ret[1][iterator] = x[id2 * dof + 1] ;
						ret[0][iterator++] = x[id3 * dof + 1] ;
					}
				}
			}
			else if( field == TWFT_CRACKS )
			{
				std::vector<DelaunayTriangle *> triangles = source->getElements2DInLayer( layer ) ;

				for( int i = 0 ; i < triangles.size() ; i++ )
				{
					if( triangles[i]->getBehaviour() && triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR &&  triangles[i]->getBehaviour()->getDamageModel() &&  triangles[i]->getBehaviour()->getDamageModel()->getState().max() > POINT_TOLERANCE_2D)
					{

						std::pair<double, double> np = triangles[i]->getBehaviour()->getFractureCriterion()->getCrackOpeningAndSlip(triangles[i]->getState()) ;
						ret[0][iterator] = np.first;
						ret[1][iterator] = np.first ;
						ret[2][iterator] = np.first ;
						ret[3][iterator] = np.second;
						ret[4][iterator] = np.second ;
						ret[5][iterator++] = np.second ;

					}
					else if ( triangles[i]->getBehaviour() && triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR )
					{
						ret[0][iterator] = 0;
						ret[1][iterator] = 0 ;
						ret[2][iterator] = 0 ;
						ret[3][iterator] = 0;
						ret[4][iterator] = 0 ;
						ret[5][iterator++] = 0 ;
					}
				}
			}
			else if( field == TWFT_DAMAGE)
			{
				std::vector<DelaunayTriangle *> triangles = source->getElements2DInLayer( layer ) ;

				for( int i = 0 ; i < triangles.size() ; i++ )
				{
					if( triangles[i]->getBehaviour() && triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR &&  triangles[i]->getBehaviour()->getDamageModel() )
					{
						double d = triangles[i]->getBehaviour()->getDamageModel()->getState().max();

						if( triangles[i]->getBehaviour()->getDamageModel()->fractured() )
							d = 1 ;

						ret[0][iterator] = d;
						ret[1][iterator] = d ;
						ret[2][iterator++] = d ;

					}
					else if ( triangles[i]->getBehaviour() && triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR )
					{
						ret[0][iterator] = 0;
						ret[1][iterator] = 0 ;
						ret[2][iterator++] = 0 ;
					}
				}
			}
			else if( field == TWFT_IMPOSED_STRESS_NORM )
			{
				std::vector<DelaunayTriangle *> triangles = source->getElements2DInLayer( layer ) ;
				Vector x( triangles.size() * 3 ) ;

				for( int i = 0 ; i < triangles.size() ; i++ )
				{
					if( triangles[i]->getBehaviour() && triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR &&  triangles[i]->getBehaviour()->hasInducedForces() && triangles[i]->getBehaviour()->getDamageModel())
					{
						Vector dv = triangles[i]->getBehaviour()->getDamageModel()->getImposedStress(triangles[i]->getCenter());

						double d = sqrt(std::inner_product(&dv[0], &dv[dv.size()], &dv[0], 0.));
						if( triangles[i]->getBehaviour()->getDamageModel()->fractured() )
							d = 1 ;

						ret[0][iterator] = d;
						ret[1][iterator] = d ;
						ret[2][iterator++] = d ;

					}
					else if ( triangles[i]->getBehaviour() && triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR )
					{
						ret[0][iterator] = 0;
						ret[1][iterator] = 0 ;
						ret[2][iterator++] = 0 ;
					}
				}
			}
			else
			{
				std::vector<DelaunayTriangle *> tri = source->getElements2DInLayer( layer ) ;

				for( size_t i = 0 ; i < tri.size() ; i++ )
				{
					std::pair<bool, std::vector<double> > val = getDoubleValue( tri[i], field ) ;

					if(  tri[i]->getBehaviour() && tri[i]->getBehaviour()->type != VOID_BEHAVIOUR )
					{
						if( val.first )
						{
							for( size_t j = 0 ; j < numberOfFields( field ) ; j++ )
								ret[j][iterator] = val.second[j] ;

							iterator++ ;
						}
					}
				}
			}
		}
	}

	return ret ;
}

std::vector<std::valarray<double> > TriangleWriter::getDoubleValues( FieldType field, int layer )
{
	std::vector<std::valarray<double> > ret ;
	int iterator = 0 ;
	
	std::vector<DelaunayTriangle *> triangles = source->getElements2DInLayer( layer ) ;
	int blocks = triangles[0]->getBehaviour()->getNumberOfDegreesOfFreedom() / 2 ;
	int size = fieldTypeElementarySize( field, SPACE_TWO_DIMENSIONAL, blocks) ;
	int pointsPerTri = triangles[0]->getBoundingPoints().size() ;
	int pointsPerPlane = pointsPerTri / triangles[0]->timePlanes() ;
	int factor = 1 ;

	if( pointsPerPlane % 6 == 0 )
		factor = 2 ;

	int time_offset = timePlane[layerTranslator[layer]] * pointsPerTri / triangles[0]->timePlanes() ;

	for( int i = 0 ; i < size*3 ; i++ )
	{
		std::valarray<double> reti( nTriangles[layerTranslator[layer]] ) ;
		ret.push_back( reti ) ;
	}
	

	for( size_t i = 0 ; i < triangles.size() ; i++ )
	{
		if(  triangles[i]->getBehaviour() && triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR )
		{
			Vector first(size) ;
			Vector second(size) ;
			Vector third(size) ;
		  
			triangles[i]->getState().getField(field,  triangles[i]->getBoundingPoint(time_offset + 0), first, false);
			triangles[i]->getState().getField(field,  triangles[i]->getBoundingPoint(time_offset + factor), second, false);
			triangles[i]->getState().getField(field,  triangles[i]->getBoundingPoint(time_offset + factor*2), third, false);
			
			for(size_t j = 0 ; j < size ; j++)
			{
				ret[j*3+0][iterator] = first[j] ;
				ret[j*3+1][iterator] = second[j] ;
				ret[j*3+2][iterator] = third[j] ;
			}
			
			iterator++ ;
			
		}
	}

	return ret ;
}

std::pair<bool, std::vector<double> > TriangleWriter::getDoubleValue( DelaunayTriangle *tri, TWFieldType field )
{
	bool found = false ;
	std::vector<double> ret( numberOfFields( field ) ) ;

	if( tri->getBehaviour() && tri->getBehaviour()->type != VOID_BEHAVIOUR )
	{
		switch( field )
		{
			case TWFT_COORDINATE:
				ret[5] = tri->getBoundingPoint(0).x ;
				ret[4] = tri->getBoundingPoint(0).y ;
				ret[3] = tri->getBoundingPoint(1).x ;
				ret[2] = tri->getBoundingPoint(1).y ;
				ret[1] = tri->getBoundingPoint(2).x ;
				ret[0] = tri->getBoundingPoint(2).y ;
				found = true ;
				break ;

			case TWFT_PRINCIPAL_ANGLE:
			{
				Vector v(0., 1) ;
				tri->getState().getField( PRINCIPAL_ANGLE_FIELD, *tri->first, v, false ) ;
				ret[2] = 180.*v[0]/M_PI ;
				tri->getState().getField( PRINCIPAL_ANGLE_FIELD, *tri->second, v, false ) ;
				ret[1] = 180.*v[0]/M_PI ;
				tri->getState().getField( PRINCIPAL_ANGLE_FIELD, *tri->third, v, false ) ;
				ret[0] = 180.*v[0]/M_PI ;

				found = true ;
				break ;
			}
			case TWFT_TRIANGLE_ANGLE:
			{
				Vector v(0., 1) ;
				double a = tri->area() ;
				double r = tri->getRadius() ;
				ret[2] = a/(M_PI*r*r) ;
				ret[1] = a/(M_PI*r*r) ;
				ret[0] = a/(M_PI*r*r) ;
// 				tri->getState().getField( PRINCIPAL_ANGLE_FIELD, *tri->first, v, false ) ;
// 				ret[2] = 180.*v[0]/M_PI ;
// 				tri->getState().getField( PRINCIPAL_ANGLE_FIELD, *tri->second, v, false ) ;
// 				ret[1] = 180.*v[0]/M_PI ;
// 				tri->getState().getField( PRINCIPAL_ANGLE_FIELD, *tri->third, v, false ) ;
// 				ret[0] = 180.*v[0]/M_PI ;

				found = true ;
				break ;
			}
			case TWFT_CRACKS:
			{
				if( tri->getBehaviour()->getDamageModel()->getState().max() > POINT_TOLERANCE_2D)
				{
					std::pair<double, double> np = tri->getBehaviour()->getFractureCriterion()->getCrackOpeningAndSlip(tri->getState()) ;
					ret[0] = np.first;
					ret[1] = np.first ;
					ret[2] = np.first ;
					ret[3] = np.second;
					ret[4] = np.second ;
					ret[5] = np.second ;

				}
				else
				{
					ret[0] = 0;
					ret[1] = 0 ;
					ret[2] = 0 ;
					ret[3] = 0;
					ret[4] = 0 ;
					ret[5] = 0 ;
				}
				found = true ;
				break ;
			}
			
			case TWFT_CRACK_ANGLE:
			{
				if( tri->getBehaviour()->getFractureCriterion() && tri->getBehaviour()->getDamageModel()->getState().max() > POINT_TOLERANCE_2D)
				{
					double angle = 180.*tri->getBehaviour()->getFractureCriterion()->smoothedCrackAngle(tri->getState())/M_PI ;
					while (angle < 0)
						angle += 90 ;
					ret[2] =  angle ;
					ret[1] =  angle ;
					ret[0] =  angle ;
				}
				else
				{
					ret[2] = 0 ;
					ret[1] = 0 ;
					ret[0] = 0 ;
				}

				found = true ;
				break ;
			}
			
			case TWFT_CRITERION:
			{
				if( tri->getBehaviour() && tri->getBehaviour()->getFractureCriterion())
				{
					double d = tri->getBehaviour()->getFractureCriterion()->getScoreAtState() ;
					ret[2] = d ;
					ret[1] = d ;
					ret[0] = d ;
				}
				else if( tri->getBehaviour() )
				{
					double d = -1 ;
					ret[2] = d ;
					ret[1] = d ;
					ret[0] = d ;
				}

				found = true ;
				break ;
			}

			case TWFT_PRINCIPAL_STRESS:
			{
				Vector v(0., 2) ;
				tri->getState().getField( PRINCIPAL_REAL_STRESS_FIELD, *tri->first, v, false) ;
				ret[5] = v[1] ;
				ret[2] = v[0] ;
				tri->getState().getField( PRINCIPAL_REAL_STRESS_FIELD, *tri->second, v, false) ;
				ret[4] = v[1] ;
				ret[1] = v[0] ;
				tri->getState().getField( PRINCIPAL_REAL_STRESS_FIELD, *tri->third, v, false) ;
				ret[3] = v[1] ;
				ret[0] = v[0] ;
				found = true ;
				break ;
			}
			case TWFT_PRINCIPAL_STRAIN:
			{
				Vector v(0., 2) ;
				tri->getState().getField( PRINCIPAL_STRAIN_FIELD, *tri->first, v, false) ;
				ret[5] = v[1] ;
				ret[2] = v[0] ;
				tri->getState().getField( PRINCIPAL_STRAIN_FIELD, *tri->second, v, false) ;
				ret[4] = v[1] ;
				ret[1] = v[0] ;
				tri->getState().getField( PRINCIPAL_STRAIN_FIELD, *tri->third, v, false) ;
				ret[3] = v[1] ;
				ret[0] = v[0] ;
				found = true ;
				break ;
			}
			case TWFT_STIFFNESS:
			{
				double t = 0 ;

				if( tri->timePlanes() > 1 )
					t = -1 + timePlane[0] * 2 / ( (int)tri->timePlanes() - 1 ) ;

				if( dynamic_cast<LinearForm *>( tri->getBehaviour() ) )
				{
					ret[2] = (dynamic_cast<LinearForm *>( tri->getBehaviour() )->getTensor( Point( 0.5, 0, 0, t ) )[0][0]+dynamic_cast<LinearForm *>( tri->getBehaviour() )->getTensor( Point( 0.5, 0, 0, t ) )[1][1])*.5 ;
					ret[1] = ret[2] ;
					ret[0] = ret[2] ;
					found = true ;
				}

				break ;
			}
			case TWFT_STIFFNESS_X:
			{
				LinearForm *b = dynamic_cast<LinearForm *>( tri->getBehaviour() ) ;
				double t = 0 ;

				if( tri->timePlanes() > 1 )
					t = -1 + timePlane[0] * 2 / ( (int)tri->timePlanes() - 1 ) ;

				if( dynamic_cast<LinearForm *>( tri->getBehaviour() ) )
				{
					ret[2] =dynamic_cast<LinearForm *>( tri->getBehaviour() )->getTensor( Point( 1, 0, 0, t ) )[0][0] ;
					ret[1] = ret[2] ;
					ret[0] = ret[2] ;
					found = true ;
				}

				break ;
			}
			case TWFT_STIFFNESS_Y:
			{
				LinearForm *b = dynamic_cast<LinearForm *>( tri->getBehaviour() ) ;
				double t = 0 ;

				if( tri->timePlanes() > 1 )
					t = -1 + timePlane[0] * 2 / ( (int)tri->timePlanes() - 1 ) ;

				if( dynamic_cast<LinearForm *>( tri->getBehaviour() ) )
				{
					ret[2] = dynamic_cast<LinearForm *>( tri->getBehaviour() )->getTensor( Point( 1, 0, 0, t ) )[1][1] ;
					ret[1] = ret[2] ;
					ret[0] = ret[2] ;
					found = true ;
				}

				break ;
			}
			case TWFT_STIFFNESS_Z:
			{
				LinearForm *b = dynamic_cast<LinearForm *>( tri->getBehaviour() ) ;
				double t = 0 ;

				if( tri->timePlanes() > 1 )
					t = -1 + timePlane[0] * 2 / ( (int)tri->timePlanes() - 1 ) ;

				if( dynamic_cast<LinearForm *>( tri->getBehaviour() ) )
				{
					ret[2] = dynamic_cast<LinearForm *>( tri->getBehaviour() )->getTensor( Point( 1, 0, 0, t ) )[2][2] ;
					ret[1] = ret[2] ;
					ret[0] = ret[2] ;
					found = true ;
				}

				break ;
			}

			case TWFT_VON_MISES:
			{
				Vector v(0., 1) ;
				tri->getState().getField( VON_MISES_REAL_STRESS_FIELD, *tri->first, v, false ) ;
				ret[2] = ret[1] = ret[0] = v[0] ;
				found = true ;
				break ;
			}

		}
	}

	return std::make_pair( found, ret ) ;

}

void TriangleWriter::writeHeader( int layer, bool append )
{
	std::fstream outstream ;

	if( append )
		outstream.open( filename.c_str(), std::ios::out | std::ios::app ) ;
	else
		outstream.open( filename.c_str(), std::ios::out ) ;

	outstream << "TRIANGLES" << std::endl ;
	outstream << ( int ) nTriangles[layerTranslator[layer]] << std::endl ;
	outstream << 3 << std::endl ;
	outstream << ( ( int ) values[layerTranslator[layer]].size() - 6 ) / 3 << std::endl ;
	outstream.close() ;
}

void BinaryTriangleWriter::writeHeader(int layer, bool append )
{
	std::fstream outstream ;

	if( append )
		outstream.open( filename.c_str(), std::ios::out | std::ios::app ) ;
	else
		outstream.open( filename.c_str(), std::ios::out ) ;

	outstream << "BIN_TRIANGLES" << std::endl ;
	outstream << ( int ) values[0].size() << std::endl ;
	outstream << 3 << std::endl ;
	outstream << ( ( int ) values.size() - 6 ) / 3 << std::endl ;
	outstream.close() ;
}

void MultiTriangleWriter::writeHeader( int layer, bool firstlayer, bool append )
{
	if(firstlayer)
	{
		nextCounter() ;
		std::fstream headstream ;

		headstream.open( head.c_str(), std::ios::out | std::ios::app ) ;
		headstream << filename << std::endl ;
		headstream.close();
	}

	std::fstream outstream ;
	if( append )
		outstream.open( filename.c_str(), std::ios::out | std::ios::app ) ;
	else
		outstream.open( filename.c_str(), std::ios::out ) ;
	outstream << "TRIANGLES" << std::endl ;
	outstream << ( int ) nTriangles[layerTranslator[layer]] << std::endl ;
	outstream << 3 << std::endl ;
	outstream << ( ( int ) values[layerTranslator[layer]].size() - 6 ) / 3 << std::endl ;
	outstream.close() ;
}



int numberOfFields( TWFieldType field )
{
	switch( field )
	{
		case TWFT_COORDINATE:
			return 6 ;
		case TWFT_DISPLACEMENTS:
			return 6 ;
		case TWFT_PRINCIPAL_ANGLE:
			return 3 ;
		case TWFT_TRIANGLE_ANGLE:
			return 3 ;
		case TWFT_CRACK_ANGLE:
			return 3 ;
		case TWFT_STIFFNESS:
			return 3 ;
		case TWFT_STIFFNESS_X:
			return 3 ;
		case TWFT_STIFFNESS_Y:
			return 3 ;
		case TWFT_STIFFNESS_Z:
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
		case TWFT_CRITERION:
			return 3 ;
		case TWFT_DAMAGE:
			return 3 ;
		case TWFT_IMPOSED_STRESS_NORM:
			return 3 ;
		case TWFT_CRACKS:
			return 6 ;
	}

	return 3 ;
}


}



























