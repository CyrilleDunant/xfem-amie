#include "composite.h"
#include "../../utilities/matrixops.h"

using namespace Mu ;

Composite::Composite( DelaunayTriangle *tri, std::vector<Feature *> feats ) : Phase( tri )
{

}

Composite::Composite( DelaunayTetrahedron *tet, std::vector<Feature *> feats ) : Phase( tet )
{

}

Composite::Composite( Phase p ) : Phase( p )
{

}

void Composite::apply()
{

}

Matrix Composite::I4( Matrix C )
{
	Matrix I4( C.numRows(), C.numCols() ) ;

	for( size_t i = 0 ; i < 3 + 3 * ( I4.size() == 36 ) ; i++ )
	{
		I4[i][i] = 1 ;
	}

	return I4 ;
}

Matrix Composite::eshelby( Matrix C )
{
	Matrix S( C.numRows(), C.numCols() ) ;
	double nu = 0.2 ;
	double a = C[0][0] ;
	double b = C[0][1] ;
//    if(C.size()==36)
	nu = b / ( a + b ) ;
//    else
//	nu = b / a ;

	if( S.size() == 9 )
	{
		double Siiii = ( 7 - 5 * nu ) / ( 15 * ( 1 - nu ) ) ;
		double Siijj = ( 5 * nu - 1 ) / ( 15 * ( 1 - nu ) ) ;
		double Sijij = ( 4 - 5 * nu ) / ( 15 * ( 1 - nu ) ) ;

		S[0][0] = Siiii ;
		S[1][1] = Siiii ;

		S[0][1] = Siijj ;
		S[1][0] = Siijj ;

		S[2][2] = Sijij ;
	}

	if( S.size() == 36 )
	{
		double Siiii = ( 7 - 5 * nu ) / ( 15 * ( 1 - nu ) ) ;
		double Siijj = ( 5 * nu - 1 ) / ( 15 * ( 1 - nu ) ) ;
		double Sijij = ( 4 - 5 * nu ) / ( 15 * ( 1 - nu ) ) ;

		S[0][0] = Siiii ;
		S[1][1] = Siiii ;
		S[2][2] = Siiii ;

		S[0][1] = Siijj ;
		S[0][2] = Siijj ;
		S[1][0] = Siijj ;
		S[1][2] = Siijj ;
		S[2][0] = Siijj ;
		S[2][1] = Siijj ;

		S[3][3] = Sijij ;
		S[4][4] = Sijij ;
		S[5][5] = Sijij ;
	}

	return S ;
}

void Composite::invertTensor( Matrix &m )
{
	Matrix m1( 3, 3 ) ;

	for( size_t i = 0 ; i < 3 ; i++ )
	{
		for( size_t j = 0 ; j < 3 ; j++ )
			m1[i][j] = m[i][j] ;
	}

	invert3x3Matrix( m1 ) ;

	for( size_t i = 0 ; i < 3 ; i++ )
	{
		for( size_t j = 0 ; j < 3 ; j++ )
			m[i][j] = m1[i][j] ;
	}

	for( size_t i = 3 ; i < m.numCols() ; i++ )
		m[i][i] = 1. / m[i][i] ;
}

MatrixInclusionComposite::MatrixInclusionComposite( DelaunayTriangle *tri, Feature *inc ) : Composite( tri, std::vector<Feature *>( 0 ) )
{
	matrix = Phase( tri ) ;
	inclusion = Phase( inc ) ;

	matrix.volume -= inclusion.volume ;
	volume = matrix.volume + inclusion.volume ;
	matrix.volume /= volume ;
	inclusion.volume /= volume ;
}

MatrixInclusionComposite::MatrixInclusionComposite( DelaunayTetrahedron *tet, Feature *inc ) : Composite( tet, std::vector<Feature *>( 0 ) )
{
	matrix = Phase( tet ) ;
	inclusion = Phase( inc ) ;

	matrix.volume -= inclusion.volume ;
	volume = matrix.volume + inclusion.volume ;
	matrix.volume /= volume ;
	inclusion.volume /= volume ;
}

MatrixInclusionComposite::MatrixInclusionComposite( Phase mat, Phase inc ) : Composite( mat )
{
	matrix = Phase( mat ) ;
	inclusion = Phase( inc ) ;

	volume = inclusion.volume ;

	inclusion.volume = 1 - matrix.volume ;

}

void MatrixInclusionComposite::apply()
{
	this->getStrainConcentrationTensor() ;
	Matrix I = Composite::I4( matrix.C ) ;

	inclusion.A = B * inclusion.volume + I * matrix.volume ;
	Composite::invertTensor( inclusion.A ) ;
	inclusion.A *= B ;

	matrix.A = ( I - inclusion.A * inclusion.volume ) / matrix.volume ;

	C = ( ( matrix.A * matrix.volume ) * matrix.C ) ;
	C += ( ( inclusion.A * inclusion.volume ) * inclusion.C ) ;
	beta.resize( matrix.beta.size() );
	beta = matrix.A * matrix.volume * matrix.beta ;
	beta += inclusion.A * inclusion.volume * inclusion.beta ;

	/*    lambda.clear() ;
	    Vector b = inclusion.beta ;
	    b -= matrix.beta ;
	    Matrix Cdiff = inclusion.C - matrix.C ;
	    if(Cdiff.size()==36)
		invert6x6Matrix(Cdiff) ;
	    else
		invert3x3Matrix(Cdiff) ;
	    Matrix G = inclusion.A - I ;
	    G *= Cdiff ;
	    Vector a = G * b ;

	    for(size_t i = 0 ; i < matrix.lambda.size() ; i++)
	    {
		Matrix S = matrix.C ;
		if(S.size()==36)
		    invert6x6Matrix(S) ;
		else
		    invert3x3Matrix(S) ;
		Vector s1 = matrix.lambda[i] - matrix.beta ;
		s1 = S * s1 ;
		s1 = s1 - a ;

		Matrix K = inclusion.A ;
		if(K.size()==36)
		    invert6x6Matrix(K) ;
		else
		    invert3x3Matrix(K) ;
		K *= C ;
		s1 = K * s1 ;
		s1 = s1 + beta ;
	    }*/

}

void MatrixInclusionComposite::getStrainConcentrationTensor()
{
	B = Matrix( C.numRows(), C.numCols() ) ;
}

DiluteMatrixInclusionComposite::DiluteMatrixInclusionComposite( DelaunayTriangle *tri, Feature *inc ) : MatrixInclusionComposite( tri, inc )
{

}

DiluteMatrixInclusionComposite::DiluteMatrixInclusionComposite( DelaunayTetrahedron *tet, Feature *inc ) : MatrixInclusionComposite( tet, inc )
{

}

DiluteMatrixInclusionComposite::DiluteMatrixInclusionComposite( Phase mat, Phase inc ) : MatrixInclusionComposite( mat, inc )
{

}

void DiluteMatrixInclusionComposite::getStrainConcentrationTensor()
{
	MoriTanakaMatrixInclusionComposite mori( matrix, inclusion ) ;
	mori.apply() ;

	Matrix I = Composite::I4( mori.B ) ;
	B = I - ( mori.B * inclusion.volume ) ;
	Composite::invertTensor( B ) ;
	B *= mori.B * matrix.volume ;
}


VoigtMatrixInclusionComposite::VoigtMatrixInclusionComposite( DelaunayTriangle *tri, Feature *inc ) : MatrixInclusionComposite( tri, inc )
{

}

VoigtMatrixInclusionComposite::VoigtMatrixInclusionComposite( DelaunayTetrahedron *tet, Feature *inc ) : MatrixInclusionComposite( tet, inc )
{

}

VoigtMatrixInclusionComposite::VoigtMatrixInclusionComposite( Phase mat, Phase inc ) : MatrixInclusionComposite( mat, inc )
{

}

void VoigtMatrixInclusionComposite::getStrainConcentrationTensor()
{
	B = Composite::I4( C ) ;
}

ReussMatrixInclusionComposite::ReussMatrixInclusionComposite( DelaunayTriangle *tri, Feature *inc ) : MatrixInclusionComposite( tri, inc )
{

}

ReussMatrixInclusionComposite::ReussMatrixInclusionComposite( DelaunayTetrahedron *tet, Feature *inc ) : MatrixInclusionComposite( tet, inc )
{

}

ReussMatrixInclusionComposite::ReussMatrixInclusionComposite( Phase mat, Phase inc ) : MatrixInclusionComposite( mat, inc )
{

}

void ReussMatrixInclusionComposite::getStrainConcentrationTensor()
{
	B = inclusion.C ;
	Composite::invertTensor( B );
	B *= matrix.C ;
}

MoriTanakaMatrixInclusionComposite::MoriTanakaMatrixInclusionComposite( DelaunayTriangle *tri, Feature *inc ) : MatrixInclusionComposite( tri, inc )
{

}

MoriTanakaMatrixInclusionComposite::MoriTanakaMatrixInclusionComposite( DelaunayTetrahedron *tet, Feature *inc ) : MatrixInclusionComposite( tet, inc )
{

}

MoriTanakaMatrixInclusionComposite::MoriTanakaMatrixInclusionComposite( Phase mat, Phase inc ) : MatrixInclusionComposite( mat, inc )
{

}

void MoriTanakaMatrixInclusionComposite::getStrainConcentrationTensor()
{
	Matrix I = Composite::I4( C ) ;
	Matrix S = Composite::eshelby( matrix.C ) ;

	B = matrix.C ;
	Composite::invertTensor( B ) ;
	B *= S * ( inclusion.C - matrix.C ) ;
	B += I ;
	Composite::invertTensor( B ) ;

}

InverseMoriTanakaMatrixInclusionComposite::InverseMoriTanakaMatrixInclusionComposite( DelaunayTriangle *tri, Feature *inc ) : MatrixInclusionComposite( tri, inc )
{

}

InverseMoriTanakaMatrixInclusionComposite::InverseMoriTanakaMatrixInclusionComposite( DelaunayTetrahedron *tet, Feature *inc ) : MatrixInclusionComposite( tet, inc )
{

}

InverseMoriTanakaMatrixInclusionComposite::InverseMoriTanakaMatrixInclusionComposite( Phase mat, Phase inc ) : MatrixInclusionComposite( mat, inc )
{

}

void InverseMoriTanakaMatrixInclusionComposite::getStrainConcentrationTensor()
{
	Matrix I = Composite::I4( C ) ;
	Matrix S = Composite::eshelby( inclusion.C ) ;

	B = inclusion.C ;
	Composite::invertTensor( B ) ;
	B *= S * ( matrix.C - inclusion.C ) ;
	B += I ;
//    Composite::invertTensor(B) ;
}

BiphasicSelfConsistentComposite::BiphasicSelfConsistentComposite( DelaunayTriangle *tri, Feature *inc ) : MatrixInclusionComposite( tri, inc )
{
	crystals.first = matrix ;
	crystals.second = inclusion ;
	VoigtMatrixInclusionComposite voigt( matrix, inclusion ) ;
	voigt.apply() ;
	fictious = voigt ;
}

BiphasicSelfConsistentComposite::BiphasicSelfConsistentComposite( DelaunayTetrahedron *tet, Feature *inc ) : MatrixInclusionComposite( tet, inc )
{
	crystals.first = matrix ;
	crystals.second = inclusion ;
	VoigtMatrixInclusionComposite voigt( matrix, inclusion ) ;
	voigt.apply() ;
	fictious = voigt ;
}

BiphasicSelfConsistentComposite::BiphasicSelfConsistentComposite( Phase mat, Phase inc ) : MatrixInclusionComposite( mat, inc )
{
	crystals.first = matrix ;
	crystals.second = inclusion ;
	VoigtMatrixInclusionComposite voigt( matrix, inclusion ) ;
	voigt.apply() ;
	fictious = static_cast<Phase>( voigt ) ;
}

void BiphasicSelfConsistentComposite::getStrainConcentrationTensor()
{
	double error = 1. ;
	Matrix S = inclusion.C ;
	Matrix I = Composite::I4( S ) ;

	while( error > POINT_TOLERANCE_2D )
	{
		MoriTanakaMatrixInclusionComposite mtFictiousFirst( fictious, crystals.first ) ;
		mtFictiousFirst.apply() ;
		MoriTanakaMatrixInclusionComposite mtFictiousSecond( fictious, crystals.second ) ;
		mtFictiousSecond.apply() ;

		Matrix G = mtFictiousFirst.inclusion.A * crystals.first.volume ;
		G += mtFictiousSecond.inclusion.A * crystals.second.volume ;
		Composite::invertTensor( G ) ;
		G *= mtFictiousSecond.inclusion.A ;
		S = G ;

		G *= ( crystals.second.C - crystals.first.C ) * crystals.second.volume ;
		G += crystals.first.C ;

		Matrix K = fictious.C - G ;
		error = 0 ;

		for( size_t i = 0 ; i < K.numCols() ; i++ )
		{
			for( size_t j = 0 ; j < K.numRows() ; j++ )
			{
				if( fictious.C[i][j] != 0 )
					K[i][j] /= fictious.C[i][j] ;

				if( std::abs( K[i][j] ) > error )
					error = std::abs( K[i][j] ) ;
			}
		}

		fictious.C = G ;
	}

	B = I - S * inclusion.volume ;
	Composite::invertTensor( B ) ;
	B *= S * matrix.volume ;
}


MatrixMultiInclusionComposite::MatrixMultiInclusionComposite( DelaunayTriangle *tri, std::vector<Feature *> inc ) : Composite( tri, inc )
{

	matrix = Phase( tri ) ;

	for( size_t i = 0 ; i < inc.size() ; i++ )
	{
		if( inclusions.empty() )
			inclusions.push_back( Phase( inc[i] ) ) ;

		else
		{
			int found = -1 ;
			Phase test( inc[i] ) ;

			for( size_t j = 0 ; j < inclusions.size() ; j++ )
			{
				if( test.C == inclusions[j].C && test.beta[0] == inclusions[j].beta[0] )
					found = ( int ) j ;
			}

			if( found == -1 )
				inclusions.push_back( test );
			else
				inclusions[found].volume += test.volume ;
		}
	}


	volume = matrix.volume ;

	for( size_t i = 0 ; i < inclusions.size() ; i++ )
		matrix.volume -= inclusions[i].volume ;

	matrix.volume /= volume ;

	for( size_t i = 0 ; i < inclusions.size() ; i++ )
		inclusions[i].volume /= volume ;

	for( size_t i = 0 ; i < inclusions.size() ; i++ )
	{
		grains.push_back( MoriTanakaMatrixInclusionComposite( matrix, inclusions[i] ) ) ;
		grains[i].volume /= 1. - matrix.volume ;
	}

}

MatrixMultiInclusionComposite::MatrixMultiInclusionComposite( DelaunayTetrahedron *tet, std::vector<Feature *> inc ) : Composite( tet, inc )
{
	matrix = Phase( tet ) ;

	for( size_t i = 0 ; i < inc.size() ; i++ )
	{
		if( inclusions.empty() )
			inclusions.push_back( Phase( inc[i] ) ) ;
		else
		{
			int found = -1 ;
			Phase test( inc[i] ) ;

			for( size_t j = 0 ; j < inclusions.size() ; j++ )
			{
				if( test.C == inclusions[j].C && test.beta[0] == inclusions[j].beta[0] )
					found = ( int ) j ;
			}

			if( found == -1 )
				inclusions.push_back( test );
			else
				inclusions[found].volume += test.volume ;
		}
	}

	for( size_t i = 0 ; i < inclusions.size() ; i++ )
		matrix.volume -= inclusions[i].volume ;

	volume = matrix.volume ;

	for( size_t i = 0 ; i < inclusions.size() ; i++ )
		volume += inclusions[i].volume ;

	matrix.volume /= volume ;

	for( size_t i = 0 ; i < inclusions.size() ; i++ )
		inclusions[i].volume /= volume ;

	for( size_t i = 0 ; i < inclusions.size() ; i++ )
	{
		grains.push_back( MoriTanakaMatrixInclusionComposite( matrix, inclusions[i] ) ) ;
		grains[i].volume /= 1. - matrix.volume ;
	}

}


void MatrixMultiInclusionComposite::apply()
{
	for( size_t i = 0 ; i < grains.size() ; i++ )
		grains[i].apply() ;

	this->getStrainLocalizationTensor() ;

	C = Matrix( C.numRows(), C.numCols() ) ;

//    std::cerr << "-----" << std::endl ;
	for( size_t i = 0 ; i < grains.size() ; i++ )
		C += grains[i].A * ( grains[i].C * grains[i].volume ) ;

//    grains[0].C.print() ;
	//  std::cerr << grains[0].volume << std::endl ;
	// grains[0].A.print() ;

	beta.resize( beta.size() );

	for( size_t i = 0 ; i < grains.size() ; i++ )
		beta += ( grains[i].A * ( grains[i].volume ) ) * grains[i].beta ;
}

void MatrixMultiInclusionComposite::getStrainLocalizationTensor()
{
	for( size_t i = 0 ; i < grains.size() ; i++ )
		grains[i].A = Matrix( C.numRows(), C.numCols() ) ;
}

VoigtMatrixMultiInclusionComposite::VoigtMatrixMultiInclusionComposite( DelaunayTriangle *tri, std::vector<Feature *> inc ) : MatrixMultiInclusionComposite( tri, inc )
{

}

VoigtMatrixMultiInclusionComposite::VoigtMatrixMultiInclusionComposite( DelaunayTetrahedron *tet, std::vector<Feature *> inc ) : MatrixMultiInclusionComposite( tet, inc )
{

}

VoigtMatrixMultiInclusionComposite::VoigtMatrixMultiInclusionComposite( std::vector<DelaunayTriangle *> tri ) : MatrixMultiInclusionComposite( tri[0], std::vector<Feature *>( 0 ) )
{
	for( size_t i = 0 ; i < tri.size() ; i++ )
		inclusions.push_back( Phase( tri[i] ) ) ;

	matrix.volume = 0 ;
	volume = 0 ;

	for( size_t i = 0 ; i < inclusions.size() ; i++ )
	{
		inclusions[i].volume = tri[i]->area() ;
		volume += inclusions[i].volume ;
	}

	for( size_t i = 0 ; i < inclusions.size() ; i++ )
		inclusions[i].volume /= volume ;

	grains.clear() ;
}

VoigtMatrixMultiInclusionComposite::VoigtMatrixMultiInclusionComposite( std::vector<DelaunayTetrahedron *> tet ) : MatrixMultiInclusionComposite( tet[0], std::vector<Feature *>( 0 ) )
{
	for( size_t i = 0 ; i < tet.size() ; i++ )
		inclusions.push_back( Phase( tet[i] ) ) ;

	matrix.volume = 0 ;
	volume = 0 ;

	for( size_t i = 0 ; i < inclusions.size() ; i++ )
	{
		inclusions[i].volume = tet[i]->volume() ;
		volume += inclusions[i].volume ;
	}

	for( size_t i = 0 ; i < inclusions.size() ; i++ )
		inclusions[i].volume /= volume ;

	grains.clear() ;
}

void VoigtMatrixMultiInclusionComposite::getStrainLocalizationTensor()
{
	for( size_t i = 0 ; i < grains.size() ; i++ )
		grains[i].A = Composite::I4( C ) ;
}

ReussMatrixMultiInclusionComposite::ReussMatrixMultiInclusionComposite( DelaunayTriangle *tri, std::vector<Feature *> inc ) : MatrixMultiInclusionComposite( tri, inc )
{

}

ReussMatrixMultiInclusionComposite::ReussMatrixMultiInclusionComposite( DelaunayTetrahedron *tet, std::vector<Feature *> inc ) : MatrixMultiInclusionComposite( tet, inc )
{

}

void ReussMatrixMultiInclusionComposite::getStrainLocalizationTensor()
{
	C = Matrix( matrix.C.numRows(), matrix.C.numCols() ) ;

	for( size_t i = 0 ; i < grains.size() ; i++ )
	{
		Matrix Cg = grains[i].C ;
		Composite::invertTensor( Cg ) ;
		C += Cg * grains[i].volume ;

	}

	Composite::invertTensor( C ) ;

	for( size_t i = 0 ; i < grains.size() ; i++ )
	{
		Matrix Cg = grains[i].C ;
		Composite::invertTensor( Cg ) ;
		grains[i].A = Cg * C ;

	}

}

GeneralizedSelfConsistentComposite::GeneralizedSelfConsistentComposite( std::vector<DelaunayTriangle *> tri ) : VoigtMatrixMultiInclusionComposite( tri )
{

}

GeneralizedSelfConsistentComposite::GeneralizedSelfConsistentComposite( std::vector<DelaunayTetrahedron *> tet ) : VoigtMatrixMultiInclusionComposite( tet )
{

}

void GeneralizedSelfConsistentComposite::apply()
{
	while( !converged() )
	{
		makeGrains() ;

		for( size_t i = 0 ; i < grains.size() ; i++ )
		{
			matrix.C += grains[i].C * grains[i].volume ;
			matrix.beta += grains[i].beta * grains[i].volume ;
		}
	}
}

void GeneralizedSelfConsistentComposite::makeGrains()
{
	previous = Phase( matrix ) ;

	if( grains.size() == 0 )
	{
		matrix.C = Matrix( inclusions[0].C ) ;

		for( size_t i = 0 ; i < inclusions.size() ; i++ )
		{
			matrix.C += inclusions[i].C * inclusions[i].volume ;
			matrix.beta += inclusions[i].beta * inclusions[i].volume ;
		}

		previous = Phase( matrix ) ;
	}

	for( size_t i = 0 ; i < inclusions.size() ; i++ )
	{
		grains.push_back( MoriTanakaMatrixInclusionComposite( matrix, inclusions[i] ) ) ;
		grains[i].inclusion.volume = inclusions[i].volume ;
		grains[i].matrix.volume = 1 - inclusions[i].volume ;
	}
}

bool GeneralizedSelfConsistentComposite::converged()
{
	if( grains.size() == 0 )
		return false ;

	Matrix epsilon = matrix.C - previous.C ;
	Matrix S = matrix.C ;
	Composite::invertTensor( S ) ;
	epsilon *= S ;

	Vector r = epsilon.array() ;
	double rmax = std::abs( r.max() ) ;
	double rmin = std::abs( r.min() ) ;
	return std::max( rmax, rmin ) < POINT_TOLERANCE_3D ;
}
