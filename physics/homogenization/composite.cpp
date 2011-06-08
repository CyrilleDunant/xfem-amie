#include "composite.h"
#include "../../utilities/matrixops.h"

using namespace Mu ;

Composite::Composite(DelaunayTriangle * tri, std::vector<Feature *> feats) : Phase(tri)
{

}

Composite::Composite(DelaunayTetrahedron * tet, std::vector<Feature *> feats) : Phase(tet)
{

}

Composite::Composite(Phase p) : Phase(p)
{

}

void Composite::apply()
{

}

Matrix Composite::I4(Matrix C)
{
    Matrix I4(C) ;
    for(size_t i = 0 ; i < 2+(I4.size()==36) ; i++)
    {
	I4[i][i] = 0 ;
    }
    return I4 ;
}

Matrix Composite::eshelby(Matrix C)
{
    Matrix S(C) ;
    double nu = 0.2 ;
    double a = C[0][0] ;
    double b = C[0][1] ;
    if(C.size()==36)
	nu = b / (a+b) ;
    else
	nu = b / a ;

    double Siiii = (7-5*nu)/(15*(1-nu)) ;
    double Siijj = (5*nu-1)/(15*(1-nu)) ;
    double Sijij = (4-5*nu)/(15*(1-nu)) ;

    if(S.size() == 9)
    {
	S[0][0] = Siiii ;
	S[1][1] = Siiii ;

	S[0][1] = Siijj ;
	S[1][0] = Siijj ;

	S[2][2] = Sijij ;
    }

    if(S.size() == 36)
    {
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

MatrixInclusionComposite::MatrixInclusionComposite(DelaunayTriangle *tri, Feature *inc) : Composite(tri, std::vector<Feature *>(0))
{
    matrix = Phase(tri) ;
    inclusion = Phase(inc) ;

    matrix.volume -= inclusion.volume ;
    volume = matrix.volume + inclusion.volume ;
    matrix.volume /= volume ;
    inclusion.volume /= volume ;
}

MatrixInclusionComposite::MatrixInclusionComposite(DelaunayTetrahedron * tet, Feature * inc) : Composite(tet, std::vector<Feature *>(0))
{
    matrix = Phase(tet) ;
    inclusion = Phase(inc) ;

    matrix.volume -= inclusion.volume ;
    volume = matrix.volume + inclusion.volume ;
    matrix.volume /= volume ;
    inclusion.volume /= volume ;
}

MatrixInclusionComposite::MatrixInclusionComposite(Phase mat, Phase inc) : Composite(mat)
{
    matrix = Phase(mat) ;
    inclusion = Phase(inc) ;

    volume = inclusion.volume ;

    inclusion.volume = 1-matrix.volume ;

}

void MatrixInclusionComposite::apply()
{
    this->getStrainConcentrationTensor() ;
    Matrix I = Composite::I4(matrix.C) ;

    inclusion.A = B * inclusion.volume + I * matrix.volume ;
    if(A.size()==36)
	invert6x6Matrix(inclusion.A) ;
    else
	invert3x3Matrix(inclusion.A) ;
    inclusion.A *= B ;

    matrix.A = (I - inclusion.A * inclusion.volume) / matrix.volume ;

    C = ((matrix.A * matrix.volume) * matrix.C) ;
    C += ((inclusion.A * inclusion.volume) * inclusion.C) ;
    beta.resize(beta.size());
    beta = matrix.A * matrix.volume * matrix.beta ;
    beta += inclusion.A * inclusion.volume * inclusion.beta ;

}

void MatrixInclusionComposite::getStrainConcentrationTensor()
{
    B = Matrix(C) ;
}

VoigtMatrixInclusionComposite::VoigtMatrixInclusionComposite(DelaunayTriangle *tri, Feature *inc) : MatrixInclusionComposite(tri, inc)
{

}

VoigtMatrixInclusionComposite::VoigtMatrixInclusionComposite(DelaunayTetrahedron *tet, Feature *inc) : MatrixInclusionComposite(tet, inc)
{

}

void VoigtMatrixInclusionComposite::getStrainConcentrationTensor()
{
    B = Composite::I4(C) ;
}

ReussMatrixInclusionComposite::ReussMatrixInclusionComposite(DelaunayTriangle *tri, Feature *inc) : MatrixInclusionComposite(tri, inc)
{

}

ReussMatrixInclusionComposite::ReussMatrixInclusionComposite(DelaunayTetrahedron *tet, Feature *inc) : MatrixInclusionComposite(tet, inc)
{

}

void ReussMatrixInclusionComposite::getStrainConcentrationTensor()
{
    B = inclusion.C ;
    if(B.size()==36)
	invert6x6Matrix(B) ;
    else
	invert3x3Matrix(B) ;
    B *= matrix.C ;
}

MoriTanakaMatrixInclusionComposite::MoriTanakaMatrixInclusionComposite(DelaunayTriangle *tri, Feature *inc) : MatrixInclusionComposite(tri, inc)
{

}

MoriTanakaMatrixInclusionComposite::MoriTanakaMatrixInclusionComposite(DelaunayTetrahedron *tet, Feature *inc) : MatrixInclusionComposite(tet, inc)
{

}

MoriTanakaMatrixInclusionComposite::MoriTanakaMatrixInclusionComposite(Phase mat, Phase inc) : MatrixInclusionComposite(mat, inc)
{

}

void MoriTanakaMatrixInclusionComposite::getStrainConcentrationTensor()
{
    Matrix I = Composite::I4(C) ;
    Matrix S = Composite::eshelby(inclusion.C) ;

    B = matrix.C ;
    if(B.size()==36)
	invert6x6Matrix(B) ;
    else
	invert6x6Matrix(B) ;
    B *= S * (inclusion.C - matrix.C) ;
    B += I ;
    if(B.size()==36)
	invert6x6Matrix(B) ;
    else
	invert3x3Matrix(B) ;

}

MatrixMultiInclusionComposite::MatrixMultiInclusionComposite(DelaunayTriangle *tri, std::vector<Feature *> inc) : Composite(tri, inc)
{
    matrix = Phase(tri) ;
    for(size_t i = 0 ; i < inc.size() ; i++)
    {
	if(inclusions.empty())
	    inclusions.push_back(Phase(inc[i])) ;
	else
	{
	    int found = -1 ;
	    Phase test(inc[i]) ;
	    for(size_t j = 0 ; j < inclusions.size() ; j++)
	    {
		if(test.C == inclusions[j].C && test.beta[0] == inclusions[j].beta[0])
		    found = (int) j ;
	    }
	    if(found == -1)
		inclusions.push_back(test);
	    else
		inclusions[found].volume += test.volume ;
	}
    }

    for(size_t i = 0 ; i < inclusions.size() ; i++)
	matrix.volume -= inclusions[i].volume ;

    volume = matrix.volume ;
    for(size_t i = 0 ; i < inclusions.size() ; i++)
	volume += inclusions[i].volume ;

    matrix.volume /= volume ;
    for(size_t i = 0 ; i < inclusions.size() ; i++)
	inclusions[i].volume /= volume ;

    for(size_t i = 0 ; i < inclusions.size() ; i++)
	grains.push_back(MoriTanakaMatrixInclusionComposite(matrix, inclusions[i])) ;

}

MatrixMultiInclusionComposite::MatrixMultiInclusionComposite(DelaunayTetrahedron *tet, std::vector<Feature *> inc) : Composite(tet, inc)
{
    matrix = Phase(tet) ;
    for(size_t i = 0 ; i < inc.size() ; i++)
    {
	if(inclusions.empty())
	    inclusions.push_back(Phase(inc[i])) ;
	else
	{
	    int found = -1 ;
	    Phase test(inc[i]) ;
	    for(size_t j = 0 ; j < inclusions.size() ; j++)
	    {
		if(test.C == inclusions[j].C && test.beta[0] == inclusions[j].beta[0])
		    found = (int) j ;
	    }
	    if(found == -1)
		inclusions.push_back(test);
	    else
		inclusions[found].volume += test.volume ;
	}
    }

    for(size_t i = 0 ; i < inclusions.size() ; i++)
	matrix.volume -= inclusions[i].volume ;

    volume = matrix.volume ;
    for(size_t i = 0 ; i < inclusions.size() ; i++)
	volume += inclusions[i].volume ;

    matrix.volume /= volume ;
    for(size_t i = 0 ; i < inclusions.size() ; i++)
	inclusions[i].volume /= volume ;

    for(size_t i = 0 ; i < inclusions.size() ; i++)
	grains.push_back(MoriTanakaMatrixInclusionComposite(matrix, inclusions[i])) ;

}


void MatrixMultiInclusionComposite::apply()
{
    for(size_t i = 0 ; i < grains.size() ; i++)
	grains[i].apply() ;

    getStrainLocalizationTensor() ;

    C = Matrix(C) ;
    for(size_t i = 0 ; i < grains.size() ; i++)
	C += grains[i].A * (grains[i].C * grains[i].volume) ;

    beta.resize(beta.size());
    for(size_t i = 0 ; i < grains.size() ; i++)
	beta += (grains[i].A * grains[i].volume) * grains[i].beta ;
}

void MatrixMultiInclusionComposite::getStrainLocalizationTensor()
{
    for(size_t i = 0 ; i < grains.size() ; i++)
	grains[i].A = Matrix(C) ;
}

VoigtMatrixMultiInclusionComposite::VoigtMatrixMultiInclusionComposite(DelaunayTriangle *tri, std::vector<Feature *> inc) : MatrixMultiInclusionComposite(tri, inc)
{

}

VoigtMatrixMultiInclusionComposite::VoigtMatrixMultiInclusionComposite(DelaunayTetrahedron * tet, std::vector<Feature *> inc) : MatrixMultiInclusionComposite(tet, inc)
{

}

VoigtMatrixMultiInclusionComposite::VoigtMatrixMultiInclusionComposite(std::vector<DelaunayTriangle *> tri) : MatrixMultiInclusionComposite(tri[0], std::vector<Feature *>(0))
{
    for(size_t i = 0 ; i < tri.size() ; i++)
	inclusions.push_back(Phase(tri[i])) ;

    matrix.volume = 0 ;
    volume = 0 ;
    for(size_t i = 0 ; i < inclusions.size() ; i++)
    {
	inclusions[i].volume = tri[i]->area() ;
	volume += inclusions[i].volume ;
    }

    for(size_t i = 0 ; i < inclusions.size() ; i++)
	inclusions[i].volume /= volume ;

    grains.clear() ;
}

VoigtMatrixMultiInclusionComposite::VoigtMatrixMultiInclusionComposite(std::vector<DelaunayTetrahedron *> tet) : MatrixMultiInclusionComposite(tet[0], std::vector<Feature *>(0))
{
    for(size_t i = 0 ; i < tet.size() ; i++)
	inclusions.push_back(Phase(tet[i])) ;

    matrix.volume = 0 ;
    volume = 0 ;
    for(size_t i = 0 ; i < inclusions.size() ; i++)
    {
	inclusions[i].volume = tet[i]->volume() ;
	volume += inclusions[i].volume ;
    }

    for(size_t i = 0 ; i < inclusions.size() ; i++)
	inclusions[i].volume /= volume ;

    grains.clear() ;
}

void VoigtMatrixMultiInclusionComposite::getStrainLocalizationTensor()
{
    for(size_t i = 0 ; i < grains.size() ; i++)
	grains[i].A = Composite::I4(C) ;
}

ReussMatrixMultiInclusionComposite::ReussMatrixMultiInclusionComposite(DelaunayTriangle *tri, std::vector<Feature *> inc) : MatrixMultiInclusionComposite(tri, inc)
{

}

ReussMatrixMultiInclusionComposite::ReussMatrixMultiInclusionComposite(DelaunayTetrahedron *tet, std::vector<Feature *> inc) : MatrixMultiInclusionComposite(tet, inc)
{

}

void ReussMatrixMultiInclusionComposite::getStrainLocalizationTensor()
{
    C = Matrix(matrix.C) ;
    for(size_t i = 0 ; i < grains.size() ; i++)
    {
	Matrix Cg = grains[i].C ;
	if(Cg.size()==36)
	    invert6x6Matrix(Cg) ;
	else
	    invert3x3Matrix(Cg) ;

	C += Cg * grains[i].volume ;

    }

    if(C.size()==36)
	invert6x6Matrix(C) ;
    else
	invert3x3Matrix(C) ;

    for(size_t i = 0 ; i < grains.size() ; i++)
    {
	Matrix Cg = grains[i].C ;
	if(Cg.size()==36)
	    invert6x6Matrix(Cg) ;
	else
	    invert3x3Matrix(Cg) ;

	grains[i].A = Cg * C ;

    }

}

GeneralizedSelfConsistentComposite::GeneralizedSelfConsistentComposite(std::vector<DelaunayTriangle *> tri) : VoigtMatrixMultiInclusionComposite(tri)
{

}

GeneralizedSelfConsistentComposite::GeneralizedSelfConsistentComposite(std::vector<DelaunayTetrahedron *> tet) : VoigtMatrixMultiInclusionComposite(tet)
{

}

void GeneralizedSelfConsistentComposite::apply()
{
    while(!converged())
    {
	makeGrains() ;
	for(size_t i = 0 ; i < grains.size() ; i++)
	{
	    matrix.C += grains[i].C * grains[i].volume ;
	    matrix.beta += grains[i].beta* grains[i].volume ;
	}
    }
}

void GeneralizedSelfConsistentComposite::makeGrains()
{
    previous = Phase(matrix) ;
    if(grains.size() == 0)
    {
	matrix.C = Matrix(inclusions[0].C) ;
	for(size_t i = 0 ; i < inclusions.size() ; i++)
	{
	    matrix.C += inclusions[i].C * inclusions[i].volume ;
	    matrix.beta += inclusions[i].beta* inclusions[i].volume ;
	}
	previous = Phase(matrix) ;
    }

    for(size_t i = 0 ; i < inclusions.size() ; i++)
    {
	grains.push_back(MoriTanakaMatrixInclusionComposite(matrix,inclusions[i])) ;
	grains[i].inclusion.volume = inclusions[i].volume ;
	grains[i].matrix.volume = 1-inclusions[i].volume ;
    }
}

bool GeneralizedSelfConsistentComposite::converged()
{
    if(grains.size() == 0)
	return false ;
    Matrix epsilon = matrix.C - previous.C ;
    Matrix S = matrix.C ;
    if(S.size()==36)
	invert6x6Matrix(S) ;
    else
	invert3x3Matrix(S) ;
    epsilon *= S ;

    Vector r = epsilon.array() ;
    double rmax =std::abs(r.max()) ;
    double rmin = std::abs(r.min()) ;
    std::cout << std::max(rmax,rmin) << std::endl ;
    return std::max(rmax,rmin) < POINT_TOLERANCE_3D ;
}
