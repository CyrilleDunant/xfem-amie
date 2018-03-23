#ifndef COMPOSITE_H
#define COMPOSITE_H

#include "phase.h"

namespace Amie
{

    
struct Composite : public Phase
{
public:
    Composite(DelaunayTriangle * tri, std::vector<Feature *> feats, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1);
    Composite(DelaunayTetrahedron * tet, std::vector<Feature *> feats, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1);
    Composite(Phase p, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;

    virtual void apply(InclusionGeometryType t, double a=1, double b=1, double c=1) ;

    static Matrix I4(Matrix C) ;
     Matrix eshelby(Matrix C) ;
     Matrix eshelbyCylinder(Matrix C, double a, double b) ;
     Matrix eshelbyEllipsoid(Matrix C, double a, double b, double c) ;
    static void invertTensor(Matrix & m) ;
} ;

struct MatrixInclusionComposite : public Composite
{
public:
    Phase matrix ;
    Phase inclusion ;
    Matrix B ;

public:
    MatrixInclusionComposite(DelaunayTriangle * tri, Feature * inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    MatrixInclusionComposite(DelaunayTetrahedron * tet, Feature * inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    MatrixInclusionComposite(Phase mat, Phase inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;

    virtual void apply(InclusionGeometryType t, double a=1, double b=1, double c=1) ;
    virtual void getStrainConcentrationTensor(InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
};

struct DiluteMatrixInclusionComposite : public MatrixInclusionComposite
{
public:
    DiluteMatrixInclusionComposite(DelaunayTriangle * tri, Feature * inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    DiluteMatrixInclusionComposite(DelaunayTetrahedron * tet, Feature * inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    DiluteMatrixInclusionComposite(Phase mat, Phase inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;

    virtual void getStrainConcentrationTensor(InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
};

struct VoigtMatrixInclusionComposite : public MatrixInclusionComposite
{
public:
    VoigtMatrixInclusionComposite(DelaunayTriangle * tri, Feature * inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    VoigtMatrixInclusionComposite(DelaunayTetrahedron * tet, Feature * inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    VoigtMatrixInclusionComposite(Phase mat, Phase inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;

    virtual void getStrainConcentrationTensor(InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
};

struct ReussMatrixInclusionComposite : public MatrixInclusionComposite
{
public:
    ReussMatrixInclusionComposite(DelaunayTriangle * tri, Feature * inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    ReussMatrixInclusionComposite(DelaunayTetrahedron * tet, Feature * inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    ReussMatrixInclusionComposite(Phase mat, Phase inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;

    virtual void getStrainConcentrationTensor(InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
};

struct MoriTanakaMatrixInclusionComposite : public MatrixInclusionComposite
{
public:
    MoriTanakaMatrixInclusionComposite(DelaunayTriangle * tri, Feature * inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    MoriTanakaMatrixInclusionComposite(DelaunayTetrahedron * tet, Feature * inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    MoriTanakaMatrixInclusionComposite(Phase mat, Phase inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;

    virtual void getStrainConcentrationTensor(InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
};

struct InverseMoriTanakaMatrixInclusionComposite : public MatrixInclusionComposite
{
public:
    InverseMoriTanakaMatrixInclusionComposite(DelaunayTriangle * tri, Feature * inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    InverseMoriTanakaMatrixInclusionComposite(DelaunayTetrahedron * tet, Feature * inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    InverseMoriTanakaMatrixInclusionComposite(Phase mat, Phase inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;

    virtual void getStrainConcentrationTensor(InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
};

struct BiphasicSelfConsistentComposite : public MatrixInclusionComposite
{
public:
	std::pair<Phase, Phase> crystals ;
	Phase fictious ;

public:
	BiphasicSelfConsistentComposite(DelaunayTriangle * tri, Feature * inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
	BiphasicSelfConsistentComposite(DelaunayTetrahedron * tet, Feature * inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
	BiphasicSelfConsistentComposite(Phase mat, Phase inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    BiphasicSelfConsistentComposite(Phase mat, Phase inc, BiphasicSelfConsistentComposite hint, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;

	virtual void getStrainConcentrationTensor(InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
};

struct MatrixMultiInclusionComposite : public Composite
{
public:
    Phase matrix ;
    std::vector<Phase> inclusions ;
    std::vector<MoriTanakaMatrixInclusionComposite> grains ;

public:
    MatrixMultiInclusionComposite(DelaunayTriangle * tri, std::vector<Feature *> inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    MatrixMultiInclusionComposite(DelaunayTetrahedron * tet, std::vector<Feature *> inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    MatrixMultiInclusionComposite( Phase m, std::vector<Phase> inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1 ) ;

    virtual void apply(InclusionGeometryType t, double a=1, double b=1, double c=1) ;
    virtual void getStrainLocalizationTensor() ;
};

struct VoigtMatrixMultiInclusionComposite : public MatrixMultiInclusionComposite
{
public:
    VoigtMatrixMultiInclusionComposite(DelaunayTriangle * tri, std::vector<Feature *> inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    VoigtMatrixMultiInclusionComposite(DelaunayTetrahedron * tet, std::vector<Feature *> inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    VoigtMatrixMultiInclusionComposite(std::vector<DelaunayTriangle *> tri, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    VoigtMatrixMultiInclusionComposite(std::vector<DelaunayTetrahedron *> tet, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    VoigtMatrixMultiInclusionComposite( Phase m, std::vector<Phase> inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1 ) : MatrixMultiInclusionComposite(m, inc, t, a, b, c) { } ;

    virtual void getStrainLocalizationTensor() ;
};

struct ReussMatrixMultiInclusionComposite : public MatrixMultiInclusionComposite
{
public:
    ReussMatrixMultiInclusionComposite(DelaunayTriangle * tri, std::vector<Feature *> inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    ReussMatrixMultiInclusionComposite(DelaunayTetrahedron * tet, std::vector<Feature *> inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    ReussMatrixMultiInclusionComposite( Phase m, std::vector<Phase> inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1 )  : MatrixMultiInclusionComposite(m, inc, t, a, b, c) { } ;

    virtual void getStrainLocalizationTensor() ;
};

struct GeneralizedSelfConsistentComposite : public VoigtMatrixMultiInclusionComposite
{
public:
    GeneralizedSelfConsistentComposite(std::vector<DelaunayTriangle *> tri, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    GeneralizedSelfConsistentComposite(std::vector<DelaunayTetrahedron *> tet, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;

    virtual void apply(InclusionGeometryType t, double a=1, double b=1, double c=1) ;

private:
    Phase previous ;
    void makeGrains(InclusionGeometryType t, double a=1, double b=1, double c=1) ;
    bool converged() ;
};

}

#endif // COMPOSITE_H
