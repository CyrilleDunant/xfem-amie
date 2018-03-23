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
    Composite(const Phase & p) ;

    virtual void apply() ;

    static Matrix I4(const Matrix & C) ;
     Matrix eshelby(const Matrix & C) ;
     Matrix eshelbyCylinder(const Matrix & C, double a, double b) ;
     Matrix eshelbyEllipsoid(const Matrix & C, double a, double b, double c) ;
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
    MatrixInclusionComposite(const Phase & mat, const Phase & inc) ;

    virtual void apply() ;
    virtual void getStrainConcentrationTensor() ;
};

struct DiluteMatrixInclusionComposite : public MatrixInclusionComposite
{
public:
    DiluteMatrixInclusionComposite(DelaunayTriangle * tri, Feature * inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    DiluteMatrixInclusionComposite(DelaunayTetrahedron * tet, Feature * inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    DiluteMatrixInclusionComposite(const Phase & mat, const Phase & inc) ;

    virtual void getStrainConcentrationTensor() ;
};

struct VoigtMatrixInclusionComposite : public MatrixInclusionComposite
{
public:
    VoigtMatrixInclusionComposite(DelaunayTriangle * tri, Feature * inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    VoigtMatrixInclusionComposite(DelaunayTetrahedron * tet, Feature * inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    VoigtMatrixInclusionComposite(const Phase &mat, const Phase &inc) ;

    virtual void getStrainConcentrationTensor() ;
};

struct ReussMatrixInclusionComposite : public MatrixInclusionComposite
{
public:
    ReussMatrixInclusionComposite(DelaunayTriangle * tri, Feature * inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    ReussMatrixInclusionComposite(DelaunayTetrahedron * tet, Feature * inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    ReussMatrixInclusionComposite(const Phase &mat, const Phase &inc) ;

    virtual void getStrainConcentrationTensor() ;
};

struct MoriTanakaMatrixInclusionComposite : public MatrixInclusionComposite
{
public:
    MoriTanakaMatrixInclusionComposite(DelaunayTriangle * tri, Feature * inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    MoriTanakaMatrixInclusionComposite(DelaunayTetrahedron * tet, Feature * inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    MoriTanakaMatrixInclusionComposite(const Phase & mat, const Phase &inc) ;

    virtual void getStrainConcentrationTensor() ;
};

struct InverseMoriTanakaMatrixInclusionComposite : public MatrixInclusionComposite
{
public:
    InverseMoriTanakaMatrixInclusionComposite(DelaunayTriangle * tri, Feature * inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    InverseMoriTanakaMatrixInclusionComposite(DelaunayTetrahedron * tet, Feature * inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    InverseMoriTanakaMatrixInclusionComposite(const Phase &mat, const Phase &inc) ;

    virtual void getStrainConcentrationTensor() ;
};

struct BiphasicSelfConsistentComposite : public MatrixInclusionComposite
{
public:
	Phase fictious ;

public:
	BiphasicSelfConsistentComposite(DelaunayTriangle * tri, Feature * inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
	BiphasicSelfConsistentComposite(DelaunayTetrahedron * tet, Feature * inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
	BiphasicSelfConsistentComposite(const Phase &mat, const Phase &inc) ;
    BiphasicSelfConsistentComposite(const Phase &mat,const  Phase &inc, const BiphasicSelfConsistentComposite & hint) ;

	virtual void getStrainConcentrationTensor() ;
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
    MatrixMultiInclusionComposite( const Phase &m, const std::vector<Phase> &inc) ;

    virtual void apply() ;
    virtual void getStrainLocalizationTensor(InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
};

struct VoigtMatrixMultiInclusionComposite : public MatrixMultiInclusionComposite
{
public:
    VoigtMatrixMultiInclusionComposite(DelaunayTriangle * tri, std::vector<Feature *> inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    VoigtMatrixMultiInclusionComposite(DelaunayTetrahedron * tet, std::vector<Feature *> inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    VoigtMatrixMultiInclusionComposite(std::vector<DelaunayTriangle *> tri, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    VoigtMatrixMultiInclusionComposite(std::vector<DelaunayTetrahedron *> tet, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    VoigtMatrixMultiInclusionComposite( const Phase &m, std::vector<Phase> &inc) : MatrixMultiInclusionComposite(m, inc) { } ;

    virtual void getStrainLocalizationTensor(InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
};

struct ReussMatrixMultiInclusionComposite : public MatrixMultiInclusionComposite
{
public:
    ReussMatrixMultiInclusionComposite(DelaunayTriangle * tri, std::vector<Feature *> inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    ReussMatrixMultiInclusionComposite(DelaunayTetrahedron * tet, std::vector<Feature *> inc, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    ReussMatrixMultiInclusionComposite( const Phase &m, std::vector<Phase> &inc)  : MatrixMultiInclusionComposite(m, inc) { } ;

    virtual void getStrainLocalizationTensor(InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
};

struct GeneralizedSelfConsistentComposite : public VoigtMatrixMultiInclusionComposite
{
public:
    GeneralizedSelfConsistentComposite(std::vector<DelaunayTriangle *> tri, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
    GeneralizedSelfConsistentComposite(std::vector<DelaunayTetrahedron *> tet, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;

    virtual void apply() ;

private:
    Phase previous ;
    void makeGrains() ;
    bool converged() ;
};

}

#endif // COMPOSITE_H
