#ifndef COMPOSITE_H
#define COMPOSITE_H

#include "phase.h"

namespace Amie
{

struct Composite : public Phase
{
public:
    Composite(DelaunayTriangle * tri, std::vector<Feature *> feats);
    Composite(DelaunayTetrahedron * tet, std::vector<Feature *> feats);
    Composite(Phase p) ;

    virtual void apply() ;

    static Matrix I4(Matrix C) ;
    static Matrix eshelby(Matrix C) ;
    static void invertTensor(Matrix & m) ;
} ;

struct MatrixInclusionComposite : public Composite
{
public:
    Phase matrix ;
    Phase inclusion ;
    Matrix B ;

public:
    MatrixInclusionComposite(DelaunayTriangle * tri, Feature * inc) ;
    MatrixInclusionComposite(DelaunayTetrahedron * tet, Feature * inc) ;
    MatrixInclusionComposite(Phase mat, Phase inc) ;

    virtual void apply() ;
    virtual void getStrainConcentrationTensor() ;
};

struct DiluteMatrixInclusionComposite : public MatrixInclusionComposite
{
public:
    DiluteMatrixInclusionComposite(DelaunayTriangle * tri, Feature * inc) ;
    DiluteMatrixInclusionComposite(DelaunayTetrahedron * tet, Feature * inc) ;
    DiluteMatrixInclusionComposite(Phase mat, Phase inc) ;

    virtual void getStrainConcentrationTensor() ;
};

struct VoigtMatrixInclusionComposite : public MatrixInclusionComposite
{
public:
    VoigtMatrixInclusionComposite(DelaunayTriangle * tri, Feature * inc) ;
    VoigtMatrixInclusionComposite(DelaunayTetrahedron * tet, Feature * inc) ;
    VoigtMatrixInclusionComposite(Phase mat, Phase inc) ;

    virtual void getStrainConcentrationTensor() ;
};

struct ReussMatrixInclusionComposite : public MatrixInclusionComposite
{
public:
    ReussMatrixInclusionComposite(DelaunayTriangle * tri, Feature * inc) ;
    ReussMatrixInclusionComposite(DelaunayTetrahedron * tet, Feature * inc) ;
    ReussMatrixInclusionComposite(Phase mat, Phase inc) ;

    virtual void getStrainConcentrationTensor() ;
};

struct MoriTanakaMatrixInclusionComposite : public MatrixInclusionComposite
{
public:
    MoriTanakaMatrixInclusionComposite(DelaunayTriangle * tri, Feature * inc) ;
    MoriTanakaMatrixInclusionComposite(DelaunayTetrahedron * tet, Feature * inc) ;
    MoriTanakaMatrixInclusionComposite(Phase mat, Phase inc) ;

    virtual void getStrainConcentrationTensor() ;
};

struct InverseMoriTanakaMatrixInclusionComposite : public MatrixInclusionComposite
{
public:
    InverseMoriTanakaMatrixInclusionComposite(DelaunayTriangle * tri, Feature * inc) ;
    InverseMoriTanakaMatrixInclusionComposite(DelaunayTetrahedron * tet, Feature * inc) ;
    InverseMoriTanakaMatrixInclusionComposite(Phase mat, Phase inc) ;

    virtual void getStrainConcentrationTensor() ;
};

struct BiphasicSelfConsistentComposite : public MatrixInclusionComposite
{
public:
	std::pair<Phase, Phase> crystals ;
	Phase fictious ;

public:
	BiphasicSelfConsistentComposite(DelaunayTriangle * tri, Feature * inc) ;
	BiphasicSelfConsistentComposite(DelaunayTetrahedron * tet, Feature * inc) ;
	BiphasicSelfConsistentComposite(Phase mat, Phase inc) ;

	virtual void getStrainConcentrationTensor() ;
};

struct MatrixMultiInclusionComposite : public Composite
{
public:
    Phase matrix ;
    std::vector<Phase> inclusions ;
    std::vector<MoriTanakaMatrixInclusionComposite> grains ;

public:
    MatrixMultiInclusionComposite(DelaunayTriangle * tri, std::vector<Feature *> inc) ;
    MatrixMultiInclusionComposite(DelaunayTetrahedron * tet, std::vector<Feature *> inc) ;

    virtual void apply() ;
    virtual void getStrainLocalizationTensor() ;
};

struct VoigtMatrixMultiInclusionComposite : public MatrixMultiInclusionComposite
{
public:
    VoigtMatrixMultiInclusionComposite(DelaunayTriangle * tri, std::vector<Feature *> inc) ;
    VoigtMatrixMultiInclusionComposite(DelaunayTetrahedron * tet, std::vector<Feature *> inc) ;
    VoigtMatrixMultiInclusionComposite(std::vector<DelaunayTriangle *> tri) ;
    VoigtMatrixMultiInclusionComposite(std::vector<DelaunayTetrahedron *> tet) ;

    virtual void getStrainLocalizationTensor() ;
};

struct ReussMatrixMultiInclusionComposite : public MatrixMultiInclusionComposite
{
public:
    ReussMatrixMultiInclusionComposite(DelaunayTriangle * tri, std::vector<Feature *> inc) ;
    ReussMatrixMultiInclusionComposite(DelaunayTetrahedron * tet, std::vector<Feature *> inc) ;

    virtual void getStrainLocalizationTensor() ;
};

struct GeneralizedSelfConsistentComposite : public VoigtMatrixMultiInclusionComposite
{
public:
    GeneralizedSelfConsistentComposite(std::vector<DelaunayTriangle *> tri) ;
    GeneralizedSelfConsistentComposite(std::vector<DelaunayTetrahedron *> tet) ;

    virtual void apply() ;

private:
    Phase previous ;
    void makeGrains() ;
    bool converged() ;
};

}

#endif // COMPOSITE_H
