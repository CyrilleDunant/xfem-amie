#include "phase.h"
#include "../stiffness_with_imposed_deformation.h"
#include "../../utilities/matrixops.h"

using namespace Mu ;

Phase::Phase()
{
    volume = 0 ;
    C = Matrix(6,6) ;
    beta = Vector(3) ;
    A = Matrix(6,6) ;
}

Phase::Phase(DelaunayTriangle * tri)
{
    behaviour = tri->getBehaviour()->getCopy() ;
    volume = tri->area() ;
    stiffnessFromBehaviour() ;
    expansionFromBehaviour() ;
    A = Matrix(C) ;
}

Phase::Phase(DelaunayTetrahedron * tet)
{
    behaviour = tet->getBehaviour()->getCopy() ;
    volume = tet->volume() ;
    stiffnessFromBehaviour() ;
    expansionFromBehaviour() ;
    A = Matrix(C) ;
}

Phase::Phase(Feature * f)
{
    behaviour = f->getBehaviour()->getCopy() ;
    if(f->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	volume = f->area() ;
    if(f->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
	volume = f->volume() ;
    stiffnessFromBehaviour() ;
    expansionFromBehaviour() ;
    A = Matrix(C) ;
}

Phase::Phase(const Phase & p)
{
    C = p.C ;
    beta = p.beta ;
    A = p.A ;
    volume = p.volume ;
}

Form * Phase::getBehaviour()
{
    Vector alpha ;
    Matrix S = C ;
    if(S.size()==36)
	invert6x6Matrix(S) ;
    else
	invert3x3Matrix(S) ;
    alpha = S * beta ;
    return new StiffnessWithImposedDeformation(C,alpha) ;
}

void Phase::stiffnessFromBehaviour()
{
    C = behaviour->getTensor(Point(1./3,1./3,1./3)) ;
}

void Phase::expansionFromBehaviour()
{
    beta = behaviour->getImposedStress(Point(1./3,1./3,1./3)) ;
}
