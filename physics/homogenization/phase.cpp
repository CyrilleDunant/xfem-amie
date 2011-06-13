#include "phase.h"
#include "../stiffness_with_imposed_deformation.h"
#include "../../utilities/matrixops.h"

using namespace Mu ;

Phase::Phase()
{
	volume = 0 ;
	C = Matrix( 6, 6 ) ;
	beta = Vector( 3 ) ;
	A = Matrix( 6, 6 ) ;
}

Phase::Phase( DelaunayTriangle *tri )
{
	behaviour = tri->getBehaviour()->getCopy() ;
	volume = tri->area() ;
	stiffnessFromBehaviour() ;
	expansionFromBehaviour() ;
	A = Matrix( C.numRows(), C.numCols() ) ;
}

Phase::Phase( DelaunayTetrahedron *tet )
{
	behaviour = tet->getBehaviour()->getCopy() ;
	volume = tet->volume() ;
	stiffnessFromBehaviour() ;
	expansionFromBehaviour() ;
	A = Matrix( C.numRows(), C.numCols() ) ;
}

Phase::Phase( Feature *f )
{
	behaviour = f->getBehaviour()->getCopy() ;

	if( f->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
		volume = f->area() ;

	if( f->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
		volume = f->volume() ;

	stiffnessFromBehaviour() ;
	expansionFromBehaviour() ;
	A = Matrix( C.numRows(), C.numCols() ) ;
}

Phase::Phase( const Phase &p )
{
	C.resize( p.C.numRows(), p.C.numCols() ) ;
	C = p.C ;
	beta.resize( p.beta.size() ) ;
	beta = p.beta ;
	A.resize( p.A.numRows(), p.A.numCols() ) ;
	A = p.A ;
	volume = p.volume ;
}

void Phase::apply()
{

}

Form *Phase::getBehaviour()
{
	this->apply() ;
	Matrix S = C ;

	if( S.size() == 36 )
		invert6x6Matrix( S ) ;
	else
		invert3x3Matrix( S ) ;

	Vector alpha = S * beta ;
	return new StiffnessWithImposedDeformation( C, alpha ) ;
}

Phase &Phase::operator =( const Phase &p )
{
	C.resize( p.C.numRows(), p.C.numCols() ) ;
	C = p.C ;
	beta.resize( p.beta.size() ) ;
	beta = p.beta ;
	A.resize( p.A.numRows(), p.A.numCols() ) ;
	A = p.A ;
	volume = p.volume ;
	return *this ;
}

void Phase::stiffnessFromBehaviour()
{
	Matrix tmp =  behaviour->getTensor( Point( 1. / 3, 1. / 3, 1. / 3 ) ) ;
	C.resize( tmp.numRows(), tmp.numCols() ) ;
	C = tmp ;
}

void Phase::expansionFromBehaviour()
{
	Vector tmp = behaviour->getImposedStress( Point( 1. / 3, 1. / 3, 1. / 3 ) ) ;
	beta.resize( tmp.size() );
	beta = tmp ;
}
