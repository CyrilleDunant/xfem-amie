#include "phase.h"
#include "../stiffness_with_imposed_deformation.h"
#include "../../utilities/matrixops.h"
#include "../fracturecriteria/mohrcoulomb.h"

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
	ruptureFromBehaviour() ;
	A = Matrix( C.numRows(), C.numCols() ) ;
}

Phase::Phase( DelaunayTetrahedron *tet )
{
	behaviour = tet->getBehaviour()->getCopy() ;
	volume = tet->volume() ;
	stiffnessFromBehaviour() ;
	expansionFromBehaviour() ;
	ruptureFromBehaviour() ;
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
	ruptureFromBehaviour() ;
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
	for(size_t i = 0 ; i < p.lambda.size() ; i++)
	{
	    Vector l(p.lambda[i].size()) ;
	    l = p.lambda[i] ;
	    lambda.push_back(l);
	}
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
	for(size_t i = 0 ; i < p.lambda.size() ; i++)
	{
	    Vector l(p.lambda[i].size()) ;
	    l = p.lambda[i] ;
	    lambda.push_back(l);
	}
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

void Phase::ruptureFromBehaviour()
{
    lambda.clear() ;
    FractureCriterion * frac = behaviour->getFractureCriterion() ;
    if(frac)
    {
	if(dynamic_cast<NonLocalMohrCoulomb *>(frac))
	{
	    double u = dynamic_cast<NonLocalMohrCoulomb *>(frac)->upVal ;
	    double d = dynamic_cast<NonLocalMohrCoulomb *>(frac)->downVal ;
	    Vector up(beta.size()) ;
	    Vector down(beta.size()) ;
	    for(size_t i = 0 ; i < beta.size() ; i++)
	    {
		up[i] = u ;
		down[i] = d ;
	    }
	    lambda.push_back(up) ;
	    lambda.push_back(down) ;
	}
	if(dynamic_cast<MohrCoulomb *>(frac))
	{
	    double u = dynamic_cast<MohrCoulomb *>(frac)->upVal ;
	    double d = dynamic_cast<MohrCoulomb *>(frac)->downVal ;
	    Vector up(beta.size()) ;
	    Vector down(beta.size()) ;
	    for(size_t i = 0 ; i < beta.size() ; i++)
	    {
		up[i] = u ;
		down[i] = d ;
	    }
	    lambda.push_back(up) ;
	    lambda.push_back(down) ;
	    return ;
	}

    }
}
