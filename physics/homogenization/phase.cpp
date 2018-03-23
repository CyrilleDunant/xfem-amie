#include "phase.h"
#include "../stiffness_with_imposed_deformation.h"
#include "../stiffness.h"
#include "../homogeneised_behaviour.h"
#include "../dual_behaviour.h"
#include "../../utilities/matrixops.h"
#include "../fracturecriteria/mohrcoulomb.h"

using namespace Amie ;

Phase::Phase(InclusionGeometryType t, double a, double b, double c): t(t), a(a), b(b), c(c)
{
	volume = 1e-6 ;
	C = Matrix( 6, 6 ) ;
	beta = Vector( 3 ) ;
	A = Matrix( 6, 6 ) ;
}

Phase::Phase( Form * f, double v, SpaceDimensionality dim,InclusionGeometryType t, double a, double b, double c): t(t), a(a), b(b), c(c)
{
	behaviour = f ;
	volume = v ;
	stiffnessFromBehaviour(dim) ;
	expansionFromBehaviour(dim) ;
	ruptureFromBehaviour(dim) ;
    volume = std::max(volume, 1e-6) ;
	A = Matrix( C.numRows(), C.numCols() ) ;
}	

Phase::Phase( DelaunayTriangle *tri ,InclusionGeometryType t, double a, double b, double c): t(t), a(a), b(b), c(c)
{
	behaviour = tri->getBehaviour() ;
	volume = tri->area() ;
	if(dynamic_cast<BimaterialInterface *>(behaviour))
	{
		int count = 0;
		int out = 0 ;
		for(double i = 0.00 ; i < 1. ; i+=0.01)
		{
			for(double j = 0.00 ; j < 1.-i ; j += 0.01)
			{
				count++ ;
				out += (int) !dynamic_cast<BimaterialInterface *>(behaviour)->insideGeometry(Point(i,j)) ;
			}
		}
		volume *= out/count ;
	}
	if(dynamic_cast<HomogeneisedBehaviour *>(behaviour))
		behaviour = dynamic_cast<HomogeneisedBehaviour *>(behaviour)->getOriginalBehaviour() ;
	stiffnessFromBehaviour() ;
	expansionFromBehaviour() ;
	ruptureFromBehaviour() ;
    volume = std::max(volume, 1e-6) ;
	A = Matrix( C.numRows(), C.numCols() ) ;
}

Phase::Phase( DelaunayTetrahedron *tet ,InclusionGeometryType t, double a, double b, double c): t(t), a(a), b(b), c(c)
{
	behaviour = tet->getBehaviour() ;
	volume = tet->volume() ;
	stiffnessFromBehaviour() ;
	expansionFromBehaviour() ;
	ruptureFromBehaviour() ;
    volume = std::max(volume, 1e-6) ;
	A = Matrix( C.numRows(), C.numCols() ) ;
}

Phase::Phase( Feature *f ,InclusionGeometryType t, double a, double b, double c): t(t), a(a), b(b), c(c)
{
	behaviour = f->getBehaviour() ;

	if( f->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
		volume = f->area() ;

	if( f->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
		volume = f->volume() ;

	stiffnessFromBehaviour() ;
	expansionFromBehaviour() ;
	ruptureFromBehaviour() ;
     volume = std::max(volume, 1e-6) ;
	A = Matrix( C.numRows(), C.numCols() ) ;
}

Phase::Phase( Feature *f , DelaunayTriangle * tri,InclusionGeometryType t, double a, double b, double c): t(t), a(a), b(b), c(c)
{
	behaviour = f->getBehaviour() ;
	
	if( f->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
	{
		volume = f->area() ;
		if(dynamic_cast<Circle *>(f))
		{
			Circle overlap(f->getRadius(), f->getCenter()) ;
			overlap.sampleSurface(50./f->getRadius(), 1.) ;
			int count = 0 ;
			for(size_t i = 0 ; i < overlap.getBoundingPoints().size() ; i++)
				count += (int) tri->in(overlap.getBoundingPoint(i)) ;
			for(size_t i = 0 ; i < overlap.getInPoints().size() ; i++)
				count += (int) tri->in(overlap.getInPoint(i)) ;
			volume *= count ;
			volume /= overlap.getBoundingPoints().size() + overlap.getInPoints().size() ;

			
		}
		
	}
	
	if( f->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
		volume = f->volume() ;
    volume = std::max(volume, 1e-6) ;
	
	stiffnessFromBehaviour() ;
	expansionFromBehaviour() ;
	ruptureFromBehaviour() ;
	A = Matrix( C.numRows(), C.numCols() ) ;
}

Phase::Phase( const Phase &p ): t(p.t), a(p.a), b(p.b), c(p.c)
{
	C.resize( p.C.numRows(), p.C.numCols() ) ;
	C = p.C ;
	beta.resize( p.beta.size() ) ;
	beta = p.beta ;
	A.resize( p.A.numRows(), p.A.numCols() ) ;
	A = p.A ;
	volume = p.volume ;
    volume = std::max(volume, 1e-6) ;
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

Form * Phase::getBehaviour()
{
	apply() ;
	Matrix S = C ;
        
//         S.print() ;
	if( S.size() == 36 )
		invert6x6Matrix( S ) ;
	else
		invert3x3Matrix( S ) ;
	
	Vector alpha = S * beta ;
//         std::cout << alpha[0] << std::endl ;
//         S.print() ;
	return new StiffnessWithImposedDeformation( C, alpha ) ;
}

Phase &Phase::operator =( const Phase &p )
{
    t = p.t ;
    a = p.a ;
    b = p.b ;
    c = p.c ;
	C.resize( p.C.numRows(), p.C.numCols() ) ;
	C = p.C ;
	beta.resize( p.beta.size() ) ;
	beta = p.beta ;
	A.resize( p.A.numRows(), p.A.numCols() ) ;
	A = p.A ;
	volume = p.volume ;
    volume = std::max(volume, 1e-6) ;
	for(size_t i = 0 ; i < p.lambda.size() ; i++)
	{
	    Vector l(p.lambda[i].size()) ;
	    l = p.lambda[i] ;
	    lambda.push_back(l);
	}
	return *this ;
}

void Phase::stiffnessFromBehaviour(SpaceDimensionality dim)
{
	if(behaviour->type == VOID_BEHAVIOUR)
	{
		C.resize( 3+3*(dim == SPACE_THREE_DIMENSIONAL), 3+3*(dim == SPACE_THREE_DIMENSIONAL)) ;
		C = 0. ;
		return ;
	}

	Matrix tmp =  behaviour->getTensor( Point(0.,0.) ) ;
	if(dynamic_cast<BimaterialInterface *>(behaviour))
	{
		tmp = Matrix(tmp.numRows(), tmp.numCols()) ;
		int count = 0;
		for(double i = 0.00 ; i < 1. ; i+=0.01)
		{
			for(double j = 0.00 ; j < 1.-i ; j += 0.01)
			{
				count++ ;
				tmp += behaviour->getTensor(Point(i,j)) ;
			}
		}
		tmp /= (count) ;
	}
	if(dynamic_cast<HomogeneisedBehaviour *>(behaviour))
	{
		tmp =  dynamic_cast<HomogeneisedBehaviour *>(behaviour)->original->getTensor( Point(0.,0.) ) ;
	}
	if(dynamic_cast<HomogeneisedBehaviour *>(behaviour))
	{
		tmp =  dynamic_cast<HomogeneisedBehaviour *>(behaviour)->original->getTensor( Point(0.,0.) ) ;
	}
	C.resize( tmp.numRows(), tmp.numCols() ) ;
	C = tmp ;
		
}

void Phase::expansionFromBehaviour(SpaceDimensionality dim)
{
	if(behaviour->type == VOID_BEHAVIOUR)
	{
		beta.resize( 3+3*(dim == SPACE_THREE_DIMENSIONAL)) ;
		beta = 0. ;
		return ;
	}

	Vector tmp = behaviour->getImposedStrain( Point( 1. / 3, 1. / 3, 1. / 3 ) ) ;
	beta.resize( tmp.size() );
	beta = tmp * C;
}

void Phase::ruptureFromBehaviour(SpaceDimensionality dim)
{
    lambda.clear() ;
    if(behaviour->type == VOID_BEHAVIOUR)
    {
	return ;
    }
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

void Phase::print()
{
	C.print() ;
	std::cout << std::endl ;
	for(size_t i = 0 ; i < beta.size() ; i++)
		std::cout << beta[i] << std::endl ;
}
