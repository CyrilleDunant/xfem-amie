#include "generalized_spacetime_viscoelasticity.h"
#include "../elements/generalized_spacetime_viscoelastic_element_state.h"

using namespace Mu ;

void placeMatrixInBlock(const Matrix & m, size_t i, size_t j, Matrix & ret)
{
	size_t r = m.numRows() ;
	size_t c = m.numCols() ;
	for(size_t k = i*r ; k < (i+1)*r ; k++)
	{
		for(size_t l = j*c ; l < (j+1)*c ; l++)
		{
			ret[k][l] = m[k-i*r][l-j*c] ;
		}
	}
}

void addMatrixInBlock(const Matrix & m, size_t i, size_t j, Matrix & ret)
{
	size_t r = m.numRows() ;
	size_t c = m.numCols() ;
	for(size_t k = i*r ; k < (i+1)*r ; k++)
	{
		for(size_t l = j*c ; l < (j+1)*c ; l++)
		{
			ret[k][l] += m[k-i*r][l-j*c] ;
		}
	}
}

void substractMatrixInBlock(const Matrix & m, size_t i, size_t j, Matrix & ret)
{
	size_t r = m.numRows() ;
	size_t c = m.numCols() ;
	for(size_t k = i*r ; k < (i+1)*r ; k++)
	{
		for(size_t l = j*c ; l < (j+1)*c ; l++)
		{
			ret[k][l] -= m[k-i*r][l-j*c] ;
		}
	}
}

void getBlockInMatrix( const Matrix & source, size_t i, size_t j, Matrix & ret)
{
	size_t r = ret.numRows() ;
	size_t c = ret.numCols() ;
	for(size_t k = i*r ; k < (i+1)*r ; k++)
	{
		for(size_t l = j*c ; l < (j+1)*c ; l++)
		{
			ret[k-i*r][l-j*c] = source[k][l] ;
		}
	}	
}


GeneralizedSpaceTimeViscoelasticity::GeneralizedSpaceTimeViscoelasticity(ViscoelasticModel m, const Matrix & rig, int n) : LinearForm(rig, false, false, (1+n)*(rig.numRows()/3+1)), model(m), blocks(1+n)
{
	v.push_back(XI);
	v.push_back(ETA);
	if(rig.size() > 9)
		v.push_back(ZETA);
	v.push_back(TIME_VARIABLE);
	
	param.resize(rig.numRows()*(1+n), rig.numCols()*(1+n)) ;
	eta.resize(rig.numRows()*(1+n), rig.numCols()*(1+n)) ;
	param.array() = 0 ;
	eta.array() = 0 ;
	
	switch(model)
	{
		case PURE_ELASTICITY:
			placeMatrixInBlock( rig, 0,0, param) ;
			break ;
		case PURE_VISCOSITY:
			placeMatrixInBlock( rig, 0,0, eta) ;
			break ;
		default:
			std::cout << "warning: wrong constructor for GeneralizedSpaceTimeViscoelasticity" << std::endl ;  
	}
} ;

GeneralizedSpaceTimeViscoelasticity::GeneralizedSpaceTimeViscoelasticity(ViscoelasticModel m, const Matrix & rig, const Matrix & e, int b, int n) : LinearForm(rig, false, false, (1+n+b+(m == MAXWELL))*(rig.numRows()/3+1)), model(m), blocks(1+n+b+(m==MAXWELL))
{ 
	v.push_back(XI);
	v.push_back(ETA);
	if(rig.size() > 9)
		v.push_back(ZETA);
	v.push_back(TIME_VARIABLE);
	
	param.resize(rig.numRows()*(1+n+b+(m == MAXWELL)), rig.numCols()*(1+n+b+(m == MAXWELL))) ;
	eta.resize(rig.numRows()*(1+n+b+(m == MAXWELL)), rig.numCols()*(1+n+b+(m == MAXWELL))) ;
	param.array() = 0 ;
	eta.array() = 0 ;
	
	switch(model)
	{
		case KELVIN_VOIGT:
			placeMatrixInBlock( rig, 0,0, param) ;
			placeMatrixInBlock( e, 0,0, eta) ;
			break ;
		case MAXWELL:
		{
			Matrix r0 = rig*(-1.) ;
			placeMatrixInBlock( rig, 0,0, param) ;
			placeMatrixInBlock( r0,  0,1+b, param) ;
			placeMatrixInBlock( r0,  1+b,0, param) ;
			placeMatrixInBlock( rig, 1+b,1+b, param) ;
			placeMatrixInBlock( e, 1+b,1+b, eta) ;
			break ;
		}
		default:
			std::cout << "warning: wrong constructor for GeneralizedSpaceTimeViscoelasticity" << std::endl ;  
	}
} ;

GeneralizedSpaceTimeViscoelasticity::GeneralizedSpaceTimeViscoelasticity(const Matrix & rig, const Matrix & e, int b, int n) : LinearForm(rig, false, false, (n+b)*((rig.numRows()/b)/3+1)), model(GENERAL_VISCOELASTICITY), blocks(n+b)
{
	v.push_back(XI);
	v.push_back(ETA);
	if(rig.size() > 9)
		v.push_back(ZETA);
	v.push_back(TIME_VARIABLE);
	
	param.resize(rig.numRows()/b*blocks, rig.numCols()/b*blocks) ;
	eta.resize(rig.numRows()/b*blocks, rig.numCols()/b*blocks) ;
	param.array() = 0 ;
	eta.array() = 0 ;
	
	switch(model)
	{
		case GENERAL_VISCOELASTICITY:
			placeMatrixInBlock( rig, 0,0, param) ;
			placeMatrixInBlock( e, 0,0, eta) ;
			break ;
		default:
			std::cout << "warning: wrong constructor for GeneralizedSpaceTimeViscoelasticity" << std::endl ;  
	}
} ;

GeneralizedSpaceTimeViscoelasticity::GeneralizedSpaceTimeViscoelasticity(ViscoelasticModel m, const Matrix & c_kv, const Matrix & e_kv, const Matrix & c_mx, const Matrix & e_mx, int b, int n) : LinearForm(c_kv, false, false, (3+n+b)*(c_kv.numRows()/3+1)), model(m), blocks(3+b+n)
{
	v.push_back(XI);
	v.push_back(ETA);
	if(c_kv.size() > 9)
		v.push_back(ZETA);
	v.push_back(TIME_VARIABLE);
	
	param.resize(c_kv.numRows()*(3+n+b), c_kv.numCols()*(3+n+b)) ;
	eta.resize(c_kv.numRows()*(3+n+b), c_kv.numCols()*(3+n+b)) ;
	param.array() = 0 ;
	eta.array() = 0 ;
	
	switch(model)
	{
		case BURGER:
		{
			placeMatrixInBlock( c_kv, 0,0, param) ;
			placeMatrixInBlock( c_kv, 1+b,1+b, param) ;
			addMatrixInBlock( c_mx, 1+b,1+b, param) ;
			placeMatrixInBlock( c_mx, 2+b,2+b, param) ;
			Matrix r_kv = c_kv*(-1) ;
			placeMatrixInBlock( r_kv, 0,1+b, param) ;
			placeMatrixInBlock( r_kv, 1+b,0, param) ;
			Matrix r_mx = c_mx*(-1) ;
			placeMatrixInBlock( r_mx, 1+b,2+b, param) ;
			placeMatrixInBlock( r_mx, 2+b,1+b, param) ;
			
			placeMatrixInBlock( e_kv, 0,0, eta) ;
			placeMatrixInBlock( e_kv, 1+b,1+b, eta) ;
			placeMatrixInBlock( e_mx, 2+b,2+b, eta) ;
			Matrix v_kv = e_kv*(-1) ;
			placeMatrixInBlock( v_kv, 0,1+b, eta) ;
			placeMatrixInBlock( v_kv, 1+b,0, eta) ;
			break ;
		}
		default:
			std::cout << "warning: wrong constructor for GeneralizedSpaceTimeViscoelasticity" << std::endl ;  
	}
} ;

GeneralizedSpaceTimeViscoelasticity::GeneralizedSpaceTimeViscoelasticity(ViscoelasticModel m, const Matrix & c0, std::vector<std::pair<Matrix, Matrix> > & branches, int b, int n) : LinearForm(c0, false, false, (1+n+b+branches.size())*(c0.numRows()/3+1)), model(m), blocks(1+n+b+branches.size())
{
	v.push_back(XI);
	v.push_back(ETA);
	if(c0.size() > 9)
		v.push_back(ZETA);
	v.push_back(TIME_VARIABLE);
	
	param.resize(c0.numRows()*(1+n+b+branches.size()), c0.numCols()*(1+n+b+branches.size())) ;
	eta.resize(c0.numRows()*(1+n+b+branches.size()), c0.numCols()*(1+n+b+branches.size())) ;
	param.array() = 0 ;
	eta.array() = 0 ;
	
	switch(model)
	{
		case GENERALIZED_KELVIN_VOIGT:
		{
			placeMatrixInBlock( c0, 0,0, param) ;
			Matrix r0 = c0*(-1) ;
			for(size_t i = 1 ; i < branches.size()+1 ; i++)
			{
				placeMatrixInBlock( r0, i+b,0, param) ;
				placeMatrixInBlock( r0, 0,i+b, param) ;
				placeMatrixInBlock( c0, i+b,i+b, param) ;
				addMatrixInBlock( branches[i-1].first, i+b,i+b, param) ;
				for(size_t j = i+1 ; j < branches.size()+1 ; j++)
				{
					placeMatrixInBlock( c0, i+b,j+b, param) ;
					placeMatrixInBlock( c0, j+b,i+b, param) ;
				}
				placeMatrixInBlock( branches[i-1].second, i+b,i+b, eta) ;
			}
			break ;
		}
		case GENERALIZED_MAXWELL:
		{
			placeMatrixInBlock( c0, 0,0, param) ;
			for(size_t i = 1 ; i < branches.size()+1 ; i++)
			{
				Matrix ri = branches[i-1].first * (-1) ;
				addMatrixInBlock( branches[i-1].first, 0,0, param) ;
				placeMatrixInBlock( branches[i-1].first, i+b,i+b, param) ;
				placeMatrixInBlock( ri, 0,i+b, param) ;
				placeMatrixInBlock( ri, i+b,0, param) ;
				placeMatrixInBlock( branches[i-1].second, i+b,i+b, eta) ;
			}
			break ;
		}
		default:
			std::cout << "warning: wrong constructor for GeneralizedSpaceTimeViscoelasticity" << std::endl ;  
	}
} ;

GeneralizedSpaceTimeViscoelasticity::GeneralizedSpaceTimeViscoelasticity(ViscoelasticModel m, const Matrix & c0, const Matrix & c1, const Matrix & e1, int b, int n) : LinearForm(c0, false, false, (2+n+b)*(c0.numRows()/3+1)), model(m), blocks(2+n+b)
{
	v.push_back(XI);
	v.push_back(ETA);
	if(c0.size() > 9)
		v.push_back(ZETA);
	v.push_back(TIME_VARIABLE);
	
	param.resize(c0.numRows()*(2+n+b), c0.numCols()*(2+n+b)) ;
	eta.resize(c0.numRows()*(2+n+b), c0.numCols()*(2+n+b)) ;
	param.array() = 0 ;
	eta.array() = 0 ;
	
	switch(model)
	{
		case GENERALIZED_KELVIN_VOIGT:
		{
			placeMatrixInBlock( c0, 0,0, param) ;
			placeMatrixInBlock( c0, 1+b,1+b, param) ;
			addMatrixInBlock( c1, 1+b,1+b, param) ;
			Matrix r0 = c0*(-1) ;
			placeMatrixInBlock( r0, 1+b,0, param) ;
			placeMatrixInBlock( r0, 0,1+b, param) ;
			addMatrixInBlock( e1, 1+b,1+b, eta) ;
			break ;
		}
		case GENERALIZED_MAXWELL:
		{
			placeMatrixInBlock( c0, 0,0, param) ;
			addMatrixInBlock( c1, 0,0, param) ;
			placeMatrixInBlock( c1, 1+b,1+b, param) ;
			Matrix r1 = c1*(-1) ;
			placeMatrixInBlock( r1, 1+b,0, param) ;
			placeMatrixInBlock( r1, 0,1+b, param) ;
			placeMatrixInBlock( e1, 1+b,1+b, eta) ;
			break ;
		}
		default:
			std::cout << "warning: wrong constructor for GeneralizedSpaceTimeViscoelasticity" << std::endl ;  
	}
} ;

GeneralizedSpaceTimeViscoelasticity::~GeneralizedSpaceTimeViscoelasticity() { } ;

ElementState * GeneralizedSpaceTimeViscoelasticity::createElementState( IntegrableEntity * e) 
{
	return new GeneralizedSpaceTimeViscoElasticElementState(e) ;  
}

void GeneralizedSpaceTimeViscoelasticity::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	Matrix a(ret.numRows()/blocks, ret.numCols()/blocks) ;
	Matrix b(ret.numRows()/blocks, ret.numCols()/blocks) ;
	
	Matrix buffer(param.numRows()/blocks, param.numCols()/blocks) ;

	switch(model)
	{
		case GENERALIZED_KELVIN_VOIGT:
		{
			// stiffness (0,0)
			getBlockInMatrix(param, 0,0, buffer) ;
			vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
			vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
			placeMatrixInBlock( a, 0,0, ret ) ;
			addMatrixInBlock( b, 0,0, ret ) ;
			for(size_t i = 1 ; i < blocks ; i++)
			{
				// first line
				substractMatrixInBlock( a, i,0, ret ) ;
				substractMatrixInBlock( b, i,0, ret ) ;
				// first column
				substractMatrixInBlock( a, 0,i, ret ) ;
				substractMatrixInBlock( b, 0,i, ret ) ;
				for(size_t j = i+1 ; j < blocks ; j++)
				{
					// upper triangle
					placeMatrixInBlock( a, i,j, ret ) ;
					addMatrixInBlock( b, i,j, ret ) ;
					// lower triangle
					placeMatrixInBlock( a, j,i, ret ) ;
					addMatrixInBlock( b, j,i, ret ) ;
				}
			}
			for(size_t i = 1 ; i < blocks ; i++)
			{
				//stiffness (diagonal)
				getBlockInMatrix(param, i,i, buffer) ;
				vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
				vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
				placeMatrixInBlock( a, i,i, ret ) ;
				addMatrixInBlock( b, i,i, ret ) ;
				// viscosity (diagonal)
				getBlockInMatrix(eta, i,i, buffer) ;
				vm->ieval(GradientDot(p_i) * buffer * GradientDot(p_j, true),    gp, Jinv,v, a) ;
				vm->ieval(GradientDotDot(p_i)    * buffer * Gradient(p_j, true), gp, Jinv,v, b) ;
				addMatrixInBlock( a, i,i, ret ) ;
				addMatrixInBlock( b, i,i, ret ) ;
			}
			return ;
		}
		
		case GENERALIZED_MAXWELL:
		{
			// stiffness (0,0)
			getBlockInMatrix(param, 0,0, buffer) ;
			vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
			vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
			placeMatrixInBlock( a, 0,0, ret ) ;
			addMatrixInBlock( b, 0,0, ret ) ;
			for(size_t i = 1 ; i < blocks ; i++)
			{
				//stiffness (diagonal)
				getBlockInMatrix(param, i,i, buffer) ;
				vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
				vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
				placeMatrixInBlock( a, i,i, ret ) ;
				addMatrixInBlock( b, i,i, ret ) ;
				// first line
				substractMatrixInBlock( a, i,0, ret ) ;
				substractMatrixInBlock( b, i,0, ret ) ;
				// first column
				substractMatrixInBlock( a, 0,i, ret ) ;
				substractMatrixInBlock( b, 0,i, ret ) ;
				
				// viscosity (diagonal)
				getBlockInMatrix(eta, i,i, buffer) ;
				vm->ieval(GradientDot(p_i) * buffer * GradientDot(p_j, true),    gp, Jinv,v, a) ;
				vm->ieval(GradientDotDot(p_i)    * buffer * Gradient(p_j, true), gp, Jinv,v, b) ;
				addMatrixInBlock( a, i,i, ret ) ;
				addMatrixInBlock( b, i,i, ret ) ;
			}
			return ;
		}
		
		case BURGER:
		{
			// stiffness KV
			getBlockInMatrix(param, 0,0, buffer) ;
			vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
			vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
			placeMatrixInBlock( a, 0,0, ret ) ;
			addMatrixInBlock( b, 0,0, ret ) ;
			placeMatrixInBlock( a, 1,1, ret ) ;
			addMatrixInBlock( b, 1,1, ret ) ;
			substractMatrixInBlock( a, 0,1, ret ) ;
			substractMatrixInBlock( b, 0,1, ret ) ;
			substractMatrixInBlock( a, 1,0, ret ) ;
			substractMatrixInBlock( b, 1,0, ret ) ;
			// stiffness Maxwell
			getBlockInMatrix(param, 2,2, buffer) ;
			vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
			vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
			addMatrixInBlock( a, 1,1, ret ) ;
			addMatrixInBlock( b, 1,1, ret ) ;
			placeMatrixInBlock( a, 2,2, ret ) ;
			addMatrixInBlock( b, 2,2, ret ) ;
			substractMatrixInBlock( a, 1,2, ret ) ;
			substractMatrixInBlock( b, 1,2, ret ) ;
			substractMatrixInBlock( a, 2,1, ret ) ;
			substractMatrixInBlock( b, 2,1, ret ) ;
			// viscosity KV
			getBlockInMatrix(eta, 0,0, buffer) ;
			vm->ieval(GradientDot(p_i)    * buffer   * GradientDot(p_j, true), gp, Jinv,v,a);
			vm->ieval(GradientDotDot(p_i) * buffer   * Gradient(p_j, true),    gp, Jinv,v,b);
			addMatrixInBlock( a, 0,0, ret ) ;
			addMatrixInBlock( b, 0,0, ret ) ;
			addMatrixInBlock( a, 1,1, ret ) ;
			addMatrixInBlock( b, 1,1, ret ) ;
			substractMatrixInBlock( a, 0,1, ret ) ;
			substractMatrixInBlock( b, 0,1, ret ) ;
			substractMatrixInBlock( a, 1,0, ret ) ;
			substractMatrixInBlock( b, 1,0, ret ) ;
			// viscosity Maxwell
			getBlockInMatrix(eta, 2,2, buffer) ;
			vm->ieval(GradientDot(p_i)    * buffer   * GradientDot(p_j, true), gp, Jinv,v,a);
			vm->ieval(GradientDotDot(p_i) * buffer   * Gradient(p_j, true),    gp, Jinv,v,b);
			addMatrixInBlock( a, 2,2, ret ) ;
			addMatrixInBlock( b, 2,2, ret ) ;
			return ;
		}
		
		case MAXWELL:
		{
			// stiffness
			getBlockInMatrix(param, 0,0, buffer) ;
			vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
			vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
			placeMatrixInBlock( a, 0,0, ret ) ;
			addMatrixInBlock( b, 0,0, ret ) ;
			placeMatrixInBlock( a, 1,1, ret ) ;
			addMatrixInBlock( b, 1,1, ret ) ;
			substractMatrixInBlock( a, 0,1, ret ) ;
			substractMatrixInBlock( b, 0,1, ret ) ;
			substractMatrixInBlock( a, 1,0, ret ) ;
			substractMatrixInBlock( b, 1,0, ret ) ;
			// viscosity
			getBlockInMatrix(eta, 1,1, buffer) ;
			vm->ieval(GradientDot(p_i)    * buffer   * GradientDot(p_j, true), gp, Jinv,v,a);
			vm->ieval(GradientDotDot(p_i) * buffer   * Gradient(p_j, true),    gp, Jinv,v,b);
			addMatrixInBlock( a, 1,1, ret ) ;
			addMatrixInBlock( b, 1,1, ret ) ;
			return ;
		}
		
		case KELVIN_VOIGT:
		{
			// stiffness
			getBlockInMatrix(param, 0,0, buffer) ;
			vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
			vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
			placeMatrixInBlock( a, 0,0, ret ) ;
			addMatrixInBlock( b, 0,0, ret ) ;
			// viscosity
			getBlockInMatrix(eta, 0,0, buffer) ;
			vm->ieval(GradientDot(p_i)    * buffer   * GradientDot(p_j, true), gp, Jinv,v,a);
			vm->ieval(GradientDotDot(p_i) * buffer   * Gradient(p_j, true),    gp, Jinv,v,b);
			addMatrixInBlock( a, 0,0, ret ) ;
			addMatrixInBlock( b, 0,0, ret ) ;
			return ;
		}
		
		case PURE_ELASTICITY:
		{
			getBlockInMatrix(param, 0,0, buffer) ;
			vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
			vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
			placeMatrixInBlock( a, 0,0, ret ) ;
			addMatrixInBlock( b, 0,0, ret ) ;
			
			
			return ;
		}
		
		case PURE_VISCOSITY:
		{
			getBlockInMatrix(eta, 0,0, buffer) ;
			vm->ieval(GradientDot(p_i)    * buffer   * GradientDot(p_j, true), gp, Jinv,v,a);
			vm->ieval(GradientDotDot(p_i) * buffer   * Gradient(p_j, true),    gp, Jinv,v,b);
			placeMatrixInBlock( a, 0,0, ret ) ;
			addMatrixInBlock( b, 0,0, ret ) ;
			
			return ;
		}
		
		default:
		{
			for(size_t i = 0 ; i < blocks ; i++)
			{
				// elastic matrix (diagonal)
				getBlockInMatrix(param, i,i, buffer) ;
				vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
				vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
				placeMatrixInBlock( a, i,i, ret ) ;
				addMatrixInBlock( b, i,i, ret ) ;
				// elastic matrix (upper-triangle)
				for(size_t j = i+1 ; j < blocks ; j++)
				{
					getBlockInMatrix(param, i,j, buffer) ;
					vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
					vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
					placeMatrixInBlock( a, i,j, ret ) ;
					addMatrixInBlock( b, i,j, ret ) ;
					// symmetry
					placeMatrixInBlock( a, j,i, ret ) ;
					addMatrixInBlock( b, j,i, ret ) ;
				}

				// viscous matrix (diagonal)
				getBlockInMatrix(eta, i,i, buffer) ;
				vm->ieval(GradientDot(p_i)    * buffer   * GradientDot(p_j, true), gp, Jinv,v,a);
				vm->ieval(GradientDotDot(p_i) * buffer   * Gradient(p_j, true),    gp, Jinv,v,b);
				addMatrixInBlock( a, i,i, ret ) ;
				addMatrixInBlock( b, i,i, ret ) ;
				// viscous matrix (upper-triangle)
				for(size_t j = i+1 ; j < blocks ; j++)
				{
					getBlockInMatrix(eta, i,j, buffer) ;
					vm->ieval(GradientDot(p_i)    * buffer   * GradientDot(p_j, true), gp, Jinv,v,a);
					vm->ieval(GradientDotDot(p_i) * buffer   * Gradient(p_j, true),    gp, Jinv,v,b);
					addMatrixInBlock( a, i,j, ret ) ;
					addMatrixInBlock( b, i,j, ret ) ;
					// symmetry
					addMatrixInBlock( a, j,i, ret ) ;
					addMatrixInBlock( b, j,i, ret ) ;
				}
			}
		}
	}
	
  
}

bool GeneralizedSpaceTimeViscoelasticity::fractured() const
{
	return false ;
}

bool GeneralizedSpaceTimeViscoelasticity::changed() const
{
	return false ;
} 

Form * GeneralizedSpaceTimeViscoelasticity::getCopy() const 
{
	return new GeneralizedSpaceTimeViscoelasticity(*this) ;
}

Vector GeneralizedSpaceTimeViscoelasticity::getForcesFromAppliedStress( Vector & data, Function & shape, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<Variable> & v) 
{
	return VirtualMachine().ieval(GradientDot( shape ) * ( data ), gp, Jinv, v) ;
}

Vector GeneralizedSpaceTimeViscoelasticity::getForcesFromAppliedStress( const Function & data, size_t index, size_t externaldofs,  Function & shape, IntegrableEntity * e,const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<Variable> & v) 
{
	VirtualMachine vm ;
	
	size_t n = e->getBoundingPoints().size() ;
	Vector field(0., n*externaldofs) ;
	for(size_t i = 0 ; i < n ; i++)
		field[ i*externaldofs + index ] = vm.eval( data, e->getBoundingPoint(i) ) ;
	
	std::vector<Vector> g(e->getGaussPoints().gaussPoints.size(), Vector(0., externaldofs)) ;
	e->getState().getExternalFieldAtGaussPoints( field, externaldofs, g) ;
	
	Vector f = vm.ieval( GradientDot( shape ) * g, e, v) ;

	field = 0 ;
	for(size_t i = 0 ; i < n ; i++)
		field[ i*externaldofs + index ] = vm.deval( data, TIME_VARIABLE, e->getBoundingPoint(i) ) ;

	e->getState().getExternalFieldAtGaussPoints( field, externaldofs, g) ;

	f += vm.ieval( Gradient( shape ) * g, e, v) ;
	
	return f ;
}