#include "viscoelasticity.h"
#include "../elements/generalized_spacetime_viscoelastic_element_state.h"

using namespace Mu ;

void Mu::placeMatrixInBlock(const Matrix & m, size_t i, size_t j, Matrix & ret)
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

void Mu::addMatrixInBlock(const Matrix & m, size_t i, size_t j, Matrix & ret)
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

void Mu::substractMatrixInBlock(const Matrix & m, size_t i, size_t j, Matrix & ret)
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

void Mu::getBlockInMatrix( const Matrix & source, size_t i, size_t j, Matrix & ret)
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


Viscoelasticity::Viscoelasticity(ViscoelasticModel m, const Matrix & rig, int n, double r) : LinearForm(rig, false, false, (1+n)*(rig.numRows()/3+1)), model(m), blocks(1+n)
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
			std::cout << "warning: wrong constructor for Viscoelasticity" << std::endl ;  
	}
} ;

Viscoelasticity::Viscoelasticity(ViscoelasticModel m, const Matrix & rig, const Matrix & e, int b, int n, double r) : LinearForm(rig, false, false, (1+n+b+(m == MAXWELL))*(rig.numRows()/3+1)), model(m), blocks(1+n+b+(m==MAXWELL))
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
			std::cout << "warning: wrong constructor for Viscoelasticity" << std::endl ;  
	}
} ;

Viscoelasticity::Viscoelasticity(const Matrix & rig, const Matrix & e, int b, int n, double r) : LinearForm(rig, false, false, (n+b)*((rig.numRows()/b)/3+1)), model(GENERAL_VISCOELASTICITY), blocks(n+b)
{
	v.push_back(XI);
	v.push_back(ETA);
	if(rig.numRows()/b > 3)
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
			std::cout << "warning: wrong constructor for Viscoelasticity" << std::endl ;  
	}
} ;

Viscoelasticity::Viscoelasticity(ViscoelasticModel m, const Matrix & c_kv, const Matrix & e_kv, const Matrix & c_mx, const Matrix & e_mx, int b, int n, double r) : LinearForm(c_kv, false, false, (3+n+b)*(c_kv.numRows()/3+1)), model(m), blocks(3+b+n)
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
			placeMatrixInBlock( c_mx, 0,0, param) ;
			placeMatrixInBlock( c_mx, 1+b,1+b, param) ;
 			placeMatrixInBlock( c_mx, 2+b,2+b, param) ;
			placeMatrixInBlock( c_mx, 1+b,2+b, param) ;
			placeMatrixInBlock( c_mx, 2+b,1+b, param) ;
			Matrix r_mx = c_mx * (-1.) ;
			placeMatrixInBlock( r_mx, 0+b,1+b, param) ;
			placeMatrixInBlock( r_mx, 0+b,2+b, param) ;
			placeMatrixInBlock( r_mx, 1+b,0+b, param) ;
 			placeMatrixInBlock( r_mx, 2+b,0+b, param) ;
			addMatrixInBlock( c_kv, 2+b,2+b, param) ;
			
			
			placeMatrixInBlock( e_mx, 1+b,1+b, eta) ;
			placeMatrixInBlock( e_kv, 2+b,2+b, eta) ;
			
// 			param.print() ;
// 			eta.print() ;
// 			exit(0) ;
			break ;
		}
		default:
			std::cout << "warning: wrong constructor for Viscoelasticity" << std::endl ;  
	}
} ;

Viscoelasticity::Viscoelasticity(ViscoelasticModel m, const Matrix & c0, std::vector<std::pair<Matrix, Matrix> > & branches, int b, int n, double r) : LinearForm(c0, false, false, (1+n+b+branches.size())*(c0.numRows()/3+1)), model(m), blocks(1+n+b+branches.size())
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
			std::cout << "warning: wrong constructor for Viscoelasticity" << std::endl ;  
	}
} ;

Viscoelasticity::Viscoelasticity(ViscoelasticModel m, const Matrix & c0, const Matrix & c1, const Matrix & e1, int b, int n, double r) : LinearForm(c0, false, false, (2+n+b)*(c0.numRows()/3+1)), model(m), blocks(2+n+b)
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
			std::cout << "warning: wrong constructor for Viscoelasticity" << std::endl ;  
	}
} ;

Viscoelasticity::~Viscoelasticity() {} ;

ElementState * Viscoelasticity::createElementState( IntegrableEntity * e) 
{
	return new GeneralizedSpaceTimeViscoElasticElementState(e) ;  
}

void Viscoelasticity::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
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
			a += b ;
			placeMatrixInBlock( a, 0,0, ret ) ;
			for(size_t i = 1 ; i < blocks ; i++)
			{
				// first line
				substractMatrixInBlock( a, i,0, ret ) ;
				// first column
				substractMatrixInBlock( a, 0,i, ret ) ;
				for(size_t j = i+1 ; j < blocks ; j++)
				{
					// upper triangle
					placeMatrixInBlock( a, i,j, ret ) ;
					// lower triangle
					placeMatrixInBlock( a, j,i, ret ) ;
				}
			}
			for(size_t i = 1 ; i < blocks ; i++)
			{
				//stiffness (diagonal)
				getBlockInMatrix(param, i,i, buffer) ;
				vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
				vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
				a += b ;
				placeMatrixInBlock( a, i,i, ret ) ;
			}
			
			return ;
		}
		
		case GENERALIZED_MAXWELL:
		{
			// stiffness (0,0)
			getBlockInMatrix(param, 0,0, buffer) ;
			vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
			vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
			a += b ;
			placeMatrixInBlock( a, 0,0, ret ) ;
			for(size_t i = 1 ; i < blocks ; i++)
			{
				//stiffness (diagonal)
				getBlockInMatrix(param, i,i, buffer) ;
				vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
				vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
				a += b ;
				placeMatrixInBlock( a, i,i, ret ) ;
				// first line
				substractMatrixInBlock( a, i,0, ret ) ;
				// first column
				substractMatrixInBlock( a, 0,i, ret ) ;				
			}
			return ;
		}
		
		case BURGER:
		{
			// stiffness Maxwell
			getBlockInMatrix(param, 0,0, buffer) ;
			vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
			vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
			a += b ;
			placeMatrixInBlock( a, 0,0, ret ) ;
			placeMatrixInBlock( a, 1,1, ret ) ;
			placeMatrixInBlock( a, 1,2, ret ) ;
			placeMatrixInBlock( a, 2,1, ret ) ;
			substractMatrixInBlock( a, 0,1, ret ) ;
 			substractMatrixInBlock( a, 1,0, ret ) ;
			substractMatrixInBlock( a, 0,2, ret ) ;
 			substractMatrixInBlock( a, 2,0, ret ) ;
			// stiffness KV
			getBlockInMatrix(param, 2,2, buffer) ;
			vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
			vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
			a += b ;
			placeMatrixInBlock( a, 2,2, ret ) ;

			return ;
		}
		
		case MAXWELL:
		{
			// stiffness
			getBlockInMatrix(param, 0,0, buffer) ;
  			vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
  			vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
  			a += b ;
			placeMatrixInBlock( a, 0,0, ret ) ;
			placeMatrixInBlock( a, 1,1, ret ) ;
			substractMatrixInBlock( a, 0,1, ret ) ;
			substractMatrixInBlock( a, 1,0, ret ) ;
			return ;
		}
		
		case KELVIN_VOIGT:
		{
			// stiffness
			getBlockInMatrix(param, 0,0, buffer) ;
			vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
			vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
  			a += b ;
			placeMatrixInBlock( a, 0,0, ret ) ;
			return ;
		}
		
		case PURE_ELASTICITY:
		{
			getBlockInMatrix(param, 0,0, buffer) ;
 			vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
 			vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
 			a += b ;
 			placeMatrixInBlock( a, 0,0, ret ) ;
			
			return ;
		}
		
		case PURE_VISCOSITY:
		{
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
				a += b ;
				placeMatrixInBlock( a, i,i, ret ) ;
				// elastic matrix (upper-triangle)
				for(size_t j = i+1 ; j < blocks ; j++)
				{
					getBlockInMatrix(param, i,j, buffer) ;
					vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
					vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
					a += b ;
					placeMatrixInBlock( a, i,j, ret ) ;
					// symmetry
					placeMatrixInBlock( a, j,i, ret ) ;
				}

			}
		}
	}
	
  
}

void Viscoelasticity::applyViscous(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	Matrix a(ret.numRows()/blocks, ret.numCols()/blocks) ;
	Matrix b(ret.numRows()/blocks, ret.numCols()/blocks) ;

	Matrix buffer(param.numRows()/blocks, param.numCols()/blocks) ;
	
	switch(model)
	{
		case GENERALIZED_KELVIN_VOIGT:
		{
			for(size_t i = 1 ; i < blocks ; i++)
			{
				// viscosity (diagonal)
				getBlockInMatrix(eta, i,i, buffer) ;
				vm->ieval(GradientDot(p_i) * buffer * GradientDot(p_j, true),    gp, Jinv,v, a) ;
				vm->ieval(GradientDotDot(p_i)    * buffer * Gradient(p_j, true), gp, Jinv,v, b) ;
				a += b ;
				placeMatrixInBlock( a, i,i, ret ) ;
			}
			
			return ;
		}
		
		case GENERALIZED_MAXWELL:
		{
			for(size_t i = 1 ; i < blocks ; i++)
			{
				// viscosity (diagonal)
				getBlockInMatrix(eta, i,i, buffer) ;
				vm->ieval(GradientDot(p_i) * buffer * GradientDot(p_j, true),    gp, Jinv,v, a) ;
				vm->ieval(GradientDotDot(p_i)    * buffer * Gradient(p_j, true), gp, Jinv,v, b) ;
				a += b ;
				placeMatrixInBlock( a, i,i, ret ) ;
			}
			return ;
		}
		
		case BURGER:
		{
			// viscosity Maxwell
			getBlockInMatrix(eta, 1,1, buffer) ;
			vm->ieval(GradientDot(p_i)    * buffer   * GradientDot(p_j, true), gp, Jinv,v,a);
// 			vm->ieval(GradientDotDot(p_i) * buffer   * Gradient(p_j, true),    gp, Jinv,v,b);
// 			a += b ;
			placeMatrixInBlock( a, 1,1, ret ) ;
			// viscosity KV
			getBlockInMatrix(eta, 2,2, buffer) ;
			vm->ieval(GradientDot(p_i)    * buffer   * GradientDot(p_j, true), gp, Jinv,v,a);
// 			vm->ieval(GradientDotDot(p_i) * buffer   * Gradient(p_j, true),    gp, Jinv,v,b);
// 			a += b ;
			placeMatrixInBlock( a, 2,2, ret ) ;
			
			
			return ;
		}
		
		case MAXWELL:
		{
			// viscosity
			getBlockInMatrix(eta, 1,1, buffer) ;
 			vm->ieval(GradientDot(p_i)    * buffer   * GradientDot(p_j, true), gp, Jinv,v,a);
 			vm->ieval(GradientDotDot(p_i) * buffer   * Gradient(p_j, true),    gp, Jinv,v,b);
 			a += b ;
			placeMatrixInBlock( a, 1,1, ret ) ;
			return ;
		}
		
		case KELVIN_VOIGT:
		{
			// viscosity
			getBlockInMatrix(eta, 0,0, buffer) ;
			vm->ieval(GradientDot(p_i)    * buffer   * GradientDot(p_j, true), gp, Jinv,v,a);
			vm->ieval(GradientDotDot(p_i) * buffer   * Gradient(p_j, true),    gp, Jinv,v,b);
			a += b ;
			placeMatrixInBlock( a, 0,0, ret ) ;
			return ;
		}
		
		case PURE_ELASTICITY:
		{
			return ;
		}
		
		case PURE_VISCOSITY:
		{
			getBlockInMatrix(eta, 0,0, buffer) ;
 			vm->ieval(GradientDot(p_i)    * buffer   * GradientDot(p_j, true), gp, Jinv,v,a);
  			vm->ieval(GradientDotDot(p_i) * buffer   * Gradient(p_j, true),    gp, Jinv,v,b);
			a+= b ;
			placeMatrixInBlock( a, 0,0, ret ) ;
			return ;
		}
		
		default:
		{
			for(size_t i = 0 ; i < blocks ; i++)
			{
				// viscous matrix (diagonal)
				getBlockInMatrix(eta, i,i, buffer) ;
				vm->ieval(GradientDot(p_i)    * buffer   * GradientDot(p_j, true), gp, Jinv,v,a);
				vm->ieval(GradientDotDot(p_i) * buffer   * Gradient(p_j, true),    gp, Jinv,v,b);
				a += b ;
				placeMatrixInBlock( a, i,i, ret ) ;
				// viscous matrix (upper-triangle)
				for(size_t j = i+1 ; j < blocks ; j++)
				{
					getBlockInMatrix(eta, i,j, buffer) ;
					vm->ieval(GradientDot(p_i)    * buffer   * GradientDot(p_j, true), gp, Jinv,v,a);
					vm->ieval(GradientDotDot(p_i) * buffer   * Gradient(p_j, true),    gp, Jinv,v,b);
					a += b ;
					addMatrixInBlock( a, i,j, ret ) ;
					// symmetry
					addMatrixInBlock( a, j,i, ret ) ;
				}
			}
		}
	}
	
  
}

bool Viscoelasticity::fractured() const
{
	return false ;
}

bool Viscoelasticity::changed() const
{
	return false ;
} 

Form * Viscoelasticity::getCopy() const 
{
	Matrix rig = param ;
	Matrix e = eta ;
	Viscoelasticity * copy = new Viscoelasticity(rig, e, blocks) ;
	copy->model = model ;
	
	if(getExtra2dMeshes())
	{
		for(size_t i = 0 ; i < getExtra2dMeshes()->size() ; i++)
			copy->addMesh((*getExtra2dMeshes())[i]);
	}
	if(getExtra3dMeshes())
	{
		for(size_t i = 0 ; i < getExtra3dMeshes()->size() ; i++)
			copy->addMesh((*getExtra3dMeshes())[i]);
	}
	return copy ; 
}

Vector Viscoelasticity::getForcesFromAppliedStress( Vector & data, Function & shape, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<Variable> & v, bool isVolumic) 
{
 	if(isVolumic)
	{
		return (VirtualMachine().ieval(GradientDot( shape ) * data, gp, Jinv, v)) ;
	}
	return data * VirtualMachine().ieval(Differential( shape, TIME_VARIABLE ), gp, Jinv, v) ;
}

Vector Viscoelasticity::getForcesFromAppliedStress( const Function & data, size_t index, size_t externaldofs,  Function & shape, IntegrableEntity * e,const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<Variable> & v, bool isVolumic) 
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

void Viscoelasticity::print() const
{
	std::cout << "I am a viscoelastic model" ;
	switch(model)
	{
	  case PURE_ELASTICITY:
	    std::cout << " (elastic only)" << std::endl ;
	    return ;
	  case PURE_VISCOSITY:
	    std::cout << " (viscous only)" << std::endl ;
	    return ;
	  case KELVIN_VOIGT:
	    std::cout << " (kelvin-voigt)" << std::endl ;
	    return ;
	  case GENERALIZED_KELVIN_VOIGT:
	    std::cout << " ( generalized kelvin-voigt)" << std::endl ;
	    return ;
	  case MAXWELL:
	    std::cout << " (maxwell)" << std::endl ;
	    return ;
	  case GENERALIZED_MAXWELL:
	    std::cout << " ( generalized maxwell)" << std::endl ;
	    return ;
	  case BURGER:
	    std::cout << " (burger)" << std::endl ;
	    return ;
	  case GENERAL_VISCOELASTICITY:
	    std::cout << " (general)" << std::endl ;
	    return ;
	}
}

