#include "viscoelasticity_and_fracture.h"
#include "viscoelasticity.h"
#include "fracturecriteria/maxstrain.h"
#include "../elements/generalized_spacetime_viscoelastic_element_state.h"

using namespace Mu ;

/*void placeMatrixInBlock(const Matrix & m, size_t i, size_t j, Matrix & ret)
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
}*/


ViscoelasticityAndFracture::ViscoelasticityAndFracture(ViscoelasticModel m, const Matrix & rig, FractureCriterion * c, DamageModel * d, int n, double r) : Viscoelasticity(m, rig, n, r), coeffsDamageElastic(blocks, blocks), coeffsDamageViscous(blocks, blocks), criterion(c), dfunc(d)
{
	for(size_t i = 0 ; i < blocks ; i++)
	{
		for(size_t j = 0 ; j < blocks ; j++)
		{
			coeffsDamageViscous[i][j] = 1. ;
			coeffsDamageElastic[i][j] = 1. ;
		}
	}
} ;

ViscoelasticityAndFracture::ViscoelasticityAndFracture(ViscoelasticModel m, const Matrix & rig, const Matrix & e, FractureCriterion * c, DamageModel * d, int b, int n, double r) : Viscoelasticity(m, rig, e, b, n, r), coeffsDamageElastic(blocks, blocks), coeffsDamageViscous(blocks, blocks), criterion(c), dfunc(d)
{ 
	for(size_t i = 0 ; i < blocks ; i++)
	{
		for(size_t j = 0 ; j < blocks ; j++)
		{
			coeffsDamageViscous[i][j] = 1. ;
			coeffsDamageElastic[i][j] = 1. ;
		}
	}
} ;

ViscoelasticityAndFracture::ViscoelasticityAndFracture(const Matrix & rig, const Matrix & e, int b, FractureCriterion * c, DamageModel * d, int n, double r) : Viscoelasticity(rig, e, b, n, r),coeffsDamageElastic(blocks, blocks), coeffsDamageViscous(blocks, blocks), criterion(c), dfunc(d)
{
	for(size_t i = 0 ; i < blocks ; i++)
	{
		for(size_t j = 0 ; j < blocks ; j++)
		{
			coeffsDamageViscous[i][j] = 1. ;
			coeffsDamageElastic[i][j] = 1. ;
		}
	}
} ;

ViscoelasticityAndFracture::ViscoelasticityAndFracture(ViscoelasticModel m, const Matrix & c_kv, const Matrix & e_kv, const Matrix & c_mx, const Matrix & e_mx,FractureCriterion * c, DamageModel * d,  int b, int n, double r) : Viscoelasticity(m, c_kv, e_kv, c_mx, e_mx, b, n, r), coeffsDamageElastic(blocks, blocks), coeffsDamageViscous(blocks, blocks), criterion(c), dfunc(d)
{
	param.resize(c_kv.numRows()*(3+n+b), c_kv.numCols()*(3+n+b)) ;
	eta.resize(c_kv.numRows()*(3+n+b), c_kv.numCols()*(3+n+b)) ;
	param.array() = 0 ;
	eta.array() = 0 ;
} ;

ViscoelasticityAndFracture::ViscoelasticityAndFracture(ViscoelasticModel m, const Matrix & c0, std::vector<std::pair<Matrix, Matrix> > & branches, FractureCriterion * c, DamageModel * d, int b, int n, double r) : Viscoelasticity(m, c0, branches, b, n, r), coeffsDamageElastic(blocks, blocks), coeffsDamageViscous(blocks, blocks), criterion(c), dfunc(d)
{
	for(size_t i = 0 ; i < blocks ; i++)
	{
		for(size_t j = 0 ; j < blocks ; j++)
		{
			coeffsDamageViscous[i][j] = 1. ;
			coeffsDamageElastic[i][j] = 1. ;
		}
	}
} ;

ViscoelasticityAndFracture::ViscoelasticityAndFracture(ViscoelasticModel m, const Matrix & c0, const Matrix & c1, const Matrix & e1, FractureCriterion * c, DamageModel * d, int b, int n, double r) : Viscoelasticity(m, c0, c1, e1, b, n, r), coeffsDamageElastic(blocks, blocks), coeffsDamageViscous(blocks, blocks), criterion(c), dfunc(d)
{
	for(size_t i = 0 ; i < blocks ; i++)
	{
		for(size_t j = 0 ; j < blocks ; j++)
		{
			coeffsDamageViscous[i][j] = 1. ;
			coeffsDamageElastic[i][j] = 1. ;
		}
	}
} ;

ViscoelasticityAndFracture::~ViscoelasticityAndFracture()
{ 
	delete dfunc ;
	delete criterion ;
} 

ElementState * ViscoelasticityAndFracture::createElementState( IntegrableEntity * e) 
{
	return new GeneralizedSpaceTimeViscoElasticElementState(e) ;  
}

FractureCriterion * ViscoelasticityAndFracture::getFractureCriterion() const
{
	return criterion ;
}

DamageModel * ViscoelasticityAndFracture::getDamageModel() const
{
	return dfunc ;
}

void ViscoelasticityAndFracture::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	Matrix a(ret.numRows()/blocks, ret.numCols()/blocks) ;
	Matrix b(ret.numRows()/blocks, ret.numCols()/blocks) ;
	if(!time_d)
	{
		
		Matrix buffer(param.numRows()/blocks, param.numCols()/blocks) ;

		switch(model)
		{
			case GENERALIZED_KELVIN_VOIGT:
			{
				// stiffness (0,0)
				getBlockInMatrix(param, 0,0, buffer) ;
				buffer = dfunc->apply(buffer) ;
				vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
				vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
				a += b ;
				placeMatrixInBlock( a * coeffsDamageElastic[0][0], 0,0, ret ) ;
				for(size_t i = 1 ; i < blocks ; i++)
				{
					// first line
					substractMatrixInBlock( a * coeffsDamageElastic[i][0], i,0, ret ) ;
					// first column
					substractMatrixInBlock( a * coeffsDamageElastic[0][i], 0,i, ret ) ;
					for(size_t j = i+1 ; j < blocks ; j++)
					{
						// upper triangle
						placeMatrixInBlock( a * coeffsDamageElastic[i][j], i,j, ret ) ;
						// lower triangle
						placeMatrixInBlock( a * coeffsDamageElastic[j][i], j,i, ret ) ;
					}
				}
				for(size_t i = 1 ; i < blocks ; i++)
				{
					//stiffness (diagonal)
					getBlockInMatrix(param, i,i, buffer) ;
					buffer = dfunc->apply(buffer) ;
					vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
					vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
					a += b ;
					placeMatrixInBlock( a * coeffsDamageElastic[i][i], i,i, ret ) ;
				}
				
				return ;
			}
			
			case GENERALIZED_MAXWELL:
			{
				// stiffness (0,0)
				getBlockInMatrix(param, 0,0, buffer) ;
				buffer = dfunc->apply(buffer) ;
				vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
				vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
				a += b ;
				placeMatrixInBlock( a * coeffsDamageElastic[0][0], 0,0, ret ) ;
				for(size_t i = 1 ; i < blocks ; i++)
				{
					//stiffness (diagonal)
					getBlockInMatrix(param, i,i, buffer) ;
					buffer = dfunc->apply(buffer) ;
					vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
					vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
					a += b ;
					placeMatrixInBlock( a * coeffsDamageElastic[i][i], i,i, ret ) ;
					// first line
					substractMatrixInBlock( a * coeffsDamageElastic[i][0], i,0, ret ) ;
					// first column
					substractMatrixInBlock( a * coeffsDamageElastic[0][i], 0,i, ret ) ;
					
				}
				return ;
			}
			
			case BURGER:
			{
				// stiffness Maxwell
				getBlockInMatrix(param, 0,0, buffer) ;
				buffer = dfunc->apply(buffer) ;
				vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
				vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
				a += b ;
				placeMatrixInBlock( a * coeffsDamageElastic[0][0], 0,0, ret ) ;
				placeMatrixInBlock( a * coeffsDamageElastic[1][1], 1,1, ret ) ;
				placeMatrixInBlock( a * coeffsDamageElastic[1][2], 1,2, ret ) ;
				placeMatrixInBlock( a * coeffsDamageElastic[2][1], 2,1, ret ) ;
				substractMatrixInBlock( a * coeffsDamageElastic[0][1], 0,1, ret ) ;
				substractMatrixInBlock( a * coeffsDamageElastic[1][0], 1,0, ret ) ;
				substractMatrixInBlock( a * coeffsDamageElastic[0][2], 0,2, ret ) ;
				substractMatrixInBlock( a * coeffsDamageElastic[2][0], 2,0, ret ) ;
				// stiffness KV
				getBlockInMatrix(param, 2,2, buffer) ;
				buffer = dfunc->apply(buffer) ;
				vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
				vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
				a += b ;
				placeMatrixInBlock( a * coeffsDamageElastic[2][2], 2,2, ret ) ;

				return ;
			}
			
			case MAXWELL:
			{
				// stiffness
				getBlockInMatrix(param, 0,0, buffer) ;
				buffer = dfunc->apply(buffer) ;
				vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
				vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
				a += b ;
				placeMatrixInBlock( a * coeffsDamageElastic[0][0], 0,0, ret ) ;
				placeMatrixInBlock( a * coeffsDamageElastic[1][1], 1,1, ret ) ;
				substractMatrixInBlock( a * coeffsDamageElastic[0][1], 0,1, ret ) ;
				substractMatrixInBlock( a * coeffsDamageElastic[1][0], 1,0, ret ) ;
				return ;
			}
			
			case KELVIN_VOIGT:
			{
				// stiffness
				getBlockInMatrix(param, 0,0, buffer) ;
				buffer = dfunc->apply(buffer) ;
				vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
				vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
				a += b ;
				placeMatrixInBlock( a * coeffsDamageElastic[0][0], 0,0, ret ) ;
				return ;
			}
			
			case PURE_ELASTICITY:
			{
				getBlockInMatrix(param, 0,0, buffer) ;
				buffer = dfunc->apply(buffer) ;
				vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
				vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
				a += b ;
				placeMatrixInBlock( a * coeffsDamageElastic[0][0], 0,0, ret ) ;
				
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
					buffer = dfunc->apply(buffer) ;
					vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
					vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
					a += b ;
					placeMatrixInBlock( a * coeffsDamageElastic[i][i], i,i, ret ) ;
					// elastic matrix (upper-triangle)
					for(size_t j = i+1 ; j < blocks ; j++)
					{
						getBlockInMatrix(param, i,j, buffer) ;
						buffer = dfunc->apply(buffer) ;
						vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
						vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
						a += b ;
						placeMatrixInBlock( a * coeffsDamageElastic[i][j], i,j, ret ) ;
						// symmetry
						placeMatrixInBlock( a * coeffsDamageElastic[j][i], j,i, ret ) ;
					}

				}
			}
		}
		
		return ;
	}
	else
	{
		std::vector<Matrix> tensor ;
		for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
			tensor.push_back(Matrix(param.numRows()/blocks, param.numCols()/blocks)) ;
		
		for(size_t i = 0 ; i < blocks ; i++)
		{
			for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
			{
				getBlockInMatrix( param, i,i, tensor[g] ) ;
				tensor[g] = dfunc->apply(tensor[g], gp.gaussPoints[g].first ) ;
			}
			
			vm->ieval(GradientDot(p_i) * tensor * Gradient(p_j, true),    gp, Jinv,v, a) ;
			vm->ieval(Gradient(p_i)    * tensor * GradientDot(p_j, true), gp, Jinv,v, b) ;
			a += b ;
			placeMatrixInBlock( a * coeffsDamageElastic[i][i], i,i, ret ) ;

			for(size_t j = i+1 ; j < blocks ; j++)
			{
				for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
				{
					getBlockInMatrix( param, i,j, tensor[g] ) ;
					tensor[g] = dfunc->apply(tensor[g], gp.gaussPoints[g].first ) ;
				}
				
				vm->ieval(GradientDot(p_i) * tensor * Gradient(p_j, true),    gp, Jinv,v, a) ;
				vm->ieval(Gradient(p_i)    * tensor * GradientDot(p_j, true), gp, Jinv,v, b) ;
				a += b ;
				placeMatrixInBlock( a * coeffsDamageElastic[i][j], i,j, ret ) ;
				placeMatrixInBlock( a * coeffsDamageElastic[j][i], j,i, ret ) ;
			}	
		}
	  
	}
	
  
}

void ViscoelasticityAndFracture::applyViscous(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
		Matrix a(ret.numRows()/blocks, ret.numCols()/blocks) ;
		Matrix b(ret.numRows()/blocks, ret.numCols()/blocks) ;
	if(!time_d)
	{

		Matrix buffer(param.numRows()/blocks, param.numCols()/blocks) ;

		switch(model)
		{
			case GENERALIZED_KELVIN_VOIGT:
			{
				for(size_t i = 1 ; i < blocks ; i++)
				{
					// viscosity (diagonal)
					getBlockInMatrix(eta, i,i, buffer) ;
					buffer = dfunc->apply(buffer) ;
					vm->ieval(GradientDot(p_i) * buffer * GradientDot(p_j, true),    gp, Jinv,v, a) ;
					vm->ieval(GradientDotDot(p_i)    * buffer * Gradient(p_j, true), gp, Jinv,v, b) ;
					a += b ;
					placeMatrixInBlock( a * coeffsDamageViscous[i][i], i,i, ret ) ;
				}
				
				return ;
			}
			
			case GENERALIZED_MAXWELL:
			{
				for(size_t i = 1 ; i < blocks ; i++)
				{
					// viscosity (diagonal)
					getBlockInMatrix(eta, i,i, buffer) ;
					buffer = dfunc->apply(buffer) ;
					vm->ieval(GradientDot(p_i) * buffer * GradientDot(p_j, true),    gp, Jinv,v, a) ;
					vm->ieval(GradientDotDot(p_i)    * buffer * Gradient(p_j, true), gp, Jinv,v, b) ;
					a += b ;
					placeMatrixInBlock( a * coeffsDamageViscous[i][i], i,i, ret ) ;
				}
				return ;
			}
			
			case BURGER:
			{
				// viscosity Maxwell
				getBlockInMatrix(eta, 1,1, buffer) ;
				buffer = dfunc->apply(buffer) ;
				vm->ieval(GradientDot(p_i)    * buffer   * GradientDot(p_j, true), gp, Jinv,v,a);
				vm->ieval(GradientDotDot(p_i) * buffer   * Gradient(p_j, true),    gp, Jinv,v,b);
				a += b ;
				placeMatrixInBlock( a * coeffsDamageViscous[1][1], 1,1, ret ) ;
				// viscosity KV
				getBlockInMatrix(eta, 2,2, buffer) ;
				vm->ieval(GradientDot(p_i)    * buffer   * GradientDot(p_j, true), gp, Jinv,v,a);
				vm->ieval(GradientDotDot(p_i) * buffer   * Gradient(p_j, true),    gp, Jinv,v,b);
				a += b ;
				placeMatrixInBlock( a * coeffsDamageViscous[2][2], 2,2, ret ) ;
				return ;
			}
			
			case MAXWELL:
			{
				// viscosity
				getBlockInMatrix(eta, 1,1, buffer) ;
				buffer = dfunc->apply(buffer) ;
				vm->ieval(GradientDot(p_i)    * buffer   * GradientDot(p_j, true), gp, Jinv,v,a);
				vm->ieval(GradientDotDot(p_i) * buffer   * Gradient(p_j, true),    gp, Jinv,v,b);
				a += b ;
				placeMatrixInBlock( a * coeffsDamageViscous[1][1], 1,1, ret ) ;
				return ;
			}
			
			case KELVIN_VOIGT:
			{
				// viscosity
				getBlockInMatrix(eta, 0,0, buffer) ;
				buffer = dfunc->apply(buffer) ;
				vm->ieval(GradientDot(p_i)    * buffer   * GradientDot(p_j, true), gp, Jinv,v,a);
				vm->ieval(GradientDotDot(p_i) * buffer   * Gradient(p_j, true),    gp, Jinv,v,b);
				addMatrixInBlock( a * coeffsDamageViscous[0][0], 0,0, ret ) ;
				return ;
			}
			
			case PURE_ELASTICITY:
			{
				return ;
			}
			
			case PURE_VISCOSITY:
			{
				getBlockInMatrix(eta, 0,0, buffer) ;
				buffer = dfunc->apply(buffer) ;
				vm->ieval(GradientDot(p_i)    * buffer   * GradientDot(p_j, true), gp, Jinv,v,a);
				vm->ieval(GradientDotDot(p_i) * buffer   * Gradient(p_j, true),    gp, Jinv,v,b);
				a += b ;
				placeMatrixInBlock( a * coeffsDamageViscous[0][0], 0,0, ret ) ;
				
				return ;
			}
			
			default:
			{
				for(size_t i = 0 ; i < blocks ; i++)
				{
					// viscous matrix (diagonal)
					getBlockInMatrix(eta, i,i, buffer) ;
					buffer = dfunc->apply(buffer) ;
					vm->ieval(GradientDot(p_i)    * buffer   * GradientDot(p_j, true), gp, Jinv,v,a);
					vm->ieval(GradientDotDot(p_i) * buffer   * Gradient(p_j, true),    gp, Jinv,v,b);
					a += b ;
					placeMatrixInBlock( a * coeffsDamageViscous[i][i], i,i, ret ) ;
					// viscous matrix (upper-triangle)
					for(size_t j = i+1 ; j < blocks ; j++)
					{
						getBlockInMatrix(eta, i,j, buffer) ;
						buffer = dfunc->apply(buffer) ;
						vm->ieval(GradientDot(p_i)    * buffer   * GradientDot(p_j, true), gp, Jinv,v,a);
						vm->ieval(GradientDotDot(p_i) * buffer   * Gradient(p_j, true),    gp, Jinv,v,b);
						a += b ;
						addMatrixInBlock( a * coeffsDamageViscous[i][j], i,j, ret ) ;
						// symmetry
						addMatrixInBlock( a * coeffsDamageViscous[j][i], j,i, ret ) ;
					}
				}
			}
		}
		
		return ;
	}
	else
	{
		std::vector<Matrix> tensor ;
		for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
			tensor.push_back(Matrix(param.numRows()/blocks, param.numCols()/blocks)) ;
		
		for(size_t i = 0 ; i < blocks ; i++)
		{
			for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
			{
				getBlockInMatrix( eta, i,i, tensor[g] ) ;
				tensor[g] = dfunc->apply(tensor[g], gp.gaussPoints[g].first ) ;
			}
			
			vm->ieval(GradientDotDot(p_i) * tensor * Gradient(p_j, true),    gp, Jinv,v, a) ;
			vm->ieval(GradientDot(p_i)    * tensor * GradientDot(p_j, true), gp, Jinv,v, b) ;
			a += b ;
			placeMatrixInBlock( a * coeffsDamageViscous[i][i], i,i, ret ) ;

			for(size_t j = i+1 ; j < blocks ; j++)
			{
				for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
				{
					getBlockInMatrix( eta, i,j, tensor[g] ) ;
					tensor[g] = dfunc->apply(tensor[g], gp.gaussPoints[g].first ) ;
				}
				
				vm->ieval(GradientDotDot(p_i) * tensor * Gradient(p_j, true),    gp, Jinv,v, a) ;
				vm->ieval(GradientDot(p_i)    * tensor * GradientDot(p_j, true), gp, Jinv,v, b) ;
				a += b ;
				placeMatrixInBlock( a * coeffsDamageViscous[i][j], i,j, ret ) ;
				placeMatrixInBlock( a * coeffsDamageViscous[j][i], j,i, ret ) ;
			}	
		}
	}
  
	
  
}

void ViscoelasticityAndFracture::step(double timestep, ElementState & currentState, double maxscore) 
{
	dfunc->step(currentState, maxscore) ;
	currentState.getParent()->behaviourUpdated = dfunc->changed() ;
	SpaceTimeNonLocalMaximumStrain * crit = dynamic_cast<SpaceTimeNonLocalMaximumStrain *>(criterion) ;
	if(crit != nullptr && dfunc->changed() && !dfunc->fractured())
	{
		crit->upVal = crit->maxstress/(param[0][0]*(1.-dfunc->getState().max())) ;
	}


}


bool ViscoelasticityAndFracture::fractured() const
{
	return dfunc->fractured() ;
}

bool ViscoelasticityAndFracture::changed() const
{
	return dfunc->changed() ;
} 

Form * ViscoelasticityAndFracture::getCopy() const 
{
  
	ViscoelasticityAndFracture * copy = new ViscoelasticityAndFracture(  param, eta, blocks, criterion->getCopy(), dfunc->getCopy())  ;
	copy->model = model ;
	copy->dfunc->getState(true).resize(dfunc->getState().size());
	copy->dfunc->getState(true) = dfunc->getState() ;
	copy->criterion->setMaterialCharacteristicRadius(criterion->getMaterialCharacteristicRadius()) ;
	copy->dfunc->setDamageDensityTolerance(dfunc->getDamageDensityTolerance());
	copy->dfunc->setThresholdDamageDensity(dfunc->getThresholdDamageDensity());
	
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

Vector ViscoelasticityAndFracture::getForcesFromAppliedStress( Vector & data, Function & shape, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<Variable> & v) 
{
	return VirtualMachine().ieval(GradientDot( shape ) * ( data ), gp, Jinv, v) ;
}

Vector ViscoelasticityAndFracture::getForcesFromAppliedStress( const Function & data, size_t index, size_t externaldofs,  Function & shape, IntegrableEntity * e,const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<Variable> & v) 
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

Matrix ViscoelasticityAndFracture::getTensor(const Point & p, IntegrableEntity * e, int g) const
{
//	return dfunc->apply(param) ; 
	Matrix tensor(param.numRows(), param.numCols()) ;
	Matrix buffer(param.numRows()/blocks, param.numCols()/blocks) ;
	Matrix tmp = dfunc->apply(param) ;
	for(size_t i = 0 ; i < blocks ; i++)
	{
		for(size_t j = 0 ; j < blocks ; j++)
		{
			getBlockInMatrix(tmp, i,j, buffer) ;
			placeMatrixInBlock( buffer * coeffsDamageElastic[i][j], i, j, tensor) ;
		}
	}
	return tensor ;
}

Matrix ViscoelasticityAndFracture::getViscousTensor(const Point & p, IntegrableEntity * e, int g) const
{
//	return dfunc->apply(eta) ; 
	Matrix tensor(param.numRows(), param.numCols()) ;
	Matrix buffer(param.numRows()/blocks, param.numCols()/blocks) ;
	Matrix tmp = dfunc->apply(eta) ;
	for(size_t i = 0 ; i < blocks ; i++)
	{
		for(size_t j = 0 ; j < blocks ; j++)
		{
			getBlockInMatrix(tmp, i,j, buffer) ;
			placeMatrixInBlock( buffer * coeffsDamageViscous[i][j], i, j, tensor) ;
		}
	}
	return tensor ;
}

void ViscoelasticityAndFracture::setFractureCriterion(FractureCriterion * frac) 
{
	if(frac)
	{
		criterion = frac ;
	}
	
}

