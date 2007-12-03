
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//

#include "matrixops.h"

using namespace Mu ;

Matrix::Matrix(size_t x, size_t y)
{
	r=x ;
	c=y ;
	v = new Vector(0., x*y) ;
}

Matrix::Matrix(size_t rl, size_t k, size_t l, const Matrix & m)
{
	r=rl ;
	c=rl ;
	v = new Vector(rl*rl) ;
	for(size_t i = 0 ; i< rl ; i++)
	{
		for(size_t j = 0 ; j< rl ; j++)
		{
			(*v)[i*rl+j] = m[k+i][l+j] ;
		}
	}
}

Matrix::Matrix(size_t i, size_t j, const Matrix & m)
{
	r=m.numRows()-1 ;
	c=m.numCols()-1 ;
	v = new Vector(r*c) ;
	for(size_t l = 0 ; l< m.numRows() ; l++)
	{
		if(l != i)
		{
			for(size_t k = 0 ; k< m.numCols() ; k++)
			{
				if(k != j)
				{
					size_t r_index = l ;
					if(l > i)
						r_index-- ;
					size_t c_index = k ;
					if(k > j)
						c_index-- ;
					
					(*this)[r_index][c_index] = m[l][k] ;
				}
			}
		}
	}
}

double& Matrix::operator()(size_t x, size_t y)
{
	return row(x)[y] ;
}

double Matrix::operator()(size_t x, size_t y) const
{
	return row(x)[y] ;
}

Matrix Matrix::transpose() const
{
	Matrix ret(c,r) ;
	for(size_t i = 0 ; i < r ; i++)
	{
		for(size_t j = 0 ; j < c ; j++)
		{
			ret[j][i] = (*this)[i][j] ;
		}
	}
	
	return ret ;
}

Matrix &Matrix::operator*=(double d)
{
	(*v) *=d ;
	return *this ;
}

Matrix Matrix::operator*(double d) const
{
	Matrix ret(*this) ;
	ret*=d ;
		
	return ret ;
}

Matrix Matrix::operator/(double d) const
{
	Matrix ret(*this) ;
	ret /=d ;
	return ret ;
}

Matrix &Matrix::operator/=(double d)
{
	(*v) /=d ;
	return *this ;
}

Matrix & Matrix::operator =(const MtM& m)
{
	(*this) = (Matrix) m ;
	
	return *this ;
}

Matrix &Matrix::operator*=(const Matrix &m)
{
	assert(m.numRows() == this->numCols()) ;
	
	Matrix ret =matrix_multiply((*this), m) ;
	
	(*v) = ret.array() ;
	r = ret.numRows() ;
	c = ret.numCols() ;
	
	return *this ;
}


Matrix &Matrix::operator +=(const Matrix &m)
{
	(*v) += (m.array()) ;
	return *this ;
}

Matrix Matrix::operator +(const  Matrix &m) const
{
	Matrix ret(*this) ;
	ret += m ; ;
	return ret ;
}

Matrix &Matrix::operator -=(const Matrix &m)
{
	(*v) -= (m.array()) ;
	return *this ;
}

Matrix Matrix::operator -(const  Matrix &m) const
{
	Matrix ret(*this) ;
	ret -= m ; ;
	return ret ;
}

bool Matrix::operator ==(const Matrix &m)
{
	if(v->size() != m.array().size())
		return false ;
	
	for(size_t i = 0 ; i < v->size() ; i++)
		if((*v)[i] != m.array()[i])
			return false ;
	
	return true ;
}

bool Matrix::operator !=(const Matrix &m)
{
	if(v->size() != m.array().size())
		return true ;
	else
		for(size_t i = 0 ; i < v->size() ; i++)
			if((*v)[i] != m.array()[i])
				return true ;
	
	return false ;
}

Matrix::Matrix(const Matrix& m) : r(m.numRows()), c( m.numCols())
{
	v = new Vector(m.array()) ;
}

Matrix &Matrix::operator =(const Matrix &m)
{
	delete v ;
	v = new Vector(m.array()) ;
	r = m.numRows() ;
	c = m.numCols() ;
	return *this ;
}

void Matrix::print() const
{
	for(size_t i = 0 ; i < numRows() ;i++)
	{
		for(size_t j = 0 ; j < numCols() ;j++)
		{
			std::cout << (*this)[i][j] << "   " << std::flush ;
		}
		
		std::cout << std::endl ;
	}
}

MtM::operator const Matrix() const
{
	return matrix_multiply(first, second) ;
}


MtV::operator const Vector()
{
	return matrix_vector_multiply(m,v) ;
}


void swapLines(size_t l0, size_t l1, Mu::Matrix * m)
{
	for(size_t i = 0 ; i < m->numCols() ; i++)
	{
		double temp = (*m)[l1][i] ;
		(*m)[l1][i] = (*m)[l0][i] ;
		(*m)[l0][i] = temp ;
// 		std::swap(m[l0][i], m[l1][i]) ;
	}
}

Mu::Matrix swapMatrix(size_t l0, size_t l1, size_t dim)
{
	
	Mu::Matrix ret(dim, dim) ;
	for(size_t i = 0 ;  i < dim ; i++)
	{
		if(i!=l0 && i!=l1 )
		{
			ret[i][i] = 1 ;
		}
		else if(i == l0)
		{
			ret[l0][l1] = 1 ;
		}
		else if(i == l1)
		{
			ret[l1][l0] = 1 ;
		}
	}
	
	return ret ;
	
}

//clever 2x2 Matrix inversion. Thanks the numerical cookbook :)
Mu::Matrix inverse2x2Matrix(const Mu::Matrix &s)
{
	Matrix s_(s) ;
	invert2x2Matrix(s_) ;
	return s_ ;
}

void invert2x2Matrix(Mu::Matrix &s)
{
	
	if(std::abs(s[0][0]) < 1e-16)
	{
		Matrix swap(2,2) ;
		swap[0][1] = 1 ;
		swap[1][0] = 1 ;
		s *= swap ;
		
		double r1 = 1./s[0][0] ;
		double r2 = s[1][0] * r1 ;
		double r3 = r1* s[0][1] ;
		double r4 = s[1][0] * r3 ;
		double r5 = r4 - s[1][1] ;
		double r6 = 1./r5 ;
		s[0][1] = r3*r6 ;
		s[1][0] = r6*r2 ;
		double r7 = r3*s[1][0] ;
		s[0][0] = r1-r7 ;
		s[1][1] = -r6 ;
		
		s = swap*s ;
		return ;
	}
	
	double r1 = 1./s[0][0] ;
	double r2 = s[1][0] * r1 ;
	double r3 = r1* s[0][1] ;
	double r4 = s[1][0] * r3 ;
	double r5 = r4 - s[1][1] ;
	double r6 = 1./r5 ;
	s[0][1] = r3*r6 ;
	s[1][0] = r6*r2 ;
	double r7 = r3*s[1][0] ;
	s[0][0] = r1-r7 ;
	s[1][1] = -r6 ;

}


Mu::Matrix inverse3x3Matrix(const Mu::Matrix & m)
{
	Mu::Matrix m_(m) ;
	invert3x3Matrix(m_) ;
	return m_ ;
}

void invert3x3Matrix(Mu::Matrix & m)
{

	
	Vector v(9) ;
	double m00 = m.array()[0] ; double m01 =  m.array()[1] ; double m02 =  m.array()[2] ; 
	double m10 = m.array()[3] ; double m11 =  m.array()[4] ; double m12 =  m.array()[5] ; 
	double m20 = m.array()[6] ; double m21 =  m.array()[7] ; double m22 =  m.array()[8] ; 
	
	double r11 = fma(m11,m22,-m12*m21) ;
	double r21 = fma(m12,m20,-m10*m22) ;
	double r31 = fma(m10,m21,-m11*m20) ;
	double det = fma(m00,(r11),fma( m01,r21, m02*r31)) ;
	
	v[0] = r11 ; 
	v[1] = fma(m02,m21, - m01*m22) ;
	v[2] = fma(m01,m12, - m02*m11) ;
	
	v[3] = r21 ; 
	v[4] = fma(m00,m22, - m02*m20) ;
	v[5] = fma(m02,m10, - m00*m12) ;
	
	v[6] = r31 ; 
	v[7] = fma(m01,m20, - m00*m21) ;
	v[8] = fma(m11,m00, - m01*m10) ;
	v /= det ;
	m.array() = v ;

	
}

Mu::Matrix inverse4x4Matrix(const Mu::Matrix &s)
{	
	Mu::Matrix ret(4,4) ;
	
	
	double deti = det(s) ;
	
	ret[0][0] =  det(Matrix(0,0,s)) ; ret[0][1] = -det(Matrix(1,0,s)) ; ret[0][2] =  det(Matrix(2,0,s)) ; ret[0][3] = -det(Matrix(3,0,s)) ;
	ret[1][0] = -det(Matrix(0,1,s)) ; ret[1][1] =  det(Matrix(1,1,s)) ; ret[1][2] = -det(Matrix(2,1,s)) ; ret[1][3] =  det(Matrix(3,1,s)) ;
	ret[2][0] =  det(Matrix(0,2,s)) ; ret[2][1] = -det(Matrix(1,2,s)) ; ret[2][2] =  det(Matrix(2,2,s)) ; ret[2][3] = -det(Matrix(3,2,s)) ;
	ret[3][0] = -det(Matrix(0,3,s)) ; ret[3][1] =  det(Matrix(1,3,s)) ; ret[3][2] = -det(Matrix(2,3,s)) ; ret[3][3] =  det(Matrix(3,3,s)) ;
	
	return ret/deti ;

}


double det(const Mu::Matrix &s) 
{
	switch(s.numCols())
	{		
		case 1:
		{
			return s[0][0];
		}
		case 2:
		{
			return fma(s[0][0],s[1][1],-s[0][1]*s[1][0]); 
		}
		case 3:
		{
			
			return  fma(s[0][0],s[1][1]*s[2][2],fma(  
							s[1][2],s[2][0]*s[0][1], fma( 
							s[0][2],s[2][1]*s[1][0], fma(- 
							s[0][0],s[1][2]*s[2][1], fma(- 
							s[0][1],s[1][0]*s[2][2], - 
							s[2][0]*s[1][1]*s[0][2])))));
		}
	
		case 4:
		{
			double ret = 0 ;
			
			for(size_t i = 0 ;  i < 1  ; i++)
			{
				for(size_t j = 0 ;  j< 4  ; j++)
				{
					double sig = 1 ;
					for(size_t k = 0 ; k < (i+j)%2 ; k++)
						sig*=-1 ;
					
					Matrix mat(i,j,s) ;
					
					ret= fma(sig*s[i][j],det(mat), ret) ;
				}
			}
			
			return ret ;
		}
	default:
		{
			return 0 ;
		}
	}
}


Mu::Matrix exp(const Mu::Matrix& m, size_t order)
{
	Matrix Identity(m.numRows(), m.numCols()) ;
	Matrix ret(m.numRows(), m.numCols()) ;
	
	for(size_t i = 0 ; i < m.numRows() ; i++)
	{
		Identity[i][i] = 1 ;
	}
	
	for(size_t i=0 ; i < order ; i++)
	{
		Matrix temp(m) ;
		
		for(size_t j = 0 ; j < i ; j++)
		{
			temp*=(m)/(j+2) ;
		}

		ret+= temp ;
	}
	
	return ret+Identity ;
}

Mu::Matrix log(const Mu::Matrix&m, size_t order)
{
	Matrix Identity(m.numRows(), m.numCols()) ;
	
	for(size_t i = 0 ; i < m.numRows() ; i++)
	{
		Identity[i][i] = 1 ;
	}
	
	Matrix ret(m-Identity) ;

	
	for(size_t i=1 ; i < order ; i++)
	{
		Matrix temp(m-Identity) ;
		
		for(size_t j = 0 ; j < i ; j++)
		{
			temp*=(m-Identity) ;
		}
		if(i%2 == 0)
			ret+= temp/(double)(i+1) ;
		else
			ret-= temp/(double)(i+1) ;
	}
	
	return ret ;
}

Vector Mu::solveSystem(const Matrix & A, const Vector & b, Vector & x)
{
	
	double epsilon = 1e-12 ;
	Vector r = b - (Vector)(A*x) ;
	Vector r_(r) ;
	double rho = std::inner_product(&r[0], &r[r.size()], &r_[0], double(0)) ;
	
	Vector invDiag(r.size()) ;
	
	for(size_t i = 0 ; i < r.size() ; i++)
	{
		invDiag[i] = 1./A[i][i] ;
// 		invDiag[i] = 1 ;
	}
	
	if(rho < epsilon)
		return x ;
	
	Vector p(r) ;
	
	Vector p_ = invDiag*p ;
	//Vector p_ = precondition(p) ;
	
	Vector v = A*p_ ;
	
	double alpha = rho/std::inner_product(&r_[0], &r_[r_.size()], &v[0], double(0)) ;
	
	Vector s = r - v*alpha ;
	
	if(std::abs(s.max()) < epsilon)
	{
		x += p_*alpha ;
		return x ;
	}
	
	Vector s_ = invDiag*s ;
	//Vector s_ = precondition(s) ;
	Vector t = A*s_ ;
	double omega = std::inner_product(&t[0], &t[t.size()], &s[0], double(0))/std::inner_product(&t[0], &t[t.size()], &t[0], double(0)) ;
	x += p_*alpha +omega*s_ ;
	r = s- t*omega ;
	double rho_ =rho ;
	
	size_t nit = 0 ;
	
	while(nit < r.size())
	{
		std::cout << "nit = "<< nit << std::endl ;
		nit++ ;
		
		rho = std::inner_product(&r[0], &r[r.size()], &r_[0], double(0)) ;
		if(std::abs(rho) < epsilon)
			return x ;
		
		double beta = (rho/rho_)*(alpha/omega) ;
		p = r + (p-v*omega)*beta ;
		p_ = invDiag*p ;
		//p_ = precondition(p) ;
		v = A*p_ ;
		alpha = rho/std::inner_product(&r_[0], &r_[r_.size()], &v[0], double(0)) ;
		s = r - v*alpha ;
		
		//s_ = precondition(s) ;
		s_ = invDiag*s ;
		t = A*s_ ;
		omega = std::inner_product(&t[0], &t[t.size()], &s[0], double(0))/std::inner_product(&t[0], &t[t.size()], &t[0], double(0)) ;
		
		if(std::abs(omega) < epsilon)
			return x ;
		
		x += p_*alpha +s_*omega ;
		r = s- t*omega ;
		rho_ = rho ;
		
	}
	
	return x ;
}

