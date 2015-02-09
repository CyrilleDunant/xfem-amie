//
// C++ Implementation: vm_function_matrix
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "vm_function_matrix.h"

namespace Amie
{

FunctionMatrix::FunctionMatrix(size_t x, size_t y)
{
	r=x ;
	c=y ;
	v = new std::valarray<Function>(x*y) ;
}

Function& FunctionMatrix::operator()(size_t x, size_t y)
{
	return row(x)[y] ;
}

Function FunctionMatrix::operator()(size_t x, size_t y) const
{
	return row(x)[y] ;
}

FunctionMatrix FunctionMatrix::transpose() const
{
	FunctionMatrix ret(c,r) ;
	for(size_t i = 0 ; i < r ; i++)
	{
		for(size_t j = 0 ; j < c ; j++)
		{
			ret[j][i] = (*this)[i][j] ;
		}
	}
	
	return ret ;
}

FunctionMatrix & FunctionMatrix::operator*=(const Function &d)
{
	(*v) *=d ;
	return *this ;
}

FunctionMatrix & FunctionMatrix::operator*=(const double d)
{
	for(size_t i = 0 ;  i < v->size() ; i++)
	{
		(*v)[i] *= d ;
	}
	return *this ;
}

FunctionMatrix FunctionMatrix::operator*(const Function & d) const
{
	FunctionMatrix ret(*this) ;
	ret*=d ;
	
	return ret ;
}

FunctionMatrix FunctionMatrix::operator*(const double d) const
{
	FunctionMatrix ret(*this) ;
	ret *=d ;
	
	return ret ;
}

FunctionMatrix FunctionMatrix::operator/(const Function & d) const
{
	FunctionMatrix ret(*this) ;
	ret /=d ;
	return ret ;
}

FunctionMatrix &FunctionMatrix::operator/=(const Function & d)
{
	(*v) /=d ;
	return *this ;
}

FunctionMatrix & FunctionMatrix::operator =(const FMtFM& m)
{
	(*this) = (FunctionMatrix) m ;
	
	return *this ;
}

FunctionMatrix & FunctionMatrix::operator =(const FMtM& m)
{
	(*this) = (FunctionMatrix) m ;
	
	return *this ;
}

FunctionMatrix &FunctionMatrix::operator*=(const FunctionMatrix &m)
{
	assert(m.numRows() == this->numCols()) ;
	
	FunctionMatrix ret =ff_matrix_multiply((*this), m) ;
	
	(*v) = ret.array() ;
	r = ret.numRows() ;
	c = ret.numCols() ;
	
	return *this ;
}

FunctionMatrix &FunctionMatrix::operator*=(const Matrix &m)
{
	assert(m.numRows() == this->numCols()) ;
	
	FunctionMatrix ret =fm_matrix_multiply((*this), m) ;
	
	(*v) = ret.array() ;
	r = ret.numRows() ;
	c = ret.numCols() ;
	
	return *this ;
}


FunctionMatrix &FunctionMatrix::operator +=(const FunctionMatrix &m)
{
	(*v) += (m.array()) ;
	return *this ;
}

FunctionMatrix &FunctionMatrix::operator +=(const Matrix &m)
{
	for(size_t i = 0 ; i < numRows() ; i++)
	{
			for(size_t j = 0 ; j < numCols() ; j++)
		{
			(*this)[i][j] += m[i][j] ;
		}
	}
	return *this ;
}

FunctionMatrix FunctionMatrix::operator +(const  FunctionMatrix &m) const
{
	FunctionMatrix ret(*this) ;
	ret += m ; ;
	return ret ;
}

FunctionMatrix &FunctionMatrix::operator -=(const FunctionMatrix &m)
{
	(*v) -= (m.array()) ;
	return *this ;
}

FunctionMatrix FunctionMatrix::operator -(const  FunctionMatrix &m) const
{
	FunctionMatrix ret(*this) ;
	ret -= m ; ;
	return ret ;
}

FunctionMatrix FunctionMatrix::d(const Variable vv) const
{
	FunctionMatrix ret(*this) ;
	
	for(size_t i = 0 ; i < ret.array().size() ; i++)
	{
		if( (*v)[i].isDifferentiable() )
			ret.array()[i]  = (*v)[i].d(vv) ;
		else
			std::cout << "\\o/" << std::flush ;
	}
	return ret ;
}

FunctionMatrix::FunctionMatrix(const FunctionMatrix& m) : r(m.numRows()), c( m.numCols())
{
	v = new std::valarray<Function>(m.size()) ;
	(*v) = m.array() ;
}

FunctionMatrix &FunctionMatrix::operator =(const FunctionMatrix &m)
{
	delete v ;
	v = new std::valarray<Function>(m.size()) ;
	(*v) = m.array() ;
	r = m.numRows() ;
	c = m.numCols() ;
	return *this ;
}

FMtM::operator const FunctionMatrix() const
{
	return fm_matrix_multiply(first, second) ;
}

MtFM::operator const FunctionMatrix() const
{
	return mf_matrix_multiply(first, second) ;
}

FMtFM::operator const FunctionMatrix() const
{
	return ff_matrix_multiply(first, second) ;
}

FMtV::operator const FunctionVector()
{
	return fmatrix_vector_multiply(m,v) ;
}

FMtFV::operator const FunctionVector()
{
	return (m*v) ;
}


FMtMtFM::operator const FunctionMatrix() const
{
	return ff_matrix_multiply(fm_matrix_multiply(first, second), third) ;
}


//clever 2x2 FunctionMatrix inversion. Thanks the numerical cookbook :)
FunctionMatrix inverse2x2FunctionMatrix(const FunctionMatrix s)
{
	
	if(s[0][0].byteCode.size() == 0 || s[1][1].byteCode.size() == 0)
	{
		Matrix swap(2,2) ; swap[0][0] = 0 ; swap[0][1] = 1 ; swap[1][0] = 1 ; swap[1][1] = 0 ;
		FunctionMatrix s_ ( s*swap );
		
		FunctionMatrix ret(2,2) ;
		Function r1 = Function("1")/s_[0][0] ;
		Function r2 = s_[1][0] * r1 ;
		Function r3 = r1* s_[0][1] ;
		Function r4 = s_[1][0] * r3 ;
		Function r5 = r4 - s_[1][1] ;
		Function r6 =  Function("1")/r5 ;
		ret[0][1] = r3*r6 ;
		ret[1][0] = r6*r2 ;
		Function r7 = r3*ret[1][0] ;
		ret[0][0] = r1-r7 ;
		ret[1][1] = r6*-1 ;
		return swap*ret ;
	}
	
	FunctionMatrix ret(2,2) ;
	Function r1 = Function("1")/s[0][0] ;
	Function r2 = s[1][0] * r1 ;
	Function r3 = r1* s[0][1] ;
	Function r4 = s[1][0] * r3 ;
	Function r5 = r4 - s[1][1] ;
	Function r6 = Function("1")/r5 ;
	ret[0][1] = r3*r6 ;
	ret[1][0] = r6*r2 ;
	Function r7 = r3*ret[1][0] ;
	ret[0][0] = r1-r7 ;
	ret[1][1] = r6*-1 ;
	
	return ret ;
}

FunctionMatrix inverse3x3FunctionMatrix(const FunctionMatrix m)
{
	FunctionMatrix m_(m) ;
	std::vector< Matrix > swap ;
	
// 	std::cout << m[0][0] << "   " << m[0][1] << "   " << m[0][2] << std::endl ;
// 	std::cout << m[1][0] << "   " << m[1][1] << "   " << m[1][2] << std::endl ;
// 	std::cout << m[2][0] << "   " << m[2][1] << "   " << m[2][2] << std::endl ;
	
	while(m_[0][0].byteCode.size() == 0 || m_[1][1].byteCode.size() == 0 || m_[2][2].byteCode.size() == 0)
	{
		
		for(size_t i = 0 ; i < 3 ; i++)
		{
			if(m_[i][i].byteCode.size() == 0)
			{
				for(size_t j = 0 ; j < 3 ; j++)
				{
					if(!m_[j][i].byteCode.size() == 0 )
					{
						Matrix temp(swapMatrix(i,j, 3)) ;
						swap.push_back(temp) ;
						m_ = m_* temp ;
					}
				}
			}
		}
	}
	
	FunctionMatrix ret(3,3) ;
	Function r11 = m_[1][1]*m_[2][2]-m_[1][2]*m_[2][1] ;
	Function r21 = m_[1][2]*m_[2][0]-m_[1][0]*m_[2][2] ;
	Function r31 = m_[1][0]*m_[2][1]-m_[1][1]*m_[2][0] ;
	Function det = m_[0][0]*(r11) + m_[0][1]*(r21) + m_[0][2]*(r31) ;
	
	ret[0][0] = r11 ; 
	ret[0][1] = m_[0][2]*m_[2][1] - m_[0][1]*m_[2][2] ;
	ret[0][2] = m_[0][1]*m_[1][2] - m_[0][2]*m_[1][1] ;
	
	ret[1][0] = r21 ; 
	ret[1][1] = m_[0][0]*m_[2][2] - m_[0][2]*m_[2][0] ;
	ret[1][2] = m_[0][2]*m_[1][0] - m_[0][0]*m_[1][2] ;
	
	ret[2][0] = r31 ; 
	ret[2][1] = m_[0][1]*m_[2][0] - m_[0][0]*m_[2][1] ;
	ret[2][2] = m_[1][1]*m_[0][0] - m_[0][1]*m_[1][0] ;
	
	for(std::vector< Matrix >::reverse_iterator i = swap.rbegin() ; i != swap.rend()  ; ++i )
	{
		ret = (*i)*ret ;
	}
	
// 	std::cout << ret[0][0] << "   " << ret[0][1] << "   " << ret[0][2] << std::endl ;
// 	std::cout << ret[1][0] << "   " << ret[1][1] << "   " << ret[1][2] << std::endl ;
// 	std::cout << ret[2][0] << "   " << ret[2][1] << "   " << ret[2][2] << std::endl ;
	
	return ret/det ;
}

inline FMtFV operator*(const FunctionMatrix& mm, const FunctionVector& v)
{
	return FMtFV(mm, v) ;
} ;

inline FMtFM operator*(const FunctionMatrix& mm, const FunctionMatrix& mmm)
{
	return FMtFM(mm, mmm) ;
} ;

inline FMtMtFM operator*(const FMtM& mm, const FunctionMatrix& mmm)
{
	return FMtMtFM(mm.first, mm.second, mmm) ;
} ;


inline FMtV operator*(const FunctionMatrix& mm, const Vector & v)
{
	return FMtV(mm, v) ;
} ;

inline FMtM operator*(const FunctionMatrix& mm, const Matrix& mmm)
{
	return FMtM(mm, mmm) ;
} ;

inline MtFM operator*(const Matrix& mm, const FunctionMatrix& mmm)
{
	return MtFM(mm, mmm) ;
} ;

inline const FunctionMatrix ff_matrix_multiply(const FunctionMatrix &m0, const FunctionMatrix &m1 )
{
	assert(m0.numCols() == m1.numRows()) ;
	
	FunctionMatrix ret(m0.numRows(), m1.numCols()) ;
	
	for(size_t i = 0 ; i < m0.numRows() ; i++)
	{
		for(size_t j = 0 ; j < m1.numCols() ; j++)
		{
			const Cslice_iter<Function>& ri = m0[i] ;
			const Cslice_iter<Function>& cj = m1.column(j) ;
			ret[i][j] = std::inner_product(&ri[0], &ri[m0.numCols()], cj, Function()) ;
		}
	}
	return ret ;
}

inline const FunctionMatrix fm_matrix_multiply(const FunctionMatrix &m0, const Matrix &m1 )
{
	assert(m0.numCols() == m1.numRows()) ;
	
	FunctionMatrix ret(m0.numRows(), m1.numCols()) ;
	
	for(size_t i = 0 ; i < m0.numRows() ; i++)
	{
		for(size_t j = 0 ; j < m1.numCols() ; j++)
		{
			const Cslice_iter<Function>& ri = m0[i] ;
			const Cslice_iter<double>& cj = m1.column(j) ;
			ret[i][j] = std::inner_product(&ri[0], &ri[m0.numCols()], cj, Function() ) ;
		}
	}
	return ret ;
}

inline const FunctionMatrix mf_matrix_multiply(const Matrix &m0, const FunctionMatrix &m1 )
{
	assert(m0.numCols() == m1.numRows()) ;
	
	FunctionMatrix ret(m0.numRows(), m1.numCols()) ;
	
	for(size_t i = 0 ; i < m0.numRows() ; i++)
	{
		for(size_t j = 0 ; j < m1.numCols() ; j++)
		{
			const Cslice_iter<double>& ri = m0[i] ;
			const Cslice_iter<Function>& cj = m1.column(j) ;
			ret[i][j] = std::inner_product(&cj[0], &cj[m0.numCols()], ri, Function() ) ;
		}
	}
	return ret ;
}


inline const FunctionVector matrix_fvector_multiply(const Matrix &m, const FunctionVector &v )
{
	assert(m.numRows() == v.size()) ;
	
	FunctionVector ret(v.size()) ;
	
	for(size_t i = 0 ; i < m.numRows() ; i++)
	{
		
		const Cslice_iter<double>& ri = m[i] ;
		ret[i] = std::inner_product(&v[0], &v[v.size()],ri, Function() ) ;
	}
	return ret ;
}

inline const FunctionVector fmatrix_vector_multiply(const FunctionMatrix &m, const Vector &v )
{
	assert(m.numRows() == v.size()) ;
	
	FunctionVector ret(v.size()) ;
	
	for(size_t i = 0 ; i < m.numRows() ; i++)
	{
		
		const Cslice_iter<Function>& ri = m[i] ;
		ret[i] = std::inner_product(ri, ri.end(), &v[0], Function() ) ;
	}
	return ret ;
}

inline const FunctionVector operator*(const FunctionVector &v , const FunctionMatrix &m )
{
	assert(m.numCols() == v.size()) ;
	
	FunctionVector ret(v.size()) ;
	
	for(size_t i = 0 ; i < m.numCols() ; i++)
	{
		
		const Cslice_iter<Function>& ri = m.column(i) ;
		ret[i] = std::inner_product(ri, ri.end(), &v[0] , Function()) ;
	}
	return ret ;
}

inline const FunctionVector operator*(const Vector &v , const FunctionMatrix &m )
{
	assert(m.numCols() == v.size()) ;
	
	FunctionVector ret(v.size()) ;
	
	for(size_t i = 0 ; i < m.numCols() ; i++)
	{
		
		const Cslice_iter<Function>& ri = m.column(i) ;
		ret[i] = std::inner_product(ri, ri.end(), &v[0] , Function()) ;
	}
	return ret ;
}

GtFMtG GtFM::operator*(const Gradient & f) const
{
	return GtFMtG(this->first, this->second, f) ;
}

GtFM operator *(const Gradient & g, const FunctionMatrix & m)
{
	return GtFM(g, m) ;
}

} ;