//
// C++ Implementation: vm_function_matrix
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "vm_function_matrix.h"

using namespace Mu ;

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

FunctionMatrix FunctionMatrix::operator*(const Function d) const
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

FunctionMatrix FunctionMatrix::operator/(const Function d) const
{
	FunctionMatrix ret(*this) ;
	ret /=d ;
	return ret ;
}

FunctionMatrix &FunctionMatrix::operator/=(Function d)
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


FunctionMatrix &FunctionMatrix::operator +=(const FunctionMatrix &m)
{
	(*v) += (m.array()) ;
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

// bool FunctionMatrix::operator ==(const FunctionMatrix &m)
// {
// 	if(v->size() != m.array().size())
// 		return false ;
// 	
// 	for(size_t i = 0 ; i < v->size() ; i++)
// 		if((*v)[i] != m.array()[i])
// 			return false ;
// 	
// 	return true ;
// }
// 
// bool FunctionMatrix::operator !=(const FunctionMatrix &m)
// {
// 	if(v->size() != m.array().size())
// 		return true ;
// 	else
// 		for(size_t i = 0 ; i < v->size() ; i++)
// 			if((*v)[i] != m.array()[i])
// 				return true ;
// 	
// 	return false ;
// }

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
Mu::FunctionMatrix inverse2x2FunctionMatrix(const Mu::FunctionMatrix s)
{
	
	if(s[0][0].getByteCode().size() == 0 || s[1][1].getByteCode().size() == 0)
	{
		Mu::Matrix swap(2,2) ; swap[0][0] = 0 ; swap[0][1] = 1 ; swap[1][0] = 1 ; swap[1][1] = 0 ;
		Mu::FunctionMatrix s_ = s*swap ;
		
		Mu::FunctionMatrix ret(2,2) ;
		Mu::Function r1 = Mu::Function("1")/s_[0][0] ;
		Mu::Function r2 = s_[1][0] * r1 ;
		Mu::Function r3 = r1* s_[0][1] ;
		Mu::Function r4 = s_[1][0] * r3 ;
		Mu::Function r5 = r4 - s_[1][1] ;
		Mu::Function r6 =  Mu::Function("1")/r5 ;
		ret[0][1] = r3*r6 ;
		ret[1][0] = r6*r2 ;
		Mu::Function r7 = r3*ret[1][0] ;
		ret[0][0] = r1-r7 ;
		ret[1][1] = r6*-1 ;
		return swap*ret ;
	}
	
	Mu::FunctionMatrix ret(2,2) ;
	Mu::Function r1 = Mu::Function("1")/s[0][0] ;
	Mu::Function r2 = s[1][0] * r1 ;
	Mu::Function r3 = r1* s[0][1] ;
	Mu::Function r4 = s[1][0] * r3 ;
	Mu::Function r5 = r4 - s[1][1] ;
	Mu::Function r6 = Mu::Function("1")/r5 ;
	ret[0][1] = r3*r6 ;
	ret[1][0] = r6*r2 ;
	Mu::Function r7 = r3*ret[1][0] ;
	ret[0][0] = r1-r7 ;
	ret[1][1] = r6*-1 ;
	
	return ret ;
}

Mu::FunctionMatrix inverse3x3FunctionMatrix(const Mu::FunctionMatrix m)
{
	Mu::FunctionMatrix m_(m) ;
	std::vector< Mu::Matrix > swap ;
	
// 	std::cout << m[0][0] << "   " << m[0][1] << "   " << m[0][2] << std::endl ;
// 	std::cout << m[1][0] << "   " << m[1][1] << "   " << m[1][2] << std::endl ;
// 	std::cout << m[2][0] << "   " << m[2][1] << "   " << m[2][2] << std::endl ;
	
	while(m_[0][0].getByteCode().size() == 0 || m_[1][1].getByteCode().size() == 0 || m_[2][2].getByteCode().size() == 0)
	{
		
		for(size_t i = 0 ; i < 3 ; i++)
		{
			if(m_[i][i].getByteCode().size() == 0)
			{
				for(size_t j = 0 ; j < 3 ; j++)
				{
					if(!m_[j][i].getByteCode().size() == 0 )
					{
						Mu::Matrix temp = swapMatrix(i,j, 3) ;
						swap.push_back(temp) ;
						m_ = m_* temp ;
					}
				}
			}
		}
	}
	
	Mu::FunctionMatrix ret(3,3) ;
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
	
	for(std::vector< Mu::Matrix >::reverse_iterator i = swap.rbegin() ; i != swap.rend()  ; ++i )
	{
		ret = (*i)*ret ;
	}
	
// 	std::cout << ret[0][0] << "   " << ret[0][1] << "   " << ret[0][2] << std::endl ;
// 	std::cout << ret[1][0] << "   " << ret[1][1] << "   " << ret[1][2] << std::endl ;
// 	std::cout << ret[2][0] << "   " << ret[2][1] << "   " << ret[2][2] << std::endl ;
	
	return ret/det ;
}

inline Mu::FMtFV operator*(const Mu::FunctionMatrix& mm, const Mu::FunctionVector& v)
{
	return Mu::FMtFV(mm, v) ;
} ;

inline Mu::FMtFM operator*(const Mu::FunctionMatrix& mm, const Mu::FunctionMatrix& mmm)
{
	return Mu::FMtFM(mm, mmm) ;
} ;

inline Mu::FMtMtFM operator*(const Mu::FMtM& mm, const Mu::FunctionMatrix& mmm)
{
	return Mu::FMtMtFM(mm.first, mm.second, mmm) ;
} ;


inline Mu::FMtV operator*(const Mu::FunctionMatrix& mm, const Vector & v)
{
	return Mu::FMtV(mm, v) ;
} ;

inline Mu::FMtM operator*(const Mu::FunctionMatrix& mm, const Mu::Matrix& mmm)
{
	return Mu::FMtM(mm, mmm) ;
} ;

inline Mu::MtFM operator*(const Mu::Matrix& mm, const Mu::FunctionMatrix& mmm)
{
	return Mu::MtFM(mm, mmm) ;
} ;

inline const Mu::FunctionMatrix ff_matrix_multiply(const Mu::FunctionMatrix &m0, const Mu::FunctionMatrix &m1 )
{
	assert(m0.numCols() == m1.numRows()) ;
	
	Mu::FunctionMatrix ret(m0.numRows(), m1.numCols()) ;
	
	for(size_t i = 0 ; i < m0.numRows() ; i++)
	{
		for(size_t j = 0 ; j < m1.numCols() ; j++)
		{
			const Mu::Cslice_iter<Mu::Function>& ri = m0.row(i) ;
			const Mu::Cslice_iter<Mu::Function>& cj = m1.column(j) ;
			ret[i][j] = std::inner_product(&ri[0], &ri[m0.numCols()], cj, Mu::Function()) ;
		}
	}
	return ret ;
}

inline const Mu::FunctionMatrix fm_matrix_multiply(const Mu::FunctionMatrix &m0, const Mu::Matrix &m1 )
{
	assert(m0.numCols() == m1.numRows()) ;
	
	Mu::FunctionMatrix ret(m0.numRows(), m1.numCols()) ;
	
	for(size_t i = 0 ; i < m0.numRows() ; i++)
	{
		for(size_t j = 0 ; j < m1.numCols() ; j++)
		{
			const Mu::Cslice_iter<Mu::Function>& ri = m0.row(i) ;
			const Mu::Cslice_iter<double>& cj = m1.column(j) ;
			ret[i][j] = std::inner_product(&ri[0], &ri[m0.numCols()], cj, Mu::Function() ) ;
		}
	}
	return ret ;
}

inline const Mu::FunctionMatrix mf_matrix_multiply(const Mu::Matrix &m0, const Mu::FunctionMatrix &m1 )
{
	assert(m0.numCols() == m1.numRows()) ;
	
	Mu::FunctionMatrix ret(m0.numRows(), m1.numCols()) ;
	
	for(size_t i = 0 ; i < m0.numRows() ; i++)
	{
		for(size_t j = 0 ; j < m1.numCols() ; j++)
		{
			const Mu::Cslice_iter<double>& ri = m0.row(i) ;
			const Mu::Cslice_iter<Mu::Function>& cj = m1.column(j) ;
			ret[i][j] = std::inner_product(&cj[0], &cj[m0.numCols()], ri, Mu::Function() ) ;
		}
	}
	return ret ;
}


inline const Mu::FunctionVector matrix_fvector_multiply(const Mu::Matrix &m, const Mu::FunctionVector &v )
{
	assert(m.numRows() == v.size()) ;
	
	Mu::FunctionVector ret(v.size()) ;
	
	for(size_t i = 0 ; i < m.numRows() ; i++)
	{
		
		const Mu::Cslice_iter<double>& ri = m.row(i) ;
		ret[i] = std::inner_product(&v[0], &v[v.size()],ri, Mu::Function() ) ;
	}
	return ret ;
}

inline const Mu::FunctionVector fmatrix_vector_multiply(const Mu::FunctionMatrix &m, const Vector &v )
{
	assert(m.numRows() == v.size()) ;
	
	Mu::FunctionVector ret(v.size()) ;
	
	for(size_t i = 0 ; i < m.numRows() ; i++)
	{
		
		const Mu::Cslice_iter<Mu::Function>& ri = m.row(i) ;
		ret[i] = std::inner_product(ri, ri.end(), &v[0], Mu::Function() ) ;
	}
	return ret ;
}

inline const Mu::FunctionVector operator*(const Mu::FunctionVector &v , const Mu::FunctionMatrix &m )
{
	assert(m.numCols() == v.size()) ;
	
	Mu::FunctionVector ret(v.size()) ;
	
	for(size_t i = 0 ; i < m.numCols() ; i++)
	{
		
		const Mu::Cslice_iter<Mu::Function>& ri = m.column(i) ;
		ret[i] = std::inner_product(ri, ri.end(), &v[0] , Mu::Function()) ;
	}
	return ret ;
}

inline const Mu::FunctionVector operator*(const Vector &v , const Mu::FunctionMatrix &m )
{
	assert(m.numCols() == v.size()) ;
	
	Mu::FunctionVector ret(v.size()) ;
	
	for(size_t i = 0 ; i < m.numCols() ; i++)
	{
		
		const Mu::Cslice_iter<Mu::Function>& ri = m.column(i) ;
		ret[i] = std::inner_product(ri, ri.end(), &v[0] , Mu::Function()) ;
	}
	return ret ;
}

GtFMtG GtFM::operator*(const Mu::Gradient & f) const
{
	return GtFMtG(this->first, this->second, f) ;
}

Mu::GtFM operator *(const Mu::Gradient & g, const Mu::FunctionMatrix & m)
{
	return Mu::GtFM(g, m) ;
}

