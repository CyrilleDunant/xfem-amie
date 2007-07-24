//
// C++ Interface: vm_function_matrix
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef VM_FUNCTION_MATRIX_H
#define VM_FUNCTION_MATRIX_H

#include <valarray>

#include "../sliceiters.h"
#include "../matrixops.h"
#include "vm_function_base.h"

namespace Mu
{


struct FMtFM ;
struct FMtM ;

typedef std::valarray<Function> FunctionVector ;

class FunctionMatrix 
{
	std::valarray< Function > *v;
	size_t r, c;
	
public:
	FunctionMatrix(size_t x, size_t y);
	FunctionMatrix()
	{
		r = 2 ;
		c = 2 ;
		v = new std::valarray< Function >(4) ;
	}
	
	FunctionMatrix(const FunctionMatrix&) ;
	
	virtual ~FunctionMatrix() {delete v; }
	
	FunctionMatrix &operator=(const FunctionMatrix &m);
	
	size_t size() const {return r*c ;}
	size_t numCols() const {return c ;}
	size_t numRows() const {return r ;}
	
	Slice_iter< Function > column(size_t i )
	{
		return Slice_iter< Function >(v, std::slice(i, r, c)) ;
	}
	Cslice_iter< Function > column(size_t i ) const 
	{
		return Cslice_iter< Function >(v, std::slice(i, r, c)) ;
	}
	
	Slice_iter< Function > row(size_t i )
	{
		return Slice_iter< Function >(v, std::slice(i*c, c, 1)) ;
	}
	Cslice_iter< Function > row(size_t i ) const
	{
		return Cslice_iter< Function >(v, std::slice(i*c, c, 1)) ;
	}
	
	FunctionMatrix transpose() const ;
	
	FunctionMatrix d(const Variable v) const ;
	
	Function& operator()(size_t x, size_t y) ;
	Function operator()(size_t x, size_t y) const;
	
	Slice_iter< Function > operator()(size_t i) {return row(i) ;}
	Cslice_iter< Function > operator()(size_t i) const {return row(i) ;}
	
	Slice_iter< Function > operator[](size_t i) {return row(i) ;}
	Cslice_iter< Function > operator[](size_t i) const {return row(i) ;}
	
	FunctionMatrix& operator =(const FMtFM& m) ;
	FunctionMatrix& operator =(const FMtM& m) ;
	
	FunctionMatrix operator *(const Function) const;
	FunctionMatrix operator /(const Function) const;
	FunctionMatrix& operator /=(const Function) ;
	FunctionMatrix& operator *=(const Function &m);
// 	const Matrix& operator *(const Matrix &m) const ;
	FunctionMatrix& operator *=(const FunctionMatrix &m);
	FunctionMatrix& operator +=(const FunctionMatrix &m) ;
	FunctionMatrix operator +(const FunctionMatrix &m) const;
	FunctionMatrix& operator -=(const FunctionMatrix &m) ;
	FunctionMatrix operator -(const FunctionMatrix &m) const;
	
	FunctionMatrix& operator *=(const double);
	FunctionMatrix operator *(const double) const;
	FunctionMatrix operator /(const double) const;
	FunctionMatrix& operator /=(const double) ;
// 	const Matrix& operator *(const Matrix &m) const ;
	FunctionMatrix& operator +=(const Matrix &m) ;
	FunctionMatrix operator +(const Matrix &m) const;
	FunctionMatrix& operator -=(const Matrix &m) ;
	FunctionMatrix operator -(const Matrix &m) const;;
	
// 	bool operator ==(const FunctionMatrix &m) ;
// 	bool operator !=(const FunctionMatrix &m) ;
	
	std::valarray< Function > &array() {return *v ;}
	std::valarray< Function > array() const {return *v ;}
} ;


struct FMtFV
{
	const FunctionMatrix &m;
	const FunctionVector &v;
	
	FMtFV(const FunctionMatrix &mm, const FunctionVector &vv) : m(mm), v(vv) { }
	
	operator const FunctionVector();
} ;

struct FMtV
{
	const FunctionMatrix &m;
	const Vector &v;
	
	FMtV(const FunctionMatrix &mm, const Vector &vv) : m(mm), v(vv) { }
	
	operator const FunctionVector();
} ;


struct FMtFM
{
	const FunctionMatrix &first;
	const FunctionMatrix &second;
	
	FMtFM(const FunctionMatrix &mm, const FunctionMatrix &mmm) : first(mm), second(mmm) { }
	
	operator const FunctionMatrix() const;
} ;

struct MtFM
{
	const Matrix &first;
	const FunctionMatrix &second;
	
MtFM(const Matrix &mm, const FunctionMatrix &mmm) : first(mm), second(mmm) { }
	
	operator const FunctionMatrix() const;
} ;

struct FMtM
{
	const FunctionMatrix &first;
	const Matrix &second;
	
	FMtM(const FunctionMatrix &mm, const Matrix &mmm) : first(mm), second(mmm) { }
	
	operator const FunctionMatrix() const;
} ;

struct FMtMtFM
{
	const FunctionMatrix &first;
	const Matrix &second;
	const FunctionMatrix &third;
	
	FMtMtFM(const FunctionMatrix &mm, const Matrix &mmm, const FunctionMatrix &mmmm) : first(mm), second(mmm), third(mmmm) { }
	
	operator const FunctionMatrix() const;
} ;


} ;

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

//clever 2x2 Matrix inversion. Thanks the numerical cookbook :)
Mu::FunctionMatrix inverse2x2FunctionMatrix(const Mu::FunctionMatrix s) ;

Mu::FunctionMatrix inverse3x3FunctionMatrix(const Mu::FunctionMatrix m) ;





#endif

