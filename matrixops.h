// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __MATRIXOPS_H__
#define __MATRIXOPS_H__

#include <valarray>
#include "sliceiters.h"
#include <iostream>
#include <vector>

#include <assert.h>

typedef std::valarray<double> Vector ;


namespace Mu
{

struct MtM ;

class Matrix 
{
	Vector *v;
	size_t r, c;

public:
	
	/** Construct an M x N matrix.
	 * 
	 * All Elements are initialised to 0.
	 * 
	 * @param x number of rows
	 * @param y number of columns
	 */
	Matrix(size_t x, size_t y);
	
		/** Construct a 2x2 matrix
	 */
	Matrix()
	{
		r = 2 ;
		c = 2 ;
		v = new Vector(0., 4) ;
	}
	
	/** Construct a square submatrix from a given source matrix.
	 * 
	 * @param rl size of the submatrix
	 * @param k staring row in the source matrix
	 * @param l staring column in the source matrix
	 * @param m source matrix
	 */
	Matrix(size_t rl, size_t k, size_t l, const Matrix & m) ;
	
	
	/** construct a matrix from a given source matrix without specified row and column.
	 * 
	 * @param i 
	 * @param j 
	 * @param m 
	 */
	Matrix(size_t i, size_t j, const Matrix & m) ;
	
	
	Matrix(const Matrix&) ;
	virtual ~Matrix() {delete v; }
	
	Matrix &operator=(const Matrix &m);
	
	size_t size() const {return r*c ;}
	size_t numCols() const {return c ;}
	size_t numRows() const {return r ;}
	
	Slice_iter< double > column(size_t i )
	{
		return Slice_iter< double >(v, std::slice(i, r, c)) ;
	}
	Cslice_iter< double > column(size_t i ) const 
	{
		return Cslice_iter< double >(v, std::slice(i, r, c)) ;
	}
	
	Slice_iter< double > row(size_t i )
	{
		return Slice_iter< double >(v, std::slice(i*c, c, 1)) ;
	}
	Cslice_iter< double > row(size_t i ) const
	{
		return Cslice_iter< double >(v, std::slice(i*c, c, 1)) ;
	}
	
	Matrix transpose() const ;

	double& operator()(size_t x, size_t y) ;
	double operator()(size_t x, size_t y) const;
	
	Slice_iter< double > operator()(size_t i) {return row(i) ;}
	Cslice_iter< double > operator()(size_t i) const {return row(i) ;}
	
	Slice_iter< double > operator[](size_t i) {return row(i) ;}
	Cslice_iter< double > operator[](size_t i) const {return row(i) ;}
	
	Matrix& operator *=(double);
	Matrix operator *(double) const;
	Matrix operator /(double) const;
	Matrix& operator /=(double) ;
	Matrix& operator *=(const Matrix &m);
// 	const Matrix& operator *(const Matrix &m) const ;
	Matrix& operator +=(const Matrix &m) ;
	Matrix operator +(const Matrix &m) const;
	Matrix& operator -=(const Matrix &m) ;
	Matrix operator -(const Matrix &m) const;
	Matrix& operator =(const MtM& m) ;
	
	bool operator ==(const Matrix &m) ;
	bool operator !=(const Matrix &m) ;
	
	Vector &array() {return *v ;}
	Vector array() const {return *v ;}
} ;


struct MtV
{
	const Matrix &m;
	const Vector &v;
	
	MtV(const Matrix &mm, const Vector &vv) : m(mm), v(vv) { }
	
	operator const Vector();
} ;


struct MtM
{
	const Matrix &first;
	const Matrix &second;
	
	MtM(const Matrix &mm, const Matrix &mmm) : first(mm), second(mmm) { }
	
	operator const Matrix() const;
} ;

}

inline Mu::MtV operator*(const Mu::Matrix& mm, const Vector& v)
{
	return Mu::MtV(mm, v) ;
} ;

inline Mu::MtM operator*(const Mu::Matrix& mm, const Mu::Matrix& mmm)
{
	return Mu::MtM(mm, mmm) ;
} ;

inline const Mu::Matrix matrix_multiply(const Mu::Matrix &m0, const Mu::Matrix &m1 )
{
	assert(m0.numCols() == m1.numRows()) ;
	
	Mu::Matrix ret(m0.numRows(), m1.numCols()) ;
	
	for(size_t i = 0 ; i < m0.numRows() ; i++)
	{
		for(size_t j = 0 ; j < m1.numCols() ; j++)
		{
			const Mu::Cslice_iter<double>& ri = m0.row(i) ;
			const Mu::Cslice_iter<double>& cj = m1.column(j) ;
			ret[i][j] = std::inner_product(&ri[0], &ri[m0.numCols()], cj, (double)(0) ) ;
		}
	}
	return ret ;
}

inline const Vector matrix_vector_multiply(const Mu::Matrix &m, const Vector &v )
{
	assert(m.numRows() == v.size()) ;
	
	Vector ret(v.size()) ;
	
	for(size_t i = 0 ; i < m.numRows() ; i++)
	{

		const Mu::Cslice_iter<double>& ri = m.row(i) ;
		ret[i] = std::inner_product(ri, ri.end(), &v[0], (double)(0) ) ;
	}
	return ret ;
}

inline const Vector operator*(const Vector &v , const Mu::Matrix &m )
{
	assert(m.numCols() == v.size()) ;
	
	Vector ret(v.size()) ;
	
	for(size_t i = 0 ; i < m.numCols() ; i++)
	{

		const Mu::Cslice_iter<double>& ri = m.column(i) ;
		ret[i] = std::inner_product(ri, ri.end(), &v[0], (double)(0) ) ;
	}
	return ret ;
}

Mu::Matrix swapMatrix(size_t l0, size_t l1, size_t dim = 2) ;
void swapLines(size_t l0, size_t l1, Mu::Matrix * m) ;

//clever 2x2 Matrix inversion. Thanks the numerical cookbook :)
Mu::Matrix inverse2x2Matrix(const Mu::Matrix &s) ;
void invert2x2Matrix(Mu::Matrix &s) ;

Mu::Matrix inverse3x3Matrix(const Mu::Matrix &m) ;
void invert3x3Matrix(Mu::Matrix & m) ;

Mu::Matrix inverse4x4Matrix(const Mu::Matrix &s) ;

double det(const Mu::Matrix &s) ;

Mu::Matrix exp(const Mu::Matrix&, size_t order = 20) ;

Mu::Matrix log(const Mu::Matrix&, size_t order) ;



# endif  // __MATRIXOPS_H__
