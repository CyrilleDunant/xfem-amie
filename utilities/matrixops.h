// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __MATRIXOPS_H__
#define __MATRIXOPS_H__

#include <valarray>
#include <iostream>
#include "sliceiters.h"
#include <vector>
#include <numeric>


#include <assert.h>


#ifdef HAVE_SSE3
#include <emmintrin.h>

typedef union 
{
	__m128d vec ;
	double val[2] ;
} vecdouble;
#endif

typedef std::valarray<double> Vector ;

namespace Mu
{

struct MtM ;
struct MtV ;
struct MtMtM ;

/** \brief Full matrix class*/
class Matrix 
{
	Vector *v;
	size_t r, c;
	bool cleanup ;
public:
	
	/** \brief Construct an M x N matrix.
	 * 
	 * All Elements are initialised to 0.
	 * 
	 * @param x number of rows
	 * @param y number of columns
	 */
	Matrix(size_t x, size_t y);

	/** \brief Special case: construct a matrix with external memory storage.*/
	Matrix(size_t x, size_t y, Vector * cache);
	
	/** Construct a 0x0 matrix
	 */
	Matrix()
	{
		cleanup = false ;
		r = 0 ;
		c = 0 ;
		v = NULL ;
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
	
	/** \brief copy-constructor */
	Matrix(const Matrix&) ;

	/** \brief constructor from a vector<vector<double> > */
	Matrix(const std::vector<std::vector<double> > & vv) ;

	virtual ~Matrix() { if(cleanup)delete v; }

	/** \brief Copy-constructor.
	* If the matrix is not yet initialised, initialise to the source matrix values. Else do 
	* a simple copy. Because this involves copying valarrays, if the matrices have
	* different dimensions, the result will be undefined.
	*/
	Matrix &operator=(const Matrix &m);
	
	/** \brief Return the number of elements in the matrix*/
	size_t size() const {return r*c ;}
	
	/** \brief Return the number of columns*/
	size_t numCols() const {return c ;}
	
	/** \brief Return the number of rows*/
	size_t numRows() const {return r ;}
	
	/** \brief Return a column iterator*/
	Slice_iter< double > column(size_t i )
	{
		return Slice_iter< double >(v, std::slice(i, r, c)) ;
	}
	
	/** \brief return a column const iterator*/
	Cslice_iter< double > column(size_t i ) const 
	{
		return Cslice_iter< double >(v, std::slice(i, r, c)) ;
	}
	
	/** \brief return a row iterator*/
	Slice_iter< double > row(size_t i )
	{
		return Slice_iter< double >(v, std::slice(i*c, c, 1)) ;
	}
	
	/** \brief return a const row iterator*/
	Cslice_iter< double > row(size_t i ) const
	{
		return Cslice_iter< double >(v, std::slice(i*c, c, 1)) ;
	}
	
	/** \brief construct a new transpose Matrix*/
	Matrix transpose() const ;

	/** \brief discard all values and resize this matrix*/
	void resize(size_t x, size_t y)
	{
		r = x ;
		c = y ;
		double zero(0);
		if(v == NULL)
			v = new Vector(zero, x*y) ;
		else
			v->resize(x*y, zero) ;

		cleanup = true ;
	}

	/** \brief check if the matrix has initialised values*/
	bool isNull() const
	{
		return v == NULL ;
	}

	/** \brief FORTRAN-style accessors*/
	double& operator()(size_t x, size_t y) ;
	
	/** \brief FORTRAN-style accessors*/
	double operator()(size_t x, size_t y) const;
	
	/** \brief FORTRAN-style accessor: return row*/
	Slice_iter< double > operator()(size_t i) {return row(i) ;}
	
	/** \brief FORTRAN-style accessor: return row*/
	Cslice_iter< double > operator()(size_t i) const {return row(i) ;}
	
	/** \brief c++-style accessor: return row*/
	Slice_iter< double > operator[](size_t i) {return row(i) ;}
	
	/** \brief c++-style accessor: return row*/
	Cslice_iter< double > operator[](size_t i) const {return row(i) ;}
	
	Matrix& operator *=(double);
	Matrix operator *(double) const;
	Matrix operator /(double) const;
	Matrix& operator /=(double) ;
	Matrix& operator *=(const Matrix &m);
// 	const Matrix& operator *(const Matrix &m) const ;
	Matrix& operator +=(const Matrix &m) ;
	Matrix& operator +=(const MtM &m) ;
	Matrix& operator +=(const MtMtM &m) ;
	Matrix operator +(const Matrix &m) const;
	Matrix& operator -=(const Matrix &m) ;
	Matrix operator -(const Matrix &m) const;
	Matrix& operator =(const MtM& m) ;
	Matrix& operator =(const MtMtM& m) ;
	
	bool operator ==(const Matrix &m) ;
	bool operator !=(const Matrix &m) ;

	/** \brief Return the array of values, stored in row-major fashion*/
	Vector &array() {return *v ;}
	
	/** \brief Return the array of values, stored in row-major fashion*/
	Vector array() const {return *v ;}
	
	void print() const ;

	

} ;

/** \brief perform the inner product of two vector of doubles
* @param v0 Pointer to the start of the first vector
* @param v1 pointer to the start of the second vector
* @param size size of the vectors
*/
double parallel_inner_product(const double * v0, const double * v1, int size) ;

/** \brief Structure used for Matrix-Vector multiplication, minimising temporaries*/
struct MtV
{
	const Matrix &m;
	const Vector &v;
	
	/** \brief constructor, initialise the references necessary for the computation*/
	MtV(const Matrix &mm, const Vector &vv) : m(mm), v(vv) { }
	
	/** \brief cast operator, compute the operation*/
	operator const Vector();
} ;

/** \brief Structure used for Matrix-Matrix-Matrix multiplication, minimising the use of temporares*/
struct MtMtM
{
	const Matrix &first;
	const Matrix &second;
	const Matrix &third;
	
	/** \brief constructor, initialise the references necessary for the computation*/
	MtMtM(const Matrix &mm, const Matrix &mmm, const Matrix &mmmm) : first(mm), second(mmm), third(mmmm) { }
	
	/** \brief cast operator, compute the operation*/
	operator const Matrix() const;
} ;

/** \brief Structure used for Matrix-Matrix multiplication, minimising temporaries*/
struct MtM
{
	const Matrix &first;
	const Matrix &second;
	
	/** \brief constructor, initialise the references necessary for the computation*/
	MtM(const Matrix &mm, const Matrix &mmm) : first(mm), second(mmm) { }
	
	/** \brief Construct a Matrix-Matrix-Matrix multiplication*/
	MtMtM operator *(const Matrix &mm) { return MtMtM(first, second, mm) ; }
	
	/** \brief cast operator, compute the operation*/
	operator const Matrix() const;
} ;

/** \brief Solve Ax = b system of equation, using a preconditionned Conjugate Gradient Algorithm*/
Vector solveSystem(const Matrix & A, const Vector & b, Vector & x) ;

}

/** \brief Construct a Matrix-Vector multiplication object*/
inline Mu::MtV operator*(const Mu::Matrix& mm, const Vector& v)
{
	return Mu::MtV(mm, v) ;
}

/** \brief Construct a Matrix-Matrix multiplication object*/
inline Mu::MtM operator*(const Mu::Matrix& mm, const Mu::Matrix& mmm)
{
	return Mu::MtM(mm, mmm) ;
}

/** \brief Perform a Matrix-Matrix multiplication */
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

/** \brief Perform a Matrix-Matrix-Matrix multiplication*/
inline const Mu::Matrix matrix_matrix_matrix_multiply(const Mu::Matrix &m0, const Mu::Matrix &m1, const Mu::Matrix &m2 )
{
	Mu::Matrix ret(m0.numRows(), m2.numCols()) ;
	
	for(size_t i = 0 ; i < m0.numRows() ; i++)
	{
		for(size_t j = 0 ; j < m1.numCols() ; j++)
		{
			const Mu::Cslice_iter<double>& ri = m0.row(i) ;
			const Mu::Cslice_iter<double>& cj = m1.column(j) ;
			double r_ij = std::inner_product(&ri[0], &ri[m0.numCols()], cj, (double)(0) ) ;
			for(size_t k = 0 ; k < m2.numCols() ; k++)
				ret[i][k] += r_ij*m2[j][k] ;
		}
	}
	
	return ret ;
}

/** \brief Perform a Matrix-Matrix-Matrix multiplication, and assign to already initialised matrix passed as further argument*/
inline void matrix_matrix_matrix_multiply_and_assign(const Mu::Matrix &m0, const Mu::Matrix &m1, const Mu::Matrix &m2, Mu::Matrix & ret )
{
	ret.array() = 0 ;
	for(size_t i = 0 ; i < m0.numRows() ; i++)
	{
		for(size_t j = 0 ; j < m1.numCols() ; j++)
		{
			const Mu::Cslice_iter<double>& ri = m0.row(i) ;
			const Mu::Cslice_iter<double>& cj = m1.column(j) ;
			double r_ij = std::inner_product(&ri[0], &ri[m0.numCols()], cj, (double)(0) ) ;
			for(size_t k = 0 ; k < m2.numCols() ; k++)
				ret[i][k] += r_ij*m2[j][k] ;
		}
	}
	
}

/** \brief Perform a Matrix-Matrix-Matrix multiplication, multiply with a factor and add to already initialised matrix passed as further argument*/
inline void matrix_matrix_matrix_multiply_and_add(const Mu::Matrix &m0, const Mu::Matrix &m1, const Mu::Matrix &m2, double f,  Mu::Matrix & ret )
{
	for(size_t i = 0 ; i < m0.numRows() ; i++)
	{
		for(size_t j = 0 ; j < m1.numCols() ; j++)
		{
			const Mu::Cslice_iter<double>& ri = m0.row(i) ;
			const Mu::Cslice_iter<double>& cj = m1.column(j) ;
			double r_ij = std::inner_product(&ri[0], &ri[m0.numCols()], cj, (double)(0) ) ;
			for(size_t k = 0 ; k < m2.numCols() ; k++)
				ret[i][k] += r_ij*m2[j][k]*f ;
		}
	}
}

/** \brief Perform a Matrix-Matrix multiplication, and assign to already initialised matrix passed as further argument*/
inline void matrix_multiply_and_assign(const Mu::Matrix &m0, const Mu::Matrix &m1,  Mu::Matrix &ret)
{
	assert(m0.numCols() == m1.numRows()) ;
// 		std::cout << ret.numRows() << ", " << ret.numCols() << ", " << m0.numRows() << ", " <<m0.numCols() << ", " << m1.numRows() << ", "  <<m1.numCols() << std::endl ;
	for(size_t i = 0 ; i < m0.numRows() ; i++)
	{
		for(size_t j = 0 ; j < m1.numCols() ; j++)
		{
			const Mu::Cslice_iter<double>& ri = m0.row(i) ;
			const Mu::Cslice_iter<double>& cj = m1.column(j) ;
			ret[i][j] = std::inner_product(&ri[0], &ri[m0.numCols()], cj, (double)(0) ) ;
		}
	}
}

/** \brief Perform a Matrix-Matrix multiplication, and add to already initialised matrix passed as further argument*/
inline void matrix_multiply_and_add(const Mu::Matrix &m0, const Mu::Matrix &m1,  Mu::Matrix &ret)
{
	assert(m0.numCols() == m1.numRows()) ;
// 	std::cout << ret.numRows() << ", " << ret.numCols() << ", " << m0.numRows() << ", " <<m0.numCols() << ", " << m1.numRows() << ", "  <<m1.numCols() << std::endl ;
	for(size_t i = 0 ; i < m0.numRows() ; i++)
	{
		for(size_t j = 0 ; j < m1.numCols() ; j++)
		{
			const Mu::Cslice_iter<double>& ri = m0.row(i) ;
			const Mu::Cslice_iter<double>& cj = m1.column(j) ;
			ret[i][j] += std::inner_product(&ri[0], &ri[m0.numCols()], cj, (double)(0) ) ;
		}
	}
}

/** \brief Perform a Matrix-Vector multiplication*/
inline const Vector matrix_vector_multiply(const Mu::Matrix &m, const Vector &v )
{
	assert(m.numRows() == v.size()) ;
	
	Vector ret(m.numRows()) ;
	
	for(size_t i = 0 ; i < m.numRows() ; i++)
	{

		const Mu::Cslice_iter<double>& ri = m.row(i) ;
		ret[i] = std::inner_product(ri, ri.end(), &v[0], (double)(0) ) ;
	}
	return ret ;
}

/** \brief Perform a Vector-Matrix multiplication*/
inline Vector operator*(const Vector &v , const Mu::Matrix &m )
{
	assert(m.numCols() == v.size()) ;
	
	Vector ret(m.numCols()) ;
	
	for(size_t i = 0 ; i < m.numCols() ; i++)
	{
		const Mu::Cslice_iter<double>& ri = m.column(i) ;
		ret[i] = std::inner_product(ri, ri.end(), &v[0], (double)(0) ) ;
	}
	return ret ;

}


/** \brief Return a matrix, the multiplication by which is equivalent to swapping two lines
*
* @param l0 index of the first line
* @param l1 index of the second line
* @param size of the matrix to produce
*/
Mu::Matrix swapMatrix(size_t l0, size_t l1, size_t dim = 2) ;
void swapLines(size_t l0, size_t l1, Mu::Matrix * m) ;

//clever 2x2 Matrix inversion. Thanks the numerical cookbook :)

/** \brief 2x2 matrix inversion*/
Mu::Matrix inverse2x2Matrix(const Mu::Matrix &s) ;

/** \brief in-place 2x2 matrix inversion*/
void invert2x2Matrix(Mu::Matrix &s) ;

/** \brief 3x3 matrix inversion*/
Mu::Matrix inverse3x3Matrix(const Mu::Matrix &m) ;

/** \brief in-place 3x3 matrix inversion*/
void invert3x3Matrix(Mu::Matrix & m) ;

/** \brief 4x4 matrix inversion*/
Mu::Matrix inverse4x4Matrix(const Mu::Matrix &s) ;

/** \brief 6x6 matrix inversion*/
Mu::Matrix inverse6x6Matrix(const Mu::Matrix &s) ;

/** \brief in-place 6x6 matrix inversion*/
void invert6x6Matrix(Mu::Matrix & m) ;

/** \brief compute the determinant of a Matrix, slow*/
double det(const Mu::Matrix &s) ;

/** \brief Compute the exp of a Matrix, using series expansion. May not converge*/
Mu::Matrix exp(const Mu::Matrix&, size_t order = 20) ;

/** \brief Compute the log of a Matrix using series Expansion, May not converge*/
Mu::Matrix log(const Mu::Matrix&, size_t order) ;



# endif  // __MATRIXOPS_H__

