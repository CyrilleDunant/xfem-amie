// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

// #define nullptr NULL


#ifndef __MATRIXOPS_H__
#define __MATRIXOPS_H__


#include <iostream>
#include <valarray>
#include "sliceiters.h"
#include <vector>
#include <numeric>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif


#include <assert.h>


// #ifdef HAVE_SSE3
// #include <emmintrin.h>
//
// typedef union
// {
// 	__m128d vec ;
// 	double val[2] ;
// } vecdouble;
// #endif

typedef std::valarray<double> Vector ;

namespace Amie
{

struct MtM ;
struct MtV ;
struct MtMtM ;

/** \brief Full matrix class*/
class Matrix final
{
    Vector *v;
    size_t r, c;
    bool noClean ;
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

    /** \brief Construct a matrix from the product of two vectors.*/
    Matrix(const Vector & v1, const Vector & v2);

    /** Construct a 0x0 matrix
     */
    Matrix()
    {
        noClean = false ;
        r = 0 ;
        c = 0 ;
        v = nullptr ;
    }

    double froebeniusNorm() const
    {
        if(isNull())
            return 0. ;

        double n = 0 ;
        for( size_t i = 0 ; i < v->size() ; i++)
        {
            n += (*v)[i]* (*v)[i] ;
        }
        return sqrt(n) ;
    }

    double squareFroebeniusNorm() const
    {
        if(isNull())
            return 0. ;
        double n = 0 ;
        for( size_t i = 0 ; i < v->size() ; i++)
        {
            n += (*v)[i]* (*v)[i] ;
        }
        return n ;
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

    virtual ~Matrix() {
        if(!noClean) delete v;
    }

    /** \brief Copy-constructor.
    * If the matrix is not yet initialised, initialise to the source matrix values. Else do
    * a simple copy. Because this involves copying valarrays, if the matrices have
    * different dimensions, the result will be undefined.
    */
    Matrix &operator=(const Matrix &m);
    Matrix& operator =(const MtM& m) ;
    Matrix& operator =(double) ;
    Matrix& operator =(const MtMtM& m) ;

    /** \brief Return the number of elements in the matrix*/
    size_t size() const {
        return r*c ;
    }

    /** \brief Return the number of columns*/
    size_t numCols() const {
        return c ;
    }

    /** \brief Return the number of rows*/
    size_t numRows() const {
        return r ;
    }

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

    /** \brief construct a new transpose Matrix*/
    Matrix transpose() const ;
    void transpose(Matrix & m) const ;

    /** \brief discard all values and resize this matrix*/
    void resize(size_t x, size_t y)
    {
        r = x ;
        c = y ;
        if(v == nullptr)
        {
            v = new Vector(0., x*y) ;
            noClean = false ;
        }
        else
            v->resize(x*y, 0.) ;
    }

    /** \brief check if the matrix has initialised values*/
    bool isNull() const
    {
        return v == nullptr ;
    }

    /** \brief c++-style accessor: return row*/
    Slice_iter< double > operator[](size_t i) {
        return Slice_iter< double >(v, std::slice(i*c, c, 1)) ;
    }

    /** \brief c++-style accessor: return row*/
    Cslice_iter< double > operator[](size_t i) const {
        return Cslice_iter< double >(v, std::slice(i*c, c, 1)) ;
    }

    /** \brief Return the array of values, stored in row-major fashion*/
    Vector &array() {
        return *v ;
    }

    /** \brief Return the array of values, stored in row-major fashion*/
    const Vector & array() const {
        return *v ;
    }

    void print() const ;



} ;

Matrix& operator *=(Matrix &m, double);
Matrix operator *(const Matrix &m, double) ;
Matrix operator /(const Matrix &m,double) ;
Matrix& operator /=(Matrix &m, double) ;
Matrix& operator *=(Matrix &m_,const Matrix &m);
Matrix& operator +=(Matrix &m_,const Matrix &m) ;
Matrix& operator +=(Matrix &m_,const MtM &m) ;
Matrix& operator +=(Matrix &m_,const MtMtM &m) ;
Matrix operator +(const Matrix &m_,const Matrix &m) ;
Matrix& operator -=(const Matrix &m_, const Matrix &m) ;
Matrix operator -(const Matrix &m_,const Matrix &m);
Matrix operator -(const Matrix & m_) ;
Matrix operator -(const Matrix &m,const Vector &v);
Matrix operator +(const Matrix &m,const Vector &v);
Matrix &operator -=(Matrix &m_,const Vector &v);
Matrix &operator +=(Matrix &m_,const Vector &v);

bool operator ==(const Matrix &m_, const Matrix &m) ;
bool operator !=(const Matrix &m_, const Matrix &m) ;

std::pair<Vector, double> getLargestEigenValueAndVector(const Amie::Matrix & m) ;
std::vector<std::pair<Vector, double> > deflate(const Amie::Matrix & m) ;

/** \brief perform the inner product of two vector of doubles
* @param v0 Pointer to the start of the first vector
* @param v1 pointer to the start of the second vector
* @param size size of the vectors
*/
inline double parallel_inner_product(const double *  v0, const double *  v1, int size)
{
    double result = 0 ;
    double compensate = 0 ;

    #pragma omp parallel for reduction(+:result,compensate)
    for(int i = 0 ; i < size ; ++i)
    {
        double y = *(v0+i)* *(v1+i)-compensate ;
        double tot = result+y ;
        compensate  = (tot-result)-y ;
        result = tot  ;
    }
    return result ;
}

/** \brief perform the inner product of two vector of doubles
* @param v0 Pointer to the start of the first vector
* @param v1 pointer to the start of the second vector
* @param size size of the vectors
*/
inline double parallel_inner_product_restricted(const double * __restrict__ v0, const double * __restrict__ v1, int size)
{
    double result = 0 ;
    double compensate = 0 ;
    #pragma omp parallel for reduction(+:result,compensate)
    for(int i = 0 ; i < size ; ++i)
    {
        double y = *(v0+i)* *(v1+i)-compensate ;
        double tot = result+y ;
        compensate  = (tot-result)-y ;
        result = tot  ;
    }
    return result ;
}


/** \brief Structure used for Matrix-Vector multiplication, minimising temporaries*/
struct MtV
{
    const Matrix &m;
    const Vector &v;

    /** \brief constructor, initialise the references necessary for the computation*/
    MtV(const Matrix &mm, const Vector &vv) : m(mm), v(vv) { }

    /** \brief cast operator, compute the operation*/
    operator  Vector() const;
} ;

/** \brief Structure used for Matrix-Matrix-Matrix multiplication, minimising the use of temporares*/
struct MtMtM
{
    const Matrix &first;
    const Matrix &second;
    const Matrix &third;

    /** \brief constructor, initialise the references necessary for the computation*/
    MtMtM(const Matrix &mm, const Matrix &mmm, const Matrix &mmmm) : first(mm), second(mmm), third(mmmm) { }

    Matrix operator *(const double & d) const ;

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
    MtMtM operator *(const Matrix &mm) const {
        return MtMtM(first, second, mm) ;
    }

    Matrix operator *(const double & d) const ;

    /** \brief cast operator, compute the operation*/
    operator const Matrix() const;
} ;

/** \brief Solve Ax = b system of equation, using a preconditionned Conjugate Gradient Algorithm*/
Vector solveSystem(const Matrix & A, const Vector & b, Vector & x) ;



/** \brief Construct a Matrix-Vector multiplication object*/
inline Amie::MtV operator*(const Amie::Matrix& mm, const Vector& v)
{
    return Amie::MtV(mm, v) ;
}

/** \brief Construct a Matrix-Matrix multiplication object*/
inline Amie::MtM operator*(const Amie::Matrix& mm, const Amie::Matrix& mmm)
{
    return Amie::MtM(mm, mmm) ;
}

/** \brief Perform a Matrix-Matrix multiplication */
inline const Amie::Matrix matrix_multiply(const Amie::Matrix &m0, const Amie::Matrix &m1 )
{
    assert(m0.numCols() == m1.numRows()) ;

    Amie::Matrix ret(m0.numRows(), m1.numCols()) ;

    for(size_t i = 0 ; i < m0.numRows() ; i++)
    {
        for(size_t j = 0 ; j < m1.numCols() ; j++)
        {
            const Amie::Cslice_iter<double>& ri = m0[i] ;
            const Amie::Cslice_iter<double>& cj = m1.column(j) ;
            ret[i][j] = std::inner_product(&ri[0], &ri[m0.numCols()], cj, (double)(0) ) ;
        }
    }
    return ret ;
}

inline const Amie::Matrix matrix_multiply(const Amie::Matrix &m0, const Amie::Matrix &m1 , const double & d)
{
    assert(m0.numCols() == m1.numRows()) ;

    Amie::Matrix ret(m0.numRows(), m1.numCols()) ;

    for(size_t i = 0 ; i < m0.numRows() ; i++)
    {
        for(size_t j = 0 ; j < m1.numCols() ; j++)
        {
            const Amie::Cslice_iter<double>& ri = m0[i] ;
            const Amie::Cslice_iter<double>& cj = m1.column(j) ;
            ret[i][j] = std::inner_product(&ri[0], &ri[m0.numCols()], cj, (double)(0) )*d ;
        }
    }
    return ret ;
}

/** \brief Perform a Matrix-Matrix-Matrix multiplication*/
inline const Amie::Matrix matrix_matrix_matrix_multiply(const Amie::Matrix &m0, const Amie::Matrix &m1, const Amie::Matrix &m2 )
{
    Amie::Matrix ret(m0.numRows(), m2.numCols()) ;

    for(size_t i = 0 ; i < m0.numRows() ; i++)
    {
        for(size_t j = 0 ; j < m1.numCols() ; j++)
        {
            const Amie::Cslice_iter<double>& ri = m0[i] ;
            const Amie::Cslice_iter<double>& cj = m1.column(j) ;
            double r_ij = std::inner_product(&ri[0], &ri[m0.numCols()], cj, (double)(0) ) ;
            for(size_t k = 0 ; k < m2.numCols() ; k++)
                ret[i][k] += r_ij*m2[j][k] ;
        }
    }

    return ret ;
}

inline const Amie::Matrix matrix_matrix_matrix_multiply(const Amie::Matrix &m0, const Amie::Matrix &m1, const Amie::Matrix &m2 , const double & d)
{
    Amie::Matrix ret(m0.numRows(), m2.numCols()) ;

    for(size_t i = 0 ; i < m0.numRows() ; i++)
    {
        for(size_t j = 0 ; j < m1.numCols() ; j++)
        {
            const Amie::Cslice_iter<double>& ri = m0[i] ;
            const Amie::Cslice_iter<double>& cj = m1.column(j) ;
            double r_ij = std::inner_product(&ri[0], &ri[m0.numCols()], cj, (double)(0) ) ;
            for(size_t k = 0 ; k < m2.numCols() ; k++)
                ret[i][k] += r_ij*m2[j][k]*d ;
        }
    }

    return ret ;
}

/** \brief Perform a Matrix-Matrix-Matrix multiplication, and assign to already initialised matrix passed as further argument*/
inline void matrix_matrix_matrix_multiply_and_assign(const Amie::Matrix &m0, const Amie::Matrix &m1, const Amie::Matrix &m2, Amie::Matrix & ret )
{
    if(&m0 == &ret || &m1 == &ret || &m2 == &ret)
    {
        Amie::Matrix r(ret.numRows(), ret.numCols()) ;
        for(size_t i = 0 ; i < m0.numRows() ; i++)
        {
            for(size_t j = 0 ; j < m1.numCols() ; j++)
            {
                const Amie::Cslice_iter<double>& ri = m0[i] ;
                const Amie::Cslice_iter<double>& cj = m1.column(j) ;
                double r_ij = std::inner_product(&ri[0], &ri[m0.numCols()], cj, (double)(0) ) ;
                for(size_t k = 0 ; k < m2.numCols() ; k++)
                    r[i][k] += r_ij*m2[j][k] ;
            }
        }
        ret = r ;
        return ;
    }

    ret.array() = 0 ;
    for(size_t i = 0 ; i < m0.numRows() ; i++)
    {
        for(size_t j = 0 ; j < m1.numCols() ; j++)
        {
            const Amie::Cslice_iter<double>& ri = m0[i] ;
            const Amie::Cslice_iter<double>& cj = m1.column(j) ;
            double r_ij = std::inner_product(&ri[0], &ri[m0.numCols()], cj, (double)(0) ) ;
            for(size_t k = 0 ; k < m2.numCols() ; k++)
                ret[i][k] += r_ij*m2[j][k] ;
        }
    }

}

/** \brief Perform a Matrix-Matrix-Matrix multiplication, multiply with a factor and add to already initialised matrix passed as further argument*/
inline void matrix_matrix_matrix_multiply_and_add(const Amie::Matrix &m0, const Amie::Matrix &m1, const Amie::Matrix &m2, double f,  Amie::Matrix & ret )
{
    if(&m0 == &ret || &m1 == &ret || &m2 == &ret)
    {
        Amie::Matrix r(ret.numRows(), ret.numCols()) ;
        Amie::Matrix c(ret.numRows(), ret.numCols()) ;
        for(size_t i = 0 ; i < m0.numRows() ; i++)
        {
            for(size_t j = 0 ; j < m1.numCols() ; j++)
            {
                const Amie::Cslice_iter<double>& ri = m0[i] ;
                const Amie::Cslice_iter<double>& cj = m1.column(j) ;
                double r_ij = std::inner_product(&ri[0], &ri[m0.numCols()], cj, (double)(0) ) ;
                for(size_t k = 0 ; k < m2.numCols() ; k++)
                {
                    double y = r_ij*m2[j][k]*f-c[i][k] ;
                    double t = r[i][k]+y ;
                    c[i][k] = (t-r[i][k])-y ;
                    r[i][k] = t ;
//                     r[i][k] += r_ij*m2[j][k]*f ;
                }
            }
        }
        ret = r ;
        return ;
    }

    Amie::Matrix c(ret.numRows(), ret.numCols()) ;
    for(size_t i = 0 ; i < m0.numRows() ; i++)
    {
        for(size_t j = 0 ; j < m1.numCols() ; j++)
        {
            const Amie::Cslice_iter<double>& ri = m0[i] ;
            const Amie::Cslice_iter<double>& cj = m1.column(j) ;
            double r_ij = std::inner_product(&ri[0], &ri[m0.numCols()], cj, (double)(0) ) ;
            for(size_t k = 0 ; k < m2.numCols() ; k++)
            {
                double y = r_ij*m2[j][k]*f -c[i][k];
                double t = ret[i][k]+y ;
                c[i][k] = (t-ret[i][k])-y ;
                ret[i][k] = t ;
//                 ret[i][k] += r_ij*m2[j][k]*f ;
            }
        }
    }
}

/** \brief Perform a Matrix-Matrix multiplication, and assign to already initialised matrix passed as further argument*/
inline void matrix_multiply_and_assign(const Amie::Matrix &m0, const Amie::Matrix &m1,  Amie::Matrix &ret) ;

/** \brief Perform a Matrix-Matrix multiplication, and add to already initialised matrix passed as further argument*/
inline void matrix_multiply_and_add(const Amie::Matrix &m0, const Amie::Matrix &m1,  Amie::Matrix &ret)
{
    assert(m0.numCols() == m1.numRows()) ;
    if(&m0 == &ret || &m1 == &ret)
    {
        Amie::Matrix r(ret.numRows(), ret.numCols()) ;
        for(size_t i = 0 ; i < m0.numRows() ; i++)
        {
            for(size_t j = 0 ; j < m1.numCols() ; j++)
            {
                const Amie::Cslice_iter<double>& ri = m0[i] ;
                const Amie::Cslice_iter<double>& cj = m1.column(j) ;
                r[i][j] += std::inner_product(&ri[0], &ri[m0.numCols()], cj, (double)(0) ) ;
            }
        }
        ret = r ;
        return ;
    }


    for(size_t i = 0 ; i < m0.numRows() ; i++)
    {
        for(size_t j = 0 ; j < m1.numCols() ; j++)
        {
            const Amie::Cslice_iter<double>& ri = m0[i] ;
            const Amie::Cslice_iter<double>& cj = m1.column(j) ;
            ret[i][j] += std::inner_product(&ri[0], &ri[m0.numCols()], cj, (double)(0) ) ;
        }
    }
}

/** \brief Perform a Matrix-Vector multiplication*/
inline const Vector matrix_vector_multiply(const Amie::Matrix &m, const Vector &v )
{
    assert(m.numRows() == v.size()) ;

    Vector ret(m.numRows()) ;

    for(size_t i = 0 ; i < m.numRows() ; i++)
    {

        const Amie::Cslice_iter<double>& ri = m[i] ;
        ret[i] = std::inner_product(ri, ri.end(), &v[0], (double)(0) ) ;
    }
    return ret ;
}

/** \brief Perform a Vector-Matrix multiplication*/
inline Vector operator*(const Vector &v , const Amie::Matrix &m )
{
    assert(m.numCols() == v.size()) ;

    Vector ret(m.numCols()) ;

    for(size_t i = 0 ; i < m.numCols() ; i++)
    {
        const Amie::Cslice_iter<double>& ri = m.column(i) ;
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
Amie::Matrix swapMatrix(size_t l0, size_t l1, size_t dim = 2) ;
void swapLines(size_t l0, size_t l1, Amie::Matrix * m) ;

double trace(const Matrix & m) ;

double secondInvariant(const Amie::Matrix & m) ;

Amie::Matrix identity(size_t i) ;

//clever 2x2 Matrix inversion. Thanks the numerical cookbook :)

/** \brief 2x2 matrix inversion*/
Amie::Matrix inverse2x2Matrix(const Amie::Matrix &s) ;

/** \brief in-place 2x2 matrix inversion*/
void invert2x2Matrix(Amie::Matrix &s) ;

/** \brief 3x3 matrix inversion*/
Amie::Matrix inverse3x3Matrix(const Amie::Matrix &m) ;

/** \brief in-place 3x3 matrix inversion*/
void invert3x3Matrix(Amie::Matrix & m) ;

void invert3x3SymetricMatrix(Amie::Matrix & m) ;

/** \brief 4x4 matrix inversion*/
Amie::Matrix inverse4x4Matrix(const Amie::Matrix &s) ;

/** \brief 6x6 matrix inversion*/
Amie::Matrix inverse6x6Matrix(const Amie::Matrix &s) ;

/** \brief in-place 6x6 matrix inversion*/
void invert6x6Matrix(Matrix & m) ;

/** \brief compute the determinant of a Matrix, slow*/
double det(const Matrix &s) ;

/** \brief Compute the exp of a Matrix, using series expansion. May not converge*/
Amie::Matrix mexp(const Amie::Matrix&, size_t order = 20) ;

/** \brief Compute the log of a Matrix using series Expansion, May not converge*/
Amie::Matrix mlog(const Amie::Matrix&, size_t order) ;

}

# endif  // __MATRIXOPS_H__

