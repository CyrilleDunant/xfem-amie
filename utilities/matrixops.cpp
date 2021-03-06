
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "matrixops.h"
#include <limits>
#include <iomanip>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

namespace Amie {

std::pair<Vector, double> getLargestEigenValueAndVector(const Matrix & m)
{
    Vector evec(1., m.numCols()) ;
    evec = m*evec ;
    double eval = sqrt(std::inner_product(&evec[0], &evec[evec.size()], &evec[0], 0.)) ;
    evec /= eval ;
    double error = eval ;
    double initialError = eval ;
    while(error  > 1e-24)
    {
        double peval = eval ;
        evec = m*evec ;
        eval = sqrt(std::inner_product(&evec[0], &evec[evec.size()], &evec[0], 0.)) ;
        evec /= eval ;
        error = std::abs(peval-eval)/initialError ;
    }
    return std::make_pair(evec, eval) ;
}

std::vector<std::pair<Vector, double> > deflate(const Matrix & m)
{
    std::vector<std::pair<Vector, double> > ret ;
    std::pair<Vector, double> largest = getLargestEigenValueAndVector(m) ;
    ret.push_back(largest);
    Matrix current = m-Matrix(largest.first, largest.first)*largest.second ;
    for(size_t j = 1 ; j < m.numCols() ; j++)
    {
        largest = getLargestEigenValueAndVector(current) ;
        ret.push_back(largest);
        current = current-Matrix(largest.first, largest.first)*largest.second ;
    }

    std::reverse(ret.begin(), ret.end());
    return ret ;
}

Matrix::Matrix(size_t x, size_t y)
{
    noClean = false ;
    r=x ;
    c=y ;
    v = new Vector(0., x*y) ;
}

Matrix::Matrix(const Vector & v1, const Vector & v2)
{
    noClean = false ;
    r=v1.size() ;
    c=v1.size() ;
    v = new Vector(0., r*c) ;
    for(size_t i = 0 ; i < r ; i++)
    {
        for(size_t j = 0 ; j < c ; j++)
        {
            (*this)[i][j] = v1[i]* v2[j] ;
        }
    }
}

Matrix::Matrix(size_t x, size_t y, Vector * cache)
{
    noClean = true ;
    r=x ;
    c=y ;
    v = cache ;
}

Matrix::Matrix(size_t rl, size_t k, size_t l, const Matrix & m)
{
    noClean = false ;
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

inline void matrix_multiply_and_assign(const Matrix &m0, const Matrix &m1,  Matrix &ret)
{
    assert(m0.numCols() == m1.numRows()) ;
    if(&m0 == &ret || &m1 == &ret)
    {
        Matrix r(ret.numRows(), ret.numCols()) ;
        for(size_t i = 0 ; i < m0.numRows() ; i++)
        {
            const Cslice_iter<double>& ri = m0[i] ;
            for(size_t j = 0 ; j < m1.numCols() ; j++)
            {
                const Cslice_iter<double>& cj = m1.column(j) ;

                r[i][j] = std::inner_product(&ri[0], &ri[m0.numCols()], cj, (double)(0) ) ;

            }
        }
        ret = r ;
        return ;
    }


    for(size_t i = 0 ; i < m0.numRows() ; i++)
    {
        const Cslice_iter<double>& ri = m0[i] ;
        for(size_t j = 0 ; j < m1.numCols() ; j++)
        {
            const Cslice_iter<double>& cj = m1.column(j) ;

            ret[i][j] = std::inner_product(&ri[0], &ri[m0.numCols()], cj, (double)(0) ) ;

        }

    }
}

Matrix::Matrix(size_t i, size_t j, const Matrix & m)
{
    noClean = false ;
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

Matrix::Matrix(const std::vector<std::vector<double> > & vv)
{
    noClean = false ;
    r = vv.size() ;
    c = vv[0].size() ;
    v = new Vector(r*c) ;
    for(size_t i = 0 ; i < vv.size() ; i++)
    {
        for(size_t j = 0 ; j < std::min(vv[0].size(), vv[i].size()) ; j++)
        {
            (*v)[i*r+j] = vv[i][j] ;
        }
    }
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

void Matrix::transpose(Matrix & m) const
{
    for(size_t i = 0 ; i < r ; i++)
    {
        for(size_t j = 0 ; j < c ; j++)
        {
            m[j][i] = (*this)[i][j] ;
        }
    }
}

Matrix identity(size_t i)
{
    Matrix I(i, i) ;
    for(size_t j = 0 ; j < i ; ++j)
        I[j][j] = 1 ;

    return I ;
}

double trace(const Matrix & m)
{
    double ret = 0 ;
    for(size_t i = 0 ; i < std::min(m.numCols(), m.numRows()) ; ++i)
        ret += m[i][i] ;

    return ret ;
}

double secondInvariant(const Matrix & m)
{
    double tr = trace(m) ;
    return tr*tr-trace(m*m) ;
}


Matrix &operator*=(Matrix &m, double d)
{
    m.array() *=d ;
    return m ;
}


Matrix operator*(const Matrix &m, double d)
{
    Matrix ret(m) ;
    ret*=d ;

    return ret ;
}

Matrix operator*( double d, const Matrix &m)
{
    Matrix ret(m) ;
    ret*=d ;

    return ret ;
}

Matrix operator/(const Matrix &m, double d)
{
    Matrix ret(m) ;
    ret /=d ;
    return ret ;
}

Matrix &operator/=( Matrix &m,double d)
{
    m.array() /=d ;
    return m ;
}

Matrix & Matrix::operator =(const MtM& m)
{
    if(isNull())
        std::cout << "AAAAARGH" << std::endl ;

    matrix_multiply_and_assign(m.first, m.second, *this) ;

    return *this ;
}

Matrix & operator +=( Matrix & m_, const MtM& m)
{
    if(m_.isNull())
        std::cout << "EEEEEERGH" << std::endl ;

    matrix_multiply_and_add(m.first, m.second, m_) ;

    return m_ ;
}

Matrix & Matrix::operator =(const MtMtM& m)
{
    if(isNull())
        std::cout << "IIIIIIK" << std::endl ;

    matrix_matrix_matrix_multiply_and_assign(m.first, m.second, m.third, *this) ;

    return *this ;
}

Matrix & operator +=(Matrix &m_ , const MtMtM& m)
{
    matrix_matrix_matrix_multiply_and_add(m.first, m.second, m.third, 1.,  m_) ;
    return m_ ;
}

Matrix & operator*=(Matrix &m_, const Matrix &m)
{
    assert(m.numRows() == m_.numCols()) ;

    Matrix ret(matrix_multiply(m_, m)) ;

    m_.resize(ret.numRows(),ret.numCols()) ;
    m_.array() = ret.array() ;

    return m_ ;
}


MtMtM::operator const Matrix() const
{
    return matrix_matrix_matrix_multiply(first, second, third) ;
}

Matrix MtMtM::operator *(const double & d) const
{
    return matrix_matrix_matrix_multiply(first, second, third, d) ;
}

Matrix &operator +=(Matrix & m_, const Matrix &m)
{
    m_.array() += m.array() ;
    return m_ ;
}

Matrix operator +(const Matrix & m_, const  Matrix &m)
{
    Matrix ret(m_) ;
    ret += m ; ;
    return ret ;
}

Matrix &operator -=( Matrix & m_, const Matrix &m)
{
    m_.array() -= m.array() ;
    return m_ ;
}

Matrix operator -(const Matrix & m_,const  Matrix &m)
{
    Matrix ret(m_) ;
    ret -= m ; ;
    return ret ;
}
Matrix operator -(const Matrix & m_)
{
    Matrix ret(m_) ;
    ret.array() = -ret.array() ;
    return ret ;
}

Matrix operator -(const Matrix & m_,const Vector &v)
{
    Matrix ret(m_) ;
    ret.array() -= v ;
    return ret ;
}

Matrix operator +(const Matrix & m_,const Vector &v)
{
    Matrix ret(m_) ;
    ret.array() += v ;
    return ret ;
}

Matrix& operator -=( Matrix & m_,const Vector &v)
{
    m_.array() -= v ;
    return m_ ;
}

Matrix& operator +=(Matrix & m_,const Vector &v)
{
    m_.array() += v ;
    return m_ ;
}

bool operator ==(const Matrix & m_,const Matrix &m)
{
    if(m_.array().size() != m.array().size())
        return false ;

    for(size_t i = 0 ; i < m_.array().size() ; i++)
        if(std::abs(m_.array()[i] - m.array()[i]) > std::numeric_limits<double>::epsilon()*std::max(m_.array()[i], m.array()[i]))
            return false ;

    return true ;
}

bool operator !=(const Matrix & m_,const Matrix &m)
{
    if(m_.array().size() != m.array().size())
        return true ;
    else
        for(size_t i = 0 ; i < m_.array().size() ; i++)
            if(m_.array()[i] != m.array()[i])
                return true ;

    return false ;
}

Matrix::Matrix(const Matrix& m) : r(m.numRows()), c( m.numCols())
{
    noClean = false ;
    if(!m.isNull())
    {
        v = new Vector(m.array()) ;
    }
    else
    {
        v = nullptr ;
    }
}

Matrix &Matrix::operator =(double a)
{
    *v = a ;
    return *this ;
}

Matrix &Matrix::operator =(const Matrix &m)
{
    
    noClean = false ;
// 	v->resize(m.array().size()) ;
// 	v = new Vector(m.array()) ;
    if(v == nullptr && !m.isNull())
    {
        r = m.numRows() ;
        c = m.numCols() ;

        v = new Vector(m.array().size()) ;
        *v = m.array() ;
    }
    else if(!m.isNull())
    {
        if(r != m.numRows() || c != m.numCols())
        {
            r = m.numRows() ;
            c = m.numCols() ;
            v->resize(m.array().size()) ;
        }
        *v = m.array() ;
    }
    else
    {
        delete v ;
        v = nullptr ;
    }
// 	r = m.numRows() ;
// 	c = m.numCols() ;
    return *this ;
}

void Matrix::print() const
{
    for(size_t i = 0 ; i < numRows() ; i++)
    {
        for(size_t j = 0 ; j < numCols() ; j++)
        {
            std::cout << std::setprecision(16) << (*this)[i][j] << "   " << std::flush ;
        }

        std::cout << std::setprecision(8) << std::endl ;
    }
}

MtM::operator const Matrix() const
{
    return matrix_multiply(first, second) ;
}

Matrix MtM::operator *(const double & d) const
{
    return matrix_multiply(first, second, d) ;
}


MtV::operator  Vector() const
{
    return matrix_vector_multiply(m,v) ;
}


void swapLines(size_t l0, size_t l1, Matrix * m)
{
    for(size_t i = 0 ; i < m->numCols() ; i++)
    {
        double temp = (*m)[l1][i] ;
        (*m)[l1][i] = (*m)[l0][i] ;
        (*m)[l0][i] = temp ;
// 		std::swap(m[l0][i], m[l1][i]) ;
    }
}

Matrix swapMatrix(size_t l0, size_t l1, size_t dim)
{

    Matrix ret(dim, dim) ;
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
Matrix inverse2x2Matrix(const Matrix &s)
{
    Matrix s_(s) ;
    invert2x2Matrix(s_) ;
    return s_ ;
}

/**
 *
 * @param s
 */
void invert2x2Matrix(Matrix &s)
{

    if(std::abs(s.array()[0]) < 1e-24)
    {
        double a = s.array()[0] ;
        double b = s.array()[1] ;
        double c = s.array()[2] ;
        double d = s.array()[3] ;
        double dt = 1./(a*d-b*c) ;
        s.array()[0] = d*dt ;
        s.array()[1] = -b*dt ;
        s.array()[2] = -c*dt ;
        s.array()[3] = a*dt ;
        return ;
    }

    double r1 = 1./s.array()[0] ;
    double r3 = r1* s.array()[1] ;
    double r6 = 1./(s.array()[2] * r3 - s.array()[3]) ;
    s.array()[1] = r3*r6 ;
    s.array()[2] = r6*s.array()[2] * r1 ;
    s.array()[0] = r1-r3*s.array()[2] ;
    s.array()[3] = -r6 ;

}

void invert6x6Matrix(Matrix &s)
{
    Matrix s00(3,0,0,s) ;
    Matrix s01(3,0,3,s) ;
    Matrix s10(3,3,0,s) ;
    Matrix s11(3,3,3,s) ;

    if(true)
    {
        Matrix ret(6,6) ;


        double deti = det(s) ;

        ret[0][0] =  det(Matrix(0,0,s)) ;
        ret[0][1] = -det(Matrix(1,0,s)) ;
        ret[0][2] =  det(Matrix(2,0,s)) ;
        ret[0][3] = -det(Matrix(3,0,s)) ;
        ret[0][4] =  det(Matrix(4,0,s)) ;
        ret[0][5] = -det(Matrix(5,0,s)) ;

        ret[1][0] = -det(Matrix(0,1,s)) ;
        ret[1][1] =  det(Matrix(1,1,s)) ;
        ret[1][2] = -det(Matrix(2,1,s)) ;
        ret[1][3] =  det(Matrix(3,1,s)) ;
        ret[1][4] = -det(Matrix(4,1,s)) ;
        ret[1][5] =  det(Matrix(5,1,s)) ;

        ret[2][0] =  det(Matrix(0,2,s)) ;
        ret[2][1] = -det(Matrix(1,2,s)) ;
        ret[2][2] =  det(Matrix(2,2,s)) ;
        ret[2][3] = -det(Matrix(3,2,s)) ;
        ret[2][4] =  det(Matrix(4,2,s)) ;
        ret[2][5] = -det(Matrix(5,2,s)) ;

        ret[3][0] = -det(Matrix(0,3,s)) ;
        ret[3][1] =  det(Matrix(1,3,s)) ;
        ret[3][2] = -det(Matrix(2,3,s)) ;
        ret[3][3] =  det(Matrix(3,3,s)) ;
        ret[3][4] = -det(Matrix(4,3,s)) ;
        ret[3][5] =  det(Matrix(5,3,s)) ;

        ret[4][0] =  det(Matrix(0,4,s)) ;
        ret[4][1] = -det(Matrix(1,4,s)) ;
        ret[4][2] =  det(Matrix(2,4,s)) ;
        ret[4][3] = -det(Matrix(3,4,s)) ;
        ret[4][4] =  det(Matrix(4,4,s)) ;
        ret[4][5] = -det(Matrix(5,4,s)) ;

        ret[5][0] = -det(Matrix(0,5,s)) ;
        ret[5][1] =  det(Matrix(1,5,s)) ;
        ret[5][2] = -det(Matrix(2,5,s)) ;
        ret[5][3] =  det(Matrix(3,5,s)) ;
        ret[5][4] = -det(Matrix(5,4,s)) ;
        ret[5][5] =  det(Matrix(5,5,s)) ;


        s = ret/deti ;

        return ;
    }

    Matrix r1 = inverse3x3Matrix(s00) ;
    Matrix r2 = s10 * r1 ;
    Matrix r3 = r1* s01 ;
    Matrix r4 = s10 * r3 ;
    Matrix r5 = r4 - s11 ;
    Matrix r6 = inverse3x3Matrix(r5) ;
    s01 = r3*r6 ;
    s10 = r6*r2 ;
    Matrix r7 = r3*s10 ;
    s00 = r1-r7 ;
    s11 = -r6 ;

    for(size_t i = 0 ; i < 3 ; i++)
    {
        for(size_t j = 0 ; j < 3 ; j++)
        {
            s[i][j] = s00[i][j] ;
        }
    }
    for(size_t i = 0 ; i < 3 ; i++)
    {
        for(size_t j = 0 ; j < 3 ; j++)
        {
            s[i+3][j] = s10[i][j] ;
        }
    }
    for(size_t i = 0 ; i < 3 ; i++)
    {
        for(size_t j = 0 ; j < 3 ; j++)
        {
            s[i][j+3] = s01[i][j] ;
        }
    }
    for(size_t i = 0 ; i < 3 ; i++)
    {
        for(size_t j = 0 ; j < 3 ; j++)
        {
            s[i+3][j+3] = s11[i][j] ;
        }
    }
}

Matrix inverse6x6Matrix(Matrix &s)
{
    Matrix ret(s) ;
    invert6x6Matrix(ret);
    return ret ;
}

Matrix inverse3x3Matrix(const Matrix & m)
{
    Matrix m_(m) ;
    invert3x3Matrix(m_) ;
    return m_ ;
}

void invert3x3Matrix(Matrix & m)
{
    double m00 = m.array()[0] ;
    double m01 =  m.array()[1] ;
    double m02 =  m.array()[2] ;
    double m10 = m.array()[3] ;
    double m11 =  m.array()[4] ;
    double m12 =  m.array()[5] ;
    double m20 = m.array()[6] ;
    double m21 =  m.array()[7] ;
    double m22 =  m.array()[8] ;

    double r11 = fma(m11,m22,-m12*m21) ;
    double r21 = fma(m12,m20,-m10*m22) ;
    double r31 = fma(m10,m21,-m11*m20) ;
    double det = 1./fma(m00,(r11),fma( m01,r21, m02*r31)) ;

    Vector v = {r11, fma(m02,m21, - m01*m22), fma(m01,m12, - m02*m11), r21, fma(m00,m22, - m02*m20), fma(m02,m10, - m00*m12), r31, fma(m01,m20, - m00*m21), fma(m11,m00, - m01*m10) } ;
    v *= det ;
    m.array() = v ;

}

void invert3x3SymetricMatrix(Matrix & m)
{
    double m00 = m.array()[0] ;
    double m01 =  m.array()[1] ;
    double m02 =  m.array()[2] ;
    double m11 =  m.array()[4] ;
    double m12 =  m.array()[5] ;
    double m22 =  m.array()[8] ;

    double r11 = fma(m11,m22,-m12*m12) ;
    double r21 = fma(m12,m02,-m01*m22) ;
    double r31 = fma(m01,m12,-m11*m02) ;
    double det = 1./fma(m00,(r11),fma( m01,r21, m02*r31)) ;

    Vector v = {r11, r21, r31, r21, fma(m00,m22, - m02*m02), fma(m02,m01, - m00*m12), r31, fma(m01,m02, - m00*m12), fma(m11,m00, - m01*m01) } ;
    v *= det ;
    m.array() = v ;

}

Matrix inverse4x4Matrix(const Matrix &s)
{

    Matrix ret(4,4) ;
//     Matrix a(2, 0, 0, s) ; Matrix b(2, 0, 2, s) ;
//     Matrix c(2, 2, 0, s) ; Matrix d(2, 2, 2, s) ;
//     Matrix r1 = inverse2x2Matrix(a) ;
//     Matrix r3 = r1*b ;
//     Matrix r6 = (Matrix)(c * r3) - d ; invert2x2Matrix(r6);
//     Matrix s00 = r1 - r3 * c ; Matrix s01 = r3 * r6 ;
//     Matrix s10 = r6 * c * r1 ;
//
//     ret[0][0] = s00[0][0]; ret[0][1] = s00[0][1]; ret[0][2] = s01[0][0]; ret[0][3] = s01[0][1];
//     ret[1][0] = s00[1][0]; ret[1][1] = s00[1][1]; ret[1][2] = s01[1][0]; ret[1][3] = s01[1][1];
//     ret[2][0] = s10[0][0];ret[2][1] =  s10[0][1]; ret[2][2] = -r6[0][0]; ret[3][3] = -r6[0][1];
//     ret[3][0] = s10[1][0];ret[3][1] =  s10[1][1]; ret[3][2] = -r6[1][0]; ret[3][3] = -r6[1][1];
//     return ret ;
//     double r1 = 1. / s.array()[0] ;
//     double r3 = r1 * s.array()[1] ;
//     double r6 = 1. / (s.array()[2] * r3 - s.array()[3]) ;
//     s.array()[1] = r3 * r6 ;
//     s.array()[2] = r6 * s.array()[2] * r1 ;
//     s.array()[0] = r1 - r3 * s.array()[2] ;
//     s.array()[3] = -r6 ;

//     Matrix a(2, 0, 0, s) ; Matrix c(2, 0, 2, s) ;
//     Matrix b(2, 2, 0, s) ; Matrix d(2, 2, 2, s) ;
// //
//     Matrix dt = inverse2x2Matrix((Matrix)(a*d)-(Matrix)(b*c)) ;
//     d *= dt ;
//     b *= dt ;
//     c *= dt ;
//     a *= dt ;
//     ret[0][0] = d[0][0] ; ret[0][1] = d[0][1] ; ret[0][2] = -c[0][0] ; ret[0][3] = -c[0][1] ;
//     ret[1][0] = d[1][0] ; ret[1][1] = d[1][1] ; ret[1][2] = -c[1][0] ; ret[1][3] = -c[1][1] ;
//     ret[2][0] = -b[0][0] ; ret[2][1] = -b[0][1] ; ret[2][2] = a[0][0] ; ret[2][3] = a[0][1] ;
//     ret[3][0] = -b[1][0] ; ret[3][1] = -b[1][1] ; ret[3][2] = a[1][0] ; ret[3][3] = a[1][1] ;
//     return ret ;
//
//         double a = s.array()[0] ;
//         double b = s.array()[1] ;
//         double c = s.array()[2] ;
//         double d = s.array()[3] ;
//         double dt = 1./(a*d-b*c) ;
//         s.array()[0] = d*dt ;
//         s.array()[1] = -b*dt ;
//         s.array()[2] = -c*dt ;
//         s.array()[3] = a*dt ;
//         return ;
//

    double deti = det(s) ;

    ret[0][0] =  det(Matrix(0,0,s)) ;
    ret[0][1] = -det(Matrix(1,0,s)) ;
    ret[0][2] =  det(Matrix(2,0,s)) ;
    ret[0][3] = -det(Matrix(3,0,s)) ;
    ret[1][0] = -det(Matrix(0,1,s)) ;
    ret[1][1] =  det(Matrix(1,1,s)) ;
    ret[1][2] = -det(Matrix(2,1,s)) ;
    ret[1][3] =  det(Matrix(3,1,s)) ;
    ret[2][0] =  det(Matrix(0,2,s)) ;
    ret[2][1] = -det(Matrix(1,2,s)) ;
    ret[2][2] =  det(Matrix(2,2,s)) ;
    ret[2][3] = -det(Matrix(3,2,s)) ;
    ret[3][0] = -det(Matrix(0,3,s)) ;
    ret[3][1] =  det(Matrix(1,3,s)) ;
    ret[3][2] = -det(Matrix(2,3,s)) ;
    ret[3][3] =  det(Matrix(3,3,s)) ;

    return ret/deti ;

}


// Mu::Matrix inverse4x4Matrix(const Mu::Matrix &s)
// {
// 	Matrix s00(2, 0, 0, s) ;
// 	Matrix s10(2, 2, 0, s) ;
// 	Matrix s01(2, 0, 2, s) ;
// 	Matrix s11(2, 2, 2, s) ;
// 	invert2x2Matrix(s00) ;
// 	Matrix r2 = s10 * s00 ;
// 	Matrix r3 = s00* s01 ;
// 	Matrix r4 = s10 * r3 ;
// 	Matrix r5 = r4 - s11 ;
// 	invert2x2Matrix(r5) ;
// 	s01 = r3*r5 ;
// 	s10 = r5*r2 ;
// 	Matrix r7 = r3*s10 ;
// 	s00 -= r7 ;
// 	s11 = -r5 ;
//
// 	Matrix ret(4,4) ;
// 	ret[0][0] = s00[0][0]; ret[0][1] = s00[0][1]; ret[0][2] = s01[0][0]; ret[0][3] = s01[0][1];
// 	ret[1][0] = s00[1][0]; ret[1][1] = s00[1][1]; ret[1][2] = s01[1][0]; ret[1][3] = s01[1][1];
// 	ret[2][0] = s10[0][0];ret[2][1] =  s10[0][1]; ret[2][2] = s11[0][0]; ret[3][3] = s11[0][1];
// 	ret[3][0] = s10[1][0];ret[3][1] =  s10[1][1]; ret[3][2] = s11[1][0]; ret[3][3] = s11[1][1];
//
// 	return ret ;
// }

double det(const Matrix &s)
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
    case 5:
    {
        double ret = 0 ;

        for(size_t i = 0 ;  i < 1  ; i++)
        {
            for(size_t j = 0 ;  j< 5  ; j++)
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
    case 6:
    {
        double ret = 0 ;

        for(size_t i = 0 ;  i < 1  ; i++)
        {
            for(size_t j = 0 ;  j< 6  ; j++)
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


Matrix mexp(const Matrix& m, size_t order)
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

Matrix mlog(const Matrix&m, size_t order)
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

// std::vector<std::vector<std::pair>> largestEigenValue(const Mu::CoordinateIndexedSparseMatrix & A, bool sym)
// {
//
// 	std::valarray<double> x_(A.row_size.size()*A.stride) ;
// 	for(size_t i = 0 ; i< x_.size() ; i++)
// 	{
// 		x_[i] = (double)rand()/RAND_MAX ;
// 	}
// 	x_/= sqrt(std::inner_product(&x_[0], &x_[x_.size()], &x_[0], double(0)))  ;
// 	std::valarray<double> x(x_) ;
//
// 	double eps = 1 ;
// 	while ( eps > 1e-6 )
// 	{
// 		x = A*x_ ;
// 		x/= sqrt(std::inner_product(&x[0], &x[x.size()], &x[0], double(0)))  ;
// 		std::valarray<double> delta = x-x_ ;
// 		eps = sqrt(std::inner_product(&delta[0], &delta[delta.size()], &delta[0], double(0))) ;
// 		x_ = x ;
// 	}
//
// 	x = A*x_ ;
// 	return std::inner_product(&x[0], &x[x.size()], &x_[0], double(0)) ;
// }

Vector solveSystem(const Matrix & A, const Vector & b, Vector & x)
{

    double epsilon = 1e-12 ;
    Vector r = b - (Vector)(A*x) ;
    Vector r_(r) ;
    double rho = std::inner_product(&r[0], &r[r.size()], &r_[0], double(0)) ;

    Vector invDiag(r.size()) ;

    for(size_t i = 0 ; i < r.size() ; i++)
    {
        invDiag[i] = 1./A[i][i] ;
    }

    if(rho < epsilon)
        return x ;

    Vector p(r) ;

    Vector p_ = invDiag*p ;

    Vector v = A*p_ ;

    double alpha = rho/std::inner_product(&r_[0], &r_[r_.size()], &v[0], double(0)) ;

    Vector s = r - v*alpha ;

    if(std::abs(s.max()) < epsilon)
    {
        x += p_*alpha ;
        return x ;
    }

    Vector s_ = invDiag*s ;
    Vector t = A*s_ ;
    double omega = std::inner_product(&t[0], &t[t.size()], &s[0], double(0))/std::inner_product(&t[0], &t[t.size()], &t[0], double(0)) ;
    x += p_*alpha +omega*s_ ;
    r = s- t*omega ;
    double rho_ =rho ;

    size_t nit = 0 ;

    while(nit < r.size())
    {
        nit++ ;

        rho = std::inner_product(&r[0], &r[r.size()], &r_[0], double(0)) ;
        if(std::abs(rho) < epsilon)
            return x ;

        double beta = (rho/rho_)*(alpha/omega) ;
        p = r + (p-v*omega)*beta ;
        p_ = invDiag*p ;
        v = A*p_ ;
        alpha = rho/std::inner_product(&r_[0], &r_[r_.size()], &v[0], double(0)) ;
        s = r - v*alpha ;

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

} 
