//
// C++ Implementation: sparse_vector
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "sparse_vector.h"
#include <string.h>
#ifdef HAVE_SSE3
#include <pmmintrin.h>
#endif

namespace Amie {


SparseVector::SparseVector(Vector & val, std::valarray<unsigned int> & idx , const size_t length , const size_t start, const size_t index, const size_t stride) : val(val), idx(idx), length(length), start(start), stride(stride), index(index)
{
    zero = 0 ;
}

double SparseVector::operator [](size_t i) const
{
    unsigned int * __start__       = &idx[start] ;
    unsigned int * __end__         = &idx[start+length] ;
    unsigned int * i_index_pointer = std::lower_bound(__start__, __end__, i/stride) ;
    unsigned int offset            = i_index_pointer - __start__ ;
    unsigned int colLength = stride+stride%2 ;
    return (std::binary_search(__start__, __end__, i/stride)) ?  val[(start+offset)*stride*colLength + (i%stride)*colLength+index%stride] : 0 ;

}

double & SparseVector::operator [](const size_t i)
{
// 	std::cout << i << "  " << start << std::endl;
    zero = 0 ;
    unsigned int * __start__       = &idx[start] ;
    unsigned int * __end__         = &idx[start+length] ;
    unsigned int * i_index_pointer = std::lower_bound(__start__, __end__, i/stride) ;
    unsigned int offset            = i_index_pointer - __start__ ;
    unsigned int colLength = stride+stride%2 ;
    return (std::binary_search(__start__, __end__, i/stride)) ?  val[(start+offset)*stride*colLength + (i%stride)*colLength+index%stride] :  zero ;

}

Vector SparseVector::operator *(const Vector &v) const
{
    int colLength = stride + stride%2 ;
    Vector ret(0., colLength) ;
#ifdef HAVE_SSE3
    const __m128d * array_iterator = (const __m128d *)&val[start*colLength*stride] ;
    for(unsigned int j = start ; j != length+start ; j++)
    {
        for(int c = 0 ; c < stride ; c++)
        {
            __m128d vval =  _mm_set1_pd(v[idx[j]*stride+c]) ;
            for(size_t i = 0 ; i != colLength ; i+=2)
            {
                _mm_store_pd(&ret[i], _mm_add_pd(_mm_load_pd(&ret[i]), _mm_mul_pd(*(array_iterator++), vval))) ;
            }
        }
    }
#else
    const double * array_iterator0 = &val[start*colLength*stride] ;
    const double * array_iterator1 = &val[start*colLength*stride] ;
    for(unsigned int j = start ; j != length+start ; j++)
    {
        for(size_t c = 0 ; c < stride ; c++)
        {
            double vval =  v[idx[j]*stride+c] ;
            for(int i = 0 ; i != colLength ; i+=2)
            {
                ret[i] +=  *(array_iterator0) * vval ;
                ret[i+1] +=  *(array_iterator1) * vval ;
                array_iterator0 += 2 ;
                array_iterator1 += 2 ;
            }
        }
    }
#endif
    return ret ;
}

Matrix SparseVector::operator *(const SparseVector &v) const
{
    Matrix ret(stride, stride) ;
    int colLength = stride+stride%2 ;
    int blocksize = stride*colLength ;

    size_t i = 0 ;
    size_t j = 0 ;
    while(i < length && j < v.length)
    {
        if(idx[start+i] > v.idx[v.start+ j])
            j++ ;
        else if(idx[start+i] < v.idx[v.start+ j])
            i++ ;
        else
        {
            int idx = 0 ;
            for(size_t k = 0 ; k < stride ; k++)
            {
                for(size_t l = 0 ; l < stride ; l++)
                {
                    ret[k][l] += val[start*blocksize+i*blocksize+colLength*k+l]*v.val[v.start*blocksize+j*blocksize + idx++] ;
                }
            }
            i++ ;
            j++ ;
        }
    }
    return ret ;

}

Vector SparseVector::operator +(const Vector &v) const
{
    Vector ret(v) ;

    for(size_t j = 0 ; j < length ; j++)
    {
        ret[idx[start + j]] += val[start + j] ;
    }

    return ret ;
}

ConstSparseVector::ConstSparseVector(const Vector & v, const std::valarray<unsigned int> & i , const size_t l , const size_t s, const size_t ind, const size_t st) : array(v), column_index(i), length(l), start(s), stride(st), index(ind)
{
}

void ConstSparseVector::parallel_product(const Vector &v, double *dest, const size_t rowstart, const size_t colstart ) const
{
    if(stride == 2)
    {
#ifdef HAVE_SSE3
        const __m128d * array_iterator = (__m128d*)&array[start*4] ;
        const double * vec_iterator = &v[column_index[start]*2] ;
        const double * dest_iterator = dest+(column_index[start+1]-column_index[start])*2 ;
        _mm_store_pd(dest,  _mm_add_pd( _mm_load_pd(dest), _mm_add_pd( _mm_mul_pd(*array_iterator,  _mm_set1_pd(*vec_iterator)), _mm_mul_pd(*(array_iterator+1), _mm_set1_pd(*(vec_iterator+1)))))) ;
        for(unsigned int j = start+1 ; j != length+start && j+1 != length+start; vec_iterator += (column_index[j+1]-column_index[j])*2, dest_iterator += (column_index[j+1]-column_index[j])*2,++j,array_iterator+=2)
        {
            __m128d tmp =  _mm_add_pd( _mm_load_pd(*dest_iterator), _mm_add_pd( _mm_mul_pd(*array_iterator,  _mm_set1_pd(*vec_iterator)), _mm_mul_pd(*(array_iterator+1), _mm_set1_pd(*(vec_iterator+1))))) ;
            _mm_store_pd(dest, _mm_add_pd(tmp, dest)) ;
            _mm_store_pd(dest_iterator, _mm_add_pd(tmp, dest_iterator)) ;
        }
#else
        const double * array_iterator0 = &array[start*4] ;
        const double * array_iterator1 = &array[start*4+1] ;
        const double * array_iterator2 = &array[start*4+2] ;
        const double * array_iterator3 = &array[start*4+3] ;
        double * dest_iterator = dest + (column_index[start+1]-column_index[start])*2;
        *dest_iterator += *array_iterator0*v[column_index[start]*2]+*array_iterator2*v[column_index[start]*2+1] ;
        *(dest_iterator+1) += *array_iterator1*v[column_index[start]*2]+*array_iterator3*v[column_index[start]*2+1] ;
        double tmp0 = *array_iterator0*v[column_index[start]*2]+*array_iterator2*v[column_index[start]*2+1] ;
        double tmp1 = *array_iterator1*v[column_index[start]*2]+*array_iterator3*v[column_index[start]*2+1] ;
        *dest += tmp0 ;
        *(dest+1) += tmp1 ;
        for(unsigned int j = start+1 ; j != length+start ; j++,array_iterator0+=4,array_iterator1+=4,array_iterator2+=4,array_iterator3+=4, dest_iterator +=(column_index[j+1]-column_index[j])*2)
        {
            double tmp0 = *array_iterator0*v[column_index[j]*2]+*array_iterator2*v[column_index[j]*2+1] ;
            double tmp1 = *array_iterator1*v[column_index[j]*2]+*array_iterator3*v[column_index[j]*2+1] ;
            *dest_iterator += tmp0 ;
            *(dest_iterator+1) += tmp1 ;
            *dest += tmp0 ;
            *(dest+1) += tmp1 ;
        }
#endif
        return ;
    }

    if(stride == 3)
    {
#ifdef HAVE_SSE3
        const __m128d * array_iterator = (__m128d*)&array[start*4*stride] ;
        const double * dest_iterator = dest + (column_index[start+1]-column_index[start])*3;
        const __m128d vval =  _mm_set1_pd(v[column_index[j]*stride]) ;
        _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*array_iterator, vval))) ;
        _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*array_iterator, vval))) ;
        vval =  _mm_set1_pd(v[column_index[j]*stride+1]) ;
        _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+2), vval))) ;
        _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+2), vval))) ;
        vval =  _mm_set1_pd(v[column_index[j]*stride+2]) ;
        _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+4), vval))) ;
        _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+4), vval))) ;
        for(unsigned int j = start+1 ; j != length+start ; ++j,array_iterator+=6, dest_iterator += (column_index[j+1]-column_index[j])*3)
        {
            const __m128d vval =  _mm_set1_pd(v[column_index[j]*stride]) ;
            _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*array_iterator, vval))) ;
            _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*array_iterator, vval))) ;
            _mm_store_pd((dest_iterator),  _mm_add_pd( _mm_load_pd((dest_iterator)), _mm_mul_pd(*array_iterator, vval))) ;
            _mm_store_pd((dest_iterator+2),  _mm_add_pd( _mm_load_pd((dest_iterator+2)), _mm_mul_pd(*array_iterator, vval))) ;
            vval =  _mm_set1_pd(v[column_index[j]*stride+1]) ;
            _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+2), vval))) ;
            _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+2), vval))) ;
            _mm_store_pd((dest_iterator),  _mm_add_pd( _mm_load_pd((dest_iterator)), _mm_mul_pd(*(array_iterator+2), vval))) ;
            _mm_store_pd((dest_iterator+2),  _mm_add_pd( _mm_load_pd((dest_iterator+2)), _mm_mul_pd(*(array_iterator+2), vval))) ;
            vval =  _mm_set1_pd(v[column_index[j]*stride+2]) ;
            _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+4), vval))) ;
            _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+4), vval))) ;
            _mm_store_pd((dest_iterator),  _mm_add_pd( _mm_load_pd((dest_iterator)), _mm_mul_pd(*(array_iterator+4), vval))) ;
            _mm_store_pd((dest_iterator+2),  _mm_add_pd( _mm_load_pd((dest_iterator+2)), _mm_mul_pd(*(array_iterator+4), vval))) ;
        }
#else
        const double * array_iterator0 = &array[start*4*stride] ;
        const double * array_iterator1 = &array[start*4*stride+1] ;
        double * dest_iterator = dest+ (column_index[start+1]-column_index[start])*stride;
        double vval =  v[column_index[start]*stride] ;
        *(dest+0) += *array_iterator0 * vval ;
        *(dest+1) += *array_iterator1 * vval ;
        *(dest+2) += *(array_iterator0+2) * vval ;
        *(dest+3) += *(array_iterator1+2) * vval ;
        vval =  v[column_index[start]*stride+1] ;
        *(dest+0) += *(array_iterator0+4) * vval ;
        *(dest+1) += *(array_iterator1+4) * vval ;
        *(dest+2) += *(array_iterator0+6) * vval ;
        *(dest+3) += *(array_iterator1+6) * vval ;
        vval =  v[column_index[start]*stride+2] ;
        *(dest+0) += *(array_iterator0+8) * vval ;
        *(dest+1) += *(array_iterator1+8) * vval ;
        *(dest+2) += *(array_iterator0+10) * vval ;
        *(dest+3) += *(array_iterator1+10) * vval ;
        for(unsigned int j = start+1 ; j != length+start ; ++j,dest_iterator += (column_index[j+1]-column_index[j])*stride)
        {
            double vval =  v[column_index[j]*stride] ;
            *(dest+0) += *array_iterator0 * vval ;
            *(dest+1) += *array_iterator1 * vval ;
            *(dest+2) += *(array_iterator0+2) * vval ;
            *(dest+3) += *(array_iterator1+2) * vval ;
            *(dest_iterator+0) += *array_iterator0 * vval ;
            *(dest_iterator+1) += *array_iterator1 * vval ;
            *(dest_iterator+2) += *(array_iterator0+2) * vval ;
            *(dest_iterator+3) += *(array_iterator1+2) * vval ;
            vval =  v[column_index[j]*stride+1] ;
            *(dest+0) += *(array_iterator0+4) * vval ;
            *(dest+1) += *(array_iterator1+4) * vval ;
            *(dest+2) += *(array_iterator0+6) * vval ;
            *(dest+3) += *(array_iterator1+6) * vval ;
            *(dest_iterator+0) += *(array_iterator0+4) * vval ;
            *(dest_iterator+1) += *(array_iterator1+4) * vval ;
            *(dest_iterator+2) += *(array_iterator0+6) * vval ;
            *(dest_iterator+3) += *(array_iterator1+6) * vval ;
            vval =  v[column_index[j]*stride+2] ;
            *(dest+0) += *(array_iterator0+8) * vval ;
            *(dest+1) += *(array_iterator1+8) * vval ;
            *(dest+2) += *(array_iterator0+10) * vval ;
            *(dest+3) += *(array_iterator1+10) * vval ;
            *(dest_iterator+0) +=  *(array_iterator0+8) * vval ;
            *(dest_iterator+1) +=  *(array_iterator1+8) * vval ;
            *(dest_iterator+2) +=  *(array_iterator0+10) * vval ;
            *(dest_iterator+3) +=  *(array_iterator1+10) * vval ;
        }
#endif
    }

    if(stride == 6)
    {
        const int colLength = stride + stride%2 ;
#ifdef HAVE_SSE3
        const __m128d * array_iterator = (__m128d*)&array[start*colLength*stride] ;
        const __m128d * dest_iterator = (__m128d*)(dest + (column_index[start+1]-column_index[start])*stride);

        __m128d vval =  _mm_set1_pd(v[column_index[j]*stride]) ;

        _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*array_iterator, vval))) ;
        _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+1), vval))) ;
        _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+2), vval))) ;

        vval =  _mm_set1_pd(v[column_index[j]*stride+1]) ;

        _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+3,) vval))) ;
        _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+4), vval))) ;
        _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+5), vval))) ;

        vval =  _mm_set1_pd(v[column_index[j]*stride+2]) ;

        _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+6), vval))) ;
        _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+7), vval))) ;
        _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+8), vval))) ;

        vval =  _mm_set1_pd(v[column_index[j]*stride+3]) ;

        _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+9), vval))) ;
        _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+10), vval))) ;
        _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+11), vval))) ;

        vval =  _mm_set1_pd(v[column_index[j]*stride+4]) ;

        _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+12), vval))) ;
        _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+13), vval))) ;
        _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+14), vval))) ;

        vval =  _mm_set1_pd(v[column_index[j]*stride+5]) ;

        _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+15), vval))) ;
        _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+16), vval))) ;
        _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+17), vval))) ;
        array_iterator +=18 ;
        for(unsigned int j = start ; j != length+start ; j++)
        {
            __m128d vval =  _mm_set1_pd(v[column_index[j]*stride]) ;

            _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*array_iterator, vval))) ;
            _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+1), vval))) ;
            _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+2), vval))) ;
            _mm_store_pd((dest_iterator),  _mm_add_pd( _mm_load_pd((dest_iterator)), _mm_mul_pd(*array_iterator, vval))) ;
            _mm_store_pd((dest_iterator+2),  _mm_add_pd( _mm_load_pd((dest_iterator+2)), _mm_mul_pd(*(array_iterator+1), vval))) ;
            _mm_store_pd((dest_iterator+4),  _mm_add_pd( _mm_load_pd((dest_iterator+4)), _mm_mul_pd(*(array_iterator+2), vval))) ;

            vval =  _mm_set1_pd(v[column_index[j]*stride+1]) ;

            _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+3,) vval))) ;
            _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+4), vval))) ;
            _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+5), vval))) ;
            _mm_store_pd((dest_iterator),  _mm_add_pd( _mm_load_pd((dest_iterator)), _mm_mul_pd(*(array_iterator+3), vval))) ;
            _mm_store_pd((dest_iterator+2),  _mm_add_pd( _mm_load_pd((dest_iterator+2)), _mm_mul_pd(*(array_iterator+4), vval))) ;
            _mm_store_pd((dest_iterator+4),  _mm_add_pd( _mm_load_pd((dest_iterator+4)), _mm_mul_pd(*(array_iterator+5), vval))) ;

            vval =  _mm_set1_pd(v[column_index[j]*stride+2]) ;

            _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+6), vval))) ;
            _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+7), vval))) ;
            _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+8), vval))) ;
            _mm_store_pd((dest_iterator),  _mm_add_pd( _mm_load_pd((dest_iterator)), _mm_mul_pd(*(array_iterator+6), vval))) ;
            _mm_store_pd((dest_iterator+2),  _mm_add_pd( _mm_load_pd((dest_iterator+2)), _mm_mul_pd(*(array_iterator+7), vval))) ;
            _mm_store_pd((dest_iterator+4),  _mm_add_pd( _mm_load_pd((dest_iterator+4)), _mm_mul_pd(*(array_iterator+8), vval))) ;

            vval =  _mm_set1_pd(v[column_index[j]*stride+3]) ;

            _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+9), vval))) ;
            _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+10), vval))) ;
            _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+11), vval))) ;
            _mm_store_pd((dest_iterator),  _mm_add_pd( _mm_load_pd((dest_iterator)), _mm_mul_pd(*(array_iterator+9), vval))) ;
            _mm_store_pd((dest_iterator+2),  _mm_add_pd( _mm_load_pd((dest_iterator+2)), _mm_mul_pd(*(array_iterator+10), vval))) ;
            _mm_store_pd((dest_iterator+4),  _mm_add_pd( _mm_load_pd((dest_iterator+4)), _mm_mul_pd(*(array_iterator+11), vval))) ;

            vval =  _mm_set1_pd(v[column_index[j]*stride+4]) ;

            _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+12), vval))) ;
            _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+13), vval))) ;
            _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+14), vval))) ;
            _mm_store_pd((dest_iterator),  _mm_add_pd( _mm_load_pd((dest_iterator)), _mm_mul_pd(*(array_iterator+12), vval))) ;
            _mm_store_pd((dest_iterator+2),  _mm_add_pd( _mm_load_pd((dest_iterator+2)), _mm_mul_pd(*(array_iterator+13), vval))) ;
            _mm_store_pd((dest_iterator+4),  _mm_add_pd( _mm_load_pd((dest_iterator+4)), _mm_mul_pd(*(array_iterator+14), vval))) ;

            vval =  _mm_set1_pd(v[column_index[j]*stride+5]) ;

            _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+15), vval))) ;
            _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+16), vval))) ;
            _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+17), vval))) ;
            _mm_store_pd((dest_iterator),  _mm_add_pd( _mm_load_pd((dest_iterator)), _mm_mul_pd(*(array_iterator+15), vval))) ;
            _mm_store_pd((dest_iterator+2),  _mm_add_pd( _mm_load_pd((dest_iterator+2)), _mm_mul_pd(*(array_iterator+16), vval))) ;
            _mm_store_pd((dest_iterator+4),  _mm_add_pd( _mm_load_pd((dest_iterator+4)), _mm_mul_pd(*(array_iterator+17), vval))) ;
            array_iterator +=18 ;
            dest_iterator += (column_index[j+1]-column_index[j])*stride ;
        }
#else
        const double * array_iterator0 = &array[start*colLength*stride] ;
        const double * array_iterator1 = &array[start*colLength*stride+1] ;
        double * dest_iterator  = dest + (column_index[start+1]-column_index[start])*stride;
        double vval =  v[column_index[start]*stride] ;

        *(dest) += *array_iterator0 * vval ;
        *(dest+1) += *array_iterator1 * vval ;

        *(dest+2) += *(array_iterator0+2) * vval ;
        *(dest+3) += *(array_iterator1+2) * vval ;

        *(dest+4) += *(array_iterator0+4) * vval ;
        *(dest+5) += *(array_iterator1+4) * vval ;

        vval =  v[column_index[start]*stride+1] ;

        *(dest) += *(array_iterator0+6) * vval ;
        *(dest+1) += *(array_iterator1+6) * vval ;

        *(dest+2) += *(array_iterator0+8) * vval ;
        *(dest+3) += *(array_iterator1+8) * vval ;

        *(dest+4) += *(array_iterator0+10) * vval ;
        *(dest+5) += *(array_iterator1+10) * vval ;

        vval =  v[column_index[start]*stride+2] ;

        *(dest) += *(array_iterator0+12) * vval ;
        *(dest+1) += *(array_iterator1+12) * vval ;

        *(dest+2) += *(array_iterator0+14) * vval ;
        *(dest+3) += *(array_iterator1+14) * vval ;

        *(dest+4) += *(array_iterator0+16) * vval ;
        *(dest+5) += *(array_iterator1+16) * vval ;

        vval =  v[column_index[start]*stride+3] ;

        *(dest) += *(array_iterator0+18) * vval ;
        *(dest+1) += *(array_iterator1+18) * vval ;

        *(dest+2) += *(array_iterator0+20) * vval ;
        *(dest+3) += *(array_iterator1+20) * vval ;

        *(dest+4) += *(array_iterator0+22) * vval ;
        *(dest+5) += *(array_iterator1+22) * vval ;

        vval =  v[column_index[start]*stride+4] ;

        *(dest) += *(array_iterator0+24) * vval ;
        *(dest+1) += *(array_iterator1+24) * vval ;

        *(dest+2) += *(array_iterator0+26) * vval ;
        *(dest+3) += *(array_iterator1+26) * vval ;

        *(dest+4) += *(array_iterator0+28) * vval ;
        *(dest+5) += *(array_iterator1+28) * vval ;

        vval =  v[column_index[start]*stride+5] ;

        *(dest) += *(array_iterator0+30) * vval ;
        *(dest+1) += *(array_iterator1+30) * vval ;

        *(dest+2) += *(array_iterator0+32) * vval ;
        *(dest+3) += *(array_iterator1+32) * vval ;

        *(dest+4) += *(array_iterator0+34) * vval ;
        *(dest+5) += *(array_iterator1+34) * vval ;
        array_iterator0 += 36 ;
        array_iterator1 += 36 ;
        for(unsigned int j = start+1 ; j != length+start ; j++)
        {
            double vval =  v[column_index[j]*stride] ;

            *(dest) += *array_iterator0 * vval ;
            *(dest+1) += *array_iterator1 * vval ;

            *(dest+2) += *(array_iterator0+2) * vval ;
            *(dest+3) += *(array_iterator1+2) * vval ;

            *(dest+4) += *(array_iterator0+4) * vval ;
            *(dest+5) += *(array_iterator1+4) * vval ;

            *(dest_iterator) += *array_iterator0 * vval ;
            *(dest_iterator+1) += *array_iterator1 * vval ;

            *(dest_iterator+2) += *(array_iterator0+2) * vval ;
            *(dest_iterator+3) += *(array_iterator1+2) * vval ;

            *(dest_iterator+4) += *(array_iterator0+4) * vval ;
            *(dest_iterator+5) += *(array_iterator1+4) * vval ;

            vval =  v[column_index[j]*stride+1] ;

            *(dest) += *(array_iterator0+6) * vval ;
            *(dest+1) += *(array_iterator1+6) * vval ;

            *(dest+2) += *(array_iterator0+8) * vval ;
            *(dest+3) += *(array_iterator1+8) * vval ;

            *(dest+4) += *(array_iterator0+10) * vval ;
            *(dest+5) += *(array_iterator1+10) * vval ;

            *(dest_iterator) += *(array_iterator0+6) * vval ;
            *(dest_iterator+1) += *(array_iterator1+6) * vval ;

            *(dest_iterator+2) += *(array_iterator0+8) * vval ;
            *(dest_iterator+3) += *(array_iterator1+8) * vval ;

            *(dest_iterator+4) += *(array_iterator0+10) * vval ;
            *(dest_iterator+5) += *(array_iterator1+10) * vval ;

            vval =  v[column_index[j]*stride+2] ;

            *(dest) += *(array_iterator0+12) * vval ;
            *(dest+1) += *(array_iterator1+12) * vval ;

            *(dest+2) += *(array_iterator0+14) * vval ;
            *(dest+3) += *(array_iterator1+14) * vval ;

            *(dest+4) += *(array_iterator0+16) * vval ;
            *(dest+5) += *(array_iterator1+16) * vval ;

            *(dest_iterator) += *(array_iterator0+12) * vval ;
            *(dest_iterator+1) += *(array_iterator1+12) * vval ;

            *(dest_iterator+2) += *(array_iterator0+14) * vval ;
            *(dest_iterator+3) += *(array_iterator1+14) * vval ;

            *(dest_iterator+4) += *(array_iterator0+16) * vval ;
            *(dest_iterator+5) += *(array_iterator1+16) * vval ;

            vval =  v[column_index[j]*stride+3] ;

            *(dest) += *(array_iterator0+18) * vval ;
            *(dest+1) += *(array_iterator1+18) * vval ;

            *(dest+2) += *(array_iterator0+20) * vval ;
            *(dest+3) += *(array_iterator1+20) * vval ;

            *(dest+4) += *(array_iterator0+22) * vval ;
            *(dest+5) += *(array_iterator1+22) * vval ;

            vval =  v[column_index[j]*stride+4] ;

            *(dest) += *(array_iterator0+24) * vval ;
            *(dest+1) += *(array_iterator1+24) * vval ;

            *(dest+2) += *(array_iterator0+26) * vval ;
            *(dest+3) += *(array_iterator1+26) * vval ;

            *(dest+4) += *(array_iterator0+28) * vval ;
            *(dest+5) += *(array_iterator1+28) * vval ;

            *(dest_iterator) += *(array_iterator0+24) * vval ;
            *(dest_iterator+1) += *(array_iterator1+24) * vval ;

            *(dest_iterator+2) += *(array_iterator0+26) * vval ;
            *(dest_iterator+3) += *(array_iterator1+26) * vval ;

            *(dest_iterator+4) += *(array_iterator0+28) * vval ;
            *(dest_iterator+5) += *(array_iterator1+28) * vval ;

            vval =  v[column_index[j]*stride+5] ;

            *(dest) += *(array_iterator0+30) * vval ;
            *(dest+1) += *(array_iterator1+30) * vval ;

            *(dest+2) += *(array_iterator0+32) * vval ;
            *(dest+3) += *(array_iterator1+32) * vval ;

            *(dest+4) += *(array_iterator0+34) * vval ;
            *(dest+5) += *(array_iterator1+34) * vval ;

            *(dest_iterator) += *(array_iterator0+30) * vval ;
            *(dest_iterator+1) += *(array_iterator1+30) * vval ;

            *(dest_iterator+2) += *(array_iterator0+32) * vval ;
            *(dest_iterator+3) += *(array_iterator1+32) * vval ;

            *(dest_iterator+4) += *(array_iterator0+34) * vval ;
            *(dest_iterator+5) += *(array_iterator1+34) * vval ;
            array_iterator0 += 36 ;
            array_iterator1 += 36 ;
            dest_iterator += (column_index[j+1]-column_index[j])*stride ;
        }
#endif
        return ;
    }

    const int colLength = stride + stride%2 ;
#ifdef HAVE_SSE3
    const __m128d * array_iterator = (__m128d*)&array[start*colLength*stride] ;
    __m128d * dest_iterator = (__m128d*)(dest +(column_index[start+1]-column_index[start])*stride);
    for(size_t c = 0 ; c != stride ; ++c)
    {
        const __m128d vval =  _mm_set1_pd(v[column_index[start]*stride+c]) ;
        for(int i = 0 ; i != colLength/2 ; ++i,++array_iterator)
        {
            _mm_store_pd((dest+i*2),  _mm_add_pd( _mm_load_pd((dest+i*2)), _mm_mul_pd(*array_iterator, vval))) ;
        }
    }

    for(unsigned int j = start+1 ; j != length+start ; j++)
    {
        for(size_t c = 0 ; c != stride ; ++c)
        {
            const __m128d vval =  _mm_set1_pd(v[column_index[j]*stride+c]) ;
            for(int i = 0 ; i != colLength/2 ; ++i,++array_iterator)
            {
                __m128d temp = _mm_mul_pd(*array_iterator, vval) ;
                _mm_store_pd((dest+i*2),  _mm_add_pd( _mm_load_pd((dest+i*2)), temp)) ;
                _mm_store_pd((dest_iterator+i*2),  _mm_add_pd( _mm_load_pd((dest+i*2)), temp)) ;
            }
        }
        dest_iterator += (column_index[j+1]-column_index[j])*stride ;
    }
#else
    const double * array_iterator0 = &array[start*colLength*stride] ;
    const double * array_iterator1 = &array[start*colLength*stride+1] ;
    double * dest_iterator = dest +(column_index[start+1]-column_index[start])*stride;

    for(size_t c = 0 ; c < stride ; c++)
    {
        const double vval =  v[column_index[start]*stride+c] ;
        for(int i = 0 ; i != colLength ; i+=2)
        {
            *(dest+i) += *array_iterator0 * vval ;
            *(dest+i+1) += *array_iterator1 * vval ;
            array_iterator0 += 2 ;
            array_iterator1 += 2 ;
        }
    }

    for(unsigned int j = start+1 ; j != length+start ; j++)
    {
        for(size_t c = 0 ; c < stride ; c++)
        {
            const double vval =  v[column_index[j]*stride+c] ;
            for(int i = 0 ; i != colLength ; i+=2)
            {
                *(dest+i) += *array_iterator0 * vval ;
                *(dest+i+1) += *array_iterator1 * vval ;
                *(dest_iterator+i) += *array_iterator0 * vval ;
                *(dest_iterator+i+1) += *array_iterator1 * vval ;
                array_iterator0 += 2 ;
                array_iterator1 += 2 ;
            }
        }
        dest_iterator += (column_index[j+1]-column_index[j])*stride;
    }
#endif
}

void ConstSparseVector::inner_product(const Vector &v, double *dest, const size_t rowstart , size_t const colstart ) const
{
    size_t mstart = start ;
    if(colstart)
    {
        const unsigned int * __start__       = &column_index[start] ;
        const unsigned int * __end__         = &column_index[start+length] ;
        const unsigned int * i_index_pointer = std::lower_bound(__start__, __end__, colstart/stride) ;
        const unsigned int offset            = i_index_pointer - __start__ ;

        mstart = (start <= start+offset && start+offset < start+length) ? (start+offset) : start+length;
    }
// 		mstart = start ;std::max(mstart,start) ;
//
//
    switch(stride)
    {
    case 1:
    {
        for(unsigned int j = mstart ; j < length+start ; j++)
        {
            *dest += v[column_index[j]]*array[j*2] ;
        }
        return ;
    }
    case 2:
    {
        const double * array_iterator = &array[mstart*2*2] ;
        const double * vec_iterator = &v[column_index[mstart]*2] ;
        for(unsigned int j = mstart ; j < length+start ; j++ )
        {
            *dest     += *array_iterator*(*vec_iterator);
            *(dest+1) += *(array_iterator+1)*(*vec_iterator) ;
            *dest     += *(array_iterator+2)*(*(vec_iterator +1)) ;
            *(dest+1) += *(array_iterator+3)*(*(vec_iterator +1)) ;
            array_iterator+=4 ;
            if(j+1 < column_index.size())
                vec_iterator += column_index[j+1]*stride-column_index[j]*stride ;
        }
        return ;
    }

    case 3:
    {
        const double * array_iterator = &array[mstart*3*4] ;
        const double * vec_iterator = &v[column_index[mstart]*3] ;
        for(unsigned int j = mstart ; j < length+start ; j++)
        {
            *dest     += *array_iterator     * (*vec_iterator);
            *(dest+1) += *(array_iterator+1) * (*vec_iterator);
            *(dest+2) += *(array_iterator+2) * (*vec_iterator);
            *dest     += *(array_iterator+4) * (*(vec_iterator + 1)) ;
            *(dest+1) += *(array_iterator+5) * (*(vec_iterator + 1)) ;
            *(dest+2) += *(array_iterator+6) * (*(vec_iterator + 1)) ;
            *dest     += *(array_iterator+8 ) * (*(vec_iterator + 2));
            *(dest+1) += *(array_iterator+9 ) * (*(vec_iterator + 2));
            *(dest+2) += *(array_iterator+10) * (*(vec_iterator + 2));
            array_iterator+=12 ;
            if(j+1 < column_index.size())
                vec_iterator += column_index[j+1]*stride-column_index[j]*stride ;
        }
        return ;
    }
    case 6 :
    {
#ifdef HAVE_SSE3
        const __m128d * array_iterator = (__m128d*)&array[mstart*36] ;
// 			#pragma omp parallel for schedule(runtime)
        for(unsigned int j = mstart ; j < length+start ; j++)
        {
            __m128d vval =  _mm_set1_pd(v[column_index[j]*stride]) ;

            _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*array_iterator, vval))) ;
            _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+1), vval))) ;
            _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+2), vval))) ;

            vval =  _mm_set1_pd(v[column_index[j]*stride+1]) ;

            _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+3,) vval))) ;
            _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+4), vval))) ;
            _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+5), vval))) ;

            vval =  _mm_set1_pd(v[column_index[j]*stride+2]) ;

            _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+6), vval))) ;
            _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+7), vval))) ;
            _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+8), vval))) ;

            vval =  _mm_set1_pd(v[column_index[j]*stride+3]) ;

            _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+9), vval))) ;
            _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+10), vval))) ;
            _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+11), vval))) ;

            vval =  _mm_set1_pd(v[column_index[j]*stride+4]) ;

            _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+12), vval))) ;
            _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+13), vval))) ;
            _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+14), vval))) ;

            vval =  _mm_set1_pd(v[column_index[j]*stride+5]) ;

            _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+15), vval))) ;
            _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+16), vval))) ;
            _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+17), vval))) ;
            array_iterator +=18 ;
        }
#else
        const double * array_iterator0 = &array[mstart*36] ;
        const double * array_iterator1 = &array[mstart*36+1] ;
// 			#pragma omp parallel for schedule(runtime)
        for(unsigned int j = mstart ; j < length+start ; j++)
        {
            double vval =  v[column_index[j]*6] ;

            *(dest) += *array_iterator0 * vval ;
            *(dest+1) += *array_iterator1 * vval ;

            *(dest+2) += *(array_iterator0+2) * vval ;
            *(dest+3) += *(array_iterator1+2) * vval ;

            *(dest+4) += *(array_iterator0+4) * vval ;
            *(dest+5) += *(array_iterator1+4) * vval ;

            vval =  v[column_index[j]*6+1] ;

            *(dest) += *(array_iterator0+6) * vval ;
            *(dest+1) += *(array_iterator1+6) * vval ;

            *(dest+2) += *(array_iterator0+8) * vval ;
            *(dest+3) += *(array_iterator1+8) * vval ;

            *(dest+4) += *(array_iterator0+10) * vval ;
            *(dest+5) += *(array_iterator1+10) * vval ;

            vval =  v[column_index[j]*6+2] ;

            *(dest) += *(array_iterator0+12) * vval ;
            *(dest+1) += *(array_iterator1+12) * vval ;

            *(dest+2) += *(array_iterator0+14) * vval ;
            *(dest+3) += *(array_iterator1+14) * vval ;

            *(dest+4) += *(array_iterator0+16) * vval ;
            *(dest+5) += *(array_iterator1+16) * vval ;

            vval =  v[column_index[j]*6+3] ;

            *(dest) += *(array_iterator0+18) * vval ;
            *(dest+1) += *(array_iterator1+18) * vval ;

            *(dest+2) += *(array_iterator0+20) * vval ;
            *(dest+3) += *(array_iterator1+20) * vval ;

            *(dest+4) += *(array_iterator0+22) * vval ;
            *(dest+5) += *(array_iterator1+22) * vval ;

            vval =  v[column_index[j]*6+4] ;

            *(dest) += *(array_iterator0+24) * vval ;
            *(dest+1) += *(array_iterator1+24) * vval ;

            *(dest+2) += *(array_iterator0+26) * vval ;
            *(dest+3) += *(array_iterator1+26) * vval ;

            *(dest+4) += *(array_iterator0+28) * vval ;
            *(dest+5) += *(array_iterator1+28) * vval ;

            vval =  v[column_index[j]*6+5] ;

            *(dest) += *(array_iterator0+30) * vval ;
            *(dest+1) += *(array_iterator1+30) * vval ;

            *(dest+2) += *(array_iterator0+32) * vval ;
            *(dest+3) += *(array_iterator1+32) * vval ;

            *(dest+4) += *(array_iterator0+34) * vval ;
            *(dest+5) += *(array_iterator1+34) * vval ;
            array_iterator0 += 36 ;
            array_iterator1 += 36 ;
        }
#endif
        return ;
    }
    default:
    {
        unsigned int colLength = stride+stride%2 ;
        const double * array_iterator0 = &array[mstart*colLength*stride] ;
        const double * array_iterator1 = &array[mstart*colLength*stride+1] ;
        for(unsigned int j = mstart ; j < length+start ; j++)
        {
            for(size_t c = 0 ; c < stride ; c++)
            {
                const double vval =  v[column_index[j]*stride+c] ;
                for(size_t i = 0 ; i != colLength ; i+=2)
                {
                    *(dest+i) += *array_iterator0 * vval ;
                    *(dest+i+1) += *array_iterator1 * vval ;
                    array_iterator0 += 2 ;
                    array_iterator1 += 2 ;
                }
            }
        }
        return ;
    }
    }
}

Vector ConstSparseVector::operator *(const Vector &v) const
{
    if(stride == 2)
    {
#ifdef HAVE_SSE3
        Vector ret(0., 2) ;
        const __m128d * array_iterator = (__m128d*)&array[start*4] ;
        const double * vec_iterator = &v[column_index[start]*2] ;
        for(unsigned int j = start ; j != length+start && j+1 != length+start; j++)
        {

            _mm_store_pd(&ret[0],  _mm_add_pd( _mm_load_pd(&ret[0]), _mm_add_pd( _mm_mul_pd(*array_iterator,  _mm_set1_pd(*vec_iterator)), _mm_mul_pd(*(array_iterator+1), _mm_set1_pd(*(vec_iterator+1)))))) ;
            array_iterator+=2 ;
            vec_iterator += column_index[j+1]*2-column_index[j]*2 ;

        }
#else
        Vector ret(0., 2) ;
        const double * array_iterator0 = &array[start*4] ;
        const double * array_iterator1 = &array[start*4+1] ;
        const double * array_iterator2 = &array[start*4+2] ;
        const double * array_iterator3 = &array[start*4+3] ;
        for(unsigned int j = start ; j != length+start ; j++)
        {
            ret[0] += *array_iterator0*v[column_index[j]*2]+*array_iterator2*v[column_index[j]*2+1] ;
            ret[1] += *array_iterator1*v[column_index[j]*2]+*array_iterator3*v[column_index[j]*2+1] ;
            array_iterator0+=4 ;
            array_iterator1+=4 ;
            array_iterator2+=4 ;
            array_iterator3+=4 ;
        }
#endif
        return ret ;
    }

    const int colLength = stride + stride%2 ;
    Vector ret(0., colLength) ;
#ifdef HAVE_SSE3
    const __m128d * array_iterator = (__m128d*)&array[start*colLength*stride] ;
    for(unsigned int j = start ; j != length+start ; j++)
    {
        for(size_t c = 0 ; c < stride ; c++)
        {
            const __m128d vval =  _mm_set1_pd(v[column_index[j]*stride+c]) ;
            for(int i = 0 ; i != colLength ; i+=2)
            {
                _mm_store_pd(&ret[i],  _mm_add_pd( _mm_load_pd(&ret[i]), _mm_mul_pd(*array_iterator, vval))) ;
                array_iterator++ ;
            }
        }
    }
#else
    const double * array_iterator0 = &array[start*colLength*stride] ;
    const double * array_iterator1 = &array[start*colLength*stride+1] ;
    for(unsigned int j = start ; j != length+start ; j++)
    {
        for(size_t c = 0 ; c < stride ; c++)
        {
            const double vval =  v[column_index[j]*stride+c] ;
            for(int i = 0 ; i != colLength ; i+=2)
            {
                ret[i] += *array_iterator0 * vval ;
                ret[i+1] += *array_iterator1 * vval ;
                array_iterator0 += 2 ;
                array_iterator1 += 2 ;
            }
        }
    }
#endif
    return ret ;

}

Matrix ConstSparseVector::operator *(const SparseVector& v) const
{
    Matrix ret(stride, stride) ;
    int colLength = stride+stride%2 ;
    int blocksize = stride*colLength ;

    size_t i = 0 ;
    size_t j = 0 ;
    while(i < length && j < v.length)
    {
        if(column_index[start+i] > v.idx[v.start+ j])
            j++ ;
        else if(column_index[start+i] < v.idx[v.start+ j])
            i++ ;
        else
        {
            int idx = 0 ;
            for(size_t k = 0 ; k < stride ; k++)
            {
                for(size_t l = 0 ; l < stride ; l++)
                {
                    ret[k][l] += array[start*blocksize+i*blocksize+colLength*k+l]*v.val[v.start*blocksize+j*blocksize + idx++] ;
                }
            }
            i++ ;
            j++ ;
        }
    }
    return ret ;

}

Matrix ConstSparseVector::operator *(const ConstSparseVector& v) const
{
    Matrix ret(stride, stride) ;
    int colLength = stride+stride%2 ;
    int blocksize = stride*colLength ;

    size_t i = 0 ;
    size_t j = 0 ;
    while(i < length && j < v.length)
    {
        if(column_index[start+i] > v.column_index[v.start+ j])
            j++ ;
        else if(column_index[start+i] < v.column_index[v.start+ j])
            i++ ;
        else
        {
            int idx = 0 ;
            for(size_t k = 0 ; k < stride ; k++)
            {
                for(size_t l = 0 ; l < stride ; l++)
                {
                    ret[k][l] += array[start*blocksize+i*blocksize+colLength*k+l]*v.array[v.start*blocksize+j*blocksize + idx++] ;
                }
            }
            i++ ;
            j++ ;
        }
    }
    return ret ;

}

void ConstSparseVector::print() const
{
    for(size_t j = 0 ; j < length ; j++)
        std::cout << column_index[start + j] << " -> " << array[start + j] << std::endl ;
}

Vector ConstSparseVector::operator +(const Vector& v) const
{
    Vector ret(v) ;

    for(size_t j = start ; j < length+start ; j++)
    {
        size_t index = column_index[j] ;
        ret[index] += array[j] ;
    }

    return ret ;

}

}

