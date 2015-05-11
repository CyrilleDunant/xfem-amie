//
// C++ Interface: sparse_vector
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//


#ifndef __SPARSE_VECTOR_H
#define __SPARSE_VECTOR_H

#include <valarray>
#include <iostream>
#include "../utilities/matrixops.h"

namespace Amie
{

/** \brief Sparse Vector implementation.
 * The sparse vector implementation is not intended to be used as a standalone class. Instead it behaves as a row of a
 * CoordinateIndexedSparseMatrix. As such a matrix is stored in blocks, the sparse vector holds references to a group
 * of rows, and an aditional member is used to indicate which should be used for access.
 *
 * This class is designed to allow modification of members.
 */
struct SparseVector
{
public:
    Vector & val ;
    std::valarray<unsigned int> & idx ;
    const size_t length ;
    const size_t start ;
    const size_t stride ;
    const size_t index ;

    double zero ;

public:
    /** \brief SparseVector constructor. Initialises the references to the data.
     @param v array containing the values
     @param idx column index of the blocks
     @param l number of blocks in the row.
     @param s block start index of the values relevant to this row
     @param index row number (actual number, not in blocks)
     @param st stride: block size
     */
    SparseVector(Vector & v, std::valarray<unsigned int> & idx , const size_t l , const size_t s, const size_t index, const size_t st) ;

    /** \brief access a value in the row
     * Return the value in the ith column of the matrix on this row. If this value is not stored, 0 is retured.
     */
    double operator [](const size_t i) const ;
    
    double * getPointer(const size_t i)  ;

    /** \brief access a value in the row
     * Return the value in the ith column of the matrix on this row. If this value is not stored, assignement has no effect.
     */
    double & operator [](const size_t) ;

    /** \brief simultaneously compute a number of dot products equal to the block size.
     */
    Vector operator *(const Vector&) const ;

    void parallel_product(const Vector &v, double *dest, const size_t rowstart = 0, const size_t colstart = 0) const
    {
        if(stride == 2)
        {
#ifdef HAVE_SSE3
            const __m128d * array_iterator = (__m128d*)&val[start*4] ;
            const double * vec_iterator = &v[idx[start]*2] ;
            const double * dest_iterator = dest+(idx[start+1]-idx[start])*2 ;
            _mm_store_pd(dest,  _mm_add_pd( _mm_load_pd(dest), _mm_add_pd( _mm_mul_pd(*array_iterator,  _mm_set1_pd(*vec_iterator)), _mm_mul_pd(*(array_iterator+1), _mm_set1_pd(*(vec_iterator+1)))))) ;
            for(unsigned int j = start+1 ; j != length+start && j+1 != length+start; vec_iterator += (idx[j+1]-idx[j])*2, dest_iterator += (idx[j+1]-idx[j])*2,++j,array_iterator+=2)
            {
                __m128d tmp =  _mm_add_pd( _mm_load_pd(*dest_iterator), _mm_add_pd( _mm_mul_pd(*array_iterator,  _mm_set1_pd(*vec_iterator)), _mm_mul_pd(*(array_iterator+1), _mm_set1_pd(*(vec_iterator+1))))) ;
                _mm_store_pd(dest, _mm_add_pd(tmp, dest)) ;
                _mm_store_pd(dest_iterator, _mm_add_pd(tmp, dest_iterator))
            }
#else
            const double * array_iterator0 = &val[start*4] ;
            const double * array_iterator1 = &val[start*4+1] ;
            const double * array_iterator2 = &val[start*4+2] ;
            const double * array_iterator3 = &val[start*4+3] ;
            double * dest_iterator = dest + (idx[start+1]-idx[start])*2;
            *dest_iterator += *array_iterator0*v[idx[start]*2]+*array_iterator2*v[idx[start]*2+1] ;
            *(dest_iterator+1) += *array_iterator1*v[idx[start]*2]+*array_iterator3*v[idx[start]*2+1] ;
            double tmp0 = *array_iterator0*v[idx[start]*2]+*array_iterator2*v[idx[start]*2+1] ;
            double tmp1 = *array_iterator1*v[idx[start]*2]+*array_iterator3*v[idx[start]*2+1] ;
            *dest += tmp0 ;
            *(dest+1) += tmp1 ;
            for(unsigned int j = start+1 ; j != length+start ; j++,array_iterator0+=4,array_iterator1+=4,array_iterator2+=4,array_iterator3+=4, dest_iterator +=(idx[j+1]-idx[j])*2)
            {
                double tmp0 = *array_iterator0*v[idx[j]*2]+*array_iterator2*v[idx[j]*2+1] ;
                double tmp1 = *array_iterator1*v[idx[j]*2]+*array_iterator3*v[idx[j]*2+1] ;
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
            const __m128d * array_iterator = (__m128d*)&val[start*4*stride] ;
            const double * dest_iterator = dest + (idx[start+1]-idx[start])*3;
            const __m128d vval =  _mm_set1_pd(v[idx[j]*stride]) ;
            _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*array_iterator, vval))) ;
            _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*array_iterator, vval))) ;
            vval =  _mm_set1_pd(v[idx[j]*stride+1]) ;
            _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+2), vval))) ;
            _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+2), vval))) ;
            vval =  _mm_set1_pd(v[idx[j]*stride+2]) ;
            _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+4), vval))) ;
            _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+4), vval))) ;
            for(unsigned int j = start+1 ; j != length+start ; ++j,array_iterator+=6, dest_iterator += (idx[j+1]-idx[j])*3)
            {
                const __m128d vval =  _mm_set1_pd(v[idx[j]*stride]) ;
                _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*array_iterator, vval))) ;
                _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*array_iterator, vval))) ;
                _mm_store_pd((dest_iterator),  _mm_add_pd( _mm_load_pd((dest_iterator)), _mm_mul_pd(*array_iterator, vval))) ;
                _mm_store_pd((dest_iterator+2),  _mm_add_pd( _mm_load_pd((dest_iterator+2)), _mm_mul_pd(*array_iterator, vval))) ;
                vval =  _mm_set1_pd(v[idx[j]*stride+1]) ;
                _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+2), vval))) ;
                _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+2), vval))) ;
                _mm_store_pd((dest_iterator),  _mm_add_pd( _mm_load_pd((dest_iterator)), _mm_mul_pd(*(array_iterator+2), vval))) ;
                _mm_store_pd((dest_iterator+2),  _mm_add_pd( _mm_load_pd((dest_iterator+2)), _mm_mul_pd(*(array_iterator+2), vval))) ;
                vval =  _mm_set1_pd(v[idx[j]*stride+2]) ;
                _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+4), vval))) ;
                _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+4), vval))) ;
                _mm_store_pd((dest_iterator),  _mm_add_pd( _mm_load_pd((dest_iterator)), _mm_mul_pd(*(array_iterator+4), vval))) ;
                _mm_store_pd((dest_iterator+2),  _mm_add_pd( _mm_load_pd((dest_iterator+2)), _mm_mul_pd(*(array_iterator+4), vval))) ;
            }
#else
            const double * array_iterator0 = &val[start*4*stride] ;
            const double * array_iterator1 = &val[start*4*stride+1] ;
            double * dest_iterator = dest+ (idx[start+1]-idx[start])*stride;
            double vval =  v[idx[start]*stride] ;
            *(dest+0) += *array_iterator0 * vval ;
            *(dest+1) += *array_iterator1 * vval ;
            *(dest+2) += *(array_iterator0+2) * vval ;
            *(dest+3) += *(array_iterator1+2) * vval ;
            vval =  v[idx[start]*stride+1] ;
            *(dest+0) += *(array_iterator0+4) * vval ;
            *(dest+1) += *(array_iterator1+4) * vval ;
            *(dest+2) += *(array_iterator0+6) * vval ;
            *(dest+3) += *(array_iterator1+6) * vval ;
            vval =  v[idx[start]*stride+2] ;
            *(dest+0) += *(array_iterator0+8) * vval ;
            *(dest+1) += *(array_iterator1+8) * vval ;
            *(dest+2) += *(array_iterator0+10) * vval ;
            *(dest+3) += *(array_iterator1+10) * vval ;
            for(unsigned int j = start+1 ; j != length+start ; ++j,dest_iterator += (idx[j+1]-idx[j])*stride)
            {
                double vval =  v[idx[j]*stride] ;
                *(dest+0) += *array_iterator0 * vval ;
                *(dest+1) += *array_iterator1 * vval ;
                *(dest+2) += *(array_iterator0+2) * vval ;
                *(dest+3) += *(array_iterator1+2) * vval ;
                *(dest_iterator+0) += *array_iterator0 * vval ;
                *(dest_iterator+1) += *array_iterator1 * vval ;
                *(dest_iterator+2) += *(array_iterator0+2) * vval ;
                *(dest_iterator+3) += *(array_iterator1+2) * vval ;
                vval =  v[idx[j]*stride+1] ;
                *(dest+0) += *(array_iterator0+4) * vval ;
                *(dest+1) += *(array_iterator1+4) * vval ;
                *(dest+2) += *(array_iterator0+6) * vval ;
                *(dest+3) += *(array_iterator1+6) * vval ;
                *(dest_iterator+0) += *(array_iterator0+4) * vval ;
                *(dest_iterator+1) += *(array_iterator1+4) * vval ;
                *(dest_iterator+2) += *(array_iterator0+6) * vval ;
                *(dest_iterator+3) += *(array_iterator1+6) * vval ;
                vval =  v[idx[j]*stride+2] ;
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
            const __m128d * array_iterator = (__m128d*)&val[start*colLength*stride] ;
            const __m128d * dest_iterator = (__m128d*)(dest + (idx[start+1]-idx[start])*stride);

            __m128d vval =  _mm_set1_pd(v[idx[j]*stride]) ;

            _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*array_iterator, vval))) ;
            _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+1), vval))) ;
            _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+2), vval))) ;

            vval =  _mm_set1_pd(v[idx[j]*stride+1]) ;

            _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+3,) vval))) ;
            _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+4), vval))) ;
            _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+5), vval))) ;

            vval =  _mm_set1_pd(v[idx[j]*stride+2]) ;

            _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+6), vval))) ;
            _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+7), vval))) ;
            _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+8), vval))) ;

            vval =  _mm_set1_pd(v[idx[j]*stride+3]) ;

            _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+9), vval))) ;
            _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+10), vval))) ;
            _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+11), vval))) ;

            vval =  _mm_set1_pd(v[idx[j]*stride+4]) ;

            _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+12), vval))) ;
            _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+13), vval))) ;
            _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+14), vval))) ;

            vval =  _mm_set1_pd(v[idx[j]*stride+5]) ;

            _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+15), vval))) ;
            _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+16), vval))) ;
            _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+17), vval))) ;
            array_iterator +=18 ;
            for(unsigned int j = start ; j != length+start ; j++)
            {
                __m128d vval =  _mm_set1_pd(v[idx[j]*stride]) ;

                _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*array_iterator, vval))) ;
                _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+1), vval))) ;
                _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+2), vval))) ;
                _mm_store_pd((dest_iterator),  _mm_add_pd( _mm_load_pd((dest_iterator)), _mm_mul_pd(*array_iterator, vval))) ;
                _mm_store_pd((dest_iterator+2),  _mm_add_pd( _mm_load_pd((dest_iterator+2)), _mm_mul_pd(*(array_iterator+1), vval))) ;
                _mm_store_pd((dest_iterator+4),  _mm_add_pd( _mm_load_pd((dest_iterator+4)), _mm_mul_pd(*(array_iterator+2), vval))) ;

                vval =  _mm_set1_pd(v[idx[j]*stride+1]) ;

                _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+3,) vval))) ;
                _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+4), vval))) ;
                _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+5), vval))) ;
                _mm_store_pd((dest_iterator),  _mm_add_pd( _mm_load_pd((dest_iterator)), _mm_mul_pd(*(array_iterator+3), vval))) ;
                _mm_store_pd((dest_iterator+2),  _mm_add_pd( _mm_load_pd((dest_iterator+2)), _mm_mul_pd(*(array_iterator+4), vval))) ;
                _mm_store_pd((dest_iterator+4),  _mm_add_pd( _mm_load_pd((dest_iterator+4)), _mm_mul_pd(*(array_iterator+5), vval))) ;

                vval =  _mm_set1_pd(v[idx[j]*stride+2]) ;

                _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+6), vval))) ;
                _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+7), vval))) ;
                _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+8), vval))) ;
                _mm_store_pd((dest_iterator),  _mm_add_pd( _mm_load_pd((dest_iterator)), _mm_mul_pd(*(array_iterator+6), vval))) ;
                _mm_store_pd((dest_iterator+2),  _mm_add_pd( _mm_load_pd((dest_iterator+2)), _mm_mul_pd(*(array_iterator+7), vval))) ;
                _mm_store_pd((dest_iterator+4),  _mm_add_pd( _mm_load_pd((dest_iterator+4)), _mm_mul_pd(*(array_iterator+8), vval))) ;

                vval =  _mm_set1_pd(v[idx[j]*stride+3]) ;

                _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+9), vval))) ;
                _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+10), vval))) ;
                _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+11), vval))) ;
                _mm_store_pd((dest_iterator),  _mm_add_pd( _mm_load_pd((dest_iterator)), _mm_mul_pd(*(array_iterator+9), vval))) ;
                _mm_store_pd((dest_iterator+2),  _mm_add_pd( _mm_load_pd((dest_iterator+2)), _mm_mul_pd(*(array_iterator+10), vval))) ;
                _mm_store_pd((dest_iterator+4),  _mm_add_pd( _mm_load_pd((dest_iterator+4)), _mm_mul_pd(*(array_iterator+11), vval))) ;

                vval =  _mm_set1_pd(v[idx[j]*stride+4]) ;

                _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+12), vval))) ;
                _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+13), vval))) ;
                _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+14), vval))) ;
                _mm_store_pd((dest_iterator),  _mm_add_pd( _mm_load_pd((dest_iterator)), _mm_mul_pd(*(array_iterator+12), vval))) ;
                _mm_store_pd((dest_iterator+2),  _mm_add_pd( _mm_load_pd((dest_iterator+2)), _mm_mul_pd(*(array_iterator+13), vval))) ;
                _mm_store_pd((dest_iterator+4),  _mm_add_pd( _mm_load_pd((dest_iterator+4)), _mm_mul_pd(*(array_iterator+14), vval))) ;

                vval =  _mm_set1_pd(v[idx[j]*stride+5]) ;

                _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+15), vval))) ;
                _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+16), vval))) ;
                _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+17), vval))) ;
                _mm_store_pd((dest_iterator),  _mm_add_pd( _mm_load_pd((dest_iterator)), _mm_mul_pd(*(array_iterator+15), vval))) ;
                _mm_store_pd((dest_iterator+2),  _mm_add_pd( _mm_load_pd((dest_iterator+2)), _mm_mul_pd(*(array_iterator+16), vval))) ;
                _mm_store_pd((dest_iterator+4),  _mm_add_pd( _mm_load_pd((dest_iterator+4)), _mm_mul_pd(*(array_iterator+17), vval))) ;
                array_iterator +=18 ;
                dest_iterator += (idx[j+1]-idx[j])*stride ;
            }
#else
            const double * array_iterator0 = &val[start*colLength*stride] ;
            const double * array_iterator1 = &val[start*colLength*stride+1] ;
            double * dest_iterator  = dest + (idx[start+1]-idx[start])*stride;
            double vval =  v[idx[start]*stride] ;

            *(dest) += *array_iterator0 * vval ;
            *(dest+1) += *array_iterator1 * vval ;

            *(dest+2) += *(array_iterator0+2) * vval ;
            *(dest+3) += *(array_iterator1+2) * vval ;

            *(dest+4) += *(array_iterator0+4) * vval ;
            *(dest+5) += *(array_iterator1+4) * vval ;

            vval =  v[idx[start]*stride+1] ;

            *(dest) += *(array_iterator0+6) * vval ;
            *(dest+1) += *(array_iterator1+6) * vval ;

            *(dest+2) += *(array_iterator0+8) * vval ;
            *(dest+3) += *(array_iterator1+8) * vval ;

            *(dest+4) += *(array_iterator0+10) * vval ;
            *(dest+5) += *(array_iterator1+10) * vval ;

            vval =  v[idx[start]*stride+2] ;

            *(dest) += *(array_iterator0+12) * vval ;
            *(dest+1) += *(array_iterator1+12) * vval ;

            *(dest+2) += *(array_iterator0+14) * vval ;
            *(dest+3) += *(array_iterator1+14) * vval ;

            *(dest+4) += *(array_iterator0+16) * vval ;
            *(dest+5) += *(array_iterator1+16) * vval ;

            vval =  v[idx[start]*stride+3] ;

            *(dest) += *(array_iterator0+18) * vval ;
            *(dest+1) += *(array_iterator1+18) * vval ;

            *(dest+2) += *(array_iterator0+20) * vval ;
            *(dest+3) += *(array_iterator1+20) * vval ;

            *(dest+4) += *(array_iterator0+22) * vval ;
            *(dest+5) += *(array_iterator1+22) * vval ;

            vval =  v[idx[start]*stride+4] ;

            *(dest) += *(array_iterator0+24) * vval ;
            *(dest+1) += *(array_iterator1+24) * vval ;

            *(dest+2) += *(array_iterator0+26) * vval ;
            *(dest+3) += *(array_iterator1+26) * vval ;

            *(dest+4) += *(array_iterator0+28) * vval ;
            *(dest+5) += *(array_iterator1+28) * vval ;

            vval =  v[idx[start]*stride+5] ;

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
                double vval =  v[idx[j]*stride] ;

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

                vval =  v[idx[j]*stride+1] ;

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

                vval =  v[idx[j]*stride+2] ;

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

                vval =  v[idx[j]*stride+3] ;

                *(dest) += *(array_iterator0+18) * vval ;
                *(dest+1) += *(array_iterator1+18) * vval ;

                *(dest+2) += *(array_iterator0+20) * vval ;
                *(dest+3) += *(array_iterator1+20) * vval ;

                *(dest+4) += *(array_iterator0+22) * vval ;
                *(dest+5) += *(array_iterator1+22) * vval ;

                vval =  v[idx[j]*stride+4] ;

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

                vval =  v[idx[j]*stride+5] ;

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
                dest_iterator += (idx[j+1]-idx[j])*stride ;
            }
#endif
            return ;
        }

        const int colLength = stride + stride%2 ;
#ifdef HAVE_SSE3
        const __m128d * array_iterator = (__m128d*)&val[start*colLength*stride] ;
        __m128d * dest_iterator = (__m128d*)(dest +(idx[start+1]-idx[start])*stride);
        for(size_t c = 0 ; c != stride ; ++c)
        {
            const __m128d vval =  _mm_set1_pd(v[idx[start]*stride+c]) ;
            for(int i = 0 ; i != colLength/2 ; ++i,++array_iterator)
            {
                _mm_store_pd((dest+i*2),  _mm_add_pd( _mm_load_pd((dest+i*2)), _mm_mul_pd(*array_iterator, vval))) ;
            }
        }

        for(unsigned int j = start+1 ; j != length+start ; j++)
        {
            for(size_t c = 0 ; c != stride ; ++c)
            {
                const __m128d vval =  _mm_set1_pd(v[idx[j]*stride+c]) ;
                for(int i = 0 ; i != colLength/2 ; ++i,++array_iterator)
                {
                    __m128d temp = _mm_mul_pd(*array_iterator, vval) ;
                    _mm_store_pd((dest+i*2),  _mm_add_pd( _mm_load_pd((dest+i*2)), temp)) ;
                    _mm_store_pd((dest_iterator+i*2),  _mm_add_pd( _mm_load_pd((dest+i*2)), temp)) ;
                }
            }
            dest_iterator += (idx[j+1]-idx[j])*stride ;
        }
#else
        const double * array_iterator0 = &val[start*colLength*stride] ;
        const double * array_iterator1 = &val[start*colLength*stride+1] ;
        double * dest_iterator = dest +(idx[start+1]-idx[start])*stride;

        for(size_t c = 0 ; c < stride ; c++)
        {
            const double vval =  v[idx[start]*stride+c] ;
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
                const double vval =  v[idx[j]*stride+c] ;
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
            dest_iterator += (idx[j+1]-idx[j])*stride;
        }
#endif
    }

    void inner_product(const Vector &v, double *dest, const size_t rowstart = 0, const size_t colstart = 0) const
    {
        if(stride == 2)
        {
#ifdef HAVE_SSE3
            const __m128d * array_iterator = (__m128d*)&val[start*4] ;
            const double * vec_iterator = &v[idx[start]*2] ;
            for(unsigned int j = start ; j != length+start && j+1 != length+start; vec_iterator += (idx[j+1]-idx[j])*2,++j,array_iterator+=2)
            {
                _mm_store_pd(dest,  _mm_add_pd( _mm_load_pd(dest), _mm_add_pd( _mm_mul_pd(*array_iterator,  _mm_set1_pd(*vec_iterator)), _mm_mul_pd(*(array_iterator+1), _mm_set1_pd(*(vec_iterator+1)))))) ;
            }
#else
            const double * array_iterator0 = &val[start*4] ;
            const double * array_iterator1 = &val[start*4+1] ;
            const double * array_iterator2 = &val[start*4+2] ;
            const double * array_iterator3 = &val[start*4+3] ;
            for(unsigned int j = start ; j != length+start ; j++,array_iterator0+=4,array_iterator1+=4,array_iterator2+=4,array_iterator3+=4)
            {
                *dest += *array_iterator0*v[idx[j]*2]+*array_iterator2*v[idx[j]*2+1] ;
                *(dest+1) += *array_iterator1*v[idx[j]*2]+*array_iterator3*v[idx[j]*2+1] ;
            }
#endif
            return ;
        }

        if(stride == 3)
        {
#ifdef HAVE_SSE3
            const __m128d * array_iterator = (__m128d*)&val[start*4*stride] ;
            for(unsigned int j = start ; j != length+start ; ++j,array_iterator+=6)
            {
                const __m128d vval =  _mm_set1_pd(v[idx[j]*stride]) ;
                _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*array_iterator, vval))) ;
                _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*array_iterator, vval))) ;
                vval =  _mm_set1_pd(v[idx[j]*stride+1]) ;
                _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+2), vval))) ;
                _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+2), vval))) ;
                vval =  _mm_set1_pd(v[idx[j]*stride+2]) ;
                _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+4), vval))) ;
                _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+4), vval))) ;
            }
#else
            const double * array_iterator0 = &val[start*4*stride] ;
            const double * array_iterator1 = &val[start*4*stride+1] ;
            for(unsigned int j = start ; j != length+start ; ++j)
            {
                double vval =  v[idx[j]*stride] ;
                *(dest+0) += *array_iterator0 * vval ;
                *(dest+1) += *array_iterator1 * vval ;
                *(dest+2) += *(array_iterator0+2) * vval ;
                *(dest+3) += *(array_iterator1+2) * vval ;
                vval =  v[idx[j]*stride+1] ;
                *(dest+0) += *(array_iterator0+4) * vval ;
                *(dest+1) += *(array_iterator1+4) * vval ;
                *(dest+2) += *(array_iterator0+6) * vval ;
                *(dest+3) += *(array_iterator1+6) * vval ;
                vval =  v[idx[j]*stride+2] ;
                *(dest+0) += *(array_iterator0+8) * vval ;
                *(dest+1) += *(array_iterator1+8) * vval ;
                *(dest+2) += *(array_iterator0+10) * vval ;
                *(dest+3) += *(array_iterator1+10) * vval ;
            }
#endif
        }

        if(stride == 6)
        {
            const int colLength = stride + stride%2 ;
#ifdef HAVE_SSE3
            const __m128d * array_iterator = (__m128d*)&val[start*colLength*stride] ;
            for(unsigned int j = start ; j != length+start ; j++)
            {
                __m128d vval =  _mm_set1_pd(v[idx[j]*stride]) ;

                _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*array_iterator, vval))) ;
                _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+1), vval))) ;
                _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+2), vval))) ;

                vval =  _mm_set1_pd(v[idx[j]*stride+1]) ;

                _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+3,) vval))) ;
                _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+4), vval))) ;
                _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+5), vval))) ;

                vval =  _mm_set1_pd(v[idx[j]*stride+2]) ;

                _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+6), vval))) ;
                _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+7), vval))) ;
                _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+8), vval))) ;

                vval =  _mm_set1_pd(v[idx[j]*stride+3]) ;

                _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+9), vval))) ;
                _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+10), vval))) ;
                _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+11), vval))) ;

                vval =  _mm_set1_pd(v[idx[j]*stride+4]) ;

                _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+12), vval))) ;
                _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+13), vval))) ;
                _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+14), vval))) ;

                vval =  _mm_set1_pd(v[idx[j]*stride+5]) ;

                _mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+15), vval))) ;
                _mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+16), vval))) ;
                _mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+17), vval))) ;
                array_iterator +=18 ;
            }
#else
            const double * array_iterator0 = &val[start*colLength*stride] ;
            const double * array_iterator1 = &val[start*colLength*stride+1] ;
            for(unsigned int j = start ; j != length+start ; j++)
            {
                double vval =  v[idx[j]*stride] ;

                *(dest) += *array_iterator0 * vval ;
                *(dest+1) += *array_iterator1 * vval ;

                *(dest+2) += *(array_iterator0+2) * vval ;
                *(dest+3) += *(array_iterator1+2) * vval ;

                *(dest+4) += *(array_iterator0+4) * vval ;
                *(dest+5) += *(array_iterator1+4) * vval ;

                vval =  v[idx[j]*stride+1] ;

                *(dest) += *(array_iterator0+6) * vval ;
                *(dest+1) += *(array_iterator1+6) * vval ;

                *(dest+2) += *(array_iterator0+8) * vval ;
                *(dest+3) += *(array_iterator1+8) * vval ;

                *(dest+4) += *(array_iterator0+10) * vval ;
                *(dest+5) += *(array_iterator1+10) * vval ;

                vval =  v[idx[j]*stride+2] ;

                *(dest) += *(array_iterator0+12) * vval ;
                *(dest+1) += *(array_iterator1+12) * vval ;

                *(dest+2) += *(array_iterator0+14) * vval ;
                *(dest+3) += *(array_iterator1+14) * vval ;

                *(dest+4) += *(array_iterator0+16) * vval ;
                *(dest+5) += *(array_iterator1+16) * vval ;

                vval =  v[idx[j]*stride+3] ;

                *(dest) += *(array_iterator0+18) * vval ;
                *(dest+1) += *(array_iterator1+18) * vval ;

                *(dest+2) += *(array_iterator0+20) * vval ;
                *(dest+3) += *(array_iterator1+20) * vval ;

                *(dest+4) += *(array_iterator0+22) * vval ;
                *(dest+5) += *(array_iterator1+22) * vval ;

                vval =  v[idx[j]*stride+4] ;

                *(dest) += *(array_iterator0+24) * vval ;
                *(dest+1) += *(array_iterator1+24) * vval ;

                *(dest+2) += *(array_iterator0+26) * vval ;
                *(dest+3) += *(array_iterator1+26) * vval ;

                *(dest+4) += *(array_iterator0+28) * vval ;
                *(dest+5) += *(array_iterator1+28) * vval ;

                vval =  v[idx[j]*stride+5] ;

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
        }

        const int colLength = stride + stride%2 ;
#ifdef HAVE_SSE3
        const __m128d * array_iterator = (__m128d*)&val[start*colLength*stride] ;
        for(unsigned int j = start ; j != length+start ; j++)
        {
            for(size_t c = 0 ; c != stride ; ++c)
            {
                const __m128d vval =  _mm_set1_pd(v[idx[j]*stride+c]) ;
                for(int i = 0 ; i != colLength/2 ; ++i,++array_iterator)
                {
                    _mm_store_pd((dest+i*2),  _mm_add_pd( _mm_load_pd((dest+i*2)), _mm_mul_pd(*array_iterator, vval))) ;
                }
            }
        }
#else
        const double * array_iterator0 = &val[start*colLength*stride] ;
        const double * array_iterator1 = &val[start*colLength*stride+1] ;
        for(unsigned int j = start ; j != length+start ; j++)
        {
            for(size_t c = 0 ; c < stride ; c++)
            {
                const double vval =  v[idx[j]*stride+c] ;
                for(int i = 0 ; i != colLength ; i+=2)
                {
                    *(dest+i) += *array_iterator0 * vval ;
                    *(dest+i+1) += *array_iterator1 * vval ;
                    array_iterator0 += 2 ;
                    array_iterator1 += 2 ;
                }
            }
        }
#endif
    }


    /** \brief simultaneously compute a number of dot products equal to the block size and assign.
     */
    inline double inner_product(const SparseVector & v0, const SparseVector & v1, const size_t end)
    {

        double ret = 0 ;

        unsigned int i = 0 ;
        unsigned int j = 0 ;
        unsigned int *i_index_pointer = std::lower_bound(&v0.idx[v0.start], &v0.idx[v0.start+v0.length], end) ;
        unsigned int *j_index_pointer = std::lower_bound(&v1.idx[v1.start], &v1.idx[v1.start+v1.length], end) ;
        unsigned int i_end = i_index_pointer-&v0.idx[v0.start] ;
        unsigned int j_end = j_index_pointer-&v1.idx[v1.start] ;

        if(v0.idx[v0.start] > v1.idx[v1.start])
        {
            j_index_pointer = std::lower_bound(&v1.idx[v1.start+j], &v1.idx[v1.start+v1.length], v0.idx[v0.start+i]) ;
            j = j_index_pointer-&v1.idx[v1.start] ;
        }
        else if(v0.idx[v0.start] < v1.idx[v1.start])
        {
            i_index_pointer = std::lower_bound(&v0.idx[v0.start+i], &v0.idx[v0.start+v0.length], v1.idx[v1.start + j]) ;
            i = i_index_pointer-&v0.idx[v0.start] ;
        }

        while(i < i_end &&  j < j_end)
        {
            ret += v0.val[v0.start +i] * v1.val[v1.start+j] ;
            i++ ;
            j++ ;

            if(v0.idx[v0.start+i] > v1.idx[v1.start+j])
            {
                j_index_pointer = std::lower_bound(&v1.idx[v1.start+j], &v1.idx[v1.start+v1.length], v0.idx[v0.start+i]) ;
                j = j_index_pointer-&v1.idx[v1.start] ;
            }
            else if(v0.idx[v0.start+i] < v1.idx[v1.start+j])
            {
                i_index_pointer = std::lower_bound(&v0.idx[v0.start+i], &v0.idx[v0.start+v0.length],v1.idx[v1.start + j]) ;
                i = i_index_pointer-&v0.idx[v0.start] ;
            }

        }
        return ret ;
    }


    /** \brief simultaneously compute a number of dot products equal to the block size, squared.
     * This effectively computes a series of matrix-matrix multiplications. The two sparse vectors
     * need not have the same sparsity pattern. They must, however have the smae block size.
     */
    Matrix operator *(const SparseVector&) const ;

    /** \brief simultaneously compute a number vector-vector additions equal to the block size.
     */
    Vector operator +(const Vector&) const ;

} ;


/** \brief Sparse Vector implementation.
 * The sparse vector implementation is not intended to be used as a standalone class. Instead it behaves as a row of a
 * CoordinateIndexedSparseMatrix. As such a matrix is stored in blocks, the sparse vector holds references to a group
 * of rows, and an aditional member is used to indicate which should be used for access.
 *
 * This class does not allow modification of members.
 */

struct ConstSparseVector
{
public:
    const Vector & array ;
    const std::valarray<unsigned int> & column_index ;
    const size_t length ;
    const size_t start ;
    const size_t stride ;
    const size_t index ;

public:
    /** \brief SparseVector constructor. Initialises the references to the data.
    @param v array containing the values
    @param idx column index of the blocks
    @param l number of blocks in the row.
    @param s block start index of the values relevant to this row
    @param index row number (actual number, not in blocks)
    @param st stride: block size
     */
    ConstSparseVector(const Vector & v,  const std::valarray<unsigned int>& idx , const size_t l, const size_t s, const size_t index, const size_t st) ;

    /** \brief access a value in the row
     * Return the value in the ith column of the matrix on this row. If this value is not stored, 0 is retured.
     */
    inline double operator [](const size_t i) const
    {
        const unsigned int * __start__       = &column_index[start] ;
        const unsigned int * __end__         = &column_index[start+length] ;
        const unsigned int * i_index_pointer = std::lower_bound(__start__, __end__, i/stride) ;
        unsigned int offset            = i_index_pointer - __start__ ;
        unsigned int colLength = stride+stride%2 ;
        return(std::binary_search(__start__, __end__, i/stride)) ?
              array[(start+offset)*stride*colLength + (i%stride)*colLength+index%stride] : 0 ;
    }

    /** \brief simultaneously compute a number of dot products equal to the block size.
     */
    Vector operator *(const Vector&) const ;

    /** \brief simultaneously compute a number of dot products equal to the block size and assign.
    */
    void inner_product(const Vector &v, double *dest, const size_t rowstart = 0, size_t const colstart = 0) const ;

    void parallel_product(const Vector &v, double *dest, const size_t rowstart = 0, const size_t colstart = 0) const ;



    /** \brief simultaneously compute a number of dot products equal to the block size, squared.
     * This effectively computes a series of matrix-matrix multiplications. The two sparse vectors
     * need not have the same sparsity pattern. They must, however have the same block size.
     */
    Matrix operator *(const SparseVector&) const ;

    /** \brief simultaneously compute a number of dot products equal to the block size, squared.
     * This effectively computes a series of matrix-matrix multiplications. The two sparse vectors
     * need not have the same sparsity pattern. They must, however have the smae block size.
     */
    Matrix operator *(const ConstSparseVector&) const ;

    /** \brief simultaneously compute a number vector-vector additions equal to the block size.
     */
    Vector operator +(const Vector&) const ;

    /** \brief Print the values of these rows. Zeros are printed where no value is stored.
     */
    void print() const ;

} ;


struct SparseMaskVector
{
public:
    std::valarray<bool> & array ;
    const std::valarray<unsigned int> & column_index ;
    const size_t length ;
    const size_t start ;
    const size_t stride ;
    const size_t index ;
    bool f ;

public:
    /** \brief SparseVector constructor. Initialises the references to the data.
    @param v array containing the values
    @param idx column index of the blocks
    @param l number of blocks in the row.
    @param s block start index of the values relevant to this row
    @param index row number (actual number, not in blocks)
    @param st stride: block size
     */
    SparseMaskVector(std::valarray<bool> & v,  const std::valarray<unsigned int>& idx , const size_t l, const size_t s, const size_t index, const size_t st) ;

    /** \brief access a value in the row
     * Return the value in the ith column of the matrix on this row. If this value is not stored, 0 is retured.
     */
    inline bool & operator [](const size_t i)
    {
        f = false ;
        const unsigned int * __start__       = &column_index[start] ;
        const unsigned int * __end__         = &column_index[start+length] ;
        const unsigned int * i_index_pointer = std::lower_bound(__start__, __end__, i/stride) ;
        unsigned int offset            = i_index_pointer - __start__ ;
        unsigned int colLength = stride+stride%2 ;
        return(std::binary_search(__start__, __end__, i/stride)) ?
              array[(start+offset)*stride*colLength + (i%stride)*colLength+index%stride] : f ;
    }


} ;



inline void reverseInnerProductAssignAndAdd(const Amie::ConstSparseVector & v0, Vector & v1, double &t,  double toAdd, size_t start)
{
    int stride = v0.stride ;
    int delta = stride+stride%2 ;
    for(size_t j = v0.length+v0.start-1 ; v0.column_index[j] > start   ; --j)
    {
        for(int i = 0 ; i < stride ; i++)
            t += v1[v0.column_index[j]*stride+i]*v0.array[j*stride*delta+i*delta] ;
    }
    t+=toAdd ;
}

inline void reverseInnerProductAssignAndAdd(const Amie::SparseVector & v0, Vector & v1, double &t,  double toAdd, size_t start)
{
    int stride = v0.stride ;
    int delta = stride+stride%2 ;
    for(size_t j = v0.length+v0.start-1 ; v0.idx[j] > start   ; --j)
    {
        for(int i = 0 ; i < stride ; i++)
            t += v1[v0.idx[j]*stride+i]*v0.val[j*stride*delta+i*delta] ;
    }
    t+=toAdd ;
}

inline void innerProductAssignAndAdd(const Amie::ConstSparseVector & v0, Vector & v1, double &t,  double toAdd, size_t end)
{
    int stride = v0.stride ;
    int delta = stride+stride%2 ;
    for(size_t j =  v0.start; v0.column_index[j] < end ; ++j)
    {
        for(int i = 0 ; i < stride ; i++)
            t += v1[v0.column_index[j]*stride+i]*v0.array[j*stride*delta+i*delta] ;
    }
    t+=toAdd ;
}

inline void innerProductAssignAndAdd(const Amie::SparseVector & v0, Vector & v1, double &t,  double toAdd, size_t end)
{
    int stride = v0.stride ;
    int delta = stride+stride%2 ;
    for(size_t j =  v0.start; v0.idx[j] < end ; ++j)
    {
        for(int i = 0 ; i < stride ; i++)
            t += v1[v0.idx[j]*stride+i]*v0.val[j*stride*delta+i*delta] ;
    }
    t+=toAdd ;
}

inline double innerProduct(const Amie::ConstSparseVector & v0, Amie::ConstSparseVector & v1, int s)
{
    Amie::Matrix ret(v0.stride, v0.stride) ;
    int colLength = v0.stride+v0.stride%2 ;
    int blocksize = v0.stride*colLength ;

    size_t i = 0 ;
    size_t j = 0 ;
    while(i < v0.length && j < v1.length)
    {
        if(v0.column_index[v0.start+i] > v1.column_index[v1.start+ j] || v0.column_index[v0.start+i] < s/v0.stride)
            j++ ;
        else if(v0.column_index[v0.start+i] < v1.column_index[v1.start+ j]|| v1.column_index[v1.start+j] < s/v0.stride)
            i++ ;
        else
        {
            int idx = 0 ;
            for(size_t k = 0 ; k < v0.stride ; k++)
            {
                for(size_t l = 0 ; l < v0.stride ; l++)
                {
                    ret[k][l] += v0.array[v0.start*blocksize+i*blocksize+colLength*k+l]*v1.array[v1.start*blocksize+j*blocksize + idx++] ;
                }
            }
            i++ ;
            j++ ;
        }
    }

    return ret[v0.index%v0.stride][v1.index%v0.stride] ;

}

}
#endif
