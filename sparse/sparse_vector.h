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
namespace Mu
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
	
	/** \brief access a value in the row
	 * Return the value in the ith column of the matrix on this row. If this value is not stored, assignement has no effect.
	 */
	double & operator [](const size_t) ;
	
	/** \brief simultaneously compute a number of dot products equal to the block size.
	 */
	Vector operator *(const Vector&) const ;
	
	void inner_product(const Vector &v, double *dest) const
	{
		if(stride == 2)
		{
		#ifdef HAVE_SSE3
			const __m128d * array_iterator = (__m128d*)&val[start*4] ;
			const double * vec_iterator = &v[idx[start]*2] ;
			for(unsigned int j = start ; j != length+start && j+1 != length+start; vec_iterator += (idx[j+1]-idx[j])*2,j++)
			{
				
				_mm_store_pd(dest,  _mm_add_pd( _mm_load_pd(dest), _mm_add_pd( _mm_mul_pd(*array_iterator,  _mm_set1_pd(*vec_iterator)), _mm_mul_pd(*(array_iterator+1), _mm_set1_pd(*(vec_iterator+1)))))) ;
				array_iterator+=2 ;
				
			}
		#else
			const double * array_iterator0 = &val[start*4] ;
			const double * array_iterator1 = &val[start*4+1] ;
			const double * array_iterator2 = &val[start*4+2] ;
			const double * array_iterator3 = &val[start*4+3] ;
			for(unsigned int j = start ; j != length+start ; j++)
			{
				*dest += *array_iterator0*v[idx[j]*2]+*array_iterator2*v[idx[j]*2+1] ;
				*(dest+1) += *array_iterator1*v[idx[j]*2]+*array_iterator3*v[idx[j]*2+1] ;
				array_iterator0+=4 ;
				array_iterator1+=4 ;
				array_iterator2+=4 ;
				array_iterator3+=4 ;
			}
		#endif
		return ;
		}
		
		const int colLength = stride + stride%2 ;
		#ifdef HAVE_SSE3
		const __m128d * array_iterator = (__m128d*)&val[start*colLength*stride] ;
		for(unsigned int j = start ; j != length+start ; j++)
		{
			for(size_t c = 0 ; c < stride ; c++)
			{
				const __m128d vval =  _mm_set1_pd(v[idx[j]*stride+c]) ;
				for(int i = 0 ; i != colLength ; i+=2)
				{
					_mm_store_pd((dest+i),  _mm_add_pd( _mm_load_pd((dest+i)), _mm_mul_pd(*array_iterator, vval))) ;
					array_iterator++ ;
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
				for(size_t i = 0 ; i != colLength ; i+=2)
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
	const Vector & val ;
	const std::valarray<unsigned int> & idx ;
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
		const unsigned int * __start__       = &idx[start] ;
		const unsigned int * __end__         = &idx[start+length] ;
		const unsigned int * i_index_pointer = std::lower_bound(__start__, __end__, i/stride) ;
		unsigned int offset            = i_index_pointer - __start__ ;
		unsigned int colLength = stride+stride%2 ;
		return(std::binary_search(__start__, __end__, i/stride)) ?
			 val[(start+offset)*stride*colLength + (i%stride)*colLength+index%stride] : 0 ;
	}
	
	/** \brief simultaneously compute a number of dot products equal to the block size.
	 */
	Vector operator *(const Vector&) const ;
	
		/** \brief simultaneously compute a number of dot products equal to the block size and assign.
	 */
	void inner_product(const Vector &v, double *dest) const
	{

		switch(stride)
		{
			case 1:
			{
				for(unsigned int j = start ; j < length+start ; j++)
				{
					*dest += v[idx[j]]*val[j*2] ;
				}
				return ;
			}
			case 2:
			{
				const double * array_iterator = &val[start*2*2] ;
				const double * vec_iterator = &v[idx[start]*2] ;
				for(unsigned int j = start ; j < length+start ; j++ )
				{
					*dest     += *array_iterator*(*vec_iterator);
					*(dest+1) += *(array_iterator+1)*(*vec_iterator) ;
					*dest     += *(array_iterator+2)*(*(vec_iterator +1)) ;
					*(dest+1) += *(array_iterator+3)*(*(vec_iterator +1)) ;
					array_iterator+=4 ;
					if(j+1 < idx.size())
						vec_iterator += idx[j+1]*stride-idx[j]*stride ;
				}
				return ;
			}
			
			case 3:
			{
				const double * array_iterator = &val[start*3*4] ;
				const double * vec_iterator = &v[idx[start]*3] ;
				for(unsigned int j = start ; j < length+start ;j++)
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
					if(j+1 < idx.size())
					 vec_iterator += idx[j+1]*stride-idx[j]*stride ;
				}
				return ;
			}
			default:
			{
				int colLength = stride + stride%2 ;
				const double * array_iterator0 = &val[start*colLength*stride] ;
				const double * array_iterator1 = &val[start*colLength*stride+1] ;
				for(unsigned int j = start ; j != length+start ; j++)
				{
					for(size_t c = 0 ; c < stride ; c++)
					{
						const double vval =  v[idx[j]*stride+c] ;
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

} ;

inline void reverseInnerProductAssignAndAdd(const Mu::ConstSparseVector & v0, Vector & v1, double &t,  double toAdd, size_t start)
{
	int stride = v0.stride ;
	int delta = stride+stride%2 ;
	for(size_t j = v0.length+v0.start-1 ; v0.idx[j] > start   ; --j)
	{
		for(int i = 0 ; i < stride ; i++)
			t += v1[v0.idx[j]*stride+i]*v0.val[j*stride*delta+i*delta] ;
	}
	t+=toAdd ;
} ;


inline void innerProductAssignAndAdd(const Mu::ConstSparseVector & v0, Vector & v1, double &t,  double toAdd, size_t end)
{
	int stride = v0.stride ;
	int delta = stride+stride%2 ;
	for(size_t j =  v0.start; v0.idx[j] < end ; ++j)
	{
		for(int i = 0 ; i < stride ; i++)
			t += v1[v0.idx[j]*stride+i]*v0.val[j*stride*delta+i*delta] ;
	}
	t+=toAdd ;
} ;

inline double innerProduct(const Mu::ConstSparseVector & v0, Mu::ConstSparseVector & v1, int s)
{
	Mu::Matrix ret(v0.stride, v0.stride) ;
	int colLength = v0.stride+v0.stride%2 ;
	int blocksize = v0.stride*colLength ;
	
	size_t i = 0 ; 
	size_t j = 0 ; 
	while(i < v0.length && j < v1.length)
	{
		if(v0.idx[v0.start+i] > v1.idx[v1.start+ j] || v0.idx[v0.start+i] < s/v0.stride)
			j++ ;
		else if(v0.idx[v0.start+i] < v1.idx[v1.start+ j]|| v1.idx[v1.start+j] < s/v0.stride)
			i++ ;
		else
		{
			int idx = 0 ;
			for(size_t k = 0 ; k < v0.stride ; k++)
			{
				for(size_t l = 0 ; l < v0.stride ; l++)
				{
					ret[k][l] += v0.val[v0.start*blocksize+i*blocksize+colLength*k+l]*v1.val[v1.start*blocksize+j*blocksize + idx++] ;
				}
			}
			i++ ;
			j++ ;
		}
	}

	return ret[v0.index%v0.stride][v1.index%v0.stride] ;
	
} ;

#endif
