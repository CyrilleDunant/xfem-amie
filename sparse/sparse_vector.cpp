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

using namespace Mu ;


SparseVector::SparseVector(Vector & v, std::valarray<unsigned int> & i , const size_t l , const size_t s, const size_t ind, const size_t st) : val(v), idx(i), length(l), start(s), stride(st), index(ind)
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
			for(size_t i = 0 ; i != colLength ; i+=2)
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


ConstSparseVector::ConstSparseVector(const Vector & v, const std::valarray<unsigned int> & i , const size_t l , const size_t s, const size_t ind, const size_t st) : val(v), idx(i), length(l), start(s), stride(st), index(ind)
{
}

void ConstSparseVector::inner_product(const Vector &v, double *dest, const size_t rowstart , size_t const colstart ) const
{
	size_t mstart = start ;
	if(colstart)
	{
		const unsigned int * __start__       = &idx[start] ;
		const unsigned int * __end__         = &idx[start+length] ;
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
				*dest += v[idx[j]]*val[j*2] ;
			}
			return ;
		}
		case 2:
		{
			const double * array_iterator = &val[mstart*2*2] ;
			const double * vec_iterator = &v[idx[mstart]*2] ;
			for(unsigned int j = mstart ; j < length+start ; j++ )
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
			const double * array_iterator = &val[mstart*3*4] ;
			const double * vec_iterator = &v[idx[mstart]*3] ;
			for(unsigned int j = mstart ; j < length+start ;j++)
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
		case 6 :
		{
			const int colLength = 6  ;
			#ifdef HAVE_SSE3
			const __m128d * array_iterator = (__m128d*)&val[mstart*36] ;
// 			#pragma omp parallel for schedule(runtime)
			for(unsigned int j = mstart ; j < length+start ; j++)
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
			const double * array_iterator0 = &val[mstart*36] ;
			const double * array_iterator1 = &val[mstart*36+1] ;
// 			#pragma omp parallel for schedule(runtime)
			for(unsigned int j = mstart ; j < length+start ; j++)
			{
				double vval =  v[idx[j]*6] ;

				*(dest) += *array_iterator0 * vval ;
				*(dest+1) += *array_iterator1 * vval ;
				
				*(dest+2) += *(array_iterator0+2) * vval ;
				*(dest+3) += *(array_iterator1+2) * vval ;
				
				*(dest+4) += *(array_iterator0+4) * vval ;
				*(dest+5) += *(array_iterator1+4) * vval ;
				
				vval =  v[idx[j]*6+1] ;

				*(dest) += *(array_iterator0+6) * vval ;
				*(dest+1) += *(array_iterator1+6) * vval ;
				
				*(dest+2) += *(array_iterator0+8) * vval ;
				*(dest+3) += *(array_iterator1+8) * vval ;
				
				*(dest+4) += *(array_iterator0+10) * vval ;
				*(dest+5) += *(array_iterator1+10) * vval ;
				
				vval =  v[idx[j]*6+2] ;

				*(dest) += *(array_iterator0+12) * vval ;
				*(dest+1) += *(array_iterator1+12) * vval ;
				
				*(dest+2) += *(array_iterator0+14) * vval ;
				*(dest+3) += *(array_iterator1+14) * vval ;
				
				*(dest+4) += *(array_iterator0+16) * vval ;
				*(dest+5) += *(array_iterator1+16) * vval ;
				
				vval =  v[idx[j]*6+3] ;

				*(dest) += *(array_iterator0+18) * vval ;
				*(dest+1) += *(array_iterator1+18) * vval ;
				
				*(dest+2) += *(array_iterator0+20) * vval ;
				*(dest+3) += *(array_iterator1+20) * vval ;
				
				*(dest+4) += *(array_iterator0+22) * vval ;
				*(dest+5) += *(array_iterator1+22) * vval ;
				
				vval =  v[idx[j]*6+4] ;

				*(dest) += *(array_iterator0+24) * vval ;
				*(dest+1) += *(array_iterator1+24) * vval ;
				
				*(dest+2) += *(array_iterator0+26) * vval ;
				*(dest+3) += *(array_iterator1+26) * vval ;
				
				*(dest+4) += *(array_iterator0+28) * vval ;
				*(dest+5) += *(array_iterator1+28) * vval ;
				
				vval =  v[idx[j]*6+5] ;

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
			const double * array_iterator0 = &val[mstart*colLength*stride] ;
			const double * array_iterator1 = &val[mstart*colLength*stride+1] ;
			for(unsigned int j = mstart ; j < length+start ; j++)
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



Vector ConstSparseVector::operator *(const Vector &v) const
{
	if(stride == 2)
	{
	#ifdef HAVE_SSE3
		Vector ret(0., 2) ;
		const __m128d * array_iterator = (__m128d*)&val[start*4] ;
		const double * vec_iterator = &v[idx[start]*2] ;
		for(unsigned int j = start ; j != length+start && j+1 != length+start; j++)
		{
			
			_mm_store_pd(&ret[0],  _mm_add_pd( _mm_load_pd(&ret[0]), _mm_add_pd( _mm_mul_pd(*array_iterator,  _mm_set1_pd(*vec_iterator)), _mm_mul_pd(*(array_iterator+1), _mm_set1_pd(*(vec_iterator+1)))))) ;
			array_iterator+=2 ;
			vec_iterator += idx[j+1]*2-idx[j]*2 ;
			
		}
	#else
		Vector ret(0., 2) ;
		const double * array_iterator0 = &val[start*4] ;
		const double * array_iterator1 = &val[start*4+1] ;
		const double * array_iterator2 = &val[start*4+2] ;
		const double * array_iterator3 = &val[start*4+3] ;
		for(unsigned int j = start ; j != length+start ; j++)
		{
			ret[0] += *array_iterator0*v[idx[j]*2]+*array_iterator2*v[idx[j]*2+1] ;
			ret[1] += *array_iterator1*v[idx[j]*2]+*array_iterator3*v[idx[j]*2+1] ;
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
	const __m128d * array_iterator = (__m128d*)&val[start*colLength*stride] ;
	for(unsigned int j = start ; j != length+start ; j++)
	{
		for(size_t c = 0 ; c < stride ; c++)
		{
			const __m128d vval =  _mm_set1_pd(v[idx[j]*stride+c]) ;
			for(int i = 0 ; i != colLength ; i+=2)
			{
				_mm_store_pd(&ret[i],  _mm_add_pd( _mm_load_pd(&ret[i]), _mm_mul_pd(*array_iterator, vval))) ;
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

Matrix ConstSparseVector::operator *(const ConstSparseVector& v) const
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

void ConstSparseVector::print() const
{
	for(size_t j = 0 ; j < length ; j++)
		std::cout << idx[start + j] << " -> " << val[start + j] << std::endl ;
}

Vector ConstSparseVector::operator +(const Vector& v) const
{
	Vector ret(v) ;
	
	for(size_t j = start ; j < length+start ; j++)
	{
		size_t index = idx[j] ;
		ret[index] += val[j] ; 
	}
	
	return ret ;

}



