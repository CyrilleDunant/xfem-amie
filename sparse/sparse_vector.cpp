//
// C++ Implementation: sparse_vector
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
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
	if(std::binary_search(__start__, __end__, i/stride))
		return val[(start+offset)*stride*colLength + (i%stride)*colLength+index%stride] ;
	
	return 0 ;

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
	if(std::binary_search(__start__, __end__, i/stride))
		return val[(start+offset)*stride*colLength + (i%stride)*colLength+index%stride] ;
	
	return zero ;

}


void SparseVector::inner_product(const Vector &v, double *dest) const
{
	if(stride == 2)
	{
	#ifdef HAVE_SSE3
		const __m128d * array_iterator = (__m128d*)&val[start*4] ;
		const double * vec_iterator = &v[idx[start]*2] ;
		for(unsigned int j = start ; j != length+start ; j++)
		{
			
			_mm_store_pd(dest,  _mm_add_pd( _mm_load_pd(dest), _mm_add_pd( _mm_mul_pd(*array_iterator,  _mm_set1_pd(*vec_iterator)), _mm_mul_pd(*(array_iterator+1), _mm_set1_pd(*(vec_iterator+1)))))) ;
			array_iterator+=2 ;
			if(j+1 == length+start)
				break ;
			vec_iterator += idx[j+1]*2-idx[j]*2 ;
			
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

void ConstSparseVector::inner_product(const Vector &v, double *dest) const
{
	if(stride == 2)
	{
	#ifdef HAVE_SSE3
		const __m128d * array_iterator = (__m128d*)&val[start*4] ;
		const double * vec_iterator = &v[idx[start]*2] ;
		for(unsigned int j = start ; j != length+start ; j++)
		{
			
			_mm_store_pd(dest,  _mm_add_pd( _mm_load_pd(dest), _mm_add_pd( _mm_mul_pd(*array_iterator,  _mm_set1_pd(*vec_iterator)), _mm_mul_pd(*(array_iterator+1), _mm_set1_pd(*(vec_iterator+1)))))) ;
			array_iterator+=2 ;
			if(j+1 == length+start)
				break ;
			vec_iterator += idx[j+1]*2-idx[j]*2 ;
			
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

double innerProduct(const SparseVector & v0, const SparseVector & v1, const size_t end)
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


Vector ConstSparseVector::operator *(const Vector &v) const
{
	if(stride == 2)
	{
	#ifdef HAVE_SSE3
		Vector ret(0., 2) ;
		const __m128d * array_iterator = (__m128d*)&val[start*4] ;
		const double * vec_iterator = &v[idx[start]*2] ;
		for(unsigned int j = start ; j != length+start ; j++)
		{
			
			_mm_store_pd(&ret[0],  _mm_add_pd( _mm_load_pd(&ret[0]), _mm_add_pd( _mm_mul_pd(*array_iterator,  _mm_set1_pd(*vec_iterator)), _mm_mul_pd(*(array_iterator+1), _mm_set1_pd(*(vec_iterator+1)))))) ;
			array_iterator+=2 ;
			if(j+1 == length+start)
				break ;
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



