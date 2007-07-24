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

using namespace Mu ;

SparseVector::SparseVector(Vector & v, std::valarray<unsigned int> & i , const size_t l , const size_t s) : val(v), idx(i), length(l), start(s)
{
	zero = 0 ;
}

double SparseVector::operator [](size_t i) const
{
	unsigned int * __start__       = &idx[start] ;
	unsigned int * __end__         = &idx[start+length] ;
	unsigned int * i_index_pointer = std::find(__start__, __end__, i) ;
	unsigned int offset            = i_index_pointer - __start__ ;
	
	if(i_index_pointer !=  __end__)
		return val[start+offset] ;
	
	return 0 ;
// 	for(size_t j = 0 ; j < length ; j++)
// 	{
// 		if((*(idx + j)) == i)
// 			return (*(val + j)) ;
// 		if((*(idx + j)) > i)
// 			return 0 ;
// 	}
// 	return 0 ;
}

double & SparseVector::operator [](const size_t i) 
{
	zero = 0 ;
	unsigned int * __start__       = &idx[start] ;
	unsigned int * __end__         = &idx[start+length] ;
	unsigned int * i_index_pointer = std::find(__start__, __end__, i) ;
	unsigned int offset            = i_index_pointer - __start__ ;
	
	if(i_index_pointer !=  __end__)
		return val[start+offset] ;
	
	return zero ;
// 	for(size_t j = 0 ; j < length ; j++)
// 	{
// 		if((*(idx + j)) == i)
// 			return (*(val + j)) ;
// 		if((*(idx + j)) > i)
// 			return zero ;
// 	}
// 	return zero ;
}


double SparseVector::operator *(const Vector &v) const
{
	double ret = 0 ;
	for(size_t j = start ; j < length+start ; j++)
	{
		size_t index = idx[j] ;
		ret += v[index]*val[j] ; 
	}
	return ret ;
}

double SparseVector::operator *(const SparseVector &v) const
{
	double ret = 0 ;
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
			ret += val[start + i] * v.val[v.start+ j] ;
			i++ ;
			j++ ;
		}
	}
	return ret ;
	
// 	for(size_t j = 0 ; j < length ; j++)
// 	{
// 		ret += v[(*(idx + j))]*(*(val + j)) ; 
// 	}
// 	return ret ;
}

double innerProduct(const SparseVector & v0, const SparseVector & v1, const size_t end)
{
	double ret = 0 ;

	unsigned int i = 0 ; 
	unsigned int j = 0 ; 
	unsigned int *i_index_pointer = std::find(&v0.idx[v0.start], &v0.idx[v0.start+v0.length], end) ;
	unsigned int *j_index_pointer = std::find(&v1.idx[v1.start], &v1.idx[v1.start+v1.length], end) ;
	unsigned int i_end = i_index_pointer-&v0.idx[v0.start] ;
	unsigned int j_end = j_index_pointer-&v1.idx[v1.start] ;
	
	if(v0.idx[v0.start] > v1.idx[v1.start])
	{
		j_index_pointer = std::find(&v1.idx[v1.start+j], &v1.idx[v1.start+v1.length], v0.idx[v0.start+i]) ;
		j = j_index_pointer-&v1.idx[v1.start] ;
	}
	else if(v0.idx[v0.start] < v1.idx[v1.start])
	{
		i_index_pointer = std::find(&v0.idx[v0.start+i], &v0.idx[v0.start+v0.length], v1.idx[v1.start + j]) ;
		i = i_index_pointer-&v0.idx[v0.start] ;
	}
	
	while(i < i_end &&  j < j_end)
	{
		ret += v0.val[v0.start +i] * v1.val[v1.start+j] ;
		i++ ;
		j++ ;
		
		if(v0.idx[v0.start+i] > v1.idx[v1.start+j])
		{
			j_index_pointer = std::find(&v1.idx[v1.start+j], &v1.idx[v1.start+v1.length], v0.idx[v0.start+i]) ;
			j = j_index_pointer-&v1.idx[v1.start] ;
		}
		else if(v0.idx[v0.start+i] < v1.idx[v1.start+j])
		{
			i_index_pointer = std::find(&v0.idx[v0.start+i], &v0.idx[v0.start+v0.length],v1.idx[v1.start + j]) ;
			i = i_index_pointer-&v0.idx[v0.start] ;
		}

	}
	return ret ;
}

double innerProduct(const SparseVector & v0, const Vector & v1, const size_t end) 
{
	double ret = 0 ;
	for(size_t j = v0.start ; v0.idx[j] < end /*&& j < v0.length*/ ; ++j)
	{
		ret += v1[v0.idx[j]]*v0.val[j] ; 
	}
	return ret ;
}

double innerProduct(const ConstSparseVector & v0, const Vector & v1, const size_t end) 
{
	double ret = 0 ;
	for(size_t j =  v0.start; v0.idx[j] < end /*&& j < v0.length*/ ; ++j)
	{
		ret += v1[v0.idx[j]]*v0.val[j] ; 
// 		ret += v1[v0.idx][j]*v0.val[j] ; 
	}
	return ret ;
}

// inline void innerProductAssignAndAdd(const Mu::ConstSparseVector & v0, Vector & v1, double &t,  double toAdd, size_t end)
// {
// 	for(size_t j =  v0.start; v0.idx[j] < end ; ++j)
// 	{
// 		t += v1[v0.idx[j]]*v0.val[j] ; 
// 	}
// 	t+=toAdd ;
// }

double reverseInnerProduct(const SparseVector & v0, const Vector & v1, const size_t s) 
{
	double ret = 0 ;
	for(size_t j = v0.length+v0.start-1 ; /*j >0 && */v0.idx[j] > s   ; --j)
	{
		ret += v1[v0.idx[j]]*v0.val[j] ; 
	}
	return ret ;
}

// inline void reverseInnerProductAssignAndAdd(const Mu::ConstSparseVector & v0, Vector & v1, double &t,  double toAdd, size_t start)
// {
// 	for(size_t j = v0.length+v0.start-1 ; v0.idx[j] > start   ; --j)
// 	{
// 		t += v1[v0.idx[j]]*v0.val[j] ; 
// 	}
// 	t+=toAdd ;
// }

double reverseInnerProduct(const ConstSparseVector & v0, const Vector & v1, const size_t s) 
{
	double ret = 0 ;
	
	for(size_t j = v0.length+v0.start-1 ; v0.idx[j] > s ; --j)
	{
		ret += v1[v0.idx[j]]*v0.val[j] ; 
// 		ret += v1[v0.idx][j]*v0.val[j] ; 
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


ConstSparseVector::ConstSparseVector(const Vector & v, const std::valarray<unsigned int> & i , const size_t l , const size_t s) : val(v), idx(i), length(l), start(s)
{
}

// double ConstSparseVector::operator [](const size_t i) const
// {
// 	const size_t *i_index_pointer = std::find(&idx[start], &idx[std::min(start+length,idx.size())], i) ;
// 	size_t offset = i_index_pointer - &idx[start] ;
// 	if(i_index_pointer != &idx[std::min(start+length,idx.size())])
// 		return (val[start+offset]) ;
// 	
// 	return 0 ;
// // 	for(size_t j = 0 ; j < length ; j++)
// // 	{
// // 		if((*(idx + j)) == i)
// // 			return (*(val + j)) ;
// // 		if((*(idx + j)) > i)
// // 			return 0 ;
// // 	}
// // 	return 0 ;
// }

double ConstSparseVector::operator *(const Vector &v) const
{
	double ret = 0 ;
	
	for(size_t j = start ; j < length+start ; j++)
	{
		size_t index = idx[j] ;
		ret += v[index]*val[j] ; 
	}
	return ret ;
// 	double ret = 0 ;
// 	for(size_t j = 0 ; j < length ; j++)
// 	{
// 		ret += v[(*(idx + j))]*(*(val + j)) ; 
// 	}
// 	return ret ;
}

double ConstSparseVector::operator *(const SparseVector& v) const
{
	double ret = 0 ;
	size_t i = start ; 
	size_t j = v.start ; 
	while(i < length+start && j < v.length+v.start)
	{

		while(idx[i] > v.idx[j])
			j++ ;

		while(idx[i] < v.idx[j])
			i++ ;
		
		ret += val[i] * v.val[j] ;
		i++ ;
		j++ ;
	}
	return ret ;
	
}

double ConstSparseVector::operator *(const ConstSparseVector& v) const
{
	double ret = 0 ;
	size_t i = start ; 
	size_t j = v.start ; 
	while(i < length+start && j < v.length+v.start)
	{
		
		while(idx[i] > v.idx[j])
			j++ ;
		
		while(idx[i] < v.idx[j])
			i++ ;
		
		ret += val[i] * v.val[j] ;
		i++ ;
		j++ ;
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
// 	Vector ret(v) ;
// 	
// 	for(size_t j = 0 ; j < length ; j++)
// 	{
// 		ret[(*(idx + j))] += (*(val + j)) ; 
// 	}
// 	
// 	return ret ;
}
	

// struct BandSparseVector
// {
// public:
// 	Vector & val ;
// 	const size_t length ;
// 	const size_t start ;
// 	
// 	double zero ;
// 	
// public:
// 	BandSparseVector(Vector & v , const size_t l , const size_t s) ;
// 	
// 	double operator [](const size_t) const ;
// 	double & operator [](const size_t) ;
// 	double operator *(const Vector&) const ;
// 	double operator *(const SparseVector&) const ;
// 	Vector operator +(const Vector&) const ;
// 	
// } ;


ConstBandSparseVector::ConstBandSparseVector(const Vector & v, const size_t l, const size_t s) : val(v), length(l), start(s) { };
	
double ConstBandSparseVector::operator [](const size_t i) const
{
	if(i < start || i >=start+length)
		return 0 ;
	else
		return val[i-start] ;
}
double ConstBandSparseVector::operator *(const Vector& v) const
{
	double ret = 0 ;
	
	for(size_t i = 0 ;  i < length  ; i++)
		ret+= v[i+start]*val[i] ;
	
	return ret ;
}
double ConstBandSparseVector::operator *(const SparseVector&v) const
{
	double ret = 0 ;
	
	size_t index = 0 ;
	size_t sparse_index = v.idx[v.start]-start ;
	
	if(sparse_index >= start)
	{
		while( (sparse_index < start+length))
		{
			ret+= val[sparse_index]*v.val[index] ;
			index++ ;
			sparse_index = v.idx[v.start + index]-start ;
		}
	}
	

	return ret ;
}

double ConstBandSparseVector::operator *(const ConstSparseVector& v) const
{
	double ret = 0 ;
	
	size_t index = 0 ;
	size_t sparse_index = v.idx[v.start] ;
	
	if(sparse_index >= start)
	{
		while( (sparse_index < start+length))
		{
			ret+= val[sparse_index-start]*v.val[index] ;
			index++ ;
			sparse_index = v.idx[v.start + index] ;
		}
	}
	
	
	return ret ;
}


double ConstBandSparseVector::operator *(const ConstBandSparseVector& v) const
{
	double ret = 0 ;
	
	if(start < v.start)
	{
		if(v.start-start < length)
			for (size_t i = v.start-start ; i< length ; i++)
				ret += val[i]*v.val[i-v.start+start] ;
	}
	else
	{
		if(start-v.start < v.length)
			for (size_t i = start-v.start ; i< v.length ; i++)
				ret += val[i-start+v.start]*v.val[i] ;
	}
	return ret ;
}

double ConstBandSparseVector::operator *(const BandSparseVector& v) const
{
	double ret = 0 ;
	
	if(start < v.start)
	{
		size_t j = 0 ;
		for (size_t i = v.start-start ; i< length ; i++ )
		{
			ret += val[i]*v.val[j] ;
			j++ ;
		}
		
	}
	else
	{
		size_t j = 0 ;
		for (size_t i = start-v.start ; i< v.length ; i++)
		{
			ret += val[j]*v.val[i] ;
			j++ ;
		}
	}
	return ret ;
}

Vector ConstBandSparseVector::operator +(const Vector& v) const
{
	Vector ret(v) ;
	
	for(size_t i = start ; i < start+length ; i++)
	{
		ret[i]+=v[i] ;
	}
	
	return ret ;
}

	


