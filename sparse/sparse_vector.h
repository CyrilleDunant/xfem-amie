//
// C++ Interface: sparse_vector
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//


#ifndef __SPARSE_VECTOR_H
#define __SPARSE_VECTOR_H

#include <valarray>
#include <iostream>
#include "../matrixops.h"
namespace Mu
{

struct SparseVector
{
public:
	Vector & val ;
	std::valarray<unsigned int> & idx ;
	const size_t length ;
	const size_t start ;
	
	double zero ;
	
public:
	SparseVector(Vector & v, std::valarray<unsigned int> & idx , const size_t l , const size_t s) ;
	
	double operator [](const size_t) const ;
	double & operator [](const size_t) ;
	double operator *(const Vector&) const ;
	double operator *(const SparseVector&) const ;
	Vector operator +(const Vector&) const ;

} ;

struct ConstSparseVector
{
public:
	const Vector & val ;
	const std::valarray<unsigned int> & idx ;
	const size_t length ;
	const size_t start ;
	
public:
	ConstSparseVector(const Vector & v,  const std::valarray<unsigned int>& idx , const size_t l, const size_t s) ;
	
	inline double operator [](const size_t i) const
	{
		const unsigned int *i_index_pointer = std::find(&idx[start], &idx[std::min(start+length,idx.size())], i) ;
		unsigned int offset = i_index_pointer - &idx[start] ;
		if(i_index_pointer != &idx[std::min(start+length,idx.size())])
			return (val[start+offset]) ;
		
		return 0 ;
	}
	double operator *(const Vector&) const ;
	double operator *(const SparseVector&) const ;
	double operator *(const ConstSparseVector&) const ;
	Vector operator +(const Vector&) const ;
	
	void print() const ;
	
} ;

struct SymetricSparseVector
{
public:
	Vector & val ;
	std::valarray<size_t> & idx ;
	const size_t length ;
	const size_t start ;
	
	double zero ;
	
public:
	SymetricSparseVector(Vector & v, std::valarray<size_t> & idx , const size_t l , const size_t s) : val(v), idx(idx), length(l), start(s)
{
	zero = 0 ;
}
	
	double operator [](const size_t) const {return 0 ;};
	double & operator [](const size_t) {return zero ;};
	double operator *(const Vector&) const {return 0 ;};
	double operator *(const SymetricSparseVector&) const {return 0 ;};
	Vector operator +(const Vector&) const {return Vector(0) ;};

} ;

struct ConstBandSparseVector ;

struct BandSparseVector
{
public:
	Vector & val ;
	const size_t length ;
	const size_t start ;
	
	double zero ;
	
public:
	BandSparseVector(Vector & v , const size_t l , const size_t s) ;
	
	double operator [](const size_t) const ;
	double & operator [](const size_t) ;
	double operator *(const Vector&) const ;
	double operator *(const SparseVector&) const ;
	double operator *(const BandSparseVector&) const ;
	double operator *(const ConstBandSparseVector&) const ;
	Vector operator +(const Vector&) const ;
	
} ;

struct ConstBandSparseVector
{
public:
	const Vector & val ;
	const size_t length ;
	const size_t start ;
	
public:
	ConstBandSparseVector(const Vector & v, const size_t l, const size_t s) ;
	
	inline double operator [](const size_t i) const ;
	double operator *(const Vector&) const ;
	double operator *(const SparseVector&) const ;
	double operator *(const ConstSparseVector&) const ;
	double operator *(const BandSparseVector&) const ;
	double operator *(const ConstBandSparseVector&) const ;
	Vector operator +(const Vector&) const ;
	
	void print() const ;
	
} ;

} ;

double innerProduct( const Mu::SparseVector & v0, const Mu::SparseVector & v1, const size_t end) ;
double innerProduct( const Mu::SparseVector & v0, const Vector & v1, const size_t end) ;
double innerProduct( const Mu::ConstSparseVector & v0, const Vector & v1, const size_t end) ;
inline void innerProductAssignAndAdd(const Mu::ConstSparseVector & v0, Vector & v1, double &t,  double toAdd, size_t end)
{
	for(size_t j =  v0.start; v0.idx[j] < end ; ++j)
	{
		t += v1[v0.idx[j]]*v0.val[j] ; 
	}
	t+=toAdd ;
} ;

double reverseInnerProduct(const Mu::SparseVector & v0, const Vector & v1, const size_t start) ;
double reverseInnerProduct(const Mu::ConstSparseVector & v0, const Vector & v1, const size_t start) ;
inline void reverseInnerProductAssignAndAdd(const Mu::ConstSparseVector & v0, Vector & v1, double &t,  double toAdd, size_t start)
{
	for(size_t j = v0.length+v0.start-1 ; v0.idx[j] > start   ; --j)
	{
		t += v1[v0.idx[j]]*v0.val[j] ; 
	}
	t+=toAdd ;
} ;

#endif
