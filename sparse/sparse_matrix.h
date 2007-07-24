//
// C++ Interface: sparse_matrix
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __SPARSE_MATRIX_H
#define __SPARSE_MATRIX_H

#include<valarray>
#include<map>
#include <complex>
#include "../matrixops.h"
#include "sparse_vector.h"

namespace Mu {

struct ConstSymetricSparseVector ;

class CoordinateIndexedIncompleteSparseMatrix
{
public:
	Vector  array ;
	std::valarray<unsigned int>  column_index ;
	std::valarray<unsigned int>  row_index ;
public:
	
	CoordinateIndexedIncompleteSparseMatrix(const std::valarray<unsigned int> &, const std::valarray<unsigned int> &) ;
	CoordinateIndexedIncompleteSparseMatrix(const std::valarray<unsigned int> &, const Vector &, const size_t i) ;
	~CoordinateIndexedIncompleteSparseMatrix() ;
	
	SparseVector operator[](const size_t i) ;
	ConstSparseVector operator[](const size_t i) const ;
	double & operator()(const size_t i, const size_t j) ;
	
	Vector operator *(const Vector v) const ;
	
} ;

class CoordinateIndexedSparseMatrix
{
public:
	Vector  array ;
	std::valarray<unsigned int>  column_index ;
	std::valarray<unsigned int>  row_size ;
	std::valarray<unsigned int>  accumulated_row_size ;
public:
	CoordinateIndexedSparseMatrix(std::map<std::pair<size_t, size_t>, double> &source) ;
	CoordinateIndexedSparseMatrix(std::map<std::pair<size_t, size_t>, Matrix> &source) ;
	CoordinateIndexedSparseMatrix(const std::valarray<unsigned int> &, const std::valarray<unsigned int> &) ;
	~CoordinateIndexedSparseMatrix() ;
	
	SparseVector operator[](const size_t i) ;
	const ConstSparseVector operator[](const size_t i) const;
	double & operator()(const size_t i, const size_t j) ;
	double  operator()(const size_t i, const size_t j) const ;
		
	Vector operator *(const Vector v) const ;
// 	CoordinateIndexedSparseMatrix operator +(const CoordinateIndexedIncompleteSparseMatrix &) const ;
	CoordinateIndexedSparseMatrix & operator +=(const CoordinateIndexedIncompleteSparseMatrix &) ;
// 	CoordinateIndexedSparseMatrix operator +(const CoordinateIndexedSparseMatrix &) const ;
	CoordinateIndexedSparseMatrix & operator=(const CoordinateIndexedSparseMatrix &) ;
	
	Vector inverseDiagonal() const ;
	Vector inverseDiagonalSquared() const ;
	Vector diagonal() const ;
	CoordinateIndexedSparseMatrix transpose() const ;
	
	double froebeniusNorm() const ;
	double infinityNorm() const ;
	
} ;

class BandSparseMatrix
{
public:
	Vector  array ;
	size_t  row_size ;
	std::valarray<size_t>  start_index ;
public:
	BandSparseMatrix(size_t num_rows, size_t row_size) ;
	~BandSparseMatrix() ;
	
	BandSparseVector operator[](const size_t i) ;
	const ConstBandSparseVector operator[](const size_t i) const;
	double & operator()(const size_t i, const size_t j) ;
	double  operator()(const size_t i, const size_t j) const ;
	
	Vector operator *(const Vector v) const ;
	CoordinateIndexedSparseMatrix operator +(const CoordinateIndexedIncompleteSparseMatrix &) const ;
	CoordinateIndexedSparseMatrix & operator +=(const CoordinateIndexedIncompleteSparseMatrix &) ;
	CoordinateIndexedSparseMatrix operator +(const CoordinateIndexedSparseMatrix &) const ;
	CoordinateIndexedSparseMatrix & operator=(const CoordinateIndexedSparseMatrix &) ;
	
	Vector inverseDiagonal() const ;
	Vector inverseDiagonalSquared() const ;
	Vector diagonal() const ;
	BandSparseMatrix transpose() const ;
	
	double froebeniusNorm() const ;
	double infinityNorm() const ;
	
} ;

class SymetricSparseMatrix
{
public:
	Vector  array ;
	std::valarray<size_t>  column_index ;
	std::valarray<size_t>  row_size ;
	std::valarray<size_t>  accumulated_row_size ;
public:
	SymetricSparseMatrix(std::map<std::pair<size_t, size_t>, double> &source) ;
	SymetricSparseMatrix(std::map<std::pair<size_t, size_t>, Matrix> &source) ;
	SymetricSparseMatrix(const std::valarray<size_t> &, const std::valarray<size_t> &) ;
	~SymetricSparseMatrix() ;
	
	SymetricSparseVector operator[](const size_t i) ;
	const ConstSymetricSparseVector operator[](const size_t i) const;
	double & operator()(const size_t i, const size_t j) ;
	double  operator()(const size_t i, const size_t j) const ;
	
	Vector operator *(const Vector v) const ;
	SymetricSparseMatrix operator +(const CoordinateIndexedIncompleteSparseMatrix &) const ;
	SymetricSparseMatrix & operator +=(const CoordinateIndexedIncompleteSparseMatrix &) ;
	SymetricSparseMatrix operator +(const CoordinateIndexedSparseMatrix &) const ;
	SymetricSparseMatrix & operator=(const CoordinateIndexedSparseMatrix &) ;
	SymetricSparseMatrix operator +(const SymetricSparseMatrix &) const ;
	
	Vector inverseDiagonal() const ;
	Vector inverseDiagonalSquared() const ;
	Vector diagonal() const ;
	const SymetricSparseMatrix & transpose() const ;
	
	double froebeniusNorm() const ;
	double infinityNorm() const ;
	
} ;


struct CompositeSparseMatrixTimesVec ;
struct CompositeSparseMatrixTimesVecPlusVec ;
struct CompositeSparseMatrixTimesVecMinusVec ;
struct CompositeSparseMatrixTimesVecMinusVecMinusVec ;

struct CompositeSparseMatrix
{
	const CoordinateIndexedSparseMatrix  & sm ;
	const CoordinateIndexedIncompleteSparseMatrix  & ism ;
	CompositeSparseMatrix(const CoordinateIndexedSparseMatrix  & sm_ , const CoordinateIndexedIncompleteSparseMatrix  & ism_)  : sm(sm_), ism(ism_) { };
	
	CompositeSparseMatrixTimesVec operator *(const Vector v) const ;

	operator CoordinateIndexedSparseMatrix() const ;
} ;

struct CompositeSparseMatrixTimesVec
{
	const CompositeSparseMatrix  & sm ;
	const Vector & ve ;
	CompositeSparseMatrixTimesVec(const CompositeSparseMatrix  &s,const Vector & v ) : sm(s), ve(v) {} ;
	
	CompositeSparseMatrixTimesVecPlusVec operator + (const Vector & v) const ;
	CompositeSparseMatrixTimesVecMinusVec operator - (const Vector & v) const ;

	operator const Vector() const ;
} ;

struct CompositeSparseMatrixTimesVecPlusVec
{
	const CompositeSparseMatrixTimesVec & co ;
	const Vector & ve ;
	CompositeSparseMatrixTimesVecPlusVec(const CompositeSparseMatrixTimesVec & c, const Vector & v) : co(c), ve(v) {} ;
	operator const Vector() const 
	{
		Vector ret(0., ve.size()) ;
		
		for (size_t i = 0 ; i < co.sm.sm.row_size.size(); i++)
		{
			ret[i] += co.sm.sm[i]*co.ve +co.sm.ism[i]*co.ve +ve[i];
		}
		return ret ;
	}
	
} ;

struct CompositeSparseMatrixTimesVecMinusVec
{
	const CompositeSparseMatrixTimesVec & co ;
	const Vector & ve ;
	CompositeSparseMatrixTimesVecMinusVec(const CompositeSparseMatrixTimesVec & c, const Vector & v) : co(c), ve(v) {} ;
	
	CompositeSparseMatrixTimesVecMinusVecMinusVec operator - (const Vector & v) const ;
	
	operator const Vector() const 
	{
		Vector ret(0., ve.size()) ;
		
		for (size_t i = 0 ; i < co.sm.sm.row_size.size(); i++)
		{
			ret[i] += co.sm.sm[i]*co.ve +co.sm.ism[i]*co.ve -ve[i];
		}
		return ret ;
	}
	
} ;

struct CompositeSparseMatrixTimesVecMinusVecMinusVec
{
	const CompositeSparseMatrixTimesVecMinusVec & co ;
	const Vector & ve ;
	CompositeSparseMatrixTimesVecMinusVecMinusVec(const CompositeSparseMatrixTimesVecMinusVec & c, const Vector & v) : co(c), ve(v) {} ;
	operator const Vector() const 
	{
		Vector ret(0., ve.size()) ;
		
		for (size_t i = 0 ; i < co.co.sm.sm.row_size.size(); i++)
		{
			ret[i] += co.co.sm.sm[i]*co.co.ve +co.co.sm.ism[i]*co.co.ve -co.ve[i] -ve[i];
		}
		return ret ;
	}
	
} ;

struct ConstSymetricSparseVector
{
public:
	const size_t id ;
	const Vector & val ;
	const std::valarray<size_t> & idx ;
	const size_t length ;
	const size_t start ;
	const SymetricSparseMatrix & m ;
public:
	ConstSymetricSparseVector(const size_t i, const Vector & v,  const std::valarray<size_t>& idx , const size_t l, const size_t s, const SymetricSparseMatrix & mm) : id(i), val(v), idx(idx), length(l), start(s), m(mm)
{
}
	
	inline double operator [](const size_t i) const
	{
		if(i >= id)
		{
			const size_t *i_index_pointer = std::find(&idx[start], &idx[std::min(start+length,idx.size())], i) ;
			size_t offset = i_index_pointer - &idx[start] ;
			if(i_index_pointer != &idx[std::min(start+length,idx.size())])
				return (val[start+offset]) ;
			
			return 0 ;
		}
		else
		{
			return m[i][id] ;
		}
	}
	
	double operator *(const Vector&) const {return 0 ;};
	double operator *(const SymetricSparseVector&) const {return 0 ;};
	double operator *(const ConstSymetricSparseVector&) const {return 0 ;};
	Vector operator +(const Vector&) const {return Vector(0) ;};
	
	void print() const ;
	
} ;


} ;

Mu::CompositeSparseMatrix operator+(const Mu::CoordinateIndexedSparseMatrix  & sm, const Mu::CoordinateIndexedIncompleteSparseMatrix  & ism) ;

Mu::CompositeSparseMatrix operator+(const Mu::CoordinateIndexedIncompleteSparseMatrix  & ism, const Mu::CoordinateIndexedSparseMatrix  & sm ) ;

#endif

