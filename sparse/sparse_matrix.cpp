//
// C++ Implementation: sparse_matrix
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "sparse_matrix.h"
#include <limits>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include <sys/time.h>
using namespace Mu ;



CoordinateIndexedSparseMatrix::CoordinateIndexedSparseMatrix(std::map<std::pair<size_t, size_t>, double> &source, size_t s): stride(s), array(source.size()*s*(s+s%2)/(s*s)), column_index(source.size()/s), row_size(source.rbegin()->first.first/s), accumulated_row_size(source.rbegin()->first.first/s)
{
	unsigned int * column_index_iterator = &column_index[0] ;
	unsigned int * row_size_iterator = &row_size[0] ;
	std::map<std::pair<size_t, size_t>, double>::const_iterator previous = source.begin() ;
	size_t r_s = 0 ;
	for(std::map<std::pair<size_t, size_t>, double>::const_iterator ij = source.begin() ; ij != source.end() ; ++ij)
	{
		if(ij->first.first/stride == previous->first.first/stride)
		{
			r_s++;
		}
		else
		{
			*row_size_iterator = r_s/stride ;
			row_size_iterator++ ;
			r_s = 1 ;	
		}
		previous = ij ;
		*column_index_iterator =ij->first.second/stride ;
		if((ij->first.second+1)%s == 0)
			column_index_iterator++ ;
	}
	*row_size_iterator = r_s/stride ;
	
	for(std::map<std::pair<size_t, size_t>, double>::const_iterator ij = source.begin() ; ij != source.end() ; ++ij)
		(*this)[ij->first.first][ij->first.second] = ij->second ;

	for(size_t i = 1 ; i < accumulated_row_size.size() ; i++)
	{
		accumulated_row_size[i] += accumulated_row_size[i-1]+row_size[i-1] ;
	}
}

CoordinateIndexedSparseMatrix::CoordinateIndexedSparseMatrix(std::map<std::pair<size_t, size_t>, Matrix> &source) : stride(source.begin()->second.numRows()),  array(source.size()*stride*(stride+stride%2)/(stride*stride)), column_index(source.size()/stride), row_size(source.rbegin()->first.first/stride), accumulated_row_size(source.rbegin()->first.first/stride)
{
	double * array_iterator = &array[0] ;
	unsigned int * column_index_iterator = &column_index[0] ;
	unsigned int * row_size_iterator = &row_size[0] ;
	std::map<std::pair<size_t, size_t>, Matrix>::const_iterator previous = source.begin() ;
	size_t r_s = 0 ;
	for(std::map<std::pair<size_t, size_t>, Matrix>::const_iterator ij = source.begin() ; ij != source.end() ; ++ij)
	{
		if(ij->first.first == previous->first.first)
		{
			r_s++;
		}
		else
		{
			*row_size_iterator = r_s ;
			row_size_iterator++ ;
			r_s = 1 ;	
		}
		previous = ij ;
		*column_index_iterator =ij->first.second ;
		column_index_iterator++ ;
	}
	*row_size_iterator = r_s ;
	array = 0 ;
	for(std::map<std::pair<size_t, size_t>, Matrix>::const_iterator ij = source.begin() ; ij != source.end() ; ++ij)
	{
		for(size_t i = 0 ; i < stride ; i++)
		{
			for(size_t j = 0 ; j < stride ; j++)
			{
				*array_iterator = ij->second[j][i] ;
				array_iterator++ ;
			}
		}
	}
	
	for(size_t i = 1 ; i < accumulated_row_size.size() ; i++)
	{
		accumulated_row_size[i] += accumulated_row_size[i-1]+row_size[i-1] ;
	}
}


CoordinateIndexedSparseMatrix::CoordinateIndexedSparseMatrix(const std::valarray<unsigned int> &rs, const std::valarray<unsigned int> &ci, size_t s) : stride(s), array(0., ci.size()*stride*(stride+stride%2)),column_index(ci),row_size(rs), accumulated_row_size(rs.size())
{
// 	accumulated_row_size[1] = row_size[0] ;
	for(size_t i = 1 ; i < accumulated_row_size.size() ; i++)
	{
		accumulated_row_size[i] += accumulated_row_size[i-1]+row_size[i-1] ;
	}
}

CoordinateIndexedSparseMatrix::~CoordinateIndexedSparseMatrix()
{
}

double CoordinateIndexedSparseMatrix::froebeniusNorm() const
{
	return sqrt(std::inner_product(&array[0], &array[array.size()], &array[0], (double)(0))) ;
}

double CoordinateIndexedSparseMatrix::infinityNorm() const
{
	return std::abs(array).max() ;
}

SparseVector CoordinateIndexedSparseMatrix::operator[](const size_t i)
{
// 	size_t start_index = std::accumulate(&row_size[0], &row_size[i+1], 0) ;
	return SparseVector(array, column_index ,row_size[i/stride], accumulated_row_size[i/stride], i, stride) ;
}

const ConstSparseVector CoordinateIndexedSparseMatrix::operator[](const size_t i) const
{
// 	size_t start_index = std::accumulate(&row_size[0], &row_size[i+1], 0) ;
	return ConstSparseVector(array, column_index ,row_size[i/stride],accumulated_row_size[i/stride], i, stride) ;
}

double & CoordinateIndexedSparseMatrix::operator()(const size_t i, const size_t j)
{
	size_t start_index = accumulated_row_size[i] ;
	for(size_t k = 0 ; k < row_size[i] ;k++)
	{
		if(column_index[start_index] == j)
			break ;
		else
			start_index++ ;
	}
	return array[start_index] ;
}

double CoordinateIndexedSparseMatrix::operator()(const size_t i, const size_t j) const 
{
	size_t start_index = accumulated_row_size[i] ;
	
	for(size_t k = 0 ; k < row_size[i] ;k++)
	{
		if(column_index[start_index] == j)
			break ;
		else
			start_index++ ;
	}
	return array[start_index] ;
}

CoordinateIndexedSparseMatrixTimesVec CoordinateIndexedSparseMatrix::operator *(const Vector & v) const
{
	return CoordinateIndexedSparseMatrixTimesVec(*this, v) ;
}

CoordinateIndexedSparseMatrixTimesVecPlusVec CoordinateIndexedSparseMatrixTimesVec::operator + (const Vector & v) const
{
	return CoordinateIndexedSparseMatrixTimesVecPlusVec(*this, v) ;
}

CoordinateIndexedSparseMatrixTimesVecMinusVec CoordinateIndexedSparseMatrixTimesVec::operator - (const Vector & v) const
{
	return CoordinateIndexedSparseMatrixTimesVecMinusVec(*this, v) ;
}


CoordinateIndexedSparseMatrixTimesVec::operator const Vector() const
{
	Vector ret(0., ve.size()) ;
	
	for (size_t i = 0 ; i < sm.row_size.size(); i++)
	{
		Vector temp = sm[i*sm.stride]*ve ;

		for(size_t j = i*sm.stride ; j < i*sm.stride+sm.stride ; j++)
			ret[j] = temp[j-i*sm.stride] ;
	}
	return ret ;
}


CompositeSparseMatrix::operator CoordinateIndexedSparseMatrix() const
{
	CoordinateIndexedSparseMatrix ret(sm) ;
	
	for(size_t i = 0 ; i < ism.array.size() ; i++)
	{
		ret[ism.row_index[i]][ism.column_index[i]] += ism.array[i] ;
	}
	
	return ret ;
}

CoordinateIndexedSparseMatrix & CoordinateIndexedSparseMatrix::operator +=(const CoordinateIndexedIncompleteSparseMatrix & M)
{
	for(size_t i = 0 ; i < M.array.size() ; i++)
	{
// 		size_t * idx = std::find(&column_index[accumulated_row_size[M.row_index[i]]], &column_index[accumulated_row_size[M.row_index[i]+1]], M.column_index[i]) ;
// 		size_t * idx_0 = &column_index[accumulated_row_size[M.row_index[i]]] ;
// 		size_t delta = idx -idx_0 ;
// 		array[column_index[accumulated_row_size[M.row_index[i]]] + delta] += M.array[i] ;
		(*this)[M.row_index[i]][M.column_index[i]] += M.array[i] ;
	}
	
	return *this ;
}


// CoordinateIndexedSparseMatrix CoordinateIndexedSparseMatrix::operator +(const CoordinateIndexedSparseMatrix & M) const
// {
// 	CoordinateIndexedSparseMatrix ret(*this) ;
// 	
// 	for(size_t i = 0 ; i < M.array.size() ; i++)
// 	{
// 		ret.array[i] += M.array[i] ;
// 	}
// 	
// 	return ret ;
// }

CompositeSparseMatrix operator+(const CoordinateIndexedSparseMatrix  & sm, const CoordinateIndexedIncompleteSparseMatrix  & ism) 
{
	return CompositeSparseMatrix(sm, ism) ;
} 

CompositeSparseMatrix operator+(const CoordinateIndexedIncompleteSparseMatrix  & ism, const CoordinateIndexedSparseMatrix  & sm )
{
	return CompositeSparseMatrix(sm, ism) ;
} ;



Vector CoordinateIndexedSparseMatrix::inverseDiagonal() const
{
	Vector ret(0., row_size.size()*stride) ;
	
	for(size_t i = 0 ; i < row_size.size()*stride ; i++)
	{
		double v = (*this)[i][i] ;
		if(std::abs(v) > 1e-6)
			ret[i] = 1./v ;
		else if(v > 0)
		{
			ret[i] = 1./1e-6 ;
		}
		else
		{
			ret[i] = -1./1e-6 ;
		}
	}
		
// 		else
// 			ret[i] = 1 ;
// 		std::cout << ret[i] << std::endl ;
// 	}
	return ret ;
}

Vector CoordinateIndexedSparseMatrix::inverseDiagonalSquared() const
{
	Vector ret(0., row_size.size()*stride) ;
	
	for(size_t i = 0 ; i < row_size.size()*stride ; i++)
	{
		ret[i] = 1./((*this)[i][i]*(*this)[i][i]) ;
		if(((*this)[i][i]*(*this)[i][i]) == 0)
			(*this)[i].print() ;
	}
	return ret ;
}


Vector CoordinateIndexedSparseMatrix::diagonal() const
{
	Vector ret(0., row_size.size()*stride) ;
	
	for(size_t i = 0 ; i < row_size.size()*stride ; i++)
	{
		ret[i] = (*this)[i][i] ;
// 		if(ret[i] == 0)
// 			std::cout << i << " is null !!!" << std::endl ;
	}
	return ret ;
}

CoordinateIndexedIncompleteSparseMatrix::CoordinateIndexedIncompleteSparseMatrix(const std::valarray<unsigned int> &ri, const std::valarray<unsigned int> &ci) : array(0.,ci.size() ), column_index(ci), row_index(ri), stride(1) { } ;
CoordinateIndexedIncompleteSparseMatrix::~CoordinateIndexedIncompleteSparseMatrix() { } ;
	

CoordinateIndexedIncompleteSparseMatrix::CoordinateIndexedIncompleteSparseMatrix(const std::valarray<unsigned int> &ri, const Vector &ci, const size_t id)
{
	
	std::map<std::pair<size_t, size_t>, double > ids ;
	
	for(size_t i  =  0 ; i < ri.size() ; i++)
	{
		ids[std::make_pair(id, ri[i])] = ci[i] ;
	}
	for(size_t i  =  0 ; i < ri.size() ; i++)
	{
		ids[std::make_pair(ri[i], id)] = ci[i] ;
	}
	
	array.resize(ids.size()) ;
	column_index.resize(ids.size()) ;
	row_index.resize(ids.size()) ;
	
	size_t idx = 0 ;
	for(std::map<std::pair<size_t, size_t>, double >::iterator i = ids.begin() ; i != ids.end() ; ++i)
	{
		row_index[idx] = (*i).first.first ;
		column_index[idx] = (*i).first.second ;
		array[idx] = (*i).second ;
		idx++ ;
	}
	
	stride = 1 ;
} ;

double & CoordinateIndexedIncompleteSparseMatrix::operator()(size_t i, size_t j)
{
	unsigned int * start_pointer = std::find(&row_index[0], &row_index[row_index.size()], i) ;
	unsigned int start_index = start_pointer - &row_index[0] ;
	unsigned int * index_pointer = std::find(&column_index[start_index], &column_index[row_index.size()], j) ;
	unsigned int index = index_pointer - &column_index[start_index] ;
	
	return array[index] ;
}
	
Vector CoordinateIndexedIncompleteSparseMatrix::operator *(const Vector v) const
{
	Vector ret(0., v.size()) ;
	
	for (unsigned int array_index = 0 ; array_index < array.size(); array_index++)
	{
		size_t r_index = row_index[array_index];
		size_t c_index = column_index[array_index];

		ret[r_index] = fma(array[array_index] , v[c_index], ret[r_index]) ;
		array_index++ ;
	}
	return ret ;
}
	
SparseVector CoordinateIndexedIncompleteSparseMatrix::operator[](const size_t i)
{
	unsigned int * start_pointer = std::find(&row_index[0], &row_index[row_index.size()], i/stride) ;
	unsigned int start_index = start_pointer - &row_index[0] ;
	unsigned int length = 0 ;
	for(unsigned int j = start_index ; j < array.size() ; j++)
	{
		if(row_index[j] != row_index[start_index])
			break ;
		else
			length++ ;
	}
	return SparseVector(array, column_index, length,start_index , i, stride) ;
}

ConstSparseVector CoordinateIndexedIncompleteSparseMatrix::operator[](const size_t i) const
{
	const unsigned int * start_pointer = std::find(&row_index[0], &row_index[row_index.size()], i/stride) ;
	unsigned int start_index = start_pointer - &row_index[0] ;
	unsigned int length = 0 ;
	for(unsigned int j = start_index ; j < array.size() ; j++)
	{
		if(row_index[j] != row_index[start_index])
			break ;
		else
			length++ ;
	}
	return ConstSparseVector(array, column_index, length, start_index, i, stride) ;
}


CoordinateIndexedSparseMatrix & CoordinateIndexedSparseMatrix::operator=(const CoordinateIndexedSparseMatrix &S)
{
	this->column_index.resize(S.column_index.size()) ;
	this->array.resize(S.array.size()) ;
	this->row_size.resize(S.row_size.size()) ;
	this->accumulated_row_size.resize(S.accumulated_row_size.size()) ;
	
	this->column_index=S.column_index ;
	this->array=S.array ;
	this->row_size=S.row_size ;
	this->accumulated_row_size = S.accumulated_row_size ;
	
	return *this ;
}

CoordinateIndexedSparseMatrix CoordinateIndexedSparseMatrix::transpose() const
{
	std::valarray<unsigned int> rs((unsigned int)0, row_size.size()) ;
	std::valarray<unsigned int> ci((unsigned int)0,column_index.size()) ;
	std::valarray<unsigned int> pos_at_row((unsigned int)0,row_size.size()) ;
	std::valarray<unsigned int> ars((unsigned int)0, accumulated_row_size.size()) ;
	
	for(size_t i = 0 ; i < column_index.size() ; i++)
	{
		rs[column_index[i]]++ ;
	}
	
	size_t imax  =  ars.size() ;
	for(size_t i = 1 ; i < imax ; i++)
	{
		ars[i] = ars[i] + ars[i-1] + rs[i-1] ;
	}
	
	
	for(size_t i = 0 ; i < row_size.size() ; i++)
	{
		for(size_t j = 0 ; j < row_size[i] ; j++)
		{
			ci[ars[column_index[accumulated_row_size[i]+j]] + pos_at_row[column_index[accumulated_row_size[i]+j]]++] = i ;
		}
	}
	
	
	CoordinateIndexedSparseMatrix T(rs, ci, stride) ;
	
	for(size_t i = 0 ; i < row_size.size() ; i++)
	{
		for(size_t j = 0 ; j < row_size[i] ; j++)
		{
			T[column_index[accumulated_row_size[i]+j]][i] = array[accumulated_row_size[i]+j] ;
		}
	}
	return T ;
}


CompositeSparseMatrixTimesVec CompositeSparseMatrix::operator *(const Vector v) const 
{
	return CompositeSparseMatrixTimesVec(*this, v) ;
}


CompositeSparseMatrixTimesVecPlusVec CompositeSparseMatrixTimesVec::operator + (const Vector & v) const 
{
	return CompositeSparseMatrixTimesVecPlusVec(*this, v) ;
}

CompositeSparseMatrixTimesVecMinusVec CompositeSparseMatrixTimesVec::operator - (const Vector & v) const 
{
	return CompositeSparseMatrixTimesVecMinusVec(*this, v) ;
}

CompositeSparseMatrixTimesVec::operator const Vector() const 
{
	Vector ret(0., ve.size()) ;
	
	for (size_t i = 0 ; i < sm.sm.row_size.size(); i++)
	{
		Vector temp = sm.sm[i*sm.sm.stride]*ve;
		for(size_t j = i*sm.sm.stride ; j < i*sm.sm.stride+sm.sm.stride ; j++)
			ret[j] = temp[j-i*sm.sm.stride]+(sm.ism[i*sm.sm.stride+j-i*sm.sm.stride]*ve)[0];
	}
	
	return ret ;
} 

CompositeSparseMatrixTimesVecMinusVecMinusVec CompositeSparseMatrixTimesVecMinusVec::operator - (const Vector & v) const
{
	return CompositeSparseMatrixTimesVecMinusVecMinusVec(*this, v) ;
}

void Mu::assign(Vector & ret, const Mu::CoordinateIndexedSparseMatrixTimesVecPlusVec & c, const int rowstart, const int colstart)
{

	size_t stride = c.co.sm.stride ;
	int end = c.co.sm.row_size.size()*stride ;
	
	ret = c.ve ;
	if(rowstart)
		ret[ std::slice(0,rowstart,1) ] = 0. ;

	#pragma omp parallel for schedule(static)
	for (int i = rowstart ; i < end ; i+=stride)
	{
		c.co.sm[i].inner_product(c.co.ve, &ret[i], rowstart,  colstart);
	}
} ;

void Mu::assign(Vector & ret, const Mu::CoordinateIndexedSparseMatrixTimesVecMinusVec & c, const int rowstart, const int colstart)
{

	size_t stride = c.co.sm.stride ;
	int end = c.co.sm.row_size.size()*stride ; 
	const Vector & ve = c.co.ve ;
	ret = -c.ve ;
	if(rowstart)
		ret[ std::slice(0,rowstart,1) ] = 0. ;
	double * begin = &ret[0] ;
	int i ;
	#pragma omp parallel for schedule(static)
	for (int i = rowstart ; i <end ; i+=stride)
	{
		c.co.sm[i].inner_product(ve, begin+i, rowstart, colstart);
	}

} ;

void Mu::assign(Vector & ret, const Mu::CoordinateIndexedSparseMatrixTimesVec & c, const int rowstart, const int colstart)
{

	size_t stride = c.sm.stride ;
	int end = c.sm.row_size.size()*stride ;
	const Vector & ve = c.ve ;
	ret = 0. ;

	double * begin = &ret[0] ;
	int i ;
	#pragma omp parallel for schedule(static)
	for (i = rowstart ; i < end; i+=stride)
	{
		c.sm[i].inner_product(ve, begin+i, rowstart, colstart);
	}
	
	if(rowstart)
		ret[ std::slice(0,rowstart,1) ] = 0. ;
} ;
