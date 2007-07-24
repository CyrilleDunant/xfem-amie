//
// C++ Implementation: sparse_matrix
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "sparse_matrix.h"

using namespace Mu ;



CoordinateIndexedSparseMatrix::CoordinateIndexedSparseMatrix(std::map<std::pair<size_t, size_t>, double> &source)
{
	std::vector<double> temp_array ;
	std::vector<unsigned int> temp_column_index ;
	std::vector<unsigned int> temp_row_size ;
	std::map<std::pair<size_t, size_t>, double>::const_iterator previous = source.begin() ;
	size_t r_s = 0 ;
	for(std::map<std::pair<size_t, size_t>, double>::const_iterator ij = source.begin() ; ij != source.end() ; ++ij)
	{
		if(ij->first.first == previous->first.first)
		{
			r_s++;
			temp_array.push_back(ij->second) ;
			previous = ij ;
		}
		else
		{
			temp_row_size.push_back(r_s) ;
			r_s = 0 ;
			previous = ij ;
			temp_array.push_back(ij->second) ;
		}
	}
	
	array.resize(temp_array.size()) ; array = Vector(temp_array.size()) ;
	std::copy(temp_array.begin(), temp_array.end(), &array[0]) ;
	
	row_size.resize(temp_row_size.size()) ; row_size = std::valarray<unsigned int>(temp_row_size.size()) ;
	std::copy(temp_row_size.begin(), temp_row_size.end(), &row_size[0]) ;
}

CoordinateIndexedSparseMatrix::CoordinateIndexedSparseMatrix(std::map<std::pair<size_t, size_t>, Matrix> &source) : 
row_size((source.rbegin()->first.first+1)*source.begin()->second.numRows())
{
	size_t ddl = source.begin()->second.numRows() ;

	std::vector<double > temp_array ;
	std::vector<unsigned int> temp_column_index ;
	std::map<std::pair<size_t, size_t>, Matrix>::const_iterator previous = source.begin() ;
	std::vector<std::vector< double > > to_linerarise(ddl);
	std::vector<std::vector< unsigned int > > col_to_linerarise(ddl);

	
	for(std::map<std::pair<size_t, size_t>, Matrix>::const_iterator ij = source.begin() ; ij != source.end() ; )
	{
		size_t offset = ij->first.first*ddl ;
		for(size_t i = offset ; i < ddl+offset ; i++)
		{
			row_size[i] += ddl ;
		}
		
		if(ij->first.first != previous->first.first)
		{
			for(size_t i = 0 ; i < ddl ; i++)
			{
				for(size_t j = 0 ; j < to_linerarise[i].size() ; j++)
				{
					temp_array.push_back(to_linerarise[i][j]) ;
					temp_column_index.push_back(col_to_linerarise[i][j]) ;
				}
				to_linerarise[i].clear();
				col_to_linerarise[i].clear();
			}
			
			for(size_t i = 0 ; i < ddl ; i++)
			{
				for(size_t j = 0 ; j < ddl ; j++)
				{
					to_linerarise[i].push_back(ij->second[i][j]) ;
					col_to_linerarise[i].push_back(ij->first.second*ddl + j) ;
				}
			}
			
			previous = ij ;
		}
		else
		{
			for(size_t i = 0 ; i < ddl ; i++)
			{
				for(size_t j = 0 ; j < ddl ; j++)
				{
					to_linerarise[i].push_back(ij->second[i][j]) ;
					col_to_linerarise[i].push_back(ij->first.second*ddl + j) ;
				}
			}
			
			previous = ij ;
		}
		++ij ;
		if(ij == source.end())
		{

			for(size_t i = 0 ; i < ddl ; i++)
			{
				for(size_t j = 0 ; j < to_linerarise[i].size() ; j++)
				{
					temp_array.push_back(to_linerarise[i][j]) ;
					temp_column_index.push_back(col_to_linerarise[i][j]) ;
				}
				to_linerarise[i].clear();
				col_to_linerarise[i].clear();
			}
		}
	}
	
	array.resize(temp_array.size()) ;
	std::copy(temp_array.begin(), temp_array.end(), &array[0]) ;
	
	column_index.resize(temp_column_index.size());
	std::copy(temp_column_index.begin(), temp_column_index.end(), &column_index[0]) ;
}


CoordinateIndexedSparseMatrix::CoordinateIndexedSparseMatrix(const std::valarray<unsigned int> &rs, const std::valarray<unsigned int> &ci) : array(0., ci.size()),column_index(ci),row_size(rs), accumulated_row_size(rs.size())
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
	return SparseVector(array, column_index ,row_size[i], accumulated_row_size[i]) ;
}

const ConstSparseVector CoordinateIndexedSparseMatrix::operator[](const size_t i) const
{
// 	size_t start_index = std::accumulate(&row_size[0], &row_size[i+1], 0) ;
	return ConstSparseVector(array, column_index ,row_size[i],accumulated_row_size[i]) ;
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

Vector CoordinateIndexedSparseMatrix::operator *(const Vector v) const
{
	Vector ret(0., v.size()) ;
	
	for (size_t i = 0 ; i < row_size.size(); i++)
	{
		ret[i] += (*this)[i]*v ;
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
	Vector ret(0., row_size.size()) ;
	
	for(size_t i = 0 ; i < row_size.size() ; i++)
	{
// 		double v = (*this)[i][i] ;
// 		if(std::abs(v) < 1e16)
			ret[i] = 1./(*this)[i][i] ;
		
// 		else
// 			ret[i] = 1 ;
// 		std::cout << ret[i] << std::endl ;
	}
	return ret ;
}

Vector CoordinateIndexedSparseMatrix::inverseDiagonalSquared() const
{
	Vector ret(0., row_size.size()) ;
	
	for(size_t i = 0 ; i < row_size.size() ; i++)
	{
		ret[i] = 1./((*this)[i]*(*this)[i]) ;
		if(((*this)[i]*(*this)[i]) == 0)
			(*this)[i].print() ;
	}
	return ret ;
}


Vector CoordinateIndexedSparseMatrix::diagonal() const
{
	Vector ret(0., row_size.size()) ;
	
	for(size_t i = 0 ; i < row_size.size() ; i++)
	{
		ret[i] = (*this)[i][i] ;
// 		if(ret[i] == 0)
// 			std::cout << i << " is null !!!" << std::endl ;
	}
	return ret ;
}

CoordinateIndexedIncompleteSparseMatrix::CoordinateIndexedIncompleteSparseMatrix(const std::valarray<unsigned int> &ri, const std::valarray<unsigned int> &ci) : array(0.,ci.size() ), column_index(ci), row_index(ri) { } ;
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

		ret[r_index] += array[array_index] * v[c_index] ;
		array_index++ ;
	}
	return ret ;
}
	
SparseVector CoordinateIndexedIncompleteSparseMatrix::operator[](const size_t i)
{
	unsigned int * start_pointer = std::find(&row_index[0], &row_index[row_index.size()], i) ;
	unsigned int start_index = start_pointer - &row_index[0] ;
	unsigned int length = 0 ;
	for(unsigned int j = start_index ; j < array.size() ; j++)
	{
		if(row_index[j] != row_index[start_index])
			break ;
		else
			length++ ;
	}
	return SparseVector(array, column_index, length,start_index ) ;
}

ConstSparseVector CoordinateIndexedIncompleteSparseMatrix::operator[](const size_t i) const
{
	const unsigned int * start_pointer = std::find(&row_index[0], &row_index[row_index.size()], i) ;
	unsigned int start_index = start_pointer - &row_index[0] ;
	unsigned int length = 0 ;
	for(unsigned int j = start_index ; j < array.size() ; j++)
	{
		if(row_index[j] != row_index[start_index])
			break ;
		else
			length++ ;
	}
	return ConstSparseVector(array, column_index, length, start_index) ;
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
	
	
	CoordinateIndexedSparseMatrix T(rs, ci) ;
	
	for(size_t i = 0 ; i < row_size.size() ; i++)
	{
		for(size_t j = 0 ; j < row_size[i] ; j++)
		{
			T[column_index[accumulated_row_size[i]+j]][i] = array[accumulated_row_size[i]+j] ;
		}
	}
	return T ;
}


SymetricSparseMatrix::SymetricSparseMatrix(std::map<std::pair<size_t, size_t>, double> &source)
{
	std::vector<double> temp_array ;
	std::vector<size_t> temp_column_index ;
	std::vector<size_t> temp_row_size ;
	std::map<std::pair<size_t, size_t>, double>::const_iterator previous = source.begin() ;
	size_t r_s = 0 ;
	for(std::map<std::pair<size_t, size_t>, double>::const_iterator ij = source.begin() ; ij != source.end() ; ++ij)
	{
		if(ij->first.first == previous->first.first && (ij->first.first >= ij->first.second))
		{
			r_s++;
			temp_array.push_back(ij->second) ;
			previous = ij ;
		}
		else if(ij->first.first >= ij->first.second)
		{
			temp_row_size.push_back(r_s) ;
			r_s = 0 ;
			previous = ij ;
			temp_array.push_back(ij->second) ;
		}
	}
	
	array.resize(temp_array.size()) ; array = Vector(temp_array.size()) ;
	std::copy(temp_array.begin(), temp_array.end(), &array[0]) ;
	
	row_size.resize(temp_row_size.size()) ; row_size = std::valarray<size_t>(temp_row_size.size()) ;
	std::copy(temp_row_size.begin(), temp_row_size.end(), &row_size[0]) ;
}
	
SymetricSparseMatrix::SymetricSparseMatrix(std::map<std::pair<size_t, size_t>, Matrix> &source)
{
	size_t ddl = source.begin()->second.numRows() ;

	std::vector<double > temp_array ;
	std::vector<size_t> temp_column_index ;
	std::map<std::pair<size_t, size_t>, Matrix>::const_iterator previous = source.begin() ;
	std::vector<std::vector< double > > to_linerarise(ddl);
	std::vector<std::vector< size_t > > col_to_linerarise(ddl);

	
	for(std::map<std::pair<size_t, size_t>, Matrix>::const_iterator ij = source.begin() ; ij != source.end() ; )
	{
		size_t offset = ij->first.first*ddl ;
		for(size_t i = offset ; i < ddl+offset ; i++)
		{
			row_size[i] += ddl ;
		}
		
		if(ij->first.first != previous->first.first && (ij->first.first >= ij->first.second))
		{
			for(size_t i = 0 ; i < ddl ; i++)
			{
				for(size_t j = i ; j < to_linerarise[i].size() ; j++)
				{
					temp_array.push_back(to_linerarise[i][j]) ;
					temp_column_index.push_back(col_to_linerarise[i][j]) ;
				}
				to_linerarise[i].clear();
				col_to_linerarise[i].clear();
			}
			
			for(size_t i = 0 ; i < ddl ; i++)
			{
				for(size_t j = i ; j < ddl ; j++)
				{
					to_linerarise[i].push_back(ij->second[i][j]) ;
					col_to_linerarise[i].push_back(ij->first.second*ddl + j) ;
				}
			}
			
			previous = ij ;
		}
		else if (ij->first.first >= ij->first.second)
		{
			for(size_t i = 0 ; i < ddl ; i++)
			{
				for(size_t j = i ; j < ddl ; j++)
				{
					to_linerarise[i].push_back(ij->second[i][j]) ;
					col_to_linerarise[i].push_back(ij->first.second*ddl + j) ;
				}
			}
			
			previous = ij ;
		}
		++ij ;
		if(ij == source.end())
		{

			for(size_t i = 0 ; i < ddl ; i++)
			{
				for(size_t j = i ; j < to_linerarise[i].size() ; j++)
				{
					temp_array.push_back(to_linerarise[i][j]) ;
					temp_column_index.push_back(col_to_linerarise[i][j]) ;
				}
				to_linerarise[i].clear();
				col_to_linerarise[i].clear();
			}
		}
	}
	
	array.resize(temp_array.size()) ;
	std::copy(temp_array.begin(), temp_array.end(), &array[0]) ;
	
	column_index.resize(temp_column_index.size());
	std::copy(temp_column_index.begin(), temp_column_index.end(), &column_index[0]) ;
}

	SymetricSparseMatrix::SymetricSparseMatrix(const std::valarray<size_t> &rs, const std::valarray<size_t> &ci) : array(0., ci.size()),column_index(ci),row_size(rs), accumulated_row_size((size_t)0, rs.size())
{
// 	accumulated_row_size[1] = row_size[0] ;
	for(size_t i = 1 ; i < accumulated_row_size.size() ; i++)
	{
		accumulated_row_size[i] += accumulated_row_size[i-1]+row_size[i-1] ;
	}
}

SymetricSparseMatrix::~SymetricSparseMatrix(){ } ;
	
SymetricSparseVector SymetricSparseMatrix::operator[](const size_t i)
{
	return SymetricSparseVector(array, column_index ,row_size[i], accumulated_row_size[i]) ;
}

const ConstSymetricSparseVector SymetricSparseMatrix::operator[](const size_t i) const
{
	return ConstSymetricSparseVector(i, array, column_index ,row_size[i], accumulated_row_size[i], *this) ;
}

double & SymetricSparseMatrix::operator()(const size_t i, const size_t j)
{
	if(i <= j)
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
	else
	{
		size_t start_index = accumulated_row_size[j] ;
		for(size_t k = 0 ; k < row_size[i] ;k++)
		{
			if(column_index[start_index] == i)
				break ;
			else
				start_index++ ;
		}
		return array[start_index] ;
	}
}

double  SymetricSparseMatrix::operator()(const size_t i, const size_t j) const
{
		if(i <= j)
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
	else
	{
		size_t start_index = accumulated_row_size[j] ;
		for(size_t k = 0 ; k < row_size[i] ;k++)
		{
			if(column_index[start_index] == i)
				break ;
			else
				start_index++ ;
		}
		return array[start_index] ;
	}
}
	
Vector SymetricSparseMatrix::operator *(const Vector v) const 
{
	Vector ret(0., v.size()) ;
	
	for (size_t i = 0 ; i < row_size.size(); i++)
	{
		ret[i] += (*this)[i]*v ;
	}
	return ret ;

}

	SymetricSparseMatrix SymetricSparseMatrix::operator +(const CoordinateIndexedIncompleteSparseMatrix &M) const 
{
	SymetricSparseMatrix ret(*this) ;
	
	for(size_t i = 0 ; i < M.array.size() ; i++)
	{
		if(M.row_index[i]>=M.column_index[i])
			ret(M.row_index[i], M.column_index[i]) += M.array[i] ;
	}
	
	return ret ;
}

	SymetricSparseMatrix & SymetricSparseMatrix::operator +=(const CoordinateIndexedIncompleteSparseMatrix &M)
{
		for(size_t i = 0 ; i < M.array.size() ; i++)
	{
		if( M.row_index[i] >= M.column_index[i] )
			(*this)[M.row_index[i]][M.column_index[i]] += M.array[i] ;
	}
	
	return *this ;
}

	SymetricSparseMatrix SymetricSparseMatrix::operator +(const CoordinateIndexedSparseMatrix &M) const
{
		SymetricSparseMatrix ret(*this) ;
	
	size_t current_row = 0 ;
	size_t current_column = 0 ;
	for(size_t i = 0 ; i < M.array.size() ; i++)
	{
		size_t temp_column = M.column_index[i] ;
		if(current_row <= temp_column)
			ret.array[i] += M.array[i] ;
		if(temp_column <= current_column)
			current_row++ ;
		current_column = temp_column ;
	}
	
	return ret ;
}

	SymetricSparseMatrix & SymetricSparseMatrix::operator=(const CoordinateIndexedSparseMatrix &M)
{

	size_t current_row = 0 ;
	size_t current_column = 0 ;
	for(size_t i = 0 ; i < M.array.size() ; i++)
	{
		size_t temp_column = M.column_index[i] ;
		if(current_row <= temp_column)
			array[i] = M.array[i] ;
		if(temp_column <= current_column)
			current_row++ ;
		current_column = temp_column ;
	}
	
	return *this ;
}

	SymetricSparseMatrix SymetricSparseMatrix::operator +(const SymetricSparseMatrix &M) const
{
	SymetricSparseMatrix ret(*this) ;
	
	for(size_t i = 0 ; i < M.array.size() ; i++)
	{
		ret.array[i] += M.array[i] ;
	}
	
	return ret ;
}
	
	Vector SymetricSparseMatrix::inverseDiagonal() const
{
	Vector ret(0., row_size.size()) ;
	
	for(size_t i = 0 ; i < row_size.size() ; i++)
	{
		ret[i] = 1./array[accumulated_row_size[i]] ;
	}
	return ret ;
}

	Vector SymetricSparseMatrix::inverseDiagonalSquared() const
{
		Vector ret(0., row_size.size()) ;
	
	for(size_t i = 0 ; i < row_size.size() ; i++)
	{
		double d = array[accumulated_row_size[i]] ;
		ret[i] = 1./(d*d) ;
	}
	return ret ;
}

	Vector SymetricSparseMatrix::diagonal() const
{
		Vector ret(0., row_size.size()) ;
	
	for(size_t i = 0 ; i < row_size.size() ; i++)
	{
		ret[i] = array[accumulated_row_size[i]] ;
	}
	return ret ;
}

	const SymetricSparseMatrix & SymetricSparseMatrix::transpose() const
{
	return *this ;
}
	
	double SymetricSparseMatrix::froebeniusNorm() const
{
	return 2.*sqrt(std::inner_product(&array[0], &array[array.size()], &array[0], (double)(0))) ;
}

double SymetricSparseMatrix::infinityNorm() const
{
	return std::abs(array).max() ;
}
	
void ConstBandSparseVector::print() const { } ;


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
		ret[i] = sm.sm[i]*ve +sm.ism[i]*ve;
	}
	
	return ret ;
} 


CompositeSparseMatrixTimesVecMinusVecMinusVec CompositeSparseMatrixTimesVecMinusVec::operator - (const Vector & v) const
{
	return CompositeSparseMatrixTimesVecMinusVecMinusVec(*this, v) ;
}
