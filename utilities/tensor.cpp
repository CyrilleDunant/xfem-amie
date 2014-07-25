#include "tensor.h"

using namespace Mu ;

Tensor::Tensor(size_t o, size_t d) : order(o), dim(d) 
{
	size_t s = 1 ;
	for(size_t i = 0 ; i < o ; i++)
		s *= dim ;
	components.resize(s) ;
	components = 0. ;
}

Tensor::Tensor(const Point & p, std::vector<Variable> v) : order(1), dim(v.size())
{
	size_t s = 1 ;
	for(size_t i = 0 ; i < order ; i++)
		s *= dim ;
	components.resize(s) ;
	components = 0. ;
	for(int i = 0 ; i < v.size() ; i++)
	{
		double d = 0 ;
		switch(v[i])
		{
		  case XI:
		    d = p.getX() ; 
		    break ;
		  case ETA:
		    d = p.getY() ; 
		    break ;
		  case ZETA:
		    d = p.getZ() ; 
		    break ;
		  case TIME_VARIABLE:
		    d = p.getT() ; 
		    break ;
		}
		components[i] = d ;
	}
}

Tensor::Tensor(const Matrix & m, bool transpose) : order(2), dim(m.numRows())
{
	components.resize(dim*dim) ;
	for(int i = 0 ; i < dim ; i++)
	{
		for(int j = 0 ; j < dim ; j++)
		{
			if(transpose)
				(*this)(j,i) = m[i][j] ;
			else
				(*this)(i,j) = m[i][j] ;
		}
	}
}

size_t Tensor::index( std::vector<int> indexes ) const
{
	if(indexes.size() == 1)
		return indexes[0] ;
	size_t i = 0 ;
	size_t s = 1 ;
	for(int k = order ; k > 0 ; k--)
	{
		i += (indexes[k-1]) * s ;
		s *= dim ;
	}
	return i ;
}

size_t Tensor::index( int i ) const 
{
	std::vector<int> idx ; idx.push_back(i) ;
	return index(idx) ;
}

size_t Tensor::index( int i, int j ) const 
{
	std::vector<int> idx ; 
	idx.push_back(i) ;
	idx.push_back(j) ;
	return index(idx) ;
}

size_t Tensor::index( int i, int j, int k ) const 
{
	std::vector<int> idx ; 
	idx.push_back(i) ;
	idx.push_back(j) ;
	idx.push_back(k) ;
	return index(idx) ;
}

size_t Tensor::index( int i, int j, int k, int l ) const 
{
	std::vector<int> idx ; 
	idx.push_back(i) ;
	idx.push_back(j) ;
	idx.push_back(k) ;
	idx.push_back(l) ;	
	return index(idx) ;
}

size_t Tensor::index( int i, int j, int k, int l, int m ) const 
{
	std::vector<int> idx ; 
	idx.push_back(i) ;
	idx.push_back(j) ;
	idx.push_back(k) ;
	idx.push_back(l) ;
	idx.push_back(m) ;
	return index(idx) ;
}

size_t Tensor::index( int i, int j, int k, int l, int m, int n ) const 
{
	std::vector<int> idx ; 
	idx.push_back(i) ;
	idx.push_back(j) ;
	idx.push_back(k) ;
	idx.push_back(l) ;
	idx.push_back(m) ;
	idx.push_back(n) ;
	return index(idx) ;
}


void Tensor::print() const
{
	for(size_t i = 0 ; i < components.size() ; i++)
	{
		std::cout << components[i] << "\t" ;
	}
	std::cout << std::endl ;
}

void Tensor::threshold(double thr) 
{
	for(size_t i = 0 ; i < components.size() ; i++)
	{
		if(std::abs(components[i]) < thr)
			  components[i] = 0 ;
	}
	
}

Matrix Tensor::toMatrix(int d1, int d2) const
{
	Matrix ret ;
	if(d1+d2 != order)
		return ret ;

	Tensor dummy1( d1, dim ) ;
	Tensor dummy2( d2, dim ) ;
	
	ret = Matrix(dummy1.size(), dummy2.size()) ;
	for(int i = 0 ; i < dummy1.size() ; i++)
	{
		for(int j = 0 ; j < dummy2.size() ; j++)
		{
			ret[i][j] = components[ i*dummy2.size() + j ] ;
		}
	}
	
	return ret ;
	
}