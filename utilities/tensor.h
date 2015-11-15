#ifndef TENSOR_H
#define TENSOR_H

#include "matrixops.h"
#include "../geometry/geometry_base.h" 
#include "../polynomial/variable.h" 

namespace Amie
{

typedef enum
{
	PLANE_STRESS,
	PLANE_STRAIN,
	PLANE_STRESS_FREE_G,
} planeType ; 

  
class Tensor
{
protected:
	size_t order ;
	size_t dim ;
	Vector components ;
	
public:
	Tensor(size_t o, size_t d) ;
	Tensor(const Point & p, std::vector<Variable> v) ;
	Tensor(const Matrix & m, bool transpose = false) ;
	
	size_t getOrder() const { return order ; }
	size_t getDimensions() const { return dim ; }

	size_t index( std::vector<int> indexes ) const ;
	size_t index( int i ) const ;
	size_t index( int i, int j ) const ;
	size_t index( int i, int j, int k ) const ;
	size_t index( int i, int j, int k, int l ) const ;
	size_t index( int i, int j, int k, int l, int m ) const ;
	size_t index( int i, int j, int k, int l, int m, int n ) const ;
	size_t size() const { return components.size() ; }
	
	double & getComponent( std::vector<int> indexes ) { return components[ index(indexes) ] ;  }
	double getComponent( std::vector<int> indexes ) const { return components[ index(indexes) ] ; }
	
	double & operator() (int i) { return components[ index(i) ] ; }
	double & operator() (int i, int j) { return components[ index(i,j) ] ; }
	double & operator() (int i, int j, int k) { return components[ index(i,j, k) ] ; }
	double & operator() (int i, int j, int k, int l) { return components[ index(i,j,k,l) ] ; }
	double & operator() (int i, int j, int k, int l, int m){ return components[ index(i,j,k,l,m) ] ; }
	double & operator() (int i, int j, int k, int l, int m, int n){ return components[ index(i,j,k,l,m,n) ] ; }

	double operator() (int i) const { return components[ index(i) ] ; }
	double operator() (int i, int j) const { return components[ index(i,j) ] ; }
	double operator() (int i, int j, int k) const { return components[ index(i,j,k) ] ; }
	double operator() (int i, int j, int k, int l) const { return components[ index(i,j,k,l) ] ; }
	double operator() (int i, int j, int k, int l, int m) const { return components[ index(i,j,k,l,m) ] ; }
	double operator() (int i, int j, int k, int l, int m, int n) const { return components[ index(i,j,k,l,m,n) ] ; }
	
	Matrix toMatrix(int d1, int d2) const ;
	
	void print() const ;
	
	void threshold(double thr) ;

	static Matrix cauchyGreen(std::pair<double,double> prop, bool hooke, SpaceDimensionality dim, planeType pt ) ;
	static Matrix cauchyGreen(double p1, double p2, bool hooke, SpaceDimensionality dim, planeType pt ) ;
	static Matrix orthotropicCauchyGreen(double E_1, double E_2, double G,  double nu, planeType pt) ;
	static Matrix orthotropicCauchyGreen(double E_1, double E_2, double G,  double nu, double angle, planeType pt) ;
	static Matrix orthotropicCauchyGreen(double E_1, double E_2, double E_3, double G_1, double G_2, double G_3,  double nu) ;
	static Vector rotate2ndOrderTensor2D( Vector & tensor, double angle ) ;
		
  
} ;

} 
#endif // TENSOR_H
