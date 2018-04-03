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

typedef enum
{
	SYMMETRY_CUBIC,
	SYMMETRY_HEXAGONAL,
	SYMMETRY_MONOCLINIC,
	SYMMETRY_ORTHORHOMBIC,
	SYMMETRY_TETRAGONAL,
	SYMMETRY_TRIGONAL,
	SYMMETRY_TRICLINIC,
} SymmetryType ;

typedef enum
{
	YOUNG_POISSON,
	BULK_SHEAR,
	YOUNG_SHEAR,
} IsotropicMaterialParameters ;

typedef enum : char {
  SINGLE_OFF_DIAGONAL_VALUES,
  DOUBLE_OFF_DIAGONAL_VALUES,  
} CompositionType ;
  
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
	std::vector<int> getIndex( size_t position ) const ;
	void swap( int i, int j ) ;

	size_t index( const std::vector<int> & indexes ) const ;
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

	static Tensor dotProduct( const Tensor &t1, const Tensor &t2, double tol = -1  ) ;
	static Tensor dotProduct( const Tensor &t1, const Tensor &t2, int i1, int i2, double tol = -1  ) ;
	static Tensor chainedDotProduct( const Tensor &tmi, const Tensor &tnj, const Tensor &tok, const Tensor & tpl, const Tensor &tijkl, double tol = -1 ) ;
	static Tensor chainedDotProduct( const Tensor &tmi, const Tensor &tnj, const Tensor &tij, double tol = -1 ) ;
	static Tensor dotdotProduct( const Tensor &t1, const Tensor &t2, double tol = -1  ) ;
	static Tensor crossProduct( const Tensor &t1, const Tensor &t2, double tol = -1  ) ;

	void print() const ;
	
	void threshold(double thr) ;

	static Matrix cauchyGreen( double p1, double p2, SpaceDimensionality = SPACE_TWO_DIMENSIONAL, planeType pt = PLANE_STRESS, IsotropicMaterialParameters param = YOUNG_POISSON ) ;
	static Matrix cauchyGreen( std::pair<double, double> props, SpaceDimensionality = SPACE_TWO_DIMENSIONAL, planeType pt = PLANE_STRESS, IsotropicMaterialParameters param = YOUNG_POISSON ) ;
	static std::pair<double, double> getIsotropicMaterialParameters( const Matrix & C, IsotropicMaterialParameters param, planeType pt ) ;
	static Matrix orthotropicCauchyGreen(double E_1, double E_2, double G,  double nu, planeType pt) ;
	static Matrix orthotropicCauchyGreen(double E_1, double E_2, double G,  double nu, double angle, planeType pt) ;
	static Matrix orthotropicCauchyGreen(double E_1, double E_2, double E_3, double G_1, double G_2, double G_3,  double nu) ;
	static Matrix orthotropicCauchyGreen( Vector data, SymmetryType sym, bool force = false ) ;
	static Matrix isotropicTransverseCauchyGreen(double E_1, double E_2, double G_12,  double nu_12, double nu_23, SpaceDimensionality dim, planeType pt) ;
	static Vector rotate2ndOrderTensor2D( const Vector & tensor, double angle ) ;
	static Vector rotate2ndOrderTensor3D( const Vector & tensor, Point angle, double tol = -1 ) ;
    
    static Matrix rotate2ndOrderTensor2D( const Matrix & tensor, double angle ) ;
    
	static Matrix rotate4thOrderTensor2D( const Matrix & tensor, double angle, double tol = -1 ) ;
	static Matrix rotate4thOrderTensor3D( const Matrix & tensor, Point angle, double tol = -1 ) ;
	static Matrix to2D( Matrix & tensor, planeType pt, Variable var = ZETA) ;
	static Matrix invert4thOrderTensor3D( const Matrix & tensor, SymmetryType sym = SYMMETRY_CUBIC ) ;
		
  
} ;

} 
#endif // TENSOR_H
