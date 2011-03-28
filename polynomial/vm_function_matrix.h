//
// C++ Interface: vm_function_matrix
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef VM_FUNCTION_MATRIX_H
#define VM_FUNCTION_MATRIX_H

#include <valarray>

#include "../utilities/sliceiters.h"
#include "../utilities/matrixops.h"
#include "vm_function_base.h"

namespace Mu
{


struct FMtFM ;
struct FMtM ;
struct Function ;
struct Gradient ;

typedef std::valarray<Function> FunctionVector ;

/** \brief Symbolic matrix with Functions as elements.
 * This class is used when symbolic matrix operations are needed.
 * It provides the usual matrix semantics: multiplication, addition and so on with .
*/
class FunctionMatrix 
{
	std::valarray< Function > *v;
	size_t r, c;
	
public:
	/** \brief Constructor.
	 * initalises the matrix with f(...) = 0
	 * @param x number of rows
	 * @param y number of columns
	 */
	FunctionMatrix(size_t x, size_t y);


	/** \brief Default constructor.
	 * creates a null two-by-two matrix.
	 */
	FunctionMatrix()
	{
		r = 2 ;
		c = 2 ;
		v = new std::valarray< Function >(4) ;
	}
	
	/** \brief Copy-constructor
	 * 
	 * @param  m source
	 */
	FunctionMatrix(const FunctionMatrix& m) ;
	
	virtual ~FunctionMatrix() {delete v; }
	

	/** \brief Assignement operator
	 * the matrix is reinitialised with the size and values of the source.
	 * @param m source
	 * @return self.
	 */
	FunctionMatrix &operator=(const FunctionMatrix &m);
	
	/** \brief return the number of elements in the matrix
	 * 
	 * @return number of elements
	 */
	size_t size() const {return r*c ;}


	/** \brief return number of columns.
	 * 
	 * @return number of columns
	 */
	size_t numCols() const {return c ;}


	/** \brief return number of rows.
	 * 
	 * @return number of rows
	 */
	size_t numRows() const {return r ;}
	

	/** \brief Return the ith column
	 * 
	 * @param i index of column to return
	 * @return valarray slice corresponding to the column
	 */
	Slice_iter< Function > column(size_t i )
	{
		return Slice_iter< Function >(v, std::slice(i, r, c)) ;
	}

	/** \brief Return the ith column
	 * 
	 * @param i index of column to return
	 * @return const valarray slice corresponding to the column
	 */
	Cslice_iter< Function > column(size_t i ) const 
	{
		return Cslice_iter< Function >(v, std::slice(i, r, c)) ;
	}
	
	/** \brief Return the ith row
	 * 
	 * @param i index of row to return
	 * @return valarray slice corresponding to the row
	 */
	Slice_iter< Function > row(size_t i )
	{
		return Slice_iter< Function >(v, std::slice(i*c, c, 1)) ;
	}

	/** \brief Return the ith row
	 * 
	 * @param i index of row to return
	 * @return valarray slice corresponding to the row
	 */
	Cslice_iter< Function > row(size_t i ) const
	{
		return Cslice_iter< Function >(v, std::slice(i*c, c, 1)) ;
	}
	
	/** \brief create a matrix transpose of the original
	 * 
	 * @return a new function matrix
	 */
	FunctionMatrix transpose() const ;
	
	/** \brief return the symbolic differential according to a variable.
	 * 
   * This methos only functions if the symbolic differential of the members of the matrix are defined, 
   * otherwise a null matrix will be returned
	 * @param v variable with respect to which the differenciation should be done
	 * @return the symbolic differential of this matrix.
	 */
	FunctionMatrix d(const Variable v) const ;
	
	/** \brief Accessor of the (x, y) element of the matrix
	 * 
	 * @param x row
	 * @param y column
	 * @return a reference to the element pointed to.
	 */
	Function& operator()(size_t x, size_t y) ;

	/** \brief Accessor of the (x, y) element of the matrix
	 * 
	 * @param x row
	 * @param y column
	 * @return a reference to the element pointed to.
	 */
	Function operator()(size_t x, size_t y) const;
	
	/** \brief Return the ith row
	 * 
	 * @param i index of row to return
	 * @return valarray slice corresponding to the row
	 */
	Slice_iter< Function > operator()(size_t i) {return row(i) ;}

	/** \brief Return the ith row
	 * 
	 * @param i index of row to return
	 * @return valarray slice corresponding to the row
	 */
	Cslice_iter< Function > operator()(size_t i) const {return row(i) ;}
	
	/** \brief Return the ith row
	 * 
	 * @param i index of row to return
	 * @return valarray slice corresponding to the row
	 */
	Slice_iter< Function > operator[](size_t i) {return row(i) ;}

	/** \brief Return the ith row
	 * 
	 * @param i index of row to return
	 * @return valarray slice corresponding to the row
	 */
	Cslice_iter< Function > operator[](size_t i) const {return row(i) ;}
	

	/** \brief Assign from the product of two FunctionMatrix s
	 * 
	 * @param m operation to perform
	 */
	FunctionMatrix& operator =(const FMtFM& m) ;


	/** \brief Assign from the product of a FunctionMatrix and a Matrix
	 * 
	 * @param m operation to perform
	 */
	FunctionMatrix& operator =(const FMtM& m) ;
	

	/** \brief Multiply by a function
	 * 
	 * @param f Function with which to multiply
	 * @return a new FunctionMatrix
	 */
	FunctionMatrix operator *(const Function & f) const;


	/** \brief Divide by a function
	 * 
	 * @param  f Function with which to divide
	 * @return a new FunctionMatrix
	 */
	FunctionMatrix operator /(const Function & f) const;

	/** \brief Divide by a function and assign
	 * 
	 * @param f Function with which to divide
	 * @return self
	 */
	FunctionMatrix& operator /=(const Function & f) ;

	/** \brief Multiply by a function and assign
	 * 
	 * @param m Function with which to multiply
	 * @return self
	 */
	FunctionMatrix& operator *=(const Function &m);

	/** \brief Multiply by a FunctionMatrix and assign
	 * 
   * There is no check on the compatibility of the respective dimension of both matrices, so the operation can fail
	 * @param m FunctionMatrix with which to multiply
	 * @return self
	 */
	FunctionMatrix& operator *=(const FunctionMatrix &m);


	/** \brief Add a FunctionMatrix and assign
	 * 
   * There is no check on the compatibility of the respective dimension of both matrices, so the operation can fail
	 * @param m FunctionMatrix to add
	 * @return self
	 */
	FunctionMatrix& operator +=(const FunctionMatrix &m) ;

	/** \brief Add a FunctionMatrix
	 * 
   * There is no check on the compatibility of the respective dimension of both matrices, so the operation can fail
	 * @param m FunctionMatrix to add
	 * @return a new FunctionMatrix
	 */
	FunctionMatrix operator +(const FunctionMatrix &m) const;

	/** \brief Substract a FunctionMatrix and assign
	 * 
   * There is no check on the compatibility of the respective dimension of both matrices, so the operation can fail
	 * @param m FunctionMatrix to substract
	 * @return self
	 */
	FunctionMatrix& operator -=(const FunctionMatrix &m) ;


	/** \brief Substract a FunctionMatrix
	 * 
   * There is no check on the compatibility of the respective dimension of both matrices, so the operation can fail
	 * @param m FunctionMatrix to substract
	 * @return A new FunctionMatrix
	 */
	FunctionMatrix operator -(const FunctionMatrix &m) const;
	

	/** \brief Multiply by a double and assign
	 * 
	 * @param a scalar multiplicator
	 * @return self
	 */
	FunctionMatrix& operator *=(const double a);

	/** \brief Multiply by a double
	 * 
	 * @param a scalar multiplicator
	 * @return A new FunctionMatrix
	 */
	FunctionMatrix operator *(const double a) const;

	/** \brief Divide by a double
	 * 
	 * @param a scalar divisor
	 * @return A new FunctionMatrix
	 */
	FunctionMatrix operator /(const double a) const;

	/** \brief Divide by a double and assign
	 * 
	 * @param a scalar divisor
	 * @return self
	 */
	FunctionMatrix& operator /=(const double a) ;

	/** \brief Add a Matrix and assign
	 * 
	 * @param m matrix to add
	 * @return self
	 */
	FunctionMatrix& operator +=(const Matrix &m) ;

	/** \brief Add a Matrix
	 * 
	 * @param m matrix to add
	 * @return A new FunctionMatrix
	 */
	FunctionMatrix operator +(const Matrix &m) const;


	/** \brief Substract a Matrix and assign
	 * 
	 * @param m matrix to substract
	 * @return self
	 */
	FunctionMatrix& operator -=(const Matrix &m) ;

	/** \brief Substract a Matrix
	 * 
	 * @param m matrix to substract
	 * @return A new FunctionMatrix
	 */
	FunctionMatrix operator -(const Matrix &m) const;
	
	/** \brief Accessor to the array of Function s used to store the elements of the FunctionMatrix
	 */
	std::valarray< Function > &array() {return *v ;}

	/** \brief Accessor to the array of Function s used to store the elements of the FunctionMatrix
	 */
	std::valarray< Function > array() const {return *v ;}
} ;

/** \brief Structure to contain the references necessary for the lazy computation of a symbolic matrix-vector multiplication
*/
struct FMtFV
{
	const FunctionMatrix &m;
	const FunctionVector &v;
	
	/** \brief Constructor. Initialises the references
	 * 
	 * @param mm FunctionMatrix reference
	 * @param vv FunctionVector reference
	 */
	FMtFV(const FunctionMatrix &mm, const FunctionVector &vv) : m(mm), v(vv) { }
	
	/** \brief Compute the result of the matrix-vector multiplication
	 * 
	 */
	operator const FunctionVector();
} ;

/** \brief Structure to contain the references necessary for the lazy computation of a symbolic matrix-vector multiplication
*/
struct FMtV
{
	const FunctionMatrix &m;
	const Vector &v;

	/** \brief Constructor. Initialises the references
	 * 
	 * @param mm FunctionMatrix reference
	 * @param mmm Vector reference
	 */
	FMtV(const FunctionMatrix &mm, const Vector &vv) : m(mm), v(vv) { }

	/** \brief Compute the result of the matrix-matrix multiplication.
	 * 
	 */
	operator const FunctionVector();
} ;

/** \brief Structure to contain the references necessary for the lazy computation of a symbolic matrix-matrix multiplication
*/
struct FMtFM
{
	const FunctionMatrix &first;
	const FunctionMatrix &second;
	
	/** \brief Constructor. Initialises the references
	 * 
	 * @param mm FunctionMatrix reference
	 * @param mmm FunctionMatrix reference
	 */
	FMtFM(const FunctionMatrix &mm, const FunctionMatrix &mmm) : first(mm), second(mmm) { }
	
	/** \brief Compute the result of the matrix-matrix multiplication.
	 * 
	 */
	operator const FunctionMatrix() const;
} ;

/** \brief Structure to contain the references necessary for the lazy computation of a  matrix - symbolic matrix multiplication
*/
struct MtFM
{
	const Matrix &first;
	const FunctionMatrix &second;
	
	/** \brief Constructor. Initialises the references
	 * 
	 * @param mm Matrix reference
	 * @param mmm FunctionMatrix reference
	 */
	MtFM(const Matrix &mm, const FunctionMatrix &mmm) : first(mm), second(mmm) { }
	
	/** \brief Compute the result of the matrix-matrix multiplication.
	 * 
	 */
	operator const FunctionMatrix() const;
} ;

/** \brief Structure to contain the references necessary for the lazy computation of a  symbolic matrix - matrix multiplication
*/
struct FMtM
{
	const FunctionMatrix &first;
	const Matrix &second;
	/** \brief Constructor. Initialises the references
	 * 
	 * @param mm  FunctionMatrix reference
	 * @param mmm  Matrix reference
	 */
	FMtM(const FunctionMatrix &mm, const Matrix &mmm) : first(mm), second(mmm) { }
	
	/** \brief Compute the result of the matrix-matrix multiplication.
	 * 
	 */
	operator const FunctionMatrix() const;
} ;

/** \brief Structure to contain the references necessary for the lazy computation of a  symbolic matrix - matrix - symbolic matrix multiplication
*/
struct FMtMtFM
{
	const FunctionMatrix &first;
	const Matrix &second;
	const FunctionMatrix &third;
	
	/** \brief Constructor. Initialises the references
	 * 
	 * @param mm  FunctionMatrix reference
	 * @param mmm  Matrix reference
	 * @param mmmm  FunctionMatrix reference
	 */
	FMtMtFM(const FunctionMatrix &mm, const Matrix &mmm, const FunctionMatrix &mmmm) : first(mm), second(mmm), third(mmmm) { }
	
	/** \brief Compute the result of the matrix-matrix-matrix multiplication.
	 * 
	 */
	operator const FunctionMatrix() const;
} ;

/** \brief Structure to contain the references necessary for the lazy computation of a  Gradient - symbolic matrix - Gradient multiplication
*/
struct GtFMtG
{
	const Gradient &first ;
	const FunctionMatrix &second ;
	const Gradient &third ;
	
	/** \brief Constructor. Initialises the references
	 * 
	 * @param g  Gradient reference
	 * @param f  FunctionMatrix reference
	 * @param g_  Gradient reference
	 */
	GtFMtG(const Gradient & g, const FunctionMatrix & f,const Gradient & g_) : first(g), second(f), third(g_) { };
	
} ;

/** \brief Structure to contain the references necessary for the lazy computation of a  Gradient - symbolic matrix multiplication
*/
struct GtFM
{
	const Gradient &first ;
	const FunctionMatrix &second ;
	
	/** \brief Constructor. Initialises the references
	 * 
	 * @param g Gradient
	 * @param f FunctionMatrix
	 */
	GtFM(const Gradient & g, const FunctionMatrix & f) : first(g), second(f) { };

	/** \brief Operator to create a Gradient - symbolic matrix - Gradient multiplication
	 * 
	 * @param f Gradient
	 */
	GtFMtG operator*(const Mu::Gradient & f) const ;
} ;


} ;


/** \brief Create a symbolic matrix - symbolic vector multiplication object of lazy evaluation
 * 
 * @param mm FunctionMatrix
 * @param v FunctionVector
 * @return FMtFV
 */
Mu::FMtFV operator*(const Mu::FunctionMatrix& mm, const Mu::FunctionVector& v);


/** \brief Create a symbolic matrix - symbolic matrix multiplication object of lazy evaluation
 * 
 * @param mm FunctionMatrix
 * @param mmm FunctionMatrix
 * @return FMtFM
 */
Mu::FMtFM operator*(const Mu::FunctionMatrix& mm, const Mu::FunctionMatrix& mmm);

/** \brief Create a symbolic matrix - matrix - symbolic matrix multiplication object of lazy evaluation
 * 
 * @param mm FMtM
 * @param mmm FunctionMatrix
 * @return FMtMtFM
 */
Mu::FMtMtFM operator*(const Mu::FMtM& mm, const Mu::FunctionMatrix& mmm);

/** \brief Create a symbolic matrix - vector multiplication object of lazy evaluation
 * 
 * @param mm FunctionMatrix
 * @param v Vector
 * @return FMtV
 */
Mu::FMtV operator*(const Mu::FunctionMatrix& mm, const Vector & v);

/** \brief Create a  symbolic matrix - matrix multiplication object of lazy evaluation
 * 
 * @param mm FunctionMatrix
 * @param mmm Matrix
 * @return FMtM
 */
Mu::FMtM operator*(const Mu::FunctionMatrix& mm, const Mu::Matrix& mmm);

/** \brief Create a  matrix - symbolic matrix multiplication object of lazy evaluation
 * 
 * @param mm Matrix
 * @param mmm FunctionMatrix
 * @return MtFM
 */
Mu::MtFM operator*(const Mu::Matrix& mm, const Mu::FunctionMatrix& mmm);

/** \brief Perform symbolic matrix - symbolic matrix multiplication
 * 
 * The respective dimensions of the matrices are not checked
 * @param m0 FunctionMatrix
 * @param m1 FunctionMatrix
 * @return FunctionMatrix
 */
const Mu::FunctionMatrix ff_matrix_multiply(const Mu::FunctionMatrix &m0, const Mu::FunctionMatrix &m1 );

/** \brief Perform symbolic matrix - matrix multiplication
 * 
 * The respective dimensions of the matrices are not checked
 * @param m0 FunctionMatrix
 * @param m1 Matrix
 * @return FunctionMatrix
 */
const Mu::FunctionMatrix fm_matrix_multiply(const Mu::FunctionMatrix &m0, const Mu::Matrix &m1 );

/** \brief Perform matrix - symbolic matrix multiplication
 * 
 * The respective dimensions of the matrices are not checked
 * @param m0 Matrix
 * @param m1 FunctionMatrix
 * @return FunctionMatrix
 */
const Mu::FunctionMatrix mf_matrix_multiply(const Mu::Matrix &m0, const Mu::FunctionMatrix &m1 );

/** \brief Perform matrix - symbolic vector multiplication
 * 
 * The respective dimensions of the objects are not checked
 * @param m Matrix
 * @param v FunctionVector
 * @return FunctionMatrix
 */
const Mu::FunctionVector matrix_fvector_multiply(const Mu::Matrix &m, const Mu::FunctionVector &v );

/** \brief Perform symbolic matrix - vector multiplication
 * The respective dimensions of the objects are not checked
 * @param m FunctionMatrix
 * @param v Vector
 * @return FunctionVector
 */
const Mu::FunctionVector fmatrix_vector_multiply(const Mu::FunctionMatrix &m, const Vector &v );

/** \brief Perform symbolic vector - symbolic matrix multiplication
 * The respective dimensions of the objects are not checked
 * @param v FunctionVector
 * @param m FunctionMatrix
 * @return FunctionVector
 */
const Mu::FunctionVector operator*(const Mu::FunctionVector &v , const Mu::FunctionMatrix &m );

/** \brief Perform symbolic vector - symbolic matrix multiplication
 * The respective dimensions of the objects are not checked
 * @param v Vector
 * @param m FunctionMatrix
 * @return FunctionVector
 */
const Mu::FunctionVector operator*(const Vector &v , const Mu::FunctionMatrix &m );

/** \brief Compute the symbolic inverse of a 2x2 symbolic matrix
 * 
 * @param s FunctionMatrix to invert
 * @return FunctionMatrix
 */
Mu::FunctionMatrix inverse2x2FunctionMatrix(const Mu::FunctionMatrix s) ;

/** \brief Compute the symbolic inverse of a 3x3 symbolic matrix
 * 
 * @param m FunctionMatrix to invert
 * @return FunctionMatrix
 */
Mu::FunctionMatrix inverse3x3FunctionMatrix(const Mu::FunctionMatrix m) ;

/** \brief Create a  Gradient - symbolic matrix multiplication object of lazy evaluation
 * 
 * @param g Gradient
 * @param m FunctionMatrix
 * @return GtFM
 */
Mu::GtFM operator *(const Mu::Gradient & g, const Mu::FunctionMatrix & m) ;



#endif

