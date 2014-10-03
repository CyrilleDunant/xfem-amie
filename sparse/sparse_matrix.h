//
// C++ Interface: sparse_matrix
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __SPARSE_MATRIX_H
#define __SPARSE_MATRIX_H

#include<valarray>
#include<map>
#include <set>
#include <complex>
#include "../utilities/matrixops.h"
#include "sparse_vector.h"

namespace Amie {

struct ConstSymetricSparseVector ;
struct CompositeSparseMatrixTimesVec ;
struct CoordinateIndexedSparseMatrixTimesVec ;
struct CoordinateIndexedSparseMatrixTimesVecPlusVec ;
struct CoordinateIndexedSparseMatrixTimesVecMinusVec ;

/** \brief Coordinate-indexed sparse matrix. 
 *
 * This sparse matrix is the most general possible: it makes no assumptions on shape or structure. Every element is indexed 
 * by its row and column indices. This type of matrix is designed ot be used in conjunction with the more numerically
 * efficient CoordinateIndexedSparseMatrix, to provided offsets in non-linear solver loops for example.
 * 
 */
class CoordinateIndexedIncompleteSparseMatrix
{
public:
	Vector  array ;
	std::valarray<unsigned int>  column_index ;
	std::valarray<unsigned int>  row_index ;

	size_t stride ;
public:
	
/** \brief Constructor using two arrays as an input. 
The arrays define the sparseness pattern of the matrix.
 * @param rows List of row indices
 * @param cols List of corresponding Column indices
 */
	CoordinateIndexedIncompleteSparseMatrix(const std::valarray<unsigned int> & rows, const std::valarray<unsigned int> &cols ) ;
	
/** \brief Form a symmetric sparse matrix, given a line of index i, and column indices and corresponding values.
 * 
 * This is useful for the formation of Lagrange multipliers.
 * @param rows List of row indices
 * @param cols List of corresponding Column indices
 */
	CoordinateIndexedIncompleteSparseMatrix(const std::valarray<unsigned int> &, const Vector &, const size_t i) ;
	~CoordinateIndexedIncompleteSparseMatrix() ;
	
	/** \brief Form and return the sparse vector corresponding to line i. The element can be modified if it is included in the sparseness pattern.*/
	SparseVector operator[](const size_t i) ;
	
	/** \brief Form and return the sparse vector corresponding to line i*/
	ConstSparseVector operator[](const size_t i) const ;
	
	/** \brief Return the element located at the given row and column. return 0 otherwise*/
	double & operator()(const size_t i, const size_t j) ;
	
	/** \brief Compute a matrix-vector multiplication. This sparse matrix is not intended for efficient linear algebra.*/
	Vector operator *(const Vector v) const ;
	
} ;

/** \brief Block-encoded sparse matrix. 
 *
 * The matrix is assumed to have more lines than rows. The stride gives the size of the blocks.
 * The matrix is stored in row-major format, meaning groups of stride lines follow each order in the memory. However, those groups
 * are in turn stored in blocks which are column major.
 * Column major order is best for Matrix-vector multiplication, though it is not the "natural" order of storage.
 * Blocks are arranged so values are 16-bit aligned, so SSE operations are possible.
 * class CoordinateIndexedSparseMatrix 
 *
 * This Matrix is the most commonly used in the solvers. Helper structures and functions allow easy expression of 
 * linear algebra operations:
 \code 
Vector result(2) ;
Vector b(2) ; b[0] = 2 ; b[1] = 3 ;
 //CoordinateIndexedSparseMatrix M  has been previously initalised

assign(result, M*b-b) ;
\endcode
*/
class CoordinateIndexedSparseMatrix
{
public:
	size_t stride ;
	Vector  array ;
	std::valarray<unsigned int>  column_index ;
	std::valarray<unsigned int>  row_size ;
	std::valarray<unsigned int>  accumulated_row_size ;
	
public:

	void reshape(std::set<std::pair<size_t, size_t>> &source, size_t stride) ;

/** \brief Constructor using a pair of valarays as input for the parseness pattern. The matrix is initialised at 0.
 * @param rowSize A valarray of unsigned ints containing the succesive length of each row
 * @param columnIndex the column indes of every value to be sotred in the matrix.
 */
	CoordinateIndexedSparseMatrix(const std::valarray<unsigned int> & rowSize, const std::valarray<unsigned int> & columnIndex, size_t stride) ;
    
    CoordinateIndexedSparseMatrix(const CoordinateIndexedSparseMatrix & src) : stride(src.stride), array(src.array), column_index(src.column_index), row_size(src.row_size), accumulated_row_size(src.accumulated_row_size) { } ;
	
	~CoordinateIndexedSparseMatrix() ;
	
	/** \brief Return the sparse row of index i. 
   *
   * All rows within the block row are returned, however this has an effect only of the returned sparse 
	 * vector is multiplied with another sparse vector. The vector is writable.*/
	SparseVector operator[](const size_t i) ;

	/** \brief  Return the sparse row of index i. 
   *
   * All rows within the block row are returned, however this has an effect only of the returned sparse vector is multiplied 
	 * with another sparse vector*/
	const ConstSparseVector operator[](const size_t i) const;
	
	/** \brief Return a reference to the element at location i,j if it exists, 0 otherwise. The matrix can be updated in this way.*/
	double & operator()(const size_t i, const size_t j) ;
	
	/** \brief Return a reference to the element at location i,j if it exists, 0 otherwise. */
	double  operator()(const size_t i, const size_t j) const ;
		
/** \brief Compute a Matrix-vector operation.*/
	CoordinateIndexedSparseMatrixTimesVec operator *(const Vector & v) const ;
// 	CoordinateIndexedSparseMatrix operator +(const CoordinateIndexedIncompleteSparseMatrix &) const ;

/** \brief Add to matrix the values stored in an incomplete sparse matrix.*/
	CoordinateIndexedSparseMatrix & operator +=(const CoordinateIndexedIncompleteSparseMatrix &) ;
// 	CoordinateIndexedSparseMatrix operator +(const CoordinateIndexedSparseMatrix &) const ;

/** \brief reinitialise matrix from a new source. all values are lost.*/
	CoordinateIndexedSparseMatrix & operator=(const CoordinateIndexedSparseMatrix &) ;
	
/** \brief Return the inverse of the diagonal. 

As there are no warranty all members of the diagonal are non-0, there is no warranty that all values in the return vector be defined. This function is useful for the initialisation of preconditionners.*/
	Vector inverseDiagonal() const ;

/** \brief Return the inverse of the diagonal, squared. 

As there are no warranty all members of the diagonal are non-0, there is no warranty that all values in the return vector be defined. */
	Vector inverseDiagonalSquared() const ;

/** \brief Return the diagonal. */
	Vector diagonal() const ;

/** \brief Construct a new coordinate indexed sparse matrix, transpose of this one.
*/
	CoordinateIndexedSparseMatrix transpose() const ;
	
/** \brief  Compute the froebenius norm of the matrix: \f$ \sqrt{\sum A_{ij}A_{ij} \forall i,j}\f$*/
	double froebeniusNorm() const ;

/** \brief Compute the infinity norm of the matrix: \f$ |\max A_{ij}|\f$*/
	double infinityNorm() const ;
	
	inline void inner_product(const double*  v, double*  dest, const size_t rowstart, const size_t colstart, const size_t index) const
	{
		size_t length = row_size[index/stride] ;
		size_t start = accumulated_row_size[index/stride] ;
		
		size_t mstart = start ;
		if(colstart)
		{
			const unsigned int * __start__       = &column_index[start] ;
			const unsigned int * __end__         = &column_index[start+length] ;
			const unsigned int * i_index_pointer = std::lower_bound(__start__, __end__, colstart/stride) ;
			const unsigned int offset            = i_index_pointer - __start__ ;
			
			mstart = (start <= start+offset && start+offset < start+length) ? (start+offset) : start+length;
		}
	// 		mstart = start ;std::max(mstart,start) ;
	// 		
	// 		
		switch(stride)
		{
			case 1:
			{
				for(unsigned int j = mstart ; j < length+start ; j++)
				{
					*dest += v[column_index[j]]*array[j*2] ;
				}
				return ;
			}
			case 2:
			{
				const double * array_iterator = &array[mstart*2*2] ;
				const double * vec_iterator = v + column_index[mstart]*2 ;
				for(unsigned int j = mstart ; j < length+start ; j++ )
				{
					*dest     += *array_iterator*(*vec_iterator);
					*(dest+1) += *(array_iterator+1)*(*vec_iterator) ;
					*dest     += *(array_iterator+2)*(*(vec_iterator +1)) ;
					*(dest+1) += *(array_iterator+3)*(*(vec_iterator +1)) ;
					array_iterator+=4 ;
					if(j+1 < column_index.size())
						vec_iterator += column_index[j+1]*stride-column_index[j]*stride ;
				}
				return ;
			}
			
			case 3:
			{
				const double * array_iterator = &array[mstart*3*4] ;
				const double * vec_iterator = v + column_index[mstart]*3 ;
				for(unsigned int j = mstart ; j < length+start ;j++)
				{
					*dest     += *array_iterator     * (*vec_iterator);
					*(dest+1) += *(array_iterator+1) * (*vec_iterator);
					*(dest+2) += *(array_iterator+2) * (*vec_iterator);
					*dest     += *(array_iterator+4) * (*(vec_iterator + 1)) ;
					*(dest+1) += *(array_iterator+5) * (*(vec_iterator + 1)) ;
					*(dest+2) += *(array_iterator+6) * (*(vec_iterator + 1)) ;
					*dest     += *(array_iterator+8 ) * (*(vec_iterator + 2));
					*(dest+1) += *(array_iterator+9 ) * (*(vec_iterator + 2));
					*(dest+2) += *(array_iterator+10) * (*(vec_iterator + 2));
					array_iterator+=12 ;
					if(j+1 < column_index.size())
						vec_iterator += column_index[j+1]*stride-column_index[j]*stride ;
				}
				return ;
			}
			case 6 :
			{
				const int colLength = 6  ;
				#ifdef HAVE_SSE3
				const __m128d * array_iterator = (__m128d*)&array[mstart*36] ;
	// 			#pragma omp parallel for schedule(runtime)
				for(unsigned int j = mstart ; j < length+start ; j++)
				{
					__m128d vval =  _mm_set1_pd(*(v+column_index[j]*stride)) ;

					_mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*array_iterator, vval))) ;
					_mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+1), vval))) ;
					_mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+2), vval))) ;
					
					vval =  _mm_set1_pd(*(v+column_index[j]*stride+1)) ;

					_mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+3,) vval))) ;
					_mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+4), vval))) ;
					_mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+5), vval))) ;
					
					vval =  _mm_set1_pd(*(v+column_index[j]*stride+2)) ;

					_mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+6), vval))) ;
					_mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+7), vval))) ;
					_mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+8), vval))) ;
					
					vval =  _mm_set1_pd(*(v+column_index[j]*stride+3)) ;

					_mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+9), vval))) ;
					_mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+10), vval))) ;
					_mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+11), vval))) ;
					
					vval =  _mm_set1_pd(*(v+column_index[j]*stride+4)) ;

					_mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+12), vval))) ;
					_mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+13), vval))) ;
					_mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+14), vval))) ;
					
					vval =  _mm_set1_pd(*(v+column_index[j]*stride+5)) ;

					_mm_store_pd((dest),  _mm_add_pd( _mm_load_pd((dest)), _mm_mul_pd(*(array_iterator+15), vval))) ;
					_mm_store_pd((dest+2),  _mm_add_pd( _mm_load_pd((dest+2)), _mm_mul_pd(*(array_iterator+16), vval))) ;
					_mm_store_pd((dest+4),  _mm_add_pd( _mm_load_pd((dest+4)), _mm_mul_pd(*(array_iterator+17), vval))) ;
					array_iterator +=18 ;
				}
				#else
				const double * array_iterator0 = &array[mstart*36] ;
				const double * array_iterator1 = &array[mstart*36+1] ;
	// 			#pragma omp parallel for schedule(runtime)
				for(unsigned int j = mstart ; j < length+start ; j++)
				{
					double vval =  *(v+column_index[j]*6) ;

					*(dest) += *array_iterator0 * vval ;
					*(dest+1) += *array_iterator1 * vval ;
					
					*(dest+2) += *(array_iterator0+2) * vval ;
					*(dest+3) += *(array_iterator1+2) * vval ;
					
					*(dest+4) += *(array_iterator0+4) * vval ;
					*(dest+5) += *(array_iterator1+4) * vval ;
					
					vval =  *(v+column_index[j]*6+1) ;

					*(dest) += *(array_iterator0+6) * vval ;
					*(dest+1) += *(array_iterator1+6) * vval ;
					
					*(dest+2) += *(array_iterator0+8) * vval ;
					*(dest+3) += *(array_iterator1+8) * vval ;
					
					*(dest+4) += *(array_iterator0+10) * vval ;
					*(dest+5) += *(array_iterator1+10) * vval ;
					
					vval =  *(v+column_index[j]*6+2) ;

					*(dest) += *(array_iterator0+12) * vval ;
					*(dest+1) += *(array_iterator1+12) * vval ;
					
					*(dest+2) += *(array_iterator0+14) * vval ;
					*(dest+3) += *(array_iterator1+14) * vval ;
					
					*(dest+4) += *(array_iterator0+16) * vval ;
					*(dest+5) += *(array_iterator1+16) * vval ;
					
					vval =  *(v+column_index[j]*6+3) ;

					*(dest) += *(array_iterator0+18) * vval ;
					*(dest+1) += *(array_iterator1+18) * vval ;
					
					*(dest+2) += *(array_iterator0+20) * vval ;
					*(dest+3) += *(array_iterator1+20) * vval ;
					
					*(dest+4) += *(array_iterator0+22) * vval ;
					*(dest+5) += *(array_iterator1+22) * vval ;
					
					vval = *(v+column_index[j]*6+4) ;

					*(dest) += *(array_iterator0+24) * vval ;
					*(dest+1) += *(array_iterator1+24) * vval ;
					
					*(dest+2) += *(array_iterator0+26) * vval ;
					*(dest+3) += *(array_iterator1+26) * vval ;
					
					*(dest+4) += *(array_iterator0+28) * vval ;
					*(dest+5) += *(array_iterator1+28) * vval ;
					
					vval =  *(v+column_index[j]*6+5) ;

					*(dest) += *(array_iterator0+30) * vval ;
					*(dest+1) += *(array_iterator1+30) * vval ;
					
					*(dest+2) += *(array_iterator0+32) * vval ;
					*(dest+3) += *(array_iterator1+32) * vval ;
					
					*(dest+4) += *(array_iterator0+34) * vval ;
					*(dest+5) += *(array_iterator1+34) * vval ;
					array_iterator0 += 36 ;
					array_iterator1 += 36 ;
				}
				#endif
				return ;
			}
			default:
			{
				unsigned int colLength = stride+stride%2 ;
				const double * array_iterator0 = &array[mstart*colLength*stride] ;
				const double * array_iterator1 = &array[mstart*colLength*stride+1] ;
				for(unsigned int j = mstart ; j < length+start ; j++)
				{
					for(size_t c = 0 ; c < stride ; c++)
					{
						const double vval =  *(v+column_index[j]*stride+c) ;
						for(size_t i = 0 ; i != colLength ; i+=2)
						{
							*(dest+i) += *array_iterator0 * vval ;
							*(dest+i+1) += *array_iterator1 * vval ;
							array_iterator0 += 2 ;
							array_iterator1 += 2 ;
						}
					}
				}
				return ;
			}
		}
	} ;
} ;

struct CompositeSparseMatrixTimesVec ;
struct CompositeSparseMatrixTimesVecPlusVec ;
struct CompositeSparseMatrixTimesVecMinusVec ;
struct CompositeSparseMatrixTimesVecMinusVecMinusVec ;

/** \brief Temporary structure used for lazy evaluation of sparse matrix sums.

his structure contains a reference to a CoordinateIndexedSparseMatrix and a CoordinateIndexedIncompleteSparseMatrix, considered to be summed.
 */
struct CompositeSparseMatrix
{
	const CoordinateIndexedSparseMatrix  & sm ;
	const CoordinateIndexedIncompleteSparseMatrix  & ism ;
	/** \brief Constructor.
	Initialises the references.
	@param sm_ CoordinateIndexedSparseMatrix to refer.
	@param ism_ CoordinateIndexedIncompleteSparseMatrix to refer.
*/
	CompositeSparseMatrix(const CoordinateIndexedSparseMatrix  & sm_ , const CoordinateIndexedIncompleteSparseMatrix  & ism_)  : sm(sm_), ism(ism_) { };
	
	/** \brief overloaded operator to multiply with a Vector
	 The result of the operation is a temporary structure used for lazy evaluation.
	 */
	CompositeSparseMatrixTimesVec operator *(const Vector v) const ;

	/** \brief Produce the CoordinateIndexedSparseMatrix corresponding to this structure.*/
	operator CoordinateIndexedSparseMatrix() const ;
} ;

/** \brief Structure contaning the elements for a (matrix + matrix) * vector operation 
 This structure allows lazy evaluation of (matrix + matrix) * vector, minimizing the use of temporaries, for better numeric performance.
 */
struct CompositeSparseMatrixTimesVec
{
	const CompositeSparseMatrix  & sm ;
	const Vector & ve ;
		/** \brief Constructor.
	Initialises the references.
	@param c CompositeSparseMatrix to refer.
	@param v Vector to refer.
*/
	CompositeSparseMatrixTimesVec(const CompositeSparseMatrix  &s,const Vector & v ) : sm(s), ve(v) {} ;
	
	/** \brief overloaded operator to add with a Vector
	 The result of the operation is a temporary structure used for lazy evaluation.
	@param v Vector to add.
	 */
	CompositeSparseMatrixTimesVecPlusVec operator + (const Vector & v) const ;
	
		/** \brief overloaded operator to substract with a Vector
	 The result of the operation is a temporary structure used for lazy evaluation.
		@param v vector to substract.
	 */
	CompositeSparseMatrixTimesVecMinusVec operator - (const Vector & v) const ;

	/** \brief Produce the Vector corresponding to this operation represented by this structure.*/
	operator const Vector() const ;
} ;

/** \brief Structure contaning the elements for a matrix * vector operation 
 This structure allows lazy evaluation of matrix * vector, minimizing the use of temporaries, for better numeric performance.
 */
struct CoordinateIndexedSparseMatrixTimesVec
{
	const CoordinateIndexedSparseMatrix  & sm ;
	const Vector & ve ;
/** \brief Constructor.
	Initialises the references.
	@param s CoordinateIndexedSparseMatrix to refer.
	@param v Vector to refer.
*/
	CoordinateIndexedSparseMatrixTimesVec(const CoordinateIndexedSparseMatrix  &s,const Vector & v ) : sm(s), ve(v) {} ;
	
/** \brief overloaded operator to add with a Vector
	 The result of the operation is a temporary structure used for lazy evaluation.
	@param v Vector to add.
	*/
	CoordinateIndexedSparseMatrixTimesVecPlusVec operator + (const Vector & v) const ;
	
/** \brief overloaded operator to substract with a Vector
	 The result of the operation is a temporary structure used for lazy evaluation.
	 @param v vector to substract.
	*/
	CoordinateIndexedSparseMatrixTimesVecMinusVec operator - (const Vector & v) const ;

	/** \brief Produce the Vector corresponding to this operation represented by this structure.*/
	operator const Vector() const ;
} ;

/** \brief Structure contaning the elements for a (matrix+matrix) * vector + vector operation 
 This structure allows lazy evaluation of (matrix+matrix) * vector, minimizing the use of temporaries, for better numeric performance.
 */
struct CompositeSparseMatrixTimesVecPlusVec
{
	const CompositeSparseMatrixTimesVec & co ;
	const Vector & ve ;
	/** \brief Constructor.
	Initialises the references.
	@param c CompositeSparseMatrixTimesVec to refer.
	@param v Vector to refer.
*/
	CompositeSparseMatrixTimesVecPlusVec(const CompositeSparseMatrixTimesVec & c, const Vector & v) : co(c), ve(v) {} ;
	
	/** \brief Produce the Vector corresponding to this operation represented by this structure.*/
	operator const Vector() const 
	{
		Vector ret(0., ve.size()) ;
		
		for (size_t i = 0 ; i < co.sm.sm.row_size.size(); i++)
		{
			Vector temp = co.sm.sm[i*co.sm.sm.stride]*ve ;
			for(size_t j = i*co.sm.sm.stride ; j < i*co.sm.sm.stride+co.sm.sm.stride ; j++)
				ret[j] = temp[j-i*co.sm.sm.stride] +ve[j];
		}
		return ret ;
	}
	
} ;

/** \brief Structure contaning the elements for a matrix * vector + vector operation 
 This structure allows lazy evaluation of matrix * vector + vector, minimizing the use of temporaries, for better numeric performance.
 */
struct CoordinateIndexedSparseMatrixTimesVecPlusVec
{
	const CoordinateIndexedSparseMatrixTimesVec & co ;
	const Vector & ve ;
	/** \brief Constructor.
	Initialises the references.
	@param c CoordinateIndexedSparseMatrixTimesVec to refer.
	@param v Vector to refer.
*/
	CoordinateIndexedSparseMatrixTimesVecPlusVec(const CoordinateIndexedSparseMatrixTimesVec & c, const Vector & v) : co(c), ve(v) {} ;
	
	/** \brief Produce the Vector corresponding to this operation represented by this structure.*/
	operator const Vector() const 
	{
		Vector ret(0., ve.size()) ;
		
		for (size_t i = 0 ; i < co.sm.row_size.size(); i++)
		{
			Vector temp = co.sm[i*co.sm.stride]*co.ve ;
			for(size_t j = i*co.sm.stride ; j < i*co.sm.stride+co.sm.stride ; j++)
				ret[j] = temp[j-i*co.sm.stride] +ve[j];
		}
		return ret ;
	}
	
} ;

/** \brief Structure contaning the elements for a matrix * vector - vector operation 
 This structure allows lazy evaluation of matrix * vector-vector, minimizing the use of temporaries, for better numeric performance.
 */
struct CoordinateIndexedSparseMatrixTimesVecMinusVec
{
	const CoordinateIndexedSparseMatrixTimesVec & co ;
	const Vector & ve ;
/** \brief Constructor.
	Initialises the references.
	@param c CoordinateIndexedSparseMatrixTimesVec to refer.
	@param v Vector to refer.
*/
	CoordinateIndexedSparseMatrixTimesVecMinusVec(const CoordinateIndexedSparseMatrixTimesVec & c, const Vector & v) : co(c), ve(v) {} ;
	
	/** \brief Produce the Vector corresponding to this operation represented by this structure.*/
	operator const Vector() const 
	{
		Vector ret(0., ve.size()) ;
// 		std::cout << ve.size() << std::endl ;
		for (size_t i = 0 ; i < co.sm.row_size.size(); i++)
		{
			Vector temp = co.sm[i*co.sm.stride]*co.ve ;
			for(size_t j = i*co.sm.stride ; j < i*co.sm.stride+co.sm.stride ; j++)
			{
// 				std::cout << j << std::endl ;
				ret[j] = temp[j-i*co.sm.stride] - ve[j];
			}
		}
		return ret ;
	}
	
} ;

/** \brief Structure contaning the elements for a (matrix+matrix) * vector - vector operation 
 This structure allows lazy evaluation of (matrix+matrix) * vector-vector, minimizing the use of temporaries, for better numeric performance.
 */
struct CompositeSparseMatrixTimesVecMinusVec
{
	const CompositeSparseMatrixTimesVec & co ;
	const Vector & ve ;
	/** \brief Constructor.
	Initialises the references.
	@param c CompositeSparseMatrixTimesVec to refer.
	@param v Vector to refer.
*/
	CompositeSparseMatrixTimesVecMinusVec(const CompositeSparseMatrixTimesVec & c, const Vector & v) : co(c), ve(v) {} ;
	
	/** \brief overloaded operator to substract with a Vector
	 The result of the operation is a temporary structure used for lazy evaluation.
	@param v vector to substract.
	*/
	CompositeSparseMatrixTimesVecMinusVecMinusVec operator - (const Vector & v) const ;
	
	/** \brief Produce the Vector corresponding to this operation represented by this structure.*/
	operator const Vector() const 
	{
		Vector ret(0., ve.size()) ;
		
		for (size_t i = 0 ; i < co.sm.sm.row_size.size(); i++)
		{
			Vector temp = co.sm.sm[i*co.sm.sm.stride]*co.ve ;
			for(size_t j = i*co.sm.sm.stride ; j < i*co.sm.sm.stride+co.sm.sm.stride ; j++)
				ret[j] = temp[j-i*co.sm.sm.stride] - ve[j];
		}
		return ret ;
	}
	
} ;

/** \brief Structure contaning the elements for a (matrix+matrix) * vector - vector - vector operation 
 This structure allows lazy evaluation of (matrix+matrix) * vector - vector - vector, minimizing the use of temporaries, for better numeric performance.
 */
struct CompositeSparseMatrixTimesVecMinusVecMinusVec
{
	const CompositeSparseMatrixTimesVecMinusVec & co ;
	const Vector & ve ;
		/** \brief Constructor.
	Initialises the references.
	@param c CompositeSparseMatrixTimesVecMinusVec to refer.
	@param v Vector to refer.
*/
	CompositeSparseMatrixTimesVecMinusVecMinusVec(const CompositeSparseMatrixTimesVecMinusVec & c, const Vector & v) : co(c), ve(v) {} ;
	
	/** \brief Produce the Vector corresponding to this operation represented by this structure.*/
	operator const Vector() const 
	{
		Vector ret(0., ve.size()) ;
		
		for (size_t i = 0 ; i < co.co.sm.sm.row_size.size(); i++)
		{
			Vector temp = co.co.sm.sm[i*co.co.sm.sm.stride]*co.co.ve ;
			for(size_t j = i*co.co.sm.sm.stride ; j < i*co.co.sm.sm.stride+co.co.sm.sm.stride ; j++)
				ret[j] = temp[j-i*co.co.sm.sm.stride] - co.ve[j] - ve[j];
		}
		return ret ;
	}
	
} ;

/** \brief parrallel compute matrix*vector + vector and assign.
 @param ret Vector to store the result
 @param c CoordinateIndexedSparseMatrixTimesVecPlusVec containing the reference necessary for the operation
 */
void assign(Vector & ret, const Amie::CoordinateIndexedSparseMatrixTimesVecPlusVec & c, const int rowstart = 0, const int colstart = 0) ;

/** \brief parrallel compute matrix*vector - vector and assign.
	@param ret Vector to store the result
	@param c CoordinateIndexedSparseMatrixTimesVecMinusVec containing the reference necessary for the operation
 */
void assign(Vector & ret, const Amie::CoordinateIndexedSparseMatrixTimesVecMinusVec & c, const int rowstart = 0, const int colstart = 0) ;

/** \brief parrallel compute matrix*vector and assign.
	@param ret Vector to store the result
@param c CoordinateIndexedSparseMatrixTimesVec containing the reference necessary for the operation
 */
void assign(Vector & ret, const Amie::CoordinateIndexedSparseMatrixTimesVec & c, const int rowstart = 0, const int colstart = 0) ;

} ;

Amie::CompositeSparseMatrix operator+(const Amie::CoordinateIndexedSparseMatrix  & sm, const Amie::CoordinateIndexedIncompleteSparseMatrix  & ism) ;

Amie::CompositeSparseMatrix operator+(const Amie::CoordinateIndexedIncompleteSparseMatrix  & ism, const Amie::CoordinateIndexedSparseMatrix  & sm ) ;

#endif

