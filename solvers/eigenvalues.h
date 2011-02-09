#ifndef EIGEN_H
#define EIGEN_H

#include <numeric>
#include <valarray>
#include "../sparse/sparse_matrix.h"

namespace Mu
{
/** \brief Utility. Return largest EigenValue of the argument*/
  double largestEigenValue(const CoordinateIndexedSparseMatrix & A, bool sym = true) ;

/** \brief Utility. Return smallest EigenValue of the argument*/
  double smallestEigenValue(const CoordinateIndexedSparseMatrix & A, bool sym = true) ;
  
}



#endif 
