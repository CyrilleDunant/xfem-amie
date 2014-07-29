// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011

#ifndef EIGEN_H
#define EIGEN_H

#include <numeric>
#include <valarray>
#include "../sparse/sparse_matrix.h"

namespace Amie
{
/** \brief Utility. Return largest EigenValue of the argument*/
  double largestEigenValue(const CoordinateIndexedSparseMatrix & A, bool sym = false) ;

/** \brief Utility. Return smallest EigenValue of the argument*/
  double smallestEigenValue(const CoordinateIndexedSparseMatrix & A, bool sym = false) ;
  
}



#endif 
