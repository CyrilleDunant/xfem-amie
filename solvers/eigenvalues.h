// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011

#ifndef EIGEN_H
#define EIGEN_H

#include <numeric>
#include <valarray>
#include "../sparse/sparse_matrix.h"

namespace Amie
{
    struct Assembly ;
/** \brief Utility. Return largest EigenValue of the argument*/
  double largestEigenValue( Amie::Assembly* a, bool sym = false ) ;

/** \brief Utility. Return smallest EigenValue of the argument*/
  double smallestEigenValue( Amie::Assembly* a, bool sym = false ) ;
  
}



#endif 
