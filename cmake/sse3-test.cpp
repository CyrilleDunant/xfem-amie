#include <emmintrin.h>
#include <iostream>
#include <limits>
#include <cmath>

typedef union 
{
	__m128d vec ;
	double val[2] ;
} vecdouble;

int main()
{
	vecdouble temp, zero ;
	temp.vec = _mm_setzero_pd() ;
	zero.vec = _mm_setzero_pd() ;
	temp.val[0] = 1 ; temp.val[1] = 2 ;
	temp.vec = _mm_sub_pd(temp.vec, zero.vec) ;
	if(std::abs(temp.val[0] -1) > std::numeric_limits<double>::epsilon())
		exit(1) ;
	std::cout << temp.val[0] << std::endl ;
	return 0 ;
}