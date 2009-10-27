#include <smmintrin.h>
#include <iostream>
#include <limits>

typedef union 
{
	__m128d vec ;
	double val[2] ;
} vecdouble;

int main()
{
	vecdouble temp ;
	temp.val[0] = 1 ; temp.val[1] = 2 ;
	temp.vec = _mm_dp_pd(temp.vec, temp.vec, 61) ;
	
	if(std::abs(temp.val[1] - 4) > std::numeric_limits<double>::epsilon())
		exit(1) ;
	std::cout << temp.val[0] << std::endl ;
	return 0 ;
}
