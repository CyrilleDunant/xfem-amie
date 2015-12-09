// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../main.h"
#include "../../utilities/tensor.h"


#include <fstream>
#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h>
#include <sys/time.h>


using namespace Amie ;

int main( int argc, char *argv[] )
{
    Matrix m3d = Tensor::isotropicTransverseCauchyGreen( 10e9, 25e9, 30e9, 0.2, 0.2, SPACE_THREE_DIMENSIONAL, PLANE_STRESS ) ;
    Matrix rot = Tensor::rotate4thOrderTensor3D( m3d, Point(-0.1,0,M_PI*0.5) ) ;
    Matrix C = Tensor::to2D(rot, PLANE_STRAIN, XI) ;
    
    std::cout << "----" << std::endl ;
    m3d.print() ;

    std::cout << "----" << std::endl ;
    rot.print() ;

    std::cout << "----" << std::endl ;
    C.print() ;

    return 0 ;
}
