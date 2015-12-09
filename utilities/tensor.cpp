#include "tensor.h"

using namespace Amie ;

Tensor::Tensor(size_t o, size_t d) : order(o), dim(d)
{
    size_t s = 1 ;
    for(size_t i = 0 ; i < o ; i++)
        s *= dim ;
    components.resize(s) ;
    components = 0. ;
}

Tensor::Tensor(const Point & p, std::vector<Variable> v) : order(1), dim(v.size())
{
    size_t s = 1 ;
    for(size_t i = 0 ; i < order ; i++)
        s *= dim ;
    components.resize(s) ;
    components = 0. ;
    for(size_t i = 0 ; i < v.size() ; i++)
    {
        double d = 0 ;
        switch(v[i])
        {
        case ONE:
            d = 1 ;
        case XI:
            d = p.getX() ;
            break ;
        case ETA:
            d = p.getY() ;
            break ;
        case ZETA:
            d = p.getZ() ;
            break ;
        case TIME_VARIABLE:
            d = p.getT() ;
            break ;
        default:
            d = 0. ;
        }
        components[i] = d ;
    }
}

Tensor::Tensor(const Matrix & m, bool transpose) : order(2), dim(m.numRows())
{
    components.resize(dim*dim) ;
    for(size_t i = 0 ; i < dim ; i++)
    {
        for(size_t j = 0 ; j < dim ; j++)
        {
            if(transpose)
                (*this)(j,i) = m[i][j] ;
            else
                (*this)(i,j) = m[i][j] ;
        }
    }
}

size_t Tensor::index( std::vector<int> indexes ) const
{
    if(indexes.size() == 1)
        return indexes[0] ;
    size_t i = 0 ;
    size_t s = 1 ;
    for(int k = order ; k > 0 ; k--)
    {
        i += (indexes[k-1]) * s ;
        s *= dim ;
    }
    return i ;
}

size_t Tensor::index( int i ) const
{
    std::vector<int> idx ;
    idx.push_back(i) ;
    return index(idx) ;
}

size_t Tensor::index( int i, int j ) const
{
    std::vector<int> idx ;
    idx.push_back(i) ;
    idx.push_back(j) ;
    return index(idx) ;
}

size_t Tensor::index( int i, int j, int k ) const
{
    std::vector<int> idx ;
    idx.push_back(i) ;
    idx.push_back(j) ;
    idx.push_back(k) ;
    return index(idx) ;
}

size_t Tensor::index( int i, int j, int k, int l ) const
{
    std::vector<int> idx ;
    idx.push_back(i) ;
    idx.push_back(j) ;
    idx.push_back(k) ;
    idx.push_back(l) ;
    return index(idx) ;
}

size_t Tensor::index( int i, int j, int k, int l, int m ) const
{
    std::vector<int> idx ;
    idx.push_back(i) ;
    idx.push_back(j) ;
    idx.push_back(k) ;
    idx.push_back(l) ;
    idx.push_back(m) ;
    return index(idx) ;
}

size_t Tensor::index( int i, int j, int k, int l, int m, int n ) const
{
    std::vector<int> idx ;
    idx.push_back(i) ;
    idx.push_back(j) ;
    idx.push_back(k) ;
    idx.push_back(l) ;
    idx.push_back(m) ;
    idx.push_back(n) ;
    return index(idx) ;
}


void Tensor::print() const
{
    for(size_t i = 0 ; i < components.size() ; i++)
    {
        std::cout << components[i] << "\t" ;
    }
    std::cout << std::endl ;
}

void Tensor::threshold(double thr)
{
    for(size_t i = 0 ; i < components.size() ; i++)
    {
        if(std::abs(components[i]) < thr)
            components[i] = 0 ;
    }

}

Matrix Tensor::toMatrix(int d1, int d2) const
{
    Matrix ret ;
    if(d1+d2 != (int)order)
        return ret ;

    Tensor dummy1( d1, dim ) ;
    Tensor dummy2( d2, dim ) ;

    ret = Matrix(dummy1.size(), dummy2.size()) ;
    for(size_t i = 0 ; i < dummy1.size() ; i++)
    {
        for(size_t j = 0 ; j < dummy2.size() ; j++)
        {
            ret[i][j] = components[ i*dummy2.size() + j ] ;
        }
    }

    return ret ;

}

Matrix Tensor::cauchyGreen(std::pair<double,double> prop, bool hooke, SpaceDimensionality dim, planeType pt)
{
    double E = prop.first ;
    double nu = prop.second ;
    if(!hooke)
    {
        double k = prop.first ;
        double mu = prop.second ;
        E = 9.*k*mu / (3.*k+mu) ;
        nu = (3.*k-2.*mu) / (6.*k+2.*mu) ;
    }
    switch(dim)
    {
    case SPACE_ONE_DIMENSIONAL:
    {
        Matrix m(1,1) ;
        m[0][0] = E ;
        return m ;
    }

    case SPACE_TWO_DIMENSIONAL:
    {
        Matrix cg(3,3) ;

        if(pt == PLANE_STRESS)
        {
            cg[0][0] = 1. ;
            cg[0][1] = nu ;
            cg[0][2] = 0 ;
            cg[1][0] = nu ;
            cg[1][1] = 1. ;
            cg[1][2] = 0 ;
            cg[2][0] = 0 ;
            cg[2][1] = 0 ;
            cg[2][2] = (1.-nu) ;
            cg *= E/(1.-nu*nu) ;
        }
        else
        {
            cg[0][0] = 1.-nu ;
            cg[0][1] = nu ;
            cg[0][2] = 0 ;
            cg[1][0] = nu ;
            cg[1][1] = 1.-nu ;
            cg[1][2] = 0 ;
            cg[2][0] = 0 ;
            cg[2][1] = 0 ;
            cg[2][2] = (1-2.*nu)  ;
            cg *= E/((1.+nu)*(1.-2.*nu)) ;
        }
        return cg ;
    }
    case SPACE_THREE_DIMENSIONAL:
    {
        Matrix cgg(6,6) ;
        cgg[0][0] = 1. - nu ;
        cgg[0][1] = nu ;
        cgg[0][2] = nu ;
        cgg[1][0] = nu ;
        cgg[1][1] = 1. - nu ;
        cgg[1][2] = nu ;
        cgg[2][0] = nu ;
        cgg[2][1] = nu ;
        cgg[2][2] = 1. - nu ;
        cgg[3][3] = (0.5 - nu)*.5 ;
        cgg[4][4] = (0.5 - nu)*.5 ;
        cgg[5][5] = (0.5 - nu)*.5 ;
        cgg *= E/((1.+nu)*(1.-2.*nu)) ;
        return cgg ;
    }
    }
    return Matrix(0,0) ;
}

Matrix Tensor::cauchyGreen(double p1, double p2, bool hooke, SpaceDimensionality dim, planeType pt)
{
    double E = p1 ;
    double nu = p2 ;
    if(!hooke)
    {
        double k = p1 ;
        double mu = p2 ;
        if(mu < POINT_TOLERANCE)
        {
            switch(dim)
            {
            case SPACE_ONE_DIMENSIONAL:
            {
                Matrix m(1,1) ;
                m[0][0] = k ;
                return m ;
            }
            case SPACE_TWO_DIMENSIONAL:
            {
                Matrix cg(3,3) ;

                cg[0][0] = k ;
                cg[0][1] = k ;
                cg[0][2] = 0 ;
                cg[1][0] = k ;
                cg[1][1] = k ;
                cg[1][2] = 0 ;
                cg[2][0] = 0 ;
                cg[2][1] = 0 ;
                cg[2][2] = 0 ;
                return cg ;
            }
            case SPACE_THREE_DIMENSIONAL:
            {
                Matrix cgg(6,6) ;
                cgg[0][0] = k ;
                cgg[0][1] = k ;
                cgg[0][2] = k ;
                cgg[1][0] = k ;
                cgg[1][1] = k ;
                cgg[1][2] = k ;
                cgg[2][0] = k ;
                cgg[2][1] = k ;
                cgg[2][2] = k ;
                cgg[3][3] = 0 ;
                cgg[4][4] = 0 ;
                cgg[5][5] = 0 ;
                return cgg ;
            }
            }
        }
        E = 9*k*mu / (3*k+mu) ;
        nu = (3*k-2*mu) / (6*k+2*mu) ;
    }
    switch(dim)
    {
    case SPACE_ONE_DIMENSIONAL:
    {
        Matrix m(1,1) ;
        m[0][0] = E ;
        return m ;
    }
    case SPACE_TWO_DIMENSIONAL:
    {
        Matrix cg(3,3) ;

        if(pt == PLANE_STRESS)
        {
            cg[0][0] = 1. ;
            cg[0][1] = nu ;
            cg[0][2] = 0 ;
            cg[1][0] = nu ;
            cg[1][1] = 1. ;
            cg[1][2] = 0 ;
            cg[2][0] = 0 ;
            cg[2][1] = 0 ;
            cg[2][2] = (1.-nu) ;
            cg *= E/(1.-nu*nu) ;
        }
        else
        {
            cg[0][0] = 1.-nu ;
            cg[0][1] = nu ;
            cg[0][2] = 0 ;
            cg[1][0] = nu ;
            cg[1][1] = 1.-nu ;
            cg[1][2] = 0 ;
            cg[2][0] = 0 ;
            cg[2][1] = 0 ;
            cg[2][2] = (1-2.*nu) ;
            cg *= E/((1.+nu)*(1.-2.*nu)) ;
        }
        return cg ;
    }
    case SPACE_THREE_DIMENSIONAL:
    {
        Matrix cgg(6,6) ;
        cgg[0][0] = 1. - nu ;
        cgg[0][1] = nu ;
        cgg[0][2] = nu ;
        cgg[1][0] = nu ;
        cgg[1][1] = 1. - nu ;
        cgg[1][2] = nu ;
        cgg[2][0] = nu ;
        cgg[2][1] = nu ;
        cgg[2][2] = 1. - nu ;
        cgg[3][3] = 0.5 - nu ;
        cgg[4][4] = 0.5 - nu ;
        cgg[5][5] = 0.5 - nu ;
        cgg *= E/((1.+nu)*(1.-2.*nu)) ;
        return cgg ;
    }
    }
    return Matrix(0,0) ;
}

Matrix Tensor::orthotropicCauchyGreen(double E_1, double E_2, double G,  double nu, double angle, planeType pt)
{
    Matrix cg = Tensor::orthotropicCauchyGreen(E_1, E_2, G, nu, pt ) ;
    Matrix transform(3,3) ;
    Matrix transformt(3,3) ;
    double c = cos(angle) ;
    double s = sin(angle) ;
    transform[0][0] =  c*c ;
    transform[0][1] = s*s ;
    transform[0][2] =  2.*s*c ;
    transform[1][0] =  s*s ;
    transform[1][1] = c*c ;
    transform[1][2] = -2.*s*c ;
    transform[2][0] = -s*c ;
    transform[2][1] = s*c ;
    transform[2][2] =     c*c - s*s ;
    transformt = transform.transpose() ;
    return (transform*cg)*transformt ;
}

Matrix Tensor::orthotropicCauchyGreen(double E_1, double E_2, double G,  double nu, planeType pt)
{

    Matrix cg(3,3) ;

    if(pt == PLANE_STRESS)
    {
        if(E_1 > POINT_TOLERANCE && E_2 > POINT_TOLERANCE)
        {
// 			Matrix A(2,2) ;
// 			A[0][0] = E_1/std::max(E_1, E_2) ; A[0][1] = -E_2/std::max(E_1, E_2) ;
// 			A[1][0] = 0 ; A[1][1] = 1. ;
// 			Vector b(2) ; b[0] = 0 ; b[1] = nu ;
// 			b = inverse2x2Matrix(A)*b ;
// 			double nu_12 = b[1];
// 			double nu_21 = b[0] ;

            //nu_12*nu_21 = nu*nu ;
            double nu_21 = (nu/E_1)*(E_1+E_2)*.5 ;
            if(nu < POINT_TOLERANCE)
                nu_21 =  0 ;
            double gamma = 1./(1.-nu*nu) ;

            cg[0][0] = E_1*gamma ;
            cg[0][1] = nu_21*E_1*gamma ;
            cg[1][0] = cg[0][1] ;
            cg[1][1] = E_2*gamma ;

            G = E_1*E_2/(E_1*(1.-nu*nu)+E_2*(1.-nu*nu)) ;

            cg[2][2] = G ;
        }
        else if(E_1 > POINT_TOLERANCE)
        {
            cg[0][0] = E_1 ;
            cg[1][1] = E_1*POINT_TOLERANCE ;
            cg[2][2] = 0 ;

        }
        else if(E_2 > POINT_TOLERANCE)
        {
            cg[0][0] = E_2*POINT_TOLERANCE ;
            cg[1][1] = E_2 ;
            cg[2][2] = 0 ;
        }
        else
            cg.array() = POINT_TOLERANCE ;
    }
    else if (pt == PLANE_STRAIN)
    {
        if(E_1 > POINT_TOLERANCE && E_2 > POINT_TOLERANCE)
        {
            Matrix A(2,2) ;

            double nu21 = (nu/E_1)*(E_1+E_2)*.5 ;
            if(nu < POINT_TOLERANCE)
                nu21 =  0 ;
            double nu12 = (nu/E_2)*(E_1+E_2)*.5 ;
            if(nu < POINT_TOLERANCE)
                nu12 =  0 ;
            double nu23 = nu ;
            double nu32 = nu ;
            double nu13 = nu ;
            double nu31 = nu ;
            double nupe = 1.-nu21*nu12-nu23*nu32-nu13*nu31-nu12*nu23*nu31-nu21*nu32*nu13 ;

            cg[0][0] = E_1*(1.-nu32*nu23)/nupe ;
            cg[0][1] = (nu21-nu23*nu31)*E_1/nupe;
            cg[1][0] = cg[0][1] ;
            cg[1][1] = E_2*(1.-nu32*nu23)/nupe ;
            cg[2][2] = E_1*E_2/(E_2*(1.+nu12)*(1.-2.*nu12)+E_1*(1.+nu21)*(1.-2.*nu21)) ;
        }
        else if(E_1 > POINT_TOLERANCE)
        {
            cg[0][0] = E_1 ;
            cg[2][2] = G ;
        }
        else if(E_2 > POINT_TOLERANCE)
        {
            cg[1][1] = E_2 ;
            cg[2][2] = G ;
        }
        else
            cg.array() = 0 ;
    }
    return cg ;
}

Matrix Tensor::orthotropicCauchyGreen(double E_1, double E_2, double E_3, double G_1, double G_2, double G_3,  double nu)
{

    Matrix cg(6,6) ;
    double nu_12 = nu ;
    double nu_13 = nu ;
    double nu_23 = nu ;
    double nu_21 = nu_12*E_2/E_1 ;
    double nu_31 = nu_13*E_3/E_1 ;
    double nu_32 = nu_23*E_3/E_2 ;
    double gamma = 1./(1.-nu_12*nu_21-nu_23*nu_32-nu_13*nu_31-2.*nu_12*nu_32*nu_13) ;
    cg[0][0] = E_1*(1.-nu_23*nu_32)*gamma ;
    cg[1][1] = E_2*(1.-nu_13*nu_31)*gamma ;
    cg[2][2] = E_3*(1.-nu_12*nu_21)*gamma ;
    cg[0][1] = E_1*(nu_21+nu_31*nu_23)*gamma ;
    cg[1][0] = cg[0][1] ;
    cg[0][2] = E_1*(nu_31+nu_21*nu_32)*gamma ;
    cg[2][0] = cg[0][2] ;
    cg[1][2] = E_2*(nu_32+nu_12*nu_31)*gamma ;
    cg[2][1] = cg[1][2] ;
    cg[3][3] = G_1*0.5 ;
    cg[4][4] = G_2*0.5 ;
    cg[5][5] = G_3*0.5 ;

    return cg ;
}

Matrix Tensor::isotropicTransverseCauchyGreen(double E_1, double E_3, double G_1,  double nu_12, double nu_23, SpaceDimensionality dim, planeType pt) 
{
     if(dim == SPACE_TWO_DIMENSIONAL)
        return orthotropicCauchyGreen(E_1, E_3, G_1, nu_12, pt) ;


    Matrix cg(6,6) ;
    double E_2 = E_1 ;
    double G_2 = G_1 ;
    double G_3 = E_3/(2.*(1.+nu_23)) ;
    double nu_13 = nu_12 ;
    double nu_21 = nu_12*E_2/E_1 ;
    double nu_31 = nu_13*E_3/E_1 ;
    double nu_32 = nu_23*E_3/E_2 ;
    double gamma = 1./(1.-nu_12*nu_21-nu_23*nu_32-nu_13*nu_31-2.*nu_12*nu_32*nu_13) ;
    cg[0][0] = E_1*(1.-nu_23*nu_32)*gamma ;
    cg[1][1] = E_2*(1.-nu_13*nu_31)*gamma ;
    cg[2][2] = E_3*(1.-nu_12*nu_21)*gamma ;

    cg[0][1] = E_1*(nu_21+nu_31*nu_23)*gamma ;
    cg[1][0] = cg[0][1] ;
    cg[0][2] = E_1*(nu_31+nu_21*nu_32)*gamma ;

    cg[2][0] = cg[0][2] ;
    cg[1][2] = E_2*(nu_32+nu_12*nu_31)*gamma ;
    cg[2][1] = cg[1][2] ;

    cg[3][3] = G_1*0.5 ;
    cg[4][4] = G_2*0.5 ;
    cg[5][5] = G_3*0.5 ;

    return cg ;
}


Vector Tensor::rotate2ndOrderTensor2D( Vector & tensor, double angle )
{
    double c = cos ( angle ) ;
    double s = sin ( angle ) ;
    Matrix nrot ( 3,3 ) ;
    nrot[0][0] = c*c ;
    nrot[0][1] = s*s ;
    nrot[0][2] = -2.*s*c ;
    nrot[1][0] = s*s ;
    nrot[1][1] = c*c ;
    nrot[1][2] = 2.*s*c ;
    nrot[2][0] = s*c ;
    nrot[2][1] = -s*c ;
    nrot[2][2] = c*c-s*s ;
    return nrot*tensor ;
}

Matrix Tensor::to2D(Matrix & tensor, planeType pt, Variable var)
{
    Matrix ret(3,3) ;
    if(tensor.numCols() != 6 || tensor.numRows() != 6)
        return ret ;

    int index = var-1 ;
    int first = 0 ;
    int second = 1 ;
    if(index == 0)
        first = 2 ;
    if(index == 1)
        second = 2 ; 

    switch(pt)
    {
        case PLANE_STRESS:
        {
            ret[0][0] = tensor[first][first]-tensor[first][index]*tensor[first][index]/tensor[index][index] ;
            ret[0][1] = tensor[first][second]-tensor[second][index]*tensor[first][index]/tensor[index][index] ;
            ret[1][0] = tensor[second][first]-tensor[first][index]*tensor[second][index]/tensor[index][index] ;
            ret[1][1] = tensor[second][second]-tensor[second][index]*tensor[second][index]/tensor[index][index] ;
            ret[2][2] = tensor[index+3][index+3]*2 ;
            ret[0][2] = 0 ; ret[1][2] = 0 ; ret[2][0] = 0 ; ret[2][1] = 0 ; 
            break;
        }
        case PLANE_STRAIN:
        {
            ret[0][0] = tensor[first][first] ;
            ret[0][1] = tensor[first][second] ;
            ret[1][0] = tensor[second][first] ;
            ret[1][1] = tensor[second][second] ;
            ret[2][2] = tensor[index+3][index+3]*2 ;
            ret[0][2] = 0 ; ret[1][2] = 0 ; ret[2][0] = 0 ; ret[2][1] = 0 ; 
            break;
        }
        default:
            ret = 0 ;
    }
    return ret ;
}

Matrix rotationMatrixX(double theta)
{
    Matrix ret(6,6) ; ret = 0 ;
    double c = cos(theta) ;
    double s = sin(theta) ;

    ret[0][0] = 1 ;
    ret[1][1] = c*c ; ret[1][2] = s*s ;
    ret[2][1] = s*s ; ret[2][2] = c*c ;
    ret[1][3] = 2*c*s ;
    ret[2][3] = -2*c*s ;
    ret[3][1] = -c*s ;
    ret[3][2] = c*s ;
    ret[3][3] = c*c-s*s ;
    ret[4][4] = c ; ret[4][5] = -s ;
    ret[5][5] = c ; ret[5][4] = s ;

    return ret ;
}

Matrix rotationMatrixY(double theta)
{
    Matrix ret(6,6) ; ret = 0 ;
    double c = cos(theta) ;
    double s = sin(theta) ;

    ret[0][0] = c*c ; ret[0][2] = s*s ;
    ret[1][1] = 1 ;
    ret[2][0] = s*s ; ret[2][2] = c*c ;
    ret[0][4] = 2*c*s ;
    ret[2][4] = -2*c*s ;
    ret[4][0] = -c*s ;
    ret[4][2] = c*s ;
    ret[3][3] = c ; ret[3][5] = -s ;
    ret[4][4] = c*c-s*s ;
    ret[5][5] = c ; ret[5][3] = s ;

    return ret ;
}

Matrix rotationMatrixZ(double theta)
{
    Matrix ret(6,6) ; ret = 0 ;
    double c = cos(theta) ;
    double s = sin(theta) ;

    ret[0][0] = c*c ; ret[0][1] = s*s ;
    ret[1][0] = s*s ; ret[1][1] = c*c ;
    ret[2][2] = 1 ;
    ret[0][5] = 2*c*s ;
    ret[1][5] = -2*c*s ;
    ret[5][0] = -c*s ;
    ret[5][1] = c*s ;
    ret[5][5] = c*c-s*s ;
    ret[3][3] = c ; ret[4][3] = -s ;
    ret[4][4] = c ; ret[3][4] = s ;

    return ret ;
}

Matrix Tensor::rotate4thOrderTensor3D( Matrix & tensor, Point angle ) 
{
    Matrix ret(6,6) ;
    if(tensor.numCols() != 6 || tensor.numRows() != 6)
        return ret ;

    Matrix K = rotationMatrixX( angle.getX() ) ; 
    K *= rotationMatrixY( angle.getY() ) ;
    K *= rotationMatrixZ( angle.getZ() ) ;

//    Matrix K(6,6) ;
/*    for(size_t i = 0 ; i < 3 ; i++)
    {
        for(size_t j = 0 ; j < 3 ; j++)
        {
            K[i+3][j+3] *= 2 ;
        }
    }*/
    Matrix Kt = K.transpose() ;

    K.print() ;

    ret = (K*tensor)*Kt ;
    for(size_t i = 0 ; i < 3 ; i++)
    {
        for(size_t j = 0 ; j < 3 ; j++)
        {
            ret[3+j][i] = 0 ; ret[i][3+j] = 0 ;
            if( j != i ) 
                ret[3+i][3+j] = 0 ;
        }
    }

    return ret ;
}



