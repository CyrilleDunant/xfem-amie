#include "tensor.h"
#include "enumeration_translator.h"
#include "../physics/homogenization/composite.h"

 namespace Amie {


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

size_t Tensor::index( const std::vector<int> & indexes ) const
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

std::vector<int> Tensor::getIndex( size_t position ) const 
{
    std::vector<int> ret( order ) ;
    size_t next = position ;
    for(size_t i = 0 ; i < order ; i++)
    {
        size_t frac = std::pow( dim, order-i-1 ) ;
        ret[i] = next / frac ;
        next = next - ret[i]*frac ;
    }
    return ret ;
}

size_t Tensor::index( int i ) const
{
    return index({i}) ;
}

size_t Tensor::index( int i, int j ) const
{

    return index({i, j}) ;
}

size_t Tensor::index( int i, int j, int k ) const
{
    return index({i, j, k}) ;
}

size_t Tensor::index( int i, int j, int k, int l ) const
{
    return index({i, j, k, l}) ;
}

size_t Tensor::index( int i, int j, int k, int l, int m ) const
{

    return index({i, j, k, l, m}) ;
}

size_t Tensor::index( int i, int j, int k, int l, int m, int n ) const
{
    return index({i, j, k, l, m, n}) ;
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
    if(thr < 0)
        return ;
    for(size_t i = 0 ; i < components.size() ; i++)
    {
        if(std::abs(components[i]) < thr)
            components[i] = 0 ;
    }

}

Tensor Tensor::crossProduct( const Tensor &t1, const Tensor &t2, double tol ) 
{
    Tensor ret( t1.getOrder()+t2.getOrder(), t1.getDimensions()) ;

    std::vector<int> index1 = t1.getIndex( 0 ) ;
    std::vector<int> index2 = t2.getIndex( 0 ) ;
    std::vector<int> index3 = ret.getIndex( 0 ) ;

    for(size_t i = 0 ; i < t1.size() ; i++)
    {
        index1 = t1.getIndex( i ) ;
        for(size_t ord = 0 ; ord < t1.getOrder() ; ord++)
            index3[ ord ] = index1[ ord ] ;

        for(size_t j = 0 ; j < t2.size() ; j++)
        {
            index2 = t2.getIndex( j ) ;
            for(size_t ord = 0 ; ord < t2.getOrder() ; ord++)
                index3[ ord + t1.getOrder() ] = index2[ ord ] ;

            ret.getComponent( index3 ) += t1.getComponent( index1 ) * t2.getComponent( index2 ) ;
        }
    }

    if(tol > 0) { ret.threshold( tol ) ; }

    return ret ;
}

void Tensor::swap( int i, int j )
{
    Vector keep = components ;
    std::vector<int> idx = getIndex(0) ;

    for(size_t pos = 0 ; pos < components.size() ; pos++)
    {
        idx = getIndex( pos ) ;
        int tmp = idx[i] ;
        idx[i] = idx[j] ;
        idx[j] = tmp ;
        components[pos] = keep[ index( idx ) ] ;
    }
}

Tensor Tensor::dotProduct( const Tensor &t1, const Tensor &t2, int i1, int i2, double tol ) 
{
    Tensor ret(t1.getOrder()+t2.getOrder()-2, t1.getDimensions()) ;
    Tensor f1( t1.getOrder()-1, t1.getDimensions()) ;
    Tensor f2( t2.getOrder()-1, t2.getDimensions()) ;


    std::vector<int> index1 = t1.getIndex( 0 ) ;
    std::vector<int> index2 = t2.getIndex( 0 ) ;
    std::vector<int> indexf1 = f1.getIndex( 0 ) ;
    std::vector<int> indexf2 = f2.getIndex( 0 ) ;
    std::vector<int> index3 = ret.getIndex( 0 ) ;

    for(size_t i = 0 ; i < f1.size() ; i++)
    {
        indexf1 = f1.getIndex( i ) ;
        bool found1 = false ;
        for(size_t k = 0 ; k < index1.size() ; k++)
        {
            if(k == (size_t) i1)
                found1 = true ;
            else
            {
                index1[k] = indexf1[k-found1] ;
                index3[k-found1] = indexf1[k-found1] ; 
            }
        }


        for(size_t j = 0 ; j < f2.size() ; j++)
        {
            indexf2 = f2.getIndex( j ) ;
            bool found2 = false ;
            for(size_t k = 0 ; k < index2.size() ; k++)
            {
                if(k == (size_t) i2)
                    found2 = true ;
                else
                {
                    index2[k] = indexf2[k-found2] ;
                    index3[k-found2+indexf1.size()] = indexf2[k-found2] ;  
                }
            }

            for( size_t k = 0 ; k < ret.getDimensions() ; k++)
            {
                index1[ i1 ] = k ;
                index2[ i2 ] = k ;
                ret.getComponent( index3 ) += t1.getComponent( index1 )*t2.getComponent( index2 ) ;
            }
        }
    }

    if(tol > 0) { ret.threshold( tol ) ; }

    return ret ;    
}

Tensor Tensor::chainedDotProduct( const Tensor &tmi, const Tensor &tnj,const  Tensor &tok, const Tensor &tpl, const Tensor &tijkl, double tol  ) 
{
    size_t dim = tijkl.getDimensions() ;
    Tensor ret(4, dim) ;

    for(size_t i = 0 ; i < dim ; i++)
    {
    for(size_t j = 0 ; j < dim ; j++)
    {
    for(size_t k = 0 ; k < dim ; k++)
    {
    for(size_t l = 0 ; l < dim ; l++)
    {
    for(size_t m = 0 ; m < dim ; m++)
    {
    for(size_t n = 0 ; n < dim ; n++)
    {
    for(size_t o = 0 ; o < dim ; o++)
    {
    for(size_t p = 0 ; p < dim ; p++)
    {
        ret(i, j, k, l) += tmi(i, m)*tnj(j, n)*tok(k, o)*tpl(l, p)*tijkl(m, n, o, p) ;
    }

    }

    }

    }

    }

    }

    }

    }


//     ret.print();
    if(tol > 0) { ret.threshold( tol ) ; }
    
    return ret ;    
}

Tensor Tensor::chainedDotProduct( const Tensor & tmi, const Tensor & tnj, const Tensor &tij, double tol  ) 
{
    size_t dim = tij.getDimensions() ;
    Tensor ret(2, dim) ;

    for(size_t i = 0 ; i < dim ; i++)
    {
    for(size_t j = 0 ; j < dim ; j++)
    {
    for(size_t m = 0 ; m < dim ; m++)
    {
    for(size_t n = 0 ; n < dim ; n++)
    {
        ret(m,n) += tmi(m,i)*tnj(n,j)*tij(i,j) ;
    }

    }

    }

    }

    if(tol > 0) { ret.threshold( tol ) ; }

    return ret ;    
}


Tensor Tensor::dotProduct( const Tensor &t1, const Tensor &t2, double tol ) 
{
    Tensor ret(t1.getOrder()+t2.getOrder()-2, t1.getDimensions()) ;
    Tensor f1( t1.getOrder()-1, t1.getDimensions()) ;
    Tensor f2( t2.getOrder()-1, t2.getDimensions()) ;


    std::vector<int> index1 = t1.getIndex( 0 ) ;
    std::vector<int> index2 = t2.getIndex( 0 ) ;
    std::vector<int> indexf1 = f1.getIndex( 0 ) ;
    std::vector<int> indexf2 = f2.getIndex( 0 ) ;
    std::vector<int> index3 = ret.getIndex( 0 ) ;

    for(size_t i = 0 ; i < f1.size() ; i++)
    {
        indexf1 = f1.getIndex( i ) ;
        for(size_t k = 0 ; k < indexf1.size() ; k++)
        {
            index1[k] = indexf1[k] ;
            index3[k] = indexf1[k] ; 
        }


        for(size_t j = 0 ; j < f2.size() ; j++)
        {
            indexf2 = f2.getIndex( j ) ;
            for(size_t k = 0 ; k < indexf2.size() ; k++)
            {
                index2[k+1] = indexf2[k] ;
                index3[k+indexf1.size()] = indexf2[k] ; 
            }

            for( size_t k = 0 ; k < ret.getDimensions() ; k++)
            {
                index1[ indexf1.size() ] = k ;
                index2[ 0 ] = k ;
                ret.getComponent( index3 ) += t1.getComponent( index1 )*t2.getComponent( index2 ) ;
            }
        }
    }

    if(tol > 0) { ret.threshold( tol ) ; }

    return ret ;    
}

Tensor Tensor::dotdotProduct( const Tensor &t1, const Tensor &t2, double tol ) 
{
    Tensor ret(t1.getOrder()+t2.getOrder()-4, t1.getDimensions()) ;
    Tensor f1( t1.getOrder()-2, t1.getDimensions()) ;
    Tensor f2( t2.getOrder()-2, t2.getDimensions()) ;


    std::vector<int> index1 = t1.getIndex( 0 ) ;
    std::vector<int> index2 = t2.getIndex( 0 ) ;
    std::vector<int> indexf1 = f1.getIndex( 0 ) ;
    std::vector<int> indexf2 = f2.getIndex( 0 ) ;
    std::vector<int> index3 = ret.getIndex( 0 ) ;

    for(size_t i = 0 ; i < f1.size() ; i++)
    {
        indexf1 = f1.getIndex( i ) ;
        for(size_t k = 0 ; k < indexf1.size() ; k++)
        {
            index1[k] = indexf1[k] ;
            index3[k] = indexf1[k] ; 
        }


        for(size_t j = 0 ; j < f2.size() ; j++)
        {
            indexf2 = f2.getIndex( j ) ;
            for(size_t k = 0 ; k < indexf2.size() ; k++)
            {
                index2[k+2] = indexf2[k] ;
                index3[k+indexf1.size()] = indexf2[k] ; 
            }

            for( size_t k = 0 ; k < ret.getDimensions() ; k++)
            {
                for( size_t l = 0 ; l < ret.getDimensions() ; l++)
                {
                    index1[ indexf1.size() ] = k ;
                    index1[ indexf1.size()+1 ] = l ;
                    index2[ 0 ] = l ;
                    index2[ 1 ] = k ;

                    ret.getComponent( index3 ) += t1.getComponent( index1 )*t2.getComponent( index2 ) ;
                }
            }
        }
    }

    if(tol > 0) { ret.threshold( tol ) ; }

    return ret ;    
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

std::pair<Matrix, Matrix> rotationMatrix2D(double angle) 
{
    Matrix transform(3,3) ;
    Matrix transformt(3,3) ;
    double c = cos(angle) ;
    double s = sin(angle) ;
    transform[0][0] =  c*c ;
    transform[0][1] = s*s ;
    transform[0][2] = 2.* s*c ;
    transform[1][0] =  s*s ;
    transform[1][1] = c*c ;
    transform[1][2] = -2.*s*c ;
    transform[2][0] = -s*c ;
    transform[2][1] = s*c ;
    transform[2][2] =     c*c - s*s ;
    transformt = transform.transpose() ;

/*    transform[0][2] *= 2. ;
    transform[1][2] *= 2. ;
    transformt[0][2] *= 2 ;
    transformt[1][2] *= 2 ;*/

    return std::make_pair( transform, transformt ) ;
}



Matrix Tensor::rotate4thOrderTensor2D( const Matrix & tensor, double angle, double tol ) 
{
    Tensor stiff(4,2) ;
    stiff(0,0,0,0) = tensor[0][0] ;
    stiff(1,1,1,1) = tensor[1][1] ;
    stiff(0,0,1,1) = tensor[0][1] ;
    stiff(1,1,0,0) = tensor[1][0] ;

    stiff(0,0,0,1) = tensor[0][2]*0.5 ;
    stiff(0,0,1,0) = tensor[0][2]*0.5 ;
    stiff(1,1,0,1) = tensor[1][2]*0.5 ;
    stiff(1,1,1,0) = tensor[1][2]*0.5 ;

    stiff(1,0,0,0) = tensor[2][0]*0.5 ;
    stiff(0,1,0,0) = tensor[2][0]*0.5 ;
    stiff(1,0,1,1) = tensor[2][1]*0.5 ;
    stiff(0,1,1,1) = tensor[2][1]*0.5 ;

    stiff(0,1,0,1) = tensor[2][2]*0.5 ;
    stiff(0,1,1,0) = tensor[2][2]*0.5 ;
    stiff(1,0,0,1) = tensor[2][2]*0.5 ;
    stiff(1,0,1,0) = tensor[2][2]*0.5 ;

    Tensor rot(2,2) ;
    rot(0,0) =  cos(angle) ;
    rot(0,1) = -sin(angle) ;
    rot(1,0) =  sin(angle) ;
    rot(1,1) =  cos(angle) ;

    Tensor rotated = Tensor::chainedDotProduct( rot, rot, rot, rot, stiff ) ;
    if(tol > 0) { rotated.threshold( tol ) ; }

    Matrix ret(3,3) ;
    ret[0][0] = rotated(0,0,0,0) ;
    ret[1][1] = rotated(1,1,1,1) ;
    ret[0][1] = rotated(0,0,1,1) ;
    ret[1][0] = rotated(1,1,0,0) ;

    ret[0][2] = (rotated(0,0,0,1)+rotated(0,0,1,0)) ;//*0.5 ;
    ret[1][2] = (rotated(1,1,0,1)+rotated(1,1,1,0)) ;//*0.5 ;
    ret[2][0] = (rotated(0,1,0,0)+rotated(1,0,0,0)) ;//*0.5 ;
    ret[2][1] = (rotated(0,1,1,1)+rotated(1,0,1,1)) ;//*0.5 ;

    ret[2][2] = (rotated(0,1,0,1)+rotated(0,1,1,0)+rotated(1,0,0,1)+rotated(1,0,1,0))*0.5 ;

    return ret ;

/*    std::pair<Matrix, Matrix> transformation = rotationMatrix2D(angle) ;

    Matrix tmp = tensor ;
    tmp[0][2] *= 0.5 ;
    tmp[1][2] *= 0.5 ;
    tmp[2][2] *= 0.5 ;
    tmp[2][1] *= 0.5 ;
    tmp[2][0] *= 0.5 ;

    Matrix ret = (transformation.first*tmp)*transformation.second ;
    for(size_t i = 0 ; i < 3 ; i++)
    {
        for(size_t j = 0 ; j < 3 ; j++)
        {
            if( std::abs( ret[i][j] ) < tol )
                ret[i][j] = 0 ;
        }
    }
    ret[0][2] /= 0.5 ;
    ret[1][2] /= 0.5 ;
    ret[2][2] /= 0.5 ;
    ret[2][1] /= 0.5 ;
    ret[2][0] /= 0.5 ;

    return ret ;*/
}

Matrix Tensor::orthotropicCauchyGreen(double E_1, double E_2, double G,  double nu, double angle, planeType pt)
{
    Matrix cg = Tensor::orthotropicCauchyGreen(E_1, E_2, G, nu, pt ) ;
    return Tensor::rotate4thOrderTensor2D( cg, angle ) ;
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
        if(E_1 > POINT_TOLERANCE*POINT_TOLERANCE && E_2 > POINT_TOLERANCE*POINT_TOLERANCE)
        {
            double nu21 = (nu/E_1)*(E_1+E_2)*.5 ;
            if(nu < POINT_TOLERANCE*POINT_TOLERANCE)
                nu21 =  0 ;
            double nu12 = (nu/E_2)*(E_1+E_2)*.5 ;
            if(nu < POINT_TOLERANCE*POINT_TOLERANCE)
                nu12 =  0 ;
            double nu23 = std::min(nu, sqrt(nu21)) ;
            double nu32 = std::min(nu, sqrt(nu21)) ;
            double nu13 = std::min(nu, sqrt(nu21)) ;
            double nu31 = std::min(nu, sqrt(nu21)) ;
            double nupe = 1.-nu21*nu12-nu23*nu32-nu13*nu31-nu12*nu23*nu31-nu21*nu32*nu13 ;
	    while(nupe < 0)
	    {
		nu21 *=.95 ;
		nu12 *=.95 ;
		nu23 = std::min(nu, sqrt(nu21)) ;
		nu32 = std::min(nu, sqrt(nu21)) ;
		nu13 = std::min(nu, sqrt(nu21)) ;
		nu31 = std::min(nu, sqrt(nu21)) ;
		nupe = 1.-nu21*nu12-nu23*nu32-nu13*nu31-nu12*nu23*nu31-nu21*nu32*nu13 ;
	    }
	    //no softening or hardening due to nu
	    double factor = nupe/((1-nu*nu)*(1.-nu32*nu23)) ;
	    

            cg[0][0] = E_1*(1.-nu32*nu23)/nupe ;
            cg[0][1] = (nu21-nu23*nu31)*E_1/nupe ;
            cg[1][0] = cg[0][1] ;
            cg[1][1] = E_2*(1.-nu32*nu23)/nupe ;
            cg[2][2] = E_1*E_2/(E_2*(1.+nu12)*(1.-2.*nu12)+E_1*(1.+nu21)*(1.-2.*nu21)) ;
	    cg *= factor ;
        }
        else if(E_1 > POINT_TOLERANCE*POINT_TOLERANCE)
        {
            cg[0][0] = E_1 ;
            cg[2][2] = G ;
        }
        else if(E_2 > POINT_TOLERANCE*POINT_TOLERANCE)
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

Matrix Tensor::orthotropicCauchyGreen( Vector data, SymmetryType sym, bool force ) 
{
    Matrix cg(6,6) ;
    if(data.size() < 2)
        return cg ;

    switch(sym)
    {
        case SYMMETRY_CUBIC:
        {
            if(data.size() != 3 && !force) { return cg ; }
            double C11 = data[0] ;
            double C44 = data[1] ;
            double C12 = data[2] ;
            cg[0][0] = C11 ; cg[0][1] = C12 ; cg[0][2] = C12 ;
            cg[1][0] = C12 ; cg[1][1] = C11 ; cg[1][2] = C12 ;
            cg[2][0] = C12 ; cg[2][1] = C12 ; cg[2][2] = C11 ;
            cg[3][3] = C44 ;
            cg[4][4] = C44 ;
            cg[5][5] = C44 ;
            break ;
        }
        case SYMMETRY_HEXAGONAL:
        {
            if(data.size() != 5  && !force) { return cg ; }
            double C11 = data[0] ;
            double C33 = data[1] ;
            double C44 = data[2] ;
            double C12 = data[3] ;
            double C13 = data[4] ;
            double C66 = (C11-C12)*0.5 ;
            cg[0][0] = C11 ; cg[0][1] = C12 ; cg[0][2] = C13 ;
            cg[1][0] = C12 ; cg[1][1] = C11 ; cg[1][2] = C13 ;
            cg[2][0] = C13 ; cg[2][1] = C13 ; cg[2][2] = C33 ;
            cg[3][3] = C44 ;
            cg[4][4] = C44 ;
            cg[5][5] = C66 ;
            break ;
        }
        case SYMMETRY_MONOCLINIC:
        {
            if(data.size() != 13 && !force) { return cg ; }
            double C11 = data[0] ;
            double C22 = data[1] ;
            double C33 = data[2] ;
            double C44 = data[3] ;
            double C55 = data[4] ;
            double C66 = data[5] ;
            double C12 = data[6] ;
            double C13 = data[7] ;
            double C16 = data[8] ;
            double C23 = data[9] ;
            double C26 = data[10] ;
            double C36 = data[11] ;
            double C45 = data[12] ;
            cg[0][0] = C11 ; cg[0][1] = C12 ; cg[0][2] = C13 ; cg[0][5] = C16 ;
            cg[1][0] = C12 ; cg[1][1] = C22 ; cg[1][2] = C23 ; cg[1][5] = C26 ;
            cg[2][0] = C13 ; cg[2][1] = C23 ; cg[2][2] = C33 ; cg[2][5] = C36 ;
            cg[3][3] = C44 ; cg[3][4] = C45 ;
            cg[4][3] = C45 ; cg[4][4] = C55 ;
            cg[5][0] = C16 ; cg[5][1] = C26 ; cg[5][2] = C36 ; cg[5][5] = C66 ;
            break ;
        }
        case SYMMETRY_ORTHORHOMBIC:	

        {
            if(data.size() != 9 && !force) { return cg ; }
            double C11 = data[0] ;
            double C22 = data[1] ;
            double C33 = data[2] ;
            double C44 = data[3] ;
            double C55 = data[4] ;
            double C66 = data[5] ;
            double C12 = data[6] ;
            double C13 = data[7] ;
            double C23 = data[8] ;
            cg[0][0] = C11 ; cg[0][1] = C12 ; cg[0][2] = C13 ;
            cg[1][0] = C12 ; cg[1][1] = C22 ; cg[1][2] = C23 ;
            cg[2][0] = C13 ; cg[2][1] = C23 ; cg[2][2] = C33 ;
            cg[3][3] = C44 ; 
            cg[4][4] = C55 ;
            cg[5][5] = C66 ;
            break ;
        }
        case SYMMETRY_TETRAGONAL:
        {
            if(data.size() != 7 && !force) { return cg ; }
            double C11 = data[0] ;
            double C33 = data[1] ;
            double C44 = data[2] ;
            double C66 = data[3] ;
            double C12 = data[4] ;
            double C13 = data[5] ;
            double C16 = data[6] ;
            if( std::abs(C66) < POINT_TOLERANCE ) { C66 = (C11-C12)/2 ; } 
            cg[0][0] = C11 ; cg[0][1] = C12 ; cg[0][2] = C13 ; cg[0][5] = C16 ;
            cg[1][0] = C12 ; cg[1][1] = C11 ; cg[1][2] = C13 ; cg[1][5] = -C16 ;
            cg[2][0] = C13 ; cg[2][1] = C13 ; cg[2][2] = C33 ;
            cg[3][3] = C44 ;
            cg[4][4] = C44 ;
            cg[5][0] = C16 ; cg[1][5] = -C16 ; cg[5][5] = C66 ;
            break ;
        }
        case SYMMETRY_TRIGONAL:
        {
            if(data.size() != 7 && !force) { return cg ; }
            double C11 = data[0] ;
            double C33 = data[1] ;
            double C44 = data[2] ;
            double C66 = data[3] ;
            double C12 = data[4] ;
            double C13 = data[5] ;
            double C14 = data[6] ;
            if( std::abs(C66) < POINT_TOLERANCE ) { C66 = (C11-C12)/2 ; } 
            cg[0][0] = C11 ; cg[0][1] = C12 ; cg[0][2] = C13 ; cg[0][3] = C14 ;
            cg[1][0] = C12 ; cg[1][1] = C11 ; cg[1][2] = C13 ; cg[1][3] = -C14 ;
            cg[2][0] = C13 ; cg[2][1] = C13 ; cg[2][2] = C33 ;
            cg[3][0] = C14 ; cg[3][1] = -C14 ; cg[3][3] = C44 ;
            cg[4][4] = C44 ; cg[4][5] = C14 ;
            cg[5][4] = C14 ; cg[5][5] = C66 ;
            break ;
        }
        case SYMMETRY_TRICLINIC:
        {
            if(data.size() != 21 && !force) { return cg ; }
            double C11 = data[0] ;
            double C22 = data[1] ;
            double C33 = data[2] ;
            double C44 = data[3] ;
            double C55 = data[4] ;
            double C66 = data[5] ;
            double C12 = data[6] ;
            double C13 = data[7] ;
            double C14 = data[8] ;
            double C15 = data[9] ;
            double C16 = data[10] ;
            double C23 = data[11] ;
            double C24 = data[12] ;
            double C25 = data[13] ;
            double C26 = data[14] ;
            double C34 = data[15] ;
            double C35 = data[16] ;
            double C36 = data[17] ;
            double C45 = data[18] ;
            double C46 = data[19] ;
            double C56 = data[20] ;
            cg[0][0] = C11 ; cg[0][1] = C12 ; cg[0][2] = C13 ; cg[0][3] = C14 ; cg[0][4] = C15 ; cg[0][5] = C16 ;
            cg[1][0] = C12 ; cg[1][1] = C22 ; cg[1][2] = C23 ; cg[1][3] = C24 ; cg[1][4] = C25 ; cg[1][5] = C26 ;
            cg[2][0] = C13 ; cg[2][1] = C23 ; cg[2][2] = C33 ; cg[2][3] = C34 ; cg[2][4] = C35 ; cg[2][5] = C36 ;
            cg[3][0] = C14 ; cg[3][1] = C24 ; cg[3][2] = C34 ; cg[3][3] = C44 ; cg[3][4] = C45 ; cg[3][5] = C46 ;
            cg[4][0] = C15 ; cg[4][1] = C25 ; cg[4][2] = C35 ; cg[4][3] = C45 ; cg[4][4] = C55 ; cg[4][5] = C56 ;
            cg[5][0] = C16 ; cg[5][1] = C26 ; cg[5][2] = C36 ; cg[5][3] = C46 ; cg[5][4] = C56 ; cg[5][5] = C66 ;
            break ;
        }
    }

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


Matrix Tensor::rotate2ndOrderTensor2D( const Matrix & tensor, double angle )
{

    Matrix r(2,2) ;
    r[0][0] = cos(angle) ; r[1][1] = r[0][0] ;
    r[1][0] = sin(angle) ; r[0][1] = -r[1][0] ;
    Matrix rt = r ;
    rt[1][0] = r[0][1] ;
    rt[0][1] = r[1][0] ;

    Matrix tmp = (r*tensor)*rt ;

    return tmp ;
}

Vector Tensor::rotate2ndOrderTensor2D( const Vector & tensor, double angle )
{
    Matrix t(2,2) ;
    t[0][0] = tensor[0] ;
    t[1][1] = tensor[1] ;
    t[0][1] = tensor[2] ;
    t[1][0] = tensor[2] ;

    Matrix r(2,2) ;
    r[0][0] = cos(angle) ; r[1][1] = r[0][0] ;
    r[1][0] = sin(angle) ; r[0][1] = -r[1][0] ;
    Matrix rt = r ;
    rt[1][0] = r[0][1] ;
    rt[0][1] = r[1][0] ;

    Matrix tmp = (r*t)*rt ;

    Vector ret(3) ;
    ret[0] = tmp[0][0] ;
    ret[1] = tmp[1][1] ;
    ret[2] = (tmp[0][1]+tmp[1][0])*0.5 ;

    return ret ;
}

Matrix Tensor::to2D(Matrix & tensor, planeType pt, Variable var)
{
    Matrix ret(3,3) ;
    if(tensor.numCols() != 6 || tensor.numRows() != 6)
        return ret ;

    Matrix moved = tensor ;
    switch(var)
    {
        case XI:
            moved = Tensor::rotate4thOrderTensor3D( tensor, Point(0,M_PI*0.5,0) ) ; 
            break ;
        case ETA:
            moved = Tensor::rotate4thOrderTensor3D( tensor, Point(M_PI*0.5,M_PI*0.5,0) ) ; 
            break ;
        case ZETA:
            break ;
        default:
            break ;
    }

    int index = 2 ;// if(index < 0) { index = 0 ; } if(index > 2) { index = 2 ; }
    int first = 0 ;
    int second = 1 ;
/*    if(index == 0)
        first = 2 ;
    if(index == 1)
        second = 2 ; */

    switch(pt)
    {
        case PLANE_STRESS:
        {
            Matrix K11(3,3) ;
            Matrix K12(3,3) ;
            Matrix K21(3,3) ;
            Matrix K22(3,3) ;

            K11[0][0] = moved[first][first] ; K11[0][1] = moved[first][second] ; K11[0][2] = moved[first][index+3] ;
            K11[1][0] = moved[second][first] ; K11[1][1] = moved[second][second] ; K11[1][2] = moved[second][index+3] ;
            K11[2][0] = moved[index+3][first] ; K11[2][1] = moved[index+3][second] ; K11[2][2] = moved[index+3][index+3] ;

            K12[0][0] = moved[first][index] ; K12[0][1] = moved[first][first+3] ; K12[0][2] = moved[first][second+3] ;
            K12[1][0] = moved[second][index] ; K12[1][1] = moved[second][first+3] ; K12[1][2] = moved[second][second+3] ;
            K12[2][0] = moved[index+3][index] ; K12[2][1] = moved[index+3][first+3] ; K12[2][2] = moved[index+3][second+3] ;
            K21 = K12.transpose() ;


            K22[0][0] = moved[index][index] ; K22[0][1] = moved[index][first+3] ; K22[0][2] = moved[index][second+3] ;
            K22[1][0] = moved[first+3][index] ; K22[1][1] = moved[first+3][first+3] ; K22[1][2] = moved[first+3][second+3] ;
            K22[2][0] = moved[second+3][index] ; K22[2][1] = moved[second+3][first+3] ; K22[2][2] = moved[second+3][second+3] ;

            invert3x3Matrix(K22) ;

            ret = K11-K12*K22*K21 ;
//            ret[2][2] *= 2 ;

/*            ret[0][0] = tensor[first][first]-tensor[first][index]*tensor[first][index]/tensor[index][index] ;
            ret[0][1] = tensor[first][second]-tensor[second][index]*tensor[first][index]/tensor[index][index] ;
            ret[1][0] = tensor[second][first]-tensor[first][index]*tensor[second][index]/tensor[index][index] ;
            ret[1][1] = tensor[second][second]-tensor[second][index]*tensor[second][index]/tensor[index][index] ;
            ret[2][2] = tensor[index+3][index+3]*2 ;
            ret[0][2] = 0 ; ret[1][2] = 0 ; ret[2][0] = 0 ; ret[2][1] = 0 ; */
            break;
        }
        case PLANE_STRAIN:
        {
            ret[0][0] = moved[first][first] ;
            ret[0][1] = moved[first][second] ;
            ret[1][0] = moved[second][first] ;
            ret[1][1] = moved[second][second] ;
            ret[2][2] = moved[index+3][index+3]*2 ;
            ret[0][2] = moved[first][index+3] ; 
            ret[1][2] = moved[second][index+3] ; 
            ret[2][0] = moved[index+3][first] ;
            ret[2][1] = moved[index+3][second] ; 
            break;
        }
        default:
            ret = 0 ;
    }

/*    ret[0][2] = 0 ;
    ret[1][2] = 0 ;
    ret[2][0] = 0 ;
    ret[2][1] = 0 ;
    ret[0][1] = std::abs(ret[0][1]) ;
    ret[1][0] = std::abs(ret[1][0]) ;
    ret[0][0] = std::abs(ret[0][0]) ;
    ret[1][1] = std::abs(ret[1][1]) ;
    ret[2][2] = std::abs(ret[2][2]) ;*/

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
    Matrix ret(3,3) ; ret = 0 ;
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
    Matrix ret(3,3) ; ret = 0 ;
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

Matrix rotationMatrixXYZ( double psi, double theta, double phi)
{
    Matrix Om(3,3) ;
    Om[0][0] =  cos(psi)*cos(theta)*cos(phi)-sin(psi)*sin(phi) ;
    Om[1][0] =  sin(psi)*cos(theta)*cos(phi)+cos(psi)*sin(phi) ;
    Om[2][0] =  -sin(theta)*cos(phi) ;

    Om[0][1] =  -cos(psi)*cos(theta)*sin(phi)-sin(psi)*cos(phi) ;
    Om[1][1] =  -sin(psi)*cos(theta)*sin(phi)+cos(psi)*cos(phi) ;
    Om[2][1] =  sin(theta)*sin(phi) ;

    Om[0][2] =  cos(psi)*sin(theta) ;
    Om[1][2] =  sin(psi)*sin(theta) ;
    Om[2][2] =  cos(theta) ;
    return Om ;
}

Matrix rotationMatrixXYZ( const Point & direction)
{
    Point dir = direction/direction.norm() ;
    Matrix Om(3,3) ;
    Point up(0, 0, 1) ;
    if((up^dir).norm() < POINT_TOLERANCE)
        up.set(0,1,0) ;
    if((up^dir).norm() < POINT_TOLERANCE)
        up.set(1,0,0) ;
    
    Point xaxis = up^dir;
    xaxis /= xaxis.norm() ;

    Point yaxis = dir^xaxis ;
    yaxis /=yaxis.norm() ;

        Om[0][0] = xaxis.x;
        Om[0][1] = yaxis.x;
        Om[0][2] = dir.x;

        Om[1][0] = xaxis.y;
        Om[1][1] = yaxis.y;
        Om[1][2] = dir.y;

        Om[2][0] = xaxis.z;
        Om[2][1] = yaxis.z;
        Om[2][2] = dir.z;
    
    return Om ;
}

Vector Tensor::rotate2ndOrderTensor3D( const Vector & tensor, Point angle, double tol ) 
{
    Matrix Om = rotationMatrixXYZ( angle.getX(), angle.getY(), angle.getZ() ) ; 

    if(false)
    {
        Tensor rot(2,3) ;
        for(size_t i = 0 ; i < 3 ; i++)
        {
            for(size_t j = 0 ; j < 3 ; j++)
            {
                rot(i,j) = Om[i][j] ;
            }
        }

        Tensor data(2,3) ;
        data(0,0) = tensor[0] ;
        data(1,1) = tensor[1] ;
        data(2,2) = tensor[2] ;
        data(0,1) = tensor[5] ;
        data(1,0) = tensor[5] ;
        data(0,2) = tensor[4] ;
        data(2,0) = tensor[4] ;
        data(2,1) = tensor[3] ;
        data(1,2) = tensor[3] ;

        Tensor rotated = Tensor::chainedDotProduct( rot, rot, data, -1. ) ;
	if(tol > 0) { rotated.threshold(tol) ; }

        Vector ret(6) ;
        ret[0] = rotated(0,0) ;
        ret[1] = rotated(1,1) ;
        ret[2] = rotated(2,2) ;
        ret[3] = (rotated(1,2)+rotated(2,1))*0.5 ;
        ret[4] = (rotated(0,2)+rotated(2,0))*0.5 ;
        ret[5] = (rotated(1,0)+rotated(0,1))*0.5 ;

        return ret ;
    }

    Matrix tr = Om.transpose() ;

    Matrix ten(3,3) ;
    ten[0][0] = tensor[0] ; ten[0][1] = tensor[5] ; ten[0][2] = tensor[4] ;
    ten[1][0] = tensor[5] ; ten[1][1] = tensor[1] ; ten[1][2] = tensor[3] ;
    ten[2][0] = tensor[4] ; ten[2][1] = tensor[3] ; ten[2][2] = tensor[2] ;

    ten = (Om*ten)*tr ;

    Vector ret(6) ;
    ret[0] = ten[0][0] ;
    ret[1] = ten[1][1] ;
    ret[2] = ten[2][2] ;
    ret[3] = (ten[1][2]+ten[2][1])*0.5 ;
    ret[4] = (ten[0][2]+ten[2][0])*0.5 ;
    ret[5] = (ten[1][0]+ten[0][1])*0.5 ;

    if(tol > 0)
    {
        for(size_t i = 0 ; i < 6 ; i++)
            if( std::abs( ret[i] ) < tol ) { ret[i] = 0 ; }
    }

    return ret ;
}

Matrix Tensor::rotate4thOrderTensor3D( const Matrix & tensor, Point angle, double tol ) 
{
    Matrix ret(6,6) ;
    if(tensor.numCols() != 6 || tensor.numRows() != 6)
        return tensor ;

//     Matrix Om = rotationMatrixXYZ(angle.getX(), angle.getY(), angle.getZ()) ; 
    Matrix Om = rotationMatrixXYZ(angle) ; 
//     Om.print();

    if(true)
    {
        Tensor rot(2,3) ;
        Tensor trs(2,3) ;
        for(size_t i = 0 ; i < 3 ; i++)
        {
            for(size_t j = 0 ; j < 3 ; j++)
            {
                rot(i,j) = Om[i][j] ;
                trs(i,j) = Om[j][i] ;
            }
        }

        Tensor stiff(4,3) ;
        double v = 0.5 ;
        double w = 0.5 ;
        double z = 0.5 ;

        stiff(0,0,0,0) = tensor[0][0] ;
        stiff(1,1,1,1) = tensor[1][1] ;
        stiff(2,2,2,2) = tensor[2][2] ;
    
        stiff(0,0,1,1) = tensor[0][1] ;
        stiff(1,1,0,0) = tensor[1][0] ;
        stiff(0,0,2,2) = tensor[0][2] ;
        stiff(2,2,0,0) = tensor[2][0] ;
        stiff(1,1,2,2) = tensor[1][2] ;
        stiff(2,2,1,1) = tensor[2][1] ;
    
        stiff(0,0,1,2) = tensor[0][3] ;//*.5
        stiff(0,0,2,1) = tensor[0][3] ;//*.5
        stiff(0,0,2,0) = tensor[0][4] ;//*.5
        stiff(0,0,0,2) = tensor[0][4] ;//*.5
        stiff(0,0,0,1) = tensor[0][5] ;//*.5
        stiff(0,0,1,0) = tensor[0][5] ;//*.5
                                       //
        stiff(1,1,1,2) = tensor[1][3] ;//*.5
        stiff(1,1,2,1) = tensor[1][3] ;//*.5
        stiff(1,1,2,0) = tensor[1][4] ;//*.5
        stiff(1,1,0,2) = tensor[1][4] ;//*.5
        stiff(1,1,0,1) = tensor[1][5] ;//*.5
        stiff(1,1,1,0) = tensor[1][5] ;//*.5
                                       //
        stiff(2,2,1,2) = tensor[2][3] ;//*.5
        stiff(2,2,2,1) = tensor[2][3] ;//*.5
        stiff(2,2,2,0) = tensor[2][4] ;//*.5
        stiff(2,2,0,2) = tensor[2][4] ;//*.5
        stiff(2,2,0,1) = tensor[2][5] ;//*.5
        stiff(2,2,1,0) = tensor[2][5] ;//*.5
                                       //
        stiff(1,2,0,0) = tensor[3][0] ;//*.5
        stiff(2,1,0,0) = tensor[3][0] ;//*.5
        stiff(2,0,0,0) = tensor[4][0] ;//*.5
        stiff(0,2,0,0) = tensor[4][0] ;//*.5
        stiff(0,1,0,0) = tensor[5][0] ;//*.5    
        stiff(1,0,0,0) = tensor[5][0] ;//*.5
                                       //
        stiff(1,2,1,1) = tensor[3][1] ;//*.5
        stiff(2,1,1,1) = tensor[3][1] ;//*.5
        stiff(2,0,1,1) = tensor[4][1] ;//*.5
        stiff(0,2,1,1) = tensor[4][1] ;//*.5    
        stiff(0,1,1,1) = tensor[5][1] ;//*.5
        stiff(1,0,1,1) = tensor[5][1] ;//*.5
                                       //
        stiff(1,2,2,2) = tensor[3][2] ;//*.5
        stiff(2,1,2,2) = tensor[3][2] ;//*.5
        stiff(2,0,2,2) = tensor[4][2] ;//*.5
        stiff(0,2,2,2) = tensor[4][2] ;//*.5
        stiff(0,1,2,2) = tensor[5][2] ;//*.5
        stiff(1,0,2,2) = tensor[5][2] ;//*.5
                                     
        stiff(1,2,1,2) = tensor[3][3]*.5 ;
        stiff(1,2,2,1) = tensor[3][3]*.5 ;
        stiff(2,1,1,2) = tensor[3][3]*.5 ;
        stiff(2,1,2,1) = tensor[3][3]*.5 ;
        stiff(1,2,2,0) = tensor[3][4]*.5 ;
        stiff(1,2,0,2) = tensor[3][4]*.5 ;
        stiff(2,1,2,0) = tensor[3][4]*.5 ;
        stiff(2,1,0,2) = tensor[3][4]*.5 ;
        stiff(1,2,0,1) = tensor[3][5]*.5 ;
        stiff(1,2,1,0) = tensor[3][5]*.5 ;
        stiff(2,1,0,1) = tensor[3][5]*.5 ;
        stiff(2,1,1,0) = tensor[3][5]*.5 ;
                                     
        stiff(2,0,1,2) = tensor[4][3]*.5 ;
        stiff(2,0,2,1) = tensor[4][3]*.5 ;
        stiff(0,2,1,2) = tensor[4][3]*.5 ;
        stiff(0,2,2,1) = tensor[4][3]*.5 ;
        stiff(2,0,2,0) = tensor[4][4]*.5 ;
        stiff(2,0,0,2) = tensor[4][4]*.5 ;
        stiff(0,2,2,0) = tensor[4][4]*.5 ;
        stiff(0,2,0,2) = tensor[4][4]*.5 ;
        stiff(2,0,0,1) = tensor[4][5]*.5 ;
        stiff(2,0,1,0) = tensor[4][5]*.5 ;
        stiff(0,2,0,1) = tensor[4][5]*.5 ;
        stiff(0,2,1,0) = tensor[4][5]*.5 ;
                                     
        stiff(0,1,1,2) = tensor[5][3]*.5 ;
        stiff(0,1,2,1) = tensor[5][3]*.5 ;
        stiff(1,0,1,2) = tensor[5][3]*.5 ;
        stiff(1,0,2,1) = tensor[5][3]*.5 ;
        stiff(0,1,2,0) = tensor[5][4]*.5 ;
        stiff(0,1,0,2) = tensor[5][4]*.5 ;
        stiff(1,0,2,0) = tensor[5][4]*.5 ;
        stiff(1,0,0,2) = tensor[5][4]*.5 ;
        stiff(0,1,0,1) = tensor[5][5]*.5 ;
        stiff(0,1,1,0) = tensor[5][5]*.5 ;
        stiff(1,0,0,1) = tensor[5][5]*.5 ;    
        stiff(1,0,1,0) = tensor[5][5]*.5 ;

        Tensor rotated = Tensor::chainedDotProduct( rot, rot, rot, rot, stiff, tol ) ;
//         stiff.print();
//         rotated.print();
        if(tol > 0) { rotated.threshold( tol ) ; }

        ret[0][0] = rotated(0,0,0,0) ;
        ret[1][1] = rotated(1,1,1,1) ;
        ret[2][2] = rotated(2,2,2,2) ;    

        ret[0][1] = rotated(0,0,1,1) ;
        ret[1][0] = rotated(1,1,0,0) ;
        ret[0][2] = rotated(0,0,2,2) ;
        ret[2][0] = rotated(2,2,0,0) ;
        ret[1][2] = rotated(1,1,2,2) ;
        ret[2][1] = rotated(2,2,1,1) ;

        ret[3][0] = (rotated(1,2,0,0)+rotated(2,1,0,0)) ;//*0.5 ;
        ret[4][0] = (rotated(2,0,0,0)+rotated(0,2,0,0)) ;//*0.5 ;
        ret[5][0] = (rotated(0,1,0,0)+rotated(1,0,0,0)) ;//*0.5 ;
        ret[0][3] = (rotated(0,0,1,2)+rotated(0,0,2,1)) ;//*0.5 ;
        ret[0][4] = (rotated(0,0,2,0)+rotated(0,0,0,2)) ;//*0.5 ;
        ret[0][5] = (rotated(0,0,0,1)+rotated(0,0,1,0)) ;//*0.5 ;
                                                       
        ret[3][1] = (rotated(1,2,1,1)+rotated(2,1,1,1)) ;//*0.5 ;
        ret[4][1] = (rotated(2,0,1,1)+rotated(0,2,1,1)) ;//*0.5 ;
        ret[5][1] = (rotated(0,1,1,1)+rotated(1,0,1,1)) ;//*0.5 ;
        ret[1][3] = (rotated(1,1,1,2)+rotated(1,1,2,1)) ;//*0.5 ;
        ret[1][4] = (rotated(1,1,2,0)+rotated(1,1,0,2)) ;//*0.5 ;/w
        ret[1][5] = (rotated(1,1,0,1)+rotated(1,1,1,0)) ;//*0.5 ;      

        ret[3][2] = (rotated(1,2,2,2)+rotated(2,1,2,2)) ;//*0.5 ;
        ret[4][2] = (rotated(2,0,2,2)+rotated(0,2,2,2)) ;//*0.5 ;
        ret[5][2] = (rotated(0,1,2,2)+rotated(1,0,2,2)) ;//*0.5 ;
        ret[2][3] = (rotated(2,2,1,2)+rotated(2,2,2,1)) ;//*0.5 ;
        ret[2][4] = (rotated(2,2,2,0)+rotated(2,2,0,2)) ;//*0.5 ;
        ret[2][5] = (rotated(2,2,0,1)+rotated(2,2,1,0)) ;//*0.5 ;

        ret[3][3] = (rotated(1,2,1,2)+rotated(1,2,2,1)+rotated(2,1,2,1)+rotated(2,1,1,2))*.5 ;//*0.25
        ret[3][4] = (rotated(1,2,2,0)+rotated(1,2,0,2)+rotated(2,1,0,2)+rotated(2,1,2,0))*.5 ;//*0.25
        ret[3][5] = (rotated(1,2,0,1)+rotated(1,2,1,0)+rotated(2,1,1,0)+rotated(2,1,0,1))*.5 ;//*0.25
        ret[4][3] = (rotated(2,0,1,2)+rotated(2,0,2,1)+rotated(0,2,2,1)+rotated(0,2,1,2))*.5 ;//*0.25    
        ret[4][4] = (rotated(2,0,2,0)+rotated(2,0,0,2)+rotated(0,2,0,2)+rotated(0,2,2,0))*.5 ;//*0.25
        ret[4][5] = (rotated(2,0,0,1)+rotated(2,0,1,0)+rotated(0,2,1,0)+rotated(0,2,0,1))*.5 ;//*0.25
        ret[5][3] = (rotated(0,1,1,2)+rotated(0,1,2,1)+rotated(1,0,2,1)+rotated(1,0,1,2))*.5 ;//*0.25
        ret[5][4] = (rotated(0,1,2,0)+rotated(0,1,0,2)+rotated(1,0,0,2)+rotated(1,0,2,0))*.5 ;//*0.25
        ret[5][5] = (rotated(0,1,0,1)+rotated(0,1,1,0)+rotated(1,0,1,0)+rotated(1,0,0,1))*.5 ;//*0.25

        return ret ;
    }
//     Matrix test = tensor ;
//     for(size_t i = 0 ; i < 3 ; i++)
//     {
//         for(size_t j = 0 ; j < 3 ; j++)
//         {
//             test[i][j+3] *= 2. ;
//             test[i+3][j] *= 2. ;
//             test[i+3][j+3] *= 2. ;
//         }
// 
//     }
// 
// 
//     Matrix K(6,6) ;
//     for(size_t i = 0 ; i < 3 ; i++)
//     {
//         for(size_t j = 0 ; j < 3 ; j++)
//         {
//             size_t i1 = (i+1)%3 ;
//             size_t i2 = (i+2)%3 ;
//             size_t j1 = (j+1)%3 ;
//             size_t j2 = (j+2)%3 ;
//             K[i][j] = Om[i][j]*Om[i][j] ;
//             K[i][j+3] = Om[i][j1]*Om[i][j2] ;
//             K[i+3][j] = Om[i1][j]*Om[i2][j] ;
//             K[i+3][j+3] = (Om[i1][j1]*Om[i2][j2] + Om[i1][j2]*Om[i2][j1]) ;
//         }
//     }
// 
//     Matrix Kt = K.transpose() ;
// 
//     for(size_t i = 0 ; i < 3 ; i++)
//     {
//         for(size_t j = 0 ; j < 3 ; j++)
//         {
//             K[i+3][j] *= 2. ;
// //            K[i][j+3] *= 0.5 ;
// //            Kt[i+3][j] *= 0.5 ;
//             Kt[i][j+3] *= 2. ;
//         }
//     }
// 
//     ret = (K)*(test*Kt) ;
// 
//     for(size_t i = 0 ; i < 3 ; i++)
//     {
//         for(size_t j = 0 ; j < 3 ; j++)
//         {
//             ret[i][j+3] /= 2. ;
//             ret[i+3][j] /= 2. ;
//             ret[i+3][j+3] /= 2. ;
//         }
// 
//     }
// 
//     for(size_t i = 0 ; i < ret.numCols() ; i++)
//     {
//         for(size_t j = 0 ; j < ret.numCols() ; j++)
//             if(std::abs(ret[i][j]) < tol) { ret[i][j] = 0 ; }
//     }
// 
//     return ret ;
}

Matrix Tensor::invert4thOrderTensor3D( const Matrix & tensor, SymmetryType sym)
{
    Matrix ret(6,6) ;
    if(tensor.numCols() != 6 || tensor.numRows() != 6)
        return ret ;

    switch(sym)
    {
        case SYMMETRY_CUBIC:
        {
            ret = tensor ;
            Composite::invertTensor(ret) ;
            break ;  
        }
        case SYMMETRY_TRIGONAL:
        {
            double c11 = tensor[0][0] ;
            double c12 = tensor[0][1] ;
            double c13 = tensor[0][2] ;
            double c33 = tensor[2][2] ;
            double c14 = tensor[0][3] ;
            double c44 = tensor[3][3] ;
            double c66 = tensor[5][5] ;

            double det4466 = c44*c66-c14*c14 ;
   
            double s44 = c66/(det4466) ;
            double s66 = c44/(det4466) ;
            double s14 = -c14/(2*det4466) ;

            double det3313 = 2*c13*c13-c33*(c11+c12) ;
            double s33 = -(c11+c12)/det3313 ;
            double s13 = c13/det3313 ;

            double det1112 = c11*c11-c12*c12 ; 
            double d1 = 1-c13*s13-c14*s14 ;
            double d0 = 0-c13*s13+c14*s14 ;
            double s11 = (c11*d1 - c12*d0)/det1112 ;
            double s12 = (-c12*d1 + c11*d0)/det1112 ;

            Vector s(7) ; 
            s[0] = s11 ;
            s[1] = s33 ;
            s[2] = s44 ;
            s[3] = s66 ;
            s[4] = s12 ;
            s[5] = s13 ;
            s[6] = s14 ;
           
            ret = Tensor::orthotropicCauchyGreen( s, sym ) ;
            break ;
        }
        default:
            std::cout << "warning: tensor inversion not implemented for specificied symmetry" << std::endl ;
            break ;
    }
    return ret ;
}

Matrix Tensor::cauchyGreen( std::pair<double, double> props, SpaceDimensionality dim, planeType pt, IsotropicMaterialParameters param ) 
{
   return Tensor::cauchyGreen( props.first, props.second, dim, pt, param ) ;
}

Matrix Tensor::cauchyGreen( double p1, double p2, SpaceDimensionality dim, planeType pt, IsotropicMaterialParameters param ) 
{
    Matrix ret( 3+3*(dim==SPACE_THREE_DIMENSIONAL), 3+3*(dim==SPACE_THREE_DIMENSIONAL)) ;
    double Ciiii = 0. ;
    double Ciijj = 0. ;
    double Cijij = 0. ;

    switch( param )
    {
        case YOUNG_POISSON:
        {
            if(dim == SPACE_THREE_DIMENSIONAL || pt == PLANE_STRAIN)
            {
                double C = p1/((1.+p2)*(1.-2.*p2)) ;
                Ciiii = C*(1.-p2) ;
                Ciijj = C*p2 ;
                Cijij = C*(1.-2.*p2) ;
            }
            else if(pt == PLANE_STRESS)
            {
                double C = p1/(1.-p2*p2) ;
                Ciiii = C ;
                Ciijj = C*p2 ;
                Cijij = C*(1.-p2) ;
            }

            break ;
        }
        case BULK_SHEAR:
        {
            if(dim == SPACE_THREE_DIMENSIONAL || pt == PLANE_STRAIN)
            {
                Ciiii = p1+4.*p2/3. ;
                Ciijj = p1-2.*p2/3. ;
                Cijij = 2.*p2 ;
            }
            else if(pt == PLANE_STRESS)
            {
                Ciiii = 4.*p2*(1.-p2/(p1+4.*p2/3)) ;
                Ciijj = 2.*p2*(1.-2.*p2/(p1+4.*p2/3.)) ;
                Cijij = 2.*p2 ;
            }

            break ;
        }
        case YOUNG_SHEAR:
        {
            if(dim == SPACE_THREE_DIMENSIONAL || pt == PLANE_STRAIN)
            {
                Ciiii = (4.*p2-p1)/(3.-p1/p2) ;
                Ciijj = (p1-2.*p2)/(3.-p1/p2) ;
                Cijij = 2.*p2 ;
            }
            else if(pt == PLANE_STRESS)
            {
                Ciiii = p2/(1.-p1/(4.*p2)) ;
                Ciijj = p1/(4.-p1/p2)-p2 ;
                Cijij = 2.*p2 ;
            }

            break ;
        }
        default:
            break ;
    }

    if(dim == SPACE_THREE_DIMENSIONAL)
    {
       ret[0][0] = Ciiii ; ret[0][1] = Ciijj ; ret[0][2] = Ciijj ;
       ret[1][0] = Ciijj ; ret[1][1] = Ciiii ; ret[1][2] = Ciijj ;
       ret[2][0] = Ciijj ; ret[2][1] = Ciijj ; ret[2][2] = Ciiii ;
       ret[3][3] = Cijij ;
       ret[4][4] = Cijij ;
       ret[5][5] = Cijij ;
    }
    else
    {
       ret[0][0] = Ciiii ; ret[0][1] = Ciijj ;
       ret[1][0] = Ciijj ; ret[1][1] = Ciiii ;
       ret[2][2] = Cijij ;
    }


    return ret ;


}


std::pair<double, double> Tensor::getIsotropicMaterialParameters( const Matrix & C, IsotropicMaterialParameters param, planeType pt ) 
{
    std::pair<double, double> ret = std::make_pair(0,0) ;
    if( pt == PLANE_STRAIN || C.numCols() == 6)
    {
        switch(param)
        {
            case YOUNG_POISSON:
            {
                double nu = C[0][1]/(C[0][0]+C[0][1]) ;
                double E = C[0][0]*(1.+nu)*(1.-2.*nu)/(1.-nu) ;
                ret.first = E ;
                ret.second = nu ;
                break ;
            }
            case BULK_SHEAR:
            {
                double G = C[ C.numRows()-1 ][ C.numRows()-1 ]*0.5 ;
                double K = C[0][0] - G*4./3. ;
                ret.first = K ;
                ret.second = G ;
                break ;
            }
            case YOUNG_SHEAR:
            {
                double G = C[ C.numRows()-1 ][ C.numRows()-1 ]*0.5 ;
                double E = (4.*G-3.*C[0][0]) / (1.-C[0][0]/G) ;
                ret.first = E ;
                ret.second = G ;
                break ; 
            }
            default:
                break ;
        }
    }
    else
    {
        switch(param)
        {
            case YOUNG_POISSON:
            {
                double nu = C[0][1]/C[0][0] ;
                double E = C[0][0]*(1.-nu*nu) ;
                ret.first = E ;
                ret.second = nu ;
                break ;
            }
            case BULK_SHEAR:
            {
                double G = C[ C.numRows()-1 ][ C.numRows()-1 ]*0.5 ;
                double K = 4.*G*G/(4.*G-C[0][0]) - G*4./3. ;
                ret.first = K ;
                ret.second = G ;
                break ;
            }
            case YOUNG_SHEAR:
            {
                double G = C[ C.numRows()-1 ][ C.numRows()-1 ]*0.5 ;
                double E = 4.*G*(1.-G/C[0][0]) ;
                ret.first = E ;
                ret.second = G ;
                break ;
            }
            default:
                break ;
        }
    }
    return ret ;
}




}

