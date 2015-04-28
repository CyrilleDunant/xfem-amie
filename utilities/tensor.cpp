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
			cg[0][0] = 1. ; cg[0][1] = nu ; cg[0][2] = 0 ;
			cg[1][0] = nu ; cg[1][1] = 1. ; cg[1][2] = 0 ;
			cg[2][0] = 0 ; cg[2][1] = 0 ; cg[2][2] = (1.-nu) ;
			cg *= E/(1.-nu*nu) ;
			}
			else
			{
			cg[0][0] = 1.-nu ; cg[0][1] = nu ; cg[0][2] = 0 ;
			cg[1][0] = nu ; cg[1][1] = 1.-nu ; cg[1][2] = 0 ;
			cg[2][0] = 0 ; cg[2][1] = 0 ; cg[2][2] = (1-2.*nu) ;
			cg *= E/((1.+nu)*(1.-2.*nu)) ;
			}
			return cg ;
		}
		case SPACE_THREE_DIMENSIONAL:
		{
			Matrix cgg(6,6) ;
			cgg[0][0] = 1. - nu ; cgg[0][1] = nu ; cgg[0][2] = nu ;
			cgg[1][0] = nu ; cgg[1][1] = 1. - nu ; cgg[1][2] = nu ;
			cgg[2][0] = nu ; cgg[2][1] = nu ; cgg[2][2] = 1. - nu ;
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
	
				cg[0][0] = k ; cg[0][1] = k ; cg[0][2] = 0 ;
				cg[1][0] = k ; cg[1][1] = k ; cg[1][2] = 0 ;
				cg[2][0] = 0 ; cg[2][1] = 0 ; cg[2][2] = 0 ;
				return cg ;
			}
			case SPACE_THREE_DIMENSIONAL:
			{
				Matrix cgg(6,6) ;
				cgg[0][0] = k ; cgg[0][1] = k ; cgg[0][2] = k ;
				cgg[1][0] = k ; cgg[1][1] = k ; cgg[1][2] = k ;
				cgg[2][0] = k ; cgg[2][1] = k ; cgg[2][2] = k ;
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
			cg[0][0] = 1. ; cg[0][1] = nu ; cg[0][2] = 0 ;
			cg[1][0] = nu ; cg[1][1] = 1. ; cg[1][2] = 0 ;
			cg[2][0] = 0 ; cg[2][1] = 0 ; cg[2][2] = (1.-nu) ;
			cg *= E/(1.-nu*nu) ;
			}
			else
			{
			cg[0][0] = 1.-nu ; cg[0][1] = nu ; cg[0][2] = 0 ;
			cg[1][0] = nu ; cg[1][1] = 1.-nu ; cg[1][2] = 0 ;
			cg[2][0] = 0 ; cg[2][1] = 0 ; cg[2][2] = (1-2.*nu) ;
			cg *= E/((1.+nu)*(1.-2.*nu)) ;
			}
			return cg ;
		}
		case SPACE_THREE_DIMENSIONAL:
		{
			Matrix cgg(6,6) ;
			cgg[0][0] = 1. - nu ; cgg[0][1] = nu ; cgg[0][2] = nu ;
			cgg[1][0] = nu ; cgg[1][1] = 1. - nu ; cgg[1][2] = nu ;
			cgg[2][0] = nu ; cgg[2][1] = nu ; cgg[2][2] = 1. - nu ;
			cgg[3][3] = 0.5 - nu ;
			cgg[4][4] = 0.5 - nu ;
			cgg[5][5] = 0.5 - nu ;
			cgg *= E/((1.+nu)*(1.-2.*nu)) ;
			return cgg ;
		}
	}
	return Matrix(0,0) ;
}

Matrix Tensor::orthothropicCauchyGreen(double E_1, double E_2, double G,  double nu, planeType pt) 
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
			double nu_21 = (nu/E_1)*sqrt(E_1*E_2) ;
// 			double nu_12 = (nu/E_2)*sqrt(E_1*E_2) ;

			double gamma = 1./(1.-nu*nu) ;

			cg[0][0] = E_1*gamma ; cg[0][1] = nu_21*E_1*gamma ;
			cg[1][0] = cg[0][1] ;  cg[1][1] = E_2*gamma ;
			
			G = E_1*E_2/(E_1*(1.-nu*nu)+E_2*(1.-nu*nu)) ;
			
			cg[2][2] = G ;
		}
		else if(E_1 > POINT_TOLERANCE)
		{
			cg[0][0] = E_1 ;
			cg[2][2] = 0 ;

		}
		else if(E_2 > POINT_TOLERANCE)
		{
			cg[1][1] = E_2 ;
			cg[2][2] = 0 ;
		}
		else
			cg.array() = 0 ;
	}
	else if (pt == PLANE_STRAIN)
	{
		if(E_1 > POINT_TOLERANCE && E_2 > POINT_TOLERANCE)
		{
			Matrix A(2,2) ;

			double nu21 = (nu/E_1)*sqrt(E_1*E_2) ;
			double nu12 = (nu/E_2)*sqrt(E_1*E_2) ;
			double nu23 = nu ;
			double nu32 = nu ;
			double nu13 = nu ;
			double nu31 = nu ;
			double nupe = 1.-nu21*nu12-nu23*nu32-nu13*nu31-nu12*nu23*nu31-nu21*nu32*nu13 ;

			cg[0][0] = E_1*(1.-nu32*nu23)/nupe ; cg[0][1] = (nu21-nu23*nu31)*E_1/nupe;
			cg[1][0] = cg[0][1] ;  cg[1][1] = E_2*(1.-nu32*nu23)/nupe ;
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

Matrix Tensor::orthothropicCauchyGreen(double E_1, double E_2, double E_3, double G_1, double G_2, double G_3,  double nu) 
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
	cg[0][1] = E_1*(nu_21+nu_31*nu_23)*gamma ; cg[1][0] = cg[0][1] ;
	cg[0][2] = E_1*(nu_31+nu_21*nu_32)*gamma ; cg[2][0] = cg[0][2] ;
	cg[1][2] = E_2*(nu_32+nu_12*nu_31)*gamma ; cg[2][1] = cg[1][2] ;
	cg[3][3] = G_1*0.5 ;
	cg[4][4] = G_2*0.5 ;
	cg[5][5] = G_3*0.5 ;
	
	return cg ;
}


