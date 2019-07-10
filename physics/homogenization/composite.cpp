#include "composite.h"
#include "../../utilities/matrixops.h"

namespace Amie {

Composite::Composite( DelaunayTriangle *tri, std::vector<Feature *> feats, InclusionGeometryType t, double a, double b, double c) : Phase( tri, t, a, b, c )
{

}

Composite::Composite( DelaunayTetrahedron *tet, std::vector<Feature *> feats , InclusionGeometryType t, double a, double b, double c) : Phase( tet, t, a, b, c )
{

}

Composite::Composite( const Phase & p ) : Phase( p)
{

}

Matrix Composite::I4( const Matrix &C )
{
    Matrix I4( C.numRows(), C.numCols() ) ;

    for( size_t i = 0 ; i < 3 + 3 * ( I4.size() == 36 ) ; i++ )
    {
        I4[i][i] = 1 ;
    }

    return I4 ;
}


Matrix makeEshelbyIsotropic(const Matrix & S)
{
//     bool isIsotropic = true ;
//     for(size_t i = 1 ; i < S.numRows() && isIsotropic; i++)
//     {
//         for(size_t j = i+1 ; j < S.numCols() && isIsotropic; j++)
//         {
//             if(std::abs(S[i][j]-S[j][i]) > 1e-12)
//                 isIsotropic = false ;
//         }
//     }
//     
//     if(isIsotropic)
//         return S;
    double r, s, t ;
    r = 0.5 ;
    s = (sqrt(5)+1)*.25 ;
    t = (sqrt(5)-1)*.25 ;
    std::valarray<Vector> directions = {{ r,  s,  t},//1
                                        { r, -s,  t},//2
                                        {-r,  s,  t},//3
                                        {-r, -s,  t},//4
                                        { t,  r,  s},//5
                                        { t, -r,  s},//6
                                        {-t,  r,  s},//7
                                        {-t, -r,  s},//8
                                        { s,  t,  r},//9
                                        { s, -t,  r},//10
                                        {-s,  t,  r},//11
                                        {-s, -t,  r},//
                                        { 1,  0,  0},//
                                        { 0,  1,  0},//
                                        { 0,  0,  1} //
    } ;
    Matrix acc(S.numRows(), S.numCols()) ;
    
    Matrix accc(S.numRows(), S.numCols()) ;
    Matrix da(S.numRows(), S.numCols()) ;
    Matrix y(S.numRows(), S.numCols()) ;
    Matrix tmp(S.numRows(), S.numCols()) ;

//     double n0 = S.froebeniusNorm() ;
    for(size_t i = 0 ; i < directions.size() ;i++)
    {
        da = Tensor::rotate4thOrderEshelbyTensor3D(S, directions[i]) ;
        y = da - accc ;
        tmp = acc+y ;
        accc = (tmp- acc) - y ;
        acc = tmp ;
    }       
//     double n1 = acc.froebeniusNorm() ;
//     std::cout << n0 << " " << n1 << std::endl;
    
    return acc/directions.size();
}

Matrix makeStiffnessIsotropic(const Matrix & S)
{

    double r, s, t ;
    r = 0.5 ;
    s = (sqrt(5)+1)*.25 ;
    t = (sqrt(5)-1)*.25 ;
    std::valarray<Vector> directions = {/*{ r,  s,  t},//1
                                        { r, -s,  t},//2
                                        {-r,  s,  t},//3
                                        {-r, -s,  t},//4
                                        { t,  r,  s},//5
                                        { t, -r,  s},//6
                                        {-t,  r,  s},//7
                                        {-t, -r,  s},//8
                                        { s,  t,  r},//9
                                        { s, -t,  r},//10
                                        {-s,  t,  r},//11
                                        {-s, -t,  r},//*/
                                        { 1,  0,  0},//
                                        { 0,  1,  0},//
                                        { 0,  0,  1}} ;
    Matrix acc(S.numRows(), S.numCols()) ;
    
    Matrix accc(S.numRows(), S.numCols()) ;
    Matrix da(S.numRows(), S.numCols()) ;
    Matrix y(S.numRows(), S.numCols()) ;
    Matrix tmp(S.numRows(), S.numCols()) ;
    

    for(size_t i = 0 ; i < directions.size() ;i++)
    {
        da = Tensor::rotate4thOrderStiffnessTensor3D(S, acos(directions[i][0]), acos(directions[i][1]), acos(directions[i][2])) ;
        y = da - accc ;
        tmp = acc+y ;
        accc = (tmp- acc) - y ;
        acc = tmp ;
    }       
    
    y = S - accc ;
    tmp = acc+y ;
    accc = (tmp- acc) - y ;
    acc = tmp ;
    
    return acc/(directions.size()+1);
}

Matrix makeIsotropic(const Matrix & S, const std::vector<Point> & points)
{
    Matrix acc(S.numRows(), S.numCols()) ;
    
    Point avg ;  
    int count = 0 ;
    for(size_t i = 0 ; i < points.size() ;i++)
    {
        avg += points[i]/points.size() ;
        if(points[i].getZ() < 0)
            continue ;
        
        acc +=Tensor::rotate4thOrderEshelbyTensor3D(S, points[i]) ;
        count++ ;
    }
    
    std::cerr << acc[5][0] << "  " << avg.getX()<< "  "<< avg.getY()<< "  "<< avg.getZ() <<std::endl ;
    
    return acc/count;
    
    
}

Matrix makeStiffnessIsotropic(const Matrix & S, const std::vector<Point> & points)
{
    Matrix acc(S.numRows(), S.numCols()) ;
    Matrix accc(S.numRows(), S.numCols()) ;
    Matrix da(S.numRows(), S.numCols()) ;
    Matrix y(S.numRows(), S.numCols()) ;
    Matrix tmp(S.numRows(), S.numCols()) ;
    
    for(size_t i = 0 ; i < points.size() ;i++)
    {
        da = Tensor::rotate4thOrderStiffnessTensor3D(S, acos(points[i].getX()), acos(points[i].getY()), acos(points[i].getZ())) ;
        y = da - accc ;
        tmp = acc+y ;
        accc = (tmp- acc) - y ;
        acc = tmp ;
    }       

    return acc/(points.size());
    
}


Matrix  Composite::eshelbyEllipsoid(const Matrix & C, double a, double b, double c)
{
    

    if (b > a)
        std::swap(a,b) ;
    if(c > a)
        std::swap(a,c) ;
    if(c > b)
        std::swap(c, b) ;
    
    double pi4 = M_PIf128*4. ;
    double aa = a*a ;
    double bb = b*b ;
    double cc = c*c ;
    
    Matrix S( C.numRows(), C.numCols() ) ;
    double nu = 0.2 ;
    double x = C[0][0] ;
    double y = C[0][1] ;
    if(std::abs(x+y) > 1e-6)
        nu = y / ( x + y ) ;
    else
        nu = 0.0 ;

    double theta = asin(sqrt((aa - cc)/aa)) ;
    double k =  (aa - bb) / (aa - cc) ;

    double F = 0 ;
    double E = 0 ;
    double Ec = 0 ;
    double Fc = 0 ;
    double dw = theta*1e-3 ;
    double g0 = -sqrt(3./5) ; 
    double g1 = -g0 ;
    for(double w = 0 ; w < theta-dw*.5 ; w += 2.*dw)
    {

        double s0 = sin(w+g0*dw+dw) ;
        double d0 = 1.-k*s0*s0 ;
        if( d0 >0 )
        {
            double itemp = sqrt(d0) ;
            if(itemp > POINT_TOLERANCE)
            {
                double df = 5./9. * dw / itemp ;
                double y = df - Fc ;
                double t = F+y ;
                Fc = (t- F) - y ;
                F = t ;

                double de = 5./9. * dw * itemp ;
                y = de - Ec ;
                t = E+y ;
                Ec = (t- E) - y ;
                E = t ;

            }
        }
        
        double s = sin(w+dw) ;
        double d = 1.-k*s*s ;
        if( d >0 )
        {
            double itemp = sqrt(d) ;
            if(itemp > POINT_TOLERANCE)
            {
                double df = 8./9. * dw / itemp ;
                double y = df - Fc ;
                double t = F+y ;
                Fc = (t- F) - y ;
                F = t ;

                double de = 8./9. * dw * itemp ;
                y = de - Ec ;
                t = E+y ;
                Ec = (t- E) - y ;
                E = t ;

            }
        }
        
        double s1 = sin(w+g1*dw+dw) ;
        double d1 = 1.-k*s1*s1 ;
        if( d1 >0 )
        {
            double itemp = sqrt(d1) ;
            if(itemp > POINT_TOLERANCE)
            {
                double df = 5./9. * dw / itemp ;
                double y = df - Fc ;
                double t = F+y ;
                Fc = (t- F) - y ;
                F = t ;

                double de = 5./9. * dw * itemp ;
                y = de - Ec ;
                t = E+y ;
                Ec = (t- E) - y ;
                E = t ;
            }
        }
    }
    
    
    if(std::abs(a-b) > 1e-6 && std::abs(c-b) > 1e-6) // general spheroid
    {
        double I1 = pi4*a*b*c / ((aa-bb)*sqrt(aa-cc))*(F-E) ;
        double I3 = pi4*a*b*c / ((bb-cc)*sqrt(aa-cc))*(b*sqrt(aa-cc)/(a*c)-E) ;
        double I2 = pi4-I1-I3 ;
        
        double I12 = (I2-I1)/(3.*(aa-bb)) ;      // +
        double I13 = (I3-I1)/(3.*(aa-cc)) ;      // +
        double I11= pi4/(3.*aa) - I12 - I13 ; // +
        
        double I23 = (I3-I2)/(3.*(bb-cc)) ;  //+
        double I21 = (I1-I2)/(3.*(bb-aa)) ;  //+
        double I22 = pi4/(3.*bb) - I21 - I23 ; // +
        
        double I31 = (I1-I3)/(3.*(cc-aa)) ; //+
        double I32 = (I2-I3)/(3.*(cc-bb)) ; //+
        double I33 = pi4/(3.*cc) - I31 - I32 ; // +

        double Q = M_1_PIf128 *3./(1.-nu)*.125 ;
        double R = M_1_PIf128*(1.-2.*nu)/(1.-nu)*.125 ;
        S[0][0] = Q*aa*I11+R*I1 ; // + 
        S[1][1] = Q*bb*I22+R*I2 ; // + 
        S[2][2] = Q*cc*I33+R*I3 ; // +

        S[0][1] = Q*bb*I12-R*I1 ; // +
        S[0][2] = Q*cc*I13-R*I1 ; // +
        S[1][0] = Q*aa*I21-R*I2 ;
        S[1][2] = Q*cc*I23-R*I2 ; // +
        S[2][0] = Q*aa*I31-R*I3 ;
        S[2][1] = Q*bb*I32-R*I3 ;

        S[3][3] = (aa+bb)*.5*Q*I12+0.5*R*(I1+I2) ; // +
        S[4][4] = (aa+cc)*.5*Q*I13+0.5*R*(I1+I3) ; // +
        S[5][5] = (bb+cc)*.5*Q*I23+0.5*R*(I2+I3) ; // +
        
        return S ;
    }
    else if(std::abs(a-b) <= 1e-6 && std::abs(c-b) > 1e-6) // oblate spheroid
    {
        
        double I1 = 2.*M_PIf128*aa*c/pow(aa-cc, 3./2.)*(acos(c/a)-(c/a)*sqrt(1.-cc/(aa))) ;
        double I2 = I1 ;
        double I3 = 4.*M_PIf128-I1-I2 ;
       
        double I13 = (I3-I1)/(3.*(aa-cc)) ;      // +
        double I12 = M_PIf128/(3.*aa)- I13 * .25  ;      // +
        double I11= pi4/(3.*aa) - I12 - I13 ; // +
        
        double I23 = (I3-I2)/(3.*(aa-cc)) ;  //+
        double I21 = M_PIf128/(3.*aa)- I23 * .25 ;  //+
        double I22 = pi4/(3.*aa) - I21 - I23 ; // +
        
        double I31 = (I1-I3)/(3.*(cc-aa)) ; //+
        double I32 = (I2-I3)/(3.*(cc-aa)) ; //+
        double I33 = pi4/(3.*cc) - I31 - I32 ; // +
        
        double Q = M_1_PIf128*3./(1.-nu)*.125 ;
        double R = M_1_PIf128*(1.-2.*nu)/(1.-nu)*.125 ;
        S[0][0] = Q*aa*I11+R*I1 ; // + 
        S[1][1] = Q*aa*I22+R*I2 ; // + 
        S[2][2] = Q*cc*I33+R*I3 ; // +

        S[0][1] = Q*aa*I12-R*I1 ; // +
        S[0][2] = Q*cc*I13-R*I1 ; // +
        S[1][0] = Q*aa*I21-R*I2 ;
        S[1][2] = Q*cc*I23-R*I2 ; // +
        S[2][0] = Q*aa*I31-R*I3 ;
        S[2][1] = Q*aa*I32-R*I3 ;

        S[3][3] = (aa+aa)*.5*Q*I12+0.5*R*(I1+I2) ; // +
        S[4][4] = (aa+cc)*.5*Q*I13+0.5*R*(I1+I3) ; // +
        S[5][5] = (aa+cc)*.5*Q*I23+0.5*R*(I2+I3) ; // +
        
        return S ;
    }
    else if(std::abs(a-b) > 1e-6 && std::abs(c-b) <= 1e-6) //prolate spheroid
    {
        double I3 = 2.*M_PIf128*a*cc/pow(aa-cc, 3./2.)*(a/c*sqrt(aa/(cc)-1)-acosh(a/c)) ;
        double I2 = I3;
        double I1 = 4.*M_PIf128-I3-I2 ;
        
       
        double I12 = (I2-I1)/(3.*(aa-bb)) ;      // +
        double I13 = (I3-I1)/(3.*(aa-bb)) ;      // +
        double I11 = 4.*M_PIf128/(3.*aa) - I12 - I13 ; // +
        
        double I21 = (I1-I2)/(3.*(bb-aa)) ;  //+
        double I23 = M_PIf128/(3.*bb)-I21*.25 ;  //+
        double I22 = 4.*M_PIf128/(3.*bb) - I21 - I23 ; // +

        double I31 = (I1-I3)/(3.*(cc-aa)) ; //+
        double I32 = M_PIf128/(3.*cc)-I31*.25 ; //+
        double I33 = 4.*M_PIf128/(3.*cc) - I31 - I32 ; // +

        double Q = M_1_PIf128*3./(1.-nu)*.125 ;
        double R = M_1_PIf128*(1.-2.*nu)/(1.-nu)*.125 ;
        S[0][0] = Q*aa*I11+R*I1 ; // + 
        S[1][1] = Q*bb*I22+R*I2 ; // + 
        S[2][2] = Q*bb*I33+R*I3 ; // +

        S[0][1] = Q*bb*I12-R*I1 ; // +
        S[0][2] = Q*bb*I13-R*I1 ; // +
        S[1][0] = Q*aa*I21-R*I2 ;
        S[1][2] = Q*bb*I23-R*I2 ; // +
        S[2][0] = Q*aa*I31-R*I3 ;
        S[2][1] = Q*bb*I32-R*I3 ;

        S[3][3] = (aa+bb)*.5*Q*I12+0.5*R*(I1+I2) ; // +
        S[4][4] = (aa+bb)*.5*Q*I13+0.5*R*(I1+I3) ; // +
        S[5][5] = (bb+bb)*.5*Q*I23+0.5*R*(I2+I3) ; // +
        return S ;
    }
    else // it's a sphere
        return eshelby(C) ;
    

}


Matrix  Composite::eshelbyCylinder(const Matrix & C, double a, double b)
{
    Matrix S( C.numRows(), C.numCols() ) ;
    double nu = 0.2 ;
    double x = C[0][0] ;
    double y = C[0][1] ;
    if(std::abs(x+y) > 1e-6)
        nu = y / ( x + y ) ;
    else
        nu = 0.01 ;
    
    double pi4 = M_PI*4. ;
    double aa = a*a ;
    double bb = b*b ;
    
    double I1 = pi4*b/(a+b) ;
    double I2 = pi4*a/(a+b);

    double I12 = pi4/(3.*(a+b)*(a+b)) ;      // +
    double I11 = pi4/(3.*aa) - I12 ; // +
    
    double I21 = pi4/(3.*(a+b)*(a+b)) ;  //+
    double I22 = pi4/(3.*bb) - I21 ; // +


    double Q = M_1_PI*3./(1.-nu)*.125 ;
    double R = M_1_PI*(1.-2.*nu)/(1.-nu)*.125 ;
    S[0][0] = Q*aa*I11+R*I1 ; // + 
    S[1][1] = Q*bb*I22+R*I2 ; // + 
    S[2][2] = 0 ; // +

    S[0][1] = Q*bb*I12-R*I1 ; // +
    S[0][2] = -R*I1 ; // +
    S[1][0] = Q*aa*I21-R*I2 ;
    S[1][2] = -R*I2 ; // +
    S[2][0] = 0 ;
    S[2][1] = 0 ;

    S[3][3] = (aa+bb)*.5*Q*I12+0.5*R*(I1+I2) ; // +
    S[4][4] = 0.5*R*I1 ; // +
    S[5][5] = 0.5*R*I2 ; // + 
    
    return S ;

}

Matrix Composite::eshelby(const Matrix & C )
{
    Matrix S( C.numRows(), C.numCols() ) ;
    double nu = 0.2 ;
    double x = C[0][0] ;
    double y = C[0][1] ;
    if(x+y > 1e-16)
        nu = y / ( x + y ) ;
    else
        nu = 0 ;

    if( S.size() == 9 )
    {
        double Siiii = ( 7. - 5. * nu ) / ( 15. * ( 1. - nu ) ) ;
        double Siijj = ( 5. * nu - 1. ) / ( 15. * ( 1. - nu ) ) ;
        double Sijij = ( 4. - 5. * nu ) / ( 15. * ( 1. - nu ) ) ;

        S[0][0] = Siiii ;
        S[1][1] = Siiii ;

        S[0][1] = Siijj ;
        S[1][0] = Siijj ;

        S[2][2] = Sijij ;
    }

    if( S.size() == 36 )
    {
        double Siiii = ( 7. - 5. * nu ) / ( 15. * ( 1. - nu ) ) ;
        double Siijj = ( 5. * nu - 1 ) / ( 15. * ( 1. - nu ) ) ;
        double Sijij = ( 4. - 5. * nu ) / ( 15. * ( 1. - nu ) ) ;

        S[0][0] = Siiii ;
        S[1][1] = Siiii ;
        S[2][2] = Siiii ;

        S[0][1] = Siijj ;
        S[0][2] = Siijj ;
        S[1][0] = Siijj ;
        S[1][2] = Siijj ;
        S[2][0] = Siijj ;
        S[2][1] = Siijj ;

        S[3][3] = Sijij ;
        S[4][4] = Sijij ;
        S[5][5] = Sijij ;

    }

    return S ;
}

void Composite::invertTensor( Matrix &m )
{
    Matrix m1( 3, 3 ) ;

    for( size_t i = 0 ; i < 3 ; i++ )
    {
        for( size_t j = 0 ; j < 3 ; j++ )
            m1[i][j] = m[i][j] ;
    }

   invert3x3Matrix( m1 ) ;

    for( size_t i = 0 ; i < 3 ; i++ )
    {
        for( size_t j = 0 ; j < 3 ; j++ )
            m[i][j] = m1[i][j] ;
    }

    for( size_t i = 3 ; i < m.numCols() ; i++ )
    {
        if(std::abs(m[i][i] ) > 1e-48)
            m[i][i] = 1. / m[i][i] ;
        else
            m[i][i] = 1e48 ;
    }
}

MatrixInclusionComposite::MatrixInclusionComposite( DelaunayTriangle *tri, Feature *inc, InclusionGeometryType t, double a, double b, double c ) : Composite( tri, std::vector<Feature *>( 0 ), t, a, b, c )
{
    matrix = Phase( tri, t, a, b, c ) ;
    inclusion = Phase( inc , tri , t, a, b, c) ;

    matrix.volume = std::max(matrix.volume, 1e-12) ;
    inclusion.volume = std::max(inclusion.volume, 1e-12) ;
    volume = matrix.volume + inclusion.volume ;
    matrix.volume /= volume ;
    inclusion.volume /= volume ;
    volume=1;

}

MatrixInclusionComposite::MatrixInclusionComposite( DelaunayTetrahedron *tet, Feature *inc, InclusionGeometryType t, double a, double b, double c) : Composite( tet, std::vector<Feature *>( 0 ), t, a, b, c )
{

    matrix = Phase( tet, t, a, b, c ) ;
    inclusion = Phase( inc, t, a, b, c ) ;

    matrix.volume = std::max(matrix.volume, 1e-12) ;
    inclusion.volume = std::max(inclusion.volume, 1e-12) ;
    volume = matrix.volume + inclusion.volume ;
    matrix.volume /= volume ;
    inclusion.volume /= volume ;
    volume=1;

}

MatrixInclusionComposite::MatrixInclusionComposite( const Phase & mat, const Phase & inc) : Composite( mat), matrix(mat), inclusion(inc)
{
    matrix.volume = std::max(matrix.volume, 1e-12) ;
    inclusion.volume = std::max(inclusion.volume, 1e-12) ;
    volume = matrix.volume + inclusion.volume ;
    matrix.volume /= volume ;
    inclusion.volume /= volume ;
    volume=1;

}


void MatrixInclusionComposite::apply()
{

    getStrainConcentrationTensor() ;
    Matrix I = Composite::I4( matrix.C ) ;

    inclusion.A = B * inclusion.volume + I * matrix.volume ;
    Composite::invertTensor( inclusion.A ) ;
    inclusion.A = inclusion.A*B ;
    
    if(matrix.volume > 1e-32)
        matrix.A = ( I - inclusion.A * inclusion.volume ) / matrix.volume ;
    else
    {
        matrix.A =  I  ;
    }


    C =  ( matrix.A * matrix.volume ) * matrix.C + ( inclusion.A * inclusion.volume ) * inclusion.C ;
    
    beta.resize( matrix.beta.size(), 0. );
    beta = matrix.A * matrix.volume * matrix.beta ;
    
//     std::cout <<" a = "<< matrix.A[0][0] <<  " * " << matrix.volume << " * " << matrix.beta[0] << std::endl ;
    beta += inclusion.A * inclusion.volume * inclusion.beta ;
//     std::cout <<" b = "<< inclusion.A[0][0] <<" * " << inclusion.volume << " * "<< inclusion.beta[0] <<" = "<<beta[0] << std::endl ;


}

void MatrixInclusionComposite::getStrainConcentrationTensor()
{
    B = Matrix( C.numRows(), C.numCols() ) ;
}

DiluteMatrixInclusionComposite::DiluteMatrixInclusionComposite( DelaunayTriangle *tri, Feature *inc, InclusionGeometryType t, double a, double b, double c) : MatrixInclusionComposite( tri, inc, t, a, b, c )
{

}

DiluteMatrixInclusionComposite::DiluteMatrixInclusionComposite( DelaunayTetrahedron *tet, Feature *inc  , InclusionGeometryType t, double a, double b, double c) : MatrixInclusionComposite( tet, inc, t, a, b, c )
{

}

DiluteMatrixInclusionComposite::DiluteMatrixInclusionComposite( const Phase & mat, const Phase & inc ) : MatrixInclusionComposite( mat, inc)
{

}

void DiluteMatrixInclusionComposite::getStrainConcentrationTensor()
{
    MoriTanakaMatrixInclusionComposite mori( matrix, inclusion ) ;
    mori.apply() ;

    Matrix I = Composite::I4( mori.B ) ;
    B = I - ( mori.B * inclusion.volume ) ;
    Composite::invertTensor( B ) ;
    B = B*mori.B * matrix.volume ;
}


VoigtMatrixInclusionComposite::VoigtMatrixInclusionComposite( DelaunayTriangle *tri, Feature *inc , InclusionGeometryType t, double a, double b, double c) : MatrixInclusionComposite( tri, inc, t, a, b, c )
{

}

VoigtMatrixInclusionComposite::VoigtMatrixInclusionComposite( DelaunayTetrahedron *tet, Feature *inc , InclusionGeometryType t, double a, double b, double c) : MatrixInclusionComposite( tet, inc, t, a, b, c )
{

}

VoigtMatrixInclusionComposite::VoigtMatrixInclusionComposite( const Phase & mat, const Phase & inc) : MatrixInclusionComposite( mat, inc)
{

}

void VoigtMatrixInclusionComposite::getStrainConcentrationTensor()
{
    B = Composite::I4( C ) ;
}

ReussMatrixInclusionComposite::ReussMatrixInclusionComposite( DelaunayTriangle *tri, Feature *inc, InclusionGeometryType t, double a, double b, double c) : MatrixInclusionComposite( tri, inc, t, a, b, c )
{

}

ReussMatrixInclusionComposite::ReussMatrixInclusionComposite( DelaunayTetrahedron *tet, Feature *inc , InclusionGeometryType t, double a, double b, double c) : MatrixInclusionComposite( tet, inc, t, a, b, c )
{

}

ReussMatrixInclusionComposite::ReussMatrixInclusionComposite( const Phase & mat, const Phase & inc) : MatrixInclusionComposite( mat, inc)
{

}

void ReussMatrixInclusionComposite::getStrainConcentrationTensor()
{
    if(inclusion.C.array().max() > 1e-16)
    {
        B = inclusion.C ;
        Composite::invertTensor( B );
        B *= matrix.C ;
    }
    else
    {
        std::cerr << "Reuss scheme not valid for empty inclusions // fall back on Voigt scheme" << std::endl ;
        B = Composite::I4( matrix.C) ;
    }
}

MoriTanakaMatrixInclusionComposite::MoriTanakaMatrixInclusionComposite( DelaunayTriangle *tri, Feature *inc, InclusionGeometryType t, double a, double b, double c ) : MatrixInclusionComposite( tri, inc, t, a, b, c )
{

}

MoriTanakaMatrixInclusionComposite::MoriTanakaMatrixInclusionComposite( DelaunayTetrahedron *tet, Feature *inc, InclusionGeometryType t, double a, double b, double c ) : MatrixInclusionComposite( tet, inc, t, a, b, c )
{

}

MoriTanakaMatrixInclusionComposite::MoriTanakaMatrixInclusionComposite( const Phase & mat, const Phase & inc ) : MatrixInclusionComposite( mat, inc )
{

}

void MoriTanakaMatrixInclusionComposite::getStrainConcentrationTensor()
{

//     if(matrix.C.array().max() > 1e-32)
//     {
        Matrix I = Composite::I4( C ) ;
        Matrix S = Composite::eshelby( matrix.C ) ;
        if(matrix.t == INCLUSION_IS_ELLIPSOID)
            S = Composite::eshelbyEllipsoid( matrix.C, matrix.a, matrix.b, matrix.c ) ;
        if(matrix.t == INCLUSION_IS_CYLINDER)
            S = Composite::eshelbyCylinder( matrix.C, matrix.a, matrix.b) ;
//         S = makeEshelbyIsotropic(S) ;
        B = matrix.C ;
        Composite::invertTensor( B ) ;
        B = S * B * ( inclusion.C - matrix.C ) ;
        B += I ;
        
        Composite::invertTensor( B ) ;
//     }
//     else
//     {
//         B = matrix.C ;
//     }


}

InverseMoriTanakaMatrixInclusionComposite::InverseMoriTanakaMatrixInclusionComposite( DelaunayTriangle *tri, Feature *inc , InclusionGeometryType t, double a, double b, double c) : MatrixInclusionComposite( tri, inc, t, a, b, c )
{

}

InverseMoriTanakaMatrixInclusionComposite::InverseMoriTanakaMatrixInclusionComposite( DelaunayTetrahedron *tet, Feature *inc , InclusionGeometryType t, double a, double b, double c) : MatrixInclusionComposite( tet, inc, t, a, b, c )
{

}

InverseMoriTanakaMatrixInclusionComposite::InverseMoriTanakaMatrixInclusionComposite( const Phase & mat, const Phase & inc ) : MatrixInclusionComposite( mat, inc)
{

}

void InverseMoriTanakaMatrixInclusionComposite::getStrainConcentrationTensor()
{
    if(matrix.C.array().max() > 1e-16)
    {
        Matrix I = Composite::I4( C ) ;
        Matrix S = Composite::eshelby( inclusion.C ) ;
        if(inclusion.t == INCLUSION_IS_ELLIPSOID)
        {
            S = Composite::eshelbyEllipsoid( inclusion.C, inclusion.a, inclusion.b, inclusion.c ) ;
        }
        if(inclusion.t == INCLUSION_IS_CYLINDER)
        {
            S = Composite::eshelbyCylinder( inclusion.C, inclusion.a, inclusion.b) ;
        }
//         S = makeEshelbyIsotropic(S) ;

        B = inclusion.C ;
        Composite::invertTensor( B ) ;
        B *= S * ( matrix.C - inclusion.C ) ;
        B += I ;
        
    }
    else
    {
        std::cerr << "Inverse Mori-Tanaka scheme not valid for empty matrix // fall back on Voigt scheme" << std::endl ;
        B = Composite::I4( matrix.C) ;
    }
}

BiphasicSelfConsistentComposite::BiphasicSelfConsistentComposite( DelaunayTriangle *tri, Feature *inc, InclusionGeometryType t, double a, double b, double c ) : MatrixInclusionComposite( tri, inc, t, a, b, c )
{

    ReussMatrixInclusionComposite voigt( matrix, inclusion  ) ;
    voigt.apply() ;
    fictious = voigt ;
    Matrix S = inclusion.C ;
    Matrix I = Composite::I4( S ) ;
    hint = fictious ;
    B = I - S * inclusion.volume ;
    Composite::invertTensor( B ) ;
    B *= S * matrix.volume ;
}

BiphasicSelfConsistentComposite::BiphasicSelfConsistentComposite( DelaunayTetrahedron *tet, Feature *inc , InclusionGeometryType t, double a, double b, double c) : MatrixInclusionComposite( tet, inc, t, a, b, c )
{

    ReussMatrixInclusionComposite voigt( matrix, inclusion ) ;
    voigt.apply() ;
    fictious = voigt ;
    Matrix S = inclusion.C ;
    Matrix I = Composite::I4( S ) ;
    hint = fictious ;
    B = I - S * inclusion.volume ;
    Composite::invertTensor( B ) ;
    B *= S * matrix.volume ;
    
}

BiphasicSelfConsistentComposite::BiphasicSelfConsistentComposite( const Phase & mat, const Phase & inc) : MatrixInclusionComposite( mat, inc)
{
    matrix.volume = std::max(matrix.volume, 1e-12) ;
    inclusion.volume = std::max(inclusion.volume, 1e-12) ;
    volume = matrix.volume + inclusion.volume ;
    matrix.volume /= volume ;
    inclusion.volume /= volume ;
    volume=1;

    ReussMatrixInclusionComposite voigt( matrix, inclusion) ;
    voigt.apply() ;
    fictious = static_cast<Phase>( voigt ) ;
    Matrix S = inclusion.C ;
    Matrix I = Composite::I4( S ) ;
    hint = fictious ;
    B = I - S * inclusion.volume ;
    Composite::invertTensor( B ) ;
    B *= S * matrix.volume ;
}


BiphasicSelfConsistentComposite::BiphasicSelfConsistentComposite( const Phase & mat, const Phase &inc,  const BiphasicSelfConsistentComposite & hint) : MatrixInclusionComposite( mat, inc ),fictious(hint), hint(hint)
{
    matrix.volume = std::max(matrix.volume, 1e-12) ;
    inclusion.volume = std::max(inclusion.volume, 1e-12) ;
    volume = matrix.volume + inclusion.volume ;
    matrix.volume /= volume ;
    inclusion.volume /= volume ;
    volume=1;
    
    Matrix S = hint.C ;
    Matrix I = Composite::I4( S ) ;
    B = I - S * inclusion.volume ;
    Composite::invertTensor( B ) ;
    B *= S * matrix.volume ;

}
void BiphasicSelfConsistentComposite::getStrainConcentrationTensor()
{
    double error = 1. ;
    Matrix Sp = inclusion.C ;
    Matrix I = Composite::I4( Sp ) ;
    double nc1 = std::inner_product(&matrix.C.array()[0], &matrix.C.array()[matrix.C.array().size()], &matrix.C.array()[0], double(0)) ;
    double nc2 = std::inner_product(&inclusion.C.array()[0], &inclusion.C.array()[inclusion.C.array().size()], &inclusion.C.array()[0], double(0)) ;
    double maxnorm = sqrt(std::max(nc1, nc2)) ;
    Matrix Gp = inclusion.C ;

//     std::cout << matrix.C.array()[0] << "  " << inclusion.C.array()[0] << std::endl ;
    double minerr = 1e9 ;    
    fictious.C = inclusion.C; //MoriTanakaMatrixInclusionComposite(Phase(matrix),Phase(inclusion)).getBehaviour()->getTensor(Point());
//     fictious.C.print();
    Matrix S = fictious.C ;
    Matrix G = fictious.C ;
    int count = 0 ;
    Matrix A0 = S ;
    Matrix A1 = S ;
    Matrix del = S ;
    
    size_t sz = 1 ;
    Vector deltav(sz) ;
    for(size_t i = 0 ;  i < deltav.size() ; i++)
        deltav[i] = ((double)i+.5)/sz ;
    Vector nu_eff(deltav.size()*deltav.size()) ;
    Vector E_eff(deltav.size()*deltav.size()) ;
    
    double mv = matrix.volume ;
    double iv = inclusion.volume ;

    count = 0 ;
    double lminerror = 1e9 ;
    
    do
    {
        fictious.volume = mv ;
        inclusion.volume = iv ;
        DiluteMatrixInclusionComposite mtFictiousSecond(fictious, inclusion) ;
        mtFictiousSecond.apply() ;
        fictious.volume = iv ;
        matrix.volume = mv ;
        DiluteMatrixInclusionComposite mtFictiousFirst(fictious, matrix) ;
        mtFictiousFirst.apply() ;
        A0 = mtFictiousFirst.MatrixInclusionComposite::inclusion.A ;
        A1 = mtFictiousSecond.MatrixInclusionComposite::inclusion.A ;
        
        G = A0 * matrix.volume + A1 * inclusion.volume ;
//         G = makeStiffnessIsotropic(G) ;
        Composite::invertTensor( G ) ;
        G = A1*G ;
//         G = makeStiffnessIsotropic(G) ;
//         G = makeEshelbyIsotropic(G) ;
        S = G ;
        G = (inclusion.C - matrix.C)*G*inclusion.volume ; 
        G += matrix.C  ;
        Vector K = fictious.C.array() - G.array() ;
        error = std::inner_product(&K[0], &K[K.size()], &K[0], double(0))/maxnorm ;
        G = makeStiffnessIsotropic(G) ;
        if(error < lminerror)
        {
            lminerror = error ;
            Gp = G ;
            Sp = S ;
        }
//         G = makeEshelbyIsotropic(G) ;

        fictious.C = G;

    } while( (++count < 1024 && lminerror > 1e-32) || count < 4) ;

    minerr = lminerror ;       
    size_t idx = 0;//i*deltav.size()+j ;
    nu_eff[idx] = Gp[1][0]/(Gp[1][1]+Gp[1][0]) ;
    E_eff[idx] = Gp[1][1]*(1.+nu_eff[idx])*(1.-2.*nu_eff[idx])/(1.-nu_eff[idx]) ;


    std::cerr << matrix.volume <<"  " <<std::flush ;
    for(size_t i = 0 ;  i < E_eff.size() ; i++)
        std::cerr << E_eff[i]*32 <<"  " <<nu_eff[i] << "  "  <<std::flush ;
    std::cerr << std::endl ;

    fictious.C = Gp; 

    B = I - Sp * inclusion.volume ;
    Composite::invertTensor( B ) ;
    B *= Sp * matrix.volume ;
        
}


MatrixMultiInclusionComposite::MatrixMultiInclusionComposite( DelaunayTriangle *tri, std::vector<Feature *> inc, InclusionGeometryType t, double a, double b, double c) : Composite( tri, inc , t, a, b, c )
{
    matrix = Phase( tri, t, a, b, c  ) ;

    for( size_t i = 0 ; i < inc.size() ; i++ )
    {
        inclusions.push_back( Phase( inc[i] , tri, t, a, b, c  ) ) ;
    }
    volume = matrix.volume ;

    for( size_t i = 0 ; i < inc.size() ; i++ )
        matrix.volume -= inc[i]->area() ;

    matrix.volume = std::max(matrix.volume, 1e-12) ;
    matrix.volume /= volume ;
    matrix.volume = std::min(matrix.volume, 1. - 1e-12) ;

    for( size_t i = 0 ; i < inclusions.size() ; i++ )
    {
        grains.push_back( MoriTanakaMatrixInclusionComposite( matrix, inclusions[i]) ) ;
        grains[i].volume = inclusions[i].volume ;
    }
//     std::cout << "vm = " << matrix.volume << ", vi = " << grains[0].volume << std::endl ;
    volume = 1. ;

}

MatrixMultiInclusionComposite::MatrixMultiInclusionComposite( const Phase & p, const std::vector<Phase> & inc) : Composite( p )
{

    matrix = p ;
    volume = matrix.volume ;

    for(size_t i = 0 ; i < inc.size() ; i++)
        inclusions.push_back( inc[i]) ;

    matrix.volume = std::max(matrix.volume, 1e-12) ;
    matrix.volume /= volume ;
    matrix.volume = std::min(matrix.volume, 1. - 1e-12) ;

    for( size_t i = 0 ; i < inclusions.size() ; i++ )
        inclusions[i].volume /= volume ;

    for( size_t i = 0 ; i < inclusions.size() ; i++ )
    {
        grains.push_back( MoriTanakaMatrixInclusionComposite( matrix, inclusions[i]) ) ;
        grains[i].volume = inclusions[i].volume ;
    }

}

MatrixMultiInclusionComposite::MatrixMultiInclusionComposite( DelaunayTetrahedron *tet, std::vector<Feature *> inc, InclusionGeometryType t, double a, double b, double c ) : Composite( tet, inc, t, a, b, c )
{
    matrix = Phase( tet, t, a, b, c  ) ;

    for( size_t i = 0 ; i < inc.size() ; i++ )
    {
        if( inclusions.empty() )
            inclusions.push_back( Phase( inc[i] ) ) ;
        else
        {
            int found = -1 ;
            Phase test( inc[i] ) ;

            for( size_t j = 0 ; j < inclusions.size() ; j++ )
            {
                if( test.C == inclusions[j].C && test.beta[0] == inclusions[j].beta[0] )
                    found = ( int ) j ;
            }

            if( found == -1 )
                inclusions.push_back( test );
            else
                inclusions[found].volume += test.volume ;
        }
    }


    volume = matrix.volume ;

    for( size_t i = 0 ; i < inclusions.size() ; i++ )
        matrix.volume -= inclusions[i].volume ;

    matrix.volume = std::max(matrix.volume, 1e-6*tet->volume()) ;

    matrix.volume /= volume ;

    for( size_t i = 0 ; i < inclusions.size() ; i++ )
        inclusions[i].volume /= volume ;

    for( size_t i = 0 ; i < inclusions.size() ; i++ )
    {
        grains.push_back( MoriTanakaMatrixInclusionComposite( matrix, inclusions[i]) ) ;
        grains[i].volume /= 1. - matrix.volume ;
    }

}

void MatrixMultiInclusionComposite::apply()
{
   
    for( size_t i = 0 ; i < grains.size() ; i++ )
        grains[i].apply() ; 
    
    this->getStrainLocalizationTensor() ;

    C = Matrix( C.numRows(), C.numCols() ) ;

//    std::cerr << "-----" << std::endl ;
    for( size_t i = 0 ; i < grains.size() ; i++ )
        C += grains[i].A * ( grains[i].C * grains[i].volume ) ;

//    grains[0].C.print() ;
    //  std::cerr << grains[0].volume << std::endl ;
    // grains[0].A.print() ;

    beta.resize( beta.size(), 0. );

    for( size_t i = 0 ; i < grains.size() ; i++ )
        beta += ( grains[i].A * grains[i].volume ) * grains[i].beta ;
    
//     std::cout << beta[0] << std::endl ;
}

void MatrixMultiInclusionComposite::getStrainLocalizationTensor()
{
    for( size_t i = 0 ; i < grains.size() ; i++ )
        grains[i].A = Matrix( C.numRows(), C.numCols() ) ;
}

VoigtMatrixMultiInclusionComposite::VoigtMatrixMultiInclusionComposite( DelaunayTriangle *tri, std::vector<Feature *> inc , InclusionGeometryType t, double a, double b, double c) : MatrixMultiInclusionComposite( tri, inc, t, a, b, c )
{

}

VoigtMatrixMultiInclusionComposite::VoigtMatrixMultiInclusionComposite( DelaunayTetrahedron *tet, std::vector<Feature *> inc, InclusionGeometryType t, double a, double b, double c ) : MatrixMultiInclusionComposite( tet, inc, t, a, b, c )
{

}

VoigtMatrixMultiInclusionComposite::VoigtMatrixMultiInclusionComposite( std::vector<DelaunayTriangle *> tri , InclusionGeometryType t, double a, double b, double c) : MatrixMultiInclusionComposite( tri[0], std::vector<Feature *>( 0 ), t, a, b, c )
{
    for( size_t i = 0 ; i < tri.size() ; i++ )
        inclusions.push_back( Phase( tri[i], t, a, b, c ) ) ;

    matrix.volume = 0 ;
    volume = 0 ;

    for( size_t i = 0 ; i < inclusions.size() ; i++ )
    {
        inclusions[i].volume = tri[i]->area() ;
        volume += inclusions[i].volume ;
    }

    for( size_t i = 0 ; i < inclusions.size() ; i++ )
        inclusions[i].volume /= volume ;

    grains.clear() ;
}

VoigtMatrixMultiInclusionComposite::VoigtMatrixMultiInclusionComposite( std::vector<DelaunayTetrahedron *> tet, InclusionGeometryType t, double a, double b, double c ) : MatrixMultiInclusionComposite( tet[0], std::vector<Feature *>( 0 ), t, a, b, c )
{
    for( size_t i = 0 ; i < tet.size() ; i++ )
        inclusions.push_back( Phase( tet[i], t, a, b, c  ) ) ;

    matrix.volume = 0 ;
    volume = 0 ;

    for( size_t i = 0 ; i < inclusions.size() ; i++ )
    {
        inclusions[i].volume = tet[i]->volume() ;
        volume += inclusions[i].volume ;
    }

    for( size_t i = 0 ; i < inclusions.size() ; i++ )
        inclusions[i].volume /= volume ;

    grains.clear() ;
}

void VoigtMatrixMultiInclusionComposite::getStrainLocalizationTensor()
{
    for( size_t i = 0 ; i < grains.size() ; i++ )
        grains[i].A = Composite::I4( C ) ;
}

ReussMatrixMultiInclusionComposite::ReussMatrixMultiInclusionComposite( DelaunayTriangle *tri, std::vector<Feature *> inc, InclusionGeometryType t, double a, double b, double c ) : MatrixMultiInclusionComposite( tri, inc, t, a, b, c )
{

}

ReussMatrixMultiInclusionComposite::ReussMatrixMultiInclusionComposite( DelaunayTetrahedron *tet, std::vector<Feature *> inc, InclusionGeometryType t, double a, double b, double c ) : MatrixMultiInclusionComposite( tet, inc, t, a, b, c )
{

}

void ReussMatrixMultiInclusionComposite::getStrainLocalizationTensor()
{
    C = Matrix( matrix.C.numRows(), matrix.C.numCols() ) ;
    bool zero = false ;

    for( size_t i = 0 ; i < grains.size() ; i++ )
    {
        if(grains[i].C.array().max() > POINT_TOLERANCE)
        {
            Matrix Cg = grains[i].C ;
            Composite::invertTensor( Cg ) ;
            C += Cg * grains[i].volume ;
        }
        else
            zero = true ;

    }

    if(!zero)
    {
        Composite::invertTensor( C ) ;

        for( size_t i = 0 ; i < grains.size() ; i++ )
        {
            Matrix Cg = grains[i].C ;
            Composite::invertTensor( Cg ) ;
            grains[i].A = Cg * C ;

        }
    }
    else
    {
        std::cout << "Reuss scheme not valid for empty inclusions // fall back on Voigt scheme" << std::endl ;
        for( size_t i = 0 ; i < grains.size() ; i++ )
            grains[i].A = Composite::I4( C ) ;
    }
}

GeneralizedSelfConsistentComposite::GeneralizedSelfConsistentComposite( std::vector<DelaunayTriangle *> tri, InclusionGeometryType t, double a, double b, double c) : VoigtMatrixMultiInclusionComposite( tri, t, a, b, c )
{

}

GeneralizedSelfConsistentComposite::GeneralizedSelfConsistentComposite( std::vector<DelaunayTetrahedron *> tet, InclusionGeometryType t, double a, double b, double c ) : VoigtMatrixMultiInclusionComposite( tet, t, a, b, c )
{

}

void GeneralizedSelfConsistentComposite::apply()
{
    while( !converged() )
    {
        makeGrains() ;

        for( size_t i = 0 ; i < grains.size() ; i++ )
        {
            matrix.C += grains[i].C * grains[i].volume ;
            matrix.beta += grains[i].beta * grains[i].volume ;
        }
    }
}

void GeneralizedSelfConsistentComposite::makeGrains()
{
    previous = Phase( matrix ) ;

    if( grains.size() == 0 )
    {
        matrix.C = Matrix( inclusions[0].C ) ;

        for( size_t i = 0 ; i < inclusions.size() ; i++ )
        {
            matrix.C += inclusions[i].C * inclusions[i].volume ;
            matrix.beta += inclusions[i].beta * inclusions[i].volume ;
        }

        previous = Phase( matrix ) ;
    }

    for( size_t i = 0 ; i < inclusions.size() ; i++ )
    {
        grains.push_back( MoriTanakaMatrixInclusionComposite( matrix, inclusions[i] ) ) ;
        grains[i].inclusion.volume = inclusions[i].volume ;
        grains[i].matrix.volume = 1. - inclusions[i].volume ;
    }
}

bool GeneralizedSelfConsistentComposite::converged()
{
    if( grains.size() == 0 )
        return false ;

    Matrix epsilon = matrix.C - previous.C ;
    if(matrix.C.array().max() > POINT_TOLERANCE)
    {
        Matrix S = matrix.C ;
        Composite::invertTensor( S ) ;
        epsilon *= S ;
    }

    Vector r = epsilon.array() ;
    double rmax = std::abs( r.max() ) ;
    double rmin = std::abs( r.min() ) ;
    return std::max( rmax, rmin ) < POINT_TOLERANCE ;
}

}
