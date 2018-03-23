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

void Composite::apply()
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

Matrix  Composite::eshelbyEllipsoid(const Matrix & C, double a, double b, double c)
{
    Matrix S( C.numRows(), C.numCols() ) ;
    double nu = 0.2 ;
    double x = C[0][0] ;
    double y = C[0][1] ;
    if(std::abs(x+y) > POINT_TOLERANCE)
        nu = y / ( x + y ) ;
    else
        nu = 0 ;

    double theta = asin(sqrt((a*a-c*c)/a*a)) ;
    double k =  (a*a-b*b) / (a*a - c*c) ;

    double F = 0 ;
    double E = 0 ;
    double Ec = 0 ;
    double Fc = 0 ;
    double dw = theta*1e-1 ;
    double g0 = -1./sqrt(3) ;
    double g1 = 1./sqrt(3) ;
    for(double w = 0 ; w < theta-dw*.5 ; w += 2.*dw)
    {

        double s0 = sin(w+g0*dw+w) ;
        double d0 = 1.-k*s0*s0 ;
        if( d0 >0 )
        {
            double itemp = sqrt(d0) ;
            if(itemp > POINT_TOLERANCE)
            {
                double df = dw / itemp ;
                double y = df - Fc ;
                double t = F+y ;
                Fc = (t- F) - y ;
                F = t ;

                double de = dw * itemp ;
                y = de - Ec ;
                t = E+y ;
                Ec = (t- E) - y ;
                E = t ;

            }
        }
        double s1 = sin(w+g1*dw+w) ;
        double d1 = 1.-k*s1*s1 ;
        if( d1 >0 )
        {
            double itemp = sqrt(d1) ;
            if(itemp > POINT_TOLERANCE)
            {
                double df = dw / itemp ;
                double y = df - Fc ;
                double t = F+y ;
                Fc = (t- F) - y ;
                F = t ;

                double de = dw * itemp ;
                y = de - Ec ;
                t = E+y ;
                Ec = (t- E) - y ;
                E = t ;
            }
        }
    }

    std::cout << F << "  " << E << std::endl ;

    double I1 = 4.*M_PI*a*b*c / ((a*a-b*b)*sqrt(a*a-c*c))*(F-E) ;
    double I3 = 4.*M_PI*a*b*c / ((b*b-c*c)*sqrt(a*a-c*c))*(b*sqrt(a*a-c*c)/(a*c)-E) ;
    double I2 = 4.*M_PI-I1-I3 ;

    double I12 = (I2-I1)/(a*a-b*b) ;
    double I11 = (3.*I1 + (c*c-b*b)*I12 - c*c*4.*M_PI/(a*a))/(3*(a*a-c*c));
    double I13 = 4.*M_PI/(a*a) - 3.*I11 - I12;

    double I23 = (I3-I2)/(b*b-c*c) ;
    double I22 = (3.*I2 + (a*a-c*c)*I23 - a*a*4.*M_PI/(b*b))/(3*(b*b-a*a));
    double I21 = 4.*M_PI/(b*b) - 3.*I22 - I23;

    double I31 = (I1-I3)/(c*c-a*a) ;
    double I33 = (3.*I3 + (b*b-a*a)*I31 - b*b*4.*M_PI/(c*c))/(3*(c*c-b*b));
    double I32 = 4.*M_PI/(c*c) - 3.*I33 - I31;

//     std::cout << I1 << "  " << I2 << "  " << I3 << "  " << I11 << "  "<< I22 << "  " << I33 << std::endl ;
//     std::cout << I12 << "  " << I12 << "  " << I23 << "  " << I21 << "  "<< I23 << "  " << I32 << std::endl ;

    double Q =  3./(8.*M_PI*(1.-nu)) ;
    double R = (1.-2.*nu)/(8*M_PI*(1.-nu)) ;
    S[0][0] = Q*a*a*I11+R*I1 ;
    S[1][1] = Q*b*b*I22+R*I2 ;
    S[2][2] = Q*c*c*I33+R*I3 ;

    S[0][1] = Q*b*b*I12-R*I1 ;
    S[0][2] = Q*c*c*I13-R*I1 ;
    S[1][0] = Q*a*a*I21-R*I2 ;
    S[1][2] = Q*c*c*I23-R*I2 ;
    S[2][0] = Q*a*a*I31-R*I3 ;
    S[2][1] = Q*b*b*I32-R*I3 ;

    S[3][3] = (a*a+b*b)*.5*Q*I12+0.5*R*(I1+I2) ;
    S[4][4] = (a*a+c*c)*.5*Q*I13+0.5*R*(I1+I3) ;
    S[5][5] = (b*b+c*c)*.5*Q*I23+0.5*R*(I2+I3) ;

//     S.print() ;

    return S ;

}

Matrix  Composite::eshelbyCylinder(const Matrix & C, double a, double b)
{
    Matrix S( C.numRows(), C.numCols() ) ;
    double nu = 0.2 ;
    double x = C[0][0] ;
    double y = C[0][1] ;
    if(std::abs(x+y) > POINT_TOLERANCE)
        nu = y / ( x + y ) ;
    else
        nu = 0 ;

    S[0][0] = 1./(2.*(1.-nu))*((b*b+2.*a*b)/((a+b)*(a+b))+(1.-2.*nu)*b/(a+b)) ;
    S[1][1] = 1./(2.*(1.-nu))*((a*a+2.*a*b)/((a+b)*(a+b))+(1.-2.*nu)*a/(a+b)) ;
    S[2][2] = 0 ;

    S[0][1] = 1./(2.*(1.-nu))*(b*b/((a+b)*(a+b))-(1.-2.*nu)*b/(a+b)) ;
    S[0][2] = nu*b/((1.-nu)*(a+b)) ;
    S[1][0] = 1./(2.*(1.-nu))*(a*a/((a+b)*(a+b))-(1.-2.*nu)*a/(a+b)) ;
    S[1][2] = nu*a/((1.-nu)*(a+b)) ;
    S[2][0] = 0 ;
    S[2][1] = 0 ;

    S[3][3] = 1./(1.-nu)*((b*b+a*a)/((a+b)*(a+b))+(1.-2.*nu))*.5 ;
    S[4][4] = b/(2.*(a+b)) ;
    S[5][5] = a/(2.*(a+b)) ;
//     double Ia = 4.*M_PI*b/(a+b) ;
//     double Ib = 4.*M_PI*a/(a+b) ;
//     double Iab = 4.*M_PI/(3.*(a+b)*(a+b)) ;
//     double Iaa =  4.*M_PI/(3.*a*a)-Iab ;
//     double Ibb =  4.*M_PI/(3.*b*b)-Iab ;
//     std::cout << S[0][1] << "  " << 3./(8.*M_PI*(1.-nu))*b*b*Iab+(1.-2.*nu)/(8.*M_PI*(1.-nu))*Ia << std::endl ;
    return S ;

}

Matrix Composite::eshelby(const Matrix & C )
{
    Matrix S( C.numRows(), C.numCols() ) ;
    double nu = 0.2 ;
    double x = C[0][0] ;
    double y = C[0][1] ;
    if(x+y > POINT_TOLERANCE)
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

//         S.print() ;
//         std::cout << std::endl ;
//         eshelbyEllipsoid( C,  a,  b,  c).print() ;
//         exit(0) ;

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
        if(std::abs(m[i][i] ) > POINT_TOLERANCE)
            m[i][i] = 1. / m[i][i] ;
    }
}

MatrixInclusionComposite::MatrixInclusionComposite( DelaunayTriangle *tri, Feature *inc, InclusionGeometryType t, double a, double b, double c ) : Composite( tri, std::vector<Feature *>( 0 ), t, a, b, c )
{


    matrix = Phase( tri, t, a, b, c ) ;
    inclusion = Phase( inc , tri , t, a, b, c) ;

    matrix.volume = std::max(matrix.volume, 1e-6) ;
    inclusion.volume = std::max(inclusion.volume, 1e-6) ;
    volume = matrix.volume + inclusion.volume ;
    matrix.volume /= volume ;
    inclusion.volume /= volume ;

}

MatrixInclusionComposite::MatrixInclusionComposite( DelaunayTetrahedron *tet, Feature *inc, InclusionGeometryType t, double a, double b, double c) : Composite( tet, std::vector<Feature *>( 0 ), t, a, b, c )
{

    matrix = Phase( tet, t, a, b, c ) ;
    inclusion = Phase( inc, t, a, b, c ) ;

    matrix.volume = std::max(matrix.volume, 1e-6) ;
    inclusion.volume = std::max(inclusion.volume, 1e-6) ;
    volume = matrix.volume + inclusion.volume ;
    matrix.volume /= volume ;
    inclusion.volume /= volume ;

}

MatrixInclusionComposite::MatrixInclusionComposite( const Phase & mat, const Phase & inc) : Composite( mat), matrix(mat), inclusion(inc)
{

    matrix.volume = std::max(matrix.volume, 1e-6) ;
    inclusion.volume = std::max(inclusion.volume, 1e-6) ;
    volume = matrix.volume + inclusion.volume ;
    matrix.volume /= volume ;
    inclusion.volume /= volume ;

}


void MatrixInclusionComposite::apply()
{
    getStrainConcentrationTensor() ;
    Matrix I = Composite::I4( matrix.C ) ;

    inclusion.A = B * inclusion.volume + I * matrix.volume ;
    Composite::invertTensor( inclusion.A ) ;
    inclusion.A *= B ;
    if(matrix.volume > 1e-6)
        matrix.A = ( I - inclusion.A * inclusion.volume ) / matrix.volume ;
    else
    {
        matrix.A = ( I - inclusion.A * 1e-6 ) / 1e-6 ;
    }
//     std::cout << "vmat = "<< matrix.volume << std::endl ;

    C = ( ( matrix.A * matrix.volume ) * matrix.C ) ;
    C += ( ( inclusion.A * inclusion.volume ) * inclusion.C ) ;
    beta.resize( matrix.beta.size(), 0. );
    beta = matrix.A * matrix.volume * matrix.beta ;
    beta += inclusion.A * inclusion.volume * inclusion.beta ;

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
    B *= mori.B * matrix.volume ;
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
    if(inclusion.C.array().max() > POINT_TOLERANCE)
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
    if(matrix.C.array().max() > POINT_TOLERANCE)
    {
        Matrix I = Composite::I4( C ) ;
        Matrix S = Composite::eshelby( matrix.C ) ;
        if(matrix.t == INCLUSION_IS_ELLIPSOID)
            S = Composite::eshelbyEllipsoid( matrix.C, matrix.a, matrix.b, matrix.c ) ;
        if(matrix.t == INCLUSION_IS_CYLINDER)
            S = Composite::eshelbyCylinder( matrix.C, matrix.a, matrix.b) ;

        B = matrix.C ;
        Composite::invertTensor( B ) ;
        B *= S * ( inclusion.C - matrix.C ) ;
        B += I ;
        Composite::invertTensor( B ) ;

    }
    else
    {
        B = matrix.C ;
    }


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
    if(matrix.C.array().max() > POINT_TOLERANCE)
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
}

BiphasicSelfConsistentComposite::BiphasicSelfConsistentComposite( DelaunayTetrahedron *tet, Feature *inc , InclusionGeometryType t, double a, double b, double c) : MatrixInclusionComposite( tet, inc, t, a, b, c )
{

    ReussMatrixInclusionComposite voigt( matrix, inclusion ) ;
    voigt.apply() ;
    fictious = voigt ;
}

BiphasicSelfConsistentComposite::BiphasicSelfConsistentComposite( const Phase & mat, const Phase & inc) : MatrixInclusionComposite( mat, inc)
{

    ReussMatrixInclusionComposite voigt( matrix, inclusion) ;
    voigt.apply() ;
    fictious = static_cast<Phase>( voigt ) ;
}


BiphasicSelfConsistentComposite::BiphasicSelfConsistentComposite( const Phase & mat, const Phase &inc,  const BiphasicSelfConsistentComposite & hint) : MatrixInclusionComposite( mat, inc )
{
    fictious = hint ;

//     std::cout << "vmat & vinc = "<< matrix.volume << "  " << inclusion.volume << std::endl ;
}
void BiphasicSelfConsistentComposite::getStrainConcentrationTensor()
{
//     std::cout << "sc composite" << std::endl ;
    double error = 1. ;
    Matrix S = inclusion.C ;
    Matrix I = Composite::I4( S ) ;

    double nc1 = std::inner_product(&matrix.C.array()[0], &matrix.C.array()[matrix.C.array().size()], &matrix.C.array()[0], double(0)) ;
    double nc2 = std::inner_product(&inclusion.C.array()[0], &inclusion.C.array()[inclusion.C.array().size()], &inclusion.C.array()[0], double(0)) ;
    double maxnorm = std::max(nc1, nc2) ;

    MoriTanakaMatrixInclusionComposite mtFictiousFirst(fictious, matrix) ;
    mtFictiousFirst.apply() ;
    MoriTanakaMatrixInclusionComposite mtFictiousSecond(fictious, inclusion) ;
    mtFictiousSecond.apply() ;

    Matrix Gp = mtFictiousFirst.inclusion.A * matrix.volume ;
    Gp += mtFictiousSecond.inclusion.A * inclusion.volume ;
    Composite::invertTensor( Gp ) ;
    Gp *= mtFictiousSecond.inclusion.A ;
    Matrix Sp = Gp ;

    Gp *= inclusion.C * inclusion.volume  + matrix.C*(1.-inclusion.volume);

    Vector K = fictious.C.array() - Gp.array() ;

    error = std::inner_product(&K[0], &K[K.size()], &K[0], double(0))/maxnorm ;
//         double n = std::inner_product(&G.array()[0], &G.array()[G.array().size()], &G.array()[0], double(0)) ;
    fictious.C = Gp ;
    double perror = error ;

    do {
        perror = error ;
        MoriTanakaMatrixInclusionComposite mtFictiousFirst(fictious, matrix) ;
        mtFictiousFirst.apply() ;
        MoriTanakaMatrixInclusionComposite mtFictiousSecond(fictious, inclusion) ;
        mtFictiousSecond.apply() ;

        Matrix G = mtFictiousFirst.inclusion.A * matrix.volume ;
        G += mtFictiousSecond.inclusion.A * inclusion.volume ;
        Composite::invertTensor( G ) ;
        G *= mtFictiousSecond.inclusion.A ;
        S = G ;

        G *= inclusion.C * inclusion.volume  + matrix.C*(1.-inclusion.volume);

        Vector K = fictious.C.array() - G.array() ;
        error = std::inner_product(&K[0], &K[K.size()], &K[0], double(0))/maxnorm ;
//         double n = std::inner_product(&G.array()[0], &G.array()[G.array().size()], &G.array()[0], double(0)) ;
        fictious.C = G ;

        if(error < perror && error > 1e-16 )
        {
            Gp = G ;
            Sp = S ;
        }
    } while( error < perror && error > 1e-16 ) ;

    fictious.C = Gp ;
    B = I - Sp * inclusion.volume ;
    Composite::invertTensor( B ) ;
    B *= Sp * matrix.volume ;
}


MatrixMultiInclusionComposite::MatrixMultiInclusionComposite( DelaunayTriangle *tri, std::vector<Feature *> inc, InclusionGeometryType t, double a, double b, double c) : Composite( tri, inc , t, a, b, c )
{


    matrix = Phase( tri, t, a, b, c  ) ;

    for( size_t i = 0 ; i < inc.size() ; i++ )
    {
        if( inclusions.empty() )
            inclusions.push_back( Phase( inc[i] , tri, t, a, b, c  ) ) ;

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
//     std::cout << "padoum " << matrix.volume << "  " << inclusions.size() << std::endl ;

    for( size_t i = 0 ; i < inclusions.size() ; i++ )
        matrix.volume -= inclusions[i].volume ;

    matrix.volume = std::max(matrix.volume, tri->area()*1e-6) ;
    matrix.volume /= volume ;
    matrix.volume = std::min(matrix.volume, 1. - 1e-6) ;

    for( size_t i = 0 ; i < inclusions.size() ; i++ )
        inclusions[i].volume /= volume ;

    for( size_t i = 0 ; i < inclusions.size() ; i++ )
    {
        grains.push_back( MoriTanakaMatrixInclusionComposite( matrix, inclusions[i]) ) ;
        grains[i].volume /= 1. - matrix.volume ;
    }

}

MatrixMultiInclusionComposite::MatrixMultiInclusionComposite( const Phase & p, const std::vector<Phase> & inc) : Composite( p )
{

    matrix = p ;

    for(size_t i = 0 ; i < inc.size() ; i++)
        inclusions.push_back( inc[i]) ;

    for( size_t i = 0 ; i < inclusions.size() ; i++ )
    {
        grains.push_back( MoriTanakaMatrixInclusionComposite( matrix, inclusions[i]) ) ;
        grains[i].volume /= 1. - matrix.volume ;
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

    getStrainLocalizationTensor() ;

    C = Matrix( C.numRows(), C.numCols() ) ;

//    std::cerr << "-----" << std::endl ;
    for( size_t i = 0 ; i < grains.size() ; i++ )
        C += grains[i].A * ( grains[i].C * grains[i].volume ) ;

//    grains[0].C.print() ;
    //  std::cerr << grains[0].volume << std::endl ;
    // grains[0].A.print() ;

    beta.resize( beta.size(), 0. );

    for( size_t i = 0 ; i < grains.size() ; i++ )
        beta += ( grains[i].A * ( grains[i].volume ) ) * grains[i].beta ;
}

void MatrixMultiInclusionComposite::getStrainLocalizationTensor(InclusionGeometryType t, double a, double b, double c)
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

void VoigtMatrixMultiInclusionComposite::getStrainLocalizationTensor(InclusionGeometryType t, double a, double b, double c)
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

void ReussMatrixMultiInclusionComposite::getStrainLocalizationTensor(InclusionGeometryType t, double a, double b, double c)
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
