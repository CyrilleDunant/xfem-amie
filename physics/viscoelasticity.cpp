#include "viscoelasticity.h"
#include "material_laws/material_laws.h"
#include "../elements/generalized_spacetime_viscoelastic_element_state.h"

using namespace Amie ;

void Amie::placeMatrixInBlock(const Matrix & m, size_t i, size_t j, Matrix & ret)
{
    size_t r = m.numRows() ;
    size_t c = m.numCols() ;
    for(size_t k = i*r ; k < (i+1)*r ; k++)
    {
        for(size_t l = j*c ; l < (j+1)*c ; l++)
        {
            ret[k][l] = m[k-i*r][l-j*c] ;
        }
    }
}

void Amie::addMatrixInBlock(const Matrix & m, size_t i, size_t j, Matrix & ret)
{
    size_t r = m.numRows() ;
    size_t c = m.numCols() ;
    for(size_t k = i*r ; k < (i+1)*r ; k++)
    {
        for(size_t l = j*c ; l < (j+1)*c ; l++)
        {
            ret[k][l] += m[k-i*r][l-j*c] ;
        }
    }
}

void Amie::substractMatrixInBlock(const Matrix & m, size_t i, size_t j, Matrix & ret)
{
    size_t r = m.numRows() ;
    size_t c = m.numCols() ;
    for(size_t k = i*r ; k < (i+1)*r ; k++)
    {
        for(size_t l = j*c ; l < (j+1)*c ; l++)
        {
            ret[k][l] -= m[k-i*r][l-j*c] ;
        }
    }
}

void Amie::getBlockInMatrix( const Matrix & source, size_t i, size_t j, Matrix & ret)
{
    size_t r = ret.numRows() ;
    size_t c = ret.numCols() ;
    for(size_t k = i*r ; k < (i+1)*r ; k++)
    {
        for(size_t l = j*c ; l < (j+1)*c ; l++)
        {
            ret[k-i*r][l-j*c] = source[k][l] ;
        }
    }
}

void Viscoelasticity::makeBlockConnectivity()
{
    connectivity.clear() ;
    connectivityViscous.clear() ;

    switch(model)
    {
    case PURE_ELASTICITY:
        connectivity.push_back( BlockConnectivity(0,0) ) ;
        break ;
    case PURE_VISCOSITY:
        connectivityViscous.push_back( BlockConnectivity(0,0) ) ;
        break ;
    case KELVIN_VOIGT:
        connectivity.push_back( BlockConnectivity(0,0) ) ;
        connectivityViscous.push_back( BlockConnectivity(0,0) ) ;
        break ;
    case MAXWELL:
   {
        BlockConnectivity elastic(0,0) ;
        elastic.add(1,1) ;
        elastic.substract(0,1) ;
        elastic.substract(1,0) ;
        connectivity.push_back( elastic ) ;
        connectivityViscous.push_back( BlockConnectivity(1,1) ) ;
        break ;
    }
    case BURGER:
    {
        BlockConnectivity elastic(0,0) ;
        elastic.add(1,2) ;
        elastic.add(2,1) ;
        elastic.substract(0,1) ;
        elastic.substract(0,2) ;
        elastic.substract(1,0) ;
        elastic.substract(2,0) ;
        connectivity.push_back( elastic ) ;
        connectivity.push_back( BlockConnectivity(1,1) ) ;
        connectivity.push_back( BlockConnectivity(2,2) ) ;
        connectivityViscous.push_back( BlockConnectivity(1,1) ) ;
        connectivityViscous.push_back( BlockConnectivity(2,2) ) ;
        break ;
    }
    case GENERALIZED_MAXWELL:
        connectivity.push_back( BlockConnectivity(0,0 ) ) ;
        for(int i = 1 ; i < blocks ; i++)
        {
            BlockConnectivity branch(i,i) ;
            branch.substract(0,i) ;
            branch.substract(i,0) ;
            connectivity.push_back(branch) ;
            connectivityViscous.push_back( BlockConnectivity(i,i) ) ;
        }
        break ;
    case GENERALIZED_KELVIN_VOIGT:
    {
        BlockConnectivity elastic(0,0) ;
        for(int i = 1 ; i < blocks ; i++)
        {
            elastic.substract(0,i) ;
            elastic.substract(i,0) ;
            for(int j = i+1 ; j < blocks ; j++)
            {
                elastic.add(i,j) ;
                elastic.add(j,i) ;
            }
            connectivity.push_back( BlockConnectivity(i,i) ) ;
            connectivityViscous.push_back( BlockConnectivity(i,i) ) ;
        }
        connectivity.push_back( elastic) ;
        break ;
    }
    default:
        break ;
    }
} 

Viscoelasticity::Viscoelasticity(ViscoelasticModel m, const Matrix & rig, int n, double r) : LinearForm(rig, false, false, (1+n)*(rig.numRows()/3+1)), model(m), blocks(1+n), effblocks(1)
{
    v.push_back(XI);
    v.push_back(ETA);
    if(rig.size() > 9)
        v.push_back(ZETA);
    v.push_back(TIME_VARIABLE);

    param.resize(rig.numRows()*(1+n), rig.numCols()*(1+n)) ;
    eta.resize(rig.numRows()*(1+n), rig.numCols()*(1+n)) ;
    param = 0 ;
    eta = 0 ;


    switch(model)
    {
    case PURE_ELASTICITY:
        placeMatrixInBlock( rig, 0,0, param) ;
        break ;
    case PURE_VISCOSITY:
        tensors.push_back(rig*0) ;
        placeMatrixInBlock( rig, 0,0, eta) ;
        break ;
    default:
        std::cout << "warning: wrong constructor for Viscoelasticity" << std::endl ;
    }

    tensors.push_back(rig) ;

    makeBlockConnectivity() ;

}

Viscoelasticity::Viscoelasticity(ViscoelasticModel m, const Matrix & rig, const Matrix & e, int b, int n, double r) : LinearForm(rig, false, false, (1+n+b+(m == MAXWELL))*(rig.numRows()/3+1)), model(m), blocks(1+n+b+(m==MAXWELL)), effblocks(1+(m==MAXWELL))
{
    v.push_back(XI);
    v.push_back(ETA);
    if(rig.size() > 9)
        v.push_back(ZETA);
    v.push_back(TIME_VARIABLE);

    param.resize(rig.numRows()*(1+n+b+(m == MAXWELL)), rig.numCols()*(1+n+b+(m == MAXWELL))) ;
    eta.resize(rig.numRows()*(1+n+b+(m == MAXWELL)), rig.numCols()*(1+n+b+(m == MAXWELL))) ;
    param.array() = 0 ;
    eta.array() = 0 ;

    switch(model)
    {
    case KELVIN_VOIGT:
        placeMatrixInBlock( rig, 0,0, param) ;
        placeMatrixInBlock( e, 0,0, eta) ;
        break ;
    case MAXWELL:
    {
        Matrix r0 = rig*(-1.) ;
        placeMatrixInBlock( rig, 0,0, param) ;
        placeMatrixInBlock( r0,  0,1+b, param) ;
        placeMatrixInBlock( r0,  1+b,0, param) ;
        placeMatrixInBlock( rig, 1+b,1+b, param) ;
        placeMatrixInBlock( e, 1+b,1+b, eta) ;

        break ;
    }
    default:
        std::cout << "warning: wrong constructor for Viscoelasticity" << std::endl ;
    }

    tensors.push_back(rig) ;
    tensors.push_back(e) ;

    makeBlockConnectivity() ;
}

Viscoelasticity::Viscoelasticity(const Matrix & rig, const Matrix & e, int b, int n, double r) : LinearForm(rig, false, false, (n+b)*((rig.numRows()/b)/3+1)), model(GENERAL_VISCOELASTICITY), blocks(n+b),effblocks(b)
{
    v.push_back(XI);
    v.push_back(ETA);
    if(rig.numCols()/(n+b) > 3)
        v.push_back(ZETA);
    v.push_back(TIME_VARIABLE);

    param.resize(rig.numRows()/b*blocks, rig.numCols()/b*blocks) ;
    eta.resize(rig.numRows()/b*blocks, rig.numCols()/b*blocks) ;
    param.array() = 0 ;
    eta.array() = 0 ;

    placeMatrixInBlock( rig, 0,0, param) ;
    placeMatrixInBlock( e, 0,0, eta) ;

    Matrix buffer( 3+3*(v.size() == 4), 3+3*(v.size() == 4)) ;
    for(int i = 0 ; i < effblocks ; i++)
    {
        for(int j = 0 ; j < effblocks ; j++)
        {
             getBlockInMatrix( param, i,j, buffer) ;
             tensors.push_back( buffer ) ;
             getBlockInMatrix( eta, i,j, buffer) ;
             tensors.push_back( buffer ) ;
        }
    }

    makeBlockConnectivity() ;

}

Viscoelasticity::Viscoelasticity(ViscoelasticModel m, const Matrix & c_kv, const Matrix & e_kv, const Matrix & c_mx, const Matrix & e_mx, int b, int n, double r) : LinearForm(c_kv, false, false, (3+n+b)*(c_kv.numRows()/3+1)), model(m), blocks(3+b+n),effblocks(3)
{
    v.push_back(XI);
    v.push_back(ETA);
    if(c_kv.size() > 9)
        v.push_back(ZETA);
    v.push_back(TIME_VARIABLE);

    param.resize(c_kv.numRows()*(3+n+b), c_kv.numCols()*(3+n+b)) ;
    eta.resize(c_kv.numRows()*(3+n+b), c_kv.numCols()*(3+n+b)) ;
    param.array() = 0 ;
    eta.array() = 0 ;

    switch(model)
    {
    case BURGER:
    {
        placeMatrixInBlock( c_mx, 0,0, param) ;
        placeMatrixInBlock( c_mx, 1+b,1+b, param) ;
        placeMatrixInBlock( c_mx, 2+b,2+b, param) ;
        placeMatrixInBlock( c_mx, 1+b,2+b, param) ;
        placeMatrixInBlock( c_mx, 2+b,1+b, param) ;
        Matrix r_mx = c_mx * (-1.) ;
        placeMatrixInBlock( r_mx, 0,1+b, param) ;
        placeMatrixInBlock( r_mx, 0,2+b, param) ;
        placeMatrixInBlock( r_mx, 1+b,0, param) ;
        placeMatrixInBlock( r_mx, 2+b,0, param) ;
        addMatrixInBlock( c_kv, 2+b,2+b, param) ;


        placeMatrixInBlock( e_mx, 1+b,1+b, eta) ;
        placeMatrixInBlock( e_kv, 2+b,2+b, eta) ;

// 			param.print() ;
// 			eta.print() ;
// 			exit(0) ;
        break ;
    }
    default:
        std::cout << "warning: wrong constructor for Viscoelasticity" << std::endl ;
    }

    tensors.push_back(c_mx) ;
    tensors.push_back(e_mx) ;
    tensors.push_back(c_kv) ;
    tensors.push_back(e_kv) ;

    makeBlockConnectivity() ;

}

Viscoelasticity::Viscoelasticity(ViscoelasticModel m, const Matrix & c0, std::vector<std::pair<Matrix, Matrix> > & branches, int b, int n, double r) : LinearForm(c0, false, false, (1+n+b+branches.size())*(c0.numRows()/3+1)), model(m), blocks(1+n+b+branches.size()), effblocks(1+branches.size())
{
    v.push_back(XI);
    v.push_back(ETA);
    if(c0.size() > 9)
        v.push_back(ZETA);
    v.push_back(TIME_VARIABLE);

    param.resize(c0.numRows()*(1+n+b+branches.size()), c0.numCols()*(1+n+b+branches.size())) ;
    eta.resize(c0.numRows()*(1+n+b+branches.size()), c0.numCols()*(1+n+b+branches.size())) ;
    param.array() = 0 ;
    eta.array() = 0 ;

    switch(model)
    {
    case GENERALIZED_KELVIN_VOIGT:
    {
        placeMatrixInBlock( c0, 0,0, param) ;
        Matrix r0 = c0*(-1) ;
        for(size_t i = 1 ; i < branches.size()+1 ; i++)
        {
            placeMatrixInBlock( r0, i+b,0, param) ;
            placeMatrixInBlock( r0, 0,i+b, param) ;
            placeMatrixInBlock( c0, i+b,i+b, param) ;
            addMatrixInBlock( branches[i-1].first, i+b,i+b, param) ;
            for(size_t j = i+1 ; j < branches.size()+1 ; j++)
            {
                placeMatrixInBlock( c0, i+b,j+b, param) ;
                placeMatrixInBlock( c0, j+b,i+b, param) ;
            }
            placeMatrixInBlock( branches[i-1].second, i+b,i+b, eta) ;
        }
        break ;
    }
    case GENERALIZED_MAXWELL:
    {
        placeMatrixInBlock( c0, 0,0, param) ;
        for(size_t i = 1 ; i < branches.size()+1 ; i++)
        {
            Matrix ri = branches[i-1].first * (-1) ;
            addMatrixInBlock( branches[i-1].first, 0,0, param) ;
            placeMatrixInBlock( branches[i-1].first, i+b,i+b, param) ;
            placeMatrixInBlock( ri, 0,i+b, param) ;
            placeMatrixInBlock( ri, i+b,0, param) ;
            placeMatrixInBlock( branches[i-1].second, i+b,i+b, eta) ;
        }
        break ;
    }
    default:
        std::cout << "warning: wrong constructor for Viscoelasticity" << std::endl ;
    }

    tensors.push_back(c0) ;
    for(size_t i = 0 ; i < branches.size() ; i++)
    {
        tensors.push_back(branches[i].first) ;
        tensors.push_back(branches[i].second) ;
    }

    makeBlockConnectivity() ;
}

Viscoelasticity::Viscoelasticity(ViscoelasticModel m, const Matrix & c0, const Matrix & c1, const Matrix & e1, int b, int n, double r) : LinearForm(c0, false, false, (2+n+b)*(c0.numRows()/3+1)), model(m), blocks(2+n+b), effblocks(2)
{
    v.push_back(XI);
    v.push_back(ETA);
    if(c0.size() > 9)
        v.push_back(ZETA);
    v.push_back(TIME_VARIABLE);

    param.resize(c0.numRows()*(2+n+b), c0.numCols()*(2+n+b)) ;
    eta.resize(c0.numRows()*(2+n+b), c0.numCols()*(2+n+b)) ;
    param.array() = 0 ;
    eta.array() = 0 ;

    switch(model)
    {
    case GENERALIZED_KELVIN_VOIGT:
    {
        placeMatrixInBlock( c0, 0,0, param) ;
        placeMatrixInBlock( c0, 1+b,1+b, param) ;
        addMatrixInBlock( c1, 1+b,1+b, param) ;
        Matrix r0 = c0*(-1) ;
        placeMatrixInBlock( r0, 1+b,0, param) ;
        placeMatrixInBlock( r0, 0,1+b, param) ;
        addMatrixInBlock( e1, 1+b,1+b, eta) ;
        break ;
    }
    case GENERALIZED_MAXWELL:
    {
        placeMatrixInBlock( c0, 0,0, param) ;
        addMatrixInBlock( c1, 0,0, param) ;
        placeMatrixInBlock( c1, 1+b,1+b, param) ;
        Matrix r1 = c1*(-1) ;
        placeMatrixInBlock( r1, 1+b,0, param) ;
        placeMatrixInBlock( r1, 0,1+b, param) ;
        placeMatrixInBlock( e1, 1+b,1+b, eta) ;
        break ;
    }
    default:
        std::cout << "warning: wrong constructor for Viscoelasticity" << std::endl ;
    }

    tensors.push_back(c0) ;
    tensors.push_back(c1) ;
    tensors.push_back(e1) ;

    makeBlockConnectivity() ;
}

Viscoelasticity::~Viscoelasticity() {}

ElementState * Viscoelasticity::createElementState( IntegrableEntity * e)
{
    return new GeneralizedSpaceTimeViscoElasticElementState(e) ;
}

void Viscoelasticity::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
    Matrix deltaParam = getTensor( Point(0,0,0,-1) ) - getTensor( Point(0,0,0,1) ) ;
    if( std::abs(deltaParam.array().max()) > POINT_TOLERANCE || std::abs(deltaParam.array().min()) > POINT_TOLERANCE )
    {
        Matrix tmpRet(ret.numRows()/blocks, ret.numCols()/blocks) ;
        std::vector<Matrix> mat(Jinv.size(), Matrix( tensors[0].numCols(), tensors[0].numRows() ) ) ;
        std::vector<Matrix> dmat(Jinv.size(), Matrix( tensors[0].numCols(), tensors[0].numRows() ) ) ;
        std::vector<std::pair<Matrix, Matrix> > paramt( Jinv.size(), std::make_pair( param, param*0 ) ) ; 
        getTensorDotAtGaussPoints( gp, Jinv, paramt, true ) ;

        for(size_t i = 0 ; i < connectivity.size() ; i++)
        {
            tmpRet = 0 ;
            for(size_t g = 0; g < Jinv.size() ; g++)
            {
                getBlockInMatrix( paramt[g].first, connectivity[i].xplus[0], connectivity[i].yplus[0], mat[g] ) ;
                getBlockInMatrix( paramt[g].second, connectivity[i].xplus[0], connectivity[i].yplus[0], dmat[g] ) ;
            }
            vm->ieval( Differential(TIME_VARIABLE)[ Gradient(p_i) * mat * Gradient(p_j, true) ], dmat, gp, Jinv,v, tmpRet) ;
            for(size_t j = 0 ; j < connectivity[i].xplus.size() ; j++)
                placeMatrixInBlock( tmpRet, connectivity[i].xplus[j], connectivity[i].yplus[j], ret ) ;
            for(size_t j = 0 ; j < connectivity[i].xminus.size() ; j++)
                placeMatrixInBlock( -tmpRet , connectivity[i].xminus[j], connectivity[i].yminus[j], ret ) ;
        }

    }
    else
    {
        Matrix tmpRet(ret.numRows()/blocks, ret.numCols()/blocks) ;
        Matrix tmpParam(param.numRows()/blocks, param.numCols()/blocks) ;
        Matrix realParam = getTensor( Point(0,0,0,1) ) ;

        for(size_t i = 0 ; i < connectivity.size() ; i++)
        {
            tmpRet = 0 ;
            getBlockInMatrix(realParam, connectivity[i].xplus[0], connectivity[i].yplus[0], tmpParam) ;
            vm->ieval( Differential(TIME_VARIABLE)[ Gradient(p_i) * tmpParam * Gradient(p_j, true) ],  gp, Jinv,v, tmpRet) ;
            for(size_t j = 0 ; j < connectivity[i].xplus.size() ; j++)
                placeMatrixInBlock( tmpRet, connectivity[i].xplus[j], connectivity[i].yplus[j], ret ) ;
            for(size_t j = 0 ; j < connectivity[i].xminus.size() ; j++)
                placeMatrixInBlock( -tmpRet , connectivity[i].xminus[j], connectivity[i].yminus[j], ret ) ;
        }
    }

}

void Viscoelasticity::applyViscous(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
    Matrix deltaParam = getViscousTensor( Point(0,0,0,-1) ) - getViscousTensor( Point(0,0,0,1) ) ;
    if( std::abs(deltaParam.array().max()) > POINT_TOLERANCE || std::abs(deltaParam.array().min()) > POINT_TOLERANCE )
    {
        Matrix tmpRet(ret.numRows()/blocks, ret.numCols()/blocks) ;
        std::vector<Matrix> mat(Jinv.size(), Matrix( tensors[0].numCols(), tensors[0].numRows() ) ) ;
        std::vector<Matrix> dmat(Jinv.size(), Matrix( tensors[0].numCols(), tensors[0].numRows() ) ) ;
        std::vector<std::pair<Matrix, Matrix> > paramt( Jinv.size(), std::make_pair( eta, eta*0 ) ) ; 
        getViscousTensorDotAtGaussPoints( gp, Jinv, paramt, true ) ;

        for(size_t i = 0 ; i < connectivityViscous.size() ; i++)
        {
            tmpRet = 0 ;
            for(size_t g = 0; g < Jinv.size() ; g++)
            {
                getBlockInMatrix( paramt[g].first, connectivityViscous[i].xplus[0], connectivityViscous[i].yplus[0], mat[g] ) ;
                getBlockInMatrix( paramt[g].second, connectivityViscous[i].xplus[0], connectivityViscous[i].yplus[0], dmat[g] ) ;
            }
            vm->ieval( Differential(TIME_VARIABLE)[ Gradient(p_i) * mat * GradientDot(p_j, true) ], dmat, gp, Jinv,v, tmpRet) ;
            for(size_t j = 0 ; j < connectivityViscous[i].xplus.size() ; j++)
                placeMatrixInBlock( tmpRet, connectivityViscous[i].xplus[j], connectivityViscous[i].yplus[j], ret ) ;
            for(size_t j = 0 ; j < connectivityViscous[i].xminus.size() ; j++)
                placeMatrixInBlock( -tmpRet , connectivityViscous[i].xminus[j], connectivityViscous[i].yminus[j], ret ) ;
        }

    }
    else
    {
        Matrix tmpRet(ret.numRows()/blocks, ret.numCols()/blocks) ;
        Matrix tmpEta(eta.numRows()/blocks, eta.numCols()/blocks) ;
        Matrix realeta = getViscousTensor( Point(0,0,0,1) ) ;

        for(size_t i = 0 ; i < connectivityViscous.size() ; i++)
        {
            tmpRet = 0 ;
            getBlockInMatrix(realeta, connectivityViscous[i].xplus[0], connectivityViscous[i].yplus[0], tmpEta) ;
            vm->ieval( Differential(TIME_VARIABLE)[ Gradient(p_i) * tmpEta * GradientDot(p_j, true) ],  gp, Jinv,v, tmpRet) ;
            for(size_t j = 0 ; j < connectivityViscous[i].xplus.size() ; j++)
                placeMatrixInBlock( tmpRet, connectivityViscous[i].xplus[j], connectivityViscous[i].yplus[j], ret ) ;
            for(size_t j = 0 ; j < connectivityViscous[i].xminus.size() ; j++)
                placeMatrixInBlock( -tmpRet , connectivityViscous[i].xminus[j], connectivityViscous[i].yminus[j], ret ) ;
        }
    }

}

void Viscoelasticity::setBlocks( int maxBlocks )
{
    if(blocks == maxBlocks)
        return ;

    blocks = maxBlocks ;

    int numCols = 3+3*(v.size() == 4) ;
    num_dof = maxBlocks * (2 + 1*(v.size() == 4)) ;
    param.resize( numCols * maxBlocks, numCols * maxBlocks ) ;
    eta.resize( numCols * maxBlocks, numCols * maxBlocks ) ;
    param = 0 ;
    eta = 0 ;
    Matrix t = tensors[0]*(-1.) ;

    switch(model)
    {
        case PURE_ELASTICITY:
            placeMatrixInBlock( tensors[0], 0,0, param ) ;
            break ;
        case PURE_VISCOSITY:
            placeMatrixInBlock( tensors[1], 0,0, eta ) ;
            break ;
        case KELVIN_VOIGT:
            placeMatrixInBlock( tensors[0], 0,0, param ) ;
            placeMatrixInBlock( tensors[1], 0,0, eta ) ;
            break ;
        case MAXWELL:
            placeMatrixInBlock( tensors[0], 0,0, param ) ;
            placeMatrixInBlock( t, 0,1, param ) ;
            placeMatrixInBlock( t, 1,0, param ) ;
            placeMatrixInBlock( tensors[0], 1,1, param ) ;
            placeMatrixInBlock( tensors[1], 1,1, eta ) ;
            break ;
        case BURGER:
            placeMatrixInBlock( tensors[0], 0,0, param ) ;
            placeMatrixInBlock( t, 0,1, param ) ;
            placeMatrixInBlock( t, 1,0, param ) ;
            placeMatrixInBlock( t, 0,2, param ) ;
            placeMatrixInBlock( t, 2,0, param ) ;
            placeMatrixInBlock( tensors[0], 1,1, param ) ;
            placeMatrixInBlock( tensors[0], 1,2, param ) ;
            placeMatrixInBlock( tensors[0], 2,1, param ) ;
            placeMatrixInBlock( tensors[0], 2,2, param ) ;
            addMatrixInBlock( tensors[2], 2,2, param ) ;
            placeMatrixInBlock( tensors[1], 1,1, eta ) ;
            placeMatrixInBlock( tensors[3], 2,2, eta ) ;
            break ;
        case GENERALIZED_MAXWELL:
            placeMatrixInBlock( tensors[0], 0,0, param) ;
            for(size_t i = 1 ; i < tensors.size() ; i+=2)
            {
                Matrix ti = tensors[i] * (-1) ;
                int b = (i+1)/2 ;
                addMatrixInBlock( tensors[i], 0,0, param) ;
                placeMatrixInBlock( tensors[i], b,b, param) ;
                placeMatrixInBlock( ti, 0,b, param) ;
                placeMatrixInBlock( ti, b,0, param) ;
                placeMatrixInBlock( tensors[i+1], b,b, eta) ;
            }
            break ;
        case GENERALIZED_KELVIN_VOIGT:
        {
            placeMatrixInBlock( tensors[0], 0,0, param) ;
            for(size_t i = 1 ; i < tensors.size() ; i+=2)
            {
                int b = (i+1)/2 ;
                placeMatrixInBlock( t, b,0, param) ;
                placeMatrixInBlock( t, 0,b, param) ;
                placeMatrixInBlock( tensors[0], b,b, param) ;
                addMatrixInBlock( tensors[i], b,b, param) ;
                for(size_t j = i+2 ; j < tensors.size()+1 ; j+=2)
                {
                    int c = (j+1)/2 ;
                    placeMatrixInBlock( tensors[0], b,c, param) ;
                    placeMatrixInBlock( tensors[0], c,b, param) ;
                }
                placeMatrixInBlock( tensors[i+1], b,b, eta) ;
            }
            break ;
        }
        case GENERAL_VISCOELASTICITY:
        {
            for(int i = 0 ; i < effblocks ; i++)
            {
                for(int j = 0 ; j < effblocks ; j++)
                {
                    placeMatrixInBlock( tensors[ i*effblocks + j ], i,j, param) ;
                }
            }
            break ;
        }

    }
}

bool Viscoelasticity::fractured() const
{
    return false ;
}

bool Viscoelasticity::changed() const
{
    return false ;
}

Form * Viscoelasticity::getCopy() const
{
    Matrix rig = param ;
    Matrix e = eta ;
    switch(model)
    {
    case PURE_ELASTICITY:
        return new Viscoelasticity( PURE_ELASTICITY, tensors[0], blocks-effblocks, rho) ;
    case PURE_VISCOSITY:
        return new Viscoelasticity( PURE_VISCOSITY, tensors[1], blocks-effblocks, rho) ;
    case MAXWELL:
        return new Viscoelasticity( MAXWELL, tensors[0], tensors[1],0, blocks-effblocks, rho) ;
    case KELVIN_VOIGT:
        return new Viscoelasticity( KELVIN_VOIGT, tensors[0], tensors[1],0, blocks-effblocks, rho) ;
    case BURGER:
        return new Viscoelasticity( BURGER, tensors[2], tensors[3], tensors[0], tensors[1],0, blocks-effblocks, rho) ;
    case GENERALIZED_MAXWELL:
    {
        std::vector<std::pair<Matrix, Matrix> > branches ;
        for(size_t i = 1 ; i < tensors.size() ; i+= 2)
        {
            branches.push_back( std::make_pair(tensors[i], tensors[i+1]) ) ;
        }
        return new Viscoelasticity( GENERALIZED_MAXWELL, tensors[0], branches,0, blocks-effblocks, rho) ;
    }
    case GENERALIZED_KELVIN_VOIGT:
    {
        std::vector<std::pair<Matrix, Matrix> > branches ;
        for(size_t i = 1 ; i < tensors.size() ; i+= 2)
        {
            branches.push_back( std::make_pair(tensors[i], tensors[i+1]) ) ;
        }
        return new Viscoelasticity( GENERALIZED_KELVIN_VOIGT, tensors[0], branches,0, blocks-effblocks, rho) ;
    }
    case GENERAL_VISCOELASTICITY:
        return new Viscoelasticity( rig, e, blocks) ;
    }

    Viscoelasticity * copy = new Viscoelasticity(rig, e, blocks) ;
    copy->model = model ;

    return copy ;
}

Vector Viscoelasticity::getForcesFromAppliedStress( const Vector & data, Function & shape, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<Variable> & v, bool isVolumic, const Vector & normal)
{
    if(isVolumic)
    {
        return VirtualMachine().ieval(GradientDot( shape ) * data, gp, Jinv, v) ;
    }
    return data * VirtualMachine().ieval(Differential( shape, TIME_VARIABLE ), gp, Jinv, v) ;
}

Vector Viscoelasticity::getForcesFromAppliedStress( const Function & data, size_t index, size_t externaldofs,  Function & shape, IntegrableEntity * e,const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<Variable> & v, bool isVolumic, const Vector & normal)
{
    VirtualMachine vm ;

    size_t n = e->getBoundingPoints().size() ;
    Vector field(0., n*externaldofs) ;
    for(size_t i = 0 ; i < n ; i++)
        field[ i*externaldofs + index ] = vm.eval( data, e->getBoundingPoint(i) ) ;

    std::vector<Vector> g(e->getGaussPoints().gaussPoints.size(), Vector(0., externaldofs)) ;
    e->getState().getExternalFieldAtGaussPoints( field, externaldofs, g) ;

    Vector f = vm.ieval( GradientDot( shape ) * g, e, v) ;

    field = 0 ;
    for(size_t i = 0 ; i < n ; i++)
        field[ i*externaldofs + index ] = vm.deval( data, TIME_VARIABLE, e->getBoundingPoint(i) ) ;

    e->getState().getExternalFieldAtGaussPoints( field, externaldofs, g) ;

    Vector f1 = vm.ieval( Gradient( shape ) * g, e, v) ;
//	std::cout << f1[0] << ";" ;
    f += f1 ;

    return f*0.5 ;
}

Vector Viscoelasticity::getImposedStress(const Point & p, IntegrableEntity * e, int g) const
{
    return Vector(double(0), getTensor(p, e).numCols()/blocks) ;
}

Vector Viscoelasticity::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const
{
    return Vector(double(0), getTensor(p).numCols()/blocks) ;
}


void Viscoelasticity::print() const
{
    std::cout << "I am a viscoelastic model" ;
    switch(model)
    {
    case PURE_ELASTICITY:
        std::cout << " (elastic only)" << std::endl ;
        return ;
    case PURE_VISCOSITY:
        std::cout << " (viscous only)" << std::endl ;
        return ;
    case KELVIN_VOIGT:
        std::cout << " (kelvin-voigt)" << std::endl ;
        return ;
    case GENERALIZED_KELVIN_VOIGT:
        std::cout << " ( generalized kelvin-voigt)" << std::endl ;
        return ;
    case MAXWELL:
        std::cout << " (maxwell)" << std::endl ;
        return ;
    case GENERALIZED_MAXWELL:
        std::cout << " ( generalized maxwell)" << std::endl ;
        return ;
    case BURGER:
        std::cout << " (burger)" << std::endl ;
        return ;
    case GENERAL_VISCOELASTICITY:
        std::cout << " (general)" << std::endl ;
        return ;
    }
}


double ViscoelasticKelvinVoigtChainGenerator::getLCoefficient( CreepComplianceModel model, double tau, std::map<std::string, double> args) 
{
    switch(model)
    {
    case LOGPOWER_CREEP:
    {
        if(args.find("n") != args.end())
        {
            double n = args["n"] ;
            double tau_n = std::pow(tau, n) ;
            return n*tau_n/(1+tau_n) ;
        }
        return exp(-1./tau) ;
    }
    case ACI_CREEP:
    {
        double a = 10 ;
        double n = 0.6 ; 
        if(args.find("a") != args.end())
            a = args["a"] ;
        if(args.find("n") != args.end())
            n = args["n"] ;
        double tau_n = std::pow(tau, n) ;
        return a * n * tau_n / ((a + tau_n)*(a + tau_n)) ;
    }
    case CEB_CREEP:
    {
        double n = 0.3 ;
        double beta = 775 ;
        if(args.find("beta") != args.end())
            beta = args["beta"] ;
        if(args.find("n") != args.end())
            n = args["n"] ;
        return tau * std::pow( tau/(tau+beta), n-1 )*beta/((beta+tau)*(beta+tau)) ;
    }
    case B3_DRYING_CREEP:
    {
        double tau_shrinkage = 0.1 ;
        double tau_cure = 28 ;
        double tau_loading = 28 ;
        double h = 1 ;
        if(args.find("tau_shrinkage") != args.end())
            tau_shrinkage = args["tau_shrinkage"] ;
        if(args.find("tau_cure") != args.end())
            tau_cure = args["tau_cure"] ;
        if(args.find("tau_loading") != args.end())
            tau_loading = args["tau_loading"] ;
        if(args.find("h") != args.end())
            h = args["h"] ;
        double b = 8*(1-h) ;
        double xhi_cure = (tau_cure - tau_loading)/tau_shrinkage ;
        if(xhi_cure < 0)
            xhi_cure = 0 ;
        double sqrttau = sqrt( tau - xhi_cure ) ;
        double expbtanhtau = exp(b*tanh(sqrttau)) ;
        double expbtanhxhi = exp(b*tanh(sqrt(-xhi_cure))) ;
        return tau*b*expbtanhtau / (4*sqrt( sqrttau*sqrttau*( expbtanhtau - expbtanhxhi)*cosh(sqrttau)*cosh(sqrttau) )) ;
    }
    case JSCE_CREEP:
    {
        double c = -0.09 ;
        double n = 0.6 ;
        if(args.find("c") != args.end())
            c = args["c"] ;
        if(args.find("n") != args.end())
            n = args["n"] ;
        double tau_n = std::pow(tau, n) ;
        return (-exp(c*tau_n))*c*tau_n ;
    }
    case FIB_CREEP:
    {
        double tau_loading = 28 ;
        if(args.find("tau_loading") != args.end())
            tau_loading = args["tau_loading"] ;
        double c = 0.035+30./tau_loading ;
        c *= c ;
        return exp(-1./(c*tau)) ;
    }
    }
    return -1 ;
}


std::vector< std::pair<Matrix, Matrix> > ViscoelasticKelvinVoigtChainGenerator::getKelvinVoigtChain( double E_creep, double nu_creep, CreepComplianceModel model, std::string args, double tau0, int b, SpaceDimensionality dim, planeType pt) 
{
    std::vector< std::pair<Matrix, Matrix> > branches ;
    double tau = tau0 ;
    std::map<std::string, double> val = parseDefaultValues(args, ',') ;
    for(int i = 0 ; i < b ; i++)
    {
        double E = E_creep/log(10) ;
        double L = ViscoelasticKelvinVoigtChainGenerator::getLCoefficient(model, tau, val) ;
        if(model != LOGPOWER_CREEP && model != FIB_CREEP)
        {
            double tauprev = tau*0.51442831563 ;
            double taunext = tau*1.94390543758 ;
            L = (ViscoelasticKelvinVoigtChainGenerator::getLCoefficient(model, tauprev, val) + ViscoelasticKelvinVoigtChainGenerator::getLCoefficient(model, taunext, val))/2 ;
        }
        if(L > POINT_TOLERANCE)
        {
            E /= L ;
            Matrix C = Tensor::cauchyGreen( E, nu_creep, true, dim, pt) ;
            Matrix eta = C*tau ;
            branches.push_back(std::make_pair(C, eta)) ;
        }
        tau *= 10 ;
    }
    return branches ;
}


