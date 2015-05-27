#include "viscoelasticity_and_fracture.h"
#include "viscoelasticity.h"
#include "fracturecriteria/maxstrain.h"
#include "../elements/generalized_spacetime_viscoelastic_element_state.h"

using namespace Amie ;

ViscoelasticityAndFracture::ViscoelasticityAndFracture(ViscoelasticModel m, const Matrix & rig, FractureCriterion * c, DamageModel * d, int n, double r) : Viscoelasticity(m, rig, n, r), dfunc(d), criterion(c)
{
    setElasticAndViscousStiffnessMatrix() ;
}

ViscoelasticityAndFracture::ViscoelasticityAndFracture(ViscoelasticModel m, const Matrix & rig, const Matrix & e, FractureCriterion * c, DamageModel * d, int b, int n, double r) : Viscoelasticity(m, rig, e, b, n, r), dfunc(d), criterion(c)
{
    setElasticAndViscousStiffnessMatrix() ;
}

ViscoelasticityAndFracture::ViscoelasticityAndFracture(const Matrix & rig, const Matrix & e, int b, FractureCriterion * c, DamageModel * d, int n, double r) : Viscoelasticity(rig, e, b, n, r), dfunc(d) , criterion(c)
{
    setElasticAndViscousStiffnessMatrix() ;
}

ViscoelasticityAndFracture::ViscoelasticityAndFracture(ViscoelasticModel m, const Matrix & c_kv, const Matrix & e_kv, const Matrix & c_mx, const Matrix & e_mx,FractureCriterion * c, DamageModel * d,  int b, int n, double r) : Viscoelasticity(m, c_kv, e_kv, c_mx, e_mx, b, n, r), dfunc(d), criterion(c)
{
    setElasticAndViscousStiffnessMatrix() ;
}

ViscoelasticityAndFracture::ViscoelasticityAndFracture(ViscoelasticModel m, const Matrix & c0, std::vector<std::pair<Matrix, Matrix> > & branches, FractureCriterion * c, DamageModel * d, int b, int n, double r) : Viscoelasticity(m, c0, branches, b, n, r), dfunc(d), criterion(c)
{
    setElasticAndViscousStiffnessMatrix() ;

}

ViscoelasticityAndFracture::ViscoelasticityAndFracture(ViscoelasticModel m, const Matrix & c0, const Matrix & c1, const Matrix & e1, FractureCriterion * c, DamageModel * d, int b, int n, double r) : Viscoelasticity(m, c0, c1, e1, b, n, r), dfunc(d), criterion(c)
{
    setElasticAndViscousStiffnessMatrix() ;
}

ViscoelasticityAndFracture::~ViscoelasticityAndFracture()
{
    delete dfunc ;
    delete criterion ;
}

ElementState * ViscoelasticityAndFracture::createElementState( IntegrableEntity * e)
{
    return new GeneralizedSpaceTimeViscoElasticElementState(e) ;
}

FractureCriterion * ViscoelasticityAndFracture::getFractureCriterion() const
{
    return criterion ;
}

DamageModel * ViscoelasticityAndFracture::getDamageModel() const
{
    return dfunc ;
}

void ViscoelasticityAndFracture::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
    Matrix a(ret.numRows()/blocks, ret.numCols()/blocks) ;
    Matrix b(ret.numRows()/blocks, ret.numCols()/blocks) ;
    if(!time_d)
    {

        Matrix buffer(param.numRows()/blocks, param.numCols()/blocks) ;

        switch(model)
        {
        case GENERALIZED_KELVIN_VOIGT:
        {
            std::vector<Matrix>  dtensors ;

            for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
            {
                dtensors.push_back(dfunc->apply(tensors[0], gp.gaussPoints[g].first));
            }

            // stiffness (0,0)
//             buffer = dfunc->apply(tensors[0], gp.gaussPoints[0].first) ;
            vm->ieval(GradientDot(p_i) * dtensors * Gradient(p_j, true),    gp, Jinv,v, a) ;
            vm->ieval(Gradient(p_i)    * dtensors * GradientDot(p_j, true), gp, Jinv,v, b) ;
            a += b ;
            placeMatrixInBlock( a  , 0,0, ret ) ;
            for(int i = 1 ; i < blocks ; i++)
            {
                // first line
                placeMatrixInBlock( -a , i,0, ret ) ;
                // first column
                placeMatrixInBlock( -a , 0,i, ret ) ;
                // diagonal
                placeMatrixInBlock( a , i,i, ret ) ;
                for(int j = i+1 ; j < blocks ; j++)
                {
                    // upper triangle
                    placeMatrixInBlock( a , i,j, ret ) ;
                    // lower triangle
                    placeMatrixInBlock( a , j,i, ret ) ;
                }
            }
            for(int i = 1 ; i < blocks ; i++)
            {

                dtensors.clear() ;
                for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
                {
                    dtensors.push_back(dfunc->applyViscous(tensors[i*2-1], gp.gaussPoints[g].first)) ;
                }
                //stiffness (diagonal)
//                 buffer = dfunc->applyViscous(tensors[i*2-1], gp.gaussPoints[0].first) ;
                vm->ieval(GradientDot(p_i) * dtensors * Gradient(p_j, true),    gp, Jinv,v, a) ;
                vm->ieval(Gradient(p_i)    * dtensors * GradientDot(p_j, true), gp, Jinv,v, b) ;
                a += b ;
                addMatrixInBlock( a , i,i, ret ) ;
            }

            return ;
        }

        case GENERALIZED_MAXWELL:
        {
            // stiffness (0,0)


            std::vector<Matrix>  dtensors ;

            for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
            {
                dtensors.push_back(dfunc->apply(tensors[0], gp.gaussPoints[g].first));
            }

            vm->ieval(GradientDot(p_i) * dtensors * Gradient(p_j, true),    gp, Jinv,v, a) ;
            vm->ieval(Gradient(p_i)    * dtensors * GradientDot(p_j, true), gp, Jinv,v, b) ;
            a += b ;
            placeMatrixInBlock( a , 0,0, ret ) ;
            for(int i = 1 ; i < blocks ; i++)
            {
                dtensors.clear();
                for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
                {
                    dtensors.push_back(dfunc->applyViscous(tensors[2*i-1], gp.gaussPoints[g].first));
                }
//                 buffer = dfunc->applyViscous(tensors[2*i-1], gp.gaussPoints[0].first) ;
                vm->ieval(GradientDot(p_i) * dtensors * Gradient(p_j, true),    gp, Jinv,v, a) ;
                vm->ieval(Gradient(p_i)    * dtensors * GradientDot(p_j, true), gp, Jinv,v, b) ;
                a += b ;
                addMatrixInBlock( a , 0,0, ret ) ;
                placeMatrixInBlock( a , i,i, ret ) ;
                // first line
                placeMatrixInBlock( -a , i,0, ret ) ;
                // first column
                placeMatrixInBlock( -a , 0,i, ret ) ;

            }
            return ;
        }

        case BURGER:
        {
            // stiffness Maxwell
            std::vector<Matrix>  dtensors ;

            for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
            {
                dtensors.push_back(dfunc->apply(tensors[0], gp.gaussPoints[g].first));
            }

//             buffer = dfunc->apply(tensors[0], gp.gaussPoints[0].first) ;
            vm->ieval(GradientDot(p_i) * dtensors * Gradient(p_j, true),    gp, Jinv,v, a) ;
            vm->ieval(Gradient(p_i)    * dtensors * GradientDot(p_j, true), gp, Jinv,v, b) ;
            a += b ;
            placeMatrixInBlock( a , 0,0, ret ) ;
            placeMatrixInBlock( a , 1,1, ret ) ;
            placeMatrixInBlock( a , 2,2, ret ) ;
            placeMatrixInBlock( a , 1,2, ret ) ;
            placeMatrixInBlock( a , 2,1, ret ) ;
            placeMatrixInBlock( -a , 0,1, ret ) ;
            placeMatrixInBlock( -a , 1,0, ret ) ;
            placeMatrixInBlock( -a , 0,2, ret ) ;
            placeMatrixInBlock( -a , 2,0, ret ) ;
            // stiffness KV

            dtensors.clear() ;

            for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
            {
                dtensors.push_back(dfunc->applyViscous(tensors[2], gp.gaussPoints[g].first));
            }
//             buffer = dfunc->applyViscous(tensors[2], gp.gaussPoints[0].first) ;
            vm->ieval(GradientDot(p_i) * dtensors * Gradient(p_j, true),    gp, Jinv,v, a) ;
            vm->ieval(Gradient(p_i)    * dtensors * GradientDot(p_j, true), gp, Jinv,v, b) ;
            a += b ;
            addMatrixInBlock( a , 2,2, ret ) ;

            return ;
        }

        case MAXWELL:
        {
            // stiffness

            std::vector<Matrix>  dtensors ;

            for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
            {
                dtensors.push_back(dfunc->apply(tensors[0], gp.gaussPoints[g].first));
            }
//             buffer = dfunc->apply(tensors[0], gp.gaussPoints[0].first) ;
            vm->ieval(GradientDot(p_i) * dtensors * Gradient(p_j, true),    gp, Jinv,v, a) ;
            vm->ieval(Gradient(p_i)    * dtensors * GradientDot(p_j, true), gp, Jinv,v, b) ;
            a += b ;
            placeMatrixInBlock( a , 0,0, ret ) ;
            placeMatrixInBlock( a , 1,1, ret ) ;
            placeMatrixInBlock( -a , 0,1, ret ) ;
            placeMatrixInBlock( -a , 1,0, ret ) ;
            return ;
        }

        case KELVIN_VOIGT:
        {
            // stiffness
            std::vector<Matrix>  dtensors ;

            for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
            {
                dtensors.push_back(dfunc->apply(tensors[0], gp.gaussPoints[g].first));
            }

//             buffer = dfunc->apply(tensors[0], gp.gaussPoints[0].first) ;
            vm->ieval(GradientDot(p_i) * dtensors * Gradient(p_j, true),    gp, Jinv,v, a) ;
            vm->ieval(Gradient(p_i)    * dtensors * GradientDot(p_j, true), gp, Jinv,v, b) ;
            a += b ;
            placeMatrixInBlock( a , 0,0, ret ) ;
            return ;
        }

        case PURE_ELASTICITY:
        {
            std::vector<Matrix>  dtensors ;

            for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
            {
                dtensors.push_back(dfunc->apply(tensors[0], gp.gaussPoints[g].first));
            }
//             buffer = dfunc->apply(tensors[0], gp.gaussPoints[0].first) ;
            vm->ieval(GradientDot(p_i) * dtensors * Gradient(p_j, true),    gp, Jinv,v, a) ;
            vm->ieval(Gradient(p_i)    * dtensors * GradientDot(p_j, true), gp, Jinv,v, b) ;
            a += b ;
            placeMatrixInBlock( a , 0,0, ret ) ;

            return ;
        }

        case PURE_VISCOSITY:
        {
            return ;
        }

        default:
        {
            for(int i = 0 ; i < blocks ; i++)
            {

                std::vector<Matrix>  dtensors ;


                // elasticParam matrix (diagonal)
                getBlockInMatrix(param, i,i, buffer) ;

                for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
                {
                    dtensors.push_back(dfunc->apply(buffer, gp.gaussPoints[g].first));
                }
//                 buffer = dfunc->apply(buffer, gp.gaussPoints[0].first) ;
                vm->ieval(GradientDot(p_i) * dtensors * Gradient(p_j, true),    gp, Jinv,v, a) ;
                vm->ieval(Gradient(p_i)    * dtensors * GradientDot(p_j, true), gp, Jinv,v, b) ;
                a += b ;
                placeMatrixInBlock( a , i,i, ret ) ;
                // elasticParam matrix (upper-triangle)
                for(int j = i+1 ; j < blocks ; j++)
                {
                    getBlockInMatrix(param, i,j, buffer) ;
                    dtensors.clear();
                    for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
                    {
                        dtensors.push_back(dfunc->apply(buffer, gp.gaussPoints[g].first));
                    }
//                     buffer = dfunc->apply(buffer, gp.gaussPoints[0].first) ;
                    vm->ieval(GradientDot(p_i) * dtensors * Gradient(p_j, true),    gp, Jinv,v, a) ;
                    vm->ieval(Gradient(p_i)    * dtensors * GradientDot(p_j, true), gp, Jinv,v, b) ;
                    a += b ;
                    placeMatrixInBlock( a , i,j, ret ) ;
                    // symmetry
                    placeMatrixInBlock( a , j,i, ret ) ;
                }

            }
        }
        }

        return ;
    }
    else
    {
        std::vector<Matrix>  dtensors ;
        for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
            dtensors.push_back(Matrix(param.numRows()/blocks, param.numCols()/blocks)) ;

        for(int i = 0 ; i < blocks ; i++)
        {
            for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
            {
                getBlockInMatrix( param, i,i,  dtensors[g] ) ;
                dtensors[g] = dfunc->apply( dtensors[g], gp.gaussPoints[g].first ) ;
            }

            vm->ieval(GradientDot(p_i) *  dtensors * Gradient(p_j, true),    gp, Jinv,v, a) ;
            vm->ieval(Gradient(p_i)    *  dtensors * GradientDot(p_j, true), gp, Jinv,v, b) ;
            a += b ;
            placeMatrixInBlock( a  , i,i, ret ) ;

            for(int j = i+1 ; j < blocks ; j++)
            {
                for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
                {
                    getBlockInMatrix( param, i,j,  dtensors[g] ) ;
                    dtensors[g] = dfunc->applyViscous(dtensors[g], gp.gaussPoints[g].first ) ;
                }

                vm->ieval(GradientDot(p_i) *  dtensors * Gradient(p_j, true),    gp, Jinv,v, a) ;
                vm->ieval(Gradient(p_i)    *  dtensors * GradientDot(p_j, true), gp, Jinv,v, b) ;
                a += b ;
                placeMatrixInBlock( a  , i,j, ret ) ;
                placeMatrixInBlock( a  , j,i, ret ) ;
            }
        }

    }


}

void ViscoelasticityAndFracture::applyViscous(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
    Matrix a(ret.numRows()/blocks, ret.numCols()/blocks) ;
    Matrix b(ret.numRows()/blocks, ret.numCols()/blocks) ;
    if(!time_d)
    {

        Matrix buffer(param.numRows()/blocks, param.numCols()/blocks) ;

        switch(model)
        {
        case GENERALIZED_KELVIN_VOIGT:
        {
            for(int i = 1 ; i < blocks ; i++)
            {
                // viscosity (diagonal)
                std::vector<Matrix>  dtensors ;

                for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
                {
                    dtensors.push_back(dfunc->applyViscous(tensors[2*i], gp.gaussPoints[g].first));
                }
//                 buffer = dfunc->applyViscous(tensors[2*i]) ;
//                a = 0 ;
                vm->ieval(GradientDot(p_i) * dtensors * GradientDot(p_j, true),    gp, Jinv,v, a) ;
                vm->ieval(GradientDotDot(p_i) * dtensors * Gradient(p_j, true), gp, Jinv,v, b) ;
                a += b ;
                placeMatrixInBlock( a  , i,i, ret ) ;
            }

            return ;
        }

        case GENERALIZED_MAXWELL:
        {
            for(int i = 1 ; i < blocks ; i++)
            {
                std::vector<Matrix>  dtensors ;

                for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
                {
                    dtensors.push_back(dfunc->applyViscous(tensors[2*i], gp.gaussPoints[g].first));
                }

                // viscosity (diagonal)
//                 buffer = dfunc->applyViscous(tensors[2*i]) ;
                vm->ieval(GradientDot(p_i) * dtensors * GradientDot(p_j, true),    gp, Jinv,v, a) ;
                vm->ieval(GradientDotDot(p_i)    * dtensors * Gradient(p_j, true), gp, Jinv,v, b) ;
                a += b ;
                placeMatrixInBlock( a  , i,i, ret ) ;
            }
            return ;
        }

        case BURGER:
        {
            // viscosity Maxwell
            std::vector<Matrix>  dtensors ;

            for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
            {
                dtensors.push_back(dfunc->applyViscous(tensors[1], gp.gaussPoints[g].first));
            }

//             buffer = dfunc->applyViscous(tensors[1]) ;
            vm->ieval(GradientDot(p_i)    * dtensors   * GradientDot(p_j, true), gp, Jinv,v,a);
            vm->ieval(GradientDotDot(p_i) * dtensors   * Gradient(p_j, true),    gp, Jinv,v,b);
            a += b ;
            placeMatrixInBlock( a  , 1,1, ret ) ;
            // viscosity KV
//             buffer = dfunc->applyViscous(tensors[3]) ;

            dtensors.clear() ;

            for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
            {
                dtensors.push_back(dfunc->applyViscous(tensors[3], gp.gaussPoints[g].first));
            }
            vm->ieval(GradientDot(p_i)    * dtensors   * GradientDot(p_j, true), gp, Jinv,v,a);
            vm->ieval(GradientDotDot(p_i) * dtensors   * Gradient(p_j, true),    gp, Jinv,v,b);
            a += b ;
            placeMatrixInBlock( a  , 2,2, ret ) ;
            return ;
        }

        case MAXWELL:
        {
            // viscosity
            std::vector<Matrix>  dtensors ;

            for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
            {
                dtensors.push_back(dfunc->applyViscous(tensors[1], gp.gaussPoints[g].first));
            }
//             buffer = dfunc->applyViscous(tensors[1]) ;
            vm->ieval(GradientDot(p_i)    * dtensors   * GradientDot(p_j, true), gp, Jinv,v,a);
            vm->ieval(GradientDotDot(p_i) * dtensors   * Gradient(p_j, true),    gp, Jinv,v,b);
            a += b ;
            placeMatrixInBlock( a  , 1,1, ret ) ;
            return ;
        }

        case KELVIN_VOIGT:
        {
            // viscosity
            std::vector<Matrix>  dtensors ;

            for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
            {
                dtensors.push_back(dfunc->applyViscous(tensors[1], gp.gaussPoints[g].first));
            }
//             buffer = dfunc->applyViscous(tensors[1]) ;
            vm->ieval(GradientDot(p_i)    * dtensors   * GradientDot(p_j, true), gp, Jinv,v,a);
            vm->ieval(GradientDotDot(p_i) * dtensors   * Gradient(p_j, true),    gp, Jinv,v,b);
            addMatrixInBlock( a  , 0,0, ret ) ;
            return ;
        }

        case PURE_ELASTICITY:
        {
            return ;
        }

        case PURE_VISCOSITY:
        {
            std::vector<Matrix>  dtensors ;

            for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
            {
                dtensors.push_back(dfunc->applyViscous(tensors[1], gp.gaussPoints[g].first));
            }

//             buffer = dfunc->applyViscous(tensors[1]) ;
            vm->ieval(GradientDot(p_i)    * dtensors   * GradientDot(p_j, true), gp, Jinv,v,a);
            vm->ieval(GradientDotDot(p_i) * dtensors   * Gradient(p_j, true),    gp, Jinv,v,b);
            a += b ;
            placeMatrixInBlock( a  , 0,0, ret ) ;

            return ;
        }

        default:
        {
            for(int i = 0 ; i < blocks ; i++)
            {
                // viscousParam matrix (diagonal)
                getBlockInMatrix(eta, i,i, buffer) ;
                std::vector<Matrix>  dtensors ;

                for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
                {
                    dtensors.push_back(dfunc->applyViscous(buffer, gp.gaussPoints[g].first));
                }

//                 buffer = dfunc->applyViscous(buffer) ;
                vm->ieval(GradientDot(p_i)    * dtensors   * GradientDot(p_j, true), gp, Jinv,v,a);
                vm->ieval(GradientDotDot(p_i) * dtensors   * Gradient(p_j, true),    gp, Jinv,v,b);
                a += b ;
                placeMatrixInBlock( a  , i,i, ret ) ;
                // viscousParam matrix (upper-triangle)
                for(int j = i+1 ; j < blocks ; j++)
                {
                    dtensors.clear() ;

                    for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
                    {
                        dtensors.push_back(dfunc->applyViscous(buffer, gp.gaussPoints[g].first));
                    }
                    getBlockInMatrix(eta, i,j, buffer) ;
                    buffer = dfunc->applyViscous(buffer) ;
                    vm->ieval(GradientDot(p_i)    * dtensors   * GradientDot(p_j, true), gp, Jinv,v,a);
                    vm->ieval(GradientDotDot(p_i) * dtensors   * Gradient(p_j, true),    gp, Jinv,v,b);
                    a += b ;
                    addMatrixInBlock( a  , i,j, ret ) ;
                    // symmetry
                    addMatrixInBlock( a  , j,i, ret ) ;
                }
            }
        }
        }

        return ;
    }
    else
    {
        std::vector<Matrix>  dtensors ;
        for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
            dtensors.push_back(Matrix(param.numRows()/blocks, param.numCols()/blocks)) ;

        for(int i = 0 ; i < blocks ; i++)
        {
            for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
            {
                getBlockInMatrix( eta, i,i,  dtensors[g] ) ;
                dtensors[g] = dfunc->applyViscous( dtensors[g], gp.gaussPoints[g].first ) ;
            }

            vm->ieval(GradientDotDot(p_i) *  dtensors * Gradient(p_j, true),    gp, Jinv,v, a) ;
            vm->ieval(GradientDot(p_i)    *  dtensors * GradientDot(p_j, true), gp, Jinv,v, b) ;
            a += b ;
            placeMatrixInBlock( a  , i,i, ret ) ;

            for(int j = i+1 ; j < blocks ; j++)
            {
                for(size_t g = 0 ; g < gp.gaussPoints.size() ; g++)
                {
                    getBlockInMatrix( eta, i,j,  dtensors[g] ) ;
                    dtensors[g] = dfunc->applyViscous( dtensors[g], gp.gaussPoints[g].first ) ;
                }

                vm->ieval(GradientDotDot(p_i) *  dtensors * Gradient(p_j, true),    gp, Jinv,v, a) ;
                vm->ieval(GradientDot(p_i)    *  dtensors * GradientDot(p_j, true), gp, Jinv,v, b) ;
                a += b ;
                placeMatrixInBlock( a  , i,j, ret ) ;
                placeMatrixInBlock( a  , j,i, ret ) ;
            }
        }
    }



}

void ViscoelasticityAndFracture::step(double timestep, ElementState & currentState, double maxscore)
{
    dfunc->step(currentState, maxscore) ;
    currentState.getParent()->behaviourUpdated = dfunc->changed() ;
}


bool ViscoelasticityAndFracture::fractured() const
{
    return dfunc->fractured() ;
}

bool ViscoelasticityAndFracture::changed() const
{
    return dfunc->changed() ;
}

Form * ViscoelasticityAndFracture::getCopy() const
{

    ViscoelasticityAndFracture * copy ;
    switch(model)
    {
    case PURE_ELASTICITY:
        copy = new ViscoelasticityAndFracture( PURE_ELASTICITY, tensors[0], criterion->getCopy(), dfunc->getCopy(), blocks-1) ;
        break ;
    case PURE_VISCOSITY:
        copy = new ViscoelasticityAndFracture( PURE_VISCOSITY, tensors[1], criterion->getCopy(), dfunc->getCopy(), blocks-1) ;
        break ;
    case MAXWELL:
        copy = new ViscoelasticityAndFracture( MAXWELL, tensors[0], tensors[1], criterion->getCopy(), dfunc->getCopy(), blocks-2) ;
        break ;
    case KELVIN_VOIGT:
        copy = new ViscoelasticityAndFracture( KELVIN_VOIGT, tensors[0], tensors[1], criterion->getCopy(), dfunc->getCopy(), blocks-1) ;
        break ;
    case GENERALIZED_MAXWELL:
    {
        std::vector<std::pair<Matrix, Matrix> > branches ;
        for(size_t i = 1 ; i < tensors.size() ; i += 2)
            branches.push_back(std::make_pair(tensors[i], tensors[i+1])) ;
        copy = new ViscoelasticityAndFracture( GENERALIZED_MAXWELL, tensors[0], branches, criterion->getCopy(), dfunc->getCopy()) ;
        break ;
    }
    case GENERALIZED_KELVIN_VOIGT:
    {
        std::vector<std::pair<Matrix, Matrix> > branches ;
        for(size_t i = 1 ; i < tensors.size() ; i += 2)
            branches.push_back(std::make_pair(tensors[i], tensors[i+1])) ;
        copy = new ViscoelasticityAndFracture( GENERALIZED_KELVIN_VOIGT, tensors[0], branches, criterion->getCopy(), dfunc->getCopy()) ;
        break ;
    }
    case BURGER:
        copy = new ViscoelasticityAndFracture( BURGER, tensors[0], tensors[1], tensors[2], tensors[3], criterion->getCopy(), dfunc->getCopy(), blocks-3) ;
        break ;
    default:
        copy = new ViscoelasticityAndFracture(  param, eta, blocks, criterion->getCopy(), dfunc->getCopy())  ;

    }
    copy->model = model ;
    copy->dfunc->getState(true).resize(dfunc->getState().size());
    copy->dfunc->getState(true) = dfunc->getState() ;
    copy->criterion->setMaterialCharacteristicRadius(criterion->getMaterialCharacteristicRadius()) ;
    copy->dfunc->setDamageDensityTolerance(dfunc->getDamageDensityTolerance());
    copy->dfunc->setThresholdDamageDensity(dfunc->getThresholdDamageDensity());

    return copy ;
}

Vector ViscoelasticityAndFracture::getForcesFromAppliedStress( const Vector & data, Function & shape, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<Variable> & v)
{
    return VirtualMachine().ieval(GradientDot( shape ) * ( data ), gp, Jinv, v) ;
}

Vector ViscoelasticityAndFracture::getForcesFromAppliedStress( const Function & data, size_t index, size_t externaldofs,  Function & shape, IntegrableEntity * e,const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<Variable> & v)
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

    f += vm.ieval( Gradient( shape ) * g, e, v) ;

    return f ;
}

void ViscoelasticityAndFracture::setElasticAndViscousStiffnessMatrix()
{
    elasticParam.resize( param.numRows(), param.numCols()) ;
    viscousParam.resize( param.numRows(), param.numCols()) ;
    switch(model)
    {
    case PURE_ELASTICITY:
        placeMatrixInBlock( tensors[0], 0,0, elasticParam) ;
        break ;
    case PURE_VISCOSITY:
        break ;
    case KELVIN_VOIGT:
        placeMatrixInBlock( tensors[0], 0,0, elasticParam) ;
        break ;
    case MAXWELL:
        placeMatrixInBlock( tensors[0], 0,0, elasticParam) ;
        placeMatrixInBlock( tensors[0], 1,0, elasticParam) ;
        placeMatrixInBlock( -tensors[0], 0,1, elasticParam) ;
        placeMatrixInBlock( -tensors[0], 1,1, elasticParam) ;
        break ;
    case BURGER:
        placeMatrixInBlock( tensors[0], 0,0, elasticParam) ;
        placeMatrixInBlock( -tensors[0], 1,0, elasticParam) ;
        placeMatrixInBlock( -tensors[0], 2,0, elasticParam) ;
        placeMatrixInBlock( -tensors[0], 0,1, elasticParam) ;
        placeMatrixInBlock( -tensors[0], 0,2, elasticParam) ;
        placeMatrixInBlock( tensors[0], 1,1, elasticParam) ;
        placeMatrixInBlock( tensors[0], 1,2, elasticParam) ;
        placeMatrixInBlock( tensors[0], 2,1, elasticParam) ;
        placeMatrixInBlock( tensors[0], 2,2, elasticParam) ;
        addMatrixInBlock( tensors[2], 2,2, viscousParam) ;
        break ;
    case GENERALIZED_MAXWELL:
        placeMatrixInBlock( tensors[0], 0,0, elasticParam) ;
        for(size_t i = 1 ; i < tensors.size() ; i+=2)
        {
            Matrix ri = tensors[i] * (-1) ;
            addMatrixInBlock( tensors[i], 0,0, elasticParam) ;
            placeMatrixInBlock( ri, 0,i/2+1, elasticParam) ;
            placeMatrixInBlock( ri, i/2+1,0, elasticParam) ;
            placeMatrixInBlock( tensors[i], i/2+1,i/2+1, viscousParam) ;
        }
        break ;
    case GENERALIZED_KELVIN_VOIGT:
        placeMatrixInBlock( tensors[0], 0,0, elasticParam) ;
        for(size_t i = 1 ; i < tensors.size() ; i+=2)
        {
            placeMatrixInBlock( tensors[0], i/2+1,i/2+1, elasticParam) ;
            placeMatrixInBlock( -tensors[0], i/2+1,0, elasticParam) ;
            placeMatrixInBlock( -tensors[0], 0,i/2+1, elasticParam) ;
            for(size_t j = 1 ; j < i/2+1 ; j++)
            {
                placeMatrixInBlock( tensors[0], i/2+1,j, elasticParam) ;
                placeMatrixInBlock( tensors[0], j,i/2+1, elasticParam) ;
            }
            placeMatrixInBlock( tensors[i], i/2+1,i/2+1, viscousParam) ;
        }
        break ;
    default:
        elasticParam=param ;
        break ;
    }

}

Matrix ViscoelasticityAndFracture::getTensor(const Point & p, IntegrableEntity * e, int g) const
{
    return  dfunc->apply(elasticParam, p) + dfunc->applyViscous(viscousParam, p) ;
}

Matrix ViscoelasticityAndFracture::getViscousTensor(const Point & p, IntegrableEntity * e, int g) const
{
    return dfunc->apply( eta, p ) ;
}

void ViscoelasticityAndFracture::setFractureCriterion(FractureCriterion * frac)
{
    if(frac)
    {
        criterion = frac ;
    }

}

