#include "viscoelasticity_and_imposed_deformation.h"
#include "../features/boundarycondition.h"
#include "../elements/generalized_spacetime_viscoelastic_element_state.h"

using namespace Amie ;

ViscoelasticityAndImposedDeformation::ViscoelasticityAndImposedDeformation( ViscoelasticModel model, const Matrix & rig, int additionnalBlocksAfter, double r ) : Viscoelasticity(model, rig, additionnalBlocksAfter, r), imposedStrain(rig.numRows())
{
    imposedStrain = 0. ;
    makeImposedStress() ;
}

ViscoelasticityAndImposedDeformation::ViscoelasticityAndImposedDeformation( ViscoelasticModel model, const Matrix & rig, Vector & imp, int additionnalBlocksAfter, double r ) : Viscoelasticity(model, rig, additionnalBlocksAfter, r), imposedStrain(imp)
{
    makeImposedStress() ;
}

ViscoelasticityAndImposedDeformation::ViscoelasticityAndImposedDeformation( ViscoelasticModel model, const Matrix & rig, const Matrix & eta,  Vector & imp, int additionnalBlocksBefore, int additionnalBlocksAfter, double r) : Viscoelasticity( model, rig, eta, additionnalBlocksBefore, additionnalBlocksAfter, r), imposedStrain(imp)
{
    makeImposedStress() ;
}

ViscoelasticityAndImposedDeformation::ViscoelasticityAndImposedDeformation( ViscoelasticModel model, const Matrix & c_kv, const Matrix & e_kv, const Matrix & c_mx, const Matrix & e_mx,  Vector & imp, int additionnalBlocksBefore , int additionnalBlocksAfter , double r) : Viscoelasticity(model, c_kv, e_kv, c_mx, e_mx, additionnalBlocksBefore, additionnalBlocksAfter, r), imposedStrain(imp)
{
    makeImposedStress() ;
}

ViscoelasticityAndImposedDeformation::ViscoelasticityAndImposedDeformation( ViscoelasticModel model, const Matrix & c_0, std::vector<std::pair<Matrix, Matrix> > & branches,  Vector & imp, int additionnalBlocksBefore, int additionnalBlocksAfter, double r)  : Viscoelasticity(model, c_0, branches, additionnalBlocksBefore, additionnalBlocksAfter, r), imposedStrain(imp)
{
    makeImposedStress() ;
}

ViscoelasticityAndImposedDeformation::ViscoelasticityAndImposedDeformation( ViscoelasticModel model, const Matrix & c_0, const Matrix & c_1, const Matrix & e_1, Vector & imp,  int additionnalBlocksBefore, int additionnalBlocksAfter, double r) : Viscoelasticity(model, c_0, c_1, e_1, additionnalBlocksBefore, additionnalBlocksAfter, r), imposedStrain(imp)
{
    makeImposedStress() ;
}

ViscoelasticityAndImposedDeformation::ViscoelasticityAndImposedDeformation( const Matrix & rig, const Matrix & eta, int blocks, Vector & imp, int additionnalBlocksAfter, double r) : Viscoelasticity(rig, eta, blocks, additionnalBlocksAfter, r), imposedStrain(imp)
{
    makeImposedStress() ;
}


void ViscoelasticityAndImposedDeformation::makeImposedStress()
{
    imposedGeneralizedStrain.resize( imposedStrain.size()*blocks ) ;
    imposedGeneralizedStrain = 0 ;
    for(size_t i = 0 ; i < imposedStrain.size() ; i++)
        imposedGeneralizedStrain[i] = imposedStrain[i] ;

    Vector imposedGeneralizedStress = param * imposedGeneralizedStrain ;

    imposedStress.resize(imposedStrain.size()) ;
    for(size_t i = 0 ; i < imposedStress.size() ; i++)
        imposedStress[i] = imposedGeneralizedStress[i] ;

}

ViscoelasticityAndImposedDeformation::~ViscoelasticityAndImposedDeformation() {}

void ViscoelasticityAndImposedDeformation::apply( const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
    Viscoelasticity::apply(p_i, p_j, gp, Jinv, ret, vm) ;
}

void ViscoelasticityAndImposedDeformation::applyViscous( const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm ) const
{
    Viscoelasticity::applyViscous(p_i, p_j, gp, Jinv, ret, vm) ;
}


Form * ViscoelasticityAndImposedDeformation::getCopy() const
{
    Matrix rig = param ;
    Matrix e = eta ;
    Vector a = imposedStrain ;
    switch(model)
    {
    case PURE_ELASTICITY:
        return new ViscoelasticityAndImposedDeformation( PURE_ELASTICITY, tensors[0], a, blocks-effblocks, rho) ;
    case PURE_VISCOSITY:
        return new ViscoelasticityAndImposedDeformation( PURE_VISCOSITY, tensors[1], a, blocks-effblocks, rho) ;
    case MAXWELL:
        return new ViscoelasticityAndImposedDeformation( MAXWELL, tensors[0], tensors[1],a, 0, blocks-effblocks, rho) ;
    case KELVIN_VOIGT:
        return new ViscoelasticityAndImposedDeformation( KELVIN_VOIGT, tensors[0], tensors[1],a, 0, blocks-effblocks, rho) ;
    case BURGER:
        return new ViscoelasticityAndImposedDeformation( BURGER, tensors[2], tensors[3], tensors[0], tensors[1],a, 0, blocks-effblocks, rho) ;
    case GENERALIZED_MAXWELL:
    {
        std::vector<std::pair<Matrix, Matrix> > branches ;
        for(size_t i = 1 ; i < tensors.size() ; i+= 2)
        {
            branches.push_back( std::make_pair(tensors[i], tensors[i+1]) ) ;
        }
        return new ViscoelasticityAndImposedDeformation( GENERALIZED_MAXWELL, tensors[0], branches,a, 0, blocks-effblocks, rho) ;
    }
    case GENERALIZED_KELVIN_VOIGT:
    {
        std::vector<std::pair<Matrix, Matrix> > branches ;
        for(size_t i = 1 ; i < tensors.size() ; i+= 2)
        {
            branches.push_back( std::make_pair(tensors[i], tensors[i+1]) ) ;
        }
        return new ViscoelasticityAndImposedDeformation( GENERALIZED_KELVIN_VOIGT, tensors[0], branches,a, 0, blocks-effblocks, rho) ;
    }
    case GENERAL_VISCOELASTICITY:
        return new ViscoelasticityAndImposedDeformation( rig, e, blocks, a) ;
    }

    ViscoelasticityAndImposedDeformation * copy = new ViscoelasticityAndImposedDeformation(rig, e, blocks, a) ;
    copy->model = model ;

    return copy ;
}

Vector ViscoelasticityAndImposedDeformation::getImposedStrain(const Point & p, IntegrableEntity * e , int g ) const
{
    return imposedStrain ;
}

Vector ViscoelasticityAndImposedDeformation::getImposedStress(const Point & p, IntegrableEntity * e , int g ) const
{
    return imposedStress ;
}

std::vector<BoundaryCondition * > ViscoelasticityAndImposedDeformation::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
    std::vector<BoundaryCondition * > ret ;
    if(v.size() == 3)
    {
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, imposedStress[0]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, imposedStress[1]));
    }
    if(v.size() == 4)
    {

        ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, imposedStress[0]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, imposedStress[1]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, imposedStress[2]));
    }
    return ret ;

}

ViscoelasticityAndVariableImposedDeformation::ViscoelasticityAndVariableImposedDeformation( ViscoelasticModel model, const Matrix & rig, const Function & f, int additionnalBlocksAfter, double r ) : ViscoelasticityAndImposedDeformation(model, rig, additionnalBlocksAfter, r), variableImposed(f)
{
}

Vector ViscoelasticityAndVariableImposedDeformation::getImposedStrain(const Point & p, IntegrableEntity * e , int g ) const
{
    Point q = p ;
    if(e)
    {
        q.getX() = VirtualMachine().eval(e->getXTransform(), p) ;
        q.getY() = VirtualMachine().eval(e->getYTransform(), p) ;
        q.getT() = VirtualMachine().eval(e->getTTransform(), p) ;
    }

    Vector v = imposedStrain ;
    for(size_t i = 0 ; i < 2+(v.size()==6) ; i++)
        v[i] = VirtualMachine().eval( variableImposed , q ) ;
    return v ;
}

Vector ViscoelasticityAndVariableImposedDeformation::getImposedStress(const Point & p, IntegrableEntity * e , int g ) const
{
    Vector strain = getImposedStrain(p,e,g) ;
    Vector genstrain( imposedGeneralizedStrain.size() ) ;
    genstrain = 0 ;
    for(size_t i = 0 ; i < imposedStrain.size() ; i++)
        genstrain[i] = strain[i] ;

    Vector genstress = param * genstrain ;

    Vector stress(strain.size()) ;
    for(size_t i = 0 ; i < imposedStress.size() ; i++)
        stress[i] = genstress[i] ;
    return stress ;
}

std::vector<BoundaryCondition * > ViscoelasticityAndVariableImposedDeformation::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
    Function f = variableImposed*(param[0][0]+param[0][1]+param[0][2]) ;
    std::vector<BoundaryCondition * > ret ;
    if(v.size() == 3)
    {
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, f));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, f));
    }
    return ret ;
}

Form * ViscoelasticityAndVariableImposedDeformation::getCopy() const
{
    Function f(variableImposed) ;
    if(blocks > 1)
    {
        Matrix rig( param.numCols()/blocks, param.numRows()/blocks) ;
        for(size_t i = 0 ; i < rig.numCols() ; i++)
        {
            for(size_t j = 0 ; j < rig.numRows() ; j++)
                rig[i][j] = param[i][j] ;
        }
//			Vector imp = imposedStrain ;
        return new ViscoelasticityAndVariableImposedDeformation( model, rig, f, blocks-1, rho) ;
    }
    Matrix rig = param ;
    return new ViscoelasticityAndVariableImposedDeformation( model, rig, f) ;
}

void ViscoelasticityAndVariableImposedDeformation::apply( const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
    Viscoelasticity::apply(p_i, p_j, gp, Jinv, ret, vm) ;
}

void ViscoelasticityAndVariableImposedDeformation::applyViscous( const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm ) const
{
    Viscoelasticity::applyViscous(p_i, p_j, gp, Jinv, ret, vm) ;
}

