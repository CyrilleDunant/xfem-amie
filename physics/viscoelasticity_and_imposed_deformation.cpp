#include "viscoelasticity_and_imposed_deformation.h"
#include "../features/boundarycondition.h"
#include "../elements/generalized_spacetime_viscoelastic_element_state.h"

using namespace Mu ;

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

ViscoelasticityAndImposedDeformation::~ViscoelasticityAndImposedDeformation() {} ;

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
	if(model == PURE_ELASTICITY)
	{
		Matrix rig( param.numCols()/blocks, param.numRows()/blocks) ;
		for(size_t i = 0 ; i < rig.numCols() ; i++)
		{
			for(size_t j = 0 ; j < rig.numRows() ; j++)
				rig[i][j] = param[i][j] ;
		}
		Vector imp = imposedStrain ;
		return new ViscoelasticityAndImposedDeformation( model, rig, imp, blocks-1, rho) ;
	}
	return new ViscoelasticityAndImposedDeformation(*this) ;
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
// 		Vector tmp = imposedStress ;
// 		tmp[17] = imposedStrain[24] ;

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
