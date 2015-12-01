
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//

#include "../physics/stiffness.h"
#include "../physics/stiffness_with_imposed_deformation.h"
#include "../physics/finite_difference_viscoelasticity.h"
#include "../physics/viscoelasticity.h"
#include "../physics/viscoelasticity_and_imposed_deformation.h"
#include "../physics/materials/aggregate_behaviour.h"
#include "../physics/materials/paste_behaviour.h"
#include "../physics/materials/gel_behaviour.h"
#include "../physics/materials/rebar_behaviour.h"
#include "../physics/materials/steel_behaviour.h"
#include "../physics/materials/concrete_behaviour.h"
#include "behaviour_converter.h"


using namespace Amie ;

bool BehaviourConverter::toSpaceTimeBehaviour( Feature * f, int maxBlocks ) 
{
     Point p ;

     Form * behaviour = f->getBehaviour() ;
     if( !behaviour )
     {
         std::cout << "warning: null-behaviour" << std::endl ;
         return true ;
     }

     // VOID behaviour, nothing to do
     if(behaviour->type == VOID_BEHAVIOUR)
     {
         return true ;
     }

     // visco-elastic behaviour
     // adjust size of blocks in elementary matrices
     else if(dynamic_cast<Viscoelasticity *>(behaviour))
     {
         dynamic_cast<Viscoelasticity *>(behaviour)->setBlocks( maxBlocks ) ;
         return true ;
     }

     // finite-difference visco-elastic behaviour
     // return equivalent space-time behaviour
     else if(dynamic_cast<FiniteDifferenceViscoelasticity *>(behaviour))
     {
         ViscoelasticModel model = dynamic_cast<FiniteDifferenceViscoelasticity *>(behaviour)->model ;
         Viscoelasticity * equivalent = nullptr ;
         switch(model)
         {
         case PURE_ELASTICITY:
             equivalent = new Viscoelasticity( PURE_ELASTICITY, dynamic_cast<FiniteDifferenceViscoelasticity *>(behaviour)->tensors[0], maxBlocks-1 ) ;
             break ;
         case PURE_VISCOSITY:
             equivalent = new Viscoelasticity( PURE_VISCOSITY, dynamic_cast<FiniteDifferenceViscoelasticity *>(behaviour)->tensors[1], maxBlocks-1 ) ;
             break ;
         case KELVIN_VOIGT:
             equivalent = new Viscoelasticity( KELVIN_VOIGT, dynamic_cast<FiniteDifferenceViscoelasticity *>(behaviour)->tensors[0], dynamic_cast<FiniteDifferenceViscoelasticity *>(behaviour)->tensors[1], maxBlocks-1 ) ;
         case MAXWELL:
             equivalent = new Viscoelasticity( MAXWELL, dynamic_cast<FiniteDifferenceViscoelasticity *>(behaviour)->tensors[0], dynamic_cast<FiniteDifferenceViscoelasticity *>(behaviour)->tensors[1], maxBlocks-2 ) ;
             break ;
         case BURGER:
             equivalent = new Viscoelasticity( BURGER, dynamic_cast<FiniteDifferenceViscoelasticity *>(behaviour)->tensors[0], dynamic_cast<FiniteDifferenceViscoelasticity *>(behaviour)->tensors[1], dynamic_cast<FiniteDifferenceViscoelasticity *>(behaviour)->tensors[2], dynamic_cast<FiniteDifferenceViscoelasticity *>(behaviour)->tensors[3], maxBlocks-3 ) ;
             break ;
         case GENERALIZED_KELVIN_VOIGT:
         {
             std::vector<std::pair<Matrix, Matrix > > branches ;
             int blocks = (dynamic_cast<FiniteDifferenceViscoelasticity *>(behaviour)->tensors.size()-1)/2 ;
             for(int i = 0 ; i < blocks ; i++)
                 branches.push_back(std::make_pair( dynamic_cast<FiniteDifferenceViscoelasticity *>(behaviour)->tensors[i*2+1], dynamic_cast<FiniteDifferenceViscoelasticity *>(behaviour)->tensors[i*2+2]) ) ;
             equivalent = new Viscoelasticity( GENERALIZED_KELVIN_VOIGT, dynamic_cast<FiniteDifferenceViscoelasticity *>(behaviour)->tensors[0], branches, maxBlocks-blocks-1 ) ;
             break ;
         }
         case GENERALIZED_MAXWELL:
         {
             std::vector<std::pair<Matrix, Matrix > > branches ;
             int blocks = (dynamic_cast<FiniteDifferenceViscoelasticity *>(behaviour)->tensors.size()-1)/2 ;
             for(int i = 0 ; i < blocks ; i++)
                 branches.push_back(std::make_pair( dynamic_cast<FiniteDifferenceViscoelasticity *>(behaviour)->tensors[i*2+1], dynamic_cast<FiniteDifferenceViscoelasticity *>(behaviour)->tensors[i*2+2]) ) ;
             equivalent = new Viscoelasticity( GENERALIZED_MAXWELL, dynamic_cast<FiniteDifferenceViscoelasticity *>(behaviour)->tensors[0], branches, maxBlocks-blocks-1 ) ;
             break ;
         }
         case GENERAL_VISCOELASTICITY:
             return false ;
         }
         if(equivalent)
         {
             f->setBehaviour(equivalent) ;
             return true ;
         }
         return false ;
     }


     // predefined aggregate behaviour
     // use predefined space-time behaviour
     else if(dynamic_cast<AggregateBehaviour *>(behaviour))
     {
         dynamic_cast<AggregateBehaviour *>(behaviour)->spaceTime = true ;
         dynamic_cast<AggregateBehaviour *>(behaviour)->freeblocks = maxBlocks-3 ;
         return true ;
     }

     // predefined cement paste behaviour
     // use predefined space-time behaviour
     else if(dynamic_cast<PasteBehaviour *>(behaviour))
     {
         dynamic_cast<PasteBehaviour *>(behaviour)->spaceTime = true ;
         dynamic_cast<PasteBehaviour *>(behaviour)->freeblocks = maxBlocks-3 ;
         return true ;
     }
     // predefined ASR gel behaviour
     // use predefined space-time behaviour
     else if(dynamic_cast<GelBehaviour *>(behaviour))
     {
         dynamic_cast<GelBehaviour *>(behaviour)->spaceTime = true ;
         dynamic_cast<GelBehaviour *>(behaviour)->freeblocks = maxBlocks-3 ;
         return true ;
     }

     // generic elastic behaviour
     // attempt space-time conversion
     else if(dynamic_cast<Stiffness *>(behaviour))
     {
         f->setBehaviour( new Viscoelasticity( PURE_ELASTICITY, behaviour->getTensor( p ), maxBlocks-1) ) ;
         return true ;
     }

     // generic elastic with imposed deformation behaviour
     // attempt space-time conversion
     else if(dynamic_cast<StiffnessWithImposedDeformation *>(behaviour))
     {
         Vector v = behaviour->getImposedStrain( p ) ;
         f->setBehaviour( new ViscoelasticityAndImposedDeformation( PURE_ELASTICITY, behaviour->getTensor( p ), v, maxBlocks-1 ) ) ;
         return true ;
     }

     return false ; 
}
