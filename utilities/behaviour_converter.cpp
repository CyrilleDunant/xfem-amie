
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//

#include "../physics/stiffness.h"
#include "../physics/stiffness_with_imposed_deformation.h"
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

     // predefined aggregate behaviour
     // use predefined space-time behaviour
     else if(dynamic_cast<AggregateBehaviour *>(behaviour))
     {
         if(dynamic_cast<ElasticOnlyAggregateBehaviour *>(behaviour))
         {
             ViscoElasticOnlyAggregateBehaviour * converted = new ViscoElasticOnlyAggregateBehaviour( 59e9, 0.3, (behaviour->getTensor( p ).numCols() == 3) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL ) ;
             converted->param = behaviour->getTensor( p ) ;
             converted->freeblocks = maxBlocks-3 ;
             f->setBehaviour( converted ) ;
             return true ;
         }
         else if(dynamic_cast<ViscoElasticOnlyAggregateBehaviour *>(behaviour))
         {
             dynamic_cast<ViscoElasticOnlyAggregateBehaviour *>(behaviour)->freeblocks = maxBlocks-3 ;
             return true ;
         }
         else if(dynamic_cast<ViscoDamageAggregateBehaviour *>(behaviour))
         {
             dynamic_cast<ViscoDamageAggregateBehaviour *>(behaviour)->freeblocks = maxBlocks-3 ;
             return true ;
         }
         ViscoDamageAggregateBehaviour * converted = new ViscoDamageAggregateBehaviour( 59e9, 0.3, dynamic_cast<AggregateBehaviour *>(behaviour)->up, 0.00008, (behaviour->getTensor( p ).numCols() == 3) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL ) ;
         converted->param = behaviour->getTensor( p ) ;
         converted->freeblocks = maxBlocks-3 ;
         f->setBehaviour( converted ) ;
         return true ;
     }

     // predefined cement paste behaviour
     // use predefined space-time behaviour
     else if(dynamic_cast<PasteBehaviour *>(behaviour))
     {
         if(dynamic_cast<ElasticOnlyPasteBehaviour *>(behaviour))
         {
             ViscoElasticOnlyPasteBehaviour * converted = new ViscoElasticOnlyPasteBehaviour( 12e9, 0.2, 0.3, 0.37, (behaviour->getTensor( p ).numCols() == 3) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL ) ;
             converted->param = behaviour->getTensor( p ) ;
             converted->freeblocks = maxBlocks-3 ;
             f->setBehaviour( converted ) ;
             return true ;
         }
         else if(dynamic_cast<ViscoElasticOnlyPasteBehaviour *>(behaviour))
         {
             dynamic_cast<ViscoElasticOnlyPasteBehaviour *>(behaviour)->freeblocks = maxBlocks-3 ;
             return true ;
         }
         else if(dynamic_cast<ViscoDamagePasteBehaviour *>(behaviour))
         {
             dynamic_cast<ViscoDamagePasteBehaviour *>(behaviour)->freeblocks = maxBlocks-3 ;
             return true ;
         }
         ViscoDamagePasteBehaviour * converted = new ViscoDamagePasteBehaviour( 12e9, 0.3, dynamic_cast<PasteBehaviour *>(behaviour)->up, 0.00015 , (behaviour->getTensor( p ).numCols() == 3) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL ) ;
         converted->param = behaviour->getTensor( p ) ;
         converted->freeblocks = maxBlocks-3 ;
         f->setBehaviour( converted ) ;
         return true ;
     }

     // predefined ASR gel paste behaviour
     // use predefined space-time behaviour
     else if(dynamic_cast<GelBehaviour *>(behaviour))
     {
         if(dynamic_cast<ViscoElasticOnlyGelBehaviour *>(behaviour))
         {
             dynamic_cast<ViscoElasticOnlyGelBehaviour *>(behaviour)->freeblocks = maxBlocks-3 ;
             return true ;
         }
         ViscoElasticOnlyGelBehaviour * converted = new ViscoElasticOnlyGelBehaviour( 22e9, 0.18, dynamic_cast<GelBehaviour*>(behaviour)->imposed[0], (behaviour->getTensor( p ).numCols() == 3) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL ) ;
         converted->param = behaviour->getTensor( p ) ;
         converted->freeblocks = maxBlocks-3 ;
         f->setBehaviour( converted ) ;
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
