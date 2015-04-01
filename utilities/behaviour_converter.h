
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//

#ifndef __BEHAVIOUR_CONVERTER_H_
#define __BEHAVIOUR_CONVERTER_H_

#include "../features/feature_base.h"
#include "../elements/integrable_entity.h"


namespace Amie
{

/** Configuration atom for the configuration tree.
 *
 */
struct BehaviourConverter
{
    static bool toSpaceTimeBehaviour( Feature * f, int maxBlocks ) ;

} ;

}



#endif // __BEHAVIOUR_CONVERTER_H_
