// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../main.h"
#include "../../features/features.h"
#include "../../physics/stiffness.h"
#include "../../physics/viscoelasticity_and_fracture.h"
#include "../../physics/logarithmic_creep_with_external_parameters.h"
#include "../../utilities/writer/triangle_writer.h"
//#include "../../physics/damagemodels/spacetimefiberbasedtensileshearlineardamage.h"
#include "../../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h"
#include "../../physics/damagemodels/spacetimebadisotropiclineardamage.h"
#include "../../physics/damagemodels/spacetimefiberbasedbilineardamage.h"
//#include "../../physics/fracturecriteria/spacetimemultisurfacefracturecriterion.h"
#include "../../physics/fracturecriteria/spacetimelimitsurfacefracturecriterion.h"
#include "../../physics/fracturecriteria/maxstrain.h"
#include "../../features/sample.h"
#include "../../features/polygonSample.h"
#include "../../features/inclusion.h"
#include "../../utilities/parser/config_parser.h"
#include "../../utilities/postprocessor.h"
#include "../../utilities/itoa.h"

#include <fstream>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include <cmath>
#include <typeinfo>
#include <limits>
#include <sys/time.h>
#define DEBUG


using namespace Amie ;

int main(int argc, char *argv[])
{
	ConfigXMLParser test("amie_sensitivity_2016_02_26.xml") ;
	test.readData() ;

	return 0 ;
}

