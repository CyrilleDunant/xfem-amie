// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../utilities/configuration.h"
#include "../utilities/parser.h"

#ifdef HAVE_OMP
#include <omp>
#endif
#include <fstream>
#include <cmath>
#include <typeinfo>
#include <limits>
#include <sys/time.h>

using namespace Amie ;

int main(int argc, char *argv[])
{
    ConfigTreeItem * problem = ConfigParser::readFile("toto.ini", nullptr, false, false) ;
    std::vector<ConfigTreeItem *> children = problem->getAllChildren() ;
    for(size_t i = 0 ; i < children.size() ; i++)
        children[i]->getExternalMaterialLaw() ;


    return 0 ;
}
