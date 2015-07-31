/* this is an auto-generated file created on 31/6/2015 at 17:27  */

#ifndef __OBJECT_TRANSLATOR_H__
#define __OBJECT_TRANSLATOR_H__

#include "../physics/material_laws/material_laws.h"

namespace Amie
{

struct Object
{

    // parsed from header file: ../physics/material_laws/material_laws.h
    static ExternalMaterialLaw * getExternalMaterialLaw(std::string type, std::map<std::string, std::string> & strings, std::map<std::string, std::vector<std::string>> & stringlists, std::map<std::string, double> & values) ;
    static void resetExternalMaterialLaw(ExternalMaterialLaw * target) ;

} ;

}

#endif // __OBJECT_TRANSLATOR_H__
