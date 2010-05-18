//
// C++ Interface: mechanical homogenization
//
// Description: 
//
//
// Author:  Alain Giorla, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef OCTAVE_MANAGER_H
#define OCTAVE_MANAGER_H

#include "../../utilities/matrixops.h"
#include <fstream>

namespace Mu
{

/**
* Helper class to write files that the can be launched in octave (or similar softwares)
*/
class OctaveManager : public std::fstream
{
protected:
	std::fstream out ;
public:
	OctaveManager(const char * filename, ios_base::openmode mode) ;

	void writeArray(std::string var, Vector val) ;
	void writePlot(std::string x, std::string y, int fig) ;
	void close() { out.close() ; } ;

} ;



} ;


#endif

