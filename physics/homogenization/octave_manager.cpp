//
// C++ Implementation: mechanical analytic homogenization
//
// Description:
//
//
// Author:  Alain Giorla, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "octave_manager.h"
#include <fstream>

using namespace Mu ;

OctaveManager::OctaveManager(const char * filename, ios_base::openmode mode)
{
	out.open(filename, mode) ;
}

void OctaveManager::writeArray(std::string var, Vector val)
{
	out << var << "= [ " ;
	for(size_t i = 0 ; i < val.size() -2 ; i++)
		out << val[i] << " , " ;
	out << val[val.size()-1] << " ] ;" << std::endl ;
}

void OctaveManager::writePlot(std::string x, std::string y, int fig) 
{
	out << "figure(" << fig << ") ;" << std::endl ;
	out << "plot(" << x << "," << y << ") ;" << std::endl ;
}

