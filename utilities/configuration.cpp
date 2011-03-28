
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//

#include "configuration.h"

std::string ConfigAtom::translate(Type t)
{
	switch(t)
	{
	case TRUNK:
		return std::string("trunk") ;
	case DOUBLE:
		return std::string("double") ;
	case INCLUSION:
		return std::string("inclusion") ;
	case CRACK:
		return std::string("crack") ;
	case TYPE_PROBLEME:
		return std::string("problemType") ;
	default:
		return std::string("unknown") ;
	}
}

ConfigAtom * ConfigAtom::getParent() 
{
	return parent ; 
}
std::string ConfigAtom::getLabel() 
{
	return this->label ; 
}

void ConfigAtom::setLabel(std::string lab)
{
	if(type == TRUNK)
		return ;
	
	label = lab ;
}

std::string ConfigAtom::getFullLabel()
{
	// first we go back to the trunk ;
	
	std::string ret = "trunk" ;
	
	if (type == TRUNK)
		return ret ;
	
	ConfigAtom * current = this ;
	std::vector<std::string> familytree ;
	while(current->getParent()!=NULL)
	{
		familytree.insert(familytree.begin(), current->getLabel()) ;
		current = current->getParent() ;
	}

	// then we reconstruct the family tree
	for (size_t i = 0  ;  i < familytree.size() ; i++)
	{
		ret += "." + familytree[i] ;
	}
	
	return ret ;
}

