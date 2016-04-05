// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "config_parser.h"
#include "../font.h"
#include "../enumeration_translator.h"
#include "../../polynomial/vm_function_extra.h"
#include <sys/stat.h>

using namespace Amie ;

int ConfigParser::getIndentLevel( std::string test )
{
    size_t i = 0 ;
    while(test[i] == '.' && i < test.size())
        i++ ;
    return i ;
}

Form * ConfigParser::getBehaviour( std::string filename, Form * def , SpaceDimensionality dim) 
{
	ConfigTreeItem * cnf = ConfigParser::readFile( filename, nullptr ) ;
	if(cnf->hasChild("behaviour"))
	{
		Form * b = cnf->getChild("behaviour")->getBehaviour(dim) ;
		if(b)
			return b ;
	}
	return def ;
}

void ConfigParser::readData()
{
    ConfigTreeItem * current = trunk ;
    int level = 1 ;
    std::string buffer ;
    std::string lastLine ;
    int lastLineNumber = 0 ;
    if(!file.fail())
    {
        for( std::string line ; getline( file, line ) ; )
        {
            if(line.length() > 0)
            {
                lastLineNumber++ ;
                lastLine = line ;
            }
            if(line[0] == '#')
                continue ;

            size_t comment = line.find('#') ;
            if(comment != std::string::npos)
            {
                line = line.substr(0, comment-1) ;
            }

            size_t found = line.find(" ") ;
            size_t quote = line.find('"') ;
            while(found != std::string::npos && found < quote)
            {
                line = line.erase(found,1) ;
                found = line.find(" ") ;
            }
            quote = line.find('"') ;
            while(quote != std::string::npos)
            {
                line = line.erase(quote,1) ;
                quote = line.find('"') ;
            }

            size_t sep = line.find("=") ;
            int l = ConfigParser::getIndentLevel(line) ;
            if(l == level)
            {
                // current = current ;
            }
            else if(l == level+1)
            {
                current = current->getLastChild() ;
                level++ ;
            }
            else if(l == 0 && line.length() == 0)
                break ;
            else if(l == 0 || l > level+1)
            {
                std::cout << "parsing error in file " << filename << std::endl ;
                current->print() ;
                exit(0) ;
            }
            else // (l between 1 and level-1)
            {
                while(level > l)
                {
                    level-- ;
                    current = current->getFather() ;
                }
            }
            if(sep == std::string::npos)
            {
                new ConfigTreeItem( current , line.substr(level) ) ;
            }
            else
            {
                std::string right = line.substr( sep+1 ) ;
                bool isDouble = (right.find_first_not_of("0123456789.e-") == std::string::npos ) ;
                std::string label = line.substr(0, sep) ;
                label = label.substr(level) ;
                if(isDouble)
                {
                    new ConfigTreeItem( current , label, atof(right.c_str()) ) ;
                }
                else
                {
                    /*					if(right.find("@") == 0 && authorizeIncludes)
                    					{
                    						right = right.substr(1) ;
                    						ConfigParser includes( right, false) ;
                    						includes.readData() ;
                    						ConfigTreeItem * cnf = includes.getData() ;
                    						ConfigTreeItem * child = cnf->getChild(label) ;
                    						current->addChild( child ) ;
                    						child->setFather( current ) ;
                    					}
                    					else if( label == "include"  && authorizeIncludes)
                    					{
                    						ConfigParser includes( right, false) ;
                    						includes.readData() ;
                    						ConfigTreeItem * cnf = includes.getData() ;
                    						std::vector<ConfigTreeItem *> children = cnf->getAllChildren() ;
                    						for(size_t i = 0 ; i < children.size() ; i++)
                    						{
                    							current->addChild( children[i] ) ;
                    							children[i]->setFather( current ) ;
                    						}
                    					}
                    					else*/
                    new ConfigTreeItem( current , label, right ) ;
                }
            }

        }

        std::cout << filename << " parsed with success!" << std::endl ;
        std::cout << "last line read: " << lastLineNumber << " >>>   " << lastLine << std::endl ;

    }
    else
        std::cout << filename << " not found!" << std::endl ;


#ifdef _WIN32
    if(trunk)
        trunk->makeWindowsPath() ;
#endif

}

std::vector<BoundaryCondition *> ConfigParser::getBoundaryConditions( std::string filename, FeatureTree * F) 
{
    std::vector<BoundaryCondition *> ret ;
    ConfigTreeItem * problem = ConfigParser::readFile(filename, nullptr, false, true) ;
    if(problem->hasChild("boundary_conditions"))
        ret = problem->getChild("boundary_conditions")->getAllBoundaryConditions(F) ;
    return ret ;
}

ConfigTreeItem * ConfigParser::readFile(std::string f, ConfigTreeItem * def, bool define, bool bind, std::vector<std::string> flags, std::string path )
{
    ConfigParser parser(f) ;
    parser.readData() ;
    ConfigTreeItem * ret = parser.getData() ;
    if(ret->hasChild("template"))
        ret = ret->getChild("template")->makeTemplate() ;
    if(define)
        ret->define(def, true) ;
    if(path.size() > 0)
        ret->definePath( path ) ;
    else if(ret->hasChildFromFullLabel("define.path"))
        ret->definePath( ret->getStringData("define.path", "./") ) ;

    if(bind)
    {
        std::string realpath = path ;
        if(path.size() == 0 && ret->hasChildFromFullLabel("define.path"))
            realpath = ret->getStringData("define.path", std::string()) ;
        std::vector<std::string> callers ;
        callers.push_back(f) ;
        ret->bindInput(callers, realpath, flags ) ;
    }

    return ret ;
}


