// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "config_parser.h"
#include "function_parser.h"
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

ConfigTreeItem * ConfigParser::readXMLFile(std::string f, ConfigTreeItem * def, bool define, bool bind, std::vector<std::string> flags, std::string path )
{
    ConfigXMLParser parser(f) ;
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

size_t findNext( std::string test, char lookup, size_t start = 0, std::string ignore = std::string() )
{
    size_t i = start ;
    bool open = false ;
    while(i < test.length() )
    {
        char current = test[i] ;
        if(ignore.length() > 0 && ignore.find(current) < std::string::npos)
            open = !open ;
        if(!open && current == lookup)
            return i ;
        i++ ;
    }
    return std::string::npos ;
}

std::vector<std::string> ConfigXMLParser::breakLine( std::string line, std::string ignore ) 
{
    std::vector<std::string> ret ;
    std::string current ;
    size_t start = line.find('<') ;
    bool quotes = false ;
    bool in = true ;
    for(size_t i = start ; i < line.length() ; i++)
    {
        current += line[i] ;
        if( ignore.length() > 0 && ignore.find(line[i]) < std::string::npos )
        {
            quotes = !quotes ;
        }
        if( in && (line[i] == '>') && !quotes )
        {
            ret.push_back(current) ;
            current = std::string() ;
            in = false ;
        }
        if( !in && (line[i] == '<') && !quotes )
        {
            ret.push_back(current.substr(0,current.length()-1)) ;
            current = "<" ;
            in = true ;
        }

    }
    ret.push_back(current) ;

    return ret ;

}

XMLTokenType ConfigXMLParser::getTokenType( std::string test ) 
{
    if( test.find("<!--") == 0 && test.rfind("-->")+3 == test.length() && test.length() > 6) { return XML_INLINE_COMMENT ; }
    if( test.find("<!--") == 0) { return XML_OPEN_COMMENT ; }
    if( test.rfind("-->")+3 == test.length() && test.length() > 3) { return XML_CLOSE_COMMENT ; }
    if( test.find("<") == 0 && test.rfind(">")+1 == test.length() && test.length() > 2)
    {
        if( test.find("?") == 1 && test.rfind("?")+2 == test.length() && test.length() > 4) { return XML_HEADER ; }
        if( test.find("/") == 1  && test.length() > 3) { return XML_CLOSE ; }
        if( test.rfind("/")+2 == test.length() && test.length() > 3) { return XML_INLINE ; }
        return XML_OPEN ;
    }
    if(test.find("<") == 0 || test.rfind(">")+1 == test.length()) { return XML_INVALID ; }
    return XML_VALUE ;
}

std::string ConfigXMLParser::cropToken(std::string token, XMLTokenType type ) 
{
    switch(type)
    {
        case XML_OPEN:
            return token.substr(1, token.length()-2) ;   
        case XML_CLOSE:
            return token.substr(2, token.length()-3) ;   
        case XML_INLINE:
            return token.substr(1, token.length()-3) ;   
        case XML_VALUE:
        {
            std::string ret ;
            size_t start = 0 ;
            while( token[start] == ' ' && start < token.size() )
                start++ ;
            if(start == token.size())
                return std::string() ;
            size_t end = token.length()-1 ;
            while( end > 0 && token[end] == ' ')
                end-- ;
            return token.substr( start, end+1-(start) ) ;
        }
        default:
            break ;
    }
    return token ;
}


void ConfigXMLParser::readData()
{
    ConfigTreeItem * current = trunk ;
    std::string buffer ;
    std::string lastLine ;
    int lastLineNumber = 0 ;
    bool comment = false ;
    if(!file.fail())
    {
        for( std::string line ; getline( file, line ) ; )
        {
            lastLineNumber++ ;
            std::vector<std::string> tokens = ConfigXMLParser::breakLine( line ) ;
            for(size_t i = 0 ; i < tokens.size() ; i++)
            {
                if( tokens[i].length() == 0) { continue ; }
                XMLTokenType type = ConfigXMLParser::getTokenType( tokens[i] ) ;
                if( type == XML_OPEN_COMMENT ) { comment = true ; }
                if(!comment)
                {
                    switch(type)
                    {
                        case XML_CLOSE_COMMENT:
                        {
                            std::cout << "XML parsing error on line " << lastLineNumber << ": extra \"-->\" found" << std::endl ;
                            exit(0) ;
                            break ;
                        }
                        case XML_INVALID:
                        {
                            std::cout << "XML parsing error on line " << lastLineNumber << ": invalid token " << tokens[i] << std::endl ;
                            exit(0) ;
                            break ;
                        }
                        case XML_OPEN:
                        {
                            current->addChild( parseXMLToken( ConfigXMLParser::cropToken( tokens[i], XML_OPEN ) ) ) ;
                            current = current->getLastChild() ;
                            break ;
                        }
                        case XML_INLINE:
                        {
                            current->addChild( parseXMLToken( ConfigXMLParser::cropToken( tokens[i], XML_INLINE ) ) ) ;
                            break ;
                        }
                        case XML_CLOSE:
                        {
                            if( current->getLabel() != ConfigXMLParser::cropToken( tokens[i], XML_CLOSE ) )
                            {
                                std::cout << "XML parsing error on line " << lastLineNumber << ": unmatching tokens " << tokens[i] << std::endl ;
                                exit(0) ;
                            }
                            current = current->getFather() ;
                            break ;
                        }
                        case XML_VALUE:
                        {
                            std::string realToken = ConfigXMLParser::cropToken( tokens[i], XML_VALUE ) ;
                            if(realToken.length() > 0)
                            {
                                if( i+1 >= tokens.size() || i==0 )
                                {
                                    std::cout << "XML parsing error on line " << lastLineNumber << ": invalid value token " << tokens[i] << std::endl ;
                                    exit(0) ;
                                }
                                bool isDouble = (realToken.find_first_not_of("0123456789.e-") == std::string::npos ) ;
                                if(isDouble)
                                    current->setData( std::atof( realToken.c_str() ) ) ;
                                else
                                    current->setStringData( realToken ) ;
                            }
                            break ;
                        }
                        default:
                            break ;
                    }
                }
                if( type == XML_CLOSE_COMMENT ) { comment = false ; }
            }
        }
    }

    trunk->printTree() ;
}

ConfigTreeItem * ConfigXMLParser::parseXMLToken( std::string xml ) 
{
    std::vector<std::string> test = FunctionParser::breakString( xml, "\"\'" ) ;
    if(test.size() == 0)
        return nullptr ;

    ConfigTreeItem * ret = new ConfigTreeItem( nullptr, test[0] ) ;
    for(size_t i = 1 ; i < test.size() ; i++)
    {
        ret->addChild(parseXMLAttribute( test[i] ) );   
    }
    return ret ;
}

ConfigTreeItem * ConfigXMLParser::parseXMLAttribute( std::string attr ) 
{
    size_t sep = attr.find('=') ;
    std::string label = attr.substr(0, sep) ;
    size_t end = attr.rfind('\"') ;
    if( end == std::string::npos )
        end = attr.rfind('\'') ;
    std::string right = attr.substr(sep+2, end-(sep+2)) ;
    bool isDouble = (right.find_first_not_of("0123456789.e-") == std::string::npos ) ;
    if(isDouble)
        return new ConfigTreeItem( nullptr, label, atof(right.c_str()) ) ;
    return new ConfigTreeItem( nullptr, label, right ) ;
}
