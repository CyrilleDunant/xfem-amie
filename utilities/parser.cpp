// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "parser.h"
#include "../polynomial/vm_function_extra.h"

void InclusionParser::readData()
{
    if (!file.fail())
    {
        circs->clear();
        size_t numberOfCircles ;
        size_t index ;
        file >> numberOfCircles ;

        for (size_t i = 0 ; i < numberOfCircles ; i++)
        {
            double centerX, centerY, radius ;
            file >> index >> radius >> centerX >> centerY ;
#ifdef DEBUG
            std::cout << "index : " << index << ", r : " << radius << ", x : " << centerX << ", y : " << centerY << std::endl ;
#endif
            circs->push_back(Circle(radius, centerX, centerY )) ;
        }
    }
}

void SegmentsParser::readData()
{
    if(!file.fail())
    {
        segments->clear() ;
        size_t numberOfSegments ;
        size_t index ;
        file >> numberOfSegments ;

        for (size_t i = 0 ; i < numberOfSegments ; i++)
        {
            if( file.eof())
                return ;

            size_t numberOfHeads, numberOfPoints ;
            file >> index >> numberOfHeads >> numberOfPoints ;
            std::valarray<Point *> m(numberOfPoints) ;

            for ( size_t j = 0 ; j < numberOfPoints ; j++)
            {
                double x, y ;
                file >> x >> y ;
                m[j] = new Point(x,y) ;

            }
            segments->push_back(SegmentedLine(m)) ;
        }
    }
}

void PointParser::readData()
{
    if(!file.fail())
    {
        pointSets->clear() ;

        size_t numberOfPoints, numberOfPointSets ;
        size_t index ;
        file >> numberOfPointSets ;

        for (size_t i = 0 ; i < numberOfPointSets ; i++)
        {
            file >> numberOfPoints ;
            PointSet * m = new PointSet(numberOfPoints) ;
            for (size_t j = 0 ;  j < numberOfPoints ; j++)
            {
                if( file.eof())
                    return ;

                file >> index >> m->getPoint(j)->getX() >> m->getPoint(j)->getY() ;
#ifdef DEBUG
                std::cout << "index : " << index << ", x : "<< m->getPoint(j)->getX() << ", y : " << m->getPoint(j)->getY() << std::endl ;
#endif
            }
            pointSets->push_back(m) ;
        }
    }
}



void PointParser3D::readData()
{
    if(!file.fail())
    {
        pointSets->clear() ;

        size_t numberOfPoints;
        size_t index ;
        size_t indexPeriodic;
        file >> numberOfPoints ;
        for (size_t i = 0 ;  i < numberOfPoints ; i++)
        {
            if( file.eof())
                return ;

            double x,y,z ;
//  			file >>  index>> x >> y >> z ;
            file >> index >> indexPeriodic >> x >> y >> z;
            Point * p = new Point(x,y,z) ;
            p->getId() = index-1 ;
            pointSets->push_back(p) ;
            this->PeriodicIds.push_back(indexPeriodic-1) ;

#ifdef DEBUG
//   			std::cout << "index : " << index << ", x : "<< p->getX() << ", y : " << p->getY() << ", z : " << p->getZ() << std::endl ;
            std::cout << "index : " << index << " indexPeriodic: " << indexPeriodic << ", x : "<< p->getX() << ", y : " << p->getY() << ", z : " << p->getZ() << std::endl ;
#endif
        }
    }

}

void PeriodParser::readData()
{
    if(!file.fail())
    {
        periodSets->clear() ;

        size_t numberOfindex;

        file >> numberOfindex ;
        std::vector<Point *>  p;

        for (size_t i = 0 ;  i < numberOfindex ; i++)
        {
            if( file.eof())
                return ;

            std::valarray<int> index_id(1) ;
            size_t index ;
            file >> index >> index_id[0] ;

            periodSets->push_back(std::pair<size_t,std::valarray<int> >(index,index_id )) ;
#ifdef DEBUG
            std::cout << "index : " << index << "id "<< index_id[0]<< std::endl ;
#endif
        }
    }

}
void TetrahedronParser::readData()
{
    if(!file.fail())
    {
        tets->clear() ;
        size_t numberOfTets ;
        file >> numberOfTets ;
        for (size_t i = 0 ;  i < numberOfTets ; i++)
        {
            if( file.eof())
                return ;

            std::valarray<int> vertex_id(4) ;
            size_t mat_index ;
            size_t index ;
            file >> index >> vertex_id[0] >> vertex_id[1] >> vertex_id[2] >> vertex_id[3] >>  mat_index;
            vertex_id[0]-- ;
            vertex_id[1]-- ;
            vertex_id[2]-- ;
            vertex_id[3]-- ;

            tets->push_back(std::pair<std::valarray<int>, size_t>(vertex_id, mat_index)) ;
#ifdef DEBUG
            std::cout << "index : " << index << ", point 0 : "<< vertex_id[0]  << ", point 1 : " << vertex_id[1]  << ", point 2 : " << vertex_id[2]  << ", point 3 :" << vertex_id[3]  << ", material : "<< mat_index<< std::endl ;
#endif
        }
    }
}
void HexahedronParser::readData()
{
    if(!file.fail())
    {
        hex->clear() ;
        size_t numberOfHex ;
        file >> numberOfHex ;
        for (size_t i = 0 ;  i < numberOfHex ; i++)
        {
            if( file.eof())
                return ;

            std::valarray<int> vertex_id(8) ;
            size_t mat_index ;
            size_t index ;
            file >> index >> vertex_id[0] >> vertex_id[1] >> vertex_id[2] >> vertex_id[3] >> vertex_id[4] >> vertex_id[5] >> vertex_id[6]  >> vertex_id[7] >> mat_index;
            vertex_id[0]-- ;
            vertex_id[1]-- ;
            vertex_id[2]-- ;
            vertex_id[3]-- ;
            vertex_id[4]-- ;
            vertex_id[5]-- ;
            vertex_id[6]-- ;
            vertex_id[7]-- ;

            hex->push_back(std::pair<std::valarray<int>, size_t>(vertex_id, mat_index)) ;
#ifdef DEBUG
            std::cout << "index : " << index << ", point 0 : "<< vertex_id[0]  << ", point 1 : " << vertex_id[1]  << ", point 2 : " << vertex_id[2]  << ", point 3 :" << vertex_id[3]  << ", point 4 :" << vertex_id[4]  << ", point 5 :" << vertex_id[5]  << ", point 6 :" << vertex_id[6]  << ", point 7 :" << vertex_id[7]  << ", material : "<< mat_index<< std::endl ;
#endif
        }
    }
}
void BoundaryParser::readData()
{
    if(!file.fail())
    {
        tri0->clear() ;
        tri1->clear() ;
        tri2->clear() ;
        tri3->clear() ;
        tri4->clear() ;
        tri5->clear() ;
        size_t numberOfTris ;
        file >> numberOfTris ;
        for (size_t i = 0 ;  i < numberOfTris ; i++)
        {
            if( file.eof())
                return ;

            std::valarray<int> vertex_id(3) ;
            size_t face_index ;
            size_t index ;
            size_t to_discard ;
            file >> index >> vertex_id[0] >> vertex_id[1] >> vertex_id[2] >>  face_index >> to_discard ;
            vertex_id[0]-- ;
            vertex_id[1]-- ;
            vertex_id[2]-- ;
            if(face_index == 1)
            {
                tri0->push_back(vertex_id) ;
            }
            else if(face_index == 2)
            {
                tri1->push_back(vertex_id) ;
            }
            else if(face_index == 3)
            {
                tri2->push_back(vertex_id) ;
            }
            else if(face_index == 4)
            {
                tri3->push_back(vertex_id) ;
            }
            else if(face_index == 5)
            {
                tri4->push_back(vertex_id) ;
            }
            else if(face_index == 6)
            {
                tri5->push_back(vertex_id) ;
            }
            else
            {
                assert(false) ;
            }
#ifdef DEBUG
            std::cout << "index : " << index << ", point 0 : "<< vertex_id[0]  << ", point 1 : " << vertex_id[1]  << ", point 2 : " << vertex_id[2]  <<  ", face : "<< face_index<< std::endl ;
#endif
        }
    }
}

std::vector<std::valarray<int> > * BoundaryParser::getData(size_t face_index)
{
    if(face_index == 0)
    {
        return tri0 ;
    }
    else if(face_index == 1)
    {
        return tri1 ;
    }
    else if(face_index == 2)
    {
        return tri2 ;
    }
    else if(face_index == 3)
    {
        return tri3 ;
    }
    else if(face_index == 4)
    {
        return tri4 ;
    }
    else if(face_index == 5)
    {
        return tri5 ;
    }
    else
    {
        assert(false) ;
        return nullptr ; //shut up the compiler ;
    }
}

void HexahedronBoundaryParser::readData()
{
    if(!file.fail())
    {
        carre0->clear() ;
        carre1->clear() ;
        carre2->clear() ;
        carre3->clear() ;
        carre4->clear() ;
        carre5->clear() ;
        size_t numberOfCarres ;
        file >> numberOfCarres ;
        for (size_t i = 0 ;  i < numberOfCarres ; i++)
        {
            if( file.eof())
                return ;

            std::valarray<int> vertex_id(4) ;
            size_t face_index ;
            size_t index ;
            size_t to_discard ;
            file >> index >> vertex_id[0] >> vertex_id[1] >> vertex_id[2] >> vertex_id[3]>>  face_index >> to_discard ;
            vertex_id[0]-- ;
            vertex_id[1]-- ;
            vertex_id[2]-- ;
            vertex_id[3]-- ;
            if(face_index == 1)
            {
                carre0->push_back(vertex_id) ;
            }
            else if(face_index == 2)
            {
                carre1->push_back(vertex_id) ;
            }
            else if(face_index == 3)
            {
                carre2->push_back(vertex_id) ;
            }
            else if(face_index == 4)
            {
                carre3->push_back(vertex_id) ;
            }
            else if(face_index == 5)
            {
                carre4->push_back(vertex_id) ;
            }
            else if(face_index == 6)
            {
                carre5->push_back(vertex_id) ;
            }
            else
            {
                assert(false) ;
            }
#ifdef DEBUG
            std::cout << "index : " << index << ", point 0 : "<< vertex_id[0]  << ", point 1 : " << vertex_id[1]  << ", point 2 : " << vertex_id[2] << ", point 3 : " << vertex_id[3] <<  ", face : "<< face_index<< std::endl ;
#endif
        }
    }
}

std::vector<std::valarray<int> > * HexahedronBoundaryParser::getData(size_t face_index)
{
    if(face_index == 0)
    {
        return carre0 ;
    }
    else if(face_index == 1)
    {
        return carre1 ;
    }
    else if(face_index == 2)
    {
        return carre2 ;
    }
    else if(face_index == 3)
    {
        return carre3 ;
    }
    else if(face_index == 4)
    {
        return carre4 ;
    }
    else if(face_index == 5)
    {
        return carre5 ;
    }
    else
    {
        assert(false) ;
        return nullptr ; //shut up the compiler ;
    }
}

int ConfigParser::getIndentLevel( std::string test )
{
    size_t i = 0 ;
    while(test[i] == '.' && i < test.size())
        i++ ;
    return i ;
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
//				std::cout << line << std::endl ;
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
//				std::cout << line.substr(level) << std::endl ;
                new ConfigTreeItem( current , line.substr(level) ) ;
//				current->printTree() ;
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

ConfigTreeItem * ConfigParser::readFile(std::string f, ConfigTreeItem * def, bool define)
{
    ConfigParser parser(f) ;
    parser.readData() ;
    ConfigTreeItem * ret = parser.getData() ;
    if(ret->hasChild("template"))
        ret = ret->getChild("template")->makeTemplate() ;
    if(define)
        ret->define(def, true) ;
    if(ret->hasChildFromFullLabel("define.path"))
        ret->definePath( ret->getStringData("define.path", "./") ) ;

    return ret ;
}

CommandLineParser::CommandLineParser(std::string d, bool c) : description(d), commandLineConfiguration(c)
{
    addFlag("--help", false, "print help") ;
    addFlag("--version", false, "print current AMIE revision") ;
    addFlag("--no-openmp", false, "disable OpenMP (equivalent to --set-num-threads 1)") ;
    addValue("--set-num-threads", -1, "set the number of threads available for OpenMP") ;
}

void CommandLineParser::parseCommandLine( int argc, char *argv[] )
{
	command = std::string(argv[0]) ;

	config = nullptr ;
	size_t i = 1 ;
	while(i < (size_t) argc)
	{
		std::string test = std::string(argv[i]) ;
		if(i-1 < arguments.size() && test[0] != '-')
		{
			arguments[i-1].str = test ;
			arguments[i-1].val = atof(test.c_str()) ;
		}
		else if( flags.find( test ) != flags.end() )
			flags[ test ] = true ;
		else if( values.find( test ) != values.end()  )
		{
			values[ test ] = atof(argv[i+1]) ;
			i++ ;
		}
		else if( strings.find( test ) != strings.end()  )
		{
			strings[ test ] = std::string(argv[i+1]) ;
			i++ ;
		}
		else if(commandLineConfiguration && test[0] == '@')
		{
			if(!config)
				config = new ConfigTreeItem(nullptr, "define") ;
			std::string testval = std::string(argv[i+1]) ;
			bool isDouble = (testval.find_first_not_of("0123456789.e-") == std::string::npos ) ;
			if(isDouble)
				new ConfigTreeItem( config, test, atof(testval.c_str()) ) ;
			else
				new ConfigTreeItem( config, test, testval ) ;
			i++ ;
		}
		else if(commandLineConfiguration && test[0] == '.')
		{
			directConfig [test.substr(1)] = std::string (argv[i+1]) ;
			i++ ;
		}
		i++ ;
	}

	if( getFlag("--version") )
		printVersion() ;
	if( getFlag("--help") )
		printHelp() ;

        if( getFlag("--no-openmp"))
		setNumThreads(1) ;
	else if( getValue("--set-num-threads") > 0)
		setNumThreads( getValue("--set-num-threads") ) ;

	if(getFlag("--help") || getFlag("--version"))
		exit(0) ;

}

void CommandLineParser::setNumThreads( int n )
{
#ifdef HAVE_OMP
	omp_set_num_threads(std::max(n,1)) ;
#endif
}

void CommandLineParser::printStatus( )
{
	for(auto f = flags.begin() ; f != flags.end() ; f++)
		std::cout << f->first << " = " << (f->second ? "TRUE" : "FALSE") << std::endl ;
	for(auto f = values.begin() ; f != values.end() ; f++)
		std::cout << f->first << " = " << f->second << std::endl ;

}

void printFormatedHelp( std::string arg, size_t max, std::string help, std::string lead)
{
	std::cout << arg ;
	for(size_t i = 0 ; i < max-arg.length() ; i++)
		std::cout << " " ;
	std::cout << lead << help << std::endl ;
}

void CommandLineParser::printHelp( )
{
	std::cout << std::endl ;

	if(description.length() > 0)
	{
		std::cout << description << std::endl ;
		std::cout << std::endl ;
	}

	std::cout << "USAGE" << std::endl ;
	std::cout << command << "  " ;
	for(size_t i = 0 ; i < arguments.size() ; i++)
		std::cout << arguments[i].name << "  " ;
	if(commandLineConfiguration)
		std::cout << "[@argument (value) ...]  " ;
	std::cout << "[options...]" << std::endl ;

	std::cout << std::endl ;
	if(arguments.size() > 0)
		std::cout << "ARGUMENTS" << std::endl ;

	size_t m = 0 ;
	for(size_t i = 0 ; i < arguments.size() ; i++)
		m = std::max(m, arguments[i].name.length() ) ;

	for(size_t i = 0 ; i < arguments.size() ; i++)
		printFormatedHelp( arguments[i].name, m, arguments[i].help, "  " ) ;
	if(arguments.size() > 0)
		std::cout << std::endl ;

	std::cout << "OPTIONS" << std::endl ;

	m = 0 ;
	for(auto h = help.begin() ; h != help.end() ; h++)
		m = std::max(m, h->first.length() ) ;

	for(auto h = help.begin() ; h != help.end() ; h++)
	{
		std::string lead = "  " ;
		if(values.find(h->first) != values.end())
			lead = "  (numeral) " ;
		if(strings.find(h->first) != strings.end())
			lead = "  (string) " ;
		printFormatedHelp( h->first, m, h->second, lead ) ;
	}

	std::cout << std::endl ;
}

void CommandLineParser::printVersion( )
{
	std::cout << "current AMIE version: svn " << std::flush ;
	std::system("svnversion ..") ;
}

bool CommandLineParser::getFlag( std::string f )
{
	if( flags.find( f ) == flags.end() )
		return false ;
	return flags[f] ;
}

double CommandLineParser::getValue( std::string f )
{
	if( values.find( f ) == values.end() )
		return 0. ;
	return values[f] ;
}

std::string CommandLineParser::getString( std::string f )
{
	if( strings.find( f ) == strings.end() )
		return std::string() ;
	return strings[f] ;
}

std::string CommandLineParser::getStringArgument( size_t i )
{
	if(i >= arguments.size())
		return std::string() ;
	return arguments[i].str ;
}

std::string CommandLineParser::getStringArgument( std::string arg )
{
	for(size_t i = 0 ; i < arguments.size() ; i++)
	{
		if(arguments[i].name == arg)
			return arguments[i].str ;
	}
	return std::string() ;
}

double CommandLineParser::getNumeralArgument( size_t i )
{
	if(i >= arguments.size())
		return 0. ;
	return arguments[i].val ;
}

double CommandLineParser::getNumeralArgument( std::string arg )
{
	for(size_t i = 0 ; i < arguments.size() ; i++)
	{
		if(arguments[i].name == arg)
			return arguments[i].val ;
	}
	return 0. ;
}

FunctionParser::FunctionParser(std::vector<std::string> f) : isLinking(false)
{
	renewExpression(f) ;
}


FunctionParser::FunctionParser(std::string f) : isLinking(false)
{
	std::string current ;
	std::vector<std::string> str ;
	size_t i = 0 ;
	while(i < f.length())
	{
		bool cut = false ;
		char test = f[i] ;
		if(test == ' ')
			cut = true ;
		else
			current += test ;
		if(cut)
		{
			if(current.length() > 0)
			{
				str.push_back( current ) ;
				current.clear() ;
			}
		}
		i++ ;
	}
	if(current.length() > 0)
		str.push_back( current ) ;
	renewExpression(str) ;
}

bool isOpenBracket( std::string s) 
{
	return s== "(" ;
}

bool isCloseBracket(std::string s) 
{
	return s==")" ;
}

void FunctionParser::renewExpression(std::vector<std::string> str) 
{
	size_t i = 0 ;
	std::vector<std::string> embeddedExpression ;
	size_t openCount = 0 ;
	
	while(i < str.size())
	{
		if(openCount == 0 && !isOpenBracket( str[i] )  )
			tokens.push_back( new FunctionParserToken( str[i] ) ) ;
		else
		{
			if( isOpenBracket( str[i] ) )
			{
				if(openCount == 0)
					embeddedExpression.clear() ;
				else
					embeddedExpression.push_back( str[i] ) ;
				openCount++ ;
			}
			else if(isCloseBracket( str[i] ) )
			{
				openCount-- ;
				if(openCount == 0)
					tokens.push_back( new FunctionParser( embeddedExpression) ) ;
				else
					embeddedExpression.push_back( str[i] ) ;
			}
			else
				embeddedExpression.push_back( str[i] ) ;
		}
		i++ ;
	}

	this->link() ;
}

void FunctionParser::link()
{
	if(isLinking || isFinalToken() )
		return ;

	isLinking=true ;

	roots.resize(tokens.size(), -1) ;

	for(size_t i = 0 ; i < tokens.size() ; i++)
	{
		if(tokens[i]->isFinalToken())
		{
			TokenOperationType elem = dynamic_cast<FunctionParserToken *>(tokens[i])->toFunctionToken() ;
			switch(elem)
			{
				case TOKEN_OPERATION_ABS:
				case TOKEN_OPERATION_COS:
				case TOKEN_OPERATION_SIN:
				case TOKEN_OPERATION_TAN:
				case TOKEN_OPERATION_COSH:
				case TOKEN_OPERATION_SINH:
				case TOKEN_OPERATION_TANH:
				case TOKEN_OPERATION_EXP:
				case TOKEN_OPERATION_SIGN:
				case TOKEN_OPERATION_POSITIVITY:
				case TOKEN_OPERATION_NEGATIVITY:
				case TOKEN_OPERATION_LOG:
				case TOKEN_OPERATION_SQRT:
					roots[i+1] = i ;
					break ;
				default:
					break ;
			}
		}
//		else
//			tokens[i]->link() ;
	}

	for(size_t i = 0 ; i < tokens.size() ; i++)
	{
		if(tokens[i]->isFinalToken())
		{
			if( dynamic_cast<FunctionParserToken *>(tokens[i])->toFunctionToken() == TOKEN_OPERATION_POWER )
			{
				size_t left = i-1 ;
				size_t c = 0 ;
				while(roots[left] > -1 && c < tokens.size())
				{
					left = roots[left] ;
					c++;
				}
				size_t right = i+1 ;
				c = 0 ;
				while(roots[right] > -1 && c < tokens.size())
				{
					right = roots[right] ;
					c++;
				}
				roots[left] = i ;
				roots[right] = i ;
			}
		}
	}

	for(size_t i = 0 ; i < tokens.size() ; i++)
	{
		if(tokens[i]->isFinalToken())
		{
			if( dynamic_cast<FunctionParserToken *>(tokens[i])->toFunctionToken() == TOKEN_OPERATION_TIMES )
			{
				size_t left = i-1 ;
				size_t c = 0 ;
				while(roots[left] > -1 && c < tokens.size())
				{
					left = roots[left] ;
					c++;
				}
				size_t right = i+1 ;
				c = 0 ;
				while(roots[right] > -1 && c < tokens.size())
				{
					right = roots[right] ;
					c++;
				}
				roots[left] = i ;
				roots[right] = i ;
			}
		}
	}

	for(size_t i = 0 ; i < tokens.size() ; i++)
	{
		if(tokens[i]->isFinalToken())
		{
			if( dynamic_cast<FunctionParserToken *>(tokens[i])->toFunctionToken() == TOKEN_OPERATION_DIVIDES )
			{
				size_t left = i-1 ;
				size_t c = 0 ;
				while(roots[left] > -1 && c < tokens.size())
				{
					left = roots[left] ;
					c++;
				}
				size_t right = i+1 ;
				c = 0 ;
				while(roots[right] > -1 && c < tokens.size())
				{
					right = roots[right] ;
					c++;
				}
				roots[left] = i ;
				roots[right] = i ;
			}
		}
	}

	for(size_t i = 0 ; i < tokens.size() ; i++)
	{
		if(tokens[i]->isFinalToken())
		{
			if( dynamic_cast<FunctionParserToken *>(tokens[i])->toFunctionToken() == TOKEN_OPERATION_MINUS )
			{
				size_t left = i-1 ;
				size_t c = 0 ;
				while(roots[left] > -1 && c < tokens.size())
				{
					left = roots[left] ;
					c++;
				}
				size_t right = i+1 ;
				c = 0 ;
				while(roots[right] > -1 && c < tokens.size())
				{
					right = roots[right] ;
					c++;
				}
				roots[left] = i ;
				roots[right] = i ;
			}
		}
	}

	for(size_t i = 0 ; i < tokens.size() ; i++)
	{
		if(tokens[i]->isFinalToken())
		{
			if( dynamic_cast<FunctionParserToken *>(tokens[i])->toFunctionToken() == TOKEN_OPERATION_PLUS )
			{
				size_t left = i-1 ;
				size_t c = 0 ;
				while(roots[left] > -1 && c < tokens.size())
				{
					left = roots[left] ;
					c++;
				}
				size_t right = i+1 ;
				c = 0 ;
				while(roots[right] > -1 && c < tokens.size())
				{
					right = roots[right] ;
					c++;
				}
				roots[left] = i ;
				roots[right] = i ;
			}
		}
	}

//	for(size_t i = 0 ; i < tokens.size() ; i++)
//		std::cout << roots[i] << std::endl ;

	isLinking = false ;

}

int FunctionParser::getLeftTokenIndex( size_t i ) const 
{
	for(size_t j = 0 ; j < roots.size() ; j++)
	{
		if(roots[j] == (int) i)
			return j ;
	}
	return -1 ;
}

int FunctionParser::getRightTokenIndex( size_t i ) const 
{
	for(size_t j = 0 ; j < roots.size() ; j++)
	{
		if(roots[roots.size()-1-j] == (int) i)
			return roots.size()-1-j ;
	}
	return -1 ;
}

Function FunctionParser::getLeftFunction( size_t i ) const 
{
	int j = getLeftTokenIndex(i) ;
	if(j >= 0)
		return getFunction(j) ;
	return Function("0") ;
}

Function FunctionParser::getRightFunction( size_t i ) const 
{
	int j = getRightTokenIndex(i) ;
	if(j >= 0)
		return getFunction(j) ;
	return Function("0") ;
}

Function FunctionParser::getFunction() const
{
	for(size_t i = 0 ; i < roots.size() ; i++)
	{
		if(roots[i] == -1)
		{
			return getFunction(i) ;
		}
	}
	return Function("0") ;
}

Function FunctionParser::getFunction( size_t i ) const 
{
	if( tokens[i]->isFinalToken())
	{
		FunctionParserToken * local = dynamic_cast<FunctionParserToken *>( tokens[i] ) ; 
		TokenOperationType type = local->toFunctionToken() ;
		switch(type)
		{
			case TOKEN_OPERATION_CONSTANT:
				return Function( (local->data).c_str() ) ;
			case TOKEN_OPERATION_X:
				return Function("x") ;
			case TOKEN_OPERATION_Y:
				return Function("y") ;
			case TOKEN_OPERATION_Z:
				return Function("z") ;
			case TOKEN_OPERATION_T:
				return Function("t") ;
			case TOKEN_OPERATION_U:
				return Function("u") ;
			case TOKEN_OPERATION_V:
				return Function("v") ;
			case TOKEN_OPERATION_W:
				return Function("w") ;
			case TOKEN_OPERATION_PLUS:
				return getLeftFunction(i)+getRightFunction(i) ;
			case TOKEN_OPERATION_MINUS:
				return getLeftFunction(i)-getRightFunction(i) ;
			case TOKEN_OPERATION_TIMES:
				return getLeftFunction(i)*getRightFunction(i) ;
			case TOKEN_OPERATION_DIVIDES:
				return getLeftFunction(i)/getRightFunction(i) ;
			case TOKEN_OPERATION_POWER:
				return getLeftFunction(i)^getRightFunction(i) ;
			case TOKEN_OPERATION_ABS:
				return f_abs( getLeftFunction(i) ) ;
			case TOKEN_OPERATION_COS:
				return f_cos( getLeftFunction(i) ) ;
			case TOKEN_OPERATION_SIN:
				return f_sin( getLeftFunction(i) ) ;
			case TOKEN_OPERATION_TAN:
				return f_tan( getLeftFunction(i) ) ;
			case TOKEN_OPERATION_COSH:
				return f_cosh( getLeftFunction(i) ) ;
			case TOKEN_OPERATION_SINH:
				return f_sinh( getLeftFunction(i) ) ;
			case TOKEN_OPERATION_TANH:
				return f_tanh( getLeftFunction(i) ) ;
			case TOKEN_OPERATION_EXP:
				return f_exp( getLeftFunction(i) ) ;
			case TOKEN_OPERATION_SIGN:
				return f_sign( getLeftFunction(i) ) ;
			case TOKEN_OPERATION_POSITIVITY:
				return f_positivity( getLeftFunction(i) ) ;
			case TOKEN_OPERATION_NEGATIVITY:
				return f_negativity( getLeftFunction(i) ) ;
			case TOKEN_OPERATION_LOG:
				return f_log( getLeftFunction(i) ) ;
			case TOKEN_OPERATION_SQRT:
				return f_sqrt( getLeftFunction(i) ) ;
			default:
				return Function("0") ;
		}
	}
	return tokens[i]->getFunction() ;
		
}

Function FunctionParser::getFunction( std::string f)
{
	FunctionParser expr(f) ;
	return expr.getFunction() ;
}

void FunctionParser::print(bool end) const
{
	std::cout << "( " ;
	for(size_t i = 0 ; i < tokens.size() ; i++)
		tokens[i]->print(false) ;
	std::cout << " ) " ;
	if(end)
		std::cout << std::endl ;
}

void FunctionParserToken::print(bool end) const
{
	std::cout << " " << data << " " ;
	if(end)
		std::cout << std::endl ;
}

TokenOperationType FunctionParserToken::toFunctionToken() const
{
	Function f ;
	std::vector<double> v ;
	return f.toToken( data, 0, v ).first ;
}
