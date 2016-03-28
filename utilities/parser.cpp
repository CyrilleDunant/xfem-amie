// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "parser.h"
#include "font.h"
#include "enumeration_translator.h"
#include "../polynomial/vm_function_extra.h"
#include <sys/stat.h>

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

void CommandLineParser::setFeatureTree( FeatureTree * f ) 
{
    if(values[std::string("--set-sampling-restriction")] > -1)
        f->setSamplingRestriction( values[std::string("--set-sampling-restriction")] ) ;
    if(values[std::string("--set-delta-time")] > 0)
        f->setDeltaTime( values[std::string("--set-delta-time")] ) ;
    if(values[std::string("--set-min-delta-time")] > 0)
        f->setMinDeltaTime( values[std::string("--set-min-delta-time")] ) ;
    if(values[std::string("--set-solver-precision")] > 0)
        f->setSolverPrecision( values[std::string("--set-solver-precision")] ) ;
    if(values[std::string("--set-max-iterations-per-step")] > 0)
        f->setMaxIterationsPerStep( values[std::string("--set-max-iterations-per-step")] ) ;
    if(values[std::string("--set-ssor-iterations")] > -1)
        f->setSSORIterations( values[std::string("--set-ssor-iterations")] ) ;
    if(values[std::string("--set-sampling-number")] > 0)
        f->setSamplingNumber( values[std::string("--set-sampling-number")] ) ;
    if(values[std::string("--set-surface-sampling-factor")] > 1-POINT_TOLERANCE)
        f->setSurfaceSamplingFactor( values[std::string("--set-surface-sampling-factor")] ) ;
    if(values[std::string("--set-min-mesh-density")] > -1)
        f->setMinimumMeshDensity( values[std::string("--set-min-mesh-density")] ) ;
    if(strings[std::string("--set-order")].size() > 0)
    {
        bool ok = false ;
        Order ord = Enum::getOrder(strings[std::string("--set-order")], &ok) ;
        if(ok)
            f->setOrder( ord ) ;
    }

}

std::string CommandLineParser::getCompleteString(std::string alias) 
{
    if(aliases.find(alias) != aliases.end())
        return aliases[alias] ;
    return alias ;
}

void CommandLineParser::disableFeatureTreeArguments() 
{
    values.erase(std::string("--set-sampling-restriction")) ;
    values.erase(std::string("--set-delta-time")) ;
    values.erase(std::string("--set-min-delta-time")) ;
    values.erase(std::string("--set-min-mesh-density")) ;
    values.erase(std::string("--set-solver-precision")) ;
    values.erase(std::string("--set-max-iterations-per-step")) ;
    values.erase(std::string("--set-sampling-number")) ;
    values.erase(std::string("--set-ssor-iterations")) ;
    values.erase(std::string("--set-surface-sampling-factor")) ;
    strings.erase(std::string("--set-order")) ;

    help.erase(std::string("--set-sampling-restriction")) ;
    help.erase(std::string("--set-delta-time")) ;
    help.erase(std::string("--set-min-delta-time")) ;
    help.erase(std::string("--set-min-mesh-density")) ;
    help.erase(std::string("--set-solver-precision")) ;
    help.erase(std::string("--set-max-iterations-per-step")) ;
    help.erase(std::string("--set-sampling-number")) ;
    help.erase(std::string("--set-ssor-iterations")) ;
    help.erase(std::string("--set-surface-sampling-factor")) ;
    help.erase(std::string("--set-order")) ;
}

CommandLineParser::CommandLineParser(std::string d, bool c, bool f) : description(d), commandLineConfiguration(c), forceUnrecognizedFlags(f)
{
    addFlag(std::string("--help"), std::string("print help"), "-h") ;
    addFlag(std::string("--version"), std::string("print current AMIE revision"), "-v") ;
    addFlag(std::string("--no-openmp"), std::string("disable OpenMP (equivalent to --set-num-threads 1)")) ;
    addFlag(std::string("--print-status"), std::string("print the list of activated command line arguments before launching the simulation")) ;
    addValue(std::string("--set-num-threads"), -1, std::string("set the number of threads available for OpenMP")) ;
    addValue(std::string("--set-sampling-restriction"), -1, std::string("set the size of the smallest inclusions meshed")) ;
    addValue(std::string("--set-delta-time"), -1, std::string("set the time step"), "-dt") ;
    addValue(std::string("--set-min-delta-time"), -1, std::string("set the minimum time increment between two damage steps")) ;
    addValue(std::string("--set-min-mesh-density"), -1, std::string("set the minimum mesh density (factor between 0 and 1)")) ;
    addValue(std::string("--set-solver-precision"), -1, std::string("set the precision of the conjugate gradient solver"), "-eps") ;
    addValue(std::string("--set-max-iterations-per-step"), -1, std::string("set the maximum number of iterations during a damage step"), "-it") ;
    addValue(std::string("--set-ssor-iterations"), -1, std::string("set the number of SSOR iterations in the solver"), "-ssor") ;
    addValue(std::string("--set-sampling-number"), -1, std::string("set the number of points on the edges of the sample")) ;
    addValue(std::string("--set-surface-sampling-factor"), -1, std::string("increases the number of points on boundaries")) ;
    addString(std::string("--set-order"), std::string(), std::string("set the order of the finite elements"), "-O") ;
    addString(std::string("--input-file"), std::string(), std::string("path to a *.ini file containing the problem description"), "-i") ;
    addString(std::string("--send-email"), std::string(), std::string("email address to send system messages"), "-m") ;
}

void CommandLineParser::setFeatureTree( FeatureTree * f, int argc, char *argv[], std::string description ) 
{
    CommandLineParser parser(description) ;
    parser.parseCommandLine( argc, argv ) ;
    parser.setFeatureTree(f) ;
}

ConfigTreeItem * CommandLineParser::parseConfigFile( std::string file, bool priority )
{
	ConfigTreeItem * tmp = ConfigParser::readFile( file, nullptr, false, false ) ;
	if(tmp->hasChild("arguments"))
	{
		int argc = 1 ;
		std::vector<std::string>  args ;
		args.push_back(file) ;
		std::string read = tmp->getStringData("arguments") ;
		size_t sep = read.find(' ') ;
		while(sep < std::string::npos)
		{
			std::string test = read.substr(0, sep) ;
			if( std::find( parsed.begin(), parsed.end(), test ) == parsed.end() || priority)
			{
				args.push_back(test) ;
				argc++ ;
			}
			read = read.substr( sep+1, std::string::npos ) ;
			sep = read.find(' ') ;
		}
		if( std::find( parsed.begin(), parsed.end(), read ) == parsed.end())
		{
			args.push_back(read) ;
			argc++ ;
		}
	
		if(args.size() > 0)
			parseCommandLine( argc, nullptr, args ) ;
	}

	input = ConfigParser::readFile( file, config, true, true, getActiveFlags(), getString("--directory")) ;
	input->configure( directConfig ) ;

	return input ;
}

ConfigTreeItem * CommandLineParser::parseCommandLine( int argc, char *argv[], std::vector<std::string> sargs )
{
	if(argv != nullptr)
	{
		command = std::string(argv[0]) ;
		commandLine = "";
		for(int i = 0 ; i < argc ; i++)
			commandLine += std::string(argv[i])+"  " ;
	}

	config = nullptr ;
	size_t i = 1 ;
	while(i < (size_t) argc)
	{
		std::string test ; 
		std::string follower ;
		if(argv != nullptr)
		{
			test = getCompleteString(std::string(argv[i])) ;
			if(i+1 < (size_t) argc)
				follower = getCompleteString(std::string(argv[i+1])) ;
		}
		else
		{
			test = getCompleteString(sargs[i]) ;
			if(i+1 < (size_t) argc)
				follower = getCompleteString(sargs[i+1]) ;
		}

		if(i-1 < arguments.size() && (test.length() < 2 || (test[0] != '-' || test[1] != '-')) )
		{
			arguments[i-1].str = test ;
			arguments[i-1].val = atof(test.c_str()) ;
		}
		else if( flags.find( test ) != flags.end() )
		{
			flags[ test ] = true ;
			parsed.push_back( test ) ;
		}
		else if( values.find( test ) != values.end()  )
		{
			values[ test ] = atof(follower.c_str()) ;
			parsed.push_back( test ) ;
			i++ ;
		}
		else if( strings.find( test ) != strings.end()  )
		{
			strings[ test ] = follower ;
			parsed.push_back( test ) ;
			i++ ;
		}
		else if(commandLineConfiguration && test[0] == '@')
		{
			if(!config)
				config = new ConfigTreeItem(nullptr, std::string("define")) ;
			std::string testval = follower ;
			bool isDouble = (testval.find_first_not_of("0123456789.e-") == std::string::npos ) ;
			if(isDouble)
				new ConfigTreeItem( config, test, atof(testval.c_str()) ) ;
			else
				new ConfigTreeItem( config, test, testval ) ;
			i++ ;
		}
		else if(commandLineConfiguration && test[0] == '.')
		{
			directConfig [test.substr(1)] = follower ;
			i++ ;
		}
		else if(forceUnrecognizedFlags && test.find(std::string("--")) == 0)
		{
			flags[test] = true ;
			parsed.push_back( test ) ;
		}
		i++ ;
	}

	if( getFlag(std::string("--version")) )
		printVersion() ;
	if( getFlag(std::string("--help")) )
		printHelp() ;
	if( getFlag(std::string("--print-status")) )
		printStatus() ;

        if( getFlag(std::string("--no-openmp")))
		setNumThreads(1) ;
	else if( getValue(std::string("--set-num-threads")) > 0)
		setNumThreads( getValue(std::string("--set-num-threads") ) );

	if(getFlag(std::string("--help")) || getFlag(std::string("--version")))
		exit(0) ;

	if(argv != nullptr && strings["--input-file"].length() > 0)
	{
		parseConfigFile( strings["--input-file"], false ) ;
	}

	if(argv != nullptr && strings["--send-email"].find("@") != std::string::npos)
	{
		sendEmail( "AMIE - Execution started", commandLine ) ;
	}

	return input ;
}

void CommandLineParser::sendEmail( std::string subject, std::string body, std::string attachment)
{
	if( strings["--send-email"].find("@") == std::string::npos)
	{
		std::cerr << "improper email address, can't send email to "+strings["--send-email"] << std::endl ;
		return ;
	}

	std::string mail = "(hostname & pwd & echo \""+body+"\") > .amie_message.tmp & mail -s \""+subject+"\" " ;
	if(attachment.length() > 0)
		mail += "-a " + attachment + " ";
	mail += strings["--send-email"]+" < .amie_message.tmp" ;

	std::system(mail.c_str()) ;
}

void CommandLineParser::setNumThreads( int n )
{
#ifdef HAVE_OMP
	omp_set_num_threads(std::max(n,1)) ;
#endif
}

void CommandLineParser::printStatus( )
{
	for(size_t i = 0 ; i < arguments.size() ; i++)
		std::cout << arguments[i].name << " =  " << arguments[i].str << std::endl ;
	for(auto f = flags.begin() ; f != flags.end() ; f++)
		std::cout << f->first << " = " << (f->second ? "TRUE" : "FALSE") << std::endl ;
	for(auto f = values.begin() ; f != values.end() ; f++)
		std::cout << f->first << " = " << f->second << std::endl ;
	for(auto f = strings.begin() ; f != strings.end() ; f++)
		std::cout << f->first << " = " << f->second << std::endl ;

}

void CommandLineParser::printFormatedHelp( std::string arg, int max, int maxalias, std::string help, std::string lead, bool printAlias)
{

	int length = arg.length() ;
	if(printAlias)
	{
		length += maxalias+1 ;
		std::string alias = getAlias(arg) ;
		std::cout << alias ;
		int alength = alias.length() ;
		for(int i = 0 ; i < maxalias-alength ; i++)
			std::cout << " " ;
	}
	std::cout << arg ;
	for(int i = 0 ; i < max-length ; i++)
		std::cout << " " ;
	std::cout << lead << help << std::endl ;
}

std::string CommandLineParser::getAlias(std::string complete)
{
	for(auto i = aliases.begin() ; i != aliases.end() ; i++)
	{
		if(aliases[i->first] == complete)
		{
			return i->first ;
		}
	}
//	std::cout << "alias for " << complete << " not found!" << std::endl ;
	return std::string() ;
}

Form * CommandLineParser::getBehaviour(std::string token, Form * b, SpaceDimensionality dim)
{
	std::string file = getString(token) ;
	for(size_t i = 0 ; i < arguments.size() ; i++)
	{
		if(arguments[i].name == token)
			file = arguments[i].str ;
	}
	if(file.size() > 4 && file.find(".ini") == file.length()-4)
		return ConfigParser::getBehaviour( file, b, dim ) ;
	return b ;
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
		std::cout << "[@variable (value) ...]  " ;
	std::cout << "[options...]" << std::endl ;

	std::cout << std::endl ;
	if(arguments.size() > 0)
		std::cout << "ARGUMENTS" << std::endl ;

	size_t m = 0 ;
	for(size_t i = 0 ; i < arguments.size() ; i++)
		m = std::max(m, arguments[i].name.length() ) ;

	for(size_t i = 0 ; i < arguments.size() ; i++)
		printFormatedHelp( arguments[i].name, m, 0, arguments[i].help, "  ", false ) ;
	if(arguments.size() > 0)
		std::cout << std::endl ;

	std::cout << "OPTIONS" << std::endl ;

	size_t a = 0 ;
	for(auto al = aliases.begin() ; al != aliases.end() ; al++)
		a = std::max(a, al->first.length()+1 ) ;

	m = 0 ;
	for(auto h = help.begin() ; h != help.end() ; h++)
		m = std::max(m, h->first.length()+a+1 ) ;


	for(auto h = help.begin() ; h != help.end() ; h++)
	{
		std::string lead = "   " ;
		if(values.find(h->first) != values.end())
			lead = "   (numeral) " ;
		if(strings.find(h->first) != strings.end())
			lead = "   (string) " ;
		printFormatedHelp( h->first, m,a, h->second, lead, true ) ;
	}

	std::cout << std::endl ;
}

void CommandLineParser::printVersion( )
{
	std::cout << "current AMIE version: svn " << std::flush ;
	std::system("svnversion ..") ;
	
	struct stat attr ;
	stat( command.c_str(), &attr ) ;
	std::cout << "last compiled on: " << ctime(&attr.st_mtime) << std::endl ;

}

bool CommandLineParser::getFlag( std::string f )
{
	if( flags.find( f ) == flags.end() )
		return false ;
	return flags[f] ;
}

std::vector<std::string> CommandLineParser::getActiveFlags( )
{
	std::vector<std::string> active ;
	for(auto f = flags.begin() ; f != flags.end() ; f++)
	{
		if( f->second )
			active.push_back( f->first ) ;
	}
	return active ;
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

FunctionParser::FunctionParser(std::vector<std::string> f, std::map<std::string, char> coordinates) : isLinking(false)
{
	renewExpression(f, coordinates) ;
}

FunctionParser::FunctionParser(std::string f, std::map<std::string, char> coordinates) : isLinking(false)
{
	renewExpression( FunctionParser::breakString(f), coordinates) ;
}

std::vector<std::string> FunctionParser::breakString(std::string f)
{
	std::string current ;
	std::vector<std::string> str ;
	size_t i = 0 ;
	while(i < f.length())
	{
		bool cut = false ;
		char test = f[i] ;
		if(test == ' ' || test == '\t')
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
	return str ;
}

bool isOpenBracket( std::string s) 
{
	return s== "(" ;
}

bool isCloseBracket(std::string s) 
{
	return s==")" ;
}

std::string FunctionParser::translateCoordinate( std::string test, std::map<std::string, char> coordinates)
{
	if(coordinates.find(test) != coordinates.end())
		return std::string(1, coordinates[test]) ;
	return test ;
}

void FunctionParser::renewExpression(std::vector<std::string> str, std::map<std::string, char> coordinates) 
{
	size_t i = 0 ;
	std::vector<std::string> embeddedExpression ;
	size_t openCount = 0 ;
	
	while(i < str.size())
	{
		std::string translated = translateCoordinate( str[i], coordinates ) ;
		if(openCount == 0 && !isOpenBracket( translated )  )
		{
			if(translated.size() > 1 && translated[0] == '-' && !FunctionParserHelper::isNumeral(translated[1]))
			{
				std::vector<std::string> tmp ;
				tmp.push_back("-1") ;
				tmp.push_back("*") ;
				tmp.push_back(translateCoordinate( str[i].substr(1), coordinates )) ;
				tokens.push_back( new FunctionParser( tmp ) ) ;
			}
			else
				tokens.push_back( new FunctionParserToken( translated ) ) ;
		}
		else
		{
			if( isOpenBracket( translated ) )
			{
				if(openCount == 0)
					embeddedExpression.clear() ;
				else
					embeddedExpression.push_back( translated ) ;
				openCount++ ;
			}
			else if(isCloseBracket( str[i] ) )
			{
				openCount-- ;
				if(openCount == 0)
					tokens.push_back( new FunctionParser( embeddedExpression) ) ;
				else
					embeddedExpression.push_back( translated ) ;
			}
			else
				embeddedExpression.push_back( translated ) ;
		}
		i++ ;
	}

	this->link() ;
}

void FunctionParser::linkLeftAndRightToken( TokenOperationType op)
{
	for(size_t i = 0 ; i < tokens.size() ; i++)
	{
		if(tokens[i]->isFinalToken())
		{
			if( dynamic_cast<FunctionParserToken *>(tokens[i])->toFunctionToken() == op )
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
				if(roots[right] == -1 && right != i)
					roots[right] = i ;
				if(roots[left] == -1 && left != i)
					roots[left] = i ;
			}
		}
	}
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
				case TOKEN_OPERATION_MIN:
				case TOKEN_OPERATION_MAX:
					roots[i+1] = i ;
					roots[i+2] = i ;
					break ;
				default:
					break ;
			}
		}
//		else
//			tokens[i]->link() ;
	}

	linkLeftAndRightToken( TOKEN_OPERATION_POWER ) ;
	linkLeftAndRightToken( TOKEN_OPERATION_TIMES ) ;
	linkLeftAndRightToken( TOKEN_OPERATION_DIVIDES ) ;
	linkLeftAndRightToken( TOKEN_OPERATION_MINUS ) ;
	linkLeftAndRightToken( TOKEN_OPERATION_PLUS ) ;

//	printRoots() ;

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
			{
				if(getLeftTokenIndex(i) == getRightTokenIndex(i))
					return getRightFunction(i) ;
				return getLeftFunction(i)+getRightFunction(i) ;
			}
			case TOKEN_OPERATION_MINUS:
			{
				if(getLeftTokenIndex(i) == getRightTokenIndex(i))
					return -getRightFunction(i) ;
				return getLeftFunction(i)-getRightFunction(i) ;
			}
			case TOKEN_OPERATION_TIMES:
			{
				if(getLeftTokenIndex(i) == getRightTokenIndex(i))
					return getRightFunction(i) ;
				return getLeftFunction(i)*getRightFunction(i) ;
			}
			case TOKEN_OPERATION_DIVIDES:
			{
				if(getLeftTokenIndex(i) == getRightTokenIndex(i))
					return getRightFunction(i)^(-1) ;
				return getLeftFunction(i)/getRightFunction(i) ;
			}
			case TOKEN_OPERATION_POWER:
			{
				if(getLeftTokenIndex(i) == getRightTokenIndex(i))
					return getRightFunction(i) ;
				return getLeftFunction(i)^getRightFunction(i) ;
			}
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
			case TOKEN_OPERATION_MIN:
				return f_min( getLeftFunction(i), getRightFunction(i) ) ;
			case TOKEN_OPERATION_MAX:
				return f_max( getLeftFunction(i), getRightFunction(i) ) ;
			default:
				return Function("0") ;
		}
	}
	return tokens[i]->getFunction() ;
		
}

Function FunctionParser::getFunction( std::string f, std::map<std::string, char> coordinates )
{
	FunctionParser expr(f, coordinates) ;
//	expr.print() ;
//	expr.printRoots() ;
	return expr.getFunction() ;
}

void FunctionParser::printRoots() const
{
	for(size_t i = 0 ; i < roots.size() ; i++)
		std::cout << roots[i] << " " ;
	std::cout << std::endl ;
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
	std::vector<double> v ;
	return FunctionParserHelper::toToken( data, 0, v ).first ;
}
