// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "command_line_parser.h"
#include "config_parser.h"
#include "parser.h"
#include "../enumeration_translator.h"
#include <sys/stat.h>

using namespace Amie ;

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
//#ifdef HAVE_OMP
        std::cout << "disabling OpenMP..." << std::endl ;
	omp_set_num_threads(std::max(n,1)) ;
//#endif
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

