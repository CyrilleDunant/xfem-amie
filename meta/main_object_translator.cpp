#include <cstdlib>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <map>
#include <algorithm>
#include <ctime>



std::vector<std::string> breakString(std::string f)
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




struct ObjectConstructor 
{
	std::string token ;
	std::string father ;
	std::vector<std::pair<std::string, std::string> > arguments ;
	std::map<std::string, std::string> preprocessor ;
	std::map<std::string, std::string> defaultValues ;
	bool reset ;

	ObjectConstructor( std::string current ) 
	{
		std::vector<std::string> command = breakString( current ) ;
		token = command[1] ;
		father = command[2] ;
		reset = false ;
		for(size_t i = 3 ; i < command.size()-1 ; i++)
		{
			if(command[i] == "--has-reset")
				reset = true ;
			if(command[i].find('@') == 0)
			{
				size_t pre = command[i].find('<') ;
				size_t post = command[i].find('>') ;
				size_t start = command[i].find('[') ;
				size_t end = command[i].find(']') ;
				std::string type = command[i].substr(1, start-1 ) ;
				std::string arg = command[i].substr(start+1, end-start-1 ) ;
				if(pre < std::string::npos)
				{
					type = command[i].substr(1, pre-1 ) ;
					std::string cast = command[i].substr( pre+1, post-pre-1 ) ;
					preprocessor[arg] = cast ;
				}
				arguments.push_back( std::make_pair( type, arg ) ) ;

				if(command[i+1].find('@') == std::string::npos && command[i+1].find("*/") == std::string::npos)
				{
					defaultValues[ arg ] = command[i+1] ;
				}
			}
		}
	}

	static std::string getMapName( std::string natural )
	{
		if(natural.find("*") < std::string::npos) 
		{
			std::string tentative ;
			std::locale loc ;
			for(size_t i = 0 ; i < natural.length() ; i++)
			{
				if(natural[i] != '*')
					tentative.append( 1, std::tolower(natural[i], loc ) ) ;
			}
			return tentative+"s" ;
		}
		return natural+"s" ;
	}

	std::string getClassName(std::string suffix)  { return token+suffix ; }

	std::string getExistenceCondition( int i )
	{
		std::string ret("            if( ") ;
		ret.append( arguments[i].first ) ;
		ret.append("s.find(\"") ;
		ret.append( arguments[i].second ) ;
		ret.append("\") == ") ;
		ret.append( arguments[i].first ) ;
		ret.append("s.end() )") ;
		return ret ;
	}

	std::string getForceExistence( int i )
	{
		std::string ret("                ") ;
		std::string type = arguments[i].first ;
		if(type == "object")
			type = preprocessor[ arguments[i].second ]+"*" ;
		ret.append(ObjectConstructor::getMapName( type )) ;
		ret.append("[\"") ;
		ret.append(arguments[i].second) ;
		ret.append("\"] = ") ;
		if( arguments[i].first == "object")
			ret.append("new ") ;
		if( arguments[i].first == "string")
			ret.append("\"") ;
		ret.append( defaultValues[arguments[i].second] ) ;
		if( arguments[i].first == "string")
			ret.append("\"") ;
		ret.append(" ;") ;
		return ret ;
	}

	std::vector<std::string> printDefaultArguments()
	{
		std::vector<std::string> ret ;	
		for(size_t i = 0 ; i < arguments.size() ; i++)
		{
			if( defaultValues.find( arguments[i].second ) != defaultValues.end() )
			{
				ret.push_back(getExistenceCondition(i)) ;
				ret.push_back(getForceExistence(i)) ;
			}
		}
		return ret ;
	}

	std::string getDynamicCast(std::string suffix, std::string target) 
	{
		std::string ret("dynamic_cast<") ;
		ret.append(getClassName(suffix)) ;
		ret.append(" *>(") ;
		ret.append(target) ;
		ret.append(")") ;
		return ret ;
	}

	std::string getConstructorArguments()
	{
		std::string ret ;
		if(arguments.size() > 0)
		{
			for(size_t i = 0 ; i < arguments.size() ; i++)
			{
				std::string type = arguments[i].first ;
				if(type == "object")
					type = preprocessor[ arguments[i].second ]+"*" ;
				if(preprocessor.find(arguments[i].second) != preprocessor.end())
				{
					if(arguments[i].first == "string")
					{
						ret.append( "Enum::get" ) ;
						ret.append( preprocessor[arguments[i].second] ) ;
					}
					else
					{
						ret.append( "dynamic_cast<" ) ;
						ret.append( preprocessor[arguments[i].second] ) ;
						ret.append( " * >" ) ;
					}
					ret.append( "(" ) ;
				}
				ret.append( ObjectConstructor::getMapName( type ) ) ;
				ret.append("[\"") ;
				ret.append( arguments[i].second ) ;
				ret.append("\"]") ;
				if(preprocessor.find(arguments[i].second) != preprocessor.end())
					ret.append( ")" ) ;
				if( i+1 < arguments.size() )
					ret.append(", ") ;
			}
		}
		return ret ;
	}

	std::string getResetLine( std::string suffix ) 
	{
		if(!reset)
			return std::string() ;
		std::string ret("        if(");
		ret.append( getDynamicCast(suffix, "target") ) ;
		ret.append( " != nullptr) { " ) ;
		ret.append( getDynamicCast(suffix, "target") ) ;
		ret.append( "->reset(" ) ;
		ret.append( getConstructorArguments() ) ;
		ret.append( ") ; }" ) ;
		return ret ;
	}

	std::string getEvaluationTypeLine()
	{
		std::string ret("        if( type == \"");
		ret.append( token ) ;
		ret.append( "\" )" ) ;
		return ret ;
	}

	std::string getNewConstructorLine(std::string suffix)
	{
		std::string ret("            return new " ) ;
		ret.append( getClassName(suffix) ) ;
		ret.append( "(" ) ;
		ret.append( getConstructorArguments() ) ;
		ret.append( ") ;" ) ;
		return ret ;
	}

	std::string getResetEvaluationTypeLine(std::string suffix)
	{
		std::string ret("        if(");
		ret.append( getDynamicCast(suffix, "target") ) ;
		ret.append( " != nullptr)" ) ;
		return ret ;
	}

	std::string getFinalResetLine(std::string suffix)
	{
		std::string ret("            " ) ;
		ret.append( getDynamicCast(suffix, "target") ) ;
		ret.append( "->reset(" ) ;
		ret.append( getConstructorArguments() ) ;
		ret.append( ") ;" ) ;
		return ret ;
	}

	std::string getConstructorLine(std::string suffix)  
	{
		std::string ret("        if( type == \"");
		ret.append( token ) ;
		ret.append( "\" ) { return new " ) ;
		ret.append( getClassName(suffix) ) ;
		ret.append( "(" ) ;
		ret.append( getConstructorArguments() ) ;
		ret.append( ") ; }" ) ;
		return ret ;
	}

	std::vector<std::string> getConstructorLines( std::string suffix)
	{
		std::vector<std::string> ret ;
		if(defaultValues.empty())
		{
			ret.push_back( getConstructorLine(suffix) ) ;
			return ret ;
		}
		ret.push_back(getEvaluationTypeLine()) ;
		ret.push_back("        { ") ;
		std::vector<std::string> def = printDefaultArguments() ;
		for(size_t i = 0 ; i < def.size() ; i++)
			ret.push_back( def[i] ) ;
		ret.push_back(getNewConstructorLine(suffix)) ;
		ret.push_back("        }") ;
		return ret ;
	}

	std::vector<std::string> getAllTypes()
	{
		std::vector<std::string> ret ;
		for(size_t i = 0 ; i < arguments.size() ; i++)
		{
			std::string type = arguments[i].first ;
			if(type == "object")
				type = preprocessor[ arguments[i].second ]+"*" ;

			if(std::find(  ret.begin(), ret.end() , type) == ret.end())
				ret.push_back( type ) ;
		}
		return ret ;
	}

	std::vector<std::string> getResetLines( std::string suffix ) 
	{
		std::vector<std::string> ret ;
		if(defaultValues.empty())
		{
			std::string reset = getResetLine(suffix) ;
			if(reset.length() > 0)
				ret.push_back(reset) ;
			return ret ;
		}
		ret.push_back( getResetEvaluationTypeLine(suffix) ) ;
		ret.push_back("        { ") ;
		std::vector<std::string> def = printDefaultArguments() ;
		for(size_t i = 0 ; i < def.size() ; i++)
			ret.push_back( def[i] ) ;
		ret.push_back(getFinalResetLine(suffix)) ;
		ret.push_back("        }") ;
		return ret ;
	}

} ;

struct ObjectHeaderFile
{

	std::string file ;
	std::vector<ObjectConstructor> construct ;

	ObjectHeaderFile(std::string f, std::string type, bool open = true) : file(f) { if(open){getAllObjectConstructor(type) ; } }

	std::string getInclude() 
	{
		std::string ret ("#include \"") ;
		ret.append( file ) ;
		ret.append("\"") ;
		return ret ;
	}

	std::string getHeaderComment()  { return "        // parsed from header file: "+file ; }

	void getAllObjectConstructor( std::string type) 
	{
		std::fstream test ;
		test.open(file.c_str(), std::ios::in) ;

		while(!test.eof())
		{
			std::string current ;
			getline( test, current ) ;
			if(current.find("/*PARSE") == 0)
			{
				if(current.find("*/") == std::string::npos)
				{
					bool found = false ;
					while(!test.eof() && !found)
					{
						std::string next ;
						getline( test, next ) ;
						found = next.find("*/") != std::string::npos ;
						size_t comment = next.find("//") ;
						current.append(next.substr(0, comment-1)+" ") ;
					}
				}
				if(current.find( type) < std::string::npos)
					construct.push_back( ObjectConstructor(current) ) ;
			}
		}
		test.close() ;
	}

	std::vector<std::string> getAllTypes()
	{
		std::vector<std::string> ret ;
		for(size_t i = 0 ; i < construct.size() ; i++)
		{
			std::vector<std::string> t = construct[i].getAllTypes() ;
			for(size_t j = 0 ; j < t.size() ; j++)
			{
				if(std::find( ret.begin(), ret.end(), t[j] ) == ret.end())
					ret.push_back( t[j] ) ;
			}
		}
		return ret ;
	}

	std::vector<std::string> getAllResetTypes()
	{
		std::vector<std::string> ret ;
		for(size_t i = 0 ; i < construct.size() ; i++)
		{
			if(!construct[i].reset)
				continue ;
			std::vector<std::string> t = construct[i].getAllTypes() ;
			for(size_t j = 0 ; j < t.size() ; j++)
			{
				if(std::find( ret.begin(), ret.end(), t[j] ) == ret.end())
					ret.push_back( t[j] ) ;
			}
		}
		return ret ;
	}

} ;


struct AMIEConstructorParser
{
	std::string className ;
	std::string suffix ;
	std::vector<std::string> types ;
	std::vector<std::string> resettypes ;
	std::vector<ObjectHeaderFile> headers ;
	std::vector<std::string> req ;

	AMIEConstructorParser( std::string c, std::string s, std::vector<std::string> r ) : className(c), suffix( s ), req(r) { } ;

	static std::string getCppType( std::string natural )
	{
		if(natural == "string") { return "std::string" ; } 
		if(natural == "stringlist") { return "std::vector<std::string>" ; } 
		if(natural == "char") { return "char" ; } 
		if(natural == "value") { return "double" ; } 
		if(natural == "vector") { return "std::vector<double>" ; } 
		if(natural.find("*") < std::string::npos) { return natural ; }
		return "double" ;
	}

	void parseFolder( std::string path ) 
	{
		DIR * dp ;
		struct dirent *dirp ;
		if((dp = opendir(path.c_str())) == NULL)
			std::cout << "test directory not found!" << std::endl ;
	
		while((dirp = readdir(dp)) != NULL)
		{
			std::string test = dirp->d_name ;
			if(test.find(".h") == test.size()-2 && test.size() > 2)
			{
				headers.push_back( ObjectHeaderFile( path+test, className ) ) ;
				std::vector<std::string> t = headers[headers.size()-1].getAllTypes() ;
				for(size_t j = 0 ; j < t.size() ; j++)
				{
					if(std::find( types.begin(), types.end(), t[j] ) == types.end())
						types.push_back( t[j] ) ;
				}
				std::vector<std::string> rt = headers[headers.size()-1].getAllResetTypes() ;
				for(size_t j = 0 ; j < rt.size() ; j++)
				{
					if(std::find( resettypes.begin(), resettypes.end(), rt[j] ) == resettypes.end())
						resettypes.push_back( rt[j] ) ;
				}
			}
		}
	}

	std::string getConstructorDeclaration(bool header)
	{
		std::string ret = "    " ;
		if(header) { ret.append("static ") ; }
		ret.append(className) ;
		if(header) { ret.append(" * get") ; }
		else { ret.append(" * Object::get") ; }		
		ret.append(className) ;
		ret.append("(std::string type") ;
		for(size_t i = 0 ; i < types.size() ; i++)
		{
			ret.append(", ") ;
			ret.append( "std::map<std::string, " ) ;
			ret.append( AMIEConstructorParser::getCppType( types[i] ) ) ;
			ret.append( "> & " ) ;
			ret.append( ObjectConstructor::getMapName( types[i] ) ) ;
		}
		ret.append(")") ;
		if(header) { ret.append(" ;") ; }
		return ret ;
	}

	std::string getResetDeclaration(bool header)
	{
		std::string ret = "    " ;
		if(header) { ret.append("static void reset") ; }
		else { ret.append("void Object::reset") ; }
		ret.append(className) ;
		ret.append("(") ;
		ret.append(className) ;
		ret.append(" * target") ;
		for(size_t i = 0 ; i < resettypes.size() ; i++)
		{
			ret.append(", ") ;
			ret.append( "std::map<std::string, " ) ;
			ret.append( AMIEConstructorParser::getCppType( resettypes[i] ) ) ;
			ret.append( "> & " ) ;
			ret.append( ObjectConstructor::getMapName( types[i] ) ) ;
		}
		ret.append(")") ;
		if(header) { ret.append(" ;") ; }
		return ret ;
	}

	std::vector<std::string> printConstructor()
	{
		std::vector<std::string> lines ;
		lines.push_back( getConstructorDeclaration(false) ) ;
		lines.push_back( "    {") ;
		for(size_t i = 0 ; i < headers.size() ; i++)
		{
			if(headers[i].construct.size() > 0)
				lines.push_back( headers[i].getHeaderComment() ) ;
			for(size_t j = 0 ; j < headers[i].construct.size() ; j++)
			{
				std::vector<std::string> l = headers[i].construct[j].getConstructorLines( suffix ) ;
				for(size_t k = 0 ; k < l.size() ; k++)
					lines.push_back( l[k] ) ;
			}
			if(headers[i].construct.size() > 0)
				lines.push_back( std::string("   ") ) ;
		}
		lines.push_back("        return nullptr ;") ;
		lines.push_back("    }") ;
		return lines ;
	}

	std::vector<std::string> printReset()
	{
		std::vector<std::string> lines ;
		lines.push_back( getResetDeclaration(false) ) ;
		lines.push_back( "    {") ;
		for(size_t i = 0 ; i < headers.size() ; i++)
		{
			if(headers[i].construct.size() > 0)
				lines.push_back( headers[i].getHeaderComment() ) ;
			for(size_t j = 0 ; j < headers[i].construct.size() ; j++)
			{
				if(!headers[i].construct[j].reset)
					continue ;
				std::vector<std::string> l = headers[i].construct[j].getResetLines( suffix ) ;
				for(size_t k = 0 ; k < l.size() ; k++)
					lines.push_back( l[k] ) ;
			}
			if(headers[i].construct.size() > 0)
				lines.push_back( std::string("   ") ) ;
		}
		lines.push_back("    }") ;
		return lines ;
	}

	std::vector<std::string> getAllRequirements( std::vector<std::string> old ) 
	{
		std::vector<std::string> ret = old ;
		for(size_t i = 0 ; i < req.size() ; i++)
		{
			if( std::find(ret.begin(), ret.end(), req[i]) == ret.end())
				ret.push_back(req[i]) ;
		}
		return ret ;
	}

} ;

int main(int argc, char *argv[])
{
	std::vector<AMIEConstructorParser> all ;

	std::vector<std::string> materialLawReq ; materialLawReq.push_back("../physics/material_laws/material_laws.h") ;
	AMIEConstructorParser materialLaw("ExternalMaterialLaw","MaterialLaw", materialLawReq) ;
	materialLaw.parseFolder("../physics/material_laws/") ;
	all.push_back(materialLaw) ;

	std::vector<std::string> req ;
	for(size_t i = 0 ; i < all.size() ; i++)
		req = all[i].getAllRequirements( req ) ;

	time_t t = std::time(0);   // get time now
	struct std::tm * now = std::localtime( & t );

	std::fstream head ;
	head.open("../utilities/object_translator.h", std::ios::out) ;
	
	head << "/* this is an auto-generated file created on " << now->tm_mday << "/" << now->tm_mon << "/" << (now->tm_year + 1900) << " at " << now->tm_hour << ":" << now->tm_min << "  */" << std::endl ;
	head << std::endl ;
	head << "#ifndef __OBJECT_TRANSLATOR_H__" << std::endl ;
	head << "#define __OBJECT_TRANSLATOR_H__" << std::endl ;
	head << std::endl ;
	for(size_t i = 0 ; i < req.size() ; i++)
	{
		ObjectHeaderFile test( req[i], "false", false ) ;
		head << test.getInclude() << std::endl ;
	}
	head << std::endl ;
	head << "namespace Amie" << std::endl ;
	head << "{" << std::endl ;
	head << std::endl ;
	head << "struct Object" << std::endl ;
	head << "{" << std::endl ;
	head << std::endl ;
	for(size_t i = 0 ; i < all.size() ; i++)
	{
		head << "    // parsed from header file: " << all[i].req[0] << std::endl ;
		head << all[i].getConstructorDeclaration( true ) << std::endl ;
		head << all[i].getResetDeclaration( true ) << std::endl ;
	}
	head << std::endl ;
	head << "} ;" << std::endl ;
	head << std::endl ;
	head << "}" << std::endl ;
	head << std::endl ;
	head << "#endif // __OBJECT_TRANSLATOR_H__" << std::endl ;


	std::fstream src ;
	src.open("../utilities/object_translator.cpp", std::ios::out) ;

	src << "/* this is an auto-generated file created on " << now->tm_mday << "/" << now->tm_mon << "/" << (now->tm_year + 1900) << " at " << now->tm_hour << ":" << now->tm_min << "  */" << std::endl ;
	src << std::endl ;
	src << "#include \"object_translator.h\"" << std::endl ;
	src << "#include \"enumeration_translator.h\"" << std::endl ;
	std::vector<std::string> done ;
	for(size_t i = 0 ; i < all.size() ; i++)
	{
		for(size_t j = 0 ; j < all[i].headers.size() ; j++)
		{
			if( std::find( done.begin(), done.end(), all[i].headers[j].file ) == done.end() )
			{
				src << all[i].headers[j].getInclude() << std::endl ;
				done.push_back( all[i].headers[j].file ) ;
			}
		}
	}
	src << std::endl ;
	src << "namespace Amie" << std::endl ;
	src << "{" << std::endl ;
	src << std::endl ;
	for(size_t i = 0 ; i < all.size() ; i++)
	{
		std::vector<std::string> cons = all[i].printConstructor() ;
		std::vector<std::string> res = all[i].printReset() ;
		for(size_t j = 0 ; j < cons.size() ; j++)
			src << cons[j] << std::endl ;
		src << std::endl ;
		for(size_t j = 0 ; j < res.size() ; j++)
			src << res[j] << std::endl ;
		src << std::endl ;
	}
	src << "}" << std::endl ;
	src << std::endl ;


	return 0 ;
}


