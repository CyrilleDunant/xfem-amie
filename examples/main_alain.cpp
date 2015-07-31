#include <cstdlib>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <map>
#include <algorithm>



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

struct EnumerationConstructor
{
	std::string type ;
	std::vector<std::string> values ;

	EnumerationConstructor( std::string current ) 
	{
		std::vector<std::string> list = breakString( current ) ;
		bool open = false ;
		bool close = false ;
		for(size_t i = 0 ; i < list.size() ; i++)
		{
			if(open && ! close && list[i].find('}') == std::string::npos)
			{
				std::string test = list[i] ;
				if(test[test.length()-1] == ',')
					test = list[i].substr(0, test.length()-1) ;
				values.push_back( test ) ;
			}
			if(!open && list[i].find('{') < std::string::npos)
				open = true ;
			if(!close && list[i].find('}') < std::string::npos)
				close = true ;
			if(close)
			{
				std::string test = list[i+1] ;
				if(test[test.length()-1] == ';')
					test = list[i+1].substr(0, test.length()-1) ;
				type = test ;
				break ;
			}
		}
	}

	std::string getConstructorLine(size_t i)
	{
		if(i+1 > values.size())
			return std::string() ;
		std::string ret("        if( type == \"") ;
		ret.append(values[i]) ;
		ret.append("\") { return ") ;
		ret.append(values[i]) ;
		ret.append(" ; }") ;		
		return ret ;
	}

	std::string getDeclarationLine()
	{
		std::string ret("    static ") ;
		ret.append(type) ;
		ret.append(" get") ;
		ret.append(type) ;
		ret.append("(std::string type)") ;
		return ret ;
	}

	std::vector<std::string> printConstructor()
	{
		std::vector<std::string> ret ;
		ret.push_back( getDeclarationLine() ) ;
		ret.push_back("    {") ;
		for(size_t i = 0 ; i < values.size() ; i++)
			ret.push_back( getConstructorLine(i)) ;
		ret.push_back("        return "+values[0]+";") ;
		ret.push_back("    }") ;		
		return ret ;
	}


} ;


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
		ret.append(arguments[i].first) ;
		ret.append("s[\"") ;
		ret.append(arguments[i].second) ;
		ret.append("\"] = ") ;
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
				if(preprocessor.find(arguments[i].second) != preprocessor.end())
				{
					ret.append( "AMIEEnumerations::get" ) ;
					ret.append( preprocessor[arguments[i].second] ) ;
					ret.append( "(" ) ;
				}
				ret.append( arguments[i].first ) ;
				ret.append("s[\"") ;
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
			if(std::find(  ret.begin(), ret.end() , arguments[i].first) == ret.end())
				ret.push_back( arguments[i].first ) ;
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

struct EnumerationHeaderFile
{
	std::string file ;
	std::vector<EnumerationConstructor> enums ;

	EnumerationHeaderFile(std::string f) : file(f) { getAllEnumerations() ; } 

	std::string getHeaderComment()  { return "    // parsed from header file: "+file ; }

	void getAllEnumerations() 
	{
		std::fstream test ;
		test.open(file.c_str(), std::ios::in) ;

		while(!test.eof())
		{
			std::string current ;
			getline( test, current ) ;
			if(current.find("typedef enum") == 0)
			{
				bool found = false ;
				while(!test.eof() && !found)
				{
					std::string next ;
					getline( test, next ) ;
					found = next.find(";") != std::string::npos ;
					size_t comment = next.find("//") ;
					current.append(next.substr(0, comment-1)+" ") ;
				}
				enums.push_back( EnumerationConstructor( current ) ) ;
			}

		}
		test.close() ;
	}

} ;

struct ObjectHeaderFile
{

	std::string file ;
	std::vector<ObjectConstructor> construct ;

	ObjectHeaderFile(std::string f, std::string type) : file(f) { getAllObjectConstructor(type) ; } 

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

struct AMIEEnumerationParser
{
	std::vector<EnumerationHeaderFile> headers ;

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
				headers.push_back( EnumerationHeaderFile( path+test ) ) ;
		}
	}

	std::vector<std::string> printEnumerations()
	{
		std::vector<std::string> ret ;
		std::vector<std::string> done ;
		for(size_t i = 0 ; i < headers.size() ; i++)
		{
			for(size_t j = 0 ; j < headers[i].enums.size() ; j++)
			{
				if( std::find(done.begin(), done.end(), headers[i].enums[j].type ) == done.end() )
				{
					ret.push_back(headers[i].getHeaderComment()) ;
					std::vector<std::string> tmp = headers[i].enums[j].printConstructor() ;
					for(size_t k = 0 ; k < tmp.size() ; k++)
						ret.push_back(tmp[k]) ;
					ret.push_back("   ") ;
					done.push_back( headers[i].enums[j].type ) ;
				}

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

	AMIEConstructorParser( std::string c, std::string s = std::string() ) : className(c), suffix( s.length() == 0 ? c : s ) { } ;

	static std::string getCppType( std::string natural )
	{
		if(natural == "string") { return "std::string" ; } 
		if(natural == "stringlist") { return "std::vector<std::string>" ; } 
		if(natural == "char") { return "char" ; } 
		if(natural == "value") { return "double" ; } 
		if(natural == "vector") { return "std::vector<double>" ; } 
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

	std::string getConstructorDeclaration()
	{
		std::string ret = "    static " ;
		ret.append(className) ;
		ret.append(" * get") ;
		ret.append(className) ;
		ret.append("(std::string type") ;
		for(size_t i = 0 ; i < types.size() ; i++)
		{
			ret.append(", ") ;
			ret.append( "std::map<std::string, " ) ;
			ret.append( AMIEConstructorParser::getCppType( types[i] ) ) ;
			ret.append( "> & " ) ;
			ret.append( types[i]+"s" ) ;
		}
		ret.append(")") ;
		return ret ;
	}

	std::string getResetDeclaration()
	{
		std::string ret = "    static void reset" ;
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
			ret.append( types[i]+"s" ) ;
		}
		ret.append(")") ;
		return ret ;
	}

	std::vector<std::string> printConstructor()
	{
		std::vector<std::string> lines ;
		lines.push_back( getConstructorDeclaration() ) ;
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
		lines.push_back("   }") ;
		return lines ;
	}

	std::vector<std::string> printReset()
	{
		std::vector<std::string> lines ;
		lines.push_back( getResetDeclaration() ) ;
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
		lines.push_back("   }") ;
		return lines ;
	}

} ;

int main(int argc, char *argv[])
{
	AMIEConstructorParser materialLaw("ExternalMaterialLaw","MaterialLaw") ;
	materialLaw.parseFolder("../physics/material_laws/") ;
	
	std::vector<std::string> c = materialLaw.printConstructor() ;
	for(size_t i = 0 ; i < c.size() ; i++)
		std::cout << c[i] << std::endl ;



	return 0 ;
}


