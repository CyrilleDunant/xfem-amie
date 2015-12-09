#include <cstdlib>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <map>
#include <ctime>
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
		bool prevEqual = false ;
		for(size_t i = 0 ; i < list.size() ; i++)
		{
			if(list[i] == "=")
				prevEqual = true ;

			if(open && ! close && !prevEqual && list[i].find('}') == std::string::npos)
			{
				std::string test = list[i] ;
				if(test[test.length()-1] == ',')
					test = list[i].substr(0, test.length()-1) ;
				if(test.find('=') != std::string::npos)
					test = list[i].substr(0, test.find('=')-1) ;
				values.push_back( test ) ;
				prevEqual = false ;
			}
			if(list[i] != "=")
				prevEqual = false ;
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

	std::string getStringConstructorLine(size_t i)
	{
		if(i+1 > values.size())
			return std::string() ;
		std::string ret("            case ") ;
		ret.append(values[i]) ;
		ret.append(": return \"") ;
		ret.append(values[i]) ;
		ret.append("\" ;") ;		
		return ret ;
	}

	std::string getDeclarationLine()
	{
		std::string ret("    static ") ;
		ret.append(type) ;
		ret.append(" get") ;
		ret.append(type) ;
		ret.append("(std::string type, bool * ok = 0)") ;
		return ret ;
	}

	std::string getStringDeclarationLine()
	{
		std::string ret("    static std::string from") ;
		ret.append(type) ;
		ret.append("(") ;
		ret.append(type) ;
		ret.append(" value)") ;
		return ret ;
	}

	std::vector<std::string> printConstructor()
	{
		std::vector<std::string> ret ;
		ret.push_back( getDeclarationLine() ) ;
		ret.push_back("    {") ;
		ret.push_back("        if(ok) { *ok = true ; }") ;
		for(size_t i = 0 ; i < values.size() ; i++)
			ret.push_back( getConstructorLine(i)) ;
		ret.push_back("        if(ok) { *ok = false ; }") ;
		ret.push_back("        return "+values[0]+" ;") ;
		ret.push_back("    }") ;		
		return ret ;
	}

	std::vector<std::string> printStringConstructor()
	{
		std::vector<std::string> ret ;
		if(type != std::string("bool"))
		{
			ret.push_back( getStringDeclarationLine() ) ;
			ret.push_back("    {") ;
			ret.push_back("        switch(value)") ;
			ret.push_back("        {") ;
			for(size_t i = 0 ; i < values.size() ; i++)
				ret.push_back( getStringConstructorLine(i)) ;
			ret.push_back("        }") ;		
			ret.push_back("        return \""+values[0]+"\" ;") ;
			ret.push_back("    }") ;
		}
		else
		{
			ret.push_back("    static std::string frombool(bool value) { return ( value ? \"true\" : \"false\" ) ; }") ;

		}
		return ret ;
	}


} ;


struct EnumerationHeaderFile
{
	std::string file ;
	std::vector<EnumerationConstructor> enums ;

	EnumerationHeaderFile(std::string f) : file(f) { getAllEnumerations() ; } 

	std::string getHeaderComment()  { return "    // parsed from header file: "+file ; }

	std::string getInclude() 
	{
		std::string ret ("#include \"") ;
		ret.append( file ) ;
		ret.append("\"") ;
		return ret ;
	}

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


struct AMIEEnumerationParser
{
	std::vector<EnumerationHeaderFile> headers ;

	void parseFolder( std::string path ) 
	{
		DIR * dp ;
		struct dirent *dirp ;
		if((dp = opendir(path.c_str())) == NULL)
		{
			return ;
		}

		while((dirp = readdir(dp)) != NULL)
		{
			std::string test = dirp->d_name ;
			if(test.find(".h") == test.size()-2 && test.size() > 2)
				headers.push_back( EnumerationHeaderFile( path+test ) ) ;
			if(test.find(".") == std::string::npos)
				parseFolder( path+test+"/" ) ;
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
					std::vector<std::string> tmpstr = headers[i].enums[j].printStringConstructor() ;
					for(size_t k = 0 ; k < tmpstr.size() ; k++)
						ret.push_back(tmpstr[k]) ;
					ret.push_back("   ") ;
					done.push_back( headers[i].enums[j].type ) ;
				}

			}
		}
		return ret ;
	}

} ;


int main(int argc, char *argv[])
{
	AMIEEnumerationParser all ;
	all.parseFolder("../elements/") ;
	all.parseFolder("../features/") ;
	all.parseFolder("../filters/") ;
	all.parseFolder("../geometry/") ;
	all.parseFolder("../mesher/") ;
	all.parseFolder("../physics/") ;
	all.parseFolder("../polynomial/") ;
	all.parseFolder("../solvers/") ;
	all.parseFolder("../sparse/") ;
	all.parseFolder("../utilities/") ;

	time_t t = std::time(0);   // get time now
	struct std::tm * now = std::localtime( & t );
	
	std::fstream out ;
	out.open("../utilities/enumeration_translator.h", std::ios::out) ;
	out << "/* this is an auto-generated file created on " << now->tm_mday << "/" << now->tm_mon << "/" << (now->tm_year + 1900) << " at " << now->tm_hour << ":" << now->tm_min << "  */" << std::endl ;
	out << std::endl ;
	out << "#ifndef __ENUMERATION_TRANSLATOR_H__" << std::endl ;
	out << "#define __ENUMERATION_TRANSLATOR_H__" << std::endl ;
	out << std::endl ;

	for(size_t i = 0 ; i < all.headers.size() ; i++)
	{
		if(all.headers[i].enums.size() > 0)
			out << all.headers[i].getInclude() << std::endl ;
	}

	out << std::endl ;
	out << "namespace Amie" << std::endl ;
	out << "{" << std::endl ;
	out << std::endl ;
	out << "struct Enum" << std::endl ;
	out << "{" << std::endl ;
	out << std::endl ;
	std::vector<std::string> lines = all.printEnumerations() ;
	for(size_t i = 0 ; i < lines.size() ; i++)
		out << lines[i] << std::endl ;
	out << std::endl ;
	out << "} ;" << std::endl ;
	out << std::endl ;
	out << "}" << std::endl ;
	out << std::endl ;
	out << "#endif // __ENUMERATION_TRANSLATOR_H__" << std::endl ;


	return 0 ;
}


