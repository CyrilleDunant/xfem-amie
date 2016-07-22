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
	std::vector<std::string> helpers ;
	std::map<std::string, std::string> preprocessor ;
	std::map<std::string, std::string> defaultValues ;
	bool reset ;
	bool noSuffix ;

	ObjectConstructor( std::string current ) 
	{
		std::vector<std::string> command = breakString( current ) ;
		token = command[1] ;
		father = command[2] ;
		reset = false ;
		noSuffix = false ;
		for(size_t i = 3 ; i < command.size()-1 ; i++)
		{
			if(command[i] == "--has-reset")
				reset = true ;
			if(command[i] == "--no-suffix")
				noSuffix = true ;
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

				if(command[i+1].find('@') == std::string::npos && command[i+1].find("*/") == std::string::npos && command[i+1] != "//")
				{
					defaultValues[ arg ] = command[i+1] ;
				}

				if(command[i+2] == "//" || command[i+1] == "//")
				{
					size_t j = i+(command[i+2] == "//" ? 3 : 2) ;
					std::string comment = command[j] ;
					while(j+1 < command.size() && command[j+1].find('@') != 0 && command[j+1].find("*/") != 0)
					{
						j++ ;
						comment.append(" "+command[j]) ;
					}
					helpers.push_back(comment) ;
				}
				else
					helpers.push_back(" ") ;
			}
		}
	}

	std::vector<std::string> printTemplate(std::string basic)
	{
		std::vector<std::string> ret ;
		ret.push_back("."+basic+" = "+token) ;
		for(size_t i = 0 ; i < arguments.size() ; i++)
		{
			std::string tmp = ".."+arguments[i].second+" = " ;
			if(defaultValues.find(arguments[i].second) != defaultValues.end() && arguments[i].first != "object")
				tmp.append(defaultValues[arguments[i].second]+" ") ;
			tmp.append("# ("+arguments[i].first+") "+helpers[i]) ;
			ret.push_back(tmp) ;
		}
		return ret ;
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

	std::string getClassName(std::string suffix)  
	{
		if(token == ".") { return suffix ; }
		if(noSuffix) { return token ; }
		return token+suffix ; 
	}

	std::string getClassKeyword()
	{
		return "      <keyword>"+token+"</keyword>" ;
	}

	std::string getExistenceCondition( int i )
	{
		std::string ret("            if( ") ;
		std::string type = arguments[i].first ;
		if(type == "object")
			type = preprocessor[ arguments[i].second ]+"*" ;
		ret.append(ObjectConstructor::getMapName( type )) ;
		ret.append(".find(\"") ;
		ret.append( arguments[i].second ) ;
		ret.append("\") == ") ;
		ret.append(ObjectConstructor::getMapName( type )) ;
		ret.append(".end() )") ;
		return ret ;
	}

	std::string getForceExistence( int i )
	{
		std::string type = arguments[i].first ;
		if(type == "object")
			type = preprocessor[ arguments[i].second ]+"*" ;
		std::string ret = ObjectConstructor::getMapName( type ) ;
		ret.append("[\"") ;
		ret.append(arguments[i].second) ;
		ret.append("\"] = ") ;
		if( arguments[i].first == "object" && defaultValues[arguments[i].second] != "nullptr")
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
				std::string def = getExistenceCondition(i) ;
				def.append(" { ") ; 
				def.append(getForceExistence(i)) ;
				def.append(" } ; ") ; 
				ret.push_back(def) ;
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
				if(preprocessor.find(arguments[i].second) != preprocessor.end() && arguments[i].first == "string")
				{
					ret.append( "Enum::get" ) ;
					ret.append( preprocessor[arguments[i].second] ) ;
					ret.append( "(" ) ;
				}
				ret.append( ObjectConstructor::getMapName( type ) ) ;
				if(arguments[i].second.length() > 0)
                                {
					ret.append("[\"") ;
					ret.append( arguments[i].second ) ;
					ret.append("\"]") ;
				}
				if(preprocessor.find(arguments[i].second) != preprocessor.end()  && arguments[i].first == "string")
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

	std::string getEvalLine()  
        {
                std::string ret("        if( type == \"") ;
                ret.append( token ) ;
                ret.append("\" ) { return true ; }") ;
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
//						size_t comment = next.find("//") ;
						current.append(next+" ") ;
					}
				}
				std::vector<std::string> test = breakString( current ) ;
				if(test[2] == type)
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
	std::string alias ;
	std::vector<std::string> types ;
	std::vector<std::string> resettypes ;
	std::vector<ObjectHeaderFile> headers ;
	std::vector<std::string> req ;

	AMIEConstructorParser( std::string c, std::string s, std::string a, std::vector<std::string> r ) : className(c), suffix( s ), alias( a ), req(r) { } ;

	static std::string getCppType( std::string natural )
	{
		if(natural == "string") { return "std::string" ; } 
		if(natural == "stringlist") { return "std::vector<std::string>" ; } 
		if(natural == "char") { return "char" ; } 
		if(natural == "bool") { return "bool" ; } 
		if(natural == "point") { return "Point" ; } 
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

	std::string getEvalDeclaration(bool header)
	{
		std::string ret = "    " ;
		if(header) { ret.append("static bool is") ; }
		else { ret.append("bool Object::is") ; }
		ret.append(className) ;
		ret.append("(std::string type)") ;
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

	std::vector<std::string> printEval()
	{
		std::vector<std::string> lines ;
		lines.push_back( getEvalDeclaration(false) ) ;
		lines.push_back( "    {") ;
		for(size_t i = 0 ; i < headers.size() ; i++)
		{
			if(headers[i].construct.size() > 0)
				lines.push_back( headers[i].getHeaderComment() ) ;
			for(size_t j = 0 ; j < headers[i].construct.size() ; j++)
			{
				lines.push_back( headers[i].construct[j].getEvalLine() ) ;
			}
			if(headers[i].construct.size() > 0)
				lines.push_back( std::string("   ") ) ;
		}
		lines.push_back("        return false ;") ;
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

	void printAllTemplates()
	{
		std::fstream out ;
		std::string file = "../examples/data/templates/"+alias+".ini" ;
		out.open(file.c_str(), std::ios::out) ;
		for(size_t i = 0 ; i < headers.size() ; i++)
		{
			for(size_t j = 0 ; j < headers[i].construct.size() ; j++)
			{
				std::vector<std::string> tmp = headers[i].construct[j].printTemplate(alias) ;
				for(size_t k = 0 ; k < tmp.size() ; k++)
					out << tmp[k] << std::endl ;
			}
		}
		out.close() ;
	}

} ;

int main(int argc, char *argv[])
{
	std::vector<AMIEConstructorParser> all ;

	std::vector<std::string> materialLawReq ; materialLawReq.push_back("../physics/material_laws/material_laws.h") ;
	AMIEConstructorParser materialLaw("ExternalMaterialLaw","MaterialLaw","material_law", materialLawReq) ;
	materialLaw.parseFolder("../physics/material_laws/") ;

	all.push_back(materialLaw) ;

	std::vector<std::string> behaviourReq ; 
	behaviourReq.push_back("../elements/integrable_entity.h") ;
	behaviourReq.push_back("../physics/material_laws/material_laws.h") ;
	behaviourReq.push_back("../physics/material_laws/logcreep_accumulator.h") ;
	behaviourReq.push_back("../physics/damagemodels/damagemodel.h") ;
	behaviourReq.push_back("../physics/fracturecriteria/fracturecriterion.h") ;
	AMIEConstructorParser behaviour("Form",std::string(),"behaviour", behaviourReq) ;
        behaviour.parseFolder("../physics/") ;
        behaviour.parseFolder("../physics/materials/") ;

	all.push_back(behaviour) ;

	std::vector<std::string> damageReq ; 
	damageReq.push_back("../physics/damagemodels/damagemodel.h") ;
	AMIEConstructorParser damage("DamageModel","LinearDamage","damage_model", damageReq) ;
        damage.parseFolder("../physics/damagemodels/") ;

	all.push_back(damage) ;

	std::vector<std::string> criterionReq ; 
	criterionReq.push_back("../physics/fracturecriteria/fracturecriterion.h") ;
	AMIEConstructorParser criterion("FractureCriterion",std::string(),"fracture_criterion", criterionReq) ;
        criterion.parseFolder("../physics/fracturecriteria/") ;

	all.push_back(criterion) ;

	std::vector<std::string> accReq ; 
	accReq.push_back("../physics/material_laws/logcreep_accumulator.h") ;
	AMIEConstructorParser acc("LogCreepAccumulator","LogCreepAccumulator","accumulator", accReq) ;
        acc.parseFolder("../physics/material_laws/") ;

	all.push_back(acc) ;

	std::vector<std::string> incGenReq ; incGenReq.push_back("../features/microstructuregenerator.h") ;
	AMIEConstructorParser incGen("InclusionGenerator","InclusionGenerator","inclusion_generator", incGenReq) ;
	incGen.parseFolder("../features/") ;

	all.push_back(incGen) ;

	std::vector<std::string> granuloReq ; granuloReq.push_back("../utilities/granulo.h") ;
	AMIEConstructorParser granulo("ParticleSizeDistribution",std::string(),"particle_size_stribution", granuloReq) ;
	granulo.parseFolder("../utilities/") ;

	all.push_back(granulo) ;

	std::vector<std::string> familyReq ; familyReq.push_back("../utilities/inclusion_family.h") ;
	AMIEConstructorParser family("InclusionFamily","InclusionFamily","family", familyReq) ;
	family.parseFolder("../utilities/") ;

	all.push_back(family) ;

	std::vector<std::string> enrichmentReq ; enrichmentReq.push_back("../features/features.h") ;
	AMIEConstructorParser enrichment("EnrichmentManager","Manager","enrichment", enrichmentReq) ;
	enrichment.parseFolder("../features/enrichmentmanagers/") ;

	all.push_back(enrichment) ;

	std::vector<std::string> samplerReq ; samplerReq.push_back("../geometry/sampler/sampler.h") ;
	AMIEConstructorParser sampler("Sampler","Sampler","sampler", samplerReq) ;
        sampler.parseFolder("../geometry/sampler/") ;

	all.push_back(sampler) ;

	std::vector<std::string> postproReq ; postproReq.push_back("../utilities/postprocessor.h") ;
	AMIEConstructorParser postpro("PostProcessor","PostProcessor","post_processor", postproReq) ;
	postpro.parseFolder("../utilities/") ;

	all.push_back(postpro) ;

/*	std::vector<std::string> stepsReq ; stepsReq.push_back("../utilities/time_step_manager.h") ;
	AMIEConstructorParser steps("TimeStepManager","TimeStepManager","time_step_manager", stepsReq) ;
	steps.parseFolder("../utilities/") ;

	all.push_back(steps) ;*/

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
		head << all[i].getEvalDeclaration( true ) << std::endl ;
		head << all[i].getResetDeclaration( true ) << std::endl ;
		head << std::endl ;
	}
	head << "} ;" << std::endl ;
	head << std::endl ;
	head << "}" << std::endl ;
	head << std::endl ;
	head << "#endif // __OBJECT_TRANSLATOR_H__" << std::endl ;
	head.close() ;


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
			if( std::find( done.begin(), done.end(), all[i].headers[j].file ) == done.end() && all[i].headers[j].construct.size() > 0)
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
		std::vector<std::string> evl = all[i].printEval() ;
		std::vector<std::string> res = all[i].printReset() ;
		for(size_t j = 0 ; j < cons.size() ; j++)
			src << cons[j] << std::endl ;
		src << std::endl ;
		for(size_t j = 0 ; j < evl.size() ; j++)
			src << evl[j] << std::endl ;
		src << std::endl ;
		for(size_t j = 0 ; j < res.size() ; j++)
			src << res[j] << std::endl ;
		src << std::endl ;
	}
	src << "}" << std::endl ;
	src << std::endl ;

	for(size_t i = 0 ; i < all.size() ; i++)
		all[i].printAllTemplates() ;

	src.close() ;


	std::fstream xml ;
	xml.open( "class_keyword.xml", std::ios::out ) ;
	for(size_t i = 0 ; i < all.size() ; i++)
	{
		for(size_t j = 0 ; j < all[i].headers.size() ; j++)
		{
			for(size_t k = 0 ; k < all[i].headers[j].construct.size() ; k++)
			{
				if( all[i].headers[j].construct[k].token != "." )
					xml << all[i].headers[j].construct[k].getClassKeyword( ) << std::endl ;
			}
		}
	} 

	return 0 ;
}


