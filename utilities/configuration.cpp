
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//

#include "configuration.h"
#include "../physics/stiffness.h"
#include "../physics/stiffness_and_fracture.h"
#include "../physics/stiffness_with_imposed_deformation.h"
#include "../physics/logarithmic_creep_with_external_parameters.h"
#include "../physics/viscoelasticity.h"
#include "../physics/viscoelasticity_and_fracture.h"
#include "../physics/viscoelasticity_and_imposed_deformation.h"
#include "../physics/material_laws/humidity_material_laws.h"
#include "../physics/material_laws/mechanical_material_laws.h"
#include "../physics/material_laws/temperature_material_laws.h"
#include "../physics/materials/aggregate_behaviour.h"
#include "../physics/materials/paste_behaviour.h"
#include "../physics/materials/c3s_behaviour.h"
#include "../physics/materials/ch_behaviour.h"
#include "../physics/materials/csh_behaviour.h"
#include "../physics/materials/gel_behaviour.h"
#include "../physics/materials/rebar_behaviour.h"
#include "../physics/materials/steel_behaviour.h"
#include "../physics/materials/concrete_behaviour.h"
#include "../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h"
#include "../physics/damagemodels/spacetimeisotropiclineardamage.h"
#include "../physics/damagemodels/isotropiclineardamage.h"
#include "../physics/damagemodels/prandtlgrauertplasticstrain.h"
#include "../physics/fracturecriteria/maxstrain.h"
#include "../physics/fracturecriteria/spacetimemultisurfacefracturecriterion.h"
#include "../physics/fracturecriteria/spacetimemultilinearsofteningfracturecriterion.h"
#include "../physics/fracturecriteria/spacetimeflowrule.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/confinedmohrcoulombwithstrain.h"
#include "../physics/fracturecriteria/confinedmohrcoulomb.h"
#include "../physics/fracturecriteria/boundedvonmises.h"
#include "../physics/fracturecriteria/limitstrains.h"
#include "../physics/fracturecriteria/mazars.h"
#include "../physics/fracturecriteria/nonlocalvonmises.h"
#include "enumeration_translator.h"
#include "object_translator.h"
#include "writer/triangle_writer.h"
#include "parser/config_parser.h"
#include "iostream"
#include "fstream"


namespace Amie {

ConfigTreeItem::ConfigTreeItem() : father(nullptr), label("trunk"), data(0.), str("")
{

}

ConfigTreeItem::ConfigTreeItem(ConfigTreeItem * f, std::string l) : father(f), label(l), data(0.), str("")
{
    if(f)
        f->addChild(this) ;
}

ConfigTreeItem::ConfigTreeItem(ConfigTreeItem * f, std::string l, double d) : father(f), label(l), data(d), str("")
{
    if(f)
        f->addChild(this) ;
}

ConfigTreeItem::ConfigTreeItem(ConfigTreeItem * f, std::string l, std::string s) : father(f), label(l), data(0.), str(s)
{
    if(f)
        f->addChild(this) ;
}

ConfigTreeItem * ConfigTreeItem::getFather() const
{
    return father ;
}

void ConfigTreeItem::configure( std::map<std::string, std::string> & cnf ) 
{
    for(auto i = cnf.begin() ; i != cnf.end() ; i++)
    {
        std::string childLabel = i->first ;
        std::string testval = i->second ;
        bool isDouble = (testval.find_first_not_of("0123456789.e-") == std::string::npos ) ;
        ConfigTreeItem * child = getChildFromFullLabel( childLabel ) ;
        if(child == nullptr)
            break ;
        if(isDouble)
            child->setData( atof( testval.c_str() ) ) ;
        else
            child->setStringData(testval) ;
    }
}


ConfigTreeItem * ConfigTreeItem::getChild(std::string childLabel)  const
{
    size_t pos_start = childLabel.find_first_of("@") ;
    if(pos_start > 0 && pos_start < std::string::npos)
    {
        std::string title = childLabel.substr( 0, pos_start) ;
        std::string remainder = childLabel.substr( pos_start+1 ) ;
        size_t index = atof(remainder.c_str()) ;
        std::vector<ConfigTreeItem *> test = getAllChildren(title) ;
        if(test.size() > index)
            return test[index] ;
        else
            return nullptr ;
    }


    for(size_t i = 0 ; i < children.size() ; i++)
    {
        if(children[children.size()-1-i]->is(childLabel))
            return children[children.size()-1-i] ;
    }
    return nullptr ;
}

ConfigTreeItem * ConfigTreeItem::getLastChild() const
{
    if(children.size() > 0)
        return children[children.size()-1] ;
    return nullptr ;
}

std::vector<ConfigTreeItem *> ConfigTreeItem::getAllChildren(std::string childLabel) const
{
    std::vector<ConfigTreeItem *> ret ;
    for(size_t i = 0 ; i < children.size() ; i++)
    {
        if(children[i]->is(childLabel))
            ret.push_back(children[i]) ;
    }
    return ret ;
}

std::vector<ConfigTreeItem *> ConfigTreeItem::getAllChildren() const
{
    return children ;
}

std::map<std::string, double> ConfigTreeItem::getDataMap( std::vector<std::string> except ) 
{
    std::map<std::string, double> ret ;
    for(size_t i = 0 ; i < children.size() ; i++)
    {
        bool good = true ;
        for(size_t j = 0 ; good && j < except.size() ; j++)
        {
             if(children[i]->is( except[j] )) { good = false ; }
        }
        if(good)
            ret[ children[i]->getLabel() ] = children[i]->getData() ;
    }
    return ret ;
}

std::map<std::string, std::string> ConfigTreeItem::getStringDataMap( std::vector<std::string> except ) 
{
    std::map<std::string, std::string> ret ;
    for(size_t i = 0 ; i < children.size() ; i++)
    {
        bool good = true ;
        for(size_t j = 0 ; good && j < except.size() ; j++)
        {
             if(children[i]->is( except[j] )) { good = false ; }
        }
        if(good)
            ret[ children[i]->getLabel() ] = children[i]->getStringData() ;
    }
    return ret ;
}

ConfigTreeItem * ConfigTreeItem::getChild(std::vector<std::string> childLabelDecomposed) const
{
    if(childLabelDecomposed.size() == 0)
        return nullptr ;
    if(childLabelDecomposed.size() == 1)
        return getChild(childLabelDecomposed[0]) ;

    ConfigTreeItem * child = getChild(childLabelDecomposed[0]) ;
    if(child)
    {
        std::vector<std::string> newDecomposition ;
        for(size_t i = 1 ; i < childLabelDecomposed.size() ; i++)
        {
            newDecomposition.push_back(childLabelDecomposed[i]) ;
        }
        return child->getChild( newDecomposition ) ;
    }
    return nullptr ;
}

ConfigTreeItem * ConfigTreeItem::getChildFromFullLabel(std::string childLabel)  const
{
    return getChild( ConfigTreeItem::decompose( childLabel ) ) ;
}

void ConfigTreeItem::removeAllChildren()
{
    for(size_t i = 0 ; i < children.size() ; i++)
        children[i]->setFather( nullptr ) ;
    children.resize(0) ;
}

void ConfigTreeItem::removeChild(std::string childLabel)
{
    if(!hasChild(childLabel))
        return ;

    std::vector<ConfigTreeItem *> next ;
    for(size_t i = 0 ; i < children.size() ; i++)
    {
        if(children[i]->is(childLabel))
            children[i]->setFather( nullptr ) ;
        else
            next.push_back( children[i] ) ;
    }

    children.resize(0) ;
    for(size_t i = 0 ; i < next.size() ; i++)
        addChild( next[i] ) ;
}

void ConfigTreeItem::removeChild(std::vector<std::string> childLabelDecomposed)
{
    if(childLabelDecomposed.size() == 1)
    {
        removeChild(childLabelDecomposed[0]) ;
        return ;
    }
    std::vector<std::string> nextChildLabelDecomposed ;
    for(size_t i = 1 ; i < childLabelDecomposed.size() ; i++)
        nextChildLabelDecomposed.push_back( childLabelDecomposed[i] ) ;
    std::vector<ConfigTreeItem *> all = getAllChildren(childLabelDecomposed[0]);
    for(size_t i = 0 ; i < all.size() ; i++)
        all[i]->removeChild( nextChildLabelDecomposed ) ;
}

void ConfigTreeItem::removeChildFromFullLabel(std::string childLabel)
{
    removeChild( ConfigTreeItem::decompose( childLabel ) ) ;
}

void ConfigTreeItem::addChild(ConfigTreeItem * c)
{
    if(c)
    {
        children.push_back(c) ;
        c->forceSetFather(this) ;
    }
}

void ConfigTreeItem::addChildren(std::vector<ConfigTreeItem *> c)
{
    for(size_t i = 0 ; i < c.size() ; i++)
        addChild(c[i]) ;
}

void ConfigTreeItem::setFather(ConfigTreeItem * f)
{
    if(father && f)
        father = f ;
}

void ConfigTreeItem::forceSetFather(ConfigTreeItem * f)
{
    father = f ;
}

ConfigTreeItem * ConfigTreeItem::getRoot() const
{
    if(isTrunk())
        return nullptr ;
    ConfigTreeItem * root = this->getFather() ;
    while(!root->isTrunk())
        root = root->getFather() ;
    return root ;
}

bool ConfigTreeItem::is(std::string l) const {
    return label == l ;
}

bool ConfigTreeItem::isTrunk() const {
    return father == nullptr ;
}

bool ConfigTreeItem::isLeaf() const {
    return children.size() == 0 ;
}

std::string ConfigTreeItem::getFullLabel() const
{
    if(isTrunk())
        return label ;

    std::string path = father->getFullLabel() ;
    path.append(".") ;
    path.append(label) ;

    return path ;
}

void ConfigTreeItem::setLabel(std::string l)
{
    if(isTrunk())
        return ;

    label = l ;
}
std::vector<std::string> ConfigTreeItem::decompose(std::string full)
{
    std::vector<std::string> decomposition ;
    std::string current = full ;
    size_t found = current.find(".") ;
    while(found != std::string::npos)
    {
        decomposition.push_back(current.substr(0, found)) ;
        current = current.substr(found+1) ;
        found = current.find(".") ;
    }
    decomposition.push_back(current) ;
    return decomposition ;
}

void ConfigTreeItem::print() const
{
    std::cout << getFullLabel() ;
    if( isLeaf() || data != 0 || str.length() > 0 )
    {
        std::cout << " = " ;
        if(str.size() > 0 && str[0] != '@')
            std::cout << str ;
        else
            std::cout << data ;
    }
    std::cout << std::endl ;
}

void ConfigTreeItem::printTree() const
{
    print() ;
    for(size_t i = 0 ; i < children.size() ; i++)
        children[i]->printTree() ;
}

bool ConfigTreeItem::hasChild( std::string l )  const
{
    size_t pos_start = l.find_first_of("@") ;
    if(pos_start > 0 && pos_start < std::string::npos)
    {
        std::string title = l.substr( 0, pos_start) ;
        std::string remainder = l.substr( pos_start+1 ) ;
        size_t index = atof(remainder.c_str()) ;
        std::vector<ConfigTreeItem *> test = getAllChildren(title) ;
        if(test.size() > index)
            return true ;
        else
            return false ;
    }

    for(size_t i = 0 ; i < children.size() ; i++)
    {
        if(children[i]->is(l))
            return true ;
    }
    return false ;
}

bool ConfigTreeItem::hasChildFromFullLabel( std::string childLabel )  const
{
    return hasChild( ConfigTreeItem::decompose(childLabel) ) ;
}

bool ConfigTreeItem::hasChild( std::vector<std::string> childLabelDecomposed )  const
{
    if(childLabelDecomposed.size() == 0)
        return false ;
    if(childLabelDecomposed.size() == 1)
        return hasChild(childLabelDecomposed[0]) ;

    if(hasChild(childLabelDecomposed[0]))
    {
        std::vector<std::string> newDecomposition ;
        for(size_t i = 1 ; i < childLabelDecomposed.size() ; i++)
        {
            newDecomposition.push_back(childLabelDecomposed[i]) ;
        }
        return getChild(childLabelDecomposed[0])->hasChild( newDecomposition ) ;

    }
    else
        return false ;

}

std::map<std::string, double> ConfigTreeItem::makeDefinition() const
{
    std::map<std::string, double> ret ;
    for(size_t i = 0 ; i < children.size() ; i++)
    {
        if(children[i]->getLabel()[0] == '@')
            ret[children[i]->getLabel()] = children[i]->getData() ;
    }
    return ret ;
}

void ConfigTreeItem::define(ConfigTreeItem * def, bool over)
{
    if(hasChild("define") && def && over)
    {
//		ConfigTreeItem * def = ret->getChild("define") ;
        ConfigTreeItem * current = getChild("define") ;
        std::vector<ConfigTreeItem *> toDefine = def->getAllChildren() ;
        for(size_t i = 0 ; i < toDefine.size() ; i++)
        {
            if(current->hasChild(toDefine[i]->getLabel()))
            {
                ConfigTreeItem * toBeDefined = current->getChild( toDefine[i]->getLabel() ) ;
                toBeDefined->setData( toDefine[i]->getData() ) ;
                toBeDefined->setStringData( toDefine[i]->getStringData() ) ;
            }
        }

    }


    if(hasChild("define"))
    {
        ConfigTreeItem * definition = getChild("define") ;
//		std::map<std::string, double> def = getChild("define")->makeDefinition() ;
        for(size_t i = 0 ; i < children.size() ; i++)
        {
            if(!children[i]->is("define"))
                children[i]->define(definition) ;
        }
    }

}

void ConfigTreeItem::define(ConfigTreeItem * definition)
{
    if( str.size() > 0 && str[0] == '@' && definition->hasChild(str))
    {
        ConfigTreeItem * local = definition->getChild(str) ;
        data = local->getData() ;
        str = local->getStringData() ;
    }
    else if( str.size() > 0 && str[0] == '@' )
    {
        std::cout << "named variable " << str << " not defined in input file!" << std::endl ;
        exit(0) ;
    }

    for(size_t i = 0 ; i < children.size() ; i++)
    {
        if(!children[i]->is("define"))
            children[i]->define(definition) ;
    }
}

void ConfigTreeItem::definePath( std::string baseDirectory)
{
    if(label == std::string("file")
            || label.find("input") == 0
            || label == std::string("file_name")
            || label == std::string("tension_file_name")
            || label == std::string("compression_name")
            || label == std::string("list_of_time_steps")
            || str.find(".dat") != std::string::npos 
            || str.find(".txt") != std::string::npos )
    {
        str = baseDirectory+"/"+str ;
    }

    for(size_t i = 0 ; i < children.size() ; i++)
        children[i]->definePath(baseDirectory) ;
}

double ConfigTreeItem::getData(std::string path, double defaultValue) const
{
    if(hasChildFromFullLabel(path))
        return getChildFromFullLabel(path)->getData() ;
    return defaultValue ;
}

std::string ConfigTreeItem::getStringData(std::string path, std::string defaultValue) const
{
    if(hasChildFromFullLabel(path))
        return getChildFromFullLabel(path)->getStringData() ;
    return defaultValue ;
}

RectangularFeature * ConfigTreeItem::getSample(std::vector<ExternalMaterialLaw *> common)
{
    double sampleWidth = getData("width",1) ;
    double sampleHeight = getData("height",1) ;
    double sampleCenterX = getData("center.x",0) ;
    double sampleCenterY = getData("center.y",0) ;
    RectangularFeature * ret = new RectangularFeature(nullptr, sampleWidth, sampleHeight, sampleCenterX, sampleCenterY) ;
    if(hasChild("behaviour"))
    {
        ret->setBehaviour( getChild("behaviour")->getBehaviour( SPACE_TWO_DIMENSIONAL, common ) ) ;
    }

    return ret ;
}

Function ConfigTreeItem::getFunction() const
{
    std::string function = getStringData() ;
    if(function.length() == 0)
        function = std::to_string( getData() ) ;
    Function f(function.c_str()) ;
    return f ;
}

Vector ConfigTreeItem::readLineAsVector(std::string line, char s, char t)
{
    std::vector<double> v ;
    std::string next = line ;
    int sep = next.find(s) ;
    size_t tab = next.find(t) ;
    if(tab != std::string::npos)
    {
        std::string left = next.substr(0, tab) ;
        size_t end = next.rfind(t) ;
        std::string right = next.substr( end+1 ) ;
        std::string middle = next.substr(tab+1, end-tab-1) ;
        double dleft = std::atof( left.c_str() ) ;
        double dright = std::atof( right.c_str() ) ;
        double incr = std::abs(std::atof( middle.c_str() )) ;
        double current = std::min(dleft, dright) ;
        v.push_back( current ) ;
        while(current < std::max( dleft, dright) )
        {
            current += incr ;
            v.push_back( current ) ;
        }
    }
    else
    {
        while(sep != (int)std::string::npos)
        {
            std::string current = next.substr(0, sep) ;
            v.push_back( atof(current.c_str() ) ) ;
            next = next.substr(sep+1) ;
            sep = next.find(s) ;
        }
        v.push_back( atof( next.c_str() ) ) ;
    }
    Vector ret(v.size()) ;
    for(size_t i = 0 ; i < v.size() ; i++)
        ret[i] = v[i] ;
    return ret ;
}

std::vector<double> ConfigTreeItem::readLineAsStdVector(std::string line, char s)
{
    std::vector<double> v ;
    std::string next = line ;
    int sep = next.find(s) ;
    while(sep != (int)std::string::npos)
    {
        std::string current = next.substr(0, sep) ;
        v.push_back( atof(current.c_str() ) ) ;
        next = next.substr(sep+1) ;
        sep = next.find(s) ;
    }
    v.push_back( atof( next.c_str() ) ) ;
    return v ;
}

bool ConfigTreeItem::bindInput( std::vector<std::string> & callers, std::string dir, std::vector<std::string> flags ) 
{
    std::vector<ConfigTreeItem *> newchildren ;
    bool mod = false ;
    for(size_t i = 0 ; i < children.size() ; i++)
    {
        if( children[i]->is( "input" ) || children[i]->getLabel().find("input[") == 0)
        {
            std::string lab = children[i]->getLabel() ;
            size_t start = lab.find("[") ;
            size_t end = lab.find("]") ;
            if(start < std::string::npos && end < std::string::npos)
            {
                std::string flag = lab.substr(start+1, end-start-1) ;
                bool invert = flag[0] == '!' ;
                if(invert) { flag = flag.substr(1) ; }
                bool found = invert ;
                for(size_t f = 0 ; f < flags.size() ; f++)
                {
                    if(invert)
                        found &= (flags[f] != flag) ;
                    else
                        found |= (flags[f] == flag) ;
                }
                if(found == invert)
                {
                    std::cout << "flag " << flag << " not activated!" << std::endl ;
                    continue ;
                }
            }

            std::string next = children[i]->getStringData() ;
            std::string nextfile = next ;
            std::vector<std::string> arg ;
            std::vector<ConfigTreeItem *> got ;
            start = next.find_first_of('[') ;
            end = next.find_last_of(']') ;
            if( start != std::string::npos)
            {
                nextfile = next.substr(0, start) ;
                size_t commas = next.find_first_of(',', start+1) ;
                while(commas < std::string::npos)
                {
                    arg.push_back( next.substr( start+1, commas-start-1 ) ) ;
                    start = commas ;
                    commas = next.find_first_of(',', start+1) ;
                }
                arg.push_back( next.substr( start+1, end-start-1) ) ;
            }
            if( std::find( callers.begin(), callers.end(), nextfile ) != callers.end() )
            {
                std::cout << "warning: recursive input of file " << nextfile << " ignored" << std::endl ;
                continue ;
            }

            std::vector<std::string> nextcallers = callers ;
            nextcallers.push_back( nextfile ) ;
            ConfigTreeItem * nextitems = ConfigParser::readFile( nextfile , nullptr, false ) ;
            if(dir.size() > 0)
                nextitems->definePath( dir ) ;
            nextitems->bindInput( nextcallers, dir, flags ) ;
            if(arg.size() > 0)
            {
                for(size_t j = 0 ; j < arg.size() ; j++)
                {
                    got.push_back( nextitems->getChild( arg[j] ) )  ;
                }
            }
            else
                got = nextitems->getAllChildren() ;

            for(size_t j = 0 ; j < got.size() ; j++)
                newchildren.push_back( got[j] ) ;

            mod = true ;


        }
        else
        {
            children[i]->bindInput( callers, dir, flags ) ;
            newchildren.push_back(children[i]) ;
        }
    }

    if(mod)
    {
        children.resize(0) ;
        addChildren( newchildren ) ;
    }

    return mod ;

}


ExternalMaterialLaw * ConfigTreeItem::getExternalMaterialLaw()
{
    if(str.length() == 0 && !FunctionParserHelper::isOperator(label[label.size()-1]) )
        return nullptr ; 

    std::map<std::string, std::string> strings ;
    std::map<std::string, std::vector<std::string>> stringlists ;
    std::map<std::string, double> values ;
    strings["output"] = label ;
    strings["operation"] = "SET" ;
    if(label[label.size()-1] == '-')
         strings["operation"] = "SUBSTRACT" ;
    if(label[label.size()-1] == '+')
         strings["operation"] = "ADD" ;
    if(label[label.size()-1] == '*')
         strings["operation"] = "MULTIPLY" ;
    if(label[label.size()-1] == '/')
         strings["operation"] = "DIVIDE" ;
    if(strings["operation"] != "SET") 
    {
         strings["output"] = label.substr(0, label.size()-1) ;
         if(str.length() == 0)
         {
             std::stringstream stream ;
             stream << std::fixed << std::setprecision(16) << data ;
             str = stream.str() ;
         }
    }

    std::string type = "Eval" ;
    if(str.find("(") < std::string::npos)
    {
        std::string test = str.substr(0, str.find("(")) ;
        if( Object::isExternalMaterialLaw( test ) )
        {
            type = test ;
            std::string args = str.substr( str.find("(")+1, str.find(")")-str.find("(")-1 ) ;
            if(args.find(",") < std::string::npos)
            {
                std::vector<std::string> params ;
                size_t start = 0 ;
                size_t next = args.find(",",start+1) ;
                while( next < std::string::npos )
                {
                    std::string tmp = args.substr( start+(start > 0), next-start-(start > 0) ) ;
                    if(tmp.find("=") == std::string::npos)
                        params.push_back( tmp ) ;
                    else
                    {
                        std::string arg = tmp.substr(0, tmp.find("=") ) ;
                        double val = atof( tmp.substr(tmp.find("=")+1 ).c_str() ) ;
                        values[arg] = val ;
                        strings[arg] = tmp.substr(tmp.find("=")+1 ) ;
                    }
                    start = next ;
                    next = args.find(",",start+1) ;
                }
                std::string tmp = args.substr( start+1 ) ;
                if(tmp.find("=") == std::string::npos)
                    params.push_back( tmp ) ;
                else
                {
                    std::string arg = tmp.substr(0, tmp.find("=") ) ;
                    double val = atof( tmp.substr(tmp.find("=")+1 ).c_str() ) ;
                    values[arg] = val ;
                    strings[arg] = tmp.substr(tmp.find("=")+1 ) ;
                }
                stringlists["parameters"] = params ;
                if(params.size() > 0)
                    strings["parameter"] = params[0] ;
            }
            else
            {
                strings["parameter"] = args ;
                std::vector<std::string> params ;
                params.push_back(args) ;
                stringlists["parameters"] = params ;
            }
        }
        else if( str.find(".txt") < std::string::npos || str.find(".dat") < std::string::npos ) 
        {
            type = "LinearInterpolated" ;
            strings["file_name"] = str.substr(0, str.find("(")) ;
            strings["input"] = str.substr(str.find("(")+1, str.find(")")-str.find("(")-1 ) ;
        }
    }
    else if( Object::isExternalMaterialLaw(str))
        type = str ;
    if(type == std::string("Eval"))
    {
        strings["function"] = str ;
    }

    for(size_t i = 0 ; i < children.size() ; i++)
    {
        if( children[i]->getStringData().length() == 0)
            values[ children[i]->getLabel() ] = children[i]->getData() ;
        else
            strings[ children[i]->getLabel() ] = children[i]->getStringData() ;
    }


    return Object::getExternalMaterialLaw( type, strings, values, stringlists ) ;
}

Form * ConfigTreeItem::getBehaviour(SpaceDimensionality dim, std::vector<ExternalMaterialLaw *> common)
{
    std::map<std::string, double> values ;
    std::map<std::string, std::string> strings ;
    std::map<std::string, FractureCriterion *> frac ;
    std::map<std::string, DamageModel *> dam ;
    std::map<std::string, ExternalMaterialLawList *> law ;
    std::map<std::string, LogCreepAccumulator *> acc ;
    ExternalMaterialLawList laws ;
    strings[ "dimension" ] = Enum::fromSpaceDimensionality( dim ) ;

    for(size_t i = 0 ; i < children.size() ; i++)
    {
       if(children[i]->getStringData().length() == 0 && !FunctionParserHelper::isOperator(children[i]->getLabel()[children[i]->getLabel().length()-1])  )
          values[ children[i]->getLabel() ] = children[i]->getData() ;
       else
       {
          strings[ children[i]->getLabel() ] = children[i]->getStringData() ;

          if(children[i]->getLabel() != "fracture_criterion" || children[i]->getLabel() != "damage_model" || children[i]->getLabel() != "accumulator" )
          {
              ExternalMaterialLaw * test = children[i]->getExternalMaterialLaw() ;
              if(test != nullptr)
		{
                  laws.push_back(test) ;
		}
          }
       }
    }
    if(laws.size() > 0)
    {
        law[ "relations" ] = &laws ;
    }


    std::string type = getStringData() ;
    if(type.length() == 0)
        type = "LogarithmicCreepWithExternalParameters" ;
    if(!Object::isForm(type))
        return new VoidForm() ;

    if(hasChild("fracture_criterion"))
    {
        std::string ctype = getStringData("fracture_criterion") ;
        if(Object::isFractureCriterion(ctype))
        {
            frac[ "fracture_criterion" ] = Object::getFractureCriterion( ctype, strings, values ) ;
            if( values.find( "material_characteristic_radius" ) != values.end() )
                frac[ "fracture_criterion" ]->setMaterialCharacteristicRadius( values["material_characteristic_radius"] ) ;
        }
        if(ctype == std::string("SpaceTimeMultiSurfaceFractureCriterion"))
        {
            std::vector<ConfigTreeItem *> crits = getChild("fracture_criterion")->getAllChildren("fracture_criterion") ;
            for(size_t i = 0 ; i < crits.size() ; i++)
            {
                if( crits[i]->is("fracture_criterion") )
                {
                    std::string cctype = crits[i]->getStringData() ;
                    if(Object::isFractureCriterion(cctype))
                    {
                         std::map<std::string, double> vals = crits[i]->getDataMap( ) ;
                         std::map<std::string, std::string> strs = crits[i]->getStringDataMap(  ) ;

                         FractureCriterion * test = Object::getFractureCriterion( cctype, strs, vals ) ;
                         if(test != nullptr)
                             dynamic_cast<SpaceTimeMultiSurfaceFractureCriterion *>(frac["fracture_criterion"])->add( test ) ;
                    }

                }
            }
        }
    }
    if(hasChild("damage_model"))
    {
        std::string dtype = getStringData("damage_model") ;
        if(Object::isDamageModel(dtype))
             dam[ "damage_model" ] = Object::getDamageModel( dtype, values, strings ) ;
            if( values.find( "residual_stiffness_fraction" ) != values.end() )
                dam[ "damage_model" ]->setResidualStiffnessFraction( values["residual_stiffness_fraction"] ) ;
    }
    if(hasChild("accumulator"))
    {
        std::string atype = getStringData("accumulator") ;
        if(Object::isLogCreepAccumulator(atype))
             acc[ "accumulator" ] = Object::getLogCreepAccumulator( atype, values ) ;
    }

    return Object::getForm( type, values, frac, dam, strings, acc, law ) ;
}

Vector ConfigTreeItem::readVectorFromFile() const
{
    std::fstream in ;
    in.open( getStringData().c_str(), std::ios::in ) ;
    if(in.fail())
    {
        std::cout << "unable to read vector from file " << getStringData() << std::endl ;
        return Vector() ;
    }
    std::vector<double> tmp ;
    while(!in.eof())
    {
        double buff ;
        in >> buff ;
        tmp.push_back(buff) ;
    }
    Vector ret(tmp.size()-1) ;
    for(size_t i = 0 ; i < ret.size() ; i++)
        ret[i] = tmp[i] ;
    return ret ;
}


ParticleSizeDistribution * ConfigTreeItem::getParticleSizeDistribution() const
{
    std::string type = str ;
    if(!Object::isParticleSizeDistribution(type))
        type = "ConstantSizeDistribution" ;

    std::map<std::string, std::string> strings ;
    std::map<std::string, double> values ;
    for(size_t i = 0 ; i < children.size() ; i++)
    {
        if(children[i]->getStringData().length() == 0)
            values[ children[i]->getLabel() ] = children[i]->getData() ;
        else
            strings[ children[i]->getLabel() ] = children[i]->getStringData() ;
    }

    return Object::getParticleSizeDistribution(type, values, strings) ;


}

InclusionGenerator * ConfigTreeItem::getInclusionGenerator() const
{
    std::string type = str ;
    if(!Object::isInclusionGenerator(type))
        type = "Circular" ;
    std::map<std::string, double> parameters ;
    std::map<std::string, std::string> strings ;
    for(size_t i = 0 ; i < children.size() ; i++)
    {
        if(children[i]->getStringData().length() == 0)
            parameters[ children[i]->getLabel() ] = children[i]->getData() ;
        else
            strings[ children[i]->getLabel() ] = children[i]->getStringData() ;
    }
    return Object::getInclusionGenerator( str, parameters ) ;
}

VoronoiGrain ConfigTreeItem::getVoronoiGrain(SpaceDimensionality dim, bool spaceTime, std::vector<ExternalMaterialLaw *> common)
{
    Form * b = new VoidForm() ;
    if(hasChild("behaviour"))
        b = getChild("behaviour")->getBehaviour( dim, common ) ;
    double r = getData("radius", 0.001) ;
    double f = getData("fraction", 0.1) ;
    double c = getData("correction_factor", 1.) ;

    return VoronoiGrain( b, r, f, c) ;

}

InclusionFamily * ConfigTreeItem::getInclusionFamily(std::string type)
{
    if(type.length() == 0 || type != ".")
        type = str ;
    if(!Object::isInclusionFamily( type ))
        type = "." ;
 
    std::map<std::string, double> values ;
    std::map<std::string, std::string> strings ;
    for(size_t i = 0 ; i < children.size() ; i++)
    {
        if(children[i]->is("particle_size_distribution") || children[i]->is("geometry") || children[i]->is("behaviour") || children[i]->is("manager"))
        {
            std::map<std::string, double> cval = children[i]->getDataMap() ;
            std::map<std::string, std::string> cstr = children[i]->getStringDataMap() ;
            values.insert( cval.begin(), cval.end() ) ;
            strings.insert( cstr.begin(), cstr.end() ) ;
        }
        else if(children[i]->getStringData().length() == 0)
            values[ children[i]->getLabel() ] = children[i]->getData() ;
        else
            strings[ children[i]->getLabel() ] = children[i]->getStringData() ;
    }
    std::map<std::string, InclusionGenerator *> geometry ;
    std::map<std::string, ParticleSizeDistribution *> distribution ;
    std::string gtype = "Circular" ;
    if(hasChild("geometry"))
    {
        std::string test = getChild("geometry")->getStringData() ;
        if(Object::isInclusionGenerator(test))
            gtype = test ;
    }
    geometry["geometry"] = Object::getInclusionGenerator( gtype, values ) ;
    std::string ptype = "Constant" ;
    if(hasChild("particle_size_distribution"))
    {
        std::string test = getChild("particle_size_distribution")->getStringData() ;
        if(Object::isParticleSizeDistribution(test))
            ptype = test ;
    }
    distribution["particle_size_distribution"] = Object::getParticleSizeDistribution( ptype, values, strings ) ;

    return Object::getInclusionFamily( type, values, distribution, geometry, strings ) ;
}

Point ConfigTreeItem::getPoint(Point def ) 
{
    return Point( getData("x", def.getX()), getData("y", def.getY()), getData("z", def.getZ()), getData("t", def.getT()) ) ;
}

InclusionFamilyTree * ConfigTreeItem::makeInclusionFamilyTree( FeatureTree * F, InclusionFamily * father )
{
    SpaceDimensionality dim = F->is2D() ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL ;

    InclusionFamilyTree * tree = new InclusionFamilyTree( father ) ;

    if(!hasChild("family"))
    {
        Rectangle * placement = dynamic_cast<Rectangle *>(F->getFeature(0)) ;
        if( hasChild("placement") )
        {
            Point c = placement->getCenter() ;
            if(hasChildFromFullLabel("placement.center"))
                c = getChildFromFullLabel("placement.center")->getPoint(c) ;
            double w = getData("placement.width", placement->width()) ;
            double h = getData("placement.height", placement->height()) ;
            placement = new Rectangle( w, h, c.getX(), c.getY() ) ;
        }
        if(hasChild("surface_fraction"))
            addChild( new ConfigTreeItem( nullptr, "surface", getData("surface_fraction")*placement->area() ) ) ; 

        InclusionFamily * inc = getInclusionFamily( str ) ;
        if(hasChild("behaviour"))
        {
            std::vector<ConfigTreeItem *> forms = getAllChildren("behaviour") ;
            if(forms.size() == 1)
                inc->setBehaviour( getChild("behaviour")->getBehaviour( dim ), getStringData("copy_grain_behaviour","FALSE") == std::string("TRUE" )) ;
            else
            {
                std::vector<Form *> behaviours ;
                for(size_t i = 0 ; i < forms.size() ; i++)
                    behaviours.push_back( forms[i]->getBehaviour( dim ) ) ;
                inc->setBehaviour( behaviours, getStringData("copy_grain_behaviour","FALSE") == std::string("TRUE" )) ;
            }
        }
	tree->insert( inc ) ;
        tree->place( placement, getData("placement.spacing", 0), getData("placement.tries", 1000), getData("placement.random_seed", 1) ) ;
        if( hasChild("sampling_factor" ) )
            inc->setSamplingFactor( getData("sampling_factor" ,1.) ) ;
        inc->addToFeatureTree(F) ;
        if(hasChild("enrichment") && Object::isEnrichmentManager( getStringData("enrichment", "NoEnrichment" ) ) )
        {
            std::map<std::string, double> values ;
            std::map<std::string, std::string> strings ;
            std::map<std::string, FeatureTree *> trees ;
            std::map<std::string, InclusionFamily *> families ;
            for(size_t i = 0 ; i < children.size() ; i++)
            {
                if(children[i]->is("particle_size_distribution") || children[i]->is("geometry") || children[i]->is("behaviour") || children[i]->is("manager"))
                    continue ;
                if(children[i]->getStringData().length() == 0)
                    values[ children[i]->getLabel() ] = children[i]->getData() ;
                else
                    strings[ children[i]->getLabel() ] = children[i]->getStringData() ;
            }
            trees["feature_tree"] = F ;
            families["zones"] = inc ;
            values["index"] = 0 ;
            EnrichmentManager * manager = Object::getEnrichmentManager( getStringData("enrichment", "NoEnrichment" ), trees, families, values, strings ) ;
            if(manager != nullptr)
               F->addManager( manager ) ;
        }
        if(hasChild("inclusions"))
        {
            InclusionFamilyTree * son = getChild("inclusions")->makeInclusionFamilyTree( F, inc ) ;
            inc->sons = son ;
        }
        return tree ;
    }

    Rectangle * placement = dynamic_cast<Rectangle *>(F->getFeature(0)) ;
    if( hasChild("placement") )
    {
        Point c = placement->getCenter() ;
        if(hasChildFromFullLabel("placement.center"))
            c = getChildFromFullLabel("placement.center")->getPoint(c) ;
        placement = new Rectangle( getData("placement.width", placement->width()), getData("placement.height", placement->height()), c.getX(), c.getY() ) ;
    }

    std::vector<ConfigTreeItem *> all = getAllChildren("family") ;
    for(size_t i = 0 ; i < all.size() ; i++)
    {
        if(all[i]->hasChild("surface_fraction"))
            all[i]->addChild( new ConfigTreeItem( nullptr, "surface", all[i]->getData("surface_fraction")*placement->area() ) ) ; 
        InclusionFamily * current = all[i]->getInclusionFamily( str ) ;
        if(all[i]->hasChild("behaviour"))
            current->setBehaviour( all[i]->getChild("behaviour")->getBehaviour( dim ), getStringData("copy_grain_behaviour","FALSE") == std::string("TRUE" )) ;
	tree->insert( current ) ;
        Sampler * sampler = nullptr ;
        if(all[i]->hasChild("sampler"))
            sampler = all[i]->getChild("sampler")->getSampler() ;
        current->setSamplingFactor( getData("sampling_factor",-1), sampler ) ;
    }

    tree->place( placement, getData("placement.spacing", 0), getData("placement.tries", 1000), getData("placement.random_seed", 1) ) ;

    for(size_t i = 0 ; i < all.size() ; i++)
    {
        tree->brothers[i]->addToFeatureTree(F) ;
        if(all[i]->hasChild("enrichment") && Object::isEnrichmentManager( all[i]->getStringData("enrichment", "NoEnrichment" ) ) )
        {
            std::map<std::string, double> values ;
            std::map<std::string, std::string> strings ;
            std::map<std::string, FeatureTree *> trees ;
            std::map<std::string, InclusionFamily *> families ;
            std::vector<ConfigTreeItem *> current = all[i]->getAllChildren() ;
            for(size_t j = 0 ; j < current.size() ; j++)
            {
                if(current[j]->is("particle_size_distribution") || current[j]->is("geometry") || current[j]->is("behaviour") || current[j]->is("manager"))
                    continue ;
                if(current[j]->getStringData().length() == 0)
                    values[ current[j]->getLabel() ] = current[j]->getData() ;
                else
                    strings[ current[j]->getLabel() ] = current[j]->getStringData() ;
            }
            trees["feature_tree"] = F ;
            families["zones"] = tree->brothers[i] ;
            values["index"] = i ;
            EnrichmentManager * manager = Object::getEnrichmentManager( getStringData("enrichment", "NoEnrichment" ), trees, families, values, strings ) ;
            if(manager != nullptr)
               F->addManager( manager ) ;
        }

        if(all[i]->hasChild("inclusions"))
        {
            InclusionFamilyTree * son = all[i]->getChild("inclusions")->makeInclusionFamilyTree( F, tree->brothers[i] ) ;
            tree->brothers[i]->sons = son ;
        }
    }

    return tree ;
}

Sampler * ConfigTreeItem::getSampler() const 
{
    std::map<std::string, Point> points ;
    std::map<std::string, double> values ;
    std::map<std::string, std::string> strings ;

    std::string type = "." ;
    if( str.length() > 0) { type = str ; }

    for(size_t i = 0 ; i < children.size() ; i++)
    {
        if(children[i]->getLabel().find("point") != std::string::npos)
            points[ children[i]->getLabel() ] = children[i]->getPoint() ;
        else if(children[i]->getStringData().length() > 0)
            strings[ children[i]->getLabel() ] = children[i]->getStringData() ;
        else
            values[ children[i]->getLabel() ] = children[i]->getData() ;
    }

    if( Object::isSampler(type) )
       return Object::getSampler( type, strings, points, values ) ;

    return nullptr ;
}

std::vector<BoundaryCondition *> ConfigTreeItem::getAllBoundaryConditions(FeatureTree * F) const 
{
    std::vector<BoundaryCondition *> cond ;
    std::vector<ConfigTreeItem *> bc = getAllChildren("boundary_condition") ;
    for(size_t j = 0 ; j < bc.size() ; j++)
    {
        cond.push_back( bc[j]->getBoundaryCondition(F) ) ;
        F->addBoundaryCondition( cond[j] ) ;
    }
    return cond ;
}

BoundaryCondition * ConfigTreeItem::getBoundaryCondition(FeatureTree * f) const
{
    BoundaryCondition * ret = nullptr ;
    if( !hasChild("rate") )
    {
        if(hasChild("geometry") && hasChild("normal"))
        {
            int index = getData("geometry.index",0) ;
            Point n( getData("normal.x", 1 ), getData("normal.y", 0 )) ;
            ret = new GeometryAndFaceDefinedSurfaceBoundaryCondition( Enum::getLagrangeMultiplierType( getStringData( "condition", "GENERAL" ) ),
                dynamic_cast<Geometry *>(f->getFeature(index)), n, getData( "value", 0 ), (int) getData( "axis", 0 ) ) ;
        }
        else if(hasChild("restriction"))
        {
            double maxx = getData("restriction.top_right.x", 1.) ;
            double maxy = getData("restriction.top_right.y", 1.) ;
            double minx = getData("restriction.bottom_left.x", -1.) ;
            double miny = getData("restriction.bottom_left.y", -1.) ;
            ret = new BoundingBoxAndRestrictionDefinedBoundaryCondition( Enum::getLagrangeMultiplierType( getStringData( "condition", "GENERAL" ) ),
                    Enum::getBoundingBoxPosition( getStringData( "position", "NOW" ) ), minx, maxx, miny, maxy,
                    getData( "value", 0 ), (int) getData( "axis", 0 ) ) ;
        }
        else if(hasChild("point"))
        {
            double x = getData("point.x", 0.) ;
            double y = getData("point.y", 0.) ;
            double z = getData("point.z", 0.) ;
            double t = getData("point.t", 0.) ;
            Point p(x,y,z,t) ;
            ret = new BoundingBoxNearestNodeDefinedBoundaryCondition( Enum::getLagrangeMultiplierType( getStringData( "condition", "GENERAL" ) ),
                    Enum::getBoundingBoxPosition( getStringData( "position", "NOW" ) ), p,
                    getData( "value", 0 ), getData( "axis", 0 ) ) ;
        }
	else
	{
            ret = new BoundingBoxDefinedBoundaryCondition( Enum::getLagrangeMultiplierType( getStringData( "condition", "GENERAL" ) ),
                Enum::getBoundingBoxPosition( getStringData( "position", "NOW" ) ),
                getData( "value", 0 ), getData( "axis", 0 ) ) ;
	}

        if( hasChild("interpolation") )
            ret->setInterpolation( getStringData("interpolation") ) ;

        if( hasChild("function") )
            ret->setData( getChild("function")->getFunction() ) ;

        return ret ;

    }
    else
    {
        Function t("t") ;
        t *= getData( "rate", 0 ) ;
    
        if(hasChild("geometry") && hasChild("normal"))
        {
            int index = getData("geometry.index",0) ;
            Point n( getData("normal.x", 1 ), getData("normal.y", 0 )) ;
            ret = new GeometryAndFaceDefinedSurfaceBoundaryCondition( Enum::getLagrangeMultiplierType( getStringData( "condition", "GENERAL" ) ),
                dynamic_cast<Geometry *>(f->getFeature(index)), n, t, (int) getData( "axis", 0 ) ) ;
        }    
        else if(hasChild("restriction"))
        {
            double maxx = getData("restriction.top_right.x", 1.) ;
            double maxy = getData("restriction.top_right.y", 1.) ;
            double minx = getData("restriction.bottom_left.x", -1.) ;
            double miny = getData("restriction.bottom_left.y", -1.) ;
            ret = new BoundingBoxAndRestrictionDefinedBoundaryCondition( Enum::getLagrangeMultiplierType( getStringData( "condition", "GENERAL" ) ),
                    Enum::getBoundingBoxPosition( getStringData( "position", "NOW" ) ), minx, maxx, miny, maxy,
                    t, getData( "axis", 0 ) ) ;
        }
        else if(hasChild("point"))
        {
            double x = getData("point.x", 0.) ;
            double y = getData("point.y", 0.) ;
            double z = getData("point.z", 0.) ;
            double t = getData("point.t", 0.) ;
            Point p(x,y,z,t) ;
            ret = new BoundingBoxNearestNodeDefinedBoundaryCondition( Enum::getLagrangeMultiplierType( getStringData( "condition", "GENERAL" ) ),
                    Enum::getBoundingBoxPosition( getStringData( "position", "NOW" ) ), p,
                    t, getData( "axis", 0 ) ) ;
        }
	else
	{
            ret = new BoundingBoxDefinedBoundaryCondition( Enum::getLagrangeMultiplierType( getStringData( "condition", "GENERAL" ) ),
                Enum::getBoundingBoxPosition( getStringData( "position", "NOW" ) ),
                t, getData( "axis", 0 ) ) ;
	}
     }
     return ret ;
}

BranchedCrack * ConfigTreeItem::getCrack() const
{
    Point * a = new Point(0,0,0,0);
    Point * b = new Point(0,0,0,0);
    if(hasChild("start"))
    {
        a->setX(getData("start.x"));
        a->setY(getData("start.y"));
        a->setZ(getData("start.z"));
    }
    if(hasChild("end"))
    {
        b->setX(getData("end.x"));
        b->setY(getData("end.y"));
        b->setZ(getData("end.z"));
    }

    BranchedCrack * crack = new BranchedCrack(a, b);

    if(hasChild("propagation_method"))
    {
        std::string propagation = getStringData("propagation_method");
        if(propagation == std::string("SCORE"))
            crack->setScorePropagationMethod();
        else if(propagation == std::string("ENERGY"))
            crack->setEnergyPropagationMethod();
    }

    return crack;
}

std::vector<BranchedCrack *> ConfigTreeItem::getAllCracks() const
{
    std::vector<BranchedCrack *> cracks;
    std::vector<ConfigTreeItem *> items = getAllChildren("crack") ;
    for(size_t i = 0; i < items.size(); i++)
        cracks.push_back(items[i]->getCrack());

    return cracks;
}


bool ConfigTreeItem::isAtTimeStep(int i, int nsteps) const
{
    if(getStringData("at","ALL") == std::string("ALL"))
        return true ;
    if(getStringData("at","ALL") == std::string("LAST"))
        return i == nsteps-1 ;
    if(getStringData("at","ALL") == std::string("FIRST"))
        return i == 1 ;
    if(getStringData("at","ALL") == std::string("REGULAR"))
        return i%((int) getData("every", 2)) == 0 ;
    return false ;
}

PostProcessor * ConfigTreeItem::getPostProcessor(double instant, std::vector<unsigned int> index) const 
{
    std::string type = label ;
    std::map<std::string, std::string> strings ;
    std::map<std::string, double> values ; values["instant"] = instant ;
    if(str.length() > 0) { strings["field"]=str ; }
    for(size_t i = 0 ; i < children.size() ; i++)
    {
        if(children[i]->getStringData().length() > 0)
            strings[ children[i]->getLabel() ] = children[i]->getStringData() ;
        else
        {
            if( children[i]->is("index") && index.size() > 0)
                values[ children[i]->getLabel() ] = index[children[i]->getData()] ;
            else
                values[ children[i]->getLabel() ] = children[i]->getData() ;
        }
    }

    if( !Object::isPostProcessor( type ) )
        type = "DoNothing" ;

    return Object::getPostProcessor( type, strings, values ) ;
}

std::vector<PostProcessor *> ConfigTreeItem::getAllPostProcessors(std::vector<unsigned int> index) const 
{
    std::vector<PostProcessor *> ret ;
    double instant = (getStringData("instant","AFTER")=="AFTER")-(getStringData("instant","AFTER")=="BEFORE") ;
    for(size_t i = 0 ; i < children.size() ; i++)
    {
        if( Object::isPostProcessor( children[i]->getLabel() ) )
            ret.push_back( children[i]->getPostProcessor( instant, index ) ) ;
    }
    return ret ;
}



ConfigTreeItem * ConfigTreeItem::makeTemplate()
{
    std::string source = getStringData("source","file_not_found") ;
    if(source == std::string("file_not_found"))
    {
        std::cout << "template without source file! exiting now" << std::endl ;
        exit(0) ;
    }
    ConfigTreeItem * ret = ConfigParser::readFile(source, nullptr, false) ;

    if(hasChild("define") && ret->hasChild("define"))
    {
        ConfigTreeItem * def = ret->getChild("define") ;
        ConfigTreeItem * current = getChild("define") ;
        std::vector<ConfigTreeItem *> toDefine = current->getAllChildren() ;
        for(size_t i = 0 ; i < toDefine.size() ; i++)
        {
            if(def->hasChild(toDefine[i]->getLabel()))
            {
                ConfigTreeItem * toBeDefined = def->getChild( toDefine[i]->getLabel() ) ;
                toBeDefined->setData( toDefine[i]->getData() ) ;
                toBeDefined->setStringData( toDefine[i]->getStringData() ) ;
            }
        }

    }

    return ret ;

}

void ConfigTreeItem::writeOutput(FeatureTree * F, int i, int nsteps, std::vector<unsigned int> cacheIndex, std::vector<std::string> flags)
{
    if(getStringData("file_name","file_not_found") == std::string("file_not_found"))
        return ;

    std::string path  = getStringData("file_name","output") ;
    for(size_t j = 0 ; j < flags.size() ; j++)
        path += flags[j] ;

    if(i == 1)
    {
        std::fstream out ;
        out.open( path.c_str(), std::ios::out) ;
        out.close() ;
    }
    if(getChild("time_step")->isAtTimeStep(i, nsteps))
    {
        std::fstream out ;
        out.open( path.c_str(), std::ios::out | std::ios::app) ;
        std::cout << F->getCurrentTime() << "\t" ;
        out << F->getCurrentTime() << "\t" ;
        std::string instant = getStringData("instant","NOW") ;
        double time = (int) (instant == std::string("AFTER")) - (int) (instant == std::string("BEFORE")) ;
        std::vector<ConfigTreeItem *> fields = getAllChildren("field") ;
        for(size_t i = 0 ; i < fields.size() ; i++)
        {
            bool isFieldType = true ;
            FieldType ft = Enum::getFieldType( fields[i]->getStringData(), &isFieldType ) ;
            if(isFieldType)
            {
                Vector f = F->getAverageField( ft, time ) ;
                for(size_t j = 0 ; j < f.size() ; j++)
                {
                    std::cout << f[j] << "\t" ;
                    out << f[j] << "\t" ;
                }
            }
            else
            {
                double f = F->get2DMesh()->getField( fields[i]->getStringData(), -1, fields[i]->getStringData("correction_factor", "correction_factor") ) ;
                std::cout << f << "\t" ;
                out << f << "\t" ;
            }
        }
        std::vector<ConfigTreeItem *> families = getAllChildren("inclusions") ;
        for(size_t i = 0 ; i < families.size() ; i++)
        {
            int index = families[i]->getData("index", 0) ;
            std::vector<ConfigTreeItem *> ffields = families[i]->getAllChildren("field") ;
            for(size_t i = 0 ; i < ffields.size() ; i++)
            {
                bool isFieldType = true ;
                FieldType ft = Enum::getFieldType( ffields[i]->getStringData(), &isFieldType ) ;
                if(isFieldType)
                {
                    Vector f = F->get2DMesh()->getField( ft, cacheIndex[index], (int) (instant == std::string("AFTER")) - (int) (instant == std::string("BEFORE")) ) ;
                    for(size_t j = 0 ; j < f.size() ; j++)
                    {
                        std::cout << f[j] << "\t" ;
                        out << f[j] << "\t" ;
                    }
                }
                else
                {
                    double f = F->get2DMesh()->getField( ffields[i]->getStringData(), cacheIndex[index], ffields[i]->getStringData("correction_factor", "correction_factor") ) ;
                    std::cout << f << "\t" ;
                    out << f << "\t" ;
                }
            }

        }

        std::vector<ConfigTreeItem *> edges = getAllChildren("edge") ;
        for(size_t i = 0 ; i < edges.size() ; i++)
        {
            BoundingBoxPosition pos = Enum::getBoundingBoxPosition(edges[i]->getStringData("position", "BOTTOM")) ;
            std::vector<ConfigTreeItem *> ffields = edges[i]->getAllChildren("field") ;
            for(size_t i = 0 ; i < ffields.size() ; i++)
            {
                bool isFieldType = true ;
                FieldType ft = Enum::getFieldType( ffields[i]->getStringData(), &isFieldType ) ;
                if(isFieldType)
                {
                    Vector f = F->getAverageFieldOnBoundary( pos, ft, (int) (instant == std::string("AFTER")) - (int) (instant == std::string("BEFORE")), getData("field_index",0) ) ;
                    for(size_t j = 0 ; j < f.size() ; j++)
                    {
                        std::cout << f[j] << "\t" ;
                        out << f[j] << "\t" ;
                    }
                }
                else
                {
                    double f = F->getAverageFieldOnBoundary( pos, ffields[i]->getStringData(), (int) (instant == std::string("AFTER")) - (int) (instant == std::string("BEFORE")) ) ;
                    std::cout << f << "\t" ;
                    out << f << "\t" ;
                }
            }

        }

        std::vector<ConfigTreeItem *> nodes = getAllChildren("point") ;
        for(size_t i = 0 ; i < nodes.size() ; i++)
        {
            Point * p = new Point( nodes[i]->getData("x") , nodes[i]->getData("y") , 0 , ((int) (instant == std::string("AFTER")) - (int) (instant == std::string("BEFORE"))) ) ;
            std::vector<ConfigTreeItem *> ffields = nodes[i]->getAllChildren("field") ;
            for(size_t i = 0 ; i < ffields.size() ; i++)
            {
                bool isFieldType = true ;
                FieldType ft = Enum::getFieldType( ffields[i]->getStringData(), &isFieldType ) ;
                if(isFieldType)
                {
                    Vector f = F->getField( ft, p , -1 , false, true) ;
                    for(size_t j = 0 ; j < f.size() ; j++)
                    {
                        std::cout << f[j] << "\t" ;
                        out << f[j] << "\t" ;
                    }
                }

            }
            delete p ;

        }

        nodes.clear() ;
        nodes = getAllChildren("point_list") ;
        for(size_t i = 0 ; i < nodes.size() ; i++)
        {
            Vector x = ConfigTreeItem::readLineAsVector( nodes[i]->getStringData("x") ) ;
            Vector y = ConfigTreeItem::readLineAsVector( nodes[i]->getStringData("y") ) ;
            for(size_t xi = 0 ; xi < x.size() ; xi++)
            {
                for(size_t yi = 0 ; yi < y.size() ; yi++)
                {
                    Point * p = new Point( x[xi] , y[yi] , 0 , ((int) (instant == std::string("AFTER")) - (int) (instant == std::string("BEFORE"))) ) ;
                    std::vector<ConfigTreeItem *> ffields = nodes[i]->getAllChildren("field") ;
                    for(size_t i = 0 ; i < ffields.size() ; i++)
                    {
                        bool isFieldType = true ;
                        FieldType ft = Enum::getFieldType( ffields[i]->getStringData(), &isFieldType ) ;
                        if(isFieldType)
                        {
                            Vector f = F->getField( ft, p , -1 , false, true) ;
                            for(size_t j = 0 ; j < f.size() ; j++)
                            {
                                std::cout << f[j] << "\t" ;
                                out << f[j] << "\t" ;
                            }
                        }

                    }
                    delete p ;
                }
            }
        }

        std::cout << std::endl ;
        out << std::endl ;
        out.close() ;
    }
}

void ConfigTreeItem::exportTriangles(FeatureTree * F, int i, int nsteps, std::vector<std::string> flags)
{
    if(getStringData("file_name","file_not_found") == std::string("file_not_found"))
        return ;

    std::string path  = getStringData("file_name","triangles") ;
    for(size_t j = 0 ; j < flags.size() ; j++)
        path += flags[j] ;
    path += (itoa(i)) ;

    if(getChild("time_step")->isAtTimeStep(i, nsteps))
    {
        std::string instant = getStringData("instant","NOW") ;
        TriangleWriter trg( path , F, (int) (instant == std::string("AFTER")) - (int) (instant == std::string("BEFORE")) ) ;
        std::vector<ConfigTreeItem *> fields = getAllChildren("field") ;
        for(size_t i = 0 ; i < fields.size() ; i++)
        {
            bool isFieldType = true ;
            FieldType f = Enum::getFieldType( fields[i]->getStringData(), &isFieldType ) ;
            bool isTWFieldType = false  ;
            TWFieldType g = Enum::getTWFieldType( fields[i]->getStringData(), &isTWFieldType ) ;
            if(isFieldType)
                trg.getField( f ) ;
            else if(isTWFieldType)
                trg.getField( g ) ;
            else
                trg.getField(  fields[i]->getStringData() ) ;
        }
        trg.write() ;
    }
}

void ConfigTreeItem::exportSvgTriangles(MultiTriangleWriter * trg, FeatureTree * F, int i, int nsteps, std::vector<std::string> flags)
{
    if(!trg)
        return ;

    if(getStringData("svg", "TRUE") == std::string("FALSE"))
    {
        exportTriangles( F, i, nsteps, flags) ;
        return ;
    }

    if(getChild("time_step")->isAtTimeStep(i, nsteps))
    {
        std::string instant = getStringData("instant","NOW") ;
        trg->reset(F, (int) (instant == std::string("AFTER")) - (int) (instant == std::string("BEFORE"))) ;
        std::vector<ConfigTreeItem *> fields = getAllChildren("field") ;
        for(size_t i = 0 ; i < fields.size() ; i++)
        {
            bool isFieldType = true ;
            FieldType f = Enum::getFieldType( fields[i]->getStringData(), &isFieldType ) ;
            bool isTWFieldType = false ;
            TWFieldType g = Enum::getTWFieldType( fields[i]->getStringData(), &isTWFieldType ) ;
            if(isFieldType)
            {
                trg->getField( f ) ;
            }
            else if(isTWFieldType)
            {
                trg->getField( g ) ;
            }
            else
                trg->getField(  fields[i]->getStringData() ) ;
        }
        trg->append() ;
        trg->writeSvg() ;
    }
}





#ifdef __WIN32
void ConfigTreeItem::makeWindowsPath()
{
    for(size_t i = 0 ; i < children.size() ; i++)
        children[i]->makeWindowsPath() ;

    if(label != "function")
    {
        std::string path = str ;
        size_t backslash = path.find("/") ;
        while(backslash < std::string::npos)
        {
            std::string left = path.substr(0, backslash) ;
            std::string right = path.substr(backslash+1) ;
            path = left+"\\"+right ;
            backslash = path.find("/") ;
        }
        str = path ;
        std::cout << path << std::endl ;
    }

}
#endif

}
