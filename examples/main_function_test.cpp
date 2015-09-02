// Author: Alain Giorla <alain.b.giorla@gmail.com>, (C) 2005-2014
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../polynomial/vm_function_base.h"
#include "../utilities/parser.h"

#include <fstream>
#include <cmath>

using namespace Amie ;

std::vector<double> range2Vector(std::string range) 
{
	std::string begin = range.substr(1, range.find_first_of(":")-1) ;
	std::string mid = range.substr(range.find_first_of(":")+1, range.find_last_of(":")-range.find_first_of(":")-1 ) ;
	std::string end = range.substr(range.find_last_of(":")+1, range.length()-range.find_last_of(":")-2 ) ;

	double start = atof(begin.c_str()) ;
	double step = std::abs(atof(mid.c_str())) ;
	double finish = atof(end.c_str()) ;

	if(step == 0)
		step = 1 ;

	std::vector<double> v ;
	v.push_back(std::min(start, finish)) ;
	for(double d = std::min(start, finish)+step ; d <= std::max(start, finish) ; d += step)
		v.push_back(d) ;

	return v ;
}

int main(int argc, char *argv[])
{
	CommandLineParser parser("Evaluate the value of an AMIE function at different coordinates and print the result in the console") ;
	parser.addArgument("\"function\"", "x", "AMIE function written in RPN (quotes are mandatory)") ;
	parser.addString("--range-x", "[0:0:0]", "range of the x coordinates formatted as [start:step:end]") ;
	parser.addString("--range-y", "[0:0:0]", "range of the y coordinates formatted as [start:step:end]") ;
	parser.addString("--range-z", "[0:0:0]", "range of the z coordinates formatted as [start:step:end]") ;
	parser.addString("--range-t", "[0:0:0]", "range of the t coordinates formatted as [start:step:end]") ;
	parser.addString("--output", "", "writes output in specified file instead of console") ;
	parser.addValue("--print-limit", 1e6, "maximum number of values printed (default 1e6)") ;
	parser.addFlag("--print-vm", false, "print instructions of the AMIE virtual machine") ;
	parser.addFlag("--print-roots", false, "print roots of token list of the expression parser") ;
	parser.addFlag("--print-parser", false, "print token list of the expression parser") ;
	parser.addFlag("--rpn", false, "parse function in Reverse Polish Notation") ;
	parser.disableFeatureTreeArguments() ;
	
	parser.parseCommandLine(argc, argv) ;
	std::vector<double> x_ = range2Vector( parser.getString("--range-x") ) ;
	std::vector<double> y_ = range2Vector( parser.getString("--range-y") ) ;
	std::vector<double> z_ = range2Vector( parser.getString("--range-z") ) ;
	std::vector<double> t_ = range2Vector( parser.getString("--range-t") ) ;
	std::string f_ = parser.getStringArgument("\"function\"") ;
	std::string file = parser.getString("--output") ;
	size_t counter = parser.getValue("--print-limit") ;
	bool instruction = parser.getFlag("--print-vm") ;
	bool rpn = parser.getFlag("--rpn") ;
	bool roots = parser.getFlag("--print-roots") ;
	bool parse = parser.getFlag("--print-parser") ;

	Function test ;
	if(rpn)
		test = Function(f_.c_str()) ;
	else
	{
		test = FunctionParser::getFunction( f_ ) ;
		if( roots || parse)
		{
			FunctionParser p( f_ ) ;
			if(parse)
				p.print() ;
			if(roots)
				p.printRoots() ;
		}
	}

	VirtualMachine vm ;

	if(instruction)
		vm.print(test) ;

	size_t c = 0 ;

	if(file.size() > 0)
	{
		std::fstream stream ;
		stream.open( file.c_str(), std::ios::out ) ;
		for( size_t x = 0 ; x<x_.size() && ++c<counter ; x++)
		{
			for( size_t y = 0 ; y<y_.size()&& ++c<counter  ; y++)
			{
				for( size_t z = 0 ; z<z_.size() && ++c<counter ; z++)
				{
					for( size_t t = 0 ; t<t_.size() && ++c<counter ; t++)
					{
						stream << x_[x] << "\t" << y_[y] << "\t" << z_[z] << "\t" << t_[t] << "\t" << vm.eval( test, x_[x], y_[y], z_[z], t_[t]) << std::endl ;
					}
				}
			}
		}
	}
	else
	{
		for( size_t x = 0 ; x<x_.size() && ++c<counter ; x++)
		{
			for( size_t y = 0 ; y<y_.size() && ++c<counter ; y++)
			{
				for( size_t z = 0 ; z<z_.size() && ++c<counter ; z++)
				{
					for( size_t t = 0 ; t<t_.size() && ++c<counter ; t++)
					{
						std::cout << x_[x] << "\t" << y_[y] << "\t" << z_[z] << "\t" << t_[t] << "\t" << vm.eval( test, x_[x], y_[y], z_[z], t_[t]) << std::endl ;
					}
				}
			}
		}
	}
	


    return 0 ;
}

