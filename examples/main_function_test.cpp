// Author: Alain Giorla <alain.b.giorla@gmail.com>, (C) 2005-2014
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../polynomial/vm_function_base.h"

#include <fstream>
#include <cmath>

using namespace Amie ;

int main(int argc, char *argv[])
{
	if(argc < 3)
	{
		std::cout << "AMIE function tester" << std::endl ;
		std::cout << "usage:    ./function_test \"function\" --range x=[min:increment:max] y=[min:increment:max] z=[min:increment:max] t=[min:increment:max] --output filename" << std::endl ;
		std::cout << "x,y,z,t coordinates can be omitted (default value = 0)" << std::endl ;
		std::cout << "output can be omitted (results will be printed in the console)" << std::endl ;
		return 0 ;
	}


	std::string f_ = std::string(argv[1]) ;
	std::vector<double> x_ ;
	std::vector<double> y_ ;
	std::vector<double> z_ ;
	std::vector<double> t_ ;
	std::string file ;
	bool function = true ;
	bool range = false ;
	bool output = false ;

	int counter = 1e6 ;

	for(int i = 2 ; i < argc ; i++)
	{
		if(function)
		{
			if(std::string(argv[i]) != std::string("--range"))
			{
				f_.append(" ") ;
				f_.append(std::string(argv[i])) ;
				std::cout << f_ << std::endl ;
			}
			else
			{
				function = false ;
				range = true ;
			}
		}
		else if(range) 
		{
			if(std::string(argv[i]) != std::string("--output"))
			{
				std::string token = std::string(argv[i]) ;
				char coord = token[0] ;
				bool exponent = (token.find("10^") == 2) ;
				double first = 0. ;
				double last = 0. ;
				double incr = 0. ;
				int left = token.find('[') ;
				int firstcolon = token.find(':') ;
				int secondcolon = token.find(':', firstcolon+1) ;
				int right = token.find(']') ;
				first = atof( token.substr(left+1, firstcolon-left ).c_str() ) ;
				incr = atof( token.substr(firstcolon+1, secondcolon-firstcolon ).c_str() ) ;
				last = atof( token.substr(secondcolon+1, right-secondcolon ).c_str() ) ;
				if(first > last)
				{
					double tmp = last ;
					last = first ;
					first = tmp ;
				} 
				if(incr < 0)
					incr = -incr ;
			
				if(!exponent)
				{
					if(coord == 'x')
					{
						for(double a = first ; a <= last ; a+= incr)
							x_.push_back(a) ;
					}
					if(coord == 'y')
					{
						for(double a = first ; a <= last ; a+= incr)
							y_.push_back(a) ;
					}
					if(coord == 'z')
					{
						for(double a = first ; a <= last ; a+= incr)
							z_.push_back(a) ;
					}
					if(coord == 't')
					{
						for(double a = first ; a <= last ; a+= incr)
							t_.push_back(a) ;
					}
				}
				else
				{
					if(coord == 'x')
					{
						for(double a = first ; a <= last ; a+= incr)
							x_.push_back(std::pow(10.,a)) ;
					}
					if(coord == 'y')
					{
						for(double a = first ; a <= last ; a+= incr)
							y_.push_back(std::pow(10.,a)) ;
					}
					if(coord == 'z')
					{
						for(double a = first ; a <= last ; a+= incr)
							z_.push_back(std::pow(10.,a)) ;
					}
					if(coord == 't')
					{
						for(double a = first ; a <= last ; a+= incr)
							t_.push_back(std::pow(10.,a)) ;
					}
				}
			}
			else
			{
				range = false ;
				output = true ;
			}
		}
		else if(output)
		{
			file = std::string(argv[i]) ;
		}

	}

	if(x_.size() == 0)
		x_.push_back(0.);
	if(y_.size() == 0)
		y_.push_back(0.);
	if(z_.size() == 0)
		z_.push_back(0.);
	if(t_.size() == 0)
		t_.push_back(0.);

	VirtualMachine vm ;
	Function test(f_.c_str()) ;

	vm.print(test) ;

	int c = 0 ;

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

