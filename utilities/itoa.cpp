// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "itoa.h"
#include <string>
#include <algorithm>
#include <cmath>


std::string itoa(int value, int base) {

	enum { kMaxDigits = 35 };
	std::string buf;
// 	buf.reserve( kMaxDigits ); // Pre-allocate enough space.

	// check that the base if valid
	if (base < 2 || base > 16) return buf;

	int quotient = value;
	
	
	// Translating number to string with base:
	do {
		buf += "0123456789abcdef"[ std::abs( quotient % base ) ];
		quotient /= base;
	} while ( quotient );
	
	// Append the negative sign for base 10
	if ( value < 0 && base == 10) buf += std::string("-");
	
	std::reverse( buf.begin(), buf.end() );
	
	return buf;
	
}

int ctoi(char c)
{
	if(c == '0')
		return 0 ;
	if(c == '1')
		return 1 ;
	if(c == '2')
		return 2 ;
	if(c == '3')
		return 3 ;
	if(c == '4')
		return 4 ;
	if(c == '5')
		return 5 ;
	if(c == '6')
		return 6 ;
	if(c == '7')
		return 7 ;
	if(c == '8')
		return 8 ;
	if(c == '9')
		return 9 ;
	return 0 ;
}
