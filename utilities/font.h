
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef FONT_H
#define FONT_H

#include <ostream>

namespace Amie
{


  typedef enum {
        DEFAULT    = 0,
        BOLD       = 1,
        DIM        = 2,
        UNDERLINED = 4,
        BLACK      = 30,
        RED        = 31,
        GREEN      = 32,
        YELLOW     = 33,
        BLUE       = 34,
        MAGENTA    = 35,
        CYAN       = 36,
        WHITE      = 97,
    } FontCode ;

    struct Font {
        std::vector<FontCode> codes ;

        Font(FontCode pCode = DEFAULT) { codes.push_back(pCode) ; }
        Font(FontCode pCode1, FontCode pCode2) { codes.push_back(pCode1) ; codes.push_back(pCode2) ; }
        Font(FontCode pCode1, FontCode pCode2, FontCode pCode3) { codes.push_back(pCode1) ; codes.push_back(pCode2) ;  codes.push_back(pCode3) ; }


        friend std::ostream&

        operator<<(std::ostream& os, const Font& mod) {
            if( mod.codes.size() == 1 )
                return os << "\033[" << mod.codes[0] << "m" ;
            if( mod.codes.size() == 2 )
                return os << "\033[" << mod.codes[0] << ";" << mod.codes[1] << "m" ;
            if( mod.codes.size() == 3 )
                return os << "\033[" << mod.codes[0] << ";" << mod.codes[2] << ";" << mod.codes[3] << "m" ;
            return os ;
        }
    };

}


#endif
