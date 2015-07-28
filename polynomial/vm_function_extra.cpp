//
// C++ Implementation: vm_function_base
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "vm_function_extra.h"
#include "vm_function_base.h"
#include "vm_token.h"

namespace Amie {

Function f_positivity(const Function &f, bool differentiate)
{
    Function ret(f) ;
    ret.byteCode.push_back(TOKEN_OPERATION_POSITIVITY);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
    ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
    ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;


    if(ret.isDifferentiable() && !differentiate)
    {
        ret.setNumberOfDerivatives(0);
    }
    else if(ret.isDifferentiable())
    {
        Function zero("0") ;

        size_t  j = 0 ;
        while(j < ret.getDerivatives().size())
        {
            if( f.isDifferentiable((const Variable) j) )
            {
                Function z("0") ;
                zero.setNumberOfDerivatives( ret.getNumberOfDerivatives() );
                for(size_t k = 0 ; k < ret.getDerivatives().size() ; k++)
                {
                    zero.setDerivative((const Variable) k, z);
                }
                j = ret.getDerivatives().size() ;
            }
        }


        for(size_t i = 0 ; i < ret.getDerivatives().size() ; i++)
            ret.setDerivative((const Variable)i,  zero) ;
    }
    return ret ;
}

Function f_negativity(const Function &f, bool differentiate)
{
    Function ret = f ;
    ret.byteCode.push_back(TOKEN_OPERATION_NEGATIVITY);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
    ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
    ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;

    if(ret.isDifferentiable() && !differentiate)
    {
        ret.setNumberOfDerivatives(0);
    }
    else if (ret.isDifferentiable())
    {
        Function zero("0") ;

        size_t  j = 0 ;
        while(j < ret.getDerivatives().size())
        {
            if( f.isDifferentiable((const Variable) j) )
            {
                Function z("0") ;
                zero.setNumberOfDerivatives( ret.getNumberOfDerivatives() );
                for(size_t k = 0 ; k < ret.getDerivatives().size() ; k++)
                {
                    zero.setDerivative((const Variable) k, z);
                }
                j = ret.getDerivatives().size() ;
            }
        }

        for(size_t i = 0 ; i < ret.getDerivatives().size() ; i++)
            ret.setDerivative((const Variable)i,  zero) ;
    }
    return ret ;
}

Function f_range(const Function &f, double min, double max, bool differentiate)
{
    return f_positivity(f-min, differentiate)*f_negativity(f-max, differentiate) ;
}


Function f_exp(const Function &f)
{
    Function ret = f ;
    ret.byteCode.push_back(TOKEN_OPERATION_EXP);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
    ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
    ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;

    return ret ;
}

Function f_abs(const Function &f, bool differentiate)
{
    Function ret = f ;
    ret.byteCode.push_back(TOKEN_OPERATION_ABS);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
    ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
    ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;


    ret.setNumberOfDerivatives(f.getNumberOfDerivatives()*differentiate);

    if(differentiate && f.isDifferentiable())
    {
        for(int i = 0 ; i < f.getNumberOfDerivatives() ; i++)
        {
            Function d = f.d((const Variable)i)*f_positivity(f,false)-f.d((const Variable)i)*f_negativity(f,false) ;
            if(f.d((const Variable) i).isDifferentiable())
            {
                d.setNumberOfDerivatives( f.d((const Variable) i).getNumberOfDerivatives() );
                for(int k = 0 ; k < d.getNumberOfDerivatives() ; k++)
                {
                    Function g = f.d((const Variable)i).d((const Variable) k)*f_positivity(f,false)-f.d((const Variable)i).d((const Variable) k)*f_negativity(f,false) ;
                    d.setDerivative((const Variable) k, g);
                }
            }
            ret.getDerivatives()[i] =  new Function(d) ;
        }
    }
    return ret ;
}

Function f_log(const Function &f)
{
    Function ret(f) ;
    ret.byteCode.push_back(TOKEN_OPERATION_LOG);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
    ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
    ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;

    return ret ;
}

Function f_sqrt(const Function &f, bool differentiate)
{
    Function ret(f) ;
    ret.byteCode.push_back(TOKEN_OPERATION_SQRT);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
    ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
    ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;


    ret.setNumberOfDerivatives(0);//f.getNumberOfDerivatives()*differentiate) ;
//
// 	if(differentiate && f.isDifferentiable())
// 	{
// 		for(size_t i = 0 ; i < f.getNumberOfDerivatives() ; i++)
// 		{
// 			Function d = f.d((const Variable)i)*0.5/f_sqrt(f,f.isDifferentiable() && f.d((const Variable)i).isDifferentiable()) ;
// 			ret.getDerivatives()[i] = new Function(d) ;
// 		}
// 	}
    return ret ;
}

Function f_atan2(const Function &f0, const Function &f1)
{
    Function ret ;
    concatenateFunctions(f0, f1, ret);
    ret.byteCode.push_back(TOKEN_OPERATION_ATAN2);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
    ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
    ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
    return ret ;
}

Function f_sin(const Function &f)
{
    Function ret = f ;
    ret.byteCode.push_back(TOKEN_OPERATION_SIN);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
    ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
    ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;

    return ret ;
}

Function f_sinh(const Function &f)
{
    Function ret = f ;
    ret.byteCode.push_back(TOKEN_OPERATION_SINH);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
    ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
    ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;

    return ret ;
}


Function f_sign(const Function &f)
{
    Function ret = f ;
    ret.byteCode.push_back(TOKEN_OPERATION_SIGN);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
    ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
    ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;

    return ret ;
}





Function f_cos(const Function &f)
{
    Function ret = f ;
    ret.byteCode.push_back(TOKEN_OPERATION_COS) ;
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
    ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
    ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;

    return ret ;
}

Function f_cosh(const Function &f)
{
    Function ret = f ;
    ret.byteCode.push_back(TOKEN_OPERATION_COSH) ;
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
    ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
    ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;

    return ret ;
}

Function f_tan(const Function &f)
{
    Function ret = f ;
    ret.byteCode.push_back(TOKEN_OPERATION_TAN) ;
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
    ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
    ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;

    return ret ;
}

Function f_tanh(const Function &f)
{
    Function ret = f ;
    ret.byteCode.push_back(TOKEN_OPERATION_TANH) ;
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a.push_back(0);
    ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
    ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
    ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;

    return ret ;
}

}




