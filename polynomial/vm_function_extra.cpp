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

#include "vm_function_base.h"
#include "vm_token.h"


Amie::Function f_positivity(const Amie::Function &f, bool differentiate)
{
	Amie::Function ret(f) ;
	ret.byteCode.push_back(Amie::TOKEN_OPERATION_POSITIVITY);
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
	
	
	if(ret.isDifferentiable() && !differentiate)
	{
		ret.setNumberOfDerivatives(0);
	}
	else if(ret.isDifferentiable())
	{
		Amie::Function zero("0") ;
		
		size_t  j = 0 ;
		while(j < ret.getDerivatives().size())
		{
			if( f.isDifferentiable((const Amie::Variable) j) )
			{
				Amie::Function z("0") ;
				zero.setNumberOfDerivatives( ret.getNumberOfDerivatives() );
				for(size_t k = 0 ; k < ret.getDerivatives().size() ; k++)
				{
					zero.setDerivative((const Amie::Variable) k, z);
				}
				j = ret.getDerivatives().size() ;
			}
		}
		
		
		for(size_t i = 0 ; i < ret.getDerivatives().size() ; i++)
			ret.setDerivative((const Amie::Variable)i,  zero) ;
	}
	return ret ;
}



Amie::Function f_negativity(const Amie::Function &f, bool differentiate)
{
	Amie::Function ret = f ;
	ret.byteCode.push_back(Amie::TOKEN_OPERATION_NEGATIVITY);
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
	
	if(ret.isDifferentiable() && !differentiate)
	{
		ret.setNumberOfDerivatives(0);
	}
	else
	{
		Amie::Function zero("0") ;

		size_t  j = 0 ;
		while(j < ret.getDerivatives().size())
		{
			if( f.isDifferentiable((const Amie::Variable) j) )
			{
				Amie::Function z("0") ;
				zero.setNumberOfDerivatives( ret.getNumberOfDerivatives() );
				for(size_t k = 0 ; k < ret.getDerivatives().size() ; k++)
				{
					zero.setDerivative((const Amie::Variable) k, z);
				}
				j = ret.getDerivatives().size() ;
			}
		}
		
		for(size_t i = 0 ; i < ret.getDerivatives().size() ; i++)
			ret.setDerivative((const Amie::Variable)i,  zero) ;
	}
	return ret ;
}

Amie::Function f_exp(const Amie::Function &f)
{
	Amie::Function ret = f ;
	ret.byteCode.push_back(Amie::TOKEN_OPERATION_EXP);
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
	
	return ret ;
}

Amie::Function f_abs(const Amie::Function &f, bool differentiate)
{
	Amie::Function ret = f ;
	ret.byteCode.push_back(Amie::TOKEN_OPERATION_ABS); 
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
	
	
	ret.setNumberOfDerivatives(f.getNumberOfDerivatives()*differentiate);
	
	if(differentiate && f.isDifferentiable())
	{
		for(int i = 0 ; i < f.getNumberOfDerivatives() ; i++)
		{
			Amie::Function d = f.d((const Amie::Variable)i)*f_positivity(f,false)-f.d((const Amie::Variable)i)*f_negativity(f,false) ;
			if(f.d((const Amie::Variable) i).isDifferentiable())
			{
				d.setNumberOfDerivatives( f.d((const Amie::Variable) i).getNumberOfDerivatives() );
				for(int k = 0 ; k < d.getNumberOfDerivatives() ; k++)
				{
					Amie::Function g = f.d((const Amie::Variable)i).d((const Amie::Variable) k)*f_positivity(f,false)-f.d((const Amie::Variable)i).d((const Amie::Variable) k)*f_negativity(f,false) ;
					d.setDerivative((const Amie::Variable) k, g);
				}
			}
			ret.getDerivatives()[i] =  new Amie::Function(d) ;
		}
	}
	return ret ;
}

Amie::Function f_log(const Amie::Function &f)
{
	Amie::Function ret(f) ;
	ret.byteCode.push_back(Amie::TOKEN_OPERATION_LOG);
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
	
	return ret ;
}

Amie::Function f_sqrt(const Amie::Function &f, bool differentiate)
{
	Amie::Function ret(f) ;
	ret.byteCode.push_back(Amie::TOKEN_OPERATION_SQRT);
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
	

 	ret.setNumberOfDerivatives(0);//f.getNumberOfDerivatives()*differentiate) ;
// 	
// 	if(differentiate && f.isDifferentiable())
// 	{
// 		for(size_t i = 0 ; i < f.getNumberOfDerivatives() ; i++)
// 		{
// 			Amie::Function d = f.d((const Variable)i)*0.5/f_sqrt(f,f.isDifferentiable() && f.d((const Variable)i).isDifferentiable()) ;
// 			ret.getDerivatives()[i] = new Amie::Function(d) ;
// 		}
// 	}
	return ret ;
}

Amie::Function f_atan2(const Amie::Function &f0, const Amie::Function &f1)
{
	Amie::Function ret ;
	Amie::concatenateFunctions(f0, f1, ret);
	ret.byteCode.push_back(Amie::TOKEN_OPERATION_ATAN2);
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
	return ret ;
}

Amie::Function f_sin(const Amie::Function &f)
{
	Amie::Function ret = f ;
	ret.byteCode.push_back(Amie::TOKEN_OPERATION_SIN);
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
	
	return ret ;
}


Amie::Function f_sign(const Amie::Function &f)
{
	Amie::Function ret = f ;
	ret.byteCode.push_back(Amie::TOKEN_OPERATION_SIGN);
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
	
	return ret ;
}





Amie::Function f_cos(const Amie::Function &f)
{
	Amie::Function ret = f ;
	ret.byteCode.push_back(Amie::TOKEN_OPERATION_COS) ;
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
	
	return ret ;
} 






