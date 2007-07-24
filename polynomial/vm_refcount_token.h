//
// C++ Interface: vm_refcount_token
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef REF_COUNTED_TOKEN_H
#define REF_COUNTED_TOKEN_H

#include<iostream>
#include<valarray>

#include "vm_token.h"

namespace Mu
{

class RefCountedTokenInternal
{
private:
	const Token * t ;
	size_t count ;
	
	RefCountedTokenInternal( const RefCountedTokenInternal & ) ;

	RefCountedTokenInternal &operator=( const RefCountedTokenInternal & )
	{
		return *this ;
	}
	
public:

	void eval(Context &context) const
	{
		t->eval(context) ;
	}
	
	RefCountedTokenInternal(const Token * tok) 
	{
		t = tok ;
		count = 1 ;
	}
	
	~RefCountedTokenInternal()
	{
		delete t ;
	}
	
	const Token * token() const 
	{
		return t ;
	}
	
	inline bool dereferenced()
	{
		count-- ;
		if (count == 0)
		{
			return true ;
		}
		else
		{
			return false ;
		}
	}
	
	inline void referenced()
	{
		count++ ;
	}
	
	size_t getCount() const
	{
		return count ;
	}
} ;

class RefCountedToken
{
protected:
	RefCountedTokenInternal * t ; 
	void referenced() 
	{
		if(t!=NULL)
			t->referenced() ;
	}
	
	void dereferenced() 
	{
		if ( t!=NULL && t->dereferenced() )
		{
			delete t;
			t = NULL ;
		}
	}
	
public:
	
	RefCountedToken() : t(NULL) { } ; 
	
	
	void eval(Context &context) const
	{
		t->eval(context) ;
	}
	
	~RefCountedToken()
	{
		dereferenced() ;
	}
	
	RefCountedToken(const Token * token)
	{
		t = new RefCountedTokenInternal(token) ;
	}
	
	RefCountedToken(const RefCountedToken & token)
	{
		t = token.t ;
		referenced() ;
	}
	
	RefCountedToken & operator = ( const RefCountedToken & token)
	{
		if(t == NULL)
		{
			t = token.t ;
			referenced() ;
		}
		else
		{
			dereferenced();
			t = token.t ;
			referenced() ;
		}
		return *this ;
	}
	
	RefCountedToken & operator = (const Token * token)
	{
		if ( t && t->getCount() == 1 ) {
			referenced() ;
		}
		else 
		{
			dereferenced();
			t = new RefCountedTokenInternal( token );
		}
		return *this;
	}
	
	bool operator== ( const RefCountedToken& p ) const
	{
		return (const Token *)(*this) == (const Token*) p;
	}
	
	bool operator!= ( const RefCountedToken & p ) const
	{
		return !( *this == p );
	}
	
	const Token * operator-> () const
	{
		return (const Token *)(t?t->token():0) ; 
	} 
	
	operator const Token * () const
	{
		return (const Token *)(t?t->token():0) ; 
	}

} ;

typedef std::valarray<RefCountedToken> ByteCode  ;

} ;

#endif
