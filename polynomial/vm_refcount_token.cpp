//
// C++ Implementation: vm_refcount_token
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "vm_refcount_token.h"

using namespace Mu ;

RefCountedTokenInternal::RefCountedTokenInternal( const RefCountedTokenInternal & )
{
// 	std::cout << "*** invalid ***" << std::endl ;
}
	
// RefCountedTokenInternal & RefCountedTokenInternal::operator=( const RefCountedTokenInternal & )
// {
// // 	std::cout << "*** invalid ***" << std::endl ;
// 	return *this ;
// }
	

// RefCountedTokenInternal::RefCountedTokenInternal(const Token * tok)
// {
// 	t = tok ;
// 	count = 1 ;
// }
	
// RefCountedTokenInternal::~RefCountedTokenInternal()
// {
// 	delete t ;
// }
// 	
// const Token * RefCountedTokenInternal::token() const 
// {
// 	return t ;
// }
// 	
// inline bool RefCountedTokenInternal::dereferenced()
// {
// 	count-- ;
// 	if (count == 0)
// 	{
// 		return true ;
// 	}
// 	else
// 	{
// 		return false ;
// 	}
// }
// 	
// inline void RefCountedTokenInternal::referenced()
// {
// 	count++ ;
// }
// 
// size_t RefCountedTokenInternal::getCount() const
// {
// 	return count ;
// }

// void RefCountedToken::referenced()
// {
// 	if(t!=nullptr)
// 		t->referenced() ;
// }
// 
// void RefCountedToken::dereferenced()
// {
// 	if ( t!=nullptr && t->dereferenced() )
// 	{
// 		delete t;
// 		t = nullptr ;
// 	}
// }
	

// RefCountedToken::RefCountedToken() : t(nullptr) { } ; 
// 	
// RefCountedToken::~RefCountedToken()
// {
// 	dereferenced() ;
// }

// RefCountedToken::RefCountedToken(const Token * token)
// {
// 	t = new RefCountedTokenInternal(token) ;
// }
// 
// RefCountedToken::RefCountedToken(const RefCountedToken & token)
// {
// 	t = token.t ;
// 	referenced() ;
// }
// 
// RefCountedToken & RefCountedToken::operator = ( const RefCountedToken & token)
// {
// 	if(t == nullptr)
// 	{
// 		t = token.t ;
// 		referenced() ;
// 	}
// 	else
// 	{
// 		dereferenced();
// 		t = token.t ;
// 		referenced() ;
// 	}
// 	return *this ;
// }
// 	
// RefCountedToken & RefCountedToken::operator = (const Token * token)
// {
// 	if ( t && t->getCount() == 1 ) {
// 		referenced() ;
// 	}
// 	else 
// 	{
// 		dereferenced();
// 		t = new RefCountedTokenInternal( token );
// 	}
// 	return *this;
// }
// 
// bool RefCountedToken::operator== ( const RefCountedToken& p ) const
// {
// 	return (const Token *)(*this) == (const Token*) p;
// }
// 
// bool RefCountedToken::operator!= ( const RefCountedToken & p ) const
// {
// 	return !( *this == p );
// }
// 
// const Token * RefCountedToken::operator-> () const
// {
// return (const Token *)(t?t->token():0) ; 
// } 
// 	
// RefCountedToken::operator const Token * () const
// {
// 	return (const Token *)(t?t->token():0) ; 
// }
