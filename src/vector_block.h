// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//----------------------------------------------------------------------------
// Created on:   28 Feb 2005
// Last update:  12 Mar 2005
//----------------------------------------------------------------------------

#ifndef TEMPLATE_CLASS_NIFS_C_VECTOR_BLOCK_H
#define TEMPLATE_CLASS_NIFS_C_VECTOR_BLOCK_H

#include "vector_n.h"

namespace NIFS{    
	
template < typename T >
class vector_block : public vector_n< T , 4 > {
	public:
		vector_block( void ) : vector_n< T , 4 >() {}
		vector_block( const vector_block<T>& vec ) : vector_n< T , 4 >( vec ) {}
		vector_block( const vector_n<T,4>& vec ) : vector_n< T , 4 >( vec ) {}
		vector_block( T r ) : vector_n< T , 4 >( r ) {}

		//-- inherited
		using vector_n< T , 4 >::operator=;
		using vector_n< T , 4 >::operator+=;
		using vector_n< T , 4 >::operator-=;
		using vector_n< T , 4 >::operator+;
		using vector_n< T , 4 >::operator-;
		using vector_n< T , 4 >::operator*;
		using vector_n< T , 4 >::mult_elem;
		using vector_n< T , 4 >::operator[];
		using vector_n< T , 4 >::operator();

		T get_pres( void ) { return (*this)(0); }	// pressure
		T get_velu( void ) { return (*this)(1); }	// velocity_u
		T get_velv( void ) { return (*this)(2); }	// velocity-v
		T get_temp( void ) { return (*this)(3); }	// temperature
	
		void set_pres( const T& pres ) { (*this)[0] = pres; }
		void set_velu( const T& velu ) { (*this)[1] = velu; }
		void set_velv( const T& velv ) { (*this)[2] = velv; }
		void set_temp( const T& temp ) { (*this)[3] = temp; }
};

} // namespace NIFS

#endif // TEMPLATE_CLASS_NIFS_C_VECTOR_BLOCK_H

