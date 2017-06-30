// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//----------------------------------------------------------------------------
// Created on:   12 Mar 2005
// Last update:  25 Mar 2005
//----------------------------------------------------------------------------

#ifndef TEMPLATE_CLASS_NIFS_MATRIX_BLOCK_H
#define TEMPLATE_CLASS_NIFS_MATRIX_BLOCK_H

#include "matrix_n.h"

namespace NIFS{    

//--- INTERFACE ---

template < typename T >
class matrix_block : public matrix_n< T , 4 >{
	public:
		matrix_block( void ) : matrix_n< T , 4 >() {}
		matrix_block( const matrix_block<T>& mat ) : matrix_n< T , 4 >( mat ) {}
		matrix_block( const matrix_n<T,4>& mat ) : matrix_n< T , 4 >( mat ) {}
		matrix_block( T t ) : matrix_n< T , 4 >( t ) {}

		//-- inherited
		using matrix_n< T , 4 >::operator=;
		using matrix_n< T , 4 >::operator+=;
		using matrix_n< T , 4 >::operator-=;
		using matrix_n< T , 4 >::operator-;
		using matrix_n< T , 4 >::operator*;
		using matrix_n< T , 4 >::operator[];
		using matrix_n< T , 4 >::operator();
};

} // namespace NIFS

#endif // TEMPLATE_CLASS_NIFS_MATRIX_BLOCK_H

