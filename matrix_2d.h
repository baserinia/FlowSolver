// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//----------------------------------------------------------------------------
// Created on: 	25 Jun 2004
// Last update: 17 Feb 2005
//-----------------------------------------------------------------------------

#ifndef CLASS_NIFS_C_MATRIX_2D_H
#define CLASS_NIFS_C_MATRIX_2D_H


#include <cmath>
#include "matrix_n.h"
#include "vector_2d.h"

namespace NIFS{    //--- Nummy Incompressible Flow Solver ---

//--- INTERFACE ---

class c_matrix_2d : public matrix_n< double , 2 > {
	public:
		c_matrix_2d( void ) : matrix_n< double , 2 >() {}
		c_matrix_2d( double r );
		c_matrix_2d( double t00, double t01, double t10, double t11 );
		c_matrix_2d( const c_vector_2d& v0, const c_vector_2d& v1 );
		c_matrix_2d( const matrix_n< double,2>& mat );
		virtual ~c_matrix_2d( void ) { }
		c_matrix_2d eig_val( void ) const; 
		c_matrix_2d eig_vec( void ) const; 
		void set_matrix( double t00, double t01, double t10, double t11);
		double determinant( void ); 
		double contraction( const c_matrix_2d& mat);
		double gen_dot_prod( const c_vector_2d& v0, const c_vector_2d& v1 ); 
		c_matrix_2d log( void ) const; 
		
		//-- inherited
		using matrix_n<double,2>::operator=;
		using matrix_n<double,2>::operator();
		
		//--- overwritten
		c_matrix_2d operator+( const c_matrix_2d& mat ) const;
		c_matrix_2d operator-( const c_matrix_2d& mat ) const;
		c_vector_2d operator*( const c_vector_2d& vec ) const;
		c_matrix_2d operator-() const;
		c_matrix_2d operator*( const double& r ) const;
		c_matrix_2d operator/( const double& r ) const;
    c_matrix_2d transpose( void ) const;
	private:
    	double sqr( double r ) const { return r * r; }
};

} // namespace NIFS

#endif // CLASS_NIFS_C_MATRIX_2D_H

