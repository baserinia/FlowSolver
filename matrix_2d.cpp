// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//-----------------------------------------------------------------------------
// Created on:	12 May 2005
// Last update:	12 May 2005
//-----------------------------------------------------------------------------

#include <cmath>
#include "matrix_2d.h"

//--- IMPLEMENTATION ---
NIFS::c_matrix_2d::c_matrix_2d( double r ) : NIFS::matrix_n< double , 2 >( r ) 
{

}

NIFS::c_matrix_2d::c_matrix_2d( double t00, double t01, double t10, double t11 ) 
{
    set_matrix( t00, t01, t10, t11 );
}  

NIFS::c_matrix_2d::c_matrix_2d( const c_vector_2d& v0, const c_vector_2d& v1 ) 
{
    (*this)[0] = v0;
    (*this)[1] = v1;
}  

NIFS::c_matrix_2d::c_matrix_2d( const matrix_n< double,2>& mat ) : matrix_n< double , 2 >( mat )
{

}

NIFS::c_matrix_2d
NIFS::c_matrix_2d::eig_val( void ) const
{
    double delta = std::sqrt( ( (*this)(0)(0) - (*this)(1)(1) ) * 
							  ( (*this)(0)(0) - (*this)(1)(1) ) + 
							  4 * (*this)(0)(1) * (*this)(1)(0) );
    double lambda_1 = ( (*this)(0)(0) + (*this)(1)(1)  - delta ) / 2.0;
	double lambda_2 = ( (*this)(0)(0) + (*this)(1)(1)  + delta ) / 2.0;
	return c_matrix_2d( lambda_1 , 0.0, 0.0, lambda_2 ); 
}  


NIFS::c_matrix_2d
NIFS::c_matrix_2d::eig_vec( void ) const
{
	c_matrix_2d ten;
	c_matrix_2d eig_val_tmp = eig_val();
	if ( std::abs( (*this)(0)(0) - eig_val_tmp(0)(0) ) < 1E-15 ) {
		ten[0][0] = 1.0;
		ten[1][0] = 0.0;
	} else 	{
		ten[1][0] = 1.0 / std::sqrt( 1.0 + sqr( (*this)(0)(1) / ( (*this)(0)(0) - eig_val_tmp(0)(0) )));
		ten[0][0] = -(*this)(0)(1) / ( (*this)(0)(0) - eig_val_tmp(0)(0) ) * ten(1)(0);
		
	}
	if ( std::abs( (*this)(1)(1) - eig_val_tmp(1)(1) ) < 1E-15 ) {
		ten[0][1] = 0.0;
		ten[1][1] = 1.0;
	} else {
		ten[0][1] = 1.0 / std::sqrt( 1.0 + sqr( (*this)(1)(0) / ( (*this)(1)(1) - eig_val_tmp(1)(1) )));
		ten[1][1] = -(*this)(1)(0) / ( (*this)(1)(1) - eig_val_tmp(1)(1) ) * ten(0)(1);
	}
	return ten;
}  


void 
NIFS::c_matrix_2d::set_matrix( double t00, double t01, double t10, double t11)
{
	(*this)[0][0] = t00;
	(*this)[0][1] = t01;
	(*this)[1][0] = t10;
	(*this)[1][1] = t11;
}  


double
NIFS::c_matrix_2d::determinant( void )
{
    return ( (*this)(0)(0) * (*this)(1)(1) - (*this)(1)(0) * (*this)(0)(1) );
}

NIFS::c_matrix_2d
NIFS::c_matrix_2d::operator+( const NIFS::c_matrix_2d& mat ) const
{
	return matrix_n<double,2>::operator+( mat ); 
}

NIFS::c_matrix_2d
NIFS::c_matrix_2d::operator-( const NIFS::c_matrix_2d& mat ) const
{
	return matrix_n<double,2>::operator-( mat ); 
}

NIFS::c_vector_2d 
NIFS::c_matrix_2d::operator*( const NIFS::c_vector_2d& vec ) const
{
	return matrix_n<double,2>::operator*( vec );
}

NIFS::c_matrix_2d
NIFS::c_matrix_2d::operator-( void ) const
{
	return matrix_n<double,2>::operator-( ); 
}

NIFS::c_matrix_2d
NIFS::c_matrix_2d::operator*( const double& r ) const
{
	return matrix_n<double,2>::operator*( r ); 
}

NIFS::c_matrix_2d
NIFS::c_matrix_2d::operator/( const double& r ) const
{
	return matrix_n<double,2>::operator/( r ); 
}

NIFS::c_matrix_2d
NIFS::c_matrix_2d::transpose( void ) const
{
	return matrix_n<double,2>::transpose( ); 
}

double
NIFS::c_matrix_2d::contraction( const c_matrix_2d& mat)
{
	return	(*this)(0)(0) * mat(0)(0) +
			(*this)(0)(1) * mat(0)(1) +
			(*this)(1)(0) * mat(1)(0) +
			(*this)(1)(1) * mat(1)(1);
}

//---( Generalized vector product )---------------------------------------------
double
NIFS::c_matrix_2d::gen_dot_prod( const c_vector_2d& v0, const c_vector_2d& v1 ) 
{
	return  (*this)(0)(0) * v0(0) * v1(0) + 
			(*this)(0)(1) * v0(0) * v1(1) + 
			(*this)(1)(0) * v0(1) * v1(0) + 
			(*this)(1)(1) * v0(1) * v1(1); 
}

NIFS::c_matrix_2d
NIFS::c_matrix_2d::log( void ) const
{
  return c_matrix_2d(0.0,0.0,0.0,0.0);                    
}  
