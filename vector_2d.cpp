// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//-----------------------------------------------------------------------------
// Created on:	29 Mar 2004
// Last update:	29 mar 2005
//-----------------------------------------------------------------------------

#include <cmath>
#include "vector_2d.h"

//--- IMPLEMENTATION ---

NIFS::c_vector_2d::c_vector_2d(  ) : NIFS::vector_n< double , 2 >( )
{
	
}

NIFS::c_vector_2d::c_vector_2d( double r ) : NIFS::vector_n< double , 2 >( r ) 
{
	
}

NIFS::c_vector_2d::c_vector_2d( const NIFS::c_vector_2d& vec ) : vector_n< double , 2 >( vec )
{
	
}

NIFS::c_vector_2d::c_vector_2d( double x, double y ) 
{ 
	set_vector( x , y ); 
}

NIFS::c_vector_2d::c_vector_2d( const NIFS::vector_n< double,2>& vec ) : 
			   vector_n< double , 2 >( vec ) 
{
	
}

NIFS::c_vector_2d
NIFS::c_vector_2d::operator+( const NIFS::c_vector_2d& vec ) const
{
	return vector_n<double,2>::operator+( vec ); 
}

NIFS::c_vector_2d
NIFS::c_vector_2d::operator-( const NIFS::c_vector_2d& vec ) const
{
	return vector_n<double,2>::operator-( vec ); 
}


NIFS::c_vector_2d
NIFS::c_vector_2d::operator-( void ) const
{
	return vector_n<double,2>::operator-( ); 
}

NIFS::c_vector_2d
NIFS::c_vector_2d::operator*( const double& r ) const
{
	return vector_n<double,2>::operator*( r ); 
}

NIFS::c_vector_2d
NIFS::c_vector_2d::operator/( const double& r ) const
{
	return vector_n<double,2>::operator/( r ); 
}


NIFS::c_vector_2d
NIFS::c_vector_2d::get_unit( void ) const
{
	return vector_n<double,2>::get_unit( ); 
}


double
NIFS::c_vector_2d::operator^( const NIFS::c_vector_2d& vec ) 
{
    return ( (*this)(0) * vec(1) - (*this)(1) * vec(0) );
}  


double
NIFS::c_vector_2d::get_angle( const NIFS::c_vector_2d& vec )
{
    return std::acos( (*this) * vec / ( get_magnit() * vec.get_magnit() ) ); 
}  


NIFS::c_vector_2d
NIFS::c_vector_2d::rotate_cw() const
{
    return c_vector_2d( (*this)(1), -(*this)(0) );
}  


NIFS::c_vector_2d
NIFS::c_vector_2d::rotate_ccw() const
{
    return c_vector_2d( -(*this)(1), (*this)(0) );
}  


NIFS::c_vector_2d
NIFS::c_vector_2d::rotate( double theta ) const
{
    return NIFS::c_vector_2d(	std::cos( theta ) * (*this)(0) - std::sin( theta ) * (*this)(1),
						std::sin( theta ) * (*this)(0) + std::cos( theta ) * (*this)(1) );
}  
