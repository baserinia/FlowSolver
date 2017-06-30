//	Copyright (C) 2006  Amir R. Baserinia
//
// 	This file is part of SMAIF (Simple Mesh Adaptor for Incompressible Flows)
//
//	SMAIF is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//
//	SMAIF is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with SMAIF; if not, write to the Free Software
//	Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

//-----------------------------------------------------------------------------
// Created on:	28 Mar 2004
// Last update:	17 Feb 2005
//-----------------------------------------------------------------------------

#ifndef CLASS_NIFS_C_VECTOR_2D_H
#define CLASS_NIFS_C_VECTOR_2D_H

#include "vector_n.h"

namespace NIFS {

//--- INTERFACE ---

class c_vector_2d : public vector_n< double , 2 > {
    public:
		c_vector_2d( void );
		c_vector_2d( double r );
		c_vector_2d( const c_vector_2d& vec );
		c_vector_2d( double x, double y );
		c_vector_2d( const vector_n< double,2>& vec );
		~c_vector_2d() { }

		//-- inherited
		using vector_n<double,2>::operator=;
		using vector_n<double,2>::operator*;

		//--- overwritten		
		c_vector_2d operator+( const c_vector_2d& vec ) const;
		c_vector_2d operator-( const c_vector_2d& vec ) const;
		c_vector_2d operator-() const;
		c_vector_2d operator*( const double& r ) const;
		c_vector_2d operator/( const double& r ) const;
		c_vector_2d get_unit() const;

		//--- specialized
		double operator^( const c_vector_2d& vec );
		double get_x() const { return (*this)(0); }
		double get_y() const { return (*this)(1); }
		double get_angle( const c_vector_2d& vec );
		void set_x( double x ) { (*this)[0] = x; }
		void set_y( double y ) { (*this)[1] = y; }
		void set_vector( double x, double y ) { set_x( x ); set_y( y ); }
		c_vector_2d rotate_cw() const;
		c_vector_2d rotate_ccw() const;
		c_vector_2d rotate( double theta ) const;        
};


} //--- namespace NIFS ---

#endif //---  CLASS_NIFS_C_VECTOR_2D_H ---

