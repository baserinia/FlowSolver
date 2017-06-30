// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//----------------------------------------------------------------------------
// Date:   29 Mar 2004
// Update: 25 Feb 2005
//----------------------------------------------------------------------------

#ifndef CLASS_NIFS_C_POINT_2D_H
#define CLASS_NIFS_C_POINT_2D_H

#include "vector_2d.h"

namespace NIFS{    

//--- INTERFACE ---

class c_point_2d {
	public:
		c_point_2d() {}
		c_point_2d( const c_point_2d& point ) { m_pv = point.m_pv; }
		c_point_2d( const c_vector_2d& vec ) { m_pv = vec; }
		c_point_2d( double x, double y ) { set_pos( x, y ); }
		~c_point_2d() {}
		void set_pos( double x, double y ) { m_pv.set_vector( x, y ); }
		void set_pos( const c_vector_2d &vec ) { m_pv = vec; }
		const c_vector_2d& get_pos() const { return m_pv; }
		void operator=( const c_point_2d &point ) { m_pv = point.m_pv; }
	private:
		c_vector_2d m_pv;  // position vector
};

} // namespace NIFS

#endif // CLASS_NIFS_POINT_2D_H

