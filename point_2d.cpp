// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//----------------------------------------------------------------------------
// Date:   29 Mar 2004
// Update: 25 Feb 2005
//----------------------------------------------------------------------------

#ifndef CLASS_NIFS_POINT_2D_H
#define CLASS_NIFS_POINT_2D_H

#include "vector_2d.h"

namespace NIFS{    

class Tpoint_2d
{
public:
    Tpoint_2d() {}
    Tpoint_2d( const Tpoint_2d& point ) { m_pv = point.mpv; }	
    Tpoint_2d( const CVector2D& vec ) { m_pv = vec; }	
    Tpoint_2d( double x, double y ) { set_pos( x, y ); }
    
	~Tpoint_2d() {}
	
    void set_pos( double x, double y ) { m_pv.set_vector( x, y ); }    	
    void set_pos( const Tvector_2d &vec ) { m_pv = vec; }

    Tvector_2d& get_pos() const { return m_pv; }	

    void operator=( const Tpoint_2d &point ) { m_pv = point.m_pv; }

private:
    Tvector_2d m_pv;  // position vector	      
};

} // namespace NIFS

#endif // CLASS_NIFS_POINT_2D_H

