// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//------------------------------------------------------------------------------
// Date:   10 Apr 2004
// Update: 28 Mar 2005
//------------------------------------------------------------------------------

#ifndef CLASS_NIFS_C_FACE_BND_H
#define CLASS_NIFS_C_FACE_BND_H

#include "face.h"

namespace NIFS{    
	
//--- INTERFACE ---

class c_face_bnd : public c_face {
	public:
		c_face_bnd() : c_face() {}
		virtual ~c_face_bnd() {}
		virtual void update_geometry( void );
		virtual unsigned  get_type() { return 1; }
		virtual double project( const c_vector_2d& pnt, c_vector_2d& prj );
				
	protected:
		virtual void calc_uv12();
		virtual void calc_s12();
		virtual void calc_coef_1() = 0;
		virtual void calc_coef_c() = 0;        
		virtual void calc_df();
		virtual void calc_mfr() = 0;
		virtual void calc_bc() = 0; 
		virtual void calc_state() = 0;
		virtual void calc_error() { m_error = 0.0; }

	private:
		virtual void calc_r2() {}		//--- N/A, only for internal faces
		virtual void calc_rc() {}		//--- N/A, only for internal faces
		virtual void calc_coef_2() {}	//--- N/A, only for internal faces
};

} //--- namespace NIFS

#endif //--- CLASS_NIFS_C_FACE_BND_H ---

