// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//------------------------------------------------------------------------------
// Date:   17 Feb 2005
// Update: 17 Feb 2005
//------------------------------------------------------------------------------

#ifndef CLASS_NIFS_C_FACE_BND_WALL_NEU_H
#define CLASS_NIFS_C_FACE_BND_WALL_NEU_H

#include "face_bnd_wall.h"

namespace NIFS{    

//--- INTERFACE ---

class c_face_bnd_wall_neu : public c_face_bnd_wall {
	public:
		c_face_bnd_wall_neu();
		virtual ~c_face_bnd_wall_neu() {}

	protected:
		virtual void calc_coef_1();
		virtual void calc_coef_c();
		virtual void calc_bc();
		virtual void calc_state();

	protected:
		double m_dt_dn;	//--- wall temperature ---
};

} //--- namespace NIFS

#endif //--- CLASS_NIFS_C_FACE_BND_WALL_NEU_H

