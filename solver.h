//	Copyright (C) 2006  Amir R. Baserinia
//
// 	This file is part of SMAIF (Simple Mesh Adaptor for Incompressible Flow)
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
//
//------------------------------------------------------------------------------
// Created on:	15 Apr 2004
// Last update:	20 Jun 2006
//------------------------------------------------------------------------------

#ifndef CLASS_NIFS_C_SOLVER_H
#define CLASS_NIFS_C_SOLVER_H

#include "mesh.h"
#include "info.h"
#include "lin_sys.h"
#include "mesh_adap.h"

namespace NIFS {    

// Application Class
class c_solver {
	public:
		c_solver() {}
		~c_solver() {}
		void initialize( c_mesh* mesh_ptr );
		void solve_lin_sys( void );
		void solution_loop( void );
		void export_results( void );
		void export_params( void );

	private:
		c_mesh*      m_mesh_ptr;
		c_mesh_adap* m_mesh_adap_ptr;
		c_info m_info;
		c_lin_sys m_lin_sys;
};

} // namespace NBCFD

#endif // CLASS_SOLVER_H

