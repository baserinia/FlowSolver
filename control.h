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
//----------------------------------------------------------------------------
// Created on:	15 Mar 2005
// Last update:	30 Mar 2005
//-----------------------------------------------------------------------------

#ifndef CLASS_NIFS_C_CONTROL_H
#define CLASS_NIFS_C_CONTROL_H

#include "info.h"
#include "args.h"
#include "mesh.h"
#include "solver.h"

namespace NIFS{

//--- INTERFACE ---

class c_control {
	public:
		c_control( ) { }
		virtual ~c_control( ) { }
		void initialize( int, char** );
		void solve( void );
	private:
		c_info	m_info;
		c_mesh	m_mesh;
		c_solver m_solver;
};

}

#endif //--- CLASS_NIFS_C_CONTROL_H ---
