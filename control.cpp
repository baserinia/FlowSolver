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
// Created on:	15 Mar 2005
// Last update:	30 Mar 2005
//------------------------------------------------------------------------------

#include <iostream>
#include <new>
#include "my_except.h"
#include "control.h"

//--- IMPLEMENTATION ---


//---( initialization )---
// initialize the global variables of the code
void 
NIFS::c_control::initialize( int argc, char* argv[] )
{
  try {
	m_info.get_args()->initialize( argc , argv );  
    m_info.get_log_file()->open( m_info.get_args()->get_log_file_name( ) );
    m_info.get_params()->load_params_file( m_info.get_args()->get_params_file_name( ) ); 
    m_info.get_bnd_cond()->load_bc_file( m_info.get_args()->get_bc_file_name( ) ); 
    m_solver.export_params();
    m_mesh.load_mesh( m_info.get_args()->get_mesh_file_name( ) );
  }
  catch ( const NIFS::my_except& e ) {
    throw;
  }
  catch ( const std::bad_alloc& error ) {
    std::cerr << "Initialization failed: " << error.what() << std::endl;
    throw;
  }
}


//---( solve )---
// solve the CFD problem
void 
NIFS::c_control::solve( void )
{
  m_solver.initialize( &m_mesh );
  m_solver.solution_loop();
  m_solver.export_results();
  m_info.release_all();
}

