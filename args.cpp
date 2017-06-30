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

//----------------------------------------------------------------------------
// Created on:  15 Mar 2005
// Last update: 15 Jun 2006
//----------------------------------------------------------------------------
    
#include <cstring>	
#include "args.h"
#include "my_except.h"

//--- IMPLEMENTATION ---


//---( c_args )---
// This is the default constructor
NIFS::c_args::c_args( void )
{
  m_arg_list.clear();
  m_flags = "";
}


//---( c_args::initialize )---
//	This method is 
void 
NIFS::c_args::initialize( int argc , char* argv[] )
{
  std::string base;

  if ( argc == 1 ) 
    // if user dosn't enter the project name, "a" is used by default
    base = "a";
  else if ( argc == 2 ) 
    // if user enters the project name, it is used.
    base = argv[1];
  else
    // if more than one argument is passed to the code, an exception is thrown.
    NIFS::my_except( "too many arguments" );

  // file names  = project name + the specific extensions	
	
  m_mesh_file_name  = base + ".msh";   // mesh
  m_param_file_name = base + ".par";   // parameters
  m_bc_file_name    = base + ".bc";    // boundary conditions
  m_log_file_name   = base + ".log";   // log file
  m_cmat_file_name  = base + ".cmat";  // coefficient matrix
  m_rhsv_file_name  = base + ".rhsv";  // right hand side vector
  m_init_file_name  = base + ".init";  // initializizer
  m_soln_file_name  = base + ".soln";  // solution in ascii format
  m_grad_file_name  = base + ".grad";  // solution gradient
  m_hess_file_name  = base + ".hess";  // solution hessian
  m_solb_file_name  = base + ".solb";  // solution in binary format
  m_geom_file_name  = base + ".geom";  // geometry informaion output
  m_rsdl_file_name  = base + ".rsdl";  // residual output
  m_omsh_file_name  = base + ".omsh";  // out mesh
}
