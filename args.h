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
// Last update:	20 Jun 2006
//-----------------------------------------------------------------------------

#ifndef CLASS_NIFS_C_ARGS_H
#define CLASS_NIFS_C_ARGS_H

#include <string>
#include <list>

namespace NIFS{
	
//--- INTERFACE ---	
	
class c_args {
  public:
    c_args( );
    virtual ~c_args( ) { }    void initialize( int , char** );
    
    std::string get_mesh_file_name( void ) 	{ return m_mesh_file_name; }
    std::string get_params_file_name( void ){ return m_param_file_name; }
    std::string get_bc_file_name( void )	{ return m_bc_file_name; }
    std::string get_log_file_name( void )	{ return m_log_file_name; }
    std::string get_cmat_file_name( void )  { return m_cmat_file_name; }
    std::string get_rhsv_file_name( void )  { return m_rhsv_file_name; }
    std::string get_init_file_name( void )  { return m_init_file_name; }
    std::string get_soln_file_name( void )  { return m_soln_file_name; }
    std::string get_grad_file_name( void )  { return m_grad_file_name; }
    std::string get_hess_file_name( void )  { return m_hess_file_name; }
    std::string get_solb_file_name( void )  { return m_solb_file_name; }
    std::string get_rsdl_file_name( void )  { return m_rsdl_file_name; }
    std::string get_omsh_file_name( void )  { return m_omsh_file_name; }
    std::string get_geom_file_name( void )  { return m_geom_file_name; }    
    
  private:
    std::string m_mesh_file_name;	// mesh file
    std::string m_param_file_name;	// parameters file
    std::string m_bc_file_name;		// boundary condition file
    std::string m_log_file_name;	// log file
    std::string m_cmat_file_name;   // coefficient matrix file
    std::string m_rhsv_file_name;   // right hand side vector file
    std::string m_init_file_name;   // initialization data file name
    std::string m_soln_file_name;   // solution file name
    std::string m_grad_file_name;   // gradient file name
    std::string m_hess_file_name;   // hessian initialization data file name		
    std::string m_solb_file_name;   // solution file name in binary format
    std::string m_geom_file_name;   // gemotry 
    std::string m_rsdl_file_name;   // residual
    std::string m_omsh_file_name;   // out mesh file name
	
    std::string m_flags;			  // flags
    std::list< std::string > m_arg_list;
};

}	//--- namespace NIFS

#endif // CLASS_NIFS_C_ARGS_H
