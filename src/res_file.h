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

//------------------------------------------------------------------------------
// Created on:	20 Apr 2005
// Last update:	26 Jun 2006
//------------------------------------------------------------------------------

#ifndef CLASS_NIFS_C_RES_FILE_H
#define CLASS_NIFS_C_RES_FILE_H

#include <string>
#include "fstream"
#include "vector_block.h"
#include "vector_2d.h"
#include "matrix_2d.h"
#include "bfstream.h"

namespace NIFS{ 

//--- INTERFACE ---

class c_res_file
{
  public:
    c_res_file() {}
    virtual ~c_res_file() {}
    void open_soln_file( std::string soln_file_name );
    void open_grad_file( std::string grad_file_name );
    void open_hess_file( std::string hess_file_name );
    void open_solb_file( std::string solb_file_name );
    void open_rsdl_file( std::string rsdl_file_name );
    void open_geom_file( std::string geom_file_name );
    void close_all();
  
    void print_line_soln( const c_vector_2d&, const vector_block< double >& );
    void print_line_grad( const c_vector_2d&, const vector_block< c_vector_2d >& );
    void print_line_hess( const c_vector_2d&, const vector_block< c_matrix_2d >& );
    void print_line_solb( const c_vector_2d&, const vector_block< double >& );
    void print_line_rsdl( const c_vector_2d&, const vector_block< double >& );
    void print_line_geom( const c_vector_2d&, double );
    
  private:
    std::ofstream m_ofs_soln; // solution
    std::ofstream m_ofs_grad; // gradient
    std::ofstream m_ofs_hess; // hessian
    std::ofstream m_ofs_geom; // geometry
    std::ofstream m_ofs_rsdl; // residual
    bofstream     m_ofs_solb; // solution in binary format
};

} //--- namespace NIFS ---

#endif // CLASS_NIFS_C_RES_FILE_H
