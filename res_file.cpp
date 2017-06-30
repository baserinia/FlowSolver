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
//	Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1401  USA

//----------------------------------------------------------------------------
// Created on:	20 Apr 2005
// Last update:	26 Jun 2006
//-----------------------------------------------------------------------------

#include <iomanip>
#include "res_file.h"
#include "my_except.h"


void 
NIFS::c_res_file::open_soln_file( std::string soln_file_name )
{
  m_ofs_soln.open( soln_file_name.c_str(), std::ios::out );	
  if ( !m_ofs_soln ) {
	throw NIFS::my_except(	"results file \"" + 
							soln_file_name + 
							"\" cannot be created." );
  }
}


void 
NIFS::c_res_file::open_grad_file( std::string grad_file_name )
{
  m_ofs_grad.open( grad_file_name.c_str(), std::ios::out );	
  if ( !m_ofs_grad ) {
	throw NIFS::my_except(	"results file \"" + 
							grad_file_name + 
							"\" cannot be created." );
  }
}


void 
NIFS::c_res_file::open_hess_file( std::string hess_file_name )
{
  m_ofs_hess.open( hess_file_name.c_str(), std::ios::out );	
  if ( !m_ofs_hess ) {
	throw NIFS::my_except(	"results file \"" + 
							hess_file_name + 
							"\" cannot be created." );
  }
}


void 
NIFS::c_res_file::open_rsdl_file( std::string rsdl_file_name )
{
  m_ofs_rsdl.open( rsdl_file_name.c_str(), std::ios::out );	
  if ( !m_ofs_rsdl ) {
	throw NIFS::my_except(	"results file \"" + 
							rsdl_file_name + 
							"\" cannot be created." );
  }
}


void 
NIFS::c_res_file::open_geom_file( std::string geom_file_name )
{
  m_ofs_geom.open( geom_file_name.c_str(), std::ios::out );	
  if ( !m_ofs_geom ) {
	throw NIFS::my_except(	"results file \"" + 
							geom_file_name + 
							"\" cannot be created." );
  }
}


void 
NIFS::c_res_file::open_solb_file( std::string solb_file_name )
{
  m_ofs_solb.open( solb_file_name.c_str(), std::ios::out | std::ios::binary );	
  if ( !m_ofs_solb ) {
	throw NIFS::my_except(	"results file \"" + 
							solb_file_name + 
							"\" cannot be created." );
  }
}


void 
NIFS::c_res_file::close_all()
{
  m_ofs_soln.close();
  m_ofs_grad.close();
  m_ofs_hess.close();
  m_ofs_solb.close();
  m_ofs_geom.close();
  m_ofs_rsdl.close();	
}


void
NIFS::c_res_file::print_line_soln( const c_vector_2d& pos,
								   const vector_block< double >& state )
{
  m_ofs_soln.setf( std::ios::scientific | std::ios::showpoint  );
  m_ofs_soln.precision( 7 ); // print up to 7 floating point digits
  m_ofs_soln.width( 14 ); 
  m_ofs_soln << std::setw(14) << pos.get_x() << "  ";
  m_ofs_soln << std::setw(14) << pos.get_y() << "  ";
  m_ofs_soln << std::setw(14) << state(0) << "  ";
  m_ofs_soln << std::setw(14) << state(1) << "  ";
  m_ofs_soln << std::setw(14) << state(2) << "  ";
  m_ofs_soln << std::setw(14) << state(3) << "  ";
  m_ofs_soln << std::endl;
}


void
NIFS::c_res_file::print_line_grad( const c_vector_2d& pos,
								   const vector_block< c_vector_2d >& grad )
{
  m_ofs_grad.setf( std::ios::scientific | std::ios::showpoint  );
  m_ofs_grad.precision( 7 ); 
  m_ofs_grad.width( 14 ); 
  m_ofs_grad << std::setw(14) << pos.get_x() << "  ";
  m_ofs_grad << std::setw(14) << pos.get_y() << "  ";
  m_ofs_grad << std::setw(14) << grad(0).get_x() << "  ";
  m_ofs_grad << std::setw(14) << grad(0).get_y() << "  ";
  m_ofs_grad << std::setw(14) << grad(1).get_x() << "  ";
  m_ofs_grad << std::setw(14) << grad(1).get_y() << "  ";
  m_ofs_grad << std::setw(14) << grad(2).get_x() << "  ";
  m_ofs_grad << std::setw(14) << grad(2).get_y() << "  ";
  m_ofs_grad << std::setw(14) << grad(3).get_x() << "  ";
  m_ofs_grad << std::setw(14) << grad(3).get_y() << "  ";
  m_ofs_grad << std::endl;
}



void 
NIFS::c_res_file::print_line_hess( const c_vector_2d& pos,
								   const vector_block< c_matrix_2d >& hess ) 
{
	m_ofs_hess.setf( std::ios::scientific | std::ios::showpoint  );
	m_ofs_hess.precision( 7 ); 
	m_ofs_hess.width( 14 ); 
	m_ofs_hess << std::setw(14) << pos.get_x() << "  ";
	m_ofs_hess << std::setw(14) << pos.get_y() << "  ";
	m_ofs_hess << std::setw(14) << hess(0)(0)(0) << "  ";
	m_ofs_hess << std::setw(14) << hess(0)(0)(1) << "  ";
	m_ofs_hess << std::setw(14) << hess(0)(1)(1) << "  ";
	m_ofs_hess << std::setw(14) << hess(1)(0)(0) << "  ";
	m_ofs_hess << std::setw(14) << hess(1)(0)(1) << "  ";
	m_ofs_hess << std::setw(14) << hess(1)(1)(1) << "  ";
	m_ofs_hess << std::setw(14) << hess(2)(0)(0) << "  ";
	m_ofs_hess << std::setw(14) << hess(2)(0)(1) << "  ";
	m_ofs_hess << std::setw(14) << hess(2)(1)(1) << "  ";
	m_ofs_hess << std::setw(14) << hess(3)(0)(0) << "  ";
	m_ofs_hess << std::setw(14) << hess(3)(0)(1) << "  ";
	m_ofs_hess << std::setw(14) << hess(3)(1)(1) << "  ";
	m_ofs_hess << std::endl;
}


void
NIFS::c_res_file::print_line_rsdl( const c_vector_2d& pos,
								   const vector_block< double >& rsdl )
{
  m_ofs_rsdl.setf( std::ios::scientific | std::ios::showpoint  );
  m_ofs_rsdl.precision( 7 ); // print up to 5 floating point digits
  m_ofs_rsdl.width( 14 ); 
  m_ofs_rsdl << std::setw(14) << pos.get_x() << "  ";
  m_ofs_rsdl << std::setw(14) << pos.get_y() << "  ";
  m_ofs_rsdl << std::setw(14) << rsdl(0) << "  ";
  m_ofs_rsdl << std::setw(14) << rsdl(1) << "  ";
  m_ofs_rsdl << std::setw(14) << rsdl(2) << "  ";
  m_ofs_rsdl << std::setw(14) << rsdl(3) << "  ";
  m_ofs_rsdl << std::endl;
}


void
NIFS::c_res_file::print_line_solb( const c_vector_2d& pos,
								   const vector_block< double >& state )
{
  m_ofs_solb << pos.get_x( );
  m_ofs_solb << pos.get_y( );
  m_ofs_solb << state(0);
  m_ofs_solb << state(1);
  m_ofs_solb << state(2);
  m_ofs_solb << state(3);
}


void
NIFS::c_res_file::print_line_geom( const c_vector_2d& pos,
								   double volume )
{
  m_ofs_geom.setf( std::ios::scientific | std::ios::showpoint  );
  m_ofs_geom.precision( 7 ); // print up to 5 floating point digits
  m_ofs_geom.width( 14 ); 
  m_ofs_geom << std::setw(14) << pos.get_x() << "  ";
  m_ofs_geom << std::setw(14) << pos.get_y() << "  ";
  m_ofs_geom << std::setw(14) << volume << "  ";
  m_ofs_geom << std::endl;
}
