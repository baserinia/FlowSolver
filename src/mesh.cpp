
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

//----------------------------------------------------------------------------
// Created on:  09 Apr 2004
// Last update:	22 Mar 2006
//-----------------------------------------------------------------------------

#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "mesh.h" 
#include "face_bnd.h"
#include "face_int.h"
#include "face_bnd_wall.h"
#include "face_bnd_wall_dir.h"
#include "face_bnd_wall_neu.h"
#include "face_bnd_wall_dir_pres.h"
#include "face_bnd_inflow.h"
#include "face_bnd_outflow.h"
#include "face_bnd_sym.h"
#include "my_except.h"
#include "log_file.h"
#include <iostream>
#include "matrix_2d.h"

//--- IMPLEMENTATION ---

NIFS::c_mesh::c_mesh( c_mesh& mesh )
{
	(*this) = mesh;
}

void
NIFS::c_mesh::release_verts( void )
{
    for ( m_vert_iter = m_verts.begin(); 
		  m_vert_iter != m_verts.end(); m_vert_iter++)
        delete (*m_vert_iter);
	m_verts.clear();
}

void
NIFS::c_mesh::release_faces( void )
{
    for ( m_face_iter = m_faces.begin(); 
		  m_face_iter != m_faces.end(); m_face_iter++)
        delete (*m_face_iter);
	m_faces.clear();
}

void
NIFS::c_mesh::release_cells( void )
{
    for ( m_cell_iter = m_cells.begin(); 
		  m_cell_iter != m_cells.end(); m_cell_iter++)
        delete (*m_cell_iter);
	m_cells.clear();
}

void
NIFS::c_mesh::release_mesh( void )
{
	release_verts();
	release_faces();
	release_cells();

}

NIFS::c_cell* 
NIFS::c_mesh::get_first_cell( void )
{
	m_cell_iter = m_cells.begin();
	return (*m_cell_iter);
}

NIFS::c_cell* 
NIFS::c_mesh::get_next_cell( void )

{
	if ( (++m_cell_iter) == m_cells.end() )
		return NULL;
	else
		return (*m_cell_iter);
}

void 
NIFS::c_mesh::reset_vert_index( void )
{
	int i = 0;
	for ( m_vert_iter  = m_verts.begin(); 
       	  m_vert_iter != m_verts.end(); ++m_vert_iter ) 
		(*m_vert_iter)->set_index( i++ );
}


void 
NIFS::c_mesh::reset_face_index( void )
{
	int i = 0;
	for ( m_face_iter = m_faces.begin(); 
       	  m_face_iter != m_faces.end(); m_face_iter++) 
		(*m_face_iter)->set_index( i++ );
}


void 
NIFS::c_mesh::reset_cell_index( void )
{
	int i = 0;
	for ( m_cell_iter = m_cells.begin(); 
       	  m_cell_iter != m_cells.end(); m_cell_iter++) 
		(*m_cell_iter)->set_index( i++ );
}


int 
NIFS::c_mesh::load_mesh( std::string file_name )
{
	unsigned int vert_num;		// number of vertices
	unsigned int face_bnd_num;	// number of boundary faces
	unsigned int face_int_num;	// number of internal faces
	unsigned int cell_num;		// number of cells
	unsigned int side_num;		// number of sides of a cell
	unsigned int i,j;
	unsigned int vert_1, vert_2;// Face vertex 1 and 2 index
	unsigned int bnd_code;
	double x;	// vertex x position
	double y;	// vertex y position
	c_vertex*	vert_ptr;
	c_face*		face_ptr;	
	c_cell*		cell_ptr;	

	// open mesh file
	std::ifstream ifs( file_name.c_str(), std::ios::in );	
	if ( !ifs ) {
		throw NIFS::my_except(	"mesh file \"" + 
								file_name + 
								"\" cannot be openned." );
	}

	m_info.get_log_file()->print_line( "mesh file loading started" );
	m_info.get_log_file()->print_elapsed_time();

	ifs >> vert_num;
	ifs	>> face_bnd_num;
	ifs	>> face_int_num;
	ifs	>> cell_num;

	// load vertices
	std::vector< c_vertex* > temp_vert( vert_num );
	for ( i = 0 ; i < vert_num ; ++i ) {
		ifs >> x >> y;
		vert_ptr = new c_vertex( x , y );
		m_verts.push_back( vert_ptr );	//  <-----
		temp_vert[i] = vert_ptr;
		vert_ptr->set_index( i );
	}

	// load boundary faces
	for ( i = 0 ; i < face_bnd_num ; ++i ) {
	    ifs >> vert_1 >> vert_2 >> bnd_code;
	    if ( bnd_code > m_info.get_bnd_cond()->get_bnd_num() ) {
			throw NIFS::my_except(	"invalid boundary number in \"" +
									file_name +
									"\" found" );
    	}
    	// dynamic construction
	    switch ( m_info.get_bnd_cond()->get_bc_type( bnd_code ) ) {     
		case 1:
			face_ptr = new c_face_bnd_inflow();		// in-flow
			break;
		case 2:
			face_ptr = new c_face_bnd_outflow();	// out-flow
			break;
		case 3:
			face_ptr = new c_face_bnd_sym();		// symmetry
			break;
		case 4:
			face_ptr = new c_face_bnd_wall_dir( );	// wall dirichlet
			break;
		case 5:
			face_ptr = new c_face_bnd_wall_neu( );	// wall neumann
			break;
		case 6:
			face_ptr = new c_face_bnd_wall_dir_pres();	// wall dirichlet & pres
			break;
        default:
			throw NIFS::my_except( "unidentified boundary condition type" );
		} // switch
		face_ptr->set_code( bnd_code );				// type of boundary
		face_ptr->push_vert( temp_vert[ vert_1 ] ); // introduce V1 to face
		face_ptr->push_vert( temp_vert[ vert_2 ] ); // introduce V2 to face
		temp_vert[ vert_1 ]->push_face( face_ptr ); // introduce face to V1
		temp_vert[ vert_2 ]->push_face( face_ptr ); // introduce face to V2
		temp_vert[ vert_1 ]->push_vert( temp_vert[ vert_2 ] );	// introduce V2 to V1
		temp_vert[ vert_2 ]->push_vert( temp_vert[ vert_1 ] );	// introduce V1 to V2
        m_faces.push_back( face_ptr );        
        face_ptr->set_index( i );
    } // for boundary faces

	// load interior faces
	for ( i = 0 ; i < face_int_num ; ++i ) {
	    ifs >> vert_1 >> vert_2;
        face_ptr = new c_face_int();
        face_ptr->push_vert( temp_vert[ vert_1 ] ); // introduce V1 to face
        face_ptr->push_vert( temp_vert[ vert_2 ] ); // introduce V2 to face
        temp_vert[ vert_1 ]->push_face( face_ptr ); // introduce face to V1
		temp_vert[ vert_2 ]->push_face( face_ptr ); // introduce face to V2
		temp_vert[ vert_1 ]->push_vert( temp_vert[ vert_2 ] ); // introduce V2 to V1
		temp_vert[ vert_2 ]->push_vert( temp_vert[ vert_1 ] ); // introduce V1 to V2
        m_faces.push_back( face_ptr );        
        face_ptr->set_index( i + face_bnd_num );
    } // for interior faces

	// load cells
	for ( i=0 ; i < cell_num ; ++i)
	{	
		cell_ptr = new c_cell();	// construct a cell
		ifs >> side_num;					// load number of faces
		for ( j = 0 ; j < side_num ; ++j ) {
			ifs >> vert_1;
			cell_ptr->push_vert( temp_vert[ vert_1 ] ); // introduce V[i] to cell
			temp_vert[ vert_1 ]->push_cell( cell_ptr );	// introduce cell to V[i]
		}
		cell_ptr->set_index( i );		// cell auxilary index
		m_cells.push_back( cell_ptr );	// push back to list
	}
	
	temp_vert.clear();	// destroy temporary vertex vector
	ifs.close();		// close input mesh file

	// extract face-cell neighborhood
	std::list< c_cell* >::iterator cell_iter;
	for ( cell_iter = m_cells.begin() ; cell_iter != m_cells.end() ; ++cell_iter ) {
		for ( i = 0 ; i < (*cell_iter)->get_vert_num() ; ++i ) {
			vert_ptr = static_cast<c_vertex*>( (*cell_iter)->get_vert( i ) );
			for ( j = 0 ; j < vert_ptr->get_face_num() ; ++j ) {
				face_ptr = static_cast<c_face*>( vert_ptr->get_face( j ) );
				if ( ( face_ptr->get_vert( 0 ) == (*cell_iter)->get_vert( i+1 ) ) ||
					 ( face_ptr->get_vert( 1 ) == (*cell_iter)->get_vert( i+1 ) ) ) {
					(*cell_iter)->push_face( face_ptr ); // introduce face to cell
					face_ptr->push_cell( *cell_iter ); // introduce cell to face
				}
			} 
		}
	}

	// extract cells neighborhood
	std::list<c_face*>::iterator face_iter;
	for ( face_iter = m_faces.begin() ; face_iter != m_faces.end(); ++face_iter ) {
		if ( (*face_iter)->get_cell_num() == 2 ) {
			static_cast<c_cell*>( (*face_iter)->get_cell(0) )->push_cell( (*face_iter)->get_cell(1) );
			static_cast<c_cell*>( (*face_iter)->get_cell(1) )->push_cell( (*face_iter)->get_cell(0) );
		}
	}
	
	

	m_info.get_log_file()->print_line( 
		"mesh file loaded and connectivities extracted." );
	m_info.get_log_file()->print_elapsed_time();
	
	load_size();
	
	return 0;
}


//---( load initial solution )--------------------------------------------------
void
NIFS::c_mesh::load_initial_sol( std::string file_name )
{
  if ( m_info.get_params()->get_initialize() == 0.0 )
	return; // the initial solution is set to zero

  double x;
  double y;
  vector_block< double > state;
  
  // open initialization dara file
  std::ifstream ifs( file_name.c_str(), std::ios::in );	
  if ( !ifs ) {
	throw NIFS::my_except(	"Solution initialization file \"" + 
							file_name + 
							"\" cannot be openned." );
  }
  
  for ( m_cell_iter = m_cells.begin() ; m_cell_iter != m_cells.end() ; ++m_cell_iter ) {
	ifs >> x >> y >> state[0] >> state[1] >> state[2] >> state[3];
	if ((x!=0.0)&&(x!=16.0)&&(y!=0.0)&&(y!=1.0)) {
	  (*m_cell_iter)->set_force_state( state );
	}
  }
  ifs.close();
}



//---( Update Solution )--------------------------------------------------------
void 
NIFS::c_mesh::update_solution( void )
{

	int gh_iter = static_cast< int >( m_info.get_params()->get_gh_iter() );
	for ( int i = 0 ; i < gh_iter ; ++i ) {
		for ( m_cell_iter = m_cells.begin() ; m_cell_iter != m_cells.end() ; ++m_cell_iter ) 
			(*m_cell_iter)->calc_gradient();
			
		for ( m_cell_iter = m_cells.begin() ; m_cell_iter != m_cells.end() ; ++m_cell_iter ) 
			(*m_cell_iter)->calc_hessian();
	}

	for ( m_face_iter = m_faces.begin() ; m_face_iter != m_faces.end() ; ++m_face_iter ) 
	{
		(*m_face_iter)->update_solution();
/*		std::cout << (*m_face_iter)->get_mfr() << "\t"
		          << (*m_face_iter)->get_fmp().get_pos().get_x() << "\t"
		          << (*m_face_iter)->get_fmp().get_pos().get_y() << std::endl; */
  }
//  std::cout << std::endl;

	m_avg_err = 0.0;
	for ( m_cell_iter = m_cells.begin() ; m_cell_iter != m_cells.end() ; ++m_cell_iter ) 
	{
		(*m_cell_iter)->calc_coefs();
		for ( int i = 0 ; i < 4 ; ++i )
			m_avg_err[i] += 
				std::pow( std::abs( (*m_cell_iter)->get_err_iso_coef()( i ) ) , 0.4 ) *
				(*m_cell_iter)->get_volume();
	}
	m_avg_err /= static_cast<double>( get_cell_num() );
	for ( int i = 0 ; i < 4 ; ++i )
		m_avg_err[i] = std::sqrt( m_avg_err( i ) );
}

//---( Update Geometry )--------------------------------------------------------
void 
NIFS::c_mesh::update_geometry( void )
{
	reset_vert_index();
	reset_face_index();
	reset_cell_index();

	// updates the physical parameters within the cells
	m_coef_num = 0;
	for ( m_cell_iter = m_cells.begin() ; m_cell_iter != m_cells.end() ; ++m_cell_iter ) {
		(*m_cell_iter)->update_geometry();
		m_coef_num += (*m_cell_iter)->get_cell_num() + 1;
	}
	
	// updates the physical parameters at the faces
	for ( m_face_iter = m_faces.begin() ; m_face_iter != m_faces.end() ; ++m_face_iter ) {
		(*m_face_iter)->update_geometry();
	}
}

//---( Residual RMS )-----------------------------------------------------------
NIFS::vector_block< double >
NIFS::c_mesh::get_residual_rms( void )
{
	vector_block< double > res( 0.0 );
	
	for ( m_cell_iter = m_cells.begin() ; m_cell_iter != m_cells.end() ; ++m_cell_iter ) {
		for ( int i = 0 ; i < 4 ; ++i )
			res[i] += sqr( (*m_cell_iter)->get_residual()(i) );
	}

	for ( int i = 0 ; i < 4 ; ++i )
		res[i] = std::sqrt( res(i) / get_cell_num() );
	
	res[0] /= m_info.get_params()->get_pres_scale();
	res[1] /= m_info.get_params()->get_vel_scale();
	res[2] /= m_info.get_params()->get_vel_scale();
	res[3] /= m_info.get_params()->get_temp_scale();

	return res;
}


NIFS::vector_block< double >
NIFS::c_mesh::get_residual_max( void )
{
	vector_block< double > res( 0.0 );
	double dummy;
	
	for ( m_cell_iter = m_cells.begin() ; m_cell_iter != m_cells.end() ; ++m_cell_iter ) {
		for ( int i = 0 ; i < 4 ; ++i )
			if ( res(i) < ( dummy = (*m_cell_iter)->get_residual()(i) ) )
				res[i] = dummy;
	}
	
	res[0] /= m_info.get_params()->get_pres_scale();
	res[1] /= m_info.get_params()->get_vel_scale();
	res[2] /= m_info.get_params()->get_vel_scale();
	res[3] /= m_info.get_params()->get_temp_scale();

	return res;
}


void 
NIFS::c_mesh::print_info( void )
{
	
	
}

int 
NIFS::c_mesh::save_mesh( )
{
  int n;

  std::ofstream ofs( m_info.get_args()->get_omsh_file_name( ).c_str() , 
                     std::ios::out );	
  ofs.setf( std::ios::scientific | 
            std::ios::showpoint  );
  ofs.precision( 5 ); 
  ofs.width( 12 ); 

  for ( m_vert_iter = m_verts.begin(); 
        m_vert_iter != m_verts.end(); 
        m_vert_iter++) {
    ofs << std::setw(12) 
        << static_cast<c_vertex*>( *m_vert_iter )->get_point().get_pos().get_x() 
        << "  "
        << std::setw(12) 
        << static_cast<c_vertex*>( *m_vert_iter )->get_point().get_pos().get_y()
        << std::endl;
  } // for

  reset_vert_index(); // re-index vertices
  for ( m_face_iter = m_faces.begin(); m_face_iter != m_faces.end(); m_face_iter++) {
    if ( static_cast<c_face*>( *m_face_iter )->get_type() == 1 ) // boundary
    	ofs << std::setw(12) 
    		<< ( static_cast<c_face*>( *m_face_iter )->get_vert( 0 ) )->get_index()
    		<< "  "
    		<< std::setw(12) 
    		<< ( static_cast<c_face*>( *m_face_iter )->get_vert( 1 ) )->get_index()
    		<< "  "
    		<< std::setw(12) 
    		<< ( static_cast<c_face*>( *m_face_iter ) )->get_code() 
    		<< std::endl;
    }

    for ( m_face_iter = m_faces.begin(); m_face_iter != m_faces.end(); m_face_iter++) {
    	if ( static_cast<c_face*>( *m_face_iter )->get_type() == 0 ) // boundary
    	ofs << std::setw(12) 
    		<< ( static_cast<c_face*>( *m_face_iter )->get_vert( 0 ) )->get_index()
    		<< "  "
    		<< std::setw(12) 
    		<< ( static_cast<c_face*>( *m_face_iter )->get_vert( 1 ) )->get_index()
    		<< std::endl;
    }

    for ( m_cell_iter = m_cells.begin(); m_cell_iter != m_cells.end(); m_cell_iter++) {
      n = static_cast<c_cell*>( *m_cell_iter )->get_vert_num();
    	ofs << n;
    	for (int i =0; i<n; i++){
    		ofs << std::setw(12) 
    		<< ( static_cast<c_cell*>( *m_cell_iter )->get_vert( i ) )->get_index()
    		<< "  ";
    	}
	   	ofs	<< std::endl;
    }
	
	return 0;
}



void
NIFS::c_mesh::load_size( void )
{
  double x,y;
	std::ifstream ifs( "cav_len.dat", std::ios::in );		
	for ( int i = 0 ; i < 256 ; ++i )
		for ( int j = 0 ; j < 256 ; ++j ) {
			ifs >> x >> y >> A11[i][j] >> A12[i][j] >> A22[i][j];
	  }
}

double 
NIFS::c_mesh::get_size( c_vector_2d vec )
{
//	int x = ( vec.get_x() == 1.0 ) ? 199 : std::floor( vec.get_x() * 200.0 );
//	int y = ( vec.get_y() == 1.0 ) ? 199 : std::floor( vec.get_y() * 200.0 );
	
//	y = 199 - y;
	
//	return 0.049 * double(spsize[y][x]) / 255.0 + 0.001;
/*
  if ( pow(vec.get_x()-0.5,2.0) + pow(vec.get_y()-0.5,2.0) < 0.01 )
    return 0.01;
  else if ( pow(vec.get_x()-0.5,2.0) + pow(vec.get_y()-0.5,2.0) > 0.04 )
  	return 0.1;
  else
    return 0.01 + (sqrt(pow(vec.get_x()-0.5,2.0) + pow(vec.get_y()-0.5,2.0))-0.1)*0.9;
*/
  return 0.1-vec.get_x()*0.05;
}

NIFS::c_matrix_2d
NIFS::c_mesh::get_size_matrix( c_vector_2d vec )
{

	int i = static_cast<int>( ( vec.get_x() >= 1.0 ) ? 255.0 : floor( vec.get_x() * 256.0 ) );
	int j = static_cast<int>( ( vec.get_y() >= 1.0 ) ? 255.0 : floor( vec.get_y() * 256.0 ) );
	double d = 40.0;
	double r,t,l1,l2;
  
  
  return c_matrix_2d( d*A11[i][j],
                      d*A12[i][j],
                      d*A12[i][j],
                      d*A22[i][j] );

/*
  return  c_matrix_2d( 55,
                       -45,
                       -45,
                       55 );
*/
                       
/*
  r = sqrt( pow( vec.get_x() - 0.5 , 2.0 ) + pow( vec.get_y() - 0.5 , 2.0 ) );
  t = atan2( vec.get_y() - 0.5 , vec.get_x() - 0.5 );
  l1 = 20.0+80.0*exp(-80*pow(r-0.25,2.0));
  l2 = 20.0;
  return c_matrix_2d( l1*pow(cos(t),2.0)+l2*pow(sin(t),2.0),
                      (l1-l2)*cos(t)*sin(t),
                      (l1-l2)*cos(t)*sin(t),
                      l1*pow(sin(t),2.0)+l2*pow(cos(t),2.0));
*/

//	reset_current_cell();
/*
	c_cell* cell_ptr = search( vec );
	double d = 1.0E10;
	double temp;
	if ( cell_ptr != NULL ){
		for ( int i = 0 ; i < 4 ; ++i ) {	
			temp = get_avg_err()( i ) / std::pow(std::abs(cell_ptr->get_err_iso_coef()(i)),0.2);
			if ( temp < d )
				d = temp; 
		}
		if ( d > 2.0 * cell_ptr->get_charac_size() )
			d = cell_ptr->get_charac_size() * 2.0;
		if ( d < 0.1 * cell_ptr->get_charac_size() )
			d = cell_ptr->get_charac_size() * 0.1;
	}
	
  return c_matrix_2d( 1.0 / d,0.0,0.0,1.0 / d  );
*/
}



// returns the desirable mesh size
/*
double 
NIFS::c_mesh::get_size( c_vector_2d vec )
{
	c_cell* cell_ptr = search( vec );
	double d = 0.0;
	double temp;
	if ( cell_ptr != NULL ){
		for ( int i = 0 ; i < 4 ; ++i ) {	
			temp = get_avg_err()( i ) / std::pow(std::abs(cell_ptr->get_err_iso_coef()(i)),0.2);
			if ( temp > d )
				d = temp; 
		}
		if ( d > 2.0 * cell_ptr->get_charac_size() )
			d = cell_ptr->get_charac_size() * 2.0;
		if ( d < 0.1 * cell_ptr->get_charac_size() )
			d = cell_ptr->get_charac_size() * 0.1;
	}

	if ( d == 0.0 )
		return 0.05;
	else
		return d;
}
*/

double 
NIFS::c_mesh::get_length_ratio_move( c_vector_2d vec1 , c_vector_2d vec2 )
{
  double t;
  double r = 0.0;
  unsigned n = 50;
  c_vector_2d v1;
  c_vector_2d v2;
  c_vector_2d vm;
  c_vector_2d vd;
  c_vector_2d vt;
  c_vector_2d t_hat = ( vec2 - vec1 ).get_unit();
  
  for ( unsigned i = 0 ; i < n ; ++i ) {
    t = static_cast<double>( i ) / static_cast<double>( n );
    v1 = vec1 + ( vec2 - vec1 ) * t ;
    t = static_cast<double>( i + 1 ) / static_cast<double>( n );
    v2 = vec1 + ( vec2 - vec1 ) * t;
    vm = ( v1 + v2 ) / 2.0;
    vd = v2 - v1;
    vt = get_size_matrix( vm ) * vd;
    r += pow( pow(vt(0),100.0)+pow(vt(1),100.0), 0.01);
  }
  return r;
}


double 
NIFS::c_mesh::get_length_ratio( c_vector_2d vec1 , c_vector_2d vec2 )
{
  double t;
  double r = 0.0;
  double theta;
  unsigned n = 50;
  c_vector_2d v1;
  c_vector_2d v2;
  c_vector_2d vm;
  c_vector_2d vd;
  c_vector_2d vt;
  c_vector_2d t_hat = ( vec2 - vec1 ).get_unit();
  c_matrix_2d eig, trans;
  c_vector_2d ev_1;
  c_vector_2d ev_2;
  double eig_1, eig_2, cndno;
  
  
  for ( unsigned i = 0 ; i < n ; ++i ) {
    t = static_cast<double>( i ) / static_cast<double>( n );
    v1 = vec1 + ( vec2 - vec1 ) * t ;
    t = static_cast<double>( i + 1 ) / static_cast<double>( n );
    v2 = vec1 + ( vec2 - vec1 ) * t;
    vm = ( v1 + v2 ) / 2.0;
    vd = v2 - v1;
    trans = get_size_matrix( vm );
    vt = trans * vd;
    eig = trans.eig_vec();
    ev_1 = eig(0).get_unit();
    ev_2 = eig(1).get_unit();
//    std::cout << ev_1(0) << ' ' << ev_1(1) << ' ' << ev_2(0) << ' ' << ev_2(1) << std::endl;  
    theta = atan2( vt * ev_2 , vt * ev_1 );  
    
    eig = trans.eig_val();
    eig_1 = std::abs( eig(0)(0) );
    eig_2 = std::abs( eig(1)(1) );
//    std::cout << eig_1 << ' ' << eig_2 << std::endl;  
    cndno = ( eig_1 < eig_2 ) ? (eig_2/eig_1) : ( eig_1/eig_2);
    if ( cndno < 4.0 )
      r += sqrt(vt(0)*vt(0)+vt(1)*vt(1));
    else
      r += sqrt(vt(0)*vt(0)+vt(1)*vt(1)) 
      / ( fabs(cos( theta )) + fabs(sin( theta )) );
  }
  return r;
}


void 
NIFS::c_mesh::operator=( c_mesh& mesh )
{
	unsigned i , j;
	c_vertex* vert_ptr = NULL;
	c_face*   face_ptr = NULL;
	c_cell*   cell_ptr = NULL;
	

	i = 0;
	std::vector< c_vertex* > temp_vert( mesh.get_vert_num() );	
	for ( m_vert_iter  = mesh.m_verts.begin(); 
		  m_vert_iter != mesh.m_verts.end(); 
		  ++m_vert_iter ) {
		vert_ptr = new c_vertex( (*m_vert_iter)->get_point().get_pos().get_x(),
								 (*m_vert_iter)->get_point().get_pos().get_y() );
		m_verts.push_back( vert_ptr );	
		temp_vert[i] = vert_ptr;
		vert_ptr->set_index( i++ );
	}

	i = 0;
	std::vector< c_face* > temp_face( mesh.get_face_num() );	
    for ( m_face_iter = mesh.m_faces.begin(); 
    	  m_face_iter != mesh.m_faces.end(); 
    	  ++m_face_iter )
    {
    	if ( (*m_face_iter)->get_type() == 0 ) {  // internal faces
    		face_ptr = new c_face_int();
    	}
    	else {   // boundary faces
		    switch ( m_info.get_bnd_cond()->get_bc_type( (*m_face_iter)->get_code() ) ) {     
			case 1:
				face_ptr = new c_face_bnd_inflow();		// in-flow
				break;
			case 2:
				face_ptr = new c_face_bnd_outflow();	// out-flow
				break;
			case 3:
				face_ptr = new c_face_bnd_sym();		// symmetry
				break;
			case 4:
				face_ptr = new c_face_bnd_wall_dir( );	// wall dirichlet
				break;
			case 5:
				face_ptr = new c_face_bnd_wall_neu( );	// wall neumann
				break;
			case 6:
				face_ptr = new c_face_bnd_wall_dir_pres();	// wall dirichlet & pres
			} // switch
			face_ptr->set_code( (*m_face_iter)->get_code() );	
	   	} // if
		m_faces.push_back( face_ptr );        
		temp_face[i] = face_ptr;            
		face_ptr->set_index( i++ );
	} // for loop
	
	// allocate memory to the new cells
	i = 0;
	std::vector< c_cell* > temp_cell( mesh.get_cell_num() );	
	for ( m_cell_iter = mesh.m_cells.begin(); 
		  m_cell_iter != mesh.m_cells.end(); 
		  ++m_cell_iter ) {
		cell_ptr = new c_cell();	
		m_cells.push_back( cell_ptr );
		temp_cell[i] = cell_ptr;
		cell_ptr->set_index( i++ );
    }

	mesh.reset_vert_index();
	mesh.reset_face_index();
	mesh.reset_cell_index();

	// construct local topology for the vertices
	i = 0;
	for ( m_vert_iter  = mesh.m_verts.begin(); 
		  m_vert_iter != mesh.m_verts.end(); 
		  ++m_vert_iter ) {
		  
		for ( j = 0 ; j < (*m_vert_iter)->get_vert_num() ; ++j ) {
			temp_vert[i]->push_vert( temp_vert[ (*m_vert_iter)->get_vert( j )->get_index() ] );
		}
		
		for ( j = 0 ; j < (*m_vert_iter)->get_face_num() ; ++j ) {
			temp_vert[i]->push_face( temp_face[ (*m_vert_iter)->get_face( j )->get_index() ] );
		}
		
		for ( j = 0 ; j < (*m_vert_iter)->get_cell_num() ; ++j ) {
			temp_vert[i]->push_cell( temp_cell[ (*m_vert_iter)->get_cell( j )->get_index() ] );
		}
		
		++i;
	}


	// construct local topology for the faces
	i = 0;
	for ( m_face_iter  = mesh.m_faces.begin(); 
		  m_face_iter != mesh.m_faces.end(); 
		  ++m_face_iter ) {
		  
		for ( j = 0 ; j < (*m_face_iter)->get_vert_num() ; ++j ) 
			temp_face[i]->push_vert( temp_vert[ (*m_face_iter)->get_vert( j )->get_index() ] );
			
		for ( j = 0 ; j < (*m_face_iter)->get_face_num() ; ++j ) 
			temp_face[i]->push_face( temp_face[ (*m_face_iter)->get_face( j )->get_index() ] );
			
		for ( j = 0 ; j < (*m_face_iter)->get_cell_num() ; ++j ) 
			temp_face[i]->push_cell( temp_cell[ (*m_face_iter)->get_cell( j )->get_index() ] );

		++i;
	}

	// construct local topology for the vertices
	i = 0;
	for ( m_cell_iter  = mesh.m_cells.begin(); 
		  m_cell_iter != mesh.m_cells.end(); 
		  ++m_cell_iter ) {
		  
		for ( j = 0 ; j < (*m_cell_iter)->get_vert_num() ; ++j ) 
			temp_cell[i]->push_vert( temp_vert[ (*m_cell_iter)->get_vert( j )->get_index() ] );
			
		for ( j = 0 ; j < (*m_cell_iter)->get_face_num() ; ++j ) 
			temp_cell[i]->push_face( temp_face[ (*m_cell_iter)->get_face( j )->get_index() ] );
			
		for ( j = 0 ; j < (*m_cell_iter)->get_cell_num() ; ++j ) 
			temp_cell[i]->push_cell( temp_cell[ (*m_cell_iter)->get_cell( j )->get_index() ] );
			

		++i;
	}
	temp_vert.clear();	
	temp_face.clear();
	temp_cell.clear();
}


//---( search )------------------------------------------
// This method looks for the cell which the vector "vec" is in.
// It uses walk search algorithm
//
NIFS::c_cell* 
NIFS::c_mesh::search( c_vector_2d vec )
{
  c_vector_2d unv;	
  c_vector_2d temp_vec;
  double max;
  c_face* face_ptr;
  c_face* face_go_thru;
  c_cell* cell_ptr = *m_cells.begin();
	
  do {
    face_go_thru = NULL;
    max = 1.0E-14;
    
    for ( unsigned i = 0 ; i < cell_ptr->get_face_num() ; ++i )
    {
      face_ptr = static_cast< c_face*>( cell_ptr->get_face( i ) );
      unv = face_ptr->get_unv();
      temp_vec = face_ptr->get_fmp().get_pos() - cell_ptr->get_node().get_pos();

      if ( unv * temp_vec < 0.0 )
        unv = -unv;

      temp_vec = ( vec - face_ptr->get_fmp().get_pos() ).get_unit();

      if ( unv * temp_vec > max ) {
        face_go_thru = face_ptr;
        max = unv * temp_vec;
      }
    }
    
    if ( ( face_go_thru != NULL ) && ( face_go_thru->get_type() == 0 ) ) {
      if ( cell_ptr == face_go_thru->get_cell( 0 ) ) 
        cell_ptr = static_cast< c_cell*>( face_go_thru->get_cell( 1 ) );
      else 
        cell_ptr = static_cast< c_cell*>( face_go_thru->get_cell( 0 ) );
    }
  } while ( ( face_go_thru != NULL ) &&
            ( face_go_thru->get_type() == 0 ) );

  return cell_ptr;
}


