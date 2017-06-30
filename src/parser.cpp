// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//----------------------------------------------------------------------------
// Created on:	07 Mar 2005
// Last update: 07 Mar 2005
//----------------------------------------------------------------------------

#include <iostream>
#include <stack>
#include "my_except.h"
#include "parser.h" 
#include <cstdlib>


NIFS::c_parser::c_parser():
  		   m_max_size( 250 ),
         m_delims(" ,;\t\n"),
			   m_numers("0123456789.")
{
    //---( binary operators )---
    func_ptr2["add"] = &NIFS::c_parser::fadd;
    func_ptr2["sub"] = &NIFS::c_parser::fsub;
    func_ptr2["mul"] = &NIFS::c_parser::fmul;
    func_ptr2["div"] = &NIFS::c_parser::fdiv;
    func_ptr2["pow"] = &NIFS::c_parser::fpow;
    
    //---( unary operators )---
    func_ptr1["pos"]   = &NIFS::c_parser::fpos;
    func_ptr1["neg"]   = &NIFS::c_parser::fneg;
    func_ptr1["abs"]   = &NIFS::c_parser::fabs;
    func_ptr1["floor"] = &NIFS::c_parser::ffloor;
    func_ptr1["ceil"]  = &NIFS::c_parser::fceil;
    func_ptr1["sqrt"]  = &NIFS::c_parser::fsqrt;
    func_ptr1["ln"]    = &NIFS::c_parser::fln;
    func_ptr1["log"]   = &NIFS::c_parser::flog;
    func_ptr1["exp"]   = &NIFS::c_parser::fexp;    
    func_ptr1["cos"]   = &NIFS::c_parser::fcos;
    func_ptr1["sin"]   = &NIFS::c_parser::fsin;
    func_ptr1["tan"]   = &NIFS::c_parser::ftan;
    func_ptr1["acos"]  = &NIFS::c_parser::facos;
    func_ptr1["asin"]  = &NIFS::c_parser::fasin;
    func_ptr1["atan"]  = &NIFS::c_parser::fatan;
    func_ptr1["erf"]   = &NIFS::c_parser::ferf;
    func_ptr1["hmz"]   = &NIFS::c_parser::fhmz;
    func_ptr1["hmzd"]   = &NIFS::c_parser::fhmzd;
    
    clear();
}


void 
NIFS::c_parser::set_postfix( std::string pfstr )
{
    if ( pfstr.length() > m_max_size )
    {
		throw NIFS::my_except( "string passed to parser is too long." );
	}
    m_pfstr = pfstr;
    return;
}


void
NIFS::c_parser::tokenize_postfix()
{
    std::string::size_type begin_index;
    std::string::size_type end_index;
    std::string::size_type diff_index;
    
    convert_lowcase( m_pfstr );
    begin_index = m_pfstr.find_first_not_of( m_delims );

    while ( begin_index != std::string::npos ) {
        end_index = m_pfstr.find_first_of( m_delims, begin_index );
		if ( end_index == std::string::npos ) {
            end_index = m_pfstr.length();
		}
		diff_index = end_index - begin_index;
		m_pflist.push_back( m_pfstr.substr( begin_index, diff_index ) );
		begin_index = m_pfstr.find_first_not_of( m_delims, end_index );
    }
}


void
NIFS::c_parser::set_param( std::string name, double value )
{
    m_param[name] = static_cast<double>(value);
}


void
NIFS::c_parser::eval_postfix_list()
{
    std::stack< double > val_stack;
	double val;
	std::string name = "";
	std::list<std::string>::iterator iter;

	if (m_pflist.empty() ) {
		return;
	}
	--( iter = m_pflist.end() ); 
    if ( ( *iter ) == "eq" ) {
		m_pflist.pop_back();
		if (m_pflist.empty() ) {
			throw NIFS::my_except( "bad assignment" );
		}
		iter = m_pflist.begin(); 
		if ( ( iter->find_first_of( m_numers ) == 0 ) ||
			 ( func_ptr1.count( *iter ) > 0 ) ||
			 ( func_ptr2.count( *iter ) > 0 ) ||
			 ( ( *iter ) == "eq" ) ) {
			throw NIFS::my_except( "bad assignment" );
		} 
		else {
			name = ( *iter );
			m_pflist.pop_front();
		}
	}

    for ( iter = m_pflist.begin(); iter != m_pflist.end(); ++iter ) 
	{
		if ( iter->find_first_of( m_numers ) == 0 ) {
		    val_stack.push( atof( iter->c_str() ) );   	 
	    } 
		else if ( func_ptr1.count( *iter ) > 0 ) {
			if ( val_stack.empty() ) 
				throw NIFS::my_except( "no operand found for \"" + (*iter) + "\"");
			val_stack.top() = (this->*func_ptr1[*iter])( val_stack.top() );
		}
		else if ( func_ptr2.count( *iter ) > 0 ) { 
			if ( val_stack.empty() ) 
				throw NIFS::my_except( "no operand found for \"" + (*iter) + "\"" );
	 	    val = val_stack.top();
			val_stack.pop();
			if ( val_stack.empty() ) 
				throw NIFS::my_except( "2nd operand not found for \"" + (*iter) + "\"" );
			val_stack.top() = (this->*func_ptr2[*iter])( val_stack.top(), val );	        
		}
		else if ( ( *iter ) == "eq" ) 
			throw NIFS::my_except( 
				"operator \"eq\" must be at the end of the expression" );
		else if ( m_param.count( *iter ) > 0 ) { 
			val_stack.push( m_param[*iter] );
		}
		else 
			throw NIFS::my_except( "undefined token \"" + (*iter) + "\"" );
    }

    if ( val_stack.size() > 1 ) 
	    throw NIFS::my_except( "missed operator" );
	else if ( !val_stack.empty() ) {
		m_answer = static_cast<double>( val_stack.top() );
		if ( name.empty() )
    		m_param["ans"] = m_answer;
    	else
    		m_param[name] = m_answer;
	}
	return;
}


void 
NIFS::c_parser::clear( void )
{
	m_pfstr.clear();
	m_pflist.clear();
	m_answer = 0.0;
}


void 
NIFS::c_parser::convert_lowcase( std::string& str )
{
	for ( unsigned i = 0; i < str.length(); ++i )
		str[i] = tolower( str[i] );
}


double 
NIFS::c_parser::get_param( std::string name )
{
	convert_lowcase( name );
	if ( m_param.count( name ) > 0 ) 
		return m_param[name];
	else 
		throw NIFS::my_except( "parameter \"" + name + "\" not found" );
}


void 
NIFS::c_parser::eval_postfix( const std::string& pfstr )
{
	clear();
	set_postfix( pfstr );
	tokenize_postfix();
	eval_postfix_list();			

	return;
}


std::string 
NIFS::c_parser::get_token( const std::string& pfstr, unsigned n )
{
	clear();
	set_postfix( pfstr );
	tokenize_postfix();
	if ( m_pflist.empty() ) return "";
	std::string::size_type m = n % m_pflist.size();
	std::list< std::string >::iterator iter = m_pflist.begin();
	for ( unsigned i = 0; i < m; ++i ) 
		++iter;
	return (*iter);
}


double
NIFS::c_parser::ferf( double v )
{
  double pi = 3.1415926535;
  double a = ( 8.0 * ( pi - 3.0 ) ) / ( 3.0 * pi * ( 4.0 - pi ) );
  double r = std::sqrt( 1.0 - std::exp( ( -v * v ) * ( 4.0 / pi + a * v * v ) / 
                                        ( 1.0 + a * v * v ) ) );

  return ( v >= 0.0 )? r : -r;
}

double
NIFS::c_parser::fhmz( double v )
{
//  return std::pow( v , 2.207 ) / ( 1.074 + std::pow( v , 1.207 ) );
  return
   -8.6697e-04 * v +
    6.1929e-01 * pow(v,2) +
   -1.6883e-01 * pow(v,3) +
   -1.9845e-03 * pow(v,4) +
    1.6413e-02 * pow(v,5) +
   -5.6835e-03 * pow(v,6) +
    1.0445e-03 * pow(v,7) +
   -1.1740e-04 * pow(v,8) +
    8.0953e-06 * pow(v,9) +
   -3.1559e-07 * pow(v,10) +
    5.3393e-09 * pow(v,11);
}

double
NIFS::c_parser::fhmzd( double v )
{
//  return std::pow( 1.0 - std::exp( -1.249 * std::pow( v , 1.361 ) ) , 1/1.361 );
  if ( v < 8 )
    return  
    1.2327e-00 * v +
   -5.0124e-01 * pow(v,2) +
    5.2818e-03 * pow(v,3) +
    5.2895e-02 * pow(v,4) +
   -9.2187e-03 * pow(v,5) +
   -4.6229e-03 * pow(v,6) +
    2.6392e-03 * pow(v,7) +
   -6.2429e-04 * pow(v,8) +
    8.5265e-05 * pow(v,9) +
   -6.9956e-06 * pow(v,10) +
    3.2189e-07 * pow(v,11) +
   -6.4124e-09 * pow(v,12);
  else
    return 1.0;
}
