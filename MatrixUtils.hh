#pragma once

// Standard library includes
#include <cfloat>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>

// ROOT includes
#include "TDecompQRH.h"
#include "TMatrixD.h"

constexpr double DEFAULT_MATRIX_INVERSION_TOLERANCE = 1e-4;

std::unique_ptr< TMatrixD > invert_matrix( const TMatrixD& mat,
  const double inversion_tolerance = DEFAULT_MATRIX_INVERSION_TOLERANCE )
{
  // Pre-scale before inversion to avoid numerical problems. Here we choose a
  // scaling factor such that the smallest nonzero entry in the original matrix
  // has an absolute value of unity. Note the use of the zero-based element
  // indices for TMatrixD.
  constexpr double BIG_DOUBLE = std::numeric_limits<double>::max();
  double min_abs = BIG_DOUBLE;
  int num_bins = mat.GetNrows();
  for ( int a = 0; a < num_bins; ++a ) {
    for ( int b = 0; b < num_bins; ++b ) {
      double element = mat( a, b );
      double abs_el = std::abs( element );
      if ( abs_el > 0. && abs_el < min_abs ) min_abs = abs_el;
    }
  }

  // If all matrix elements are zero, then this scaling won't work
  // and something is wrong. Complain if this is the case.
  if ( min_abs == BIG_DOUBLE ) {
    throw std::runtime_error( "Cannot invert a null matrix" );
  }

  // We're ready. Do the actual pre-scaling here. Keep the first TMatrixD
  // and invert a copy. This allows us to check that the inversion was
  // successful below.
  auto inverse_matrix = std::make_unique< TMatrixD >( mat );
  double scaling_factor = 1. / min_abs;
  inverse_matrix->operator*=( scaling_factor );

  // Do the inversion. Here we try the QR method. Other options involve using
  // TDecompLU, TDecompBK, TDecompChol, TDecompQRH and TDecompSVD. The
  // Invert() member function of TMatrixDSym uses a TDecompBK object to do
  // the same thing internally.
  //inverse_matrix->Invert();
  TDecompQRH qr_decomp( *inverse_matrix, DBL_EPSILON );
  qr_decomp.Invert( *inverse_matrix );

  // Undo the scaling by re-applying it to the inverse matrix
  inverse_matrix->operator*=( scaling_factor );

  // Double-check that we get a unit matrix by multiplying the
  // original by its inverse
  TMatrixD unit_mat( mat, TMatrixD::kMult, *inverse_matrix );
  for ( int a = 0; a < num_bins; ++a ) {
    for ( int b = 0; b < num_bins; ++b ) {
      double element = unit_mat( a, b );
      double expected_element = 0.;
      if ( a == b ) expected_element = 1.;
      double abs_diff = std::abs( element - expected_element );
      if ( abs_diff > inversion_tolerance ) {
        throw std::runtime_error( "Matrix inversion failed" );
      }
    }
  }

  return inverse_matrix;
}

void dump_text_matrix( const std::string& output_file_name,
  const TMatrixD& matrix )
{
  // Open the output file and set up the output stream so that full numerical
  // precision is preserved in the ascii text representation
  std::ofstream out_file( output_file_name );
  out_file << std::scientific
    << std::setprecision( std::numeric_limits<double>::max_digits10 );

  int num_x_bins = matrix.GetNrows();
  int num_y_bins = matrix.GetNcols();

  out_file << "numXbins " << num_x_bins << '\n';
  out_file << "numYbins " << num_y_bins << '\n';
  out_file << "xbin  ybin  z\n";

  for ( int xb = 0; xb < num_x_bins; ++xb ) {
    for ( int yb = 0; yb < num_y_bins; ++yb ) {

      double z = matrix( xb, yb );

      // Use zero-based bin indices in the dump
      out_file << xb << "  " <<  yb << "  " << z << '\n';

    } // loop over columns (y bins)
  } // loop over rows (x bins)
}

void dump_text_column_vector( const std::string& output_file_name,
  const TMatrixD& matrix )
{
  // Open the output file and set up the output stream so that full numerical
  // precision is preserved in the ascii text representation
  std::ofstream out_file( output_file_name );
  out_file << std::scientific
    << std::setprecision( std::numeric_limits<double>::max_digits10 );

  int num_x_bins = matrix.GetNrows();
  int num_y_bins = matrix.GetNcols();

  if ( num_y_bins != 1 ) {
    throw std::runtime_error( "Input matrix is not a column vector" );
  }

  out_file << "numXbins " << num_x_bins;

  for ( int xb = 0; xb < num_x_bins; ++xb ) {

      double z = matrix( xb, 0 );

      // Use zero-based bin indices in the dump
      out_file << '\n' << xb << "  " << z;

  } // loop over rows (x bins)
}

// Load a TMatrixD object saved in a text file by a previous call to
// dump_text_matrix() or dump_text_column_vector()
TMatrixD load_matrix( const std::string& input_file_name ) {

  // Get the table of matrix element values
  std::ifstream matrix_table_file( input_file_name );

  // Peek at the file contents to decide whether we're working with
  // a matrix or a column vector
  std::string dummy;
  matrix_table_file >> dummy >> dummy >> dummy;
  bool is_matrix = ( dummy == "numYbins" );

  // Return to the beginning of the file for parsing
  matrix_table_file.seekg( 0 );

  // Get the matrix or vector dimensions from the header line(s)
  int num_x_bins, num_y_bins;

  matrix_table_file >> dummy >> num_x_bins;
  if ( is_matrix ) {
    matrix_table_file >> dummy >> num_y_bins;

    // Skip the next header line which contains the data column names
    std::getline( matrix_table_file, dummy );
  }
  else {
    num_y_bins = 1;
  }

  // Create a TMatrixD with the correct dimensions
  TMatrixD matrix( num_x_bins, num_y_bins );

  // Parse its contents from the remaining lines
  std::string line;
  while ( std::getline(matrix_table_file, line) ) {
    int bin1, bin2;
    double element;

    std::stringstream temp_ss( line );
    temp_ss >> bin1;
    if ( is_matrix ) {
      temp_ss >> bin2;
    }
    else {
      bin2 = 0;
    }
    temp_ss >> element;

    if ( bin1 < num_x_bins && bin2 < num_y_bins ) {
      matrix( bin1, bin2 ) = element;
    }
  }

  return matrix;
}

// Compute the direct sum of a vector of input TMatrixD objects
TMatrixD direct_sum( const std::vector< const TMatrixD* >& matrices ) {
  // Determine the dimensions of the direct sum of the input matrices
  int num_rows = 0;
  int num_cols = 0;
  for ( const auto& mat : matrices ) {
    num_rows += mat->GetNrows();
    num_cols += mat->GetNcols();
  }

  // Create storage for the direct sum matrix, filling with zeros to start
  TMatrixD dir_sum( num_rows, num_cols );
  dir_sum.Zero();

  // Fill it with the appropriate elements of the input matrices
  int start_row = 0;
  int start_col = 0;
  for ( const auto& mat : matrices ) {
    int mat_rows = mat->GetNrows();
    int mat_cols = mat->GetNcols();
    for ( int r = 0; r < mat_rows; ++r ) {
      for ( int c = 0; c < mat_cols; ++c ) {
        dir_sum( start_row + r, start_col + c ) = mat->operator()( r, c );
      }
    }

    start_row += mat_rows;
    start_col += mat_cols;
  }

  return dir_sum;
}

// Overloaded version for a pair of input matrices
TMatrixD direct_sum( const TMatrixD& m1, const TMatrixD& m2 ) {
  std::vector< const TMatrixD* > matrices = { &m1, &m2 };
  return direct_sum( matrices );
}
