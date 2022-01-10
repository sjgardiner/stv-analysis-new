#pragma once

// Standard library includes
#include <cfloat>
#include <cmath>
#include <limits>
#include <memory>

// ROOT includes
#include "TDecompQRH.h"
#include "TMatrixD.h"

std::unique_ptr< TMatrixD > invert_matrix( const TMatrixD& mat ) {

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
  constexpr double INVERSION_TOLERANCE = 1e-4;
  for ( int a = 0; a < num_bins; ++a ) {
    for ( int b = 0; b < num_bins; ++b ) {
      double element = unit_mat( a, b );
      double expected_element = 0.;
      if ( a == b ) expected_element = 1.;
      double abs_diff = std::abs( element - expected_element );
      if ( abs_diff > INVERSION_TOLERANCE ) throw std::runtime_error(
        "Matrix inversion failed" );
    }
  }

  return inverse_matrix;
}
