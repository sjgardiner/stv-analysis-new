#pragma once

// Standard library includes
#include <memory>
#include <set>
#include <stdexcept>

// Simple container for the output of Unfolder::unfold()
struct UnfoldedMeasurement {
  UnfoldedMeasurement( TMatrixD* unfolded_signal, TMatrixD* cov_matrix,
    TMatrixD* unfolding_matrix, TMatrixD* err_prop_matrix,
    TMatrixD* add_smear_matrix ) : unfolded_signal_( unfolded_signal ),
    cov_matrix_( cov_matrix ), unfolding_matrix_( unfolding_matrix ),
    err_prop_matrix_( err_prop_matrix ), add_smear_matrix_( add_smear_matrix )
    {}

  std::unique_ptr< TMatrixD > unfolded_signal_;
  std::unique_ptr< TMatrixD > cov_matrix_;
  std::unique_ptr< TMatrixD > unfolding_matrix_;
  std::unique_ptr< TMatrixD > err_prop_matrix_;
  std::unique_ptr< TMatrixD > add_smear_matrix_;
};

// Abstract base class for objects that implement an algorithm for unfolding
// measured background-subtracted event counts from reco space to
// true space, possibly with regularization.
class Unfolder {

  public:

    Unfolder() {}

    // Function that actually implements a specific unfolding algorithm
    virtual UnfoldedMeasurement unfold( const TMatrixD& data_signal,
      const TMatrixD& data_covmat, const TMatrixD& smearcept,
      const TMatrixD& prior_true_signal ) const = 0;

  protected:

    // Helper function that does some sanity checks on the dimensions of the
    // input matrices passed to unfold()
    static void check_matrices( const TMatrixD& data_signal,
      const TMatrixD& data_covmat, const TMatrixD& smearcept,
      const TMatrixD& prior_true_signal );
};

void Unfolder::check_matrices( const TMatrixD& data_signal,
  const TMatrixD& data_covmat, const TMatrixD& smearcept,
  const TMatrixD& prior_true_signal )
{
  // Check the matrix dimensions for sanity
  int num_ordinary_reco_bins = smearcept.GetNrows();
  int num_true_signal_bins = smearcept.GetNcols();

  if ( data_signal.GetNcols() != 1 ) {
    throw std::runtime_error( "The background-subtracted data event counts"
      " must be expressed as a column vector" );
  }

  if ( prior_true_signal.GetNcols() != 1 ) {
    throw std::runtime_error( "The prior true signal event counts must be"
      " expressed as a column vector" );
  }

  if ( data_signal.GetNrows() != num_ordinary_reco_bins ) {
    throw std::runtime_error( "Reco bin mismatch between background-"
      "subtracted data and the smearceptance matrix" );
  }

  if ( data_covmat.GetNrows() != num_ordinary_reco_bins
    || data_covmat.GetNcols() != num_ordinary_reco_bins )
  {
    throw std::runtime_error( "Dimension mismatch between data covariance"
      " matrix and the smearceptance matrix" );
  }

  if ( prior_true_signal.GetNrows() != num_true_signal_bins ) {
    throw std::runtime_error( "Dimension mismatch between prior true signal"
      " event counts and the smearceptance matrix" );
  }

}
