#pragma once

// ROOT includes
#include "TVectorD.h"

// STV analysis includes
#include "Unfolder.hh"

// Implementation of the iterative D'Agostini unfolding method
// G. D'Agostini, Nucl. Instrum. Methods Phys. Res. A 362, 487-498 (1995)
// https://hep.physics.utoronto.ca/~orr/wwwroot/Unfolding/d-agostini.pdf.
//
// Uncertainties on the measurement are propagated through the unfolding
// procedure following a corrected expression given in Eq. (20) from
// https://arxiv.org/abs/1806.03350. The original paper by D'Agostini
// propagates the measurement uncertainties as if the unfolding matrix is
// independent of the measured data points, but this is only true for the first
// iteration.
class DAgostiniUnfolder : public Unfolder {

  public:

    DAgostiniUnfolder( unsigned int default_iterations = 1u )
      : Unfolder(), num_iterations_( default_iterations ) {}

    // Trick taken from https://stackoverflow.com/a/18100999
    using Unfolder::unfold;

    virtual UnfoldedMeasurement unfold( const TMatrixD& data_signal,
      const TMatrixD& data_covmat, const TMatrixD& smearcept,
      const TMatrixD& prior_true_signal ) const override;

    inline unsigned int get_iterations() const { return num_iterations_; }
    inline void set_iterations( unsigned int iters )
      { num_iterations_ = iters; }

  protected:

    unsigned int num_iterations_;
};

UnfoldedMeasurement DAgostiniUnfolder::unfold( const TMatrixD& data_signal,
  const TMatrixD& data_covmat, const TMatrixD& smearcept,
  const TMatrixD& prior_true_signal ) const
{
  // Check input matrix dimensions for sanity
  this->check_matrices( data_signal, data_covmat,
    smearcept, prior_true_signal );

  int num_ordinary_reco_bins = smearcept.GetNrows();
  int num_true_signal_bins = smearcept.GetNcols();

  // Copy the measured data into a TVectorD for easy use with the
  // TMatrixD::NormByColumn() function in the unfolding loop
  TVectorD data_signal_vec( num_ordinary_reco_bins,
    data_signal.GetMatrixArray() );

  // Start the iterations by copying the prior on the true signal
  auto* true_signal = new TMatrixD( prior_true_signal );

  // Precompute the efficiency in each true bin by summing over reco bins in
  // the smearceptance matrix. Note that the smearceptance matrix remains
  // constant across iterations, and thus the efficiency does as well. Note
  // that the precomputed efficiencies are expressed as a column vector.
  TMatrixD eff_vec( num_true_signal_bins, 1 );
  for ( int t = 0; t < num_true_signal_bins; ++t ) {
    double efficiency = 0.;
    for ( int r = 0; r < num_ordinary_reco_bins; ++r ) {
      efficiency += smearcept( r, t );
    }
    eff_vec( t, 0 ) = efficiency;
  }

  // Create the unfolding matrix (always applied to the original
  // background-subtracted data). Note that its rows correspond to true signal
  // bins and its columns correspond to ordinary reco bins. This matrix will be
  // updated upon every iteration.
  auto unfold_mat = std::make_unique< TMatrixD >( num_true_signal_bins,
    num_ordinary_reco_bins );

  // Also create a matrix that we'll need to propagate the measurement
  // uncertainties through the unfolding procedure. The iterations introduce
  // non-trivial correlations which cause this matrix to differ from the
  // unfolding matrix itself.
  TMatrixD err_prop_mat( num_true_signal_bins, num_ordinary_reco_bins );

  err_prop_mat.Zero(); // Zero out the elements (just in case)

  // Start the iterations for the D'Agostini method
  for ( int it = 0; it < num_iterations_; ++it ) {

    // Compute the column vector of expected reco-space signal event counts
    // given the (fixed) smearceptance matrix and the current estimate of the
    // signal events in true-space
    TMatrixD reco_expected( smearcept, TMatrixD::EMatrixCreatorsOp2::kMult,
      *true_signal );

    // Update the unfolding matrix for the current iteration
    for ( int t = 0; t < num_true_signal_bins; ++t ) {
      for ( int r = 0; r < num_ordinary_reco_bins; ++r ) {
        double smearcept_element = smearcept( r, t );
        double true_signal_element = true_signal->operator()( t, 0 );
        double efficiency = eff_vec( t, 0 );
        double reco_expected_element = reco_expected( r, 0 );

        double unfold_mat_element = ( smearcept_element * true_signal_element )
          / ( efficiency * reco_expected_element );
        unfold_mat->operator()( t, r ) = unfold_mat_element;
      }
    }

    // Temporarily keep the true signal event counts from the previous
    // iteration (needed to update the error propagation matrices below)
    TMatrixD old_true_signal( *true_signal );

    // Update the estimated true signal event counts by applying the current
    // unfolding matrix to the background-subtracted data
    delete true_signal;
    true_signal = new TMatrixD( *unfold_mat,
      TMatrixD::EMatrixCreatorsOp2::kMult, data_signal );

    // Update the matrix used to propagate the uncertainties on the measured
    // data points through the unfolding procedure. We used to do this
    // element-by-element. The faster (and equivalent) procedure below was
    // taken from the RooUnfold implementation (http://roounfold.web.cern.ch/)

    // Initialize some vectors of factors that will be applied to rows and
    // columns of matrices below
    TVectorD minus_eff_over_old_iter( num_true_signal_bins );
    TVectorD new_iter_over_old_iter( num_true_signal_bins );

    for ( int t = 0; t < num_true_signal_bins; ++t ) {
      double ots = old_true_signal( t, 0 );
      if ( ots <= 0. ) {
        minus_eff_over_old_iter( t ) = 0.;
        new_iter_over_old_iter( t ) = 0.;
        continue;
      }

      double ts = true_signal->operator()( t, 0 );
      double eff = eff_vec( t, 0 );

      minus_eff_over_old_iter( t ) = -eff / ots;
      new_iter_over_old_iter( t ) = ts / ots;
    }

    // Duplicate the existing error propagation matrix so that we can
    // update the existing version in place
    TMatrixD temp_mat1( err_prop_mat );

    // If the second argument passed to TMatrixD::NormByColumn() or
    // TMatrixD::NormByRow() is "D", then the matrix elements will be
    // divided by the input TVectorD elements rather than multiplied
    // by them. Following RooUnfold, I use "M" here to stress that
    // multiplication is the intended behavior.
    temp_mat1.NormByColumn( new_iter_over_old_iter, "M" );

    TMatrixD temp_mat2(
      TMatrixD::EMatrixCreatorsOp1::kTransposed, *unfold_mat );

    temp_mat2.NormByColumn( data_signal_vec, "M" );
    temp_mat2.NormByRow( minus_eff_over_old_iter, "M" );

    TMatrixD temp_mat3( temp_mat2,
      TMatrixD::EMatrixCreatorsOp2::kMult, err_prop_mat );

    // We're ready. Update the measurement error propagation matrix.
    err_prop_mat.Mult( *unfold_mat, temp_mat3 );
    err_prop_mat += *unfold_mat;
    err_prop_mat += temp_mat1;

  } // D'Agostini method iterations

  // Now that we're finished with the iterations, we can also transform the
  // data covariance matrix to the unfolded true space using the error
  // propagation matrix.
  TMatrixD err_prop_mat_tr( TMatrixD::kTransposed, err_prop_mat );
  TMatrixD temp_mat( data_covmat, TMatrixD::EMatrixCreatorsOp2::kMult,
    err_prop_mat_tr );

  // TODO: add propagation of MC errors on the smearceptance matrix elements
  // here

  auto* true_signal_covmat = new TMatrixD( err_prop_mat,
    TMatrixD::EMatrixCreatorsOp2::kMult, temp_mat );

  UnfoldedMeasurement result( true_signal, true_signal_covmat );
  return result;
}
