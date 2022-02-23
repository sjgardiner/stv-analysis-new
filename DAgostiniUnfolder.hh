#pragma once

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
    // data points through the unfolding procedure. Do it element-by-element.
    // Make a copy of the old error propagation matrix first since we'll
    // be updating the original in the loop below.
    TMatrixD old_err_prop_mat( err_prop_mat );

    for ( int t = 0; t < num_true_signal_bins; ++t ) {
      for ( int r = 0; r < num_ordinary_reco_bins; ++r ) {
        // Start with the first two terms of the updated matrix element. These
        // can be evaluated without another for loop.
        double temp_el = unfold_mat->operator()( t, r );
        double ots = old_true_signal( t, 0 );
        if ( ots > 0. ) {
          temp_el += old_err_prop_mat( t, r )
            * true_signal->operator()( t, 0 ) / ots;
        }

        // The last term involves a sum that we evaluate using for loops
        double last_term = 0.;
        for ( int t2 = 0; t2 < num_true_signal_bins; ++t2 ) {
          double ots2 = old_true_signal( t2, 0 );
          if ( ots2 <= 0. ) continue;

          for ( int r2 = 0; r2 < num_ordinary_reco_bins; ++r2 ) {
            last_term += eff_vec( t2, 0 ) * data_signal( r2, 0 )
              * unfold_mat->operator()( t, r2 )
              * unfold_mat->operator()( t2, r2 )
              * old_err_prop_mat( t2, r ) / ots2;
          }
        }

        // Note that the last term has an overall minus sign (we do the
        // subtraction on the line below)
        temp_el -= last_term;

        // We're ready. Update the current element of the error propagation
        // matrix
        err_prop_mat( t, r ) = temp_el;
      }
    }

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
