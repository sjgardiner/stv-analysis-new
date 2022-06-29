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

    // Sets the value of the flag which determines whether MC statistical
    // uncertainties on the response (smearceptance) matrix elements should
    // be propagated through the unfolding procedure. This is more correct,
    // but it is often a small effect and relatively slow to compute.
    inline void set_include_respmat_covariance( bool do_it )
      { include_respmat_covariance_ = do_it; }

  protected:

    unsigned int num_iterations_;
    bool include_respmat_covariance_ = false;
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
  auto* err_prop_mat = new TMatrixD( num_true_signal_bins,
    num_ordinary_reco_bins );

  err_prop_mat->Zero(); // Zero out the elements (just in case)

  // We need a 3D tensor to do the propagation of MC uncertainties. It is
  // convenient in this case to represent it as a vector of TMatrixD objects.
  auto err_prop_mc_vec = std::vector< TMatrixD >();

  // Only populate the 3D tensor if it is needed (i.e., because we're including
  // the MC uncertainties on the final result)
  if ( include_respmat_covariance_ ) {
    for ( int s = 0; s < num_true_signal_bins; ++s ) {
      // Using emplace_back() here returns a reference to the new vector
      // element. We will use it below to explicitly zero out the matrix
      // contents (just in case). Use the same matrix element indexing scheme
      // as the smearceptance matrix (i.e., rows are ordinary reco bins,
      // columns are true signal bins).
      auto& mat_ref = err_prop_mc_vec.emplace_back( num_ordinary_reco_bins,
        num_true_signal_bins );

      mat_ref.Zero();
    }
  }

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
    TMatrixD temp_mat1( *err_prop_mat );

    // If the second argument passed to TMatrixD::NormByColumn() or
    // TMatrixD::NormByRow() is "D", then the matrix elements will be
    // divided by the input TVectorD elements rather than multiplied
    // by them. Following RooUnfold, I use "M" here to stress that
    // multiplication is the intended behavior.
    temp_mat1.NormByColumn( new_iter_over_old_iter, "M" );

    TMatrixD temp_mat2(
      TMatrixD::EMatrixCreatorsOp1::kTransposed, *unfold_mat );

    // The measured data points already exist as a one-column TMatrixD
    // rather than a TVectorD. To avoid an unnecessary copy, we scale
    // the columns of temp_mat2 here manually rather than calling
    // TMatrixD::NormByColumn()
    for ( int r2 = 0; r2 < num_ordinary_reco_bins; ++r2 ) {
      for ( int t2 = 0; t2 < num_true_signal_bins; ++t2 ) {
        temp_mat2( r2, t2 ) = temp_mat2( r2, t2 ) * data_signal( r2, 0 );
      }
    }

    temp_mat2.NormByRow( minus_eff_over_old_iter, "M" );

    TMatrixD temp_mat3( temp_mat2,
      TMatrixD::EMatrixCreatorsOp2::kMult, *err_prop_mat );

    // We're ready. Update the measurement error propagation matrix.
    err_prop_mat->Mult( *unfold_mat, temp_mat3 );
    err_prop_mat->operator+=( *unfold_mat );
    err_prop_mat->operator+=( temp_mat1 );

    // Also update the 3D tensor needed to propagate the MC statistical
    // uncertainties on the smearceptance matrix through the unfolding
    // procedure. As was done for the other error propagation matrix above,
    // first make a full copy of the tensor, then update the original.
    auto old_err_prop_mc_vec = std::vector< TMatrixD >();

    // Only update the 3D tensor if it is needed (i.e., because we're including
    // the MC uncertainties on the final result)
    if ( include_respmat_covariance_ ) {

      for ( int s = 0; s < num_true_signal_bins; ++s ) {
        old_err_prop_mc_vec.emplace_back( err_prop_mc_vec.at(s) );
      }

      // The outer loop over true signal bins runs over the bins used to report
      // the unfolded measurement. The inner loops run over the elements of the
      // smearceptance matrix.
      for ( int t = 0; t < num_true_signal_bins; ++t ) {
        for ( int r = 0; r < num_ordinary_reco_bins; ++r ) {
          for ( int t2 = 0; t2 < num_true_signal_bins; ++t2 ) {
            // Handle the Kronecker delta in the first term using an if
            // statement
            double temp_el = 0.;
            if ( t == t2 ) {
              double aux1 = ( data_signal(r, 0) * old_true_signal(t, 0)
                / reco_expected(r, 0) ) - true_signal->operator()( t, 0 );
              temp_el += aux1 / eff_vec( t, 0 );
            }

            temp_el -= data_signal( r, 0 ) * old_true_signal( t2, 0 )
              * unfold_mat->operator()( t, r ) / reco_expected( r, 0 );

            temp_el += true_signal->operator()( t, 0 )
              * old_err_prop_mc_vec.at( t )( r, t2 ) / old_true_signal( t, 0 );

            // These for loops take care of the sum in the last term
            double last_term = 0.;
            for ( int t3 = 0; t3 < num_true_signal_bins; ++t3 ) {
              for ( int r2 = 0; r2 < num_ordinary_reco_bins; ++r2 ) {
                last_term += data_signal( r2, 0 ) * eff_vec( t3, 0 )
                  * unfold_mat->operator()( t3, r2 )
                  * unfold_mat->operator()( t, r2 )
                  * old_err_prop_mc_vec.at( t3 )( r, t2 )
                  / old_true_signal( t3, 0 );
              }
            }

            // We account for the overall minus sign on the last term by
            // subtracting on the line below
            temp_el -= last_term;

            // We're ready. Update the MC error propagation matrix.
            err_prop_mc_vec.at( t )( r, t2 ) = temp_el;
          }
        }
      }

    }

  } // D'Agostini method iterations

  // Now that we're finished with the iterations, we can also transform the
  // data covariance matrix to the unfolded true space using the error
  // propagation matrix.
  TMatrixD err_prop_mat_tr( TMatrixD::kTransposed, *err_prop_mat );
  TMatrixD temp_mat( data_covmat, TMatrixD::EMatrixCreatorsOp2::kMult,
    err_prop_mat_tr );

  auto* true_signal_covmat = new TMatrixD( *err_prop_mat,
    TMatrixD::EMatrixCreatorsOp2::kMult, temp_mat );


  if ( include_respmat_covariance_ ) {

    // Here we also calculate a contribution to the covariance matrix on the
    // unfolded result that comes from the MC statistical uncertainty on the
    // smearceptance matrix elements
    TMatrixD mc_covmat( num_true_signal_bins, num_true_signal_bins );

    for ( int t = 0; t < num_true_signal_bins; ++t ) {

      const TMatrixD& prop_mat1 = err_prop_mc_vec.at( t );

      for ( int t2 = 0; t2 < num_true_signal_bins; ++t2 ) {

        const TMatrixD& prop_mat2 = err_prop_mc_vec.at( t2 );

        double temp_elem = 0.;

        for ( int t3 = 0; t3 < num_true_signal_bins; ++t3 ) {

          double prior_sig = prior_true_signal( t3, 0 );
          if ( prior_sig <= 0. ) continue;

          for ( int r = 0; r < num_ordinary_reco_bins; ++r ) {

            double smear1 = smearcept( r, t3 );

            for ( int r2 = 0; r2 < num_ordinary_reco_bins; ++r2 ) {

              double smear2 = smearcept( r2, t3 );

              // Calculate the covariance matrix element for the smearceptance
              // matrix elements. Assume independent multinomial distributions
              // for each true bin (as D'Agostini does)
              double covariance = 0.;
              // TODO: Account for effective statistics when using
              // weighted events
              // TODO: Account for situations in which the prior
              // differs from the true event counts used to compute
              // the smearceptance matrix elements
              if ( r == r2 ) {
                covariance = smear1 * ( 1. - smear1 ) / prior_sig;
              }
              else {
                covariance = -1. * smear1 * smear2 / prior_sig;
              }

              temp_elem += prop_mat1( r, t3 ) * covariance
                * prop_mat2( r2, t3 );
            }
          }
        }

        mc_covmat( t, t2 ) = temp_elem;
      }
    }

    // Add the MC statistical uncertainty to the other uncertainties to obtain
    // the final covariance matrix on the unfolded measurement
    true_signal_covmat->operator+=( mc_covmat );
  }

  // Compute the additional smearing matrix
  auto* add_smear = new TMatrixD( *unfold_mat,
    TMatrixD::EMatrixCreatorsOp2::kMult, smearcept );

  UnfoldedMeasurement result( true_signal, true_signal_covmat,
    unfold_mat.release(), err_prop_mat, add_smear );
  return result;
}
