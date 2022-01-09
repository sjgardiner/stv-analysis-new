#pragma once

// Standard library includes
#include <memory>

// STV analysis includes
#include "SystematicsCalculator.hh"

// Simple container for the output of Unfolder::unfold()
struct UnfoldedMeasurement {
  UnfoldedMeasurement( TMatrixD* unfolded_signal, TMatrixD* cov_matrix )
    : unfolded_signal_( unfolded_signal ), cov_matrix_( cov_matrix ) {}

  std::unique_ptr< TMatrixD > unfolded_signal_;
  std::unique_ptr< TMatrixD > cov_matrix_;
};

// Abstract base class for objects that implement an algorithm for unfolding
// measured background-subtracted event counts from reco space to
// true space, possibly with regularization.
class Unfolder {

  public:

    Unfolder() {}

    virtual UnfoldedMeasurement unfold( const TMatrixD& data_signal,
      const TMatrixD& data_covmat, const TMatrixD& smearcept,
      const TMatrixD& prior_true_signal ) const = 0;

    virtual UnfoldedMeasurement unfold(
      const SystematicsCalculator& syst_calc ) const final;
};

UnfoldedMeasurement Unfolder::unfold(
  const SystematicsCalculator& syst_calc ) const
{
  auto smearcept = syst_calc.get_cv_smearceptance_matrix();
  auto true_signal = syst_calc.get_cv_true_signal();
  auto meas = syst_calc.get_measured_events();
  const auto& data_signal = meas.reco_signal_;
  const auto& data_covmat = meas.cov_matrix_;

  return this->unfold( *meas.reco_signal_, *meas.cov_matrix_,
    *smearcept, *true_signal );
}
