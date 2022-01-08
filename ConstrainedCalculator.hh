#pragma once

// Standard library includes
#include <algorithm>

// STV analysis includes
#include "SystematicsCalculator.hh"

// Calculates covariance matrices describing the uncertainty on the reco-space
// event counts. Distinguishes between "signal region" reco bins and "sideband"
// reco bins. For the former, covariances are calculated for both the total
// event count (signal + background) and background only. For the latter,
// covariances are calculated on just the total event count. Using the
// concept of "conditional covariance" (see MicroBooNE DocDB #32672), this
// class will calculate a constrained signal + background prediction
// together with its uncertainty. A constrained CV background prediction
// is also obtained that should be subtracted from the measured data points
// and the constrained signal + background prediction before unfolding.
class ConstrainedCalculator : public SystematicsCalculator {

  public:

    ConstrainedCalculator( const std::string& input_respmat_file_name,
      const std::string& syst_cfg_file_name = "",
      const std::string& respmat_tdirectoryfile_name = "" );

    virtual double evaluate_observable( const Universe& univ, int cm_bin,
      int flux_universe_index = -1 ) const override;

    virtual double evaluate_mc_stat_unc( const Universe& univ,
      int cm_bin ) const override;

    virtual double evaluate_data_stat_unc( int cm_bin,
      bool use_ext ) const override;

    // This class uses a dimension for the covariance matrix that is different
    // than just the number of reco bins, so we need to override this virtual
    // function.
    virtual size_t get_covariance_matrix_size() const override;

  protected:

    // Helper type used to decide how to handle the covariance calculation
    enum ConstrainedCalculatorBinType {
      // The bin is an "ordinary" one, and the observable of interest is the
      // total number of reconstructed events (signal + background)
      kOrdinaryRecoBinAll,

      // The bin is an "ordinary" one, and the observable of interest is the
      // number of reconstructed background events
      kOrdinaryRecoBinBkgd,

      // The bin is a "sideband" one, and the observable of interest is the
      // total number of reconstructed events (signal + background)
      kSidebandRecoBinAll
    };

    // Helper function for linking a covariance matrix bin number to
    // information about how to compute observables within it
    int get_reco_bin_and_type( int cm_bin,
      ConstrainedCalculatorBinType& bin_type ) const;

    // Number of "ordinary" reco bins (assumed to come first in the
    // ResponseMatrixMaker configuration file)
    size_t num_ordinary_reco_bins_ = 0u;

    // Number of "sideband" reco bins
    size_t num_sideband_reco_bins_ = 0u;
};

ConstrainedCalculator::ConstrainedCalculator(
  const std::string& input_respmat_file_name,
  const std::string& syst_cfg_file_name,
  const std::string& respmat_tdirectoryfile_name )
  : SystematicsCalculator( input_respmat_file_name,
  syst_cfg_file_name, respmat_tdirectoryfile_name )
{
  num_ordinary_reco_bins_ = 0u;
  num_sideband_reco_bins_ = 0u;
  bool found_first_sideband_bin = false;
  for ( const auto& rbin : reco_bins_ ) {
    if ( rbin.type_ == kOrdinaryRecoBin ) {
      ++num_ordinary_reco_bins_;
      if ( found_first_sideband_bin ) throw std::runtime_error( "Ordinary"
        " reco bins must precede sideband ones in the ResponseMatrixMaker"
        " configuration file." );
    }
    if ( rbin.type_ == kSidebandRecoBin ) {
      found_first_sideband_bin = true;
      ++num_sideband_reco_bins_;
    }
  }
}

double ConstrainedCalculator::evaluate_observable( const Universe& univ,
  int cm_bin, int flux_universe_index ) const
{
  // For the ConstrainedCalculator class, the observable of interest is the
  // total number of events (either signal + background or background only) in
  // the current bin in reco space
  double reco_bin_events = 0.;

  bool use_detVar_CV = this->is_detvar_universe( univ );

  // Get access to the CV universe. We need it regardless of the input universe
  // so that we can use it in the denominator of the smearceptance matrix
  // element. Note that we should use the detVarCV universe as the CV when the
  // input universe corresponds to a detector variation (or is the detVarCV
  // universe itself). Based on the check above, we assign a pointer to
  // either the regular or detVar CV here as appropriate.
  const Universe* cv_univ = nullptr;
  if ( use_detVar_CV ) {
    cv_univ = detvar_universes_.at( NFT::kDetVarMCCV ).get();
  }
  else {
    cv_univ = &this->cv_universe();
  }

  ConstrainedCalculatorBinType bin_type;
  int reco_bin = this->get_reco_bin_and_type( cm_bin, bin_type );

  size_t num_true_bins = true_bins_.size();

  // We need to sum the contributions of the various true bins,
  // so loop over them while checking whether each one is associated
  // with either signal or background
  for ( size_t tb = 0u; tb < num_true_bins; ++tb ) {
    const auto& tbin = true_bins_.at( tb );

    if ( tbin.type_ == kSignalTrueBin ) {

      // Get the CV event count for the current true bin
      double denom_CV = cv_univ->hist_true_->GetBinContent( tb + 1 );

      // For the systematic variation universes, we want to assess
      // uncertainties on the signal only through the smearceptance
      // matrix. We therefore compute the smearceptance matrix element
      // here and then apply it to the CV expected event count in
      // each true bin.
      // NOTE: ROOT histogram bin numbers are one-based (bin zero is always
      // the underflow bin). Our bin indices therefore need to be offset by
      // +1 in all cases here.
      double numer = univ.hist_2d_->GetBinContent( tb + 1, reco_bin + 1 );
      double denom = univ.hist_true_->GetBinContent( tb + 1 );

      // I plan to extract the flux-averaged cross sections in terms of the
      // *nominal* flux model (as opposed to the real flux). I therefore
      // vary the numerator of the smearceptance matrix for these while
      // keeping the denominator equal to the CV expectation under the
      // nominal flux model. This is the same strategy as is used in the
      // Wire-Cell CC inclusive analysis.
      if ( flux_universe_index >= 0 ) {
        denom = denom_CV;
      }

      // If the denominator is nonzero actually calculate the fraction.
      // Otherwise, just leave it zeroed out.
      // TODO: revisit this, think about MC statistical uncertainties
      // on the empty bins
      double smearcept = 0.;
      if ( denom > 0. ) smearcept = numer / denom;

      // Compute the expected signal events in this universe
      // by multiplying the varied smearceptance matrix element
      // by the unaltered CV prediction in the current true bin.
      double expected_signal = smearcept * denom_CV;

      // For the ordinary background reco bins, don't include any signal
      // contribution
      if ( bin_type == kOrdinaryRecoBinBkgd ) {
        expected_signal = 0.;
      }
      //// TODO: REVISIT THIS
      //// For the sideband bins, go ahead and vary the signal prediction
      //// just like the background one
      //else if ( bin_type == kSidebandRecoBinAll ) {
      //  expected_signal = numer;
      //}

      // Compute the expected signal events in the current reco bin
      // with the varied smearceptance matrix (and, for flux universes,
      // the varied integrated flux)
      reco_bin_events += expected_signal;
    }
    else if ( tbin.type_ == kBackgroundTrueBin ) {
      // For background events, we can use the same procedure regardless
      // of whether we're in the CV universe or not
      double background = univ.hist_2d_->GetBinContent( tb + 1, reco_bin + 1 );
      reco_bin_events += background;
    }
  } // true bins

  return reco_bin_events;
}

double ConstrainedCalculator::evaluate_mc_stat_unc( const Universe& univ,
  int cm_bin ) const
{
  ConstrainedCalculatorBinType bin_type;
  int reco_bin = this->get_reco_bin_and_type( cm_bin, bin_type );

  if ( bin_type != kOrdinaryRecoBinBkgd ) {
    // ROOT histograms use one-based bin indices, so I correct for that here
    double err = univ.hist_reco_->GetBinError( reco_bin + 1 );
    return err;
  }

  // Include only background bins when evaluating the MC stat uncertainty
  // for the "ordinary background" bin type
  size_t num_true_bins = true_bins_.size();
  double err2 = 0.;
  for ( size_t tb = 0u; tb < num_true_bins; ++tb ) {
    const auto& tbin = true_bins_.at( tb );
    if ( tbin.type_ == kBackgroundTrueBin ) {
      double bkgd_err = univ.hist_2d_->GetBinError( tb + 1, reco_bin + 1 );
      err2 += bkgd_err * bkgd_err;
    }
  }
  // We added the errors in quadrature over the true bins above, so take the
  // square root here to get the uncertainty
  double err = std::sqrt( std::max(0., err2) );
  return err;
}

double ConstrainedCalculator::evaluate_data_stat_unc( int cm_bin,
  bool use_ext ) const
{
  ConstrainedCalculatorBinType bin_type;
  int reco_bin = this->get_reco_bin_and_type( cm_bin, bin_type );

  const TH1D* d_hist = nullptr;
  if ( use_ext ) d_hist = data_hists_.at( NFT::kExtBNB ).get(); // EXT data
  else d_hist = data_hists_.at( NFT::kOnBNB ).get(); // BNB data
  // ROOT histograms use one-based bin indices, so I correct for that here
  double err = d_hist->GetBinError( reco_bin + 1 );

  // Don't evaluate the data statistical uncertainty for the "ordinary
  // background" bins unless we're considering EXT events
  if ( bin_type == kOrdinaryRecoBinBkgd && !use_ext ) {
    return 0.;
  }

  return err;
}

size_t ConstrainedCalculator::get_covariance_matrix_size() const {
  size_t cm_size = 2u * num_ordinary_reco_bins_;
  cm_size += num_sideband_reco_bins_;
  return cm_size;
}

// Note that zero-based bin indices are assumed by this function. The
// return value is the reco bin index to use for computing covariances.
// The bin_type variable will be loaded with information about the kind
// of bin that should be assumed by the ConstrainedCalculator class.
int ConstrainedCalculator::get_reco_bin_and_type( int cm_bin,
  ConstrainedCalculatorBinType& bin_type ) const
{
  if ( cm_bin < 0 ) throw std::runtime_error( "Negative bin index encountered"
    " in ConstrainedCalculator::get_reco_bin_and_type()" );

  int bin = 0;

  // Bin numbers within the first set of ordinary reco bins do not need
  // to be remapped
  if ( cm_bin < num_ordinary_reco_bins_ ) {
    bin = cm_bin;
    bin_type = kOrdinaryRecoBinAll;
  }
  // A second copy of the ordinary reco bins follows immediately after.
  // All bins after the copies are assumed to be sideband reco bins.
  else {
    bin = cm_bin - num_ordinary_reco_bins_;
    if ( cm_bin < 2*num_ordinary_reco_bins_ ) bin_type = kOrdinaryRecoBinBkgd;
    else bin_type = kSidebandRecoBinAll;
  }

  return bin;
}
