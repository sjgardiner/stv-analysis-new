#pragma once

// Standard library includes
#include <algorithm>

// STV analysis includes
#include "SystematicsCalculator.hh"

// Calculates covariance matrices describing the uncertainty on the reco-space
// event counts. Uses a default "recipe" (represented by the default enum value
// MCC9SystematicsCalculator::SystMode::ForXSec) appropriate for unfolding
// via the Wiener-SVD or D'Agostini methods (e.g., cross-section systematic
// uncertainties are evaluated fully on backgrounds but enter only via the
// response matrix for signal events).
class MCC9SystematicsCalculator : public SystematicsCalculator {

  public:

    enum SystMode { ForXSec, VaryOnlyBackground, VaryOnlySignalResponse,
      VaryOnlySignal, VaryBackgroundAndSignalDirectly };

    MCC9SystematicsCalculator( const std::string& input_respmat_file_name,
      const std::string& syst_cfg_file_name = "",
      const std::string& respmat_tdirectoryfile_name = "" );

    virtual double evaluate_observable( const Universe& univ, int reco_bin,
      int flux_universe_index = -1 ) const override;

    virtual double evaluate_mc_stat_covariance( const Universe& univ,
      int reco_bin_a, int reco_bin_b ) const override;

    virtual double evaluate_data_stat_covariance( int reco_bin_a,
      int reco_bin_b, bool use_ext ) const override;

    inline void set_syst_mode( SystMode mode )
      { syst_mode_ = mode; }

  protected:

    // Default to using the standard recipe of calculating covariance matrix
    // elements for a cross-section measurement
    SystMode syst_mode_ = SystMode::ForXSec;

};

MCC9SystematicsCalculator::MCC9SystematicsCalculator(
  const std::string& input_respmat_file_name,
  const std::string& syst_cfg_file_name,
  const std::string& respmat_tdirectoryfile_name )
  : SystematicsCalculator( input_respmat_file_name,
  syst_cfg_file_name, respmat_tdirectoryfile_name )
{

}

double MCC9SystematicsCalculator::evaluate_observable( const Universe& univ,
  int reco_bin, int flux_universe_index ) const
{
  // For the MCC9SystematicsCalculator class, the observable of interest is the
  // total number of events (signal + background) in the current bin in reco
  // space
  double reco_bin_events = 0.;

  // Look up the requested reco bin so that we can determine its block index.
  // Contributions from signal true bins outside of its block will be ignored
  // to avoid double-counting.
  const auto& rbin = reco_bins_.at( reco_bin );
  int reco_block_index = rbin.block_index_;

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

  size_t num_true_bins = true_bins_.size();

  // We need to sum the contributions of the various true bins,
  // so loop over them while checking whether each one is associated
  // with either signal or background
  for ( size_t tb = 0u; tb < num_true_bins; ++tb ) {
    const auto& tbin = true_bins_.at( tb );

    if ( tbin.type_ == kSignalTrueBin ) {

      // Ignore contributions from signal true bins outside of the block for
      // the current reco bin. This avoids issues with double-counting.
      int true_block_index = tbin.block_index_;
      if ( reco_block_index != true_block_index ) continue;

      // Get the CV event count for the current true bin
      double denom_CV = cv_univ->hist_true_->GetBinContent( tb + 1 );

      // Get the CV expectation for the number of signal events from the
      // current true bin that fall into the current reco bin
      double numer_CV = cv_univ->hist_2d_->GetBinContent( tb + 1,
        reco_bin + 1 );

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

      double expected_signal = 0.;

      if ( syst_mode_ == SystMode::ForXSec
        || syst_mode_ == SystMode::VaryOnlySignalResponse )
      {
        // Compute the expected signal events in this universe
        // by multiplying the varied smearceptance matrix element
        // by the unaltered CV prediction in the current true bin.
        expected_signal = smearcept * denom_CV;
      }
      else if ( syst_mode_ == SystMode::VaryOnlySignal
             || syst_mode_ == SystMode::VaryBackgroundAndSignalDirectly )
      {
        // Use the current universe's expectation for the number of
        // reconstructed signal events from the current true bin that fall into
        // the current reco bin
        expected_signal = numer;
      }
      else if ( syst_mode_ == SystMode::VaryOnlyBackground ) {
        // Use the CV expectation for the number of reconstructed signal events
        // from the current true bin that fall into the current reco bin,
        // ignoring any systematic variation
        expected_signal = numer_CV;
      }
      else throw std::runtime_error( "Unrecognized SystMode enum value"
        " in MCC9SystematicsCalculator::evaluate_observable()" );

      // Compute the expected signal events in the current reco bin
      // with the varied smearceptance matrix (and, for flux universes,
      // the varied integrated flux)
      reco_bin_events += expected_signal;
    }
    else if ( tbin.type_ == kBackgroundTrueBin ) {

      double bkg = univ.hist_2d_->GetBinContent( tb + 1, reco_bin + 1 );
      double bkg_CV = cv_univ->hist_2d_->GetBinContent( tb + 1, reco_bin + 1 );

      if ( syst_mode_ == SystMode::ForXSec
        || syst_mode_ == SystMode::VaryOnlyBackground
        || syst_mode_ == SystMode::VaryBackgroundAndSignalDirectly )
      {
        // Use the current universe's expectation for the number of
        // background events from the current true bin that fall into
        // the current reco bin
        reco_bin_events += bkg;
      }
      else if ( syst_mode_ == SystMode::VaryOnlySignal
        || syst_mode_ == SystMode::VaryOnlySignalResponse )
      {
        // Use the CV universe's expectation for the number of
        // background events from the current true bin that fall into
        // the current reco bin, ignoring any systematic variation
        reco_bin_events += bkg_CV;
      }
      else throw std::runtime_error( "Unrecognized SystMode enum value"
        " in MCC9SystematicsCalculator::evaluate_observable()" );
    }
  } // true bins

  return reco_bin_events;
}

double MCC9SystematicsCalculator::evaluate_mc_stat_covariance(
  const Universe& univ, int reco_bin_a, int reco_bin_b ) const
{
  // ROOT histograms use one-based bin indices, so I correct for that here.
  // Note that using the bin error (rather than the bin contents) enables a
  // correct treatment for weighted events provided TH1::Sumw2() was called
  // before filling the histogram.
  double err = univ.hist_reco2d_->GetBinError( reco_bin_a + 1, reco_bin_b + 1 );
  double err2 = err * err;
  return err2;
}

double MCC9SystematicsCalculator::evaluate_data_stat_covariance( int reco_bin_a,
  int reco_bin_b, bool use_ext ) const
{
  const TH2D* d_hist = nullptr;
  if ( use_ext ) d_hist = data_hists2d_.at( NFT::kExtBNB ).get(); // EXT data
  else d_hist = data_hists2d_.at( NFT::kOnBNB ).get(); // BNB data
  // ROOT histograms use one-based bin indices, so I correct for that here.
  // Note that using the bin error (rather than the bin contents) enables a
  // correct treatment for weighted events provided TH1::Sumw2() was called
  // before filling the histogram.
  double err = d_hist->GetBinError( reco_bin_a + 1, reco_bin_b + 1 );
  double err2 = err * err;
  return err2;
}
