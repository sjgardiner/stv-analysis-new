#pragma once

// Standard library includes
#include <algorithm>

// STV analysis includes
#include "SystematicsCalculator.hh"

// Calculates covariance matrices describing the uncertainty on the central
// value prediction for the true event counts
class TruthSystematicsCalculator : public SystematicsCalculator {

  public:

    TruthSystematicsCalculator( const std::string& input_respmat_file_name,
      const std::string& syst_cfg_file_name = "",
      const std::string& respmat_tdirectoryfile_name = "" );

    virtual double evaluate_observable( const Universe& univ, int true_bin,
      int flux_universe_index = -1 ) const override;

    virtual double evaluate_mc_stat_covariance( const Universe& univ,
      int true_bin_a, int true_bin_b ) const override;

    virtual double evaluate_data_stat_covariance( int true_bin_a,
      int true_bin_b, bool use_ext ) const override;

    inline virtual size_t get_covariance_matrix_size() const override
      { return true_bins_.size(); }

};

TruthSystematicsCalculator::TruthSystematicsCalculator(
  const std::string& input_respmat_file_name,
  const std::string& syst_cfg_file_name,
  const std::string& respmat_tdirectoryfile_name )
  : SystematicsCalculator( input_respmat_file_name,
  syst_cfg_file_name, respmat_tdirectoryfile_name )
{

}

double TruthSystematicsCalculator::evaluate_observable( const Universe& univ,
  int true_bin, int /*flux_universe_index*/ ) const
{
  double true_events = univ.hist_true_->GetBinContent( true_bin + 1 );
  return true_events;
}

double TruthSystematicsCalculator::evaluate_mc_stat_covariance(
  const Universe& univ, int true_bin_a, int true_bin_b ) const
{
  // ROOT histograms use one-based bin indices, so I correct for that here.
  // Note that using the bin error (rather than the bin contents) enables a
  // correct treatment for weighted events provided TH1::Sumw2() was called
  // before filling the histogram.
  double err = univ.hist_true2d_->GetBinError( true_bin_a + 1, true_bin_b + 1 );
  double err2 = err*err;
  return err2;
}

double TruthSystematicsCalculator::evaluate_data_stat_covariance(
  int /*true_bin_a*/, int /*true_bin_b*/, bool /*use_ext*/ ) const
{
  // The predicted true event counts have no data stat uncertainty
  return 0.;
}
