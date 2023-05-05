#pragma once

// STV analysis includes
#include "FiducialVolume.hh"
#include "IntegratedFluxUniverseManager.hh"
#include "SystematicsCalculator.hh"

// Calculates covariance matrices describing the uncertainty on the
// MCC8-style forward-folded differential cross section
class MCC8ForwardFolder : public SystematicsCalculator {

  public:

    using NFT = NtupleFileType;

    MCC8ForwardFolder( const std::string& input_respmat_file_name,
      const std::string& syst_cfg_file_name = "",
      const std::string& respmat_tdirectoryfile_name = "" );

    virtual double evaluate_observable( const Universe& univ, int reco_bin,
      int flux_universe_index = -1 ) const override
    {
      // Note that the forward_folded_xsec() member function takes a
      // one-based reco bin index (to match the ROOT histograms), while
      // SystematicsCalculator::evaluate_observable() uses a zero-based
      // index (to match the owned bin definitions in that class). I
      // correct for this here in the call to forward_folded_xsec().
      return this->forward_folded_xsec( univ, reco_bin + 1,
        flux_universe_index );
    }

    virtual double evaluate_data_stat_covariance( int reco_bin_a,
      int reco_bin_b, bool use_ext ) const override;

    virtual double evaluate_mc_stat_covariance( const Universe& univ,
      int reco_bin_a, int reco_bin_b ) const override;

    // Helper functions for comparison to CCNp result
    double effective_efficiency( const Universe& univ, int reco_bin ) const;

    double xsec_scale_factor( int reco_bin,
      int flux_universe_index = -1 ) const;

    double expected_mc_background( const Universe& univ, int reco_bin,
      bool stat_var = false ) const;

    double forward_folded_xsec( const Universe& univ, int reco_bin,
      int flux_universe_index = -1 ) const;
};

MCC8ForwardFolder::MCC8ForwardFolder(
  const std::string& input_respmat_file_name,
  const std::string& syst_cfg_file_name,
  const std::string& respmat_tdirectoryfile_name )
  : SystematicsCalculator( input_respmat_file_name,
  syst_cfg_file_name, respmat_tdirectoryfile_name )
{

}

// NOTE: this uses a one-based reco bin index
double MCC8ForwardFolder::effective_efficiency(
  const Universe& univ, int reco_bin ) const
{
  int num_true_bins = true_bins_.size();
  int num_reco_bins = reco_bins_.size();

  if ( reco_bin > num_reco_bins ) {
    throw std::runtime_error( "Invalid reco bin" );
  }

  // Correct for the one-based reco bin index here. Access the reco bin
  // so that we can check its block index.
  const auto& rbin = reco_bins_.at( reco_bin - 1 );
  int reco_block_index = rbin.block_index_;

  double numerator = 0.;
  double denominator = 0.;
  for ( int tb = 1; tb <= num_true_bins; ++tb ) {
    // Note that the true bin definitions have a zero-based index
    auto& tbin = true_bins_.at( tb - 1 );
    // If this isn't a signal true bin, just skip it
    if ( tbin.type_ != kSignalTrueBin ) continue;

    // If the signal true bin is not in the same block as the requested
    // reco bin, also skip it. This avoids issues with double-counting.
    int true_block_index = tbin.block_index_;
    if ( reco_block_index != true_block_index ) continue;

    double num_j_gen = univ.hist_true_->GetBinContent( tb );
    double num_ij = univ.hist_2d_->GetBinContent( tb, reco_bin );

    double num_j_sel = 0.;
    for ( int rb = 1; rb <= num_reco_bins; ++rb ) {
      num_j_sel += univ.hist_2d_->GetBinContent( tb, rb );
    }

    double denom_term = num_ij * num_j_gen;
    if ( num_j_sel > 0. ) denom_term /= num_j_sel;
    else denom_term = 0.;

    numerator += num_ij;
    denominator += denom_term;
  }

  double eff = 0.;
  if ( denominator > 0. ) eff = numerator / denominator;

  return eff;
}

// NOTE: this uses a one-based reco bin index
double MCC8ForwardFolder::expected_mc_background(
  const Universe& univ, int reco_bin, bool stat_var ) const
{
  int num_true_bins = true_bins_.size();
  int num_reco_bins = reco_bins_.size();

  if ( reco_bin > num_reco_bins ) {
    throw std::runtime_error( "Invalid reco bin" );
  }

  double Bi = 0.;
  for ( int tb = 1; tb <= num_true_bins; ++tb ) {
    // Note that the true bin definitions have a zero-based index
    auto& tbin = true_bins_.at( tb - 1 );
    // If this isn't a background true bin, just skip it
    if ( tbin.type_ != kBackgroundTrueBin ) continue;

    // If this flag is set, then calculate the MC statistical variance
    // on the prediction instead of the event count
    if ( stat_var ) {
      double stat_err = univ.hist_2d_->GetBinError( tb, reco_bin );
      Bi += std::pow( stat_err, 2 );
    }
    else {
      Bi += univ.hist_2d_->GetBinContent( tb, reco_bin );
    }
  }

  return Bi;
}

// NOTE: this uses a one-based reco bin index
double MCC8ForwardFolder::xsec_scale_factor( int reco_bin,
  int flux_universe_index ) const
{
  const auto& ifum = IntegratedFluxUniverseManager::Instance();
  double numu_flux = ifum.cv_numu_integrated_flux(); // numu / POT / cm^2

  if ( flux_universe_index >= 0 ) {
    numu_flux *= ifum.flux_factor( flux_universe_index );
  }

  // Convert to numu / cm^2 by multiplying by the total analyzed beam data POT
  numu_flux *= total_bnb_data_pot_;

  const auto& cv_univ = this->cv_universe();
  // TODO: fix this. The universe histograms are expressed in terms of
  // reco bin number, so the widths are always trivially one. You need
  // to make this code aware of the physics units on the x-axis.
  double reco_bin_width = cv_univ.hist_reco_->GetBinWidth( reco_bin );

  double num_Ar_targets = num_Ar_targets_in_FV();

  double factor = 1. / ( numu_flux * num_Ar_targets * reco_bin_width );
  return factor;
}

// NOTE: this uses a one-based reco bin index
double MCC8ForwardFolder::forward_folded_xsec( const Universe& univ,
  int reco_bin, int flux_universe_index ) const
{
  const TH1D* bnb_hist = data_hists_.at( NFT::kOnBNB ).get();
  double data_counts = bnb_hist->GetBinContent( reco_bin );

  const TH1D* ext_hist = data_hists_.at( NFT::kExtBNB ).get();
  double ext_counts = ext_hist->GetBinContent( reco_bin );

  double mc_bkgd_counts = this->expected_mc_background( univ, reco_bin );

  double eff = this->effective_efficiency( univ, reco_bin );
  double scaling = this->xsec_scale_factor( reco_bin, flux_universe_index );

  double xsec = ( data_counts - ext_counts - mc_bkgd_counts ) * scaling / eff;
  return xsec;
}

double MCC8ForwardFolder::evaluate_mc_stat_covariance( const Universe& univ,
  int reco_bin_a, int reco_bin_b ) const
{
  // TODO: refactor to correctly treat possible correlations between
  // reco bins
  if ( reco_bin_a != reco_bin_b ) return 0.;

  // The input reco bin index is zero-based, but the helper functions defined
  // below use a one-based index. Correct for this here before proceeding.
  int rb_a = reco_bin_a + 1;

  double err2_a = this->expected_mc_background( univ, rb_a, true );

  double eff_a = this->effective_efficiency( univ, rb_a );
  double scaling_a = this->xsec_scale_factor( rb_a );

  if ( eff_a <= 0. ) err2_a = 0.;
  else {
    err2_a *= std::pow( scaling_a / eff_a, 2 );
  }

  return err2_a;
}

double MCC8ForwardFolder::evaluate_data_stat_covariance( int reco_bin_a,
  int reco_bin_b, bool use_ext ) const
{
  const TH2D* d_hist = nullptr;
  if ( use_ext ) d_hist = data_hists2d_.at( NFT::kExtBNB ).get(); // EXT data
  else d_hist = data_hists2d_.at( NFT::kOnBNB ).get(); // BNB data

  const auto& cv_univ = this->cv_universe();

  // The input reco bin indices is zero-based, but the ROOT histograms are
  // one-based, so I correct for that here. The helper functions below
  // also work with the one-based indices.
  int rb_a = reco_bin_a + 1;
  int rb_b = reco_bin_b + 1;

  double err = d_hist->GetBinError( rb_a, rb_b );
  double err2 = err * err;

  double eff_a = this->effective_efficiency( cv_univ, rb_a );
  double eff_b = this->effective_efficiency( cv_univ, rb_b );

  double scaling_a = this->xsec_scale_factor( rb_a );
  double scaling_b = this->xsec_scale_factor( rb_b );

  if ( eff_a <= 0. || eff_b <= 0. ) err2 = 0.;
  else err2 *= ( scaling_a / eff_a ) * ( scaling_b / eff_b );

  return err2;
}
