#pragma once

// STV analysis includes
#include "IntegratedFluxUniverseManager.hh"
#include "SystematicsCalculator.hh"

// Calculates covariance matrices describing the uncertainty on the
// MCC8-style forward-folded differential cross section
class MCC8ForwardFolder : public SystematicsCalculator {

  public:

    using NFT = NtupleFileType;

    MCC8ForwardFolder( const std::string& input_respmat_file_name,
      const std::string& respmat_tdirectoryfile_name = "" );

    virtual void calculate_covariances( const std::string& cov_type,
      CovMatrix& cov_mat, std::istream& config_file ) const override;

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
  const std::string& respmat_tdirectoryfile_name )
  : SystematicsCalculator( input_respmat_file_name,
  respmat_tdirectoryfile_name )
{

}

void MCC8ForwardFolder::calculate_covariances(
  const std::string& cov_type, CovMatrix& cov_mat,
  std::istream& config_file ) const
{
  if ( cov_type == "MCstat" ) {

    int num_reco_bins = reco_bins_.size();
    const auto& cv_univ = this->cv_universe();

    // TODO: optimize so that you can use the Sumw2 array entries instead.
    // This avoids needing to take many square roots only to square them
    // again.
    // NOTE: The reco bin index used here is one-based since ROOT histograms
    // always include an underflow bin
    for ( int rb = 1; rb <= num_reco_bins; ++rb ) {

      double err2 = this->expected_mc_background( cv_univ, rb, true );

      double eff = this->effective_efficiency( cv_univ, rb );
      double scaling = this->xsec_scale_factor( rb );

      if ( eff <= 0. ) err2 = 0.;
      else {
        err2 *= std::pow( scaling / eff, 2 );
      }

      cov_mat.cov_matrix_->SetBinContent( rb, rb, err2 );

    } // reco bins

  } // MCstat type

  else if ( cov_type == "BNBstat" || cov_type == "EXTstat" ) {

    const TH1D* d_hist = nullptr;
    if ( cov_type == "BNBstat" ) d_hist = data_hists_.at( NFT::kOnBNB ).get();
    else d_hist = data_hists_.at( NFT::kExtBNB ).get(); // EXTstat

    int num_reco_bins = d_hist->GetNbinsX();

    const auto& cv_univ = this->cv_universe();

    // Note the one-based bin numbering convention for TH1D
    for ( int rb = 1; rb <= num_reco_bins; ++rb ) {
      double err = d_hist->GetBinError( rb );
      double err2 = err * err;

      // TODO: reduce code duplication with the MCstat type
      double eff = this->effective_efficiency( cv_univ, rb );
      double scaling = this->xsec_scale_factor( rb );

      if ( eff <= 0. ) err2 = 0.;
      else {
        err2 *= std::pow( scaling / eff, 2 );
      }

      cov_mat.cov_matrix_->SetBinContent( rb, rb, err2 );
    } // reco bins

  } // BNBstat and EXTstat types

  else {
    throw std::runtime_error( "Unrecognized covariance matrix type \""
      + cov_type + '\"' );
  }
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

  double numerator = 0.;
  double denominator = 0.;
  for ( int tb = 1; tb <= num_true_bins; ++tb ) {
    // Note that the true bin definitions have a zero-based index
    auto& tbin = true_bins_.at( tb - 1 );
    // If this isn't a signal true bin, just skip it
    if ( tbin.type_ != kSignalTrueBin ) continue;

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

  constexpr double NUM_TARGETS_ACTIVE_VOL = 1.0068e30;

  double factor = 1. / ( numu_flux * NUM_TARGETS_ACTIVE_VOL * reco_bin_width );
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
