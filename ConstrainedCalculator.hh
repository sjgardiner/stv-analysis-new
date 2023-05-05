#pragma once

// Standard library includes
#include <algorithm>
#include <stdexcept>

// STV analysis includes
#include "MatrixUtils.hh"
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

    virtual double evaluate_mc_stat_covariance( const Universe& univ,
      int cm_bin_a, int cm_bin_b ) const override;

    virtual double evaluate_data_stat_covariance( int cm_bin_a, int cm_bin_b,
      bool use_ext ) const override;

    // This class uses a dimension for the covariance matrix that is different
    // than just the number of reco bins, so we need to override this virtual
    // function.
    virtual size_t get_covariance_matrix_size() const override;

    // Apply a sideband constraint before subtracting the EXT+MC background.
    // NOTE: this function assumes that the ordinary reco bins are all listed
    // before any sideband reco bins
    virtual MeasuredEvents get_measured_events() const override;

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
  num_sideband_reco_bins_ = 0u;
  bool found_first_sideband_bin = false;
  for ( const auto& rbin : reco_bins_ ) {
    if ( rbin.type_ == kOrdinaryRecoBin ) {
      if ( found_first_sideband_bin ) throw std::runtime_error( "Ordinary"
        " reco bins must precede sideband ones in the UniverseMaker"
        " configuration file." );
    }
    else if ( rbin.type_ == kSidebandRecoBin ) {
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

  // Look up the block index for the current reco bin. We will use this
  // to avoid double-counting when summing over true bins below.
  const auto& rbin = reco_bins_.at( reco_bin );
  int reco_block_index = rbin.block_index_;

  size_t num_true_bins = true_bins_.size();

  // We need to sum the contributions of the various true bins,
  // so loop over them while checking whether each one is associated
  // with either signal or background
  for ( size_t tb = 0u; tb < num_true_bins; ++tb ) {
    const auto& tbin = true_bins_.at( tb );

    if ( tbin.type_ == kSignalTrueBin ) {

      // Ignore signal true bins outside of the same block as the current
      // reco bin. This avoids issues with double-counting.
      int true_block_index = tbin.block_index_;
      if ( reco_block_index != true_block_index ) continue;

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

double ConstrainedCalculator::evaluate_mc_stat_covariance( const Universe& univ,
  int cm_bin_a, int cm_bin_b ) const
{
  ConstrainedCalculatorBinType bin_type_a, bin_type_b;
  int reco_bin_a = this->get_reco_bin_and_type( cm_bin_a, bin_type_a );
  int reco_bin_b = this->get_reco_bin_and_type( cm_bin_b, bin_type_b );

  if ( bin_type_a != kOrdinaryRecoBinBkgd
    && bin_type_b != kOrdinaryRecoBinBkgd )
  {
    // ROOT histograms use one-based bin indices, so I correct for that here
    double err = univ.hist_reco2d_->GetBinError( reco_bin_a + 1,
      reco_bin_b + 1 );
    double err2 = err * err;
    return err2;
  }
  else {
    // TODO: add proper handling of MC stat correlations between "ordinary
    // background" bins
    if ( cm_bin_a != cm_bin_b ) return 0.;
  }

  // Include only background bins when evaluating the MC stat uncertainty
  // for the "ordinary background" bin type
  size_t num_true_bins = true_bins_.size();
  double err2 = 0.;
  for ( size_t tb = 0u; tb < num_true_bins; ++tb ) {
    const auto& tbin = true_bins_.at( tb );
    if ( tbin.type_ == kBackgroundTrueBin ) {
      double bkgd_err = univ.hist_2d_->GetBinError( tb + 1, reco_bin_a + 1 );
      err2 += bkgd_err * bkgd_err;
    }
  }

  return err2;
}

double ConstrainedCalculator::evaluate_data_stat_covariance( int cm_bin_a,
  int cm_bin_b, bool use_ext ) const
{
  ConstrainedCalculatorBinType bin_type_a, bin_type_b;
  int reco_bin_a = this->get_reco_bin_and_type( cm_bin_a, bin_type_a );
  int reco_bin_b = this->get_reco_bin_and_type( cm_bin_b, bin_type_b );

  const TH2D* d_hist = nullptr;
  if ( use_ext ) d_hist = data_hists2d_.at( NFT::kExtBNB ).get(); // EXT data
  else d_hist = data_hists2d_.at( NFT::kOnBNB ).get(); // BNB data
  // ROOT histograms use one-based bin indices, so I correct for that here
  double err = d_hist->GetBinError( reco_bin_a + 1, reco_bin_b + 1 );

  bool orbb = bin_type_a == kOrdinaryRecoBinBkgd
    || bin_type_b == kOrdinaryRecoBinBkgd;

  if ( orbb ) {
    // Don't evaluate the data statistical uncertainty for the "ordinary
    // background" bins unless we're considering EXT events
    if ( !use_ext ) return 0.;
  }

  double err2 = err * err;
  return err2;
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

MeasuredEvents ConstrainedCalculator::get_measured_events() const
{
  const int two_times_ord_bins = 2*num_ordinary_reco_bins_;
  const auto& cv_univ = this->cv_universe();

  // First create vectors of the measured event counts and central-value
  // prediction in each of the sideband bins. Note here and elsewhere in this
  // function that ROOT histogram bin indices are one-based to allow for
  // underflow. The TMatrixD element indices, on the other hand, are zero-based.
  TMatrixD sideband_data( num_sideband_reco_bins_, 1 );
  TMatrixD sideband_mc_plus_ext( num_sideband_reco_bins_, 1 );

  TH1D* d_hist = data_hists_.at( NFT::kOnBNB ).get(); // BNB data
  TH1D* ext_hist = data_hists_.at( NFT::kExtBNB ).get(); // EXT data
  for ( int s = 0; s < num_sideband_reco_bins_; ++s ) {
    // Zero-based reco bin index
    int r = s + num_ordinary_reco_bins_;

    // Switch to using the one-based TH1D index when retrieving these values
    double bnb_events = d_hist->GetBinContent( r + 1 );
    double ext_events = ext_hist->GetBinContent( r + 1 );
    double cv_mc_events = cv_univ.hist_reco_->GetBinContent( r + 1 );

    sideband_data( s, 0 ) = bnb_events;
    sideband_mc_plus_ext( s, 0 ) = cv_mc_events + ext_events;
  }

  // Create the vector of measured event counts in the ordinary reco bins
  TMatrixD ordinary_data( num_ordinary_reco_bins_, 1 );
  for ( int r = 0; r < num_ordinary_reco_bins_; ++r ) {
    // Switch to using the one-based TH1D index when retrieving these values
    double bnb_events = d_hist->GetBinContent( r + 1 );

    ordinary_data( r, 0 ) = bnb_events;
  }

  // Now create the vector which stores the unconstrained prediction for both
  // signal+background and background-only bins
  TMatrixD cv_pred_vec( two_times_ord_bins, 1 );
  for ( int r = 0; r < two_times_ord_bins; ++r ) {
    // This will automatically handle the signal+background versus
    // background-only bin definitions correctly. Recall that this
    // function takes a zero-based index.
    double cv_mc_events = this->evaluate_observable( cv_univ, r );

    // Also get the EXT event count for the reco bin of interest
    ConstrainedCalculatorBinType dummy_bin_type;
    int reco_bin_index = this->get_reco_bin_and_type( r, dummy_bin_type );

    // We need to use a one-based bin index to retrieve this value from the TH1D
    double ext_events = ext_hist->GetBinContent( reco_bin_index + 1 );

    cv_pred_vec( r, 0 ) = cv_mc_events + ext_events;
  }

  // All we need now to apply the sideband constraint are submatrices of
  // the total covariance matrix. Build it and pull out the blocks that
  // we need.
  auto cov_map_ptr = this->get_covariances();
  auto tot_cov_mat = cov_map_ptr->at( "total" ).get_matrix();

  // Zero-based indices for the covariance matrix elements describing the
  // ordinary reco bins (ob) and sideband reco bins (sb)
  int first_ob_cm_idx = 0;
  int last_ob_cm_idx = two_times_ord_bins - 1;
  int first_sb_cm_idx = two_times_ord_bins;
  int last_sb_cm_idx = two_times_ord_bins + num_sideband_reco_bins_ - 1;

  // Covariance matrix block that describes the ordinary reco bins
  TMatrixD ordinary_cov_mat = tot_cov_mat->GetSub( first_ob_cm_idx,
    last_ob_cm_idx, first_ob_cm_idx, last_ob_cm_idx );

  // Block that describes the sideband reco bins
  TMatrixD sideband_cov_mat = tot_cov_mat->GetSub( first_sb_cm_idx,
    last_sb_cm_idx, first_sb_cm_idx, last_sb_cm_idx );

  // Block that describes correlations between the sideband and ordinary bins.
  // This version uses rows for sideband bins and columns for ordinary bins.
  TMatrixD s_o_cov_mat = tot_cov_mat->GetSub( first_sb_cm_idx, last_sb_cm_idx,
    first_ob_cm_idx, last_ob_cm_idx );

  // Invert the sideband covariance matrix in preparation for applying
  // the sideband constraint
  auto inverse_sideband_cov_mat = invert_matrix( sideband_cov_mat );

  // We're ready. Apply the sideband constraint to the prediction vector first.
  TMatrixD sideband_data_mc_diff( sideband_data,
    TMatrixD::EMatrixCreatorsOp2::kMinus, sideband_mc_plus_ext );
  TMatrixD temp1( *inverse_sideband_cov_mat,
    TMatrixD::EMatrixCreatorsOp2::kMult, sideband_data_mc_diff );
  TMatrixD add_to( s_o_cov_mat,
    TMatrixD::EMatrixCreatorsOp2::kTransposeMult, temp1 );

  TMatrixD constr_cv_pred_vec( cv_pred_vec,
    TMatrixD::EMatrixCreatorsOp2::kPlus, add_to );

  // Now get the corresponding updated covariance matrix
  TMatrixD temp2( *inverse_sideband_cov_mat,
    TMatrixD::EMatrixCreatorsOp2::kMult, s_o_cov_mat );
  TMatrixD subtract_from( s_o_cov_mat,
    TMatrixD::EMatrixCreatorsOp2::kTransposeMult, temp2 );

  TMatrixD constr_ordinary_cov_mat( ordinary_cov_mat,
    TMatrixD::EMatrixCreatorsOp2::kMinus, subtract_from );

  // Pull out the block of the constrained covariance matrix that describes
  // the uncertainty on signal+background only
  auto* sig_plus_bkgd_cov_mat = new TMatrixD(
    constr_ordinary_cov_mat.GetSub( 0, num_ordinary_reco_bins_ - 1,
      0, num_ordinary_reco_bins_ - 1 )
  );

  // Get the constrained background prediction column vector
  auto* reco_bkgd_vec = new TMatrixD(
    constr_cv_pred_vec.GetSub( num_ordinary_reco_bins_,
      two_times_ord_bins - 1, 0, 0 )
  );

  // Get the constrained signal+background prediction column vector
  TMatrixD* mc_plus_ext_vec = new TMatrixD(
    constr_cv_pred_vec.GetSub( 0, num_ordinary_reco_bins_ - 1, 0, 0 )
  );

  // Get the ordinary reco bin data column vector after subtracting the
  // constrained background prediction
  auto* reco_data_minus_bkgd = new TMatrixD( ordinary_data,
    TMatrixD::EMatrixCreatorsOp2::kMinus, *reco_bkgd_vec );

  MeasuredEvents result( reco_data_minus_bkgd, reco_bkgd_vec,
    mc_plus_ext_vec, sig_plus_bkgd_cov_mat );
  return result;
}
