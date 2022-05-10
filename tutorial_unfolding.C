// Standard library includes
#include <iomanip>
#include <iostream>
#include <sstream>

// ROOT includes
#include "TCanvas.h"
#include "TLegend.h"

// STV analysis includes
#include "DAgostiniUnfolder.hh"
#include "FiducialVolume.hh"
#include "MCC9SystematicsCalculator.hh"
#include "SliceBinning.hh"
#include "SliceHistogram.hh"
#include "WienerSVDUnfolder.hh"

#define USE_FAKE_DATA ""

constexpr int NUM_DAGOSTINI_ITERATIONS = 6;

void tutorial_unfolding() {

  #ifdef USE_FAKE_DATA
    // Initialize the FilePropertiesManager and tell it to treat the NuWro
    // MC ntuples as if they were data
    auto& fpm = FilePropertiesManager::Instance();
    fpm.load_file_properties( "nuwro_file_properties.txt" );
  #endif

  // Do the systematics calculations in preparation for unfolding
  auto* mcc9 = new MCC9SystematicsCalculator(
    "/uboone/data/users/gardiner/tutorial_univmake_output.root",
    "systcalc.conf" );

  // Count the number of true bins in which signal events are stored
  int num_true_signal_bins = 0;
  for ( int t = 0; t < mcc9->true_bins_.size(); ++t ) {
    const auto& tbin = mcc9->true_bins_.at( t );
    if ( tbin.type_ == kSignalTrueBin ) ++num_true_signal_bins;
  }

  // Get the tuned GENIE CV prediction in each true bin (including the
  // background true bins)
  TH1D* genie_cv_truth = mcc9->cv_universe().hist_true_.get();
  int num_true_bins = genie_cv_truth->GetNbinsX();

  // While we're at it, clone the histogram and zero it out. We'll fill this
  // one with our unfolded result for easy comparison
  TH1D* unfolded_events = dynamic_cast< TH1D* >(
    genie_cv_truth->Clone("unfolded_events") );
  unfolded_events->Reset();

  // If present, then get the fake data event counts in each true bin
  // (including the background true bins). We hope to approximately reproduce
  // these event counts in the signal true bins via unfolding the fake data.
  const auto& fake_data_univ = mcc9->fake_data_universe();
  TH1D* fake_data_truth_hist = nullptr;

  bool using_fake_data = false;
  if ( fake_data_univ ) {
    using_fake_data = true;
    fake_data_truth_hist = fake_data_univ->hist_true_.get();
  }

  // Get a map from covariance matrix name (std::string, as defined in systcalc.conf) to
  // to the CovMatrix objects produced by the MCC9SystematicsCalculator class
  auto& cov_mat_map = *mcc9->get_covariances().release();

  // Get a TH2D containing the total covariance matrix for the expected number of events
  // in reconstructed space
  auto* cov_mat = cov_mat_map.at( "total" ).cov_matrix_.get();

  // Initialize an object derived from the Unfolder base class
  std::unique_ptr< Unfolder > unfolder (
    new DAgostiniUnfolder( NUM_DAGOSTINI_ITERATIONS )
    //new WienerSVDUnfolder( true,
    //WienerSVDUnfolder::RegularizationMatrixType::kSecondDeriv )
  );

  // Apply the unfolding algorithm to the measured data and propagate the
  // uncertainties analytically
  auto result = unfolder->unfold( *mcc9 );

  // Set the event counts in each bin of the histogram that displays the
  // unfolded result.
  for ( int t = 0; t < num_true_bins; ++t ) {
    double evts = 0.;
    double error = 0.;
    if ( t < num_true_signal_bins ) {
      evts = result.unfolded_signal_->operator()( t, 0 );
      error = std::sqrt( std::max(0., result.cov_matrix_->operator()( t, t )) );
    }

    // We need to use one-based indices while working with TH1D bins, while the
    // indices for TMatrixD are zero-based
    unfolded_events->SetBinContent( t + 1, evts );
    unfolded_events->SetBinError( t + 1, error );
  }

  unfolded_events->SetStats( false );
  unfolded_events->SetLineColor( kBlack );
  unfolded_events->SetLineWidth( 3 );
  // Show only the true bin counts in the signal region of interest
  unfolded_events->GetXaxis()->SetRangeUser( 0, num_true_signal_bins );

  genie_cv_truth->SetStats( false );
  genie_cv_truth->SetLineColor( kRed );
  genie_cv_truth->SetLineWidth( 3 );
  genie_cv_truth->SetLineStyle( 9 );

  unfolded_events->Draw( "e" );
  genie_cv_truth->Draw( "hist same" );

  if ( using_fake_data ) {
    fake_data_truth_hist->SetStats( false );
    fake_data_truth_hist->SetLineColor( kBlue );
    fake_data_truth_hist->SetLineWidth( 3 );
    fake_data_truth_hist->SetLineStyle( 2 );
    fake_data_truth_hist->Draw( "hist same" );
  }

  TLegend* lg = new TLegend( 0.15, 0.7, 0.3, 0.85 );
  lg->AddEntry( unfolded_events, "unfolded", "l" );
  lg->AddEntry( genie_cv_truth, "uB tune", "l" );
  if ( using_fake_data ) {
    lg->AddEntry( fake_data_truth_hist, "truth", "l" );
  }

  lg->Draw( "same" );

  // Plot slices of the unfolded result
  auto* sb = new SliceBinning( "tutorial_true_slice_config.txt" );

  // Get the factors needed to convert to cross-section units
  double total_pot = mcc9->total_bnb_data_pot_;
  double integ_flux = integrated_numu_flux_in_FV( total_pot );
  double num_Ar = num_Ar_targets_in_FV();

  // Retrieve the true-space expected event counts from NUISANCE output files
  // for each available generator model
  double conv_factor = num_Ar * integ_flux;

  // Add the fake data truth using a column vector of event counts
  TMatrixD fake_data_truth( num_true_signal_bins, 1 );
  if ( using_fake_data ) {
    for ( int b = 0; b < num_true_signal_bins; ++b ) {
      double true_evts = fake_data_truth_hist->GetBinContent( b + 1 );
      fake_data_truth( b, 0 ) = true_evts;
    }
  }

  // Add the GENIE CV model using a column vector of event counts
  TMatrixD genie_cv_truth_vec( num_true_signal_bins, 1 );
  for ( int b = 0; b < num_true_signal_bins; ++b ) {
    double true_evts = genie_cv_truth->GetBinContent( b + 1 );
    genie_cv_truth_vec( b, 0 ) = true_evts;
  }

  // If we're working with Wiener-SVD unfolding, then apply the additional
  // smearing matrix to all true-space distributions (except for the unfolded
  // data itself)
  WienerSVDUnfolder* wsvd_ptr = dynamic_cast< WienerSVDUnfolder* >(
    unfolder.get() );

  if ( wsvd_ptr ) {

    // Get access to the additional smearing matrix
    const TMatrixD& A_C = wsvd_ptr->additional_smearing_matrix();

    // Start with the fake data truth if present
    if ( using_fake_data ) {
      TMatrixD ac_truth( A_C, TMatrixD::kMult, fake_data_truth );
      fake_data_truth = ac_truth;
    }

    // Also transform the GENIE CV model
    TMatrixD genie_cv_temp( A_C, TMatrixD::kMult, genie_cv_truth_vec );
    genie_cv_truth_vec = genie_cv_temp;
  }

  for ( size_t sl_idx = 0u; sl_idx < sb->slices_.size(); ++sl_idx ) {

    const auto& slice = sb->slices_.at( sl_idx );

    // Make a histogram showing the unfolded true event counts in the current
    // slice
    SliceHistogram* slice_unf = SliceHistogram::make_slice_histogram(
      *result.unfolded_signal_, slice, result.cov_matrix_.get() );

    // Also use the GENIE CV model to do the same
    SliceHistogram* slice_cv = SliceHistogram::make_slice_histogram(
      genie_cv_truth_vec, slice, nullptr );

    // If present, also use the truth information from the fake data to do the
    // same
    SliceHistogram* slice_truth = nullptr;
    if ( using_fake_data ) {
      slice_truth = SliceHistogram::make_slice_histogram( fake_data_truth,
        slice, nullptr );
    }

    // Keys are legend labels, values are SliceHistogram objects containing
    // true-space predictions from the corresponding generator models
    auto slice_gen_map = std::map< std::string, SliceHistogram* >();

    slice_gen_map[ "unfolded data" ] = slice_unf;
    if ( using_fake_data ) {
      slice_gen_map[ "truth" ] = slice_truth;
    }
    slice_gen_map[ "MicroBooNE Tune" ] = slice_cv;

    int var_count = 0;
    std::string diff_xsec_denom;
    std::string diff_xsec_units_denom;
    double other_var_width = 1.;
    for ( const auto& ov_spec : slice.other_vars_ ) {
      double high = ov_spec.high_bin_edge_;
      double low = ov_spec.low_bin_edge_;
      const auto& var_spec = sb->slice_vars_.at( ov_spec.var_index_ );
      if ( high != low ) {
        ++var_count;
        other_var_width *= ( high - low );
        diff_xsec_denom += 'd' + var_spec.name_;
        const std::string& temp_units = var_spec.units_;
        if ( !temp_units.empty() ) {
          diff_xsec_units_denom += " / " + temp_units;
        }
      }
    }

    for ( size_t av_idx : slice.active_var_indices_ ) {
      const auto& var_spec = sb->slice_vars_.at( av_idx );
      const std::string& temp_name = var_spec.name_;
      if ( temp_name != "true bin number" ) {
        var_count += slice.active_var_indices_.size();
        diff_xsec_denom += 'd' + var_spec.name_;
        diff_xsec_units_denom += " / " + var_spec.units_;
      }
    }

    // NOTE: This currently assumes that each slice is a 1D histogram
    // TODO: revisit as needed
    int num_slice_bins = slice_unf->hist_->GetNbinsX();
    TMatrixD trans_mat( num_slice_bins, num_slice_bins );
    for ( int b = 0; b < num_slice_bins; ++b ) {
      double width = slice_unf->hist_->GetBinWidth( b + 1 );
      width *= other_var_width;
      trans_mat( b, b ) = 1e37 / ( width * integ_flux * num_Ar );
    }

    // Convert all slice histograms from true event counts to differential
    // cross-section units
    for ( auto& pair : slice_gen_map ) {
      auto* slice_h = pair.second;
      slice_h->transform( trans_mat );

      std::string slice_y_title;
      if ( var_count > 0 ) {
        slice_y_title += "d";
        if ( var_count > 1 ) slice_y_title += "^{" + std::to_string( var_count ) + "}";
        slice_y_title += "#sigma/" + diff_xsec_denom;
      }
      else {
        slice_y_title += "#sigma";
      }
      slice_y_title += " (10^{-37} cm^{2}" + diff_xsec_units_denom + " / Ar)";

      slice_h->hist_->GetYaxis()->SetTitle( slice_y_title.c_str() );
    }

    // Keys are generator legend labels, values are the results of a chi^2
    // test compared to the unfolded data (or, in the case of the unfolded
    // data, to the fake data truth)
    std::map< std::string, SliceHistogram::Chi2Result > chi2_map;
    std::cout << '\n';
    for ( const auto& pair : slice_gen_map ) {
      const auto& name = pair.first;
      const auto* slice_h = pair.second;

      // Decide what other slice histogram should be compared to this one,
      // then calculate chi^2
      SliceHistogram* other = nullptr;
      // We don't need to compare the unfolded data to itself, so just skip to
      // the next SliceHistogram and leave a dummy Chi2Result object in the map
      if ( name == "unfolded data" ) {
        chi2_map[ name ] = SliceHistogram::Chi2Result();
        continue;
      }
      // Compare all other distributions to the unfolded data
      else {
        other = slice_gen_map.at( "unfolded data" );
      }

      // Store the chi^2 results in the map
      const auto& chi2_result = chi2_map[ name ] = slice_h->get_chi2( *other );

      std::cout << "Slice " << sl_idx << ", " << name << ": \u03C7\u00b2 = "
        << chi2_result.chi2_ << '/' << chi2_result.num_bins_ << " bin";
      if ( chi2_result.num_bins_ > 1 ) std::cout << 's';
      std::cout << ", p-value = " << chi2_result.p_value_ << '\n';
    }

    TCanvas* c1 = new TCanvas;
    slice_unf->hist_->SetLineColor( kBlack );
    slice_unf->hist_->SetLineWidth( 3 );
    slice_unf->hist_->SetMarkerStyle( kFullCircle );
    slice_unf->hist_->SetMarkerSize( 0.8 );
    slice_unf->hist_->SetStats( false );

    double ymax = -DBL_MAX;
    slice_unf->hist_->Draw( "e" );
    for ( const auto& pair : slice_gen_map ) {
      const auto& name = pair.first;
      const auto* slice_h = pair.second;

      double max = slice_h->hist_->GetMaximum();
      if ( max > ymax ) ymax = max;

      if ( name == "unfolded data" || name == "truth"
        || name == "MicroBooNE Tune" ) continue;

      slice_h->hist_->SetLineWidth( 4 );

      slice_h->hist_->Draw( "hist same" );
    }

    slice_cv->hist_->SetStats( false );
    slice_cv->hist_->SetLineColor( kAzure - 7 );
    slice_cv->hist_->SetLineWidth( 5 );
    slice_cv->hist_->SetLineStyle( 5 );
    slice_cv->hist_->Draw( "hist same" );

    if ( using_fake_data ) {
      slice_truth->hist_->SetStats( false );
      slice_truth->hist_->SetLineColor( kOrange );
      slice_truth->hist_->SetLineWidth( 5 );
      slice_truth->hist_->Draw( "hist same" );
    }

    slice_unf->hist_->GetYaxis()->SetRangeUser( 0., ymax*1.07 );
    slice_unf->hist_->Draw( "e same" );

    TLegend* lg = new TLegend( 0.15, 0.6, 0.5, 0.88 );
    for ( const auto& pair : slice_gen_map ) {
      const auto& name = pair.first;
      const auto* slice_h = pair.second;

      std::string label = name;
      std::ostringstream oss;
      const auto& chi2_result = chi2_map.at( name );
      oss << std::setprecision( 3 ) << chi2_result.chi2_ << " / "
        << chi2_result.num_bins_ << " bin";
      if ( chi2_result.num_bins_ > 1 ) oss << 's';

      if ( name != "unfolded data" ) {
        label += ": #chi^{2} = " + oss.str();
      }

      lg->AddEntry( slice_h->hist_.get(), label.c_str(), "l" );
    }

    lg->Draw( "same" );

  } // slices

}

int main() {
  tutorial_unfolding();
  return 0;
}
