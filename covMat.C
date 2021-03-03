// Makes covariance matrices for the STV analysis

// Standard library includes
#include <iostream>
#include <memory>
#include <vector>

// STV analysis includes
#include "ResponseMatrixMaker.hh"

const std::string DETVAR_CV_NAME = "prodgenie_bnb_nu_overlay_DetVar_CV_reco2"
  "_v08_00_00_38_run3b_reco2_reco2";
const std::vector< std::string > DETVAR_SAMPLE_NAMES = {
  "prodgenie_bnb_nu_overlay_DetVar_LYAttenuation_v08_00_00_38_run3b_reco2_reco2",
  "prodgenie_bnb_nu_overlay_DetVar_LYDown_v08_00_00_37_v2_run3b_reco2_reco2",
  "prodgenie_bnb_nu_overlay_DetVar_LYRayleigh_v08_00_00_37_run3b_reco2_reco2",
  "prodgenie_bnb_nu_overlay_DetVar_Recomb2_reco2_v08_00_00_39_run3b_reco2_reco2",
  "prodgenie_bnb_nu_overlay_DetVar_SCE_reco2_v08_00_00_38_run3b_reco2_reco2",
  "prodgenie_bnb_nu_overlay_DetVar_WireModAngleXZ_v08_00_00_38_exe_run3b_reco2_reco2",
  "prodgenie_bnb_nu_overlay_DetVar_WireModAngleYZ_v08_00_00_38_exe_run3b_reco2_reco2",
  "prodgenie_bnb_nu_overlay_DetVar_wiremod_ScaleX_v08_00_00_38_run3b_reco2_reco2",
  "prodgenie_bnb_nu_overlay_DetVar_wiremod_ScaleYZ_v08_00_00_38_run3b_reco2_reco2"
};

std::string detvar_sample_to_label( const std::string& sample_name ) {
  // Extract just the "DetVar" part of the file name to use as a label for
  // the associated covariance matrix
  unsigned first = sample_name.find( "DetVar_" );
  unsigned last = sample_name.find( "_v0" );
  std::string label = sample_name.substr( first, last - first );

  // If we have a trailing "_reco2" after doing so, also get rid of it
  const std::string bad_suffix( "_reco2" );
  if ( label.find(bad_suffix) != std::string::npos ) {
    label.erase( label.length() - bad_suffix.length() );
  }
  return label;
}


std::string sample_to_respmat_file_name( const std::string& sample_name,
  const std::string& respmat_folder )
{
  std::string file_name = respmat_folder + "/respmat-stv-" + sample_name + ".root";
  return file_name;
}

// Note: the caller takes ownership. The use of bare pointers is slightly
// dangerous here.
TH2D* make_covariance_matrix_histogram( const std::string& hist_name,
  const std::string& title, size_t num_reco_bins )
{
  TH2D* hist = new TH2D( hist_name.c_str(),
    (title + "; reco bin; reco bin; covariance").c_str(),
    num_reco_bins, 0., num_reco_bins, num_reco_bins, 0., num_reco_bins );
  hist->SetDirectory( nullptr );
  hist->SetStats( false );
  return hist;
}

void covMat() {

  // The response matrices in each universe have already been made. However,
  // we'll re-create the object used to make them. It provides an easy
  // interface for interpreting the bin structure that was used.
  // TODO: consider refactoring this to have a "read-only" version. Also, be
  // careful. With the current approach, if you use the wrong configuration file
  // here, all of your covariance matrices could be invalid.
  ResponseMatrixMaker resp_mat( "myconfig.txt" );

  // Get the total number of true and reco bins for later reference
  size_t num_true_bins = resp_mat.true_bins().size();
  size_t num_reco_bins = resp_mat.reco_bins().size();

  const std::string respmat_folder = "/uboone/data/users/gardiner/"
    "old-ntuples-stv/resp-all-backup";

  // Load detVar_CV histograms
  std::string detVar_CV_file_name = sample_to_respmat_file_name( DETVAR_CV_NAME,
    respmat_folder );

  std::unique_ptr< TFile > cv_file( new TFile(detVar_CV_file_name.c_str(), "read") );

  std::unique_ptr< TH1D > cv_true_hist( dynamic_cast<TH1D*>(
    cv_file->Get("unweighted_0_true")) );

  std::unique_ptr< TH2D > cv_2d_hist( dynamic_cast<TH2D*>(
    cv_file->Get("unweighted_0_2d")) );

  // Prepare the total covariance matrices for detector variation systematics
  // on the signal and on the backgrounds
  TH2D* covMat_DetVar_total_signal = make_covariance_matrix_histogram(
    "covMat_DetVar_total_signal", "DetVar_total_signal", num_reco_bins );

  TH2D* covMat_DetVar_total_bkgd = make_covariance_matrix_histogram(
    "covMat_DetVar_total_bkgd", "DetVar_total_bkgd", num_reco_bins );

  // Create the covariance matrices for detector variation systematics
  for ( const auto& sample_name : DETVAR_SAMPLE_NAMES ) {

    // Load the input histograms for the current detector variation sample
    std::string sample_file_name = sample_to_respmat_file_name( sample_name,
      respmat_folder );

    std::unique_ptr< TFile > sample_file( new TFile(sample_file_name.c_str(),
      "read") );

    std::unique_ptr< TH1D > sample_true_hist( dynamic_cast<TH1D*>(
      sample_file->Get("unweighted_0_true")) );

    std::unique_ptr< TH2D > sample_2d_hist( dynamic_cast<TH2D*>(
      sample_file->Get("unweighted_0_2d")) );

    // Create new covariance matrices for signal and background
    std::string var_name = detvar_sample_to_label( sample_name );

    TH2D* covMat_sample_signal = make_covariance_matrix_histogram(
      "covMat_" + var_name + "_signal", var_name + "_signal", num_reco_bins );

    TH2D* covMat_sample_bkgd = make_covariance_matrix_histogram(
      "covMat_" + var_name + "_bkgd", var_name + "_bkgd", num_reco_bins );

    // Get the expected signal and background event counts in each
    // reco bin in the CV universe
    std::vector< double > cv_signal_events( num_reco_bins, 0. );
    std::vector< double > cv_bkgd_events( num_reco_bins, 0. );

    for ( size_t rb = 0u; rb < num_reco_bins; ++rb ) {
      // We need to sum the contributions of the various true bins,
      // so loop over them while checking whether each one is associated
      // with either signal or background
      for ( size_t tb = 0u; tb < num_true_bins; ++tb ) {
        const auto& tbin = resp_mat.true_bins().at( tb );

        // For the CV universe, we don't have to worry about
        // computing the smearceptance matrix element explicitly
        // because the denominator cancels out when multiplying by
        // the CV prediction in the current true bin. We can therefore
        // use a similar recipe in this loop for both the signal and
        // background predictions in each reco bin.
        if ( tbin.type_ == kSignalTrueBin ) {
          // Note that ROOT histogram bin numbers are one-based (bin zero is
          // always the underflow bin). Our zero-based bin indices therefore
          // need to be offset by +1 in all cases here.
          double signal = cv_2d_hist->GetBinContent( tb + 1, rb + 1 );
          cv_signal_events.at( rb ) += signal;
        }
        else if ( tbin.type_ == kBackgroundTrueBin ) {
          double background = cv_2d_hist->GetBinContent( tb + 1, rb + 1 );
          cv_bkgd_events.at( rb ) += background;
        }
      } // true bins
    } // reco bins


    // Now get the expected signal and background event counts in
    // each reco bin in each universe.
    // TODO: loop over universes
    std::vector< double > sample_signal_events( num_reco_bins, 0. );
    std::vector< double > sample_bkgd_events( num_reco_bins, 0. );

    for ( size_t rb = 0u; rb < num_reco_bins; ++rb ) {
      // We need to sum the contributions of the various true bins,
      // so loop over them while checking whether each one is associated
      // with either signal or background
      for ( size_t tb = 0u; tb < num_true_bins; ++tb ) {
        const auto& tbin = resp_mat.true_bins().at( tb );

        if ( tbin.type_ == kSignalTrueBin ) {
          // For the systematic variation universes, we want to assess
          // uncertainties on the signal only through the smearceptance
          // matrix. We therefore compute the smearceptance matrix element
          // here and then apply it to the CV expected event count in
          // each true bin.
          // NOTE: ROOT histogram bin numbers are one-based (bin zero is always
          // the underflow bin). Our bin indices therefore need to be offset by
          // +1 in all cases here.
          double numer = sample_2d_hist->GetBinContent( tb + 1, rb + 1 );
          double denom = sample_true_hist->GetBinContent( tb + 1 );
          // If the denominator is nonzero actually calculate the fraction.
          // Otherwise, just leave it zeroed out.
          // TODO: revisit this, think about MC statistical uncertainties
          // on the empty bins
          double smearcept = 0.;
          if ( denom > 0. ) smearcept = numer / denom;

          // Get the CV event count for the current true bin
          double denom_CV = cv_true_hist->GetBinContent( tb + 1 );

          // Compute the expected signal events in the current reco bin
          // with the varied smearceptance matrix
          sample_signal_events.at( rb ) += smearcept * denom_CV;
        }
        else if ( tbin.type_ == kBackgroundTrueBin ) {
          // For background events, we can use the same procedure as
          // in the CV universe
          double background = sample_2d_hist->GetBinContent( tb + 1, rb + 1 );
          sample_bkgd_events.at( rb ) += background;
        }
      } // true bins
    } // reco bins

    // We have all the needed ingredients to get the contribution of this
    // universe to the covariance matrices. Loop over each pair of reco bins
    // and fill the corresponding covariance matrix elements.
    // TODO: the covariance matrices are symmetric by definition. You can
    // therefore make this more efficient by calculating only the subset of
    // elements that you need.
    for ( size_t a = 0u; a < num_reco_bins; ++a ) {

      double Sa_CV = cv_signal_events.at( a );
      double Ba_CV = cv_bkgd_events.at( a );

      double Sa_sample = sample_signal_events.at( a );
      double Ba_sample = sample_bkgd_events.at( a );

      for ( size_t b = 0u; b < num_reco_bins; ++b ) {

        double Sb_CV = cv_signal_events.at( b );
        double Bb_CV = cv_bkgd_events.at( b );

        double Sb_sample = sample_signal_events.at( b );
        double Bb_sample = sample_bkgd_events.at( b );

        double cov_signal = ( Sa_CV - Sa_sample ) * ( Sb_CV - Sb_sample );
        double cov_bkgd   = ( Ba_CV - Ba_sample ) * ( Bb_CV - Bb_sample );

        // Renormalize to get the *fractional* covariance matrix
        // TODO: make this conditional
        if ( Sa_CV <= 0. || Sb_CV <= 0. ) cov_signal = 0.;
        else cov_signal /= ( Sa_CV * Sb_CV );

        if ( Ba_CV <= 0. || Bb_CV <= 0. ) cov_bkgd = 0.;
        else cov_bkgd /= ( Ba_CV * Bb_CV );

        // We cheat here by noting that the lower bound of each covariance
        // matrix TH2D bin is the bin index. Filling using the zero-based bin
        // indices and the covariance as the weight yields the desired behavior
        // (increment the existing element by the current covariance value) in
        // an easy-to-read (if slightly evil) way.
        covMat_sample_signal->Fill( a, b, cov_signal );
        covMat_sample_bkgd->Fill( a, b, cov_bkgd );
      } // reco bin index b
    } // reco bin index a

    //TCanvas* c1 = new TCanvas;
    //covMat_sample_signal->Draw( "colz" );
    //TCanvas* c2 = new TCanvas;
    //covMat_sample_bkgd->Draw( "colz" );

    covMat_DetVar_total_signal->Add( covMat_sample_signal );
    covMat_DetVar_total_bkgd->Add( covMat_sample_bkgd );

  } // detVar sample

  TCanvas* c1 = new TCanvas;
  covMat_DetVar_total_signal->Draw( "colz" );
  TCanvas* c2 = new TCanvas;
  covMat_DetVar_total_bkgd->Draw( "colz" );

  //KEY: TH1D     unweighted_0_reco;1
  //KEY: TH1D     unweighted_0_true;1
  //KEY: TH2D     unweighted_0_2d;1
}
