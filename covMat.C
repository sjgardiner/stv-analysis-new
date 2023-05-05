// Makes covariance matrices for the STV analysis

// Standard library includes
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

// ROOT includes
#include "TH1D.h"
#include "TH2D.h"
#include "TParameter.h"

// STV analysis includes
#include "CovMatUtils.hh"
#include "FilePropertiesManager.hh"
#include "PlotUtils.hh"
#include "UniverseMaker.hh"

using NFT = NtupleFileType;

// Simple container for the name of a subdirectory of the main response matrix
// TDirectoryFile (dir) paired with the namecycle of a universe histogram
// (namecycle)
struct DirNamecycle {

  DirNamecycle( const std::string& d, const std::string n )
    : dir_( d ), namecycle_( n ) {}

  std::string dir_;
  std::string namecycle_;
};

// Detector systematic variations are currently only available for Run 3
constexpr int DETVAR_RUN = 3;

// Namecycle to use when prepraring universe specifiers for the
// detector variation samples (including the detVar CV)
const std::string DETVAR_NAMECYCLE = "unweighted_0";

// Central value universe for reweightable systematics
const std::string CV_WEIGHT_NAMECYCLE = "weight_TunedCentralValue_UBGenie_0";

// Helper function that produces labels to use when naming covariance
// matrix histograms for detector systematic uncertainties
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

// Stores "settings" for calling the make_cov_mat function for a single
// systematic variation of interest
struct SystInfo {
  SystInfo( const std::string& prefix, size_t count, bool frac, bool avg,
    bool is_flux ) : wgt_prefix_( prefix ), univ_count_( count ),
    needs_fractional_( frac ), average_over_universes_( avg ),
    is_flux_variation_( is_flux ) {}

  std::string wgt_prefix_;
  size_t univ_count_;
  bool needs_fractional_;
  bool average_over_universes_;
  bool is_flux_variation_;
};

const std::map< std::string, SystInfo > SYSTEMATICS_TO_USE {

  { "flux", {"weight_flux_all_", 1000u, false, true, true} },

  { "reint", {"weight_reint_all_", 1000u, false, true, false} },

  { "xsec_multi", {"weight_All_UBGenie_", 500u, false, true, false} },

  // TODO: Double-check the unisims (I wrote this in a hurry)
  { "xsec_AxFFCCQEshape", {"weight_AxFFCCQEshape_UBGenie_", 2u,
    false, false, false} },

  { "xsec_DecayAngMEC", {"weight_DecayAngMEC_UBGenie_", 2u,
    false, false, false} },

  { "xsec_NormCCCOH", {"weight_NormCCCOH_UBGenie_", 2u,
    false, false, false} },

  { "xsec_NormNCCOH", {"weight_NormNCCOH_UBGenie_", 2u,
    false, false, false} },

  { "xsec_RPA_CCQE", {"weight_RPA_CCQE_UBGenie_", 2u,
    false, true, false} },

  { "xsec_ThetaDelta2NRad", {"weight_ThetaDelta2NRad_UBGenie_", 2u,
    false, false, false} },

  { "xsec_Theta_Delta2Npi", {"weight_Theta_Delta2Npi_UBGenie_", 2u,
    false, false, false} },

  { "xsec_VecFFCCQEshape", {"weight_VecFFCCQEshape_UBGenie_", 2u,
    false, false, false} },

  { "xsec_XSecShape_CCMEC", {"weight_XSecShape_CCMEC_UBGenie_", 2u,
    false, false, false} },

  { "xsec_xsr_scc_Fa3_SCC", {"weight_xsr_scc_Fa3_SCC_", 10u,
    false, true, false} },

  { "xsec_xsr_scc_Fv3_SCC", {"weight_xsr_scc_Fv3_SCC_", 10u,
    false, true, false} },

};

CovMatResults make_cov_mat( const std::string& cov_mat_name,
  TDirectoryFile& main_dir_file, const UniverseMaker& resp_mat,
  const DirNamecycle& cv_spec,
  const std::vector< DirNamecycle >& universe_spec, bool fractional,
  bool average_over_universes, bool is_flux_variation )
{
  // Get the total number of true and reco bins for later reference
  size_t num_true_bins = resp_mat.true_bins().size();
  size_t num_reco_bins = resp_mat.reco_bins().size();

  std::cout << "PROCESSING subdirectory " << cv_spec.dir_ << " for CV "
    << cv_spec.namecycle_ << '\n';

  std::string cv_true_hist_name = cv_spec.dir_
    + '/' + cv_spec.namecycle_ + "_true";

  std::unique_ptr< TH1D > cv_true_hist( dynamic_cast< TH1D* >(
    main_dir_file.Get( cv_true_hist_name.c_str() ) )
  );

  std::string cv_2d_hist_name = cv_spec.dir_ + '/' + cv_spec.namecycle_
    + "_2d";

  std::unique_ptr< TH2D > cv_2d_hist( dynamic_cast<TH2D*>(
    main_dir_file.Get( cv_2d_hist_name.c_str() ) )
  );

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

  // Prepare the covariance matrices for systematic variations on the signal
  // and on the backgrounds
  TH2D* covMat_signal = make_covariance_matrix_histogram(
    (cov_mat_name + "_signal").c_str(), "signal covariance", num_reco_bins );

  TH2D* covMat_bkgd = make_covariance_matrix_histogram(
    (cov_mat_name + "_bkgd").c_str(), "background covariance", num_reco_bins );

  // Loop over universes
  size_t num_universes = universe_spec.size();
  for ( size_t univ_idx = 0u; univ_idx < num_universes; ++univ_idx ) {

    const auto& univ_info = universe_spec.at( univ_idx );

    std::string univ_true_hist_name = univ_info.dir_
      + '/' + univ_info.namecycle_ + "_true";

    std::unique_ptr< TH1D > univ_true_hist( dynamic_cast<TH1D*>(
      main_dir_file.Get( univ_true_hist_name.c_str() ) )
    );

    std::string univ_2d_hist_name = univ_info.dir_
      + '/' + univ_info.namecycle_ + "_2d";

    std::unique_ptr< TH2D > univ_2d_hist( dynamic_cast<TH2D*>(
      main_dir_file.Get( univ_2d_hist_name.c_str() ) )
    );

    // Get the expected signal and background event counts in
    // each reco bin in the current universe.
    std::vector< double > univ_signal_events( num_reco_bins, 0. );
    std::vector< double > univ_bkgd_events( num_reco_bins, 0. );

    for ( size_t rb = 0u; rb < num_reco_bins; ++rb ) {
      // We need to sum the contributions of the various true bins,
      // so loop over them while checking whether each one is associated
      // with either signal or background
      for ( size_t tb = 0u; tb < num_true_bins; ++tb ) {
        const auto& tbin = resp_mat.true_bins().at( tb );

        if ( tbin.type_ == kSignalTrueBin ) {

          // Get the CV event count for the current true bin
          double denom_CV = cv_true_hist->GetBinContent( tb + 1 );

          // For the systematic variation universes, we want to assess
          // uncertainties on the signal only through the smearceptance
          // matrix. We therefore compute the smearceptance matrix element
          // here and then apply it to the CV expected event count in
          // each true bin.
          // NOTE: ROOT histogram bin numbers are one-based (bin zero is always
          // the underflow bin). Our bin indices therefore need to be offset by
          // +1 in all cases here.
          double numer = univ_2d_hist->GetBinContent( tb + 1, rb + 1 );
          double denom = univ_true_hist->GetBinContent( tb + 1 );

          // I plan to extract the flux-averaged cross sections in terms of the
          // *nominal* flux model (as opposed to the real flux). I therefore
          // vary the numerator of the smearceptance matrix for these while
          // keeping the denominator equal to the CV expectation under the
          // nominal flux model. This is the same strategy as is used in the
          // Wire-Cell CC inclusive analysis.
          if ( is_flux_variation ) {
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
          double expected_CV = smearcept * denom_CV;

          // Compute the expected signal events in the current reco bin
          // with the varied smearceptance matrix (and, for flux universes,
          // the varied integrated flux)
          univ_signal_events.at( rb ) += expected_CV;
        }
        else if ( tbin.type_ == kBackgroundTrueBin ) {
          // For background events, we can use the same procedure as
          // in the CV universe
          double background = univ_2d_hist->GetBinContent( tb + 1, rb + 1 );
          univ_bkgd_events.at( rb ) += background;
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

      double Sa_univ = univ_signal_events.at( a );
      double Ba_univ = univ_bkgd_events.at( a );

      for ( size_t b = 0u; b < num_reco_bins; ++b ) {

        double Sb_CV = cv_signal_events.at( b );
        double Bb_CV = cv_bkgd_events.at( b );

        double Sb_univ = univ_signal_events.at( b );
        double Bb_univ = univ_bkgd_events.at( b );

        double cov_signal = ( Sa_CV - Sa_univ ) * ( Sb_CV - Sb_univ );
        double cov_bkgd   = ( Ba_CV - Ba_univ ) * ( Bb_CV - Bb_univ );

        // Renormalize to get a fractional covariance matrix if the user
        // requested one
        if ( fractional ) {
          if ( Sa_CV <= 0. || Sb_CV <= 0. ) cov_signal = 0.;
          else cov_signal /= ( Sa_CV * Sb_CV );

          if ( Ba_CV <= 0. || Bb_CV <= 0. ) cov_bkgd = 0.;
          else cov_bkgd /= ( Ba_CV * Bb_CV );
        }

        // We cheat here by noting that the lower bound of each covariance
        // matrix TH2D bin is the bin index. Filling using the zero-based bin
        // indices and the covariance as the weight yields the desired behavior
        // (increment the existing element by the current covariance value) in
        // an easy-to-read (if slightly evil) way.
        covMat_signal->Fill( a, b, cov_signal );
        covMat_bkgd->Fill( a, b, cov_bkgd );
      } // reco bin index b
    } // reco bin index a

  } // universe

  // If requested, average the final covariance matrix elements over all
  // universes
  if ( average_over_universes ) {
    covMat_signal->Scale( 1. / num_universes );
    covMat_bkgd->Scale( 1. / num_universes );
  }

  // Prepare histograms of the CV expected counts for signal and background in
  // each reco bin
  TH1D* reco_cv_signal = new TH1D( ("reco_cv_signal_" + cov_mat_name).c_str(),
    "reco signal event counts; reco bin; events", num_reco_bins, 0.,
    num_reco_bins );
  reco_cv_signal->SetDirectory( nullptr );
  reco_cv_signal->SetStats( false );

  TH1D* reco_cv_bkgd = new TH1D( ("reco_cv_bkgd_" + cov_mat_name).c_str(),
    "reco background event counts; reco bin; events", num_reco_bins, 0.,
    num_reco_bins );
  reco_cv_bkgd->SetDirectory( nullptr );
  reco_cv_bkgd->SetStats( false );

  for ( size_t rb = 0u; rb < num_reco_bins; ++rb ) {
    double S_CV = cv_signal_events.at( rb );
    double B_CV = cv_bkgd_events.at( rb );

    reco_cv_signal->SetBinContent( rb + 1, S_CV );
    reco_cv_bkgd->SetBinContent( rb + 1, B_CV );

    // Set error bars on the CV histograms based on the diagonal covariance
    // matrix elements (i.e., the variances)
    double signal_var = covMat_signal->GetBinContent( rb + 1, rb + 1 );
    double bkgd_var = covMat_bkgd->GetBinContent( rb + 1, rb + 1 );

    // If we've generated a fractional covariance matrix, then correct
    // for this
    if ( fractional ) {
      signal_var *= std::pow( S_CV, 2 );
      bkgd_var *= std::pow( B_CV, 2 );
    }

    reco_cv_signal->SetBinError( rb + 1, std::sqrt(signal_var) );
    reco_cv_bkgd->SetBinError( rb + 1, std::sqrt(bkgd_var) );
  }

  CovMatResults result( covMat_signal, covMat_bkgd, reco_cv_signal,
    reco_cv_bkgd, fractional );

  return result;

  //KEY: TH1D     unweighted_0_reco;1
  //KEY: TH1D     unweighted_0_true;1
  //KEY: TH2D     unweighted_0_2d;1
}

// Overloaded version that takes settings from a SystInfo object
CovMatResults make_cov_mat( const std::string& cov_mat_name,
  TDirectoryFile& main_dir_file, const UniverseMaker& resp_mat,
  const DirNamecycle& cv_spec, const SystInfo& info )
{
  std::vector< DirNamecycle > universe_spec;
  for ( size_t u = 0u; u < info.univ_count_; ++u ) {
    universe_spec.emplace_back( cv_spec.dir_,
      info.wgt_prefix_ + std::to_string(u) );
  }

  return make_cov_mat( cov_mat_name, main_dir_file, resp_mat, cv_spec,
    universe_spec, info.needs_fractional_, info.average_over_universes_,
    info.is_flux_variation_ );
}

void covMat( const std::string& input_respmat_file_name,
  const std::string& config_file_name,
  const std::string& root_output_file = "" )
{
  // Map to hold all the covariance matrices that we'll need.
  // Outer keys are stv ntuple file names, inner keys are systematic
  // type labels ("detVar", etc.), and values are CovMatResults objects
  // containing the finished matrices.
  MatrixMap matrix_map;

  // The response matrices in each universe have already been made. However,
  // we'll re-create the object used to make them. It provides an easy
  // interface for interpreting the bin structure that was used.
  // TODO: consider refactoring this to have a "read-only" version. Also, be
  // careful. With the current approach, if you use the wrong configuration
  // file here, all of your covariance matrices could be invalid.
  UniverseMaker rmm( config_file_name );

  // Get access to the TDirectoryFile within the input ROOT file that holds the
  // response matrix histograms for the various analysis ntuples
  TFile in_tfile( input_respmat_file_name.c_str(), "read" );

  TDirectoryFile* respmat_dir = nullptr;
  in_tfile.GetObject( rmm.dir_name().c_str(), respmat_dir );

  if ( !respmat_dir ) {
    std::cout << "ERROR: Could not open the TDirectoryFile \""
      << rmm.dir_name() << "\" within the ROOT file "
      << input_respmat_file_name << '\n';
    return;
  }

  // If we've got to this point, then the TDirectoryFile containing the response
  // matrix histograms is available. Set it to be the active directory for ROOT.
  respmat_dir->cd();

  // Use the FilePropertiesManager to look up all of the detVar systematics
  // samples. At the moment, these are only available for Run 3. The fractional
  // uncertainty calculated for Run 3 is then applied to the other runs.
  // TODO: update this part when/if new detVar samples become available.
  const auto& fpm = FilePropertiesManager::Instance();

  // Get access to the files for Run 3
  const auto& detvar_run_file_map = fpm.ntuple_file_map().at( DETVAR_RUN );

  // Get the central-value detVar sample. I assume here (somewhat dangerously)
  // that there is only one nutple file in each Run 3 detector systematics
  // sample. This is currently true, but it could cause problems later.
  // TODO: revisit how you handle looking up the files to address this potential
  // problem.
  auto detvar_cv_iter = detvar_run_file_map.at( NFT::kDetVarMCCV ).cbegin();

  // Get the name of the detVarCV subdirectory of the main response matrix
  // TDirectoryFile
  std::string detvar_cv_subdirectory = ntuple_subfolder_from_file_name(
    *detvar_cv_iter );

  DirNamecycle detvar_cv_spec( detvar_cv_subdirectory, DETVAR_NAMECYCLE );

  // Specify the detVar samples (other than the CV) that we want to consider.
  // Note that a wireMod dE/dx sample is available, but we exclude it because it
  // is deprecated in favor of Recomb2. Using both together would be
  // double-counting.
  constexpr std::array< NtupleFileType, 9 > detVar_labels = {
    NFT::kDetVarMCLYatten, NFT::kDetVarMCLYdown, NFT::kDetVarMCLYrayl,
    NFT::kDetVarMCRecomb2, NFT::kDetVarMCSCE, NFT::kDetVarMCWMAngleXZ,
    NFT::kDetVarMCWMAngleYZ, NFT::kDetVarMCWMX, NFT::kDetVarMCWMYZ
  };

  std::vector< DirNamecycle > detvar_universes;

  for ( const auto& label : detVar_labels ) {

    // As I did for the CV above, I assume here that there is only one ntuple in
    // each of the Run 3 detector variation samples.
    // TODO: revisit this
    auto temp_iter = detvar_run_file_map.at( label ).cbegin();

    // Get the name of the subdirectory of the main response matrix
    // TDirectoryFile for the current detector variation sample
    std::string subdir = ntuple_subfolder_from_file_name(
      *temp_iter );

    detvar_universes.emplace_back( subdir, DETVAR_NAMECYCLE );
  }

  CovMatResults detVarResults = make_cov_mat( "detVar_total", *respmat_dir,
    rmm, detvar_cv_spec, detvar_universes, true, false, false );

  // The detector variation uncertainties will be applied globally using the
  // fractional covariance matrices for Run 3 calculated above. For all other
  // systematics, we go file-by-file across all MC samples with weights. The
  // ntuple file type labels we need to find these using the
  // FilePropertiesManager are listed here.
  constexpr std::array< NtupleFileType, 3 > rw_sample_labels = { NFT::kNumuMC,
    NFT::kIntrinsicNueMC, NFT::kDirtMC };

  // Here we loop over all central-value (i.e., non-detVar) MC samples in *all*
  // runs, not just Run 3 (like we did for the detVar samples above). For code
  // readability below, we'll precompute a vector of all of the ntuple file
  // names before the main loop for reweightable systematics.

  std::vector< std::string > cv_mc_file_names;

  for ( const auto& run_and_type_pair : fpm.ntuple_file_map() ) {

    const auto& type_map = run_and_type_pair.second;

    for ( const auto& type : rw_sample_labels ) {

      const auto& file_set = type_map.at( type );

      // Note that, since we iterate over the whole set of ntuple files here for
      // every run/type combination, the reweightable systematics automatically
      // accomodate an arbitrary number of input ntuple files.
      for ( const std::string& file_name : file_set ) {
        cv_mc_file_names.push_back( file_name );
      }

    } // ntuple file type loop

  } // run loop

  // Main loop over the central-value MC samples for calculating reweightable
  // systematic uncertainties
  for ( const auto& file_name : cv_mc_file_names ) {

    // Create the universe specification for the CV prediction first
    std::string file_subdir = ntuple_subfolder_from_file_name( file_name );

    DirNamecycle cv_spec( file_subdir, CV_WEIGHT_NAMECYCLE );

    for ( const auto& pair : SYSTEMATICS_TO_USE ) {
      std::string syst_name = pair.first;
      const auto& info = pair.second;

      CovMatResults temp_results = make_cov_mat( syst_name, *respmat_dir,
        rmm, cv_spec, info );

      if ( matrix_map.find(file_name) == matrix_map.end() ) {
        matrix_map[ file_name ] = std::map< std::string, CovMatResults >();
      }

      matrix_map[ file_name ][ syst_name ] = std::move( temp_results );

    } // loop over reweightable systematics

  } // loop over central-value MC input files

  // The map holding the covariance matrices is now complete for the high-stats
  // CV MC samples with one exception: we need to add in properly normalized
  // detector systematics. Do so below using the fractional covariance matrix
  // calculated above using the Run 3 detector variations.
  size_t num_reco_bins = 0u;
  for ( auto& pair : matrix_map ) {

    //std::string stv_file_name = pair.first;
    auto& syst_map = pair.second;

    // Retrieve the covariance matrix results for the cross-section multisim
    // systematics
    auto& xsec_multi_results = syst_map.at( "xsec_multi" );

    // Clone the CV predictions for event counts in each reco bin
    TH1D* reco_signal_cv = dynamic_cast< TH1D* >(
      xsec_multi_results.reco_signal_cv_->Clone()
    );
    reco_signal_cv->SetDirectory( nullptr );

    TH1D* reco_bkgd_cv = dynamic_cast< TH1D* >(
      xsec_multi_results.reco_bkgd_cv_->Clone()
    );
    reco_bkgd_cv->SetDirectory( nullptr );

    // Create empty covariance matrix histograms with the right number of bins
    num_reco_bins = reco_signal_cv->GetNbinsX();
    TH2D* signal_cov_mat = make_covariance_matrix_histogram( "detVar_total",
      "detector variations", num_reco_bins );
    TH2D* bkgd_cov_mat = make_covariance_matrix_histogram( "detVar_total",
      "detector variations", num_reco_bins );

    // Get access to the fractional covariance matrices for detector
    // systematics
    const auto& signal_frac_cm = detVarResults.signal_cov_mat_;
    const auto& bkgd_frac_cm = detVarResults.bkgd_cov_mat_;

    // Loop over each pair of reco bins and fill the new covariance matrices
    // with elements computed using the fractional covariance matrix and
    // the CV predictions in reco space.
    for ( size_t a = 0u; a < num_reco_bins; ++a ) {

      double Sa_CV = reco_signal_cv->GetBinContent( a + 1 );
      double Ba_CV = reco_bkgd_cv->GetBinContent( a + 1 );

      for ( size_t b = 0u; b < num_reco_bins; ++b ) {

        double Sb_CV = reco_signal_cv->GetBinContent( b + 1 );
        double Bb_CV = reco_bkgd_cv->GetBinContent( b + 1 );

        double frac_cov_signal = signal_frac_cm->GetBinContent( a + 1, b + 1 );
        double frac_cov_bkgd = bkgd_frac_cm->GetBinContent( a + 1, b + 1 );

        double cov_signal = Sa_CV * Sb_CV * frac_cov_signal;
        double cov_bkgd   = Ba_CV * Bb_CV * frac_cov_bkgd;

        signal_cov_mat->SetBinContent( a + 1, b + 1, cov_signal );
        bkgd_cov_mat->SetBinContent( a + 1, b + 1, cov_bkgd );
      } // reco bin b
    } // reco bin a

    // We're done preparing the new detector systematics covariance matrices.
    // Wrap them and the cloned CV predictions in a CovMatResults object
    // which will take ownership of the histograms.
    CovMatResults cmr( signal_cov_mat, bkgd_cov_mat, reco_signal_cv,
      reco_bkgd_cv, false );

    // Add the completed information to the map
    syst_map[ "detVar" ] = std::move( cmr );

  } // loop over high-stats CV MC STV ntuple files

  // Prepare a map of total beam exposures (POT) for each of the CV MC samples.
  // These will be used to normalize the results to the beam-on data below.
  // NOTE: we could have computed this in the loop immediately above this one.
  // However, separating this task probably makes the code more readable.

  // Keys are STV ntuple file names, values are beam exposures (POT)
  // TODO: consider changing the UniverseMaker class so that the
  // POT information is saved together with the response matrices themselves.
  std::map< std::string, float > pot_map;
  for ( const auto& pair : matrix_map ) {
    const std::string& stv_file_name = pair.first;
    TFile stv_file( stv_file_name.c_str(), "read" );

    TParameter<float>* summed_pot = nullptr;
    stv_file.GetObject( "summed_pot", summed_pot );

    pot_map[ stv_file_name ] = summed_pot->GetVal();
  }

  // We have what we need for the MC samples. Now accumulate beam-on and
  // beam-off data events in 1D histograms in reco space.
  TH1D* reco_bnb_hist = new TH1D( "reco_bnb_hist", "; reco bin; events",
    num_reco_bins, 0., num_reco_bins );
  reco_bnb_hist->Sumw2();

  TH1D* reco_ext_hist = new TH1D( "reco_ext_hist", "; reco bin; events",
    num_reco_bins, 0., num_reco_bins );
  reco_ext_hist->Sumw2();

  // Tallies needed for normalizing the final results
  double bnb_pot = 0.;
  double bnb_trigs = 0.;
  double ext_trigs = 0.;

  const auto& norm_map = fpm.data_norm_map();
  const auto& file_map = fpm.ntuple_file_map();

  std::map< int, double > bnb_pot_per_run;

  for ( const auto& pair : file_map ) {
    int run = pair.first;

    bnb_pot_per_run[ run ] = 0.;

    const auto& type_to_file_set = pair.second;

    const auto& bnb_set = type_to_file_set.at( NtupleFileType::kOnBNB );
    const auto& ext_set = type_to_file_set.at( NtupleFileType::kExtBNB );

    for ( const std::string& bnb_file : bnb_set ) {

      // Look up the correct TDirectoryFile subdirectory name to use for the
      // current beam-on data ntuple
      std::string bnb_subdir = ntuple_subfolder_from_file_name( bnb_file );

      // Retrieve the histogram containing the reco event counts
      std::string bnb_hist_name = bnb_subdir + "/unweighted_0_reco";

      TH1D* reco_evts = nullptr;
      respmat_dir->GetObject( bnb_hist_name.c_str(), reco_evts );

      // Add the current event counts to the total across all beam-on samples
      reco_bnb_hist->Add( reco_evts );

      const auto& trigs_and_pot = norm_map.at( bnb_file );
      bnb_pot += trigs_and_pot.pot_;
      bnb_trigs += trigs_and_pot.trigger_count_;

      bnb_pot_per_run.at( run ) += trigs_and_pot.pot_;

    } // loop over beam-on data files

    for ( const std::string& ext_file : ext_set ) {

      // Look up the correct TDirectoryFile subdirectory name to use for the
      // current extBNB (i.e., beam-off) data ntuple
      std::string ext_subdir = ntuple_subfolder_from_file_name( ext_file );

      // Retrieve the histogram containing the reco event counts
      std::string ext_hist_name = ext_subdir + "/unweighted_0_reco";

      TH1D* reco_evts = nullptr;
      respmat_dir->GetObject( ext_hist_name.c_str(), reco_evts );

      // Add the current event counts to the total across all extBNB samples
      reco_ext_hist->Add( reco_evts );

      const auto& trigs_and_pot = norm_map.at( ext_file );
      ext_trigs += trigs_and_pot.trigger_count_;

    } // loop over extBNB (i.e., beam-off) data files

  } // loop over runs

  // We've accumulated all of the measurements. Now scale the EXT events to
  // BNB based on the ratio of the recorded triggers.
  reco_ext_hist->Scale( bnb_trigs / ext_trigs );

  // Build the final event counts and covariance matrix for the MC + EXT
  // prediction ("pred"). Start by adding in the ready-to-go EXT events.
  TH1D* reco_pred_hist = new TH1D( "reco_pred_hist", "; reco bin; events",
    num_reco_bins, 0., num_reco_bins );
  reco_pred_hist->Sumw2();

  TH2D* pred_cov_mat = make_covariance_matrix_histogram( "pred_cov",
    "; reco bin; reco bin; covariance", num_reco_bins );

  // Add in the EXT event counts
  reco_pred_hist->Add( reco_ext_hist );

  // For EXT events, the covariance matrix is just the statistical variances
  // along the diagonal.
  for ( size_t a = 0u; a < num_reco_bins; ++a ) {
    pred_cov_mat->SetBinContent( a + 1, a + 1,
      std::pow( reco_ext_hist->GetBinError(a + 1), 2 )
    );
  }

  // Before entering the loop below, create a 2D histogram to accumulate
  // the POT-scaled total number of events predicted by the central-value MC
  // model in each reco/true bin combination.
  size_t num_true_bins = rmm.true_bins().size();

  TH2D* total_mc_hist_cv_2D = new TH2D( "total_mc_hist_cv_2D",
    "; true bin number; reco bin number; events", num_true_bins, 0.,
    num_true_bins, num_reco_bins, 0., num_reco_bins );
  total_mc_hist_cv_2D->Sumw2();
  total_mc_hist_cv_2D->SetStats( false );

  // Also create a new map of CovMatResults indexed by systematic label.
  // We will fill these with the POT-scaled total covariance matrices
  // for each individual source of uncertainty. Iterate over the systematic
  // labels used for the first ntuple file in the matrix_map to ensure that
  // you don't miss any
  std::map< std::string, CovMatResults > syst_to_total_cov_mat;
  for ( const auto& label_pair : matrix_map.cbegin()->second ) {

    const std::string& label = label_pair.first;

    TH1D* temp_signal_hist = new TH1D( "total_mc_cv_signal",
      "; reco bin; events", num_reco_bins, 0., num_reco_bins );
    temp_signal_hist->SetStats( false );
    temp_signal_hist->SetDirectory( nullptr );

    TH1D* temp_bkgd_hist = new TH1D( "total_mc_cv_bkgd",
      "; reco bin; events", num_reco_bins, 0., num_reco_bins );
    temp_bkgd_hist->SetStats( false );
    temp_bkgd_hist->SetDirectory( nullptr );

    TH2D* temp_signal_cov = make_covariance_matrix_histogram(
      (label + "_total_signal").c_str(), "total covariance matrix",
      num_reco_bins );

    TH2D* temp_bkgd_cov = make_covariance_matrix_histogram(
      (label + "_total_bkgd").c_str(), "total covariance matrix",
      num_reco_bins );

    CovMatResults my_temp_results( temp_signal_cov,
      temp_bkgd_cov, temp_signal_hist, temp_bkgd_hist, false );

    syst_to_total_cov_mat[ label ] = std::move( my_temp_results );

  } // loop over systematic labels

  // All that remains is to sum the MC contributions while scaling to
  // the correct POT
  for ( auto& pair : matrix_map ) {

    std::string stv_file_name = pair.first;
    auto& syst_map = pair.second;

    double file_pot = pot_map.at( stv_file_name );

    // TODO: super hacky. Please change this. It can easily break.
    int current_run = 0;
    if ( stv_file_name.find("run1") != std::string::npos ) {
      current_run = 1;
    }
    else if ( stv_file_name.find("run2") != std::string::npos ) {
      current_run = 2;
    }
    else if ( stv_file_name.find("run3") != std::string::npos ) {
      current_run = 3;
    }

    double run_bnb_pot = bnb_pot_per_run.at( current_run );

    double pot_scale = run_bnb_pot / file_pot;

    // Avoid double-counting by only adding the CV reco space predictions once
    const auto& temp_results = syst_map.at( "xsec_multi" );
    reco_pred_hist->Add( temp_results.reco_signal_cv_.get(), pot_scale );
    reco_pred_hist->Add( temp_results.reco_bkgd_cv_.get(), pot_scale );

    // Add each of the relevant covariance matrices to the total. Note that
    // the POT scaling factor needs to be squared in this case.
    for ( const auto& syst : syst_map ) {

      const std::string& syst_label = syst.first;
      const auto& syst_results = syst.second;

      double pot_scale2 = std::pow( pot_scale, 2 );
      pred_cov_mat->Add( syst_results.signal_cov_mat_.get(), pot_scale2 );
      pred_cov_mat->Add( syst_results.bkgd_cov_mat_.get(), pot_scale2 );

      // Also add the current covariance matrices to the running totals
      // for the appropriate individual systematic category
      auto& category_results = syst_to_total_cov_mat.at( syst_label );

      category_results.signal_cov_mat_->Add(
        syst_results.signal_cov_mat_.get(), pot_scale2 );

      category_results.bkgd_cov_mat_->Add(
        syst_results.bkgd_cov_mat_.get(), pot_scale2 );

      // We increment the category-specific central-value predictions here
      // because we only see each systematic category once in the loop. We thus
      // avoid double-counting and can conveniently do everything in the same
      // loop.
      category_results.reco_signal_cv_->Add(
        temp_results.reco_signal_cv_.get(), pot_scale );

      category_results.reco_bkgd_cv_->Add(
        temp_results.reco_bkgd_cv_.get(), pot_scale );

    } // systematic categories

    // While we're at it, retrieve the 2D central-value MC prediction for the
    // current ntuple file. Add it to the total 2D CV result with the
    // appropriate POT scaling.

    std::string subdir_2d = ntuple_subfolder_from_file_name( stv_file_name );
    std::string cv_2d_name = subdir_2d + '/' + CV_WEIGHT_NAMECYCLE + "_2d";

    TH2D* temp_cv_2d_hist( dynamic_cast<TH2D*>(
      respmat_dir->Get( cv_2d_name.c_str() ) )
    );

    total_mc_hist_cv_2D->Add( temp_cv_2d_hist, pot_scale );

  } // loop over the matrix map


  // For the "total_mc" results, we can sum subcategories of the systematics
  // here as we like. Currently, this is done only to (1) combine the various
  // xsec unisims into a summary, and (2) add that summary to the multisims
  // in order to get a combined xsec covariance matrix.

  const std::vector< std::string > extra_syst_labels = {
    "xsec_unisim", "xsec_all"
  };

  for ( const auto& label : extra_syst_labels ) {

    // For the CV predictions, just clone the histograms from
    // one of the pre-existing sets of results
    auto& results_to_clone = syst_to_total_cov_mat.begin()->second;

    TH1D* temp_signal_hist = dynamic_cast< TH1D* >(
      results_to_clone.reco_signal_cv_->Clone(
        (label + "_reco_signal_cv").c_str() )
    );
    temp_signal_hist->SetDirectory( nullptr );

    TH1D* temp_bkgd_hist = dynamic_cast< TH1D* >(
      results_to_clone.reco_bkgd_cv_->Clone(
        (label + "_reco_bkgd_cv").c_str() )
    );
    temp_bkgd_hist->SetDirectory( nullptr );

    TH2D* temp_signal_cov = make_covariance_matrix_histogram(
      (label + "_signal_cov_mat").c_str(), "total covariance matrix",
      num_reco_bins );

    TH2D* temp_bkgd_cov = make_covariance_matrix_histogram(
      (label + "_bkgd_cov_mat").c_str(), "total covariance matrix",
      num_reco_bins );

    CovMatResults my_temp_results( temp_signal_cov,
      temp_bkgd_cov, temp_signal_hist, temp_bkgd_hist, false );

    syst_to_total_cov_mat[ label ] = std::move( my_temp_results );

  }

  // All right, add up the xsec unisim covariance matrices in this loop.
  const std::array< std::string, 11 > xsec_unisim_labels = {
    "xsec_AxFFCCQEshape", "xsec_DecayAngMEC", "xsec_NormCCCOH",
    "xsec_NormNCCOH", "xsec_RPA_CCQE", "xsec_ThetaDelta2NRad",
    "xsec_Theta_Delta2Npi", "xsec_VecFFCCQEshape", "xsec_XSecShape_CCMEC",
    "xsec_xsr_scc_Fa3_SCC", "xsec_xsr_scc_Fv3_SCC"
  };

  auto& unisim_results = syst_to_total_cov_mat.at( "xsec_unisim" );

  for ( const auto& label : xsec_unisim_labels ) {

    auto& temp_results = syst_to_total_cov_mat.at( label );
    unisim_results.signal_cov_mat_->Add( temp_results.signal_cov_mat_.get() );
    unisim_results.bkgd_cov_mat_->Add( temp_results.bkgd_cov_mat_.get() );

  }

  auto& xsec_all_results = syst_to_total_cov_mat.at( "xsec_all" );

  const std::array< std::string, 2 > temp_labels = { "xsec_unisim",
    "xsec_multi" };

  for ( const auto& label : temp_labels ) {

    auto& temp_results = syst_to_total_cov_mat.at( label );
    xsec_all_results.signal_cov_mat_->Add( temp_results.signal_cov_mat_.get() );
    xsec_all_results.bkgd_cov_mat_->Add( temp_results.bkgd_cov_mat_.get() );

  }

  // Now that we don't need to do any more iterations over analysis ntuple
  // files, go ahead and add the map of total covariances (summed over the POT
  // for all non-detVar MC samples) to the main matrix_map. We can get it
  // into the map efficiently via move semantics, but note that the
  // syst_to_total_cov_mat map should not be used afterwards.
  matrix_map[ "total_mc" ] = std::move( syst_to_total_cov_mat );

  // We should be done now. Set the uncertainties on the final MC prediction
  // to be equal to the diagonal elements of the total covariance matrix.
  for ( size_t a = 0u; a < num_reco_bins; ++a ) {
    double cov = pred_cov_mat->GetBinContent( a + 1, a + 1 );
    reco_pred_hist->SetBinError( a + 1, std::sqrt(std::max(0., cov)) );
  }

  TCanvas* c1 = new TCanvas;
  pred_cov_mat->Draw( "colz" );

  TCanvas* c2 = new TCanvas;
  reco_bnb_hist->SetLineColor( kBlack );
  reco_bnb_hist->SetLineWidth( 3 );
  reco_bnb_hist->SetStats( false );
  reco_bnb_hist->Draw( "e" );
  reco_pred_hist->Draw( "same hist e" );

  reco_ext_hist->SetLineColor( kRed );
  reco_ext_hist->Draw( "same hist e" );

  TLegend* lg = new TLegend( 0.7, 0.7, 0.9, 0.9 );

  //std::string legend_title = get_legend_title( bnb_pot );
  //lg->SetHeader( legend_title.c_str(), "C" );

  lg->AddEntry( reco_bnb_hist, "BNB data", "l" );
  lg->AddEntry( reco_pred_hist, "MC (stat+syst)", "l" );
  lg->AddEntry( reco_ext_hist, "EXT BNB (stat)", "l" );
  lg->Draw( "same" );

  std::cout << "DATA POT = " << bnb_pot << '\n';

  // If the user has requested to save the results, do so now
  if ( !root_output_file.empty() ) {
    TFile out_tfile( root_output_file.c_str(), "recreate" );

    save_matrix_map( matrix_map, out_tfile );

    out_tfile.cd();
    out_tfile.WriteObject( &pot_map, "pot_map" );

    reco_bnb_hist->Write();
    reco_ext_hist->Write();
    reco_pred_hist->Write();

    pred_cov_mat->Write();

    total_mc_hist_cv_2D->Write();

    TParameter<double> temp_pot_param( "bnb_pot", bnb_pot );
    temp_pot_param.Write();
  }

}

  //// Plot some stuff
  //auto& inner_map = matrix_map.rbegin()->second;
  //for ( auto& pair : inner_map ) {
  //  auto& testResults = pair.second;

  //  TH2D* signal_cm = testResults.signal_cov_mat_.release();
  //  TH1D* signal_reco = testResults.reco_signal_cv_.release();

  //  TH2D* bkgd_cm = testResults.bkgd_cov_mat_.release();
  //  TH1D* bkgd_reco = testResults.reco_bkgd_cv_.release();

  //  TCanvas* c1 = new TCanvas;
  //  signal_cm->Draw( "colz" );
  //  TCanvas* c2 = new TCanvas;
  //  signal_reco->Draw( "hist e" );
  //  bkgd_reco->Draw( "hist e same" );

  //  TCanvas* c3 = new TCanvas;
  //  bkgd_cm->Draw( "colz" );
  //}
