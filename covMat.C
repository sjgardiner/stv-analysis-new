// Makes covariance matrices for the STV analysis

// Standard library includes
#include <iostream>
#include <map>
#include <memory>
#include <vector>

// STV analysis includes
#include "FilePropertiesManager.hh"
#include "IntegratedFluxUniverseManager.hh"
#include "ResponseMatrixMaker.hh"

struct FileNamecycle {

  FileNamecycle( const std::string& f, const std::string n )
    : file_( f ), namecycle_( n ) {}

  std::string file_;
  std::string namecycle_;
};

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

struct CovMatResults {

  CovMatResults() {}

  CovMatResults( TH2D* signal_cm, TH2D* bkgd_cm, TH1D* rc_signal_cv,
    TH1D* rc_bkgd_cv, bool frac ) : signal_cov_mat_( signal_cm ),
    bkgd_cov_mat_( bkgd_cm ), reco_signal_cv_( rc_signal_cv ),
    reco_bkgd_cv_( rc_bkgd_cv ), fractional_( frac ) {}

  std::unique_ptr< TH2D > signal_cov_mat_;
  std::unique_ptr< TH2D > bkgd_cov_mat_;

  std::unique_ptr< TH1D > reco_signal_cv_;
  std::unique_ptr< TH1D > reco_bkgd_cv_;

  bool fractional_ = false;

};


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
};

std::string stv_file_to_respmat_file( const std::string& stv_file_name,
  const std::string& respmat_folder )
{
  // Get the basename of the input file using the trick described here:
  // https://stackoverflow.com/a/24386991
  std::string stv_basename = stv_file_name.substr(
    stv_file_name.find_last_of('/') + 1 );

  std::string resp_mat_file_name = respmat_folder + "/respmat-" + stv_basename;
  return resp_mat_file_name;
}

CovMatResults make_cov_mat( const std::string& cov_mat_name,
  const ResponseMatrixMaker& resp_mat, const FileNamecycle& cv_spec,
  const std::vector< FileNamecycle >& universe_spec,
  bool fractional, bool average_over_universes, bool is_flux_variation )
{
  // Get the total number of true and reco bins for later reference
  size_t num_true_bins = resp_mat.true_bins().size();
  size_t num_reco_bins = resp_mat.reco_bins().size();

  std::cout << "PROCESSING file " << cv_spec.file_ << " for "
    << cv_spec.namecycle_ << '\n';

  std::unique_ptr< TFile > cv_file( new TFile(cv_spec.file_.c_str(), "read") );

  std::unique_ptr< TH1D > cv_true_hist( dynamic_cast<TH1D*>(
    cv_file->Get((cv_spec.namecycle_ + "_true").c_str()) )
  );

  std::unique_ptr< TH2D > cv_2d_hist( dynamic_cast<TH2D*>(
    cv_file->Get((cv_spec.namecycle_ + "_2d").c_str()) )
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

  // Prepare a TFile smart pointer to use to access the histograms in each
  // universe
  std::unique_ptr< TFile > univ_file;

  // Loop over universes
  size_t num_universes = universe_spec.size();
  for ( size_t univ_idx = 0u; univ_idx < num_universes; ++univ_idx ) {
    const auto& univ_info = universe_spec.at( univ_idx );

    // Check whether we need to open a new file for this universe. This
    // can happen either if we haven't opened one previously or if the
    // current universe's file name differs from the previous one.
    if ( !univ_file || univ_file->GetName() != univ_info.file_ ) {
      univ_file.reset( new TFile(univ_info.file_.c_str(), "read") );
    }

    std::unique_ptr< TH1D > univ_true_hist( dynamic_cast<TH1D*>(
      univ_file->Get((univ_info.namecycle_ + "_true").c_str()) )
    );

    std::unique_ptr< TH2D > univ_2d_hist( dynamic_cast<TH2D*>(
      univ_file->Get((univ_info.namecycle_ + "_2d").c_str()) )
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
          // If the denominator is nonzero actually calculate the fraction.
          // Otherwise, just leave it zeroed out.
          // TODO: revisit this, think about MC statistical uncertainties
          // on the empty bins
          double smearcept = 0.;
          if ( denom > 0. ) smearcept = numer / denom;

          // Get the CV event count for the current true bin
          double denom_CV = cv_true_hist->GetBinContent( tb + 1 );

          // Compute the expected signal events in this universe
          // by multiplying the varied smearceptance matrix element
          // by the unaltered CV prediction in the current true bin.
          double expected_CV = smearcept * denom_CV;

          // If we're working with flux variations, then we also need
          // to account for the (correlated) changes to the integrated
          // flux. We can do this by multiplying by the scaling factor
          // given here, which is the ratio of the integrated numu flux
          // in the current universe to the CV numu flux.
          if ( is_flux_variation ) {
            const auto& ifum = IntegratedFluxUniverseManager::Instance();
            double scale_factor = ifum.flux_factor( univ_idx );
            expected_CV *= scale_factor;
          }

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
  const ResponseMatrixMaker& resp_mat, const FileNamecycle& cv_spec,
  const SystInfo& info )
{
  std::vector< FileNamecycle > universe_spec;
  for ( size_t u = 0u; u < info.univ_count_; ++u ) {
    universe_spec.emplace_back( cv_spec.file_,
      info.wgt_prefix_ + std::to_string(u) );
  }

  return make_cov_mat( cov_mat_name, resp_mat, cv_spec, universe_spec,
    info.needs_fractional_, info.average_over_universes_,
    info.is_flux_variation_ );
}

void covMat() {

  // Map to hold all the covariance matrices that we'll need.
  // Outer keys are stv ntuple file names, inner keys are systematic
  // type labels ("detVar", etc.), and values are CovMatResults objects
  // containing the finished matrices.
  std::map< std::string, std::map<std::string, CovMatResults> > matrix_map;

  // The response matrices in each universe have already been made. However,
  // we'll re-create the object used to make them. It provides an easy
  // interface for interpreting the bin structure that was used.
  // TODO: consider refactoring this to have a "read-only" version. Also, be
  // careful. With the current approach, if you use the wrong configuration
  // file here, all of your covariance matrices could be invalid.
  ResponseMatrixMaker rmm( "myconfig.txt" );

  const std::string respmat_folder = "/uboone/data/users/gardiner/"
    "old-ntuples-stv/resp-all-backup";

  std::string detvar_cv_file_name = sample_to_respmat_file_name(
    DETVAR_CV_NAME, respmat_folder );

  FileNamecycle detvar_cv_spec( detvar_cv_file_name, "unweighted_0" );

  std::vector< FileNamecycle > detvar_universes;
  for ( const auto& sample : DETVAR_SAMPLE_NAMES ) {
    std::string sample_file_name = sample_to_respmat_file_name( sample,
      respmat_folder );
    detvar_universes.emplace_back( sample_file_name, "unweighted_0" );
  }

  CovMatResults detVarResults = make_cov_mat( "detVar_total", rmm,
    detvar_cv_spec, detvar_universes, true, false, false );

  // CV for reweightable systematics
  const std::string cv_weight_name = "weight_TunedCentralValue_UBGenie_0";

  // The detector variation uncertainties will be applied globally using
  // the fractional covariance matrices calculated above. For all other
  // systematics, we go file-by-file across all MC samples with weights.
  std::ifstream list_file( "foo3" );
  std::string file_name;
  while ( std::getline(list_file, file_name, '\n') ) {

    std::string respmat_file = stv_file_to_respmat_file(
      file_name, respmat_folder );

    FileNamecycle cv_spec( respmat_file, cv_weight_name );

    for ( const auto& pair : SYSTEMATICS_TO_USE ) {
      std::string syst_name = pair.first;
      const auto& info = pair.second;

      CovMatResults temp_results = make_cov_mat( syst_name, rmm,
        cv_spec, info );

      if ( matrix_map.find(file_name) == matrix_map.end() ) {
        matrix_map[ file_name ] = std::map< std::string, CovMatResults >();
      }
      matrix_map[ file_name ][ syst_name ] = std::move( temp_results );
    }

  } // input files for reweightable systematics

  // The map holding the covariance matrices is now complete for the high-stats
  // CV MC samples with one exception: we need to add in properly normalized
  // detector systematics. Do so below using the fractional covariance matrix
  // calculated above using the Run 3 detector variations.
  size_t num_reco_bins = 0u;
  for ( auto& pair : matrix_map ) {
    std::string stv_file_name = pair.first;
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
  std::map< std::string, float > pot_map;
  for ( const auto& pair : matrix_map ) {
    std::string stv_file_name = pair.first;
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

  const auto& fpm = FilePropertiesManager::Instance();
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
      std::string rf = stv_file_to_respmat_file( bnb_file, respmat_folder );
      TFile temp_file( rf.c_str(), "read" );

      TH1D* reco_evts = dynamic_cast< TH1D* >(
        temp_file.Get( "unweighted_0_reco" )
      );

      reco_bnb_hist->Add( reco_evts );

      const auto& trigs_and_pot = norm_map.at( bnb_file );
      bnb_pot += trigs_and_pot.pot_;
      bnb_trigs += trigs_and_pot.trigger_count_;

      bnb_pot_per_run.at( run ) += trigs_and_pot.pot_;
    }

    for ( const std::string& ext_file : ext_set ) {
      std::string rf = stv_file_to_respmat_file( ext_file, respmat_folder );
      TFile temp_file( rf.c_str(), "read" );

      TH1D* reco_evts = dynamic_cast< TH1D* >(
        temp_file.Get( "unweighted_0_reco" )
      );

      reco_ext_hist->Add( reco_evts );

      const auto& trigs_and_pot = norm_map.at( ext_file );
      ext_trigs += trigs_and_pot.trigger_count_;
    }

  }

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

  // All that remains is to sum the MC contributions while scaling to
  // the correct POT
  for ( auto& pair : matrix_map ) {

    std::string stv_file_name = pair.first;
    auto& syst_map = pair.second;

    double file_pot = pot_map.at( stv_file_name );

    // TODO: super hacky. Please change this.
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
      const auto& syst_results = syst.second;

      double pot_scale2 = std::pow( pot_scale, 2 );
      pred_cov_mat->Add( syst_results.signal_cov_mat_.get(), pot_scale2 );
      pred_cov_mat->Add( syst_results.bkgd_cov_mat_.get(), pot_scale2 );
    }
  }

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
  lg->AddEntry( reco_bnb_hist, "BNB data", "l" );
  lg->AddEntry( reco_pred_hist, "MC (stat+syst)", "l" );
  lg->AddEntry( reco_ext_hist, "EXT BNB (stat)", "l" );
  lg->Draw( "same" );

  std::cout << "DATA POT = " << bnb_pot << '\n';
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

// std::string var_name = detvar_sample_to_label( sample_name );
