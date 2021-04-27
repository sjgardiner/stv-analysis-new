// ROOT includes
#include "TFile.h"

// STV analysis includes
#include "FilePropertiesManager.hh"
#include "SystematicsCalculator.hh"

using NFT = NtupleFileType;

// Helper function for testing the start of a std::string. Taken from
// https://stackoverflow.com/a/40441240/4081973.
bool beginsWith( const std::string& str, const std::string& prefix ) {
  int result = str.rfind( prefix, 0 );
  if ( result == 0 ) return true;
  return false;
}

struct CovMatResults {

  CovMatResults() {}

  CovMatResults( TH2D* signal_cm, TH2D* bkgd_cm, TH2D* total_cm,
    bool fractional = false )
    : signal_cov_mat_( signal_cm ), bkgd_cov_mat_( bkgd_cm ),
    total_cov_mat_( total_cm ), fractional_( fractional ) {}

  std::unique_ptr< TH2D > signal_cov_mat_;
  std::unique_ptr< TH2D > bkgd_cov_mat_;
  std::unique_ptr< TH2D > total_cov_mat_;

  bool fractional_ = false;

  // Helper function for operator+=
  void add_or_clone( std::unique_ptr<TH2D>& mine, TH2D* other ) {
    if ( !other ) return;
    if ( mine ) mine->Add( other );
    else {
      TH2D* temp_clone = dynamic_cast< TH2D* >(
        other->Clone( "temp_clone" )
      );
      temp_clone->SetDirectory( nullptr );
      temp_clone->SetStats( false );
      mine.reset( temp_clone );
    }
  }

  CovMatResults& operator+=( const CovMatResults& other ) {
    if ( this->fractional_ || other.fractional_ ) throw std::runtime_error(
      "UNIMPLEMENTED!" );

    add_or_clone( signal_cov_mat_, other.signal_cov_mat_.get() );
    add_or_clone( bkgd_cov_mat_, other.bkgd_cov_mat_.get() );
    add_or_clone( total_cov_mat_, other.total_cov_mat_.get() );

    return *this;
  }

};

// Stores "settings" for calling the make_cov_mat function for a single
// systematic variation of interest
struct SystInfo {
  SystInfo( const std::string& name, bool avg, bool is_flux )
    : wgt_name_( name ), average_over_universes_( avg ),
    is_flux_variation_( is_flux ) {}

  std::string wgt_name_;
  bool average_over_universes_;
  bool is_flux_variation_;
};

const std::map< std::string, SystInfo > SYSTEMATICS_TO_USE {

  { "flux", {"weight_flux_all", true, true} },

  { "reint", {"weight_reint_all", true, false} },

  { "xsec_multi", {"weight_All_UBGenie", true, false} },

  // TODO: Double-check the unisims (I wrote this in a hurry)
  { "xsec_AxFFCCQEshape", {"weight_AxFFCCQEshape_UBGenie", false, false} },

  { "xsec_DecayAngMEC", {"weight_DecayAngMEC_UBGenie", false, false} },

  { "xsec_NormCCCOH", {"weight_NormCCCOH_UBGenie", false, false} },

  { "xsec_NormNCCOH", {"weight_NormNCCOH_UBGenie", false, false} },

  { "xsec_RPA_CCQE", {"weight_RPA_CCQE_UBGenie", true, false} },

  { "xsec_ThetaDelta2NRad", {"weight_ThetaDelta2NRad_UBGenie", false, false} },

  { "xsec_Theta_Delta2Npi", {"weight_Theta_Delta2Npi_UBGenie", false, false} },

  { "xsec_VecFFCCQEshape", {"weight_VecFFCCQEshape_UBGenie", false, false} },

  { "xsec_XSecShape_CCMEC", {"weight_XSecShape_CCMEC_UBGenie", false, false} },

  { "xsec_xsr_scc_Fa3_SCC", {"weight_xsr_scc_Fa3_SCC", true, false} },

  { "xsec_xsr_scc_Fv3_SCC", {"weight_xsr_scc_Fv3_SCC", true, false} },

};



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

template < class UniverseContainer >
  CovMatResults make_cov_mat( const std::string& cov_mat_name,
  const ResponseMatrixMaker& resp_mat, const Universe& cv_univ,
  const UniverseContainer& universes,
  bool average_over_universes, bool is_flux_variation,
  bool fractional = false )
{
  // Get the total number of true and reco bins for later reference
  size_t num_true_bins = resp_mat.true_bins().size();
  size_t num_reco_bins = resp_mat.reco_bins().size();

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
        double signal = cv_univ.hist_2d_->GetBinContent( tb + 1, rb + 1 );
        cv_signal_events.at( rb ) += signal;
      }
      else if ( tbin.type_ == kBackgroundTrueBin ) {
        double background = cv_univ.hist_2d_->GetBinContent( tb + 1, rb + 1 );
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

  TH2D* covMat_total = make_covariance_matrix_histogram(
    (cov_mat_name + "_total").c_str(), "total covariance", num_reco_bins );

  // Loop over universes
  size_t num_universes = universes.size();
  for ( const auto& univ : universes ) {

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
          double denom_CV = cv_univ.hist_true_->GetBinContent( tb + 1 );

          // For the systematic variation universes, we want to assess
          // uncertainties on the signal only through the smearceptance
          // matrix. We therefore compute the smearceptance matrix element
          // here and then apply it to the CV expected event count in
          // each true bin.
          // NOTE: ROOT histogram bin numbers are one-based (bin zero is always
          // the underflow bin). Our bin indices therefore need to be offset by
          // +1 in all cases here.
          double numer = univ->hist_2d_->GetBinContent( tb + 1, rb + 1 );
          double denom = univ->hist_true_->GetBinContent( tb + 1 );

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
          double background = univ->hist_2d_->GetBinContent( tb + 1, rb + 1 );
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

      double Ta_CV = Sa_CV + Ba_CV;

      double Sa_univ = univ_signal_events.at( a );
      double Ba_univ = univ_bkgd_events.at( a );

      double Ta_univ = Sa_univ + Ba_univ;

      for ( size_t b = 0u; b < num_reco_bins; ++b ) {

        double Sb_CV = cv_signal_events.at( b );
        double Bb_CV = cv_bkgd_events.at( b );

        double Tb_CV = Sb_CV + Bb_CV;

        double Sb_univ = univ_signal_events.at( b );
        double Bb_univ = univ_bkgd_events.at( b );

        double Tb_univ = Sb_univ + Bb_univ;

        double cov_signal = ( Sa_CV - Sa_univ ) * ( Sb_CV - Sb_univ );
        double cov_bkgd   = ( Ba_CV - Ba_univ ) * ( Bb_CV - Bb_univ );
        double cov_total  = ( Ta_CV - Ta_univ ) * ( Tb_CV - Tb_univ );

        // Renormalize to get a fractional covariance matrix if the user
        // requested one
        if ( fractional ) {
          if ( Sa_CV <= 0. || Sb_CV <= 0. ) cov_signal = 0.;
          else cov_signal /= ( Sa_CV * Sb_CV );

          if ( Ba_CV <= 0. || Bb_CV <= 0. ) cov_bkgd = 0.;
          else cov_bkgd /= ( Ba_CV * Bb_CV );

          if ( Ta_CV <= 0. || Tb_CV <= 0. ) cov_total = 0.;
          else cov_total /= ( Ta_CV * Tb_CV );
        }

        // We cheat here by noting that the lower bound of each covariance
        // matrix TH2D bin is the bin index. Filling using the zero-based bin
        // indices and the covariance as the weight yields the desired behavior
        // (increment the existing element by the current covariance value) in
        // an easy-to-read (if slightly evil) way.
        covMat_signal->Fill( a, b, cov_signal );
        covMat_bkgd->Fill( a, b, cov_bkgd );
        covMat_total->Fill( a, b, cov_total );
      } // reco bin index b
    } // reco bin index a

  } // universe

  // If requested, average the final covariance matrix elements over all
  // universes
  if ( average_over_universes ) {
    covMat_signal->Scale( 1. / num_universes );
    covMat_bkgd->Scale( 1. / num_universes );
    covMat_total->Scale( 1. / num_universes );
  }

  CovMatResults result( covMat_signal, covMat_bkgd, covMat_total, fractional );
  return result;
}

// Overloaded version that takes a single alternate universe wrapped in a
// std::unique_ptr
CovMatResults make_cov_mat( const std::string& cov_mat_name,
  const ResponseMatrixMaker& resp_mat, const Universe& cv_univ,
  const std::unique_ptr<Universe>& alt_univ,
  bool average_over_universes = false, bool is_flux_variation = false )
{
  // Evil, but we're going to temporarily move this around, so it's needed
  // TODO: do something far better than this crap!
  auto& my_univ = const_cast< std::unique_ptr<Universe>& >( alt_univ );
  std::vector< std::unique_ptr<Universe> > temp_univ_vec;

  temp_univ_vec.emplace_back( std::unique_ptr<Universe>(my_univ.release()) );

  auto result = make_cov_mat( cov_mat_name, resp_mat, cv_univ, temp_univ_vec,
    average_over_universes, is_flux_variation );

  // Put back the pointer
  my_univ.reset( temp_univ_vec.back().release() );

  return result;
}

// Overloaded version that takes settings from a SystInfo object and uses
// a SystematicsCalculator object to get the vector of universes
CovMatResults make_cov_mat( const std::string& cov_mat_name,
  const ResponseMatrixMaker& resp_mat, const SystInfo& info,
  const SystematicsCalculator& syst )
{
  const auto& univ_vec = syst.rw_universes_.at( info.wgt_name_ );
  return make_cov_mat( cov_mat_name, resp_mat, syst.cv_universe(),
    univ_vec, info.average_over_universes_, info.is_flux_variation_ );
}

// Overloaded version that takes a detVar NtupleFileType and uses a
// SystematicsCalculator object to get the CV and alternate universe
CovMatResults make_cov_mat( const ResponseMatrixMaker& resp_mat,
  const SystematicsCalculator& syst, NFT type )
{
  bool is_not_detVar = !ntuple_type_is_detVar( type );
  if ( is_not_detVar ) throw std::runtime_error( "Invalid NtupleFileType!" );

  const auto& detVar_cv_u = syst.detvar_universes_.at( NFT::kDetVarMCCV );
  const auto& detVar_alt_u = syst.detvar_universes_.at( type );

  const auto& fpm = FilePropertiesManager::Instance();
  std::string cov_mat_name = fpm.ntuple_type_to_string( type );

  return make_cov_mat( cov_mat_name, resp_mat, *detVar_cv_u,
    detVar_alt_u, false, false );
}

// Overloaded version that uses a single alternate universe scaled by a
// fully-correlated fractional uncertainty from the CV universe. This version
// is used to easily assess POT and "number of targets" normalization
// uncertainties.
CovMatResults make_cov_mat( const std::string& cov_mat_name,
  const ResponseMatrixMaker& resp_mat, const Universe& cv_univ,
  double fully_correlated_fractional_uncertainty )
{
  // TODO: implement signal and background contributions here!
  const double frac2 = std::pow( fully_correlated_fractional_uncertainty, 2 );

  size_t num_reco_bins = resp_mat.reco_bins().size();

  CovMatResults results;

  TH2D* signal_cm = make_covariance_matrix_histogram( "scaled_signal",
    "signal covariance", num_reco_bins );

  TH2D* bkgd_cm = make_covariance_matrix_histogram( "scaled_bkgd",
    "background covariance", num_reco_bins );

  TH2D* total_cm = make_covariance_matrix_histogram( "scaled_total",
    "total covariance", num_reco_bins );

  results.signal_cov_mat_.reset( signal_cm );
  results.bkgd_cov_mat_.reset( bkgd_cm );
  results.total_cov_mat_.reset( total_cm );

  for ( int a = 1; a <= num_reco_bins; ++a ) {

    double Ta_CV = cv_univ.hist_reco_->GetBinContent( a );

    for ( int b = 1; b <= num_reco_bins; ++b ) {

      double Tb_CV = cv_univ.hist_reco_->GetBinContent( b );

      double total_cov = Ta_CV * Tb_CV * frac2;

      //results.signal_cov_mat_->SetBinContent( a, b, signal_cov );
      //results.bkgd_cov_mat_->SetBinContent( a, b, bkgd_cov );
      results.total_cov_mat_->SetBinContent( a, b, total_cov );

    } // reco bin b

  } // reco bin a

  return results;
}

// Use the MC statistical uncertainties from the CV universe to make a diagonal
// covariance matrix.
CovMatResults make_MC_stats_cov_mat( const ResponseMatrixMaker& resp_mat,
  const Universe& cv_univ )
{
  // Get the bin structure so that we can distinguish between signal and
  // background when making the covariance matrices
  const auto& true_bins = resp_mat.true_bins();
  const auto& reco_bins = resp_mat.reco_bins();

  size_t num_true_bins = true_bins.size();
  size_t num_reco_bins = reco_bins.size();

  CovMatResults results;

  TH2D* signal_cm = make_covariance_matrix_histogram( "MCstats_signal",
    "signal covariance", num_reco_bins );

  TH2D* bkgd_cm = make_covariance_matrix_histogram( "MCstats_bkgd",
    "background covariance", num_reco_bins );

  TH2D* total_cm = make_covariance_matrix_histogram( "MCstats_total",
    "total covariance", num_reco_bins );

  results.signal_cov_mat_.reset( signal_cm );
  results.bkgd_cov_mat_.reset( bkgd_cm );
  results.total_cov_mat_.reset( total_cm );

  // TODO: optimize so that you can use the Sumw2 array entries instead.
  // This avoids needing to take many square roots only to square them again.
  for ( size_t rb = 0u; rb < num_reco_bins; ++rb ) {

    double signal_err2 = 0.;
    double bkgd_err2 = 0.;

    for ( size_t tb = 0u; tb < num_true_bins; ++tb ) {
      // To account for the underflow bin, we need to increment our
      // zero-based indices here by one
      double err2 = cv_univ.hist_2d_->GetBinError( tb + 1, rb + 1 );
      err2 *= err2;

      const auto& tbin = true_bins.at( tb );
      if ( tbin.type_ == kSignalTrueBin ) {
        signal_err2 += err2;
      }
      else {
        bkgd_err2 += err2;
      }

    } // true bins

    double total_err2 = signal_err2 + bkgd_err2;

    results.signal_cov_mat_->SetBinContent( rb + 1, rb + 1, signal_err2 );
    results.bkgd_cov_mat_->SetBinContent( rb + 1, rb + 1, bkgd_err2 );
    results.total_cov_mat_->SetBinContent( rb + 1, rb + 1, total_err2 );

  } // reco bins

  return results;
}

// Use the EXT data statistical uncertainties a diagonal covariance matrix.
CovMatResults make_ext_stats_cov_mat( const SystematicsCalculator& syst )
{
  const TH1D* ext_hist = syst.data_hists_.at( NFT::kExtBNB ).get();
  size_t num_reco_bins = ext_hist->GetNbinsX();

  CovMatResults results;

  TH2D* signal_cm = make_covariance_matrix_histogram( "EXTstats_signal",
    "signal covariance", num_reco_bins );

  TH2D* bkgd_cm = make_covariance_matrix_histogram( "EXTstats_bkgd",
    "background covariance", num_reco_bins );

  TH2D* total_cm = make_covariance_matrix_histogram( "EXTstats_total",
    "total covariance", num_reco_bins );

  results.signal_cov_mat_.reset( signal_cm );
  results.bkgd_cov_mat_.reset( bkgd_cm );
  results.total_cov_mat_.reset( total_cm );

  // By definition, EXT events are background. So we'll have zero signal
  // covariance everywhere. The total is just a copy of the background one.
  // Note the one-based bin numbering convention for TH1D.
  for ( int rb = 1; rb <= num_reco_bins; ++rb ) {
    double err2 = ext_hist->GetBinError( rb );
    err2 *= err2;

    results.bkgd_cov_mat_->SetBinContent( rb, rb, err2 );
    results.total_cov_mat_->SetBinContent( rb, rb, err2 );
  }

  return results;
}


void covMat( const std::string& input_respmat_file_name,
  const std::string& config_file_name,
  const std::string& root_output_file = "" )
{
  // Keys are covariance matrix types, values are CovMatResults containing the
  // corresponding matrices
  auto* matrix_map_ptr = new std::map< std::string, CovMatResults >;
  auto& matrix_map = *matrix_map_ptr;

  // The response matrices in each universe have already been made. However,
  // we'll re-create the object used to make them. It provides an easy
  // interface for interpreting the bin structure that was used.
  // TODO: Be careful. With the current approach, if you use the wrong
  // configuration file here, all of your covariance matrices could be invalid.
  ResponseMatrixMaker rmm( config_file_name );

  auto* syst_ptr = new SystematicsCalculator( input_respmat_file_name,
    rmm.dir_name() );
  auto& syst = *syst_ptr;

  const auto& fpm = FilePropertiesManager::Instance();

  for ( const auto& pair : syst.detvar_universes_ ) {
    NFT type = pair.first;
    if ( type == NFT::kDetVarMCCV ) continue;

    std::string type_str = fpm.ntuple_type_to_string( type );
    CovMatResults temp_results = make_cov_mat( rmm, syst, type );

    std::cout << type_str << '\n';
    matrix_map[ type_str ] = std::move( temp_results );
  }

  for ( const auto& pair : SYSTEMATICS_TO_USE ) {
    const std::string& cov_mat_name = pair.first;
    const auto& info = pair.second;

    std::cout << cov_mat_name << '\n';
    CovMatResults temp_results = make_cov_mat( cov_mat_name, rmm, info, syst );

    matrix_map[ cov_mat_name ] = std::move( temp_results );
  }

  // Compute a few "subtotal" covariance matrices. Store them in the map.
  CovMatResults detVar_total;
  CovMatResults xsec_unisim;
  for ( const auto& pair : matrix_map ) {
    const std::string& name = pair.first;
    const CovMatResults& cmr = pair.second;

    // Add the covariance matrices together for detVar and xsec universes
    if ( beginsWith(name, "xsec") && name != "xsec_multi" ) {
      xsec_unisim += cmr;
    }
    else if ( beginsWith(name, "detVar") ) {
      detVar_total += cmr;
    }

  }

  matrix_map[ "xsec_unisim" ] = std::move( xsec_unisim );
  matrix_map[ "detVar_total" ] = std::move( detVar_total );

  CovMatResults xsec_total;
  xsec_total += matrix_map.at( "xsec_multi" );
  xsec_total += matrix_map.at( "xsec_unisim" );

  matrix_map[ "xsec_total" ] = std::move( xsec_total );

  // POT fractional uncertainty is 2% (fully correlated)
  CovMatResults pot_covMat = make_cov_mat( "POT", rmm,
    syst.cv_universe(), 0.02 );

  matrix_map[ "POT" ] = std::move( pot_covMat );

  // Number-of-targets uncertainty is 1% (fully correlated)
  CovMatResults numTargets_covMat = make_cov_mat( "numTargets", rmm,
    syst.cv_universe(), 0.01 );

  matrix_map[ "numTargets" ] = std::move( numTargets_covMat );

  // MC stats uncertainty is based on CV universe histograms
  CovMatResults mc_stats_covMat = make_MC_stats_cov_mat(
    rmm, syst.cv_universe() );

  matrix_map[ "MCstats" ] = std::move( mc_stats_covMat );

  CovMatResults ext_stats_covMat = make_ext_stats_cov_mat( syst );

  matrix_map[ "EXTstats" ] = std::move( ext_stats_covMat );

////////////////////////////////////////////////

  // Build the final event counts and total covariance matrix for the MC + EXT
  // prediction ("pred").
  int num_reco_bins = rmm.reco_bins().size();
  TH1D* reco_pred_hist = new TH1D( "reco_pred_hist", "; reco bin; events",
    num_reco_bins, 0., num_reco_bins );
  reco_pred_hist->Sumw2();

  TH1D* reco_bnb_hist = syst.data_hists_.at( NFT::kOnBNB ).get();
  TH1D* reco_ext_hist = syst.data_hists_.at( NFT::kExtBNB ).get();

  // Add in the EXT event counts
  reco_pred_hist->Add( reco_ext_hist );

  // Add in the CV MC prediction
  reco_pred_hist->Add( syst.cv_universe().hist_reco_.get() );

  // Terms needed for the total covariance matrix
  const std::vector< std::string > total_cov_mat_keys = { "detVar_total",
    "flux", "reint", "xsec_total", "POT", "numTargets", "MCstats", "EXTstats"
  };

  // Keys are labels, values are fractional uncertainty histograms
  auto* fr_unc_hists = new std::map< std::string, TH1D* >();
  auto& frac_uncertainty_hists = *fr_unc_hists;

  CovMatResults total_covMat;
  int color = 1;
  for ( const auto& key : total_cov_mat_keys ) {

    const auto& temp_results = matrix_map.at( key );
    total_covMat += temp_results;

    TH1D* temp_hist = new TH1D( ("myfrac_temp_hist_" + key).c_str(),
      "; reco bin; fractional uncertainty", num_reco_bins, 0., num_reco_bins );

    temp_hist->SetStats( false );
    //temp_hist->SetDirectory( nullptr );

    for ( int rb = 1; rb <= num_reco_bins; ++rb ) {
      double err2 = temp_results.total_cov_mat_->GetBinContent( rb, rb );
      double err = std::sqrt( std::max(0., err2) );
      double cv = reco_pred_hist->GetBinContent( rb );
      if ( cv > 0. ) err /= cv;
      else err = 0.;

      temp_hist->SetBinContent( rb, err );

      frac_uncertainty_hists[ key ] = temp_hist;
    }

    if ( color <= 9 ) ++color;
    if ( color == 5 ) ++color;
    if ( color >= 10 ) color += 10;

    temp_hist->SetLineColor( color );
    temp_hist->SetLineWidth( 3 );
    temp_hist->GetYaxis()->SetRangeUser( 0., 1. );
  }

  matrix_map[ "total" ] = std::move( total_covMat );

  // We should be done now. Set the uncertainties on the final MC prediction
  // to be equal to the diagonal elements of the total covariance matrix.
  const auto& total_cov_matrix = matrix_map.at( "total" ).total_cov_mat_;
  for ( size_t a = 1u; a <= num_reco_bins; ++a ) {
    double cov = total_cov_matrix->GetBinContent( a, a );
    reco_pred_hist->SetBinError( a, std::sqrt(std::max(0., cov)) );
  }

  TH1D* total_frac_err_hist = new TH1D( "total_frac_err_hist",
    "; reco bin; events", num_reco_bins, 0., num_reco_bins );
  for ( size_t a = 1u; a <= num_reco_bins; ++a ) {
    double cv = reco_pred_hist->GetBinContent( a );
    double err = reco_pred_hist->GetBinError( a );
    if ( cv > 0. ) err /= cv;
    else err = 0.;
    total_frac_err_hist->SetBinContent( a, err );
  }

  TCanvas* c1 = new TCanvas;
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

  TCanvas* c2 = new TCanvas;
  TLegend* lg2 = new TLegend( 0.7, 0.7, 0.9, 0.9 );

  total_frac_err_hist->SetStats( false );
  total_frac_err_hist->GetYaxis()->SetRangeUser( 0.,
    total_frac_err_hist->GetMaximum() * 1.05 );
  total_frac_err_hist->SetLineColor( kBlack );
  total_frac_err_hist->SetLineWidth( 3 );
  total_frac_err_hist->Draw( "hist" );

  lg2->AddEntry( total_frac_err_hist, "total", "l" );

  for ( auto& pair : frac_uncertainty_hists ) {
    const auto& name = pair.first;
    TH1D* hist = pair.second;

    lg2->AddEntry( hist, name.c_str(), "l" );
    hist->Draw( "same hist" );
  }

  lg2->Draw( "same" );

}

void norm() {

  covMat( "/uboone/data/users/gardiner/respmat-myconfig_delta_pTx.root",
    "myconfig_mcc9_delta_pTx.txt" );

  //// Write the final histograms to the output TFile
  //TFile out_tfile( "myout.root", "recreate" );
  //sc.save_universes( out_tfile );

}
