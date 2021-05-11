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

struct CovMatrix {

  CovMatrix() {}

  CovMatrix( TH2D* cov_mat ) : cov_matrix_( cov_mat ) {}

  std::unique_ptr< TH2D > cov_matrix_;

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

  CovMatrix& operator+=( const CovMatrix& other ) {

    add_or_clone( cov_matrix_, other.cov_matrix_.get() );

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

template < class UniversePointerContainer >
  CovMatrix make_cov_mat( const std::string& cov_mat_name,
  const ResponseMatrixMaker& resp_mat, const Universe& cv_univ,
  const UniversePointerContainer& universes,
  bool average_over_universes, bool is_flux_variation )
{
  // Get the total number of true and reco bins for later reference
  size_t num_true_bins = resp_mat.true_bins().size();
  size_t num_reco_bins = resp_mat.reco_bins().size();

  // Get the expected event counts in each reco bin in the CV universe
  std::vector< double > cv_reco_events( num_reco_bins, 0. );

  for ( size_t rb = 0u; rb < num_reco_bins; ++rb ) {
    cv_reco_events.at( rb ) = cv_univ.hist_reco_->GetBinContent( rb + 1 );
  }

  // Prepare the covariance matrices for systematic variations on the signal
  // and on the backgrounds
  TH2D* cov_mat = make_covariance_matrix_histogram(
    cov_mat_name.c_str(), "covariance", num_reco_bins );

  // Loop over universes
  size_t num_universes = universes.size();
  for ( const auto& univ : universes ) {

    // Get the expected event counts in each reco bin in the current universe.
    std::vector< double > univ_reco_events( num_reco_bins, 0. );

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
          univ_reco_events.at( rb ) += expected_CV;
        }
        else if ( tbin.type_ == kBackgroundTrueBin ) {
          // For background events, we can use the same procedure as
          // in the CV universe
          double background = univ->hist_2d_->GetBinContent( tb + 1, rb + 1 );
          univ_reco_events.at( rb ) += background;
        }
      } // true bins
    } // reco bins

    // We have all the needed ingredients to get the contribution of this
    // universe to the covariance matrix. Loop over each pair of reco bins and
    // fill the corresponding covariance matrix elements.
    // TODO: the covariance matrix are symmetric by definition. You can
    // therefore make this more efficient by calculating only the subset of
    // elements that you need.
    for ( size_t a = 0u; a < num_reco_bins; ++a ) {

      double cv_a = cv_reco_events.at( a );
      double univ_a = univ_reco_events.at( a );

      for ( size_t b = 0u; b < num_reco_bins; ++b ) {

        double cv_b = cv_reco_events.at( b );
        double univ_b = univ_reco_events.at( b );

        double covariance  = ( cv_a - univ_a ) * ( cv_b - univ_b );

        // We cheat here by noting that the lower bound of each covariance
        // matrix TH2D bin is the bin index. Filling using the zero-based bin
        // indices and the covariance as the weight yields the desired behavior
        // (increment the existing element by the current covariance value) in
        // an easy-to-read (if slightly evil) way.
        cov_mat->Fill( a, b, covariance );
      } // reco bin index b
    } // reco bin index a

  } // universe

  // If requested, average the final covariance matrix elements over all
  // universes
  if ( average_over_universes ) {
    cov_mat->Scale( 1. / num_universes );
  }

  CovMatrix result( cov_mat );
  return result;
}

// Overloaded version that takes a single alternate universe wrapped in a
// std::unique_ptr
CovMatrix make_cov_mat( const std::string& cov_mat_name,
  const ResponseMatrixMaker& resp_mat, const Universe& cv_univ,
  const std::unique_ptr<Universe>& alt_univ,
  bool average_over_universes = false, bool is_flux_variation = false )
{
  std::vector< const Universe* > temp_univ_vec;

  temp_univ_vec.emplace_back( alt_univ.get() );

  auto result = make_cov_mat( cov_mat_name, resp_mat, cv_univ, temp_univ_vec,
    average_over_universes, is_flux_variation );

  return result;
}

// Overloaded version that takes settings from a SystInfo object and uses
// a SystematicsCalculator object to get the vector of universes
CovMatrix make_cov_mat( const std::string& cov_mat_name,
  const ResponseMatrixMaker& resp_mat, const SystInfo& info,
  const SystematicsCalculator& syst )
{
  const auto& univ_vec = syst.rw_universes_.at( info.wgt_name_ );
  return make_cov_mat( cov_mat_name, resp_mat, syst.cv_universe(),
    univ_vec, info.average_over_universes_, info.is_flux_variation_ );
}

// Overloaded version that takes a detVar NtupleFileType and uses a
// SystematicsCalculator object to get the CV and alternate universe
CovMatrix make_cov_mat( const ResponseMatrixMaker& resp_mat,
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
CovMatrix make_cov_mat( const std::string& cov_mat_name,
  const ResponseMatrixMaker& resp_mat, const Universe& cv_univ,
  double fully_correlated_fractional_uncertainty )
{
  // TODO: implement signal and background contributions here!
  const double frac2 = std::pow( fully_correlated_fractional_uncertainty, 2 );

  size_t num_reco_bins = resp_mat.reco_bins().size();

  TH2D* cov_mat = make_covariance_matrix_histogram( "cov_mat",
    "covariance", num_reco_bins );

  cov_mat->SetDirectory( nullptr );
  CovMatrix results( cov_mat );

  for ( int a = 1; a <= num_reco_bins; ++a ) {

    double cv_a = cv_univ.hist_reco_->GetBinContent( a );

    for ( int b = 1; b <= num_reco_bins; ++b ) {

      double cv_b = cv_univ.hist_reco_->GetBinContent( b );

      double covariance = cv_a * cv_b * frac2;

      results.cov_matrix_->SetBinContent( a, b, covariance );

    } // reco bin b

  } // reco bin a

  return results;
}

// Use the MC statistical uncertainties from the CV universe to make a diagonal
// covariance matrix.
CovMatrix make_MC_stats_cov_mat( const ResponseMatrixMaker& resp_mat,
  const Universe& cv_univ )
{
  // Get the bin structure so that we can distinguish between signal and
  // background when making the covariance matrices
  const auto& reco_bins = resp_mat.reco_bins();
  size_t num_reco_bins = reco_bins.size();

  TH2D* cov_mat = make_covariance_matrix_histogram( "MCstats",
    "covariance", num_reco_bins );

  CovMatrix results( cov_mat );

  // TODO: optimize so that you can use the Sumw2 array entries instead.
  // This avoids needing to take many square roots only to square them again.
  for ( size_t rb = 0u; rb < num_reco_bins; ++rb ) {

    // To account for the underflow bin, we need to increment the zero-based
    // index here by one
    double err2 = cv_univ.hist_reco_->GetBinError( rb + 1 );
    err2 *= err2;

    results.cov_matrix_->SetBinContent( rb + 1, rb + 1, err2 );

  } // reco bins

  return results;
}

// Use the EXT data statistical uncertainties a diagonal covariance matrix.
CovMatrix make_ext_stats_cov_mat( const SystematicsCalculator& syst )
{
  const TH1D* ext_hist = syst.data_hists_.at( NFT::kExtBNB ).get();
  size_t num_reco_bins = ext_hist->GetNbinsX();

  TH2D* cov_mat = make_covariance_matrix_histogram( "EXTstats",
    "covariance", num_reco_bins );

  CovMatrix results( cov_mat );

  // Note the one-based bin numbering convention for TH1D
  for ( int rb = 1; rb <= num_reco_bins; ++rb ) {
    double err2 = ext_hist->GetBinError( rb );
    err2 *= err2;

    results.cov_matrix_->SetBinContent( rb, rb, err2 );
  }

  return results;
}


void covMat( const std::string& input_respmat_file_name,
  const std::string& config_file_name,
  const std::string& root_output_file = "" )
{
  // Keys are covariance matrix types, values are CovMatrix containing the
  // corresponding matrices
  auto* matrix_map_ptr = new std::map< std::string, CovMatrix >;
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
    CovMatrix temp_results = make_cov_mat( rmm, syst, type );

    std::cout << type_str << '\n';
    matrix_map[ type_str ] = std::move( temp_results );
  }

  for ( const auto& pair : SYSTEMATICS_TO_USE ) {
    const std::string& cov_mat_name = pair.first;
    const auto& info = pair.second;

    std::cout << cov_mat_name << '\n';
    CovMatrix temp_results = make_cov_mat( cov_mat_name, rmm, info, syst );

    matrix_map[ cov_mat_name ] = std::move( temp_results );
  }

  // Compute a few "subtotal" covariance matrices. Store them in the map.
  CovMatrix detVar_total;
  CovMatrix xsec_unisim;
  for ( const auto& pair : matrix_map ) {
    const std::string& name = pair.first;
    const CovMatrix& cmr = pair.second;

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

  CovMatrix xsec_total;
  xsec_total += matrix_map.at( "xsec_multi" );
  xsec_total += matrix_map.at( "xsec_unisim" );

  matrix_map[ "xsec_total" ] = std::move( xsec_total );

  // POT fractional uncertainty is 2% (fully correlated)
  CovMatrix pot_covMat = make_cov_mat( "POT", rmm,
    syst.cv_universe(), 0.02 );

  matrix_map[ "POT" ] = std::move( pot_covMat );

  // Number-of-targets uncertainty is 1% (fully correlated)
  CovMatrix numTargets_covMat = make_cov_mat( "numTargets", rmm,
    syst.cv_universe(), 0.01 );

  matrix_map[ "numTargets" ] = std::move( numTargets_covMat );

  // MC stats uncertainty is based on CV universe histograms
  CovMatrix mc_stats_covMat = make_MC_stats_cov_mat(
    rmm, syst.cv_universe() );

  matrix_map[ "MCstats" ] = std::move( mc_stats_covMat );

  CovMatrix ext_stats_covMat = make_ext_stats_cov_mat( syst );

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

  // Also make an EXT + CV bkgd plot
  TH1D* reco_ext_plus_bkgd_hist
    = (TH1D*)reco_ext_hist->Clone( "ext_plus_bkgd" );

  int num_true_bins = rmm.true_bins().size();
  TH1D* reco_bkgd_hist = syst.cv_universe().hist_2d_->ProjectionY(
    "reco_bkgd", num_reco_bins + 1, num_true_bins );

  reco_ext_plus_bkgd_hist->Add( reco_bkgd_hist );

  // Terms needed for the total covariance matrix
  const std::vector< std::string > total_cov_mat_keys = { "detVar_total",
    "flux", "reint", "xsec_total", "POT", "numTargets", "MCstats", "EXTstats"
  };

  // Keys are labels, values are fractional uncertainty histograms
  auto* fr_unc_hists = new std::map< std::string, TH1D* >();
  auto& frac_uncertainty_hists = *fr_unc_hists;

  CovMatrix total_covMat;
  int color = 1;
  for ( const auto& key : total_cov_mat_keys ) {

    const auto& temp_results = matrix_map.at( key );
    total_covMat += temp_results;

    TH1D* temp_hist = new TH1D( ("myfrac_temp_hist_" + key).c_str(),
      "; reco bin; fractional uncertainty", num_reco_bins, 0., num_reco_bins );

    temp_hist->SetStats( false );
    //temp_hist->SetDirectory( nullptr );

    for ( int rb = 1; rb <= num_reco_bins; ++rb ) {
      double err2 = temp_results.cov_matrix_->GetBinContent( rb, rb );
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
  const auto& total_cov_matrix = matrix_map.at( "total" ).cov_matrix_;
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

  //reco_pred_hist->SetLineWidth( 3 );
  reco_pred_hist->Draw( "same hist e" );

  reco_bnb_hist->GetYaxis()->SetRangeUser( 0., 40e3 );
  reco_bnb_hist->Draw( "same e" );

  reco_ext_hist->SetLineColor( kRed );
  reco_ext_hist->SetLineWidth( 2 );
  reco_ext_hist->Draw( "same hist e" );

  reco_ext_plus_bkgd_hist->SetLineColor( 28 );
  reco_ext_plus_bkgd_hist->SetLineWidth( 2 );
  reco_ext_plus_bkgd_hist->Draw( "same hist e" );

  TLegend* lg = new TLegend( 0.7, 0.7, 0.9, 0.9 );

  //std::string legend_title = get_legend_title( bnb_pot );
  //lg->SetHeader( legend_title.c_str(), "C" );

  lg->AddEntry( reco_bnb_hist, "BNB data", "l" );
  lg->AddEntry( reco_pred_hist, "MC (stat+syst)", "l" );
  lg->AddEntry( reco_ext_hist, "EXT BNB (stat)", "l" );
  lg->AddEntry( reco_ext_plus_bkgd_hist, "Background (stat)", "l" );
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

    std::cout << name << " frac err in bin #1 = "
      << hist->GetBinContent( 1 )*100. << "%\n";
  }

  lg2->Draw( "same" );

  std::cout << "Total frac error in bin #1 = "
    << total_frac_err_hist->GetBinContent( 1 )*100. << "%\n";

}

void norm() {

  covMat( "/uboone/data/users/gardiner/respmat-myconfig_one_bin.root",
    "myconfig_mcc9_one_bin.txt" );

  //covMat( "/uboone/data/users/gardiner/respmat-test-muon2D.root",
  //  "myconfig_muon2D.txt" );

  //covMat( "/uboone/data/users/gardiner/respmat-myconfig_delta_phiT.root",
  //  "myconfig_mcc9_delta_phiT.txt" );

  //covMat( "/uboone/data/users/gardiner/respmat-myconfig_pn.root",
  //  "myconfig_mcc9_pn.txt" );

  //covMat( "/uboone/data/users/gardiner/respmat-myconfig_cth_p.root",
  //  "myconfig_mcc9_cth_p.txt" );

  //// Write the final histograms to the output TFile
  //TFile out_tfile( "myout.root", "recreate" );
  //sc.save_universes( out_tfile );

}
