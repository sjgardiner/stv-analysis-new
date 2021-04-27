// ROOT includes
#include "TFile.h"

// STV analysis includes
#include "FilePropertiesManager.hh"
#include "SystematicsCalculator.hh"

using NFT = NtupleFileType;

struct CovMatResults {

  CovMatResults() {}

  CovMatResults( TH2D* signal_cm, TH2D* bkgd_cm, TH2D* total_cm,
    bool fractional = false )
    : signal_cov_mat_( signal_cm ), bkgd_cov_mat_( bkgd_cm ),
    total_cov_mat_( total_cm ), fractional_( fractional ) {}

  std::unique_ptr< TH2D > signal_cov_mat_;
  std::unique_ptr< TH2D > bkgd_cov_mat_;
  std::unique_ptr< TH2D > total_cov_mat_;

  bool fractional_;

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

  SystematicsCalculator syst( input_respmat_file_name, rmm.dir_name() );

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

  for ( const auto& pair : matrix_map ) {
    TCanvas* c1 = new TCanvas;
    pair.second.total_cov_mat_->Draw( "colz" );
    break;
  }

  //// Build the final event counts and covariance matrix for the MC + EXT
  //// prediction ("pred"). Start by adding in the ready-to-go EXT events.
  //TH1D* reco_pred_hist = new TH1D( "reco_pred_hist", "; reco bin; events",
  //  num_reco_bins, 0., num_reco_bins );
  //reco_pred_hist->Sumw2();

  //TH2D* pred_cov_mat = make_covariance_matrix_histogram( "pred_cov",
  //  "; reco bin; reco bin; covariance", num_reco_bins );

  //// Add in the EXT event counts
  //reco_pred_hist->Add( reco_ext_hist );

  //// For EXT events, the covariance matrix is just the statistical variances
  //// along the diagonal.
  //for ( size_t a = 0u; a < num_reco_bins; ++a ) {
  //  pred_cov_mat->SetBinContent( a + 1, a + 1,
  //    std::pow( reco_ext_hist->GetBinError(a + 1), 2 )
  //  );
  //}

  //// Before entering the loop below, create a 2D histogram to accumulate
  //// the POT-scaled total number of events predicted by the central-value MC
  //// model in each reco/true bin combination.
  //size_t num_true_bins = rmm.true_bins().size();

  //TH2D* total_mc_hist_cv_2D = new TH2D( "total_mc_hist_cv_2D",
  //  "; true bin number; reco bin number; events", num_true_bins, 0.,
  //  num_true_bins, num_reco_bins, 0., num_reco_bins );
  //total_mc_hist_cv_2D->Sumw2();
  //total_mc_hist_cv_2D->SetStats( false );

  //// Also create a new map of CovMatResults indexed by systematic label.
  //// We will fill these with the POT-scaled total covariance matrices
  //// for each individual source of uncertainty. Iterate over the systematic
  //// labels used for the first ntuple file in the matrix_map to ensure that
  //// you don't miss any
  //std::map< std::string, CovMatResults > syst_to_total_cov_mat;
  //for ( const auto& label_pair : matrix_map.cbegin()->second ) {

  //  const std::string& label = label_pair.first;

  //  TH1D* temp_signal_hist = new TH1D( "total_mc_cv_signal",
  //    "; reco bin; events", num_reco_bins, 0., num_reco_bins );
  //  temp_signal_hist->SetStats( false );
  //  temp_signal_hist->SetDirectory( nullptr );

  //  TH1D* temp_bkgd_hist = new TH1D( "total_mc_cv_bkgd",
  //    "; reco bin; events", num_reco_bins, 0., num_reco_bins );
  //  temp_bkgd_hist->SetStats( false );
  //  temp_bkgd_hist->SetDirectory( nullptr );

  //  TH2D* temp_signal_cov = make_covariance_matrix_histogram(
  //    (label + "_total_signal").c_str(), "total covariance matrix",
  //    num_reco_bins );

  //  TH2D* temp_bkgd_cov = make_covariance_matrix_histogram(
  //    (label + "_total_bkgd").c_str(), "total covariance matrix",
  //    num_reco_bins );

  //  CovMatResults my_temp_results( temp_signal_cov,
  //    temp_bkgd_cov, temp_signal_hist, temp_bkgd_hist, false );

  //  syst_to_total_cov_mat[ label ] = std::move( my_temp_results );

  //} // loop over systematic labels

  //// All that remains is to sum the MC contributions while scaling to
  //// the correct POT
  //for ( auto& pair : matrix_map ) {

  //  std::string stv_file_name = pair.first;
  //  auto& syst_map = pair.second;

  //  double file_pot = pot_map.at( stv_file_name );

  //  // TODO: super hacky. Please change this. It can easily break.
  //  int current_run = 0;
  //  if ( stv_file_name.find("run1") != std::string::npos ) {
  //    current_run = 1;
  //  }
  //  else if ( stv_file_name.find("run2") != std::string::npos ) {
  //    current_run = 2;
  //  }
  //  else if ( stv_file_name.find("run3") != std::string::npos ) {
  //    current_run = 3;
  //  }

  //  double run_bnb_pot = bnb_pot_per_run.at( current_run );

  //  double pot_scale = run_bnb_pot / file_pot;

  //  // Avoid double-counting by only adding the CV reco space predictions once
  //  const auto& temp_results = syst_map.at( "xsec_multi" );
  //  reco_pred_hist->Add( temp_results.reco_signal_cv_.get(), pot_scale );
  //  reco_pred_hist->Add( temp_results.reco_bkgd_cv_.get(), pot_scale );

  //  // Add each of the relevant covariance matrices to the total. Note that
  //  // the POT scaling factor needs to be squared in this case.
  //  for ( const auto& syst : syst_map ) {

  //    const std::string& syst_label = syst.first;
  //    const auto& syst_results = syst.second;

  //    double pot_scale2 = std::pow( pot_scale, 2 );
  //    pred_cov_mat->Add( syst_results.signal_cov_mat_.get(), pot_scale2 );
  //    pred_cov_mat->Add( syst_results.bkgd_cov_mat_.get(), pot_scale2 );

  //    // Also add the current covariance matrices to the running totals
  //    // for the appropriate individual systematic category
  //    auto& category_results = syst_to_total_cov_mat.at( syst_label );

  //    category_results.signal_cov_mat_->Add(
  //      syst_results.signal_cov_mat_.get(), pot_scale2 );

  //    category_results.bkgd_cov_mat_->Add(
  //      syst_results.bkgd_cov_mat_.get(), pot_scale2 );

  //    // We increment the category-specific central-value predictions here
  //    // because we only see each systematic category once in the loop. We thus
  //    // avoid double-counting and can conveniently do everything in the same
  //    // loop.
  //    category_results.reco_signal_cv_->Add(
  //      temp_results.reco_signal_cv_.get(), pot_scale );

  //    category_results.reco_bkgd_cv_->Add(
  //      temp_results.reco_bkgd_cv_.get(), pot_scale );

  //  } // systematic categories

  //  // While we're at it, retrieve the 2D central-value MC prediction for the
  //  // current ntuple file. Add it to the total 2D CV result with the
  //  // appropriate POT scaling.

  //  std::string subdir_2d = ntuple_subfolder_from_file_name( stv_file_name );
  //  std::string cv_2d_name = subdir_2d + '/' + CV_WEIGHT_NAMECYCLE + "_2d";

  //  TH2D* temp_cv_2d_hist( dynamic_cast<TH2D*>(
  //    respmat_dir->Get( cv_2d_name.c_str() ) )
  //  );

  //  total_mc_hist_cv_2D->Add( temp_cv_2d_hist, pot_scale );

  //} // loop over the matrix map


  //// For the "total_mc" results, we can sum subcategories of the systematics
  //// here as we like. Currently, this is done only to (1) combine the various
  //// xsec unisims into a summary, and (2) add that summary to the multisims
  //// in order to get a combined xsec covariance matrix.

  //const std::vector< std::string > extra_syst_labels = {
  //  "xsec_unisim", "xsec_all"
  //};

  //for ( const auto& label : extra_syst_labels ) {

  //  // For the CV predictions, just clone the histograms from
  //  // one of the pre-existing sets of results
  //  auto& results_to_clone = syst_to_total_cov_mat.begin()->second;

  //  TH1D* temp_signal_hist = dynamic_cast< TH1D* >(
  //    results_to_clone.reco_signal_cv_->Clone(
  //      (label + "_reco_signal_cv").c_str() )
  //  );
  //  temp_signal_hist->SetDirectory( nullptr );

  //  TH1D* temp_bkgd_hist = dynamic_cast< TH1D* >(
  //    results_to_clone.reco_bkgd_cv_->Clone(
  //      (label + "_reco_bkgd_cv").c_str() )
  //  );
  //  temp_bkgd_hist->SetDirectory( nullptr );

  //  TH2D* temp_signal_cov = make_covariance_matrix_histogram(
  //    (label + "_signal_cov_mat").c_str(), "total covariance matrix",
  //    num_reco_bins );

  //  TH2D* temp_bkgd_cov = make_covariance_matrix_histogram(
  //    (label + "_bkgd_cov_mat").c_str(), "total covariance matrix",
  //    num_reco_bins );

  //  CovMatResults my_temp_results( temp_signal_cov,
  //    temp_bkgd_cov, temp_signal_hist, temp_bkgd_hist, false );

  //  syst_to_total_cov_mat[ label ] = std::move( my_temp_results );

  //}

  //// All right, add up the xsec unisim covariance matrices in this loop.
  //const std::array< std::string, 11 > xsec_unisim_labels = {
  //  "xsec_AxFFCCQEshape", "xsec_DecayAngMEC", "xsec_NormCCCOH",
  //  "xsec_NormNCCOH", "xsec_RPA_CCQE", "xsec_ThetaDelta2NRad",
  //  "xsec_Theta_Delta2Npi", "xsec_VecFFCCQEshape", "xsec_XSecShape_CCMEC",
  //  "xsec_xsr_scc_Fa3_SCC", "xsec_xsr_scc_Fv3_SCC"
  //};

  //auto& unisim_results = syst_to_total_cov_mat.at( "xsec_unisim" );

  //for ( const auto& label : xsec_unisim_labels ) {

  //  auto& temp_results = syst_to_total_cov_mat.at( label );
  //  unisim_results.signal_cov_mat_->Add( temp_results.signal_cov_mat_.get() );
  //  unisim_results.bkgd_cov_mat_->Add( temp_results.bkgd_cov_mat_.get() );

  //}

  //auto& xsec_all_results = syst_to_total_cov_mat.at( "xsec_all" );

  //const std::array< std::string, 2 > temp_labels = { "xsec_unisim",
  //  "xsec_multi" };

  //for ( const auto& label : temp_labels ) {

  //  auto& temp_results = syst_to_total_cov_mat.at( label );
  //  xsec_all_results.signal_cov_mat_->Add( temp_results.signal_cov_mat_.get() );
  //  xsec_all_results.bkgd_cov_mat_->Add( temp_results.bkgd_cov_mat_.get() );

  //}

  //// Now that we don't need to do any more iterations over analysis ntuple
  //// files, go ahead and add the map of total covariances (summed over the POT
  //// for all non-detVar MC samples) to the main matrix_map. We can get it
  //// into the map efficiently via move semantics, but note that the
  //// syst_to_total_cov_mat map should not be used afterwards.
  //matrix_map[ "total_mc" ] = std::move( syst_to_total_cov_mat );

  //// We should be done now. Set the uncertainties on the final MC prediction
  //// to be equal to the diagonal elements of the total covariance matrix.
  //for ( size_t a = 0u; a < num_reco_bins; ++a ) {
  //  double cov = pred_cov_mat->GetBinContent( a + 1, a + 1 );
  //  reco_pred_hist->SetBinError( a + 1, std::sqrt(std::max(0., cov)) );
  //}

  //TCanvas* c1 = new TCanvas;
  //pred_cov_mat->Draw( "colz" );

  //TCanvas* c2 = new TCanvas;
  //reco_bnb_hist->SetLineColor( kBlack );
  //reco_bnb_hist->SetLineWidth( 3 );
  //reco_bnb_hist->SetStats( false );
  //reco_bnb_hist->Draw( "e" );
  //reco_pred_hist->Draw( "same hist e" );

  //reco_ext_hist->SetLineColor( kRed );
  //reco_ext_hist->Draw( "same hist e" );

  //TLegend* lg = new TLegend( 0.7, 0.7, 0.9, 0.9 );

  ////std::string legend_title = get_legend_title( bnb_pot );
  ////lg->SetHeader( legend_title.c_str(), "C" );

  //lg->AddEntry( reco_bnb_hist, "BNB data", "l" );
  //lg->AddEntry( reco_pred_hist, "MC (stat+syst)", "l" );
  //lg->AddEntry( reco_ext_hist, "EXT BNB (stat)", "l" );
  //lg->Draw( "same" );

  //std::cout << "DATA POT = " << bnb_pot << '\n';

  //// If the user has requested to save the results, do so now
  //if ( !root_output_file.empty() ) {
  //  TFile out_tfile( root_output_file.c_str(), "recreate" );

  //  save_matrix_map( matrix_map, out_tfile );

  //  out_tfile.cd();
  //  out_tfile.WriteObject( &pot_map, "pot_map" );

  //  reco_bnb_hist->Write();
  //  reco_ext_hist->Write();
  //  reco_pred_hist->Write();

  //  pred_cov_mat->Write();

  //  total_mc_hist_cv_2D->Write();

  //  TParameter<double> temp_pot_param( "bnb_pot", bnb_pot );
  //  temp_pot_param.Write();
  //}

}

void norm() {

  covMat( "/uboone/data/users/gardiner/respmat-myconfig_delta_pTx.root",
    "myconfig_mcc9_delta_pTx.txt" );

  //// Write the final histograms to the output TFile
  //TFile out_tfile( "myout.root", "recreate" );
  //sc.save_universes( out_tfile );

}
