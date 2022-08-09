// Standard library includes
#include <iomanip>
#include <iostream>
#include <sstream>

// ROOT includes
#include "TCanvas.h"
#include "TLegend.h"

// STV analysis includes
#include "../DAgostiniUnfolder.hh"
#include "../FiducialVolume.hh"
#include "../MCC9SystematicsCalculator.hh"
#include "../SliceBinning.hh"
#include "../SliceHistogram.hh"
#include "../WienerSVDUnfolder.hh"

// Wiener-SVD includes
//#include "svd/include/WienerSVD.h"

// RooUnfold includes
//#include "RooUnfold/src/RooUnfoldResponse.h"
//#include "RooUnfold/src/RooUnfoldBayes.h"

//std::unique_ptr< RooUnfoldResponse > get_test_response(
//  const SystematicsCalculator& sc )
//{
//  // Make a TH2D containing the response (or "smearceptance") matrix elements
//  // in the format expected by RooUnfoldResponse.
//  // Note that RooUnfoldResponse uses the axis convention reco <-> x,
//  // true <-> y for defining the response matrix. In the TMatrixD
//  // representation that I use here, reco <-> row and true <-> column.
//  // Note also that the RooUnfoldResponse constructor clones the input
//  // histograms (rather than taking ownership). In light of that behavior, we
//  // will use temporary histograms that will go out of scope when this function
//  // exits.
//  auto smearcept = sc.get_cv_smearceptance_matrix();
//  auto true_signal = sc.get_cv_true_signal();
//
//  int num_ordinary_reco_bins = smearcept->GetNrows();
//  int num_true_signal_bins = smearcept->GetNcols();
//
//  TH2D resp( "resp", "response matrix", num_ordinary_reco_bins, 0.,
//    num_ordinary_reco_bins, num_true_signal_bins, 0., num_true_signal_bins );
//
//  for ( int r = 0; r < num_ordinary_reco_bins; ++r ) {
//    for ( int t = 0; t < num_true_signal_bins; ++t ) {
//      double elem = smearcept->operator()( r, t );
//      // RooUnfold expects a response matrix normalized in terms of event
//      // counts (as opposed to event fractions). We multiply here by the
//      // CV true signal prediction to convert from one to the other.
//      elem *= true_signal->operator()( t, 0 );
//      // Note that TMatrixD objects have zero-based indices while TH2D objects
//      // have one-based indices. We correct for that explicitly here.
//      resp.SetBinContent( r + 1, t + 1, elem );
//    }
//  }
//
//  // Now prepare TH1D objects for the prior on the true events and the measured
//  // (background-subtracted) reco events
//  TH1D sig( "sig", "true signal", num_true_signal_bins, 0.,
//    num_true_signal_bins );
//
//  for ( int t = 0; t < num_true_signal_bins; ++t ) {
//    double elem = true_signal->operator()( t, 0 );
//    // Note that TMatrixD objects have zero-based indices while TH1D objects
//    // have one-based indices. We correct for that explicitly here.
//    sig.SetBinContent( t + 1, elem );
//  }
//
//  // Get the background-subtracted measurement
//  auto meas = sc.get_measured_events();
//  const auto& reco_signal = meas.reco_signal_;
//
//  TH1D reco( "reco", "background-subtracted data", num_ordinary_reco_bins, 0.,
//    num_ordinary_reco_bins );
//
//  for ( int r = 0; r < num_ordinary_reco_bins; ++r ) {
//    double elem = reco_signal->operator()( r, 0 );
//    // Note that TMatrixD objects have zero-based indices while TH1D objects
//    // have one-based indices. We correct for that explicitly here.
//    reco.SetBinContent( r + 1, elem );
//  }
//
//  auto result = std::make_unique< RooUnfoldResponse >( &reco, &sig, &resp,
//    "MyTestResponse", "test response" );
//
//  return result;
//}

constexpr double BIG_DOUBLE = 1e300;

void multiply_1d_hist_by_matrix( TMatrixD* mat, TH1* hist ) {
  // Copy the histogram contents into a column vector
  int num_bins = mat->GetNcols();
  TMatrixD hist_mat( num_bins, 1 );
  for ( int r = 0; r < num_bins; ++r ) {
    hist_mat( r, 0 ) = hist->GetBinContent( r + 1 );
  }

  // Multiply the column vector by the input matrix
  // TODO: add error handling here related to matrix dimensions
  TMatrixD hist_mat_transformed( *mat, TMatrixD::EMatrixCreatorsOp2::kMult,
    hist_mat );

  // Update the input histogram contents with the new values
  for ( int r = 0; r < num_bins; ++r ) {
    double val = hist_mat_transformed( r, 0 );
    hist->SetBinContent( r + 1, val );
  }

}

//const std::string SAMPLE_NAME = "MicroBooNE_CC1MuNp_XSec_2D_PmuCosmu_nu_MC";
const std::string SAMPLE_NAME = "MicroBooNE_CC1MuNp_XSec_2D_PpCosp_nu_MC";

struct TruthFileInfo {
  TruthFileInfo() {}
  TruthFileInfo( const std::string& file_name, int color, int style )
    : file_name_( file_name ), color_( color ), style_( style ) {}

  std::string file_name_;
  int color_;
  int style_;
};

// Keys are generator names and versions, values are TruthFileInfo objects
// describing nuiscomp output files containing the differential cross-section
// predictions in each true bin
std::map< std::string, TruthFileInfo > truth_file_map = {
  { "GENIE 2.12.10",
    {"/uboone/app/users/gardiner/temp-gen/BuildEventGenerators/ubmc/comp-tki/comp-gv2.root", kBlue, 1 } },
  { "GENIE 3.0.6",
    {"/uboone/app/users/gardiner/temp-gen/BuildEventGenerators/ubmc/comp-tki/comp-gv3.root", kBlack, 2} },
  { "NEUT 5.4.0.1",
    {"/uboone/app/users/gardiner/temp-gen/BuildEventGenerators/ubmc/comp-tki/comp-neut.root", kRed, 9} },
  { "NuWro 19.02.1",
    {"/uboone/app/users/gardiner/temp-gen/BuildEventGenerators/ubmc/comp-tki/comp-nuwro.root", kViolet, 7} },
 //{ "GiBUU 2019",
 //  {"/uboone/app/users/gardiner/stv/mc/comp_cc1muNp_gibuu.root", kGreen, 10} },

};

struct SampleInfo {

  SampleInfo() {}

  SampleInfo( const std::string& resp, const std::string& width,
    const std::string& sb ) : respmat_file_( resp ), widths_file_( width ),
    sb_file_( sb ) {}

  // File containing the output of respmat.C with the same true binning
  // as was used in NUISANCE to make theoretical predictions
  std::string respmat_file_;

  // NUISANCE data file in which the bin widths (along each relevant axis) are
  // stored. These files are used by the preliminary NUISANCE implementation of
  // this cross-section measurement.
  std::string widths_file_;

  // SliceBinning configuration file (needs to be consistent with the bin
  // definitions used in the other two files)
  std::string sb_file_;
};

// Keys are NUISANCE sample names, values are SampleInfo objects that give
// the corresponding input file paths needed for this script
std::map< std::string, SampleInfo > sample_info_map = {

  // 2D muon measurement
  { "MicroBooNE_CC1MuNp_XSec_2D_PmuCosmu_nu_MC", { "/uboone/data/users"
    "/gardiner/respmat-test-new-Muon2D.root", "/uboone/app/users/gardiner"
    "/stv/mc/nuisance/data/MicroBooNE/mybins2Dmuon.txt",
    "mybins_mcc9_2D_muon.txt" } },

  // 2D proton measurement
  { "MicroBooNE_CC1MuNp_XSec_2D_PpCosp_nu_MC", { "/uboone/data/users"
    "/gardiner/respmat-test-new-Proton2D.root", "/uboone/app/users/gardiner"
    "/stv/mc/nuisance/data/MicroBooNE/mybins2Dproton.txt",
    "mybins_mcc9_2D_proton.txt" } },

};

// Multiplying by the conversion factor conv_factor should change a total cross
// section value (in cm^2 / Ar) into an expected number of true event counts.
// The value of the factor can be obtained by multiplying the number of Ar
// targets in the fiducial volume by the integrated numu flux (for the measured
// POT exposure) in the fiducial volume.
//
// This function returns a map in which the keys are legend labels for each
// generator model. The values are TMatrixD column vectors containing the
// expected true event counts (directly comparable to unfolded event counts
// from the data).
std::map< std::string, TMatrixD* > get_true_events_nuisance(
  const SampleInfo& info, double conv_factor )
{
  // We'll create one prediction per generator model in the truth_file_map
  std::map< std::string, TMatrixD* > truth_counts_map;
  for ( const auto& pair : truth_file_map ) {
    // First retrieve the raw NUISANCE histogram. It is expressed as a
    // total cross section with true bin number along the x-axis
    std::string generator_label = pair.first;
    const auto& file_info = pair.second;
    std::string nuisance_file = file_info.file_name_;
    TFile temp_in_file( nuisance_file.c_str(), "read" );
    TH1D* temp_hist = nullptr;
    temp_in_file.GetObject( SAMPLE_NAME.c_str(), temp_hist );

    // Set the associated directory to a null pointer. That way this histogram
    // won't be auto-deleted when the corresponding TFile object goes out of
    // scope
    temp_hist->SetDirectory( nullptr );

    // Set the style of the histogram to match the configuration in the
    // TruthFileInfo object
    temp_hist->SetLineColor( file_info.color_ );
    temp_hist->SetLineStyle( file_info.style_ );

    // Disable displaying the stats box
    temp_hist->SetStats( false );

    // Convert the content (and error) of each bin to an expected true event
    // count. Do this using the input conversion factor (integrated numu
    // flux * number of Ar targets in the fiducial volume) and the 2D bin
    // width.
    size_t num_bins = temp_hist->GetNbinsX();
    for ( size_t b = 0u; b < num_bins; ++b ) {
      double xsec = temp_hist->GetBinContent( b + 1 );
      double err = temp_hist->GetBinError( b + 1 );
      temp_hist->SetBinContent( b + 1, xsec * conv_factor );
      temp_hist->SetBinError( b + 1, err * conv_factor );
    }

    // Now change the TH1D into a TMatrixD column vector
    TMatrixD* temp_mat = new TMatrixD( num_bins, 1 );
    for ( size_t b = 0u; b < num_bins; ++b ) {
      temp_mat->operator()( b, 0 ) = temp_hist->GetBinContent( b + 1 );
    }

    // The conversion is done, so add the finished true event counts histogram
    // to the map
    truth_counts_map[ generator_label ] = temp_mat;
  }

  return truth_counts_map;
}

void test_unfolding() {

  // Initialize the FilePropertiesManager and tell it to treat the NuWro
  // MC ntuples as if they were data
  auto& fpm = FilePropertiesManager::Instance();
  fpm.load_file_properties( "../nuwro_file_properties.txt" );

  const auto& sample_info = sample_info_map.at( SAMPLE_NAME );
  //const auto& respmat_file_name = sample_info.respmat_file_;

  const std::string respmat_file_name( "/uboone/data/users/gardiner/"
    "afro-tki-more.root" );
    //"afro-new-both2D.root" );

  // Do the systematics calculations in preparation for unfolding
  auto* syst_ptr = new MCC9SystematicsCalculator( respmat_file_name, "../systcalc_unfold_fd.conf" );
  //auto* syst_ptr = new MCC9SystematicsCalculator( respmat_file_name, "../systcalc.conf" );
  auto& syst = *syst_ptr;

  // Get the tuned GENIE CV prediction in each true bin (including the
  // background true bins)
  TH1D* genie_cv_truth = syst.cv_universe().hist_true_.get();
  int num_true_bins = genie_cv_truth->GetNbinsX();

  // While we're at it, clone the histogram and zero it out. We'll fill this
  // one with our unfolded result for easy comparison
  TH1D* unfolded_events = dynamic_cast< TH1D* >(
    genie_cv_truth->Clone("unfolded_events") );
  unfolded_events->Reset();

  // If present, then get the fake data event counts in each true bin
  // (including the background true bins). We hope to approximately reproduce
  // these event counts in the signal true bins via unfolding the fake data.
  const auto& fake_data_univ = syst.fake_data_universe();
  TH1D* fake_data_truth_hist = nullptr;

  bool using_fake_data = false;
  if ( fake_data_univ ) {
    using_fake_data = true;
    fake_data_truth_hist = fake_data_univ->hist_true_.get();
  }

  int num_ordinary_reco_bins = 0;
  int num_sideband_reco_bins = 0;
  for ( int b = 0; b < syst.reco_bins_.size(); ++b ) {
    const auto& rbin = syst.reco_bins_.at( b );
    if ( rbin.type_ == kSidebandRecoBin ) ++num_sideband_reco_bins;
    else ++num_ordinary_reco_bins;
  }

  int num_true_signal_bins = 0;
  for ( int t = 0; t < syst.true_bins_.size(); ++t ) {
    const auto& tbin = syst.true_bins_.at( t );
    if ( tbin.type_ == kSignalTrueBin ) ++num_true_signal_bins;
  }

  std::cout << "NUM ORDINARY RECO BINS = " << num_ordinary_reco_bins << '\n';
  std::cout << "NUM TRUE SIGNAL BINS = " << num_true_signal_bins << '\n';

  auto* matrix_map_ptr = syst.get_covariances().release();
  auto& matrix_map = *matrix_map_ptr;

  auto* cov_mat = matrix_map.at( "total" ).cov_matrix_.get();

  constexpr int NUM_DAGOSTINI_ITERATIONS = 2;
  constexpr bool USE_ADD_SMEAR = true;

  std::unique_ptr< Unfolder > unfolder (
    //new DAgostiniUnfolder( NUM_DAGOSTINI_ITERATIONS )
    new DAgostiniUnfolder( DAgostiniUnfolder::ConvergenceCriterion
      ::FigureOfMerit, 0.025 )
    //new WienerSVDUnfolder( true,
    //  WienerSVDUnfolder::RegularizationMatrixType::kSecondDeriv )
  );

  UnfoldedMeasurement result = unfolder->unfold( syst );

  //// Test against RooUnfold implementation of the D'Agostini method
  //auto test_response = get_test_response( mcc9 );
  //RooUnfoldBayes roounfold_unfolder( test_response.get(),
  //  test_response->Hmeasured(), NUM_DAGOSTINI_ITERATIONS );
  //auto measured = mcc9.get_measured_events();
  //roounfold_unfolder.SetMeasuredCov( *measured.cov_matrix_ );
  //roounfold_unfolder.IncludeSystematics( 1 );

  //auto* roo_hist = roounfold_unfolder.Hreco(
  //  RooUnfold::ErrorTreatment::kCovariance );

  //auto roo_cov = roounfold_unfolder.Ereco(
  //  RooUnfold::ErrorTreatment::kCovariance );

  //for ( int t = 0; t < num_true_signal_bins; ++t ) {
  //  double mine = result.unfolded_signal_->operator()( t, 0 );
  //  double my_err = std::sqrt( std::max(0.,
  //    result.cov_matrix_->operator()( t, t )) );
  //  double theirs = roo_hist->GetBinContent( t + 1 );
  //  double their_err = std::sqrt( std::max(0., roo_cov(t, t)) );

  //  std::cout << "t = " << t << ' ' << mine << " ± " << my_err
  //    << "  " << theirs << " ± " << their_err << "  "
  //    << (mine - theirs) / theirs << " ± "
  //    << (my_err - their_err) / their_err << '\n';
  //}

  //auto smearcept = mcc9.get_cv_smearceptance_matrix();
  //auto true_signal = syst.get_cv_true_signal();

  //// Sanity check (these were equal to each other as expected)
  //// int num_reco_bins = smearcept->GetNrows();
  ////TMatrixD reco_signal( *smearcept, TMatrixD::kMult, *true_signal );
  ////auto data_sig = mcc9.get_cv_ordinary_reco_signal();
  ////for ( int r = 0; r < num_reco_bins; ++r ) {
  ////  double r1 = reco_signal( r, 0 );
  ////  double r2 = data_sig->operator()( r, 0 );
  ////  std::cout << "r = " << r << ", r1 = " << r1 << ", r2 = " << r2
  ////    << ", diff = " << r1 - r2 << '\n';
  ////}

  //auto meas = mcc9.get_measured_events();
  //// ASIMOV TEST: USE CV RECO SIGNAL PREDICTION INSTEAD OF
  //// BACKGROUND-SUBTRACTED DATA
  ////const auto& data_signal = meas.reco_signal_;
  //auto data_signal = mcc9.get_cv_ordinary_reco_signal();
  //const auto& data_covmat = meas.cov_matrix_;

  //auto result = unfolder->unfold( *data_signal, *data_covmat,
  //  *smearcept, *true_signal );

  ////// Test with original implementation
  ////TVectorD true_sig_vec( num_true_bins );
  ////for ( int t = 0; t < num_true_bins; ++t ) {
  ////  true_sig_vec( t ) = true_signal->operator()( t, 0 );
  ////}

  ////TVectorD data_sig_vec( num_ordinary_reco_bins );
  ////for ( int r = 0; r < num_ordinary_reco_bins; ++r ) {
  ////  data_sig_vec( r ) = data_signal->operator()( r, 0 );
  ////}

  ////TMatrixD A_C_svd( num_true_bins, num_true_bins );
  ////TVectorD WF_svd( num_true_bins );
  ////TMatrixD UnfoldCov_svd( num_true_bins, num_true_bins );

  ////TVectorD wsvd_unfolded_signal = WienerSVD( *smearcept, true_sig_vec,
  ////  data_sig_vec, *data_covmat, 0, 0, A_C_svd, WF_svd, UnfoldCov_svd, 0. );

  //for ( int t = 0; t < num_true_bins; ++t ) {
  //  double t_evts = true_signal->operator()( t, 0 );
  //  double evts = result.unfolded_signal_->operator()( t, 0 );
  //  double error = std::sqrt( std::max(0.,
  //    result.cov_matrix_->operator()( t, t )) );
  //  //double evts_svd = wsvd_unfolded_signal( t );
  //  //double error_svd = std::sqrt( std::max(0.,
  //  //  UnfoldCov_svd( t, t )) );
  //  std::cout << "t = " << t << ", t_evts = " << t_evts
  //    << ", evts = " << evts << " ± " << error
  //    << ", diff = " << t_evts - evts << '\n';
  //}

  // Set the event counts in each bin of the histogram that displays the
  // unfolded result. Note that we don't care about the background true bins
  // (which are assumed to follow all of the signal true bins) since we've
  // subtracted out an estimate of the background before unfolding.
  for ( int t = 0; t < num_true_bins; ++t ) {
    double evts = 0.;
    double error = 0.;
    if ( t < num_true_signal_bins ) {
      evts = result.unfolded_signal_->operator()( t, 0 );
      error = std::sqrt( std::max(0., result.cov_matrix_->operator()( t, t )) );
    }

    // We need to use one-based indices while working with TH1D bins
    unfolded_events->SetBinContent( t + 1, evts );
    unfolded_events->SetBinError( t + 1, error );
  }

  unfolded_events->SetStats( false );
  unfolded_events->SetLineColor( kBlack );
  unfolded_events->SetLineWidth( 3 );
  unfolded_events->GetXaxis()->SetRangeUser( 0, num_true_signal_bins );

  // Multiply the truth-level GENIE prediction by the additional smearing
  // matrix
  TMatrixD* A_C = result.add_smear_matrix_.get();
  multiply_1d_hist_by_matrix( A_C, genie_cv_truth );

  genie_cv_truth->SetStats( false );
  genie_cv_truth->SetLineColor( kRed );
  genie_cv_truth->SetLineWidth( 3 );
  genie_cv_truth->SetLineStyle( 9 );

  unfolded_events->Draw( "e" );
  genie_cv_truth->Draw( "hist same" );

  if ( using_fake_data ) {

    // Multiply the fake data truth by the additional smearing matrix
    multiply_1d_hist_by_matrix( A_C, fake_data_truth_hist );

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
  auto* sb_ptr = new SliceBinning( "../mybins_tki.txt" );
  //auto* sb_ptr = new SliceBinning( "../mybins_mcc9_2D_both-new.txt" );
  auto& sb = *sb_ptr;

  // Get the factors needed to convert to cross-section units
  double total_pot = syst.total_bnb_data_pot_;
  double integ_flux = integrated_numu_flux_in_FV( total_pot );
  double num_Ar = num_Ar_targets_in_FV();

  // Retrieve the true-space expected event counts from NUISANCE output files
  // for each available generator model
  double conv_factor = num_Ar * integ_flux;
  auto generator_truth_map = get_true_events_nuisance( sample_info,
    conv_factor );

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

  if ( USE_ADD_SMEAR ) {

    // Get access to the additional smearing matrix
    const TMatrixD& A_C = *result.add_smear_matrix_;

    // Start with the fake data truth if present
    if ( using_fake_data ) {
      TMatrixD ac_truth( A_C, TMatrixD::kMult, fake_data_truth );
      fake_data_truth = ac_truth;
    }

    // Also transform the GENIE CV model
    TMatrixD genie_cv_temp( A_C, TMatrixD::kMult, genie_cv_truth_vec );
    genie_cv_truth_vec = genie_cv_temp;

    // Now do the other generator predictions
    for ( const auto& pair : generator_truth_map ) {
      const auto& model_name = pair.first;
      TMatrixD* truth_mat = pair.second;

      TMatrixD ac_temp( A_C, TMatrixD::kMult, *truth_mat );
      *truth_mat = ac_temp;
    }
  }

  for ( size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx ) {

    //if ( sl_idx == 14 || sl_idx == 16 | sl_idx == 22 ) ++sl_idx;
    const auto& slice = sb.slices_.at( sl_idx );

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
    auto* slice_gen_map_ptr = new std::map< std::string, SliceHistogram* >();
    auto& slice_gen_map = *slice_gen_map_ptr;

    slice_gen_map[ "unfolded data" ] = slice_unf;
    if ( using_fake_data ) {
      slice_gen_map[ "truth" ] = slice_truth;
    }
    slice_gen_map[ "MicroBooNE Tune" ] = slice_cv;

    for ( const auto& pair : generator_truth_map ) {
      const auto& model_name = pair.first;
      TMatrixD* truth_mat = pair.second;

      SliceHistogram* temp_slice = SliceHistogram::make_slice_histogram(
        *truth_mat, slice, nullptr );

      slice_gen_map[ model_name ] = temp_slice;
    }

    int var_count = 0;
    std::string diff_xsec_denom;
    std::string diff_xsec_units_denom;
    double other_var_width = 1.;
    for ( const auto& ov_spec : slice.other_vars_ ) {
      double high = ov_spec.high_bin_edge_;
      double low = ov_spec.low_bin_edge_;
      const auto& var_spec = sb.slice_vars_.at( ov_spec.var_index_ );
      if ( high != low && std::abs(high - low) < BIG_DOUBLE ) {
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
      const auto& var_spec = sb.slice_vars_.at( av_idx );
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
        slice_y_title += "d^{" + std::to_string( var_count ) + "}#sigma";
        slice_y_title += '/' + diff_xsec_denom;
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

      const auto& file_info = truth_file_map.at( name );
      slice_h->hist_->SetLineColor( file_info.color_ );
      slice_h->hist_->SetLineStyle( file_info.style_ );
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
  return;

  // ******* Also look at reco-space results
  TH1D* reco_data_hist = dynamic_cast< TH1D* >(
    syst.data_hists_.at( NFT::kOnBNB )->Clone( "reco_data_hist" )
  );
  TH1D* reco_ext_hist = syst.data_hists_.at( NFT::kExtBNB ).get();
  const auto& cv_univ = syst.cv_universe();
  int num_reco_bins = reco_data_hist->GetNbinsX();

  // Clone the reco data hist twice. We will fill the clones with the CV
  // MC+EXT prediction and the constrained one
  TH1D* reco_mc_and_ext_hist = dynamic_cast< TH1D* >(
    reco_data_hist->Clone( "reco_mc_and_ext_hist" )
  );
  reco_mc_and_ext_hist->Reset();
  reco_mc_and_ext_hist->Add( reco_ext_hist );
  reco_mc_and_ext_hist->Add( cv_univ.hist_reco_.get() );

  TH1D* reco_constrained_hist = dynamic_cast< TH1D* >(
    reco_data_hist->Clone( "reco_constrained_hist" )
  );
  reco_constrained_hist->Reset();

  // Get the post-constraint event counts and covariance matrix in the
  // signal region
  auto meas = syst.get_measured_events();

  for ( int rb = 0; rb < num_reco_bins; ++rb ) {

    double mcc9_err = std::sqrt(
      std::max( 0., cov_mat->GetBinContent(rb + 1, rb + 1) )
    );
    reco_mc_and_ext_hist->SetBinError( rb + 1, mcc9_err );

    if ( rb >= num_ordinary_reco_bins ) {
      double data_evts = reco_data_hist->GetBinContent( rb + 1 );
      reco_constrained_hist->SetBinContent( rb + 1, data_evts );
      reco_constrained_hist->SetBinError( rb + 1, 0. );
    }
    else {
      double constr_pred = meas.reco_mc_plus_ext_->operator()( rb, 0 );
      double constr_err = std::sqrt(
        std::max( 0., meas.cov_matrix_->operator()(rb, rb) )
      );

      reco_constrained_hist->SetBinContent( rb + 1, constr_pred );
      reco_constrained_hist->SetBinError( rb + 1, constr_err );
    }

  }

  TCanvas* c2 = new TCanvas;

  reco_data_hist->SetLineColor( kBlack );
  reco_data_hist->SetLineWidth( 5 );

  reco_mc_and_ext_hist->SetLineColor( kRed );
  reco_mc_and_ext_hist->SetLineStyle( 2 );
  reco_mc_and_ext_hist->SetLineWidth( 4 );

  reco_constrained_hist->SetLineColor( kBlue );
  reco_constrained_hist->SetLineStyle( 9 );
  reco_constrained_hist->SetLineWidth( 4 );

  reco_data_hist->Draw( "e" );
  reco_mc_and_ext_hist->Draw( "same hist e" );
  reco_constrained_hist->Draw( "same hist e" );

  reco_data_hist->Draw( "same e" );

  TLegend* lg2 = new TLegend( 0.15, 0.7, 0.3, 0.85 );
  lg2->AddEntry( reco_data_hist, using_fake_data ? "fake data" : "data",
    "l" );
  lg2->AddEntry( reco_mc_and_ext_hist, "uB tune + EXT", "l" );
  lg2->AddEntry( reco_constrained_hist, "post-constraint", "l" );

  lg2->Draw( "same" );

}

int main() {
  test_unfolding();
  return 0;
}
