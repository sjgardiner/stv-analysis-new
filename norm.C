// ROOT includes
#include "TFile.h"

// STV analysis includes
#include "FilePropertiesManager.hh"
#include "PlotUtils.hh"
#include "SystematicsCalculator.hh"

using NFT = NtupleFileType;

void covMat( const std::string& input_respmat_file_name ) {

  auto* syst_ptr = new SystematicsCalculator( input_respmat_file_name );
  auto& syst = *syst_ptr;

  // Keys are covariance matrix types, values are CovMatrix objects that
  // represent the corresponding matrices
  auto* matrix_map_ptr = syst.get_covariances().release();
  auto& matrix_map = *matrix_map_ptr;

  // Build the final event counts and total covariance matrix for the MC + EXT
  // prediction ("pred").
  int num_reco_bins = syst.reco_bins_.size();
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

  int num_true_bins = syst.true_bins_.size();
  TH1D* reco_bkgd_hist = syst.cv_universe().hist_2d_->ProjectionY(
    "reco_bkgd", num_reco_bins + 1, num_true_bins );

  reco_ext_plus_bkgd_hist->Add( reco_bkgd_hist );

  // Keys are labels, values are fractional uncertainty histograms
  auto* fr_unc_hists = new std::map< std::string, TH1D* >();
  auto& frac_uncertainty_hists = *fr_unc_hists;

  // Terms needed for the total covariance matrix
  const std::vector< std::string > total_cov_mat_keys = { "detVar_total",
    "flux", "reint", "xsec_total", "POT", "numTargets", "MCstats", "EXTstats"
  };

  int color = 1;
  for ( const auto& key : total_cov_mat_keys ) {

    const auto& temp_results = matrix_map.at( key );

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

// The covariance matrix passed to this function should represent the total
// uncertainty on the difference between the data and prediction
double get_chi2( const TH1D& data_hist, const TH1D& pred_hist,
  const CovMatrix& cov_mat )
{
  int num_reco_bins = data_hist.GetNbinsX();
  if ( pred_hist.GetNbinsX() != num_reco_bins ) {
    throw std::runtime_error( "Incompatible vector sizes in chi^2"
      " calculation" );
  }

  // Evaluate the inverse of the covariance matrix
  auto inverse_cov_mat = cov_mat.get_matrix();

  // Pre-scale before inversion to avoid numerical problems
  constexpr double BIG_SCALING_FACTOR = 1e76;
  inverse_cov_mat->operator*=( BIG_SCALING_FACTOR );

  // Do the inversion
  inverse_cov_mat->Invert();

  // Undo the scaling by re-applying it to the inverse matrix
  inverse_cov_mat->operator*=( BIG_SCALING_FACTOR );

  // Double-check that we get a unit matrix by multiplying the
  // original by its inverse
  //TMatrixD unit_mat( *sym_mat, TMatrixD::kMult, *inverse_cov_mat );
  //unit_mat.Print();

  // Create a 1D vector containing the difference between the data
  // and the prediction in each reco bin
  TMatrixD diff_vec( 1, num_reco_bins );
  for ( int a = 0; a < num_reco_bins; ++a ) {
    double data_counts = data_hist.GetBinContent( a + 1 );
    double pred_counts = pred_hist.GetBinContent( a + 1 );
    diff_vec( 0, a ) = data_counts - pred_counts;
  }

  // Multiply diff * covMat^{-1} * diff^{T} to get chi-squared
  TMatrixD temp1( diff_vec, TMatrixD::kMult, *inverse_cov_mat );
  TMatrixD temp2( temp1, TMatrixD::kMult, diff_vec.T() );

  // We'll now have a 1x1 matrix containing the chi-squared value
  double chi2 = temp2( 0, 0 );
  return chi2;
}

void compare_mcc8_mcc9( const std::string& input_respmat_file_name,
  const std::string& mcc8_hist_suffix,
  const std::string& mcc8_hist_title )
{
  TFile* mcc8_file = new TFile( "CCNp_data_MC_cov_dataRelease.root", "read" );

  auto* syst_ptr = new SystematicsCalculator( input_respmat_file_name );
  auto& syst = *syst_ptr;

  // Keys are covariance matrix types, values are CovMatrix objects that
  // represent the corresponding matrices
  auto* matrix_map_ptr = syst.get_covariances().release();
  auto& matrix_map = *matrix_map_ptr;

  // Build the final measurement histogram from the CV forward-folded
  // differential cross sections in each reco bin
  int num_reco_bins = syst.reco_bins_.size();
  TH1D* result_hist = new TH1D( "result_hist", "; reco bin; differential xsec"
    " (cm^{2} / Ar / x-axis unit)", num_reco_bins, 0., num_reco_bins );
  result_hist->Sumw2();
  result_hist->SetDirectory( nullptr );

  auto& total_cov_matrix = matrix_map.at( "total" );
  TH2D* total_cov_matrix_hist = total_cov_matrix.cov_matrix_.get();

  const auto& cv_univ = syst.cv_universe();
  for ( int rb = 1; rb <= num_reco_bins; ++rb ) {
    double xsec = syst.forward_folded_xsec( cv_univ, rb );
    double err2 = total_cov_matrix_hist->GetBinContent( rb, rb );

    double err = std::sqrt( std::max(0., err2) );

    result_hist->SetBinContent( rb, xsec );
    result_hist->SetBinError( rb, err );
  }

  // Get the MCC8 data and covariance matrix
  TH1D* mcc8_data = nullptr;
  mcc8_file->GetObject( ("DataXsec_" + mcc8_hist_suffix).c_str(), mcc8_data );

  mcc8_data->SetTitle( mcc8_hist_title.c_str()  );

  TH2D* mcc8_cov_hist = nullptr;
  mcc8_file->GetObject( ("CovarianceMatrix_" + mcc8_hist_suffix).c_str(),
    mcc8_cov_hist );

  // Convert the MCC8 units from 10^{-38} cm^2 / nucleon to cm^2 / Ar
  constexpr double mcc8_scale_factor = 1e-38 * 40;
  mcc8_data->Scale( mcc8_scale_factor );
  mcc8_cov_hist->Scale( std::pow(mcc8_scale_factor, 2) );

  // Correct the MCC9 results for the missing bin width factors in the
  // denominator. Since we used the MCC8 binning, we can get this information
  // from the published results.
  for ( int a = 1; a <= num_reco_bins; ++a ) {

    double width_a = mcc8_data->GetBinWidth( a );
    double old_counts = result_hist->GetBinContent( a );
    double old_err = result_hist->GetBinError( a );

    result_hist->SetBinContent( a, old_counts / width_a );
    result_hist->SetBinError( a, old_err / width_a );

    for ( int b = 1; b <= num_reco_bins; ++b ) {
      double width_b = mcc8_data->GetBinWidth( b );

      double old_cov = total_cov_matrix_hist->GetBinContent( a, b );
      double new_cov = old_cov / ( width_a * width_b );
      total_cov_matrix_hist->SetBinContent( a, b, new_cov );
    }
  }

  // Create an MCC9 histogram with "physics binning" to match the MCC8 data
  TH1D* mcc9_data = dynamic_cast< TH1D* >( mcc8_data->Clone("mcc9_data") );
  mcc9_data->SetDirectory( nullptr );
  for ( int a = 1; a <= num_reco_bins; ++a ) {
    double xsec = result_hist->GetBinContent( a );
    double err = result_hist->GetBinError( a );
    mcc9_data->SetBinContent( a, xsec );
    mcc9_data->SetBinError( a, err );
  }

  // Create a CovMatrix object from the MCC8 covariance matrix histogram.
  // Add it to ours to obtain the total covariance on the difference between
  // the two measurements
  mcc8_cov_hist->SetDirectory( nullptr );
  CovMatrix mcc8_cov_mat( mcc8_cov_hist );

  // If we're comparing the muon momentum distributions, there's a small
  // discrepancy: I included an overflow bin while MCC8 did not. This
  // makes the covariance matrices have different dimensions, which causes
  // problems when computing chi-squared. I'll manually fix this by dropping
  // the overflow bin here.
  if ( mcc8_hist_suffix == "mumom" ) {
    // Create a new MCC9 covariance matrix using only the reco bins considered
    // in MCC8. For convenience, clone the MCC8 histogram to start (we get
    // the binning right "for free" and just have to overwrite the contents).
    TH2D* mcc9_truncated_cov_hist = dynamic_cast< TH2D* >(
      mcc8_cov_hist->Clone("mcc9_truncated_cov_hist") );
    mcc9_truncated_cov_hist->SetDirectory( nullptr );

    int num_mcc8_reco_bins = mcc8_data->GetNbinsX();
    for ( int a = 1; a <= num_mcc8_reco_bins; ++a ) {
      for ( int b = 1; b <= num_mcc8_reco_bins; ++b ) {
        double covariance_mcc9 = total_cov_matrix_hist->GetBinContent( a, b );
        mcc9_truncated_cov_hist->SetBinContent( a, b, covariance_mcc9 );
      }
    }

    // We've made the truncated MCC9 covariance matrix. Replace the old one
    // with it and move on.
    total_cov_matrix.cov_matrix_.reset( mcc9_truncated_cov_hist );

  } // muon momentum comparison

  CovMatrix comp_cov_mat;
  comp_cov_mat += total_cov_matrix;
  comp_cov_mat += mcc8_cov_mat;

  // Compute a chi-squared value for the comparison of the measurements
  double chi2 = get_chi2( *mcc9_data, *mcc8_data, comp_cov_mat );
  std::cout << mcc8_hist_suffix << ": \u03C7\u00b2 = " << chi2
    << " / " << mcc8_data->GetNbinsX() << " bins\n";

  // Draw the comparison plot
  TCanvas* c1 = new TCanvas;
  c1->SetLeftMargin( 0.12 );
  c1->SetRightMargin( 0.05 );

  mcc8_data->GetXaxis()->SetLabelSize( 0.045 );
  mcc8_data->GetYaxis()->SetLabelSize( 0.045 );

  mcc8_data->SetLineColor( kAzure - 1 );
  mcc8_data->SetMarkerColor( kAzure - 1 );
  mcc8_data->SetLineWidth( 3 );
  mcc8_data->SetStats( false );

  mcc9_data->SetLineColor( kBlack );
  mcc9_data->SetLineWidth( 3 );
  mcc9_data->SetStats( false );

  mcc8_data->Draw( "e" );
  mcc9_data->Draw( "e same" );

  // Get a reasonable plot range by plotting the
  // cross section with the greater maximum bin value first
  double mcc8_max = mcc8_data->GetMaximum();
  double mcc9_max = mcc9_data->GetMaximum();
  if ( mcc8_max > mcc9_max ) {
    mcc8_data->Draw( "e" );
    mcc9_data->Draw( "e same" );
  }
  else {
    mcc9_data->Draw( "e" );
    mcc8_data->Draw( "e same" );
    // We plotted MCC9 first, so redraw it on top of the MCC8
    // result to ensure that it will be visible
    mcc9_data->Draw( "e same" );
  }

  TLegend* lg = new TLegend( 0.15, 0.65, 0.45, 0.88 );

  lg->AddEntry( mcc9_data, "MCC9", "le" );
  lg->AddEntry( mcc8_data, "PRD 102, 112013", "le" );

  std::ostringstream lg_oss;
  lg_oss << std::setprecision(3) << "#chi^{2} = " << chi2
    << " / " << mcc8_data->GetNbinsX() << " bins";
  TObject* dummy_ptr = nullptr;
  lg->AddEntry( dummy_ptr, lg_oss.str().c_str(), "" );
  lg->Draw( "same" );

  std::string uboone_label = get_legend_title( syst.total_bnb_data_pot_ );
  TLatex* ltx = new TLatex( 0.38, 0.92, uboone_label.c_str() );
  ltx->SetTextSize( 0.045 );
  ltx->SetNDC( true ); // Use the pad coordinate system, not the axes
  ltx->Draw( "same" );

}

void norm() {

  compare_mcc8_mcc9( "/uboone/data/users/gardiner/respmat_mcc8-cth_mu.root",
    "muangle", "; cos#theta_{#mu}; d#sigma/dcos#theta_{#mu} (cm^{2} / Ar)" );

  compare_mcc8_mcc9( "/uboone/data/users/gardiner/respmat_mcc8-cth_p.root",
    "pangle", "; cos#theta_{p}; d#sigma/dcos#theta_{p} (cm^{2} / Ar)" );

  compare_mcc8_mcc9( "/uboone/data/users/gardiner/respmat_mcc8-p_mu.root",
    "mumom", "; p_{#mu} (GeV); d#sigma/dp_{#mu} (cm^{2} / GeV / Ar)" );

  compare_mcc8_mcc9( "/uboone/data/users/gardiner/respmat_mcc8-p_p.root",
    "pmom", "; p_{p} (GeV); d#sigma/dp_{p} (cm^{2} / GeV / Ar)" );

  compare_mcc8_mcc9( "/uboone/data/users/gardiner/respmat_mcc8-th_mu_p.root",
    "thetamup", "; #theta_{#mu-p} (radian); d#sigma/d#theta_{#mu-p}"
    " (cm^{2} / radian / Ar)" );

}
