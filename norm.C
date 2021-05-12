// ROOT includes
#include "TFile.h"

// STV analysis includes
#include "FilePropertiesManager.hh"
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

void norm() {

  std::string input_respmat_file_name( "respmat-myconfig_one_bin.root" );

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
    " (cm^2 / Ar / x-axis unit)", num_reco_bins, 0., num_reco_bins );
  result_hist->Sumw2();

  TH2D* total_cov_matrix = matrix_map.at( "total" ).cov_matrix_.get();

  const auto& cv_univ = syst.cv_universe();
  for ( int rb = 1; rb <= num_reco_bins; ++rb ) {
    double xsec = syst.forward_folded_xsec( cv_univ, rb );
    double err2 = total_cov_matrix->GetBinContent( rb, rb );

    double err = std::sqrt( std::max(0., err2) );

    result_hist->SetBinContent( rb, xsec );
    result_hist->SetBinError( rb, err );
  }

  TCanvas* c1 = new TCanvas;
  result_hist->SetLineColor( kBlack );
  result_hist->SetLineWidth( 3 );
  result_hist->SetStats( false );
  result_hist->Draw( "e" );

  TCanvas* c2 = new TCanvas;
  total_cov_matrix->Draw( "colz" );

}
