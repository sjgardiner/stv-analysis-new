#pragma once

// STV analysis includes
#include "UniverseMaker.hh"

struct CovMatResults {

  CovMatResults() {}

  CovMatResults( TH2D* signal_cm, TH2D* bkgd_cm, TH1D* rc_signal_cv,
    TH1D* rc_bkgd_cv, bool frac ) : signal_cov_mat_( signal_cm ),
    bkgd_cov_mat_( bkgd_cm ), reco_signal_cv_( rc_signal_cv ),
    reco_bkgd_cv_( rc_bkgd_cv ), fractional_( frac ) {}

  // Returns the number of reco bins (i.e., the number of matrix elements along
  // either axis)
  int num_reco_bins() const { return signal_cov_mat_->GetNbinsX(); }

  // Returns the fractional covariance matrix element for the signal only.
  // NOTE: bin numbering in this interface is one-based (not zero-based)
  // to match the conventions of ROOT histograms.
  double frac_covariance_signal( int bin_a, int bin_b ) const {
    double covar = signal_cov_mat_->GetBinContent( bin_a, bin_b );
    if ( fractional_ ) return covar;
    else {
      double cv_a = reco_signal_cv_->GetBinContent( bin_a );
      double cv_b = reco_signal_cv_->GetBinContent( bin_b );
      double cv_prod = cv_a * cv_b;
      if ( cv_prod == 0. ) return 0.;
      covar /= cv_prod;
    }
    return covar;
  }

  // TODO: reduce code duplication here
  // Returns the fractional covariance matrix element for the background only.
  double frac_covariance_bkgd( int bin_a, int bin_b ) const {
    double covar = bkgd_cov_mat_->GetBinContent( bin_a, bin_b );
    if ( fractional_ ) return covar;
    else {
      double cv_a = reco_bkgd_cv_->GetBinContent( bin_a );
      double cv_b = reco_bkgd_cv_->GetBinContent( bin_b );
      double cv_prod = cv_a * cv_b;
      if ( cv_prod == 0. ) return 0.;
      covar /= cv_prod;
    }
    return covar;
  }

  // TODO: reduce code duplication here
  // Returns the fractional covariance matrix element the full MC prediction
  // (signal plus background)
  double frac_covariance_total( int bin_a, int bin_b ) const {
    double covar_bkgd = bkgd_cov_mat_->GetBinContent( bin_a, bin_b );
    double covar_signal = signal_cov_mat_->GetBinContent( bin_a, bin_b );

    double cv_a_bkgd = reco_bkgd_cv_->GetBinContent( bin_a );
    double cv_b_bkgd = reco_bkgd_cv_->GetBinContent( bin_b );

    double cv_a_signal = reco_signal_cv_->GetBinContent( bin_a );
    double cv_b_signal = reco_signal_cv_->GetBinContent( bin_b );

    // If these are fractional covariances already, we need to scale them
    // back into regular covariances before summing signal and background.
    if ( fractional_ ) {
      covar_bkgd *= cv_a_bkgd * cv_b_bkgd;
      covar_signal *= cv_a_signal * cv_b_signal;
    }

    // OK, now add them and divide by the total CV prediction to get a total
    // fractional covariance
    double covar = covar_bkgd + covar_signal;

    double cv_a = cv_a_bkgd + cv_a_signal;
    double cv_b = cv_b_bkgd + cv_b_signal;

    double cv_prod = cv_a * cv_b;

    if ( cv_prod == 0. ) return 0.;
    covar /= cv_prod;

    return covar;
  }

  std::unique_ptr< TH2D > signal_cov_mat_;
  std::unique_ptr< TH2D > bkgd_cov_mat_;

  std::unique_ptr< TH1D > reco_signal_cv_;
  std::unique_ptr< TH1D > reco_bkgd_cv_;

  bool fractional_ = false;

};

using MatrixMap = std::map< std::string,
  std::map<std::string, CovMatResults> >;

// Helper function that saves the contents of the map of covariance
// matrices to an output ROOT file
void save_matrix_map( const MatrixMap& matrix_map, TFile& out_tfile )
{
  TDirectoryFile* root_tdir = new TDirectoryFile( "covMat",
    "covariance matrices", "", &out_tfile );

  root_tdir->cd();

  // Save a set of the covariance matrix labels to assist in easily
  // reading them back in
  // TODO: revisit this and consider a simpler technique
  std::set< std::string > cov_mat_labels;

  for ( const auto& pair : matrix_map ) {

    const std::string& ntuple_file = pair.first;
    const auto& results_map = pair.second;

    std::string subdir = ntuple_subfolder_from_file_name( ntuple_file );

    TDirectoryFile* ntuple_tdir = new TDirectoryFile( subdir.c_str(),
      "covariance matrices", "", root_tdir );

    ntuple_tdir->cd();

    for ( const auto& results_pair : results_map ) {
      const std::string& label = results_pair.first;
      const CovMatResults& results = results_pair.second;

      cov_mat_labels.insert( label );

      results.signal_cov_mat_->Write( (label + "_signal_cov_mat").c_str() );
      results.bkgd_cov_mat_->Write( (label + "_bkgd_cov_mat").c_str() );

      results.reco_signal_cv_->Write( (label + "_reco_signal_cv").c_str() );
      results.reco_bkgd_cv_->Write( (label + "_reco_bkgd_cv").c_str() );

      TParameter<bool> frac( (label + "_fractional").c_str(),
        results.fractional_ );
      frac.Write();

    } // covariance matrix categories

  } // ntuple files

  out_tfile.WriteObject( &cov_mat_labels, "cov_mat_labels" );

}

// Helper function that reinstantiates a map of covariance
// matrices from an input ROOT file created with save_matrix_map()
MatrixMap load_matrix_map( TFile& in_tfile ) {

  MatrixMap retrieved_map;

  // Here we cheat a bit. Since the pot_map is also saved to the file and has
  // the same ntuple file names as keys, we can iterate over it instead of
  // subdirectories of the root TDirectoryFile. I'm sure there's a better way
  // to do this with the TDirectoryFile itself, but for now I use this quick
  // hack.
  // TODO: revisit this
  std::map< std::string, float >* pot_map = nullptr;

  in_tfile.GetObject( "pot_map", pot_map );

  if ( !pot_map ) {
    throw std::runtime_error( "Missing POT map!" );
  }

  std::set< std::string >* cov_mat_labels;
  in_tfile.GetObject( "cov_mat_labels", cov_mat_labels );

  if ( !cov_mat_labels ) {
    throw std::runtime_error( "Missing covMat labels!" );
  }

  TDirectoryFile* root_tdir = nullptr;

  in_tfile.GetObject( "covMat", root_tdir );
  if ( !root_tdir ) {
    throw std::runtime_error( "Could not find covMat TDirectoryFile" );
  }

  root_tdir->cd();

  // Build a vector of file names from the POT map. Then add "total_mc"
  // which is also included for the POT-summed total central-value MC
  // results.
  std::vector< std::string > ntuple_files;
  for ( const auto& pair : *pot_map ) {
    const std::string& ntuple_file = pair.first;
    //float pot = pair.second;
    ntuple_files.push_back( ntuple_file );
  }

  ntuple_files.push_back( "total_mc" );

  // Now loop over the file names and reconstruct the full matrix_map
  for ( const auto& ntuple_file : ntuple_files ) {

    if ( !retrieved_map.count(ntuple_file) ) {
      retrieved_map[ ntuple_file ] = std::map< std::string, CovMatResults >();
    }

    std::string subdir = ntuple_subfolder_from_file_name( ntuple_file );

    TDirectoryFile* ntuple_tdir = nullptr;

    root_tdir->GetObject( subdir.c_str(), ntuple_tdir );

    if ( !ntuple_tdir ) {
      throw std::runtime_error( "Missing covMat subdirectory " + subdir );
    }

    ntuple_tdir->cd();

    auto& ntuple_submap = retrieved_map.at( ntuple_file );

    for ( const auto& label : *cov_mat_labels ) {

      // Skip the special summed covariance matrix labels that only
      // exist for the "total_mc" TDirectoryFile
      if ( ntuple_file != "total_mc"
        && (label == "xsec_all" || label == "xsec_unisim") ) continue;

      // TODO: add error handling here for missing objects
      TH2D* signal_cov = nullptr;
      TH2D* bkgd_cov = nullptr;

      TH1D* signal_cv = nullptr;
      TH1D* bkgd_cv = nullptr;

      ntuple_tdir->GetObject( (label + "_signal_cov_mat").c_str(), signal_cov );
      signal_cov->SetDirectory( nullptr );

      ntuple_tdir->GetObject( (label + "_bkgd_cov_mat").c_str(), bkgd_cov );
      bkgd_cov->SetDirectory( nullptr );

      ntuple_tdir->GetObject( (label + "_reco_signal_cv").c_str(), signal_cv );
      signal_cv->SetDirectory( nullptr );

      ntuple_tdir->GetObject( (label + "_reco_bkgd_cv").c_str(), bkgd_cv );
      bkgd_cv->SetDirectory( nullptr );

      TParameter<bool>* frac = nullptr;
      ntuple_tdir->GetObject( (label + "_fractional").c_str(), frac );

      CovMatResults temp_results( signal_cov, bkgd_cov, signal_cv, bkgd_cv,
        frac->GetVal() );

      ntuple_submap[ label ] = std::move( temp_results );

    } // covariance matrix categories

  } // ntuple files

  return retrieved_map;
}
