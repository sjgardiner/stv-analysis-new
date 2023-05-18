#pragma once

// Standard library includes
#include <memory>

// ROOT includes
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMatrixD.h"
#include "TParameter.h"

// STV analysis includes
#include "FilePropertiesManager.hh"
#include "UniverseMaker.hh"

// Helper function template that retrieves an object from a TDirectoryFile
// and loads a pointer to it into a std::unique_ptr of the correct type
template< typename T > std::unique_ptr< T > get_object_unique_ptr(
  const std::string& namecycle, TDirectory& df )
{
  T* temp_ptr = nullptr;
  df.GetObject( namecycle.c_str(), temp_ptr );

  // Set the directory to nullptr in the case of ROOT histograms. This will
  // avoid extra deletion attempts
  TH1* temp_hist_ptr = dynamic_cast< TH1* >( temp_ptr );
  if ( temp_hist_ptr ) temp_hist_ptr->SetDirectory( nullptr );

  return std::unique_ptr< T >( temp_ptr );
}

// Helper function to use on a new Universe object's histograms. Prevents
// auto-deletion problems by disassociating the Universe object's histograms
// from any TDirectory. Also turns off the stats box when plotting for
// convenience.
void set_stats_and_dir( Universe& univ ) {
  univ.hist_reco_->SetStats( false );
  univ.hist_reco_->SetDirectory( nullptr );

  univ.hist_true_->SetStats( false );
  univ.hist_true_->SetDirectory( nullptr );

  univ.hist_2d_->SetStats( false );
  univ.hist_2d_->SetDirectory( nullptr );

  univ.hist_categ_->SetStats( false );
  univ.hist_categ_->SetDirectory( nullptr );

  univ.hist_reco2d_->SetStats( false );
  univ.hist_reco2d_->SetDirectory( nullptr );
}

// Tests whether a string ends with another string. Taken from
// https://stackoverflow.com/a/874160/4081973. In C++20, you could use
// std::string::ends_with() instead.
bool has_ending( const std::string& fullString, const std::string& ending ) {
  if ( fullString.length() >= ending.length() ) {
    return ( 0 == fullString.compare(
      fullString.length() - ending.length(), ending.length(), ending)
    );
  }

  return false;
}

// Simple container for a TH2D that represents a covariance matrix
struct CovMatrix {

  CovMatrix() {}

  CovMatrix( TH2D* cov_mat ) : cov_matrix_( cov_mat ) {}

  CovMatrix( const TMatrixD& matrix ) {
    int num_bins = matrix.GetNrows();
    if ( matrix.GetNcols() != num_bins ) throw std::runtime_error( "Non-square"
      " TMatrixD passed to the constructor of CovMatrix" );

    TH2D* temp_hist = new TH2D( "temp_hist", "covariance; bin;"
      " bin; covariance", num_bins, 0., num_bins, num_bins, 0., num_bins );
    temp_hist->SetDirectory( nullptr );
    temp_hist->SetStats( false );

    for ( int r = 0; r < num_bins; ++r ) {
      for ( int c = 0; c < num_bins; ++c ) {
        double cov = matrix( r, c );
        // Note that TMatrixD element indices are zero-based while TH2D bin
        // indices are one-based. We therefore correct for this here.
        temp_hist->SetBinContent( r + 1, c + 1, cov );
      }
    }

    cov_matrix_.reset( temp_hist );
  }

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

  std::unique_ptr< TMatrixD > get_matrix() const {
    // Note that ROOT histogram bin indices are one-based to allow for
    // underflow. The TMatrixDSym element indices, on the other hand,
    // are zero-based.
    int num_cm_bins = cov_matrix_->GetNbinsX();
    auto result = std::make_unique< TMatrixD >( num_cm_bins, num_cm_bins );
    // TODO: consider doing something more efficient than setting each
    // element manually
    for ( int a = 0; a < num_cm_bins; ++a ) {
      for ( int b = 0; b < num_cm_bins; ++b ) {
        result->operator()( a, b ) = cov_matrix_->GetBinContent( a + 1, b + 1 );
      }
    }
    return result;
  }

};

// Container that holds the results of subtracting the EXT+MC background from
// the measured data event counts in ordinary reco bins
struct MeasuredEvents {

  MeasuredEvents() {}

  MeasuredEvents( TMatrixD* bkgd_subtracted_data, TMatrixD* bkgd,
    TMatrixD* mc_plus_ext, TMatrixD* cov_mat )
    : reco_signal_( bkgd_subtracted_data ), reco_bkgd_( bkgd ),
    reco_mc_plus_ext_( mc_plus_ext ), cov_matrix_( cov_mat )
  {
    if ( bkgd_subtracted_data->GetNcols() != 1 ) throw std::runtime_error(
      "Non-row-vector background-subtracted signal passed to MeasuredEvents" );

    if ( bkgd->GetNcols() != 1 ) throw std::runtime_error( "Non-row-vector"
      " background prediction passed to MeasuredEvents" );

    if ( mc_plus_ext->GetNcols() != 1 ) throw std::runtime_error(
      "Non-row-vector MC+EXT prediction passed to MeasuredEvents" );

    int num_ordinary_reco_bins = bkgd_subtracted_data->GetNrows();
    if ( cov_mat->GetNcols() != num_ordinary_reco_bins
      || cov_mat->GetNrows() != num_ordinary_reco_bins )
    {
      throw std::runtime_error( "Bad covariance matrix dimensions passed to"
        " MeasuredEvents" );
    }
  }

  // Background-subtracted data event counts in the ordinary reco bins
  std::unique_ptr< TMatrixD > reco_signal_;

  // Background that was subtracted from each reco bin to form the signal
  // measurement
  std::unique_ptr< TMatrixD > reco_bkgd_;

  // Total MC+EXT prediction in each reco bin (used to estimate the covariance
  // matrix on the data)
  std::unique_ptr< TMatrixD > reco_mc_plus_ext_;

  // Covariance matrix for the background-subtracted data
  std::unique_ptr< TMatrixD > cov_matrix_;

};

using CovMatrixMap = std::map< std::string, CovMatrix >;
using NFT = NtupleFileType;

class SystematicsCalculator {

  public:

    SystematicsCalculator( const std::string& input_respmat_file_name,
      const std::string& syst_cfg_file_name = "",
      const std::string& respmat_tdirectoryfile_name = "" );

    void load_universes( TDirectoryFile& total_subdir );

    void build_universes( TDirectoryFile& root_tdir );

    void save_universes( TDirectoryFile& out_tdf );

    const Universe& cv_universe() const {
      return *rw_universes_.at( CV_UNIV_NAME ).front();
    }

    const std::unique_ptr< Universe >& fake_data_universe() const {
      return fake_data_universe_;
    }

    std::unique_ptr< CovMatrixMap > get_covariances() const;

    // Returns a background-subtracted measurement in all ordinary reco bins
    // with the total covariance matrix and the background event counts that
    // were subtracted.
    // NOTE: this function assumes that the ordinary reco bins are all listed
    // before any sideband reco bins
    virtual MeasuredEvents get_measured_events() const;

    // Utility functions to help with unfolding
    inline std::unique_ptr< TMatrixD > get_cv_smearceptance_matrix() const {
      const auto& cv_univ = this->cv_universe();
      return this->get_smearceptance_matrix( cv_univ );
    }

    std::unique_ptr< TMatrixD > get_smearceptance_matrix(
      const Universe& univ ) const;

    std::unique_ptr< TMatrixD > get_cv_true_signal() const;

    // Returns the expected background in each ordinary reco bin (including
    // both EXT and the central-value MC prediction for beam-correlated
    // backgrounds)
    // NOTE: This function assumes that all "ordinary" reco bins are listed
    // before the sideband ones.
    inline std::unique_ptr< TMatrixD > get_cv_ordinary_reco_bkgd() const
      { return this->get_cv_ordinary_reco_helper( true ); }

    // Returns the expected signal event counts in each ordinary reco bin
    // NOTE: This function assumes that all "ordinary" reco bins are listed
    // before the sideband ones.
    inline std::unique_ptr< TMatrixD > get_cv_ordinary_reco_signal() const
      { return this->get_cv_ordinary_reco_helper( false ); }

    inline size_t get_num_signal_true_bins() const
      { return num_signal_true_bins_; }

  //protected:

    // Implements both get_cv_ordinary_reco_bkgd() and
    // get_cv_ordinary_reco_signal() in order to reduce code duplication. If
    // return_bkgd is false (true), then the background (signal) event counts
    // in each ordinary reco bin will be returned as a column vector.
    std::unique_ptr< TMatrixD > get_cv_ordinary_reco_helper(
      bool return_bkgd ) const;

    // Returns true if a given Universe represents a detector variation or
    // false otherwise
    bool is_detvar_universe( const Universe& univ ) const;

    // Overload for special cases in which the N*N covariance matrix does not
    // have dimension parameter N equal to the number of reco bins
    inline virtual size_t get_covariance_matrix_size() const
      { return reco_bins_.size(); }

    CovMatrix make_covariance_matrix( const std::string& hist_name ) const;

    // Evaluate the observable described by the covariance matrices in
    // a given universe and reco-space bin. NOTE: the reco bin index given
    // as an argument to this function is zero-based.
    virtual double evaluate_observable( const Universe& univ, int reco_bin,
      int flux_universe_index = -1 ) const = 0;

    // Evaluate a covariance matrix element for the data statistical
    // uncertainty on the observable of interest for a given pair of reco bins.
    // In cases where every event falls into a unique reco bin, only the
    // diagonal covariance matrix elements are non-vanishing. Do the
    // calculation either for BNB data (use_ext = false) or for EXT data
    // (use_ext = true). NOTE: the reco bin indices consumed by this function
    // are zero-based.
    virtual double evaluate_data_stat_covariance( int reco_bin_a,
      int reco_bin_b, bool use_ext ) const = 0;

    // Evaluate a covariance matrix element for the MC statistical uncertainty
    // on the observable of interest (including contributions from both signal
    // and background events) for a given pair of reco bins within a particular
    // universe. Typically the CV universe should be used, but MC statistical
    // uncertainties are tracked for all Universe objects in case they are
    // needed. In cases where every event falls into a unique reco bin, only
    // the diagonal covariance matrix elements are non-vanishing. NOTE: the
    // reco bin indices consumed by this function are zero-based.
    virtual double evaluate_mc_stat_covariance( const Universe& univ,
      int reco_bin_a, int reco_bin_b ) const = 0;

    // Central value universe name
    const std::string CV_UNIV_NAME = "weight_TunedCentralValue_UBGenie";

    // Beginning of the subdirectory name for the TDirectoryFile containing the
    // POT-summed histograms for the various universes across all analysis
    // ntuples. The full name is formed from this prefix and the name of the
    // FilePropertiesManager configuration file that is currently active.
    const std::string TOTAL_SUBFOLDER_NAME_PREFIX = "total_";

    // Holds reco-space histograms for data (BNB and EXT) bin counts
    std::map< NFT, std::unique_ptr<TH1D> > data_hists_;
    std::map< NFT, std::unique_ptr<TH2D> > data_hists2d_;

    // Holds universe objects for reweightable systematics
    std::map< std::string, std::vector< std::unique_ptr<Universe> > >
      rw_universes_;

    // Detector systematic universes (and the detVar CV) will be indexed using
    // ntuple file type values. We're currently scaling one set of ntuples
    // (from Run 3b) to the full dataset.
    // TODO: revisit this procedure if new detVar samples become available
    std::map< NFT, std::unique_ptr<Universe> > detvar_universes_;

    // If we are working with fake data, then this will point to a Universe
    // object containing full reco and truth information for the MC portion
    // (as opposed to the EXT contribution which is added to the reco MC
    // counts in the BNB "data" histogram). If we are working with real data,
    // then this will be a null pointer.
    std::unique_ptr< Universe > fake_data_universe_ = nullptr;

    // "Alternate CV" universes for assessing unisim systematics related to
    // interaction modeling
    std::map< NFT, std::unique_ptr<Universe> > alt_cv_universes_;

    // True bin configuration that was used to compute the universes
    std::vector< TrueBin > true_bins_;

    // Reco bin configuration that was used to compute the universes
    std::vector< RecoBin > reco_bins_;

    // Total POT exposure for the analyzed BNB data
    double total_bnb_data_pot_ = 0.;

    // Name of the systematics configuration file that should be used
    // when computing covariance matrices
    std::string syst_config_file_name_;

    // Number of entries in the reco_bins_ vector that are "ordinary" bins
    size_t num_ordinary_reco_bins_ = 0u;

    // Number of entries in the true_bins_ vector that are "signal" bins
    size_t num_signal_true_bins_ = 0u;
};

SystematicsCalculator::SystematicsCalculator(
  const std::string& input_respmat_file_name,
  const std::string& syst_cfg_file_name,
  const std::string& respmat_tdirectoryfile_name )
  : syst_config_file_name_( syst_cfg_file_name )
{
  // Get access to the FilePropertiesManager singleton class
  const auto& fpm = FilePropertiesManager::Instance();

  // If the user didn't specify a particular systematics configuration file
  // to use when calculating covariance matrices, then use the default one
  if ( syst_config_file_name_.empty() ) {
    // Look up the location of the default configuration file using the
    // FilePropertiesManager to get the directory name
    syst_config_file_name_ = fpm.analysis_path() + "/systcalc.conf";
  }

  // Open in "update" mode so that we can save POT-summed histograms
  // for the combination of all analysis ntuples. Otherwise, we won't
  // write to the file.
  // TODO: consider adjusting this to be less dangerous
  TFile in_tfile( input_respmat_file_name.c_str(), "update" );

  TDirectoryFile* root_tdir = nullptr;

  // If we haven't been handed a root TDirectoryFile name explicitly, then
  // just grab the first key from the input TFile and assume it's the right
  // one to use.
  std::string tdf_name = respmat_tdirectoryfile_name;
  if ( tdf_name.empty() ) {
    tdf_name = in_tfile.GetListOfKeys()->At( 0 )->GetName();
  }

  in_tfile.GetObject( tdf_name.c_str(), root_tdir );
  if ( !root_tdir ) {
    throw std::runtime_error( "Invalid root TDirectoryFile!" );
  }

  // Construct the "total subfolder name" using the constant prefix
  // and the name of the FilePropertiesManager configuration file that
  // is currently active. Replace any '/' characters in the latter to
  // avoid TDirectoryFile path problems.
  std::string fpm_config_file = fpm.config_file_name();
  // Do the '/' replacement here in the same way as is done for
  // TDirectoryFile subfolders by the UniverseMaker class
  fpm_config_file = ntuple_subfolder_from_file_name( fpm_config_file );

  std::string total_subfolder_name = TOTAL_SUBFOLDER_NAME_PREFIX
    + fpm_config_file;

  // Check whether a set of POT-summed histograms for each universe
  // is already present in the input response matrix file. This is
  // signalled by a TDirectoryFile with a name matching the string
  // total_subfolder_name.
  TDirectoryFile* total_subdir = nullptr;
  root_tdir->GetObject( total_subfolder_name.c_str(), total_subdir );

  if ( !total_subdir ) {

    // We couldn't find the pre-computed POT-summed universe histograms,
    // so make them "on the fly" and store them in this object
    this->build_universes( *root_tdir );

    // Create a new TDirectoryFile as a subfolder to hold the POT-summed
    // universe histograms
    total_subdir = new TDirectoryFile( total_subfolder_name.c_str(),
      "universes", "", root_tdir );

    // Write the universes to the new subfolder for faster loading
    // later
    this->save_universes( *total_subdir );
  }
  else {
    // Retrieve the POT-summed universe histograms that were built
    // previously
    this->load_universes( *total_subdir );
  }

  // Also load the configuration of true and reco bins used to create the
  // universes
  std::string* true_bin_spec = nullptr;
  std::string* reco_bin_spec = nullptr;

  root_tdir->GetObject( TRUE_BIN_SPEC_NAME.c_str(), true_bin_spec );
  root_tdir->GetObject( RECO_BIN_SPEC_NAME.c_str(), reco_bin_spec );

  if ( !true_bin_spec || !reco_bin_spec ) {
    throw std::runtime_error( "Failed to load bin specifications" );
  }

  num_signal_true_bins_ = 0u;
  std::istringstream iss_true( *true_bin_spec );
  TrueBin temp_true_bin;
  while ( iss_true >> temp_true_bin ) {
    if ( temp_true_bin.type_ == TrueBinType::kSignalTrueBin ) {
      ++num_signal_true_bins_;
    }
    true_bins_.push_back( temp_true_bin );
  }

  num_ordinary_reco_bins_ = 0u;
  std::istringstream iss_reco( *reco_bin_spec );
  RecoBin temp_reco_bin;
  while ( iss_reco >> temp_reco_bin ) {
    if ( temp_reco_bin.type_ == RecoBinType::kOrdinaryRecoBin ) {
      ++num_ordinary_reco_bins_;
    }
    reco_bins_.push_back( temp_reco_bin );
  }

}

void SystematicsCalculator::load_universes( TDirectoryFile& total_subdir ) {

  const auto& fpm = FilePropertiesManager::Instance();

  // TODO: reduce code duplication between this function
  // and SystematicsCalculator::build_universes()
  TList* universe_key_list = total_subdir.GetListOfKeys();
  int num_keys = universe_key_list->GetEntries();

  // Loop over the keys in the TDirectoryFile. Build a universe object
  // for each key ending in "_2d" and store it in the rw_universes_
  // map (for reweightable systematic universes), the detvar_universes_
  // map (for detector systematic universes), or the alt_cv_universes_
  // map (for alternate CV model universes)
  for ( int k = 0; k < num_keys; ++k ) {
    // To avoid double-counting universes, search only for the 2D event
    // count histograms
    std::string key = universe_key_list->At( k )->GetName();
    bool is_not_2d_hist = !has_ending( key, "_2d" );
    if ( is_not_2d_hist ) continue;

    // Get rid of the trailing "_2d" by deleting the last three
    // characters from the current key
    key.erase( key.length() - 3u );

    // The last underscore separates the universe name from its
    // index. Split the key into these two parts.
    size_t temp_idx = key.find_last_of( '_' );

    std::string univ_name = key.substr( 0, temp_idx );
    std::string univ_index_str = key.substr( temp_idx + 1u );

    int univ_index = std::stoi( univ_index_str );

    TH1D* hist_true = nullptr;
    TH1D* hist_reco = nullptr;
    TH2D* hist_2d = nullptr;
    TH2D* hist_categ = nullptr;
    TH2D* hist_reco2d = nullptr;

    total_subdir.GetObject( (key + "_true").c_str(), hist_true );
    total_subdir.GetObject( (key + "_reco").c_str(), hist_reco );
    total_subdir.GetObject( (key + "_2d").c_str(), hist_2d );
    total_subdir.GetObject( (key + "_categ").c_str(), hist_categ );
    total_subdir.GetObject( (key + "_reco2d").c_str(), hist_reco2d );

    if ( !hist_true || !hist_reco || !hist_2d || !hist_categ || !hist_reco2d ) {
      throw std::runtime_error( "Failed to retrieve histograms for the "
        + key + " universe" );
    }

    // Reconstruct the Universe object from the retrieved histograms
    auto temp_univ = std::make_unique< Universe >( univ_name, univ_index,
      hist_true, hist_reco, hist_2d, hist_categ, hist_reco2d );

    // Determine whether the current universe represents a detector
    // variation or a reweightable variation. We'll use this information to
    // decide where it should be stored.
    NFT temp_type = fpm.string_to_ntuple_type( univ_name );
    if ( temp_type != NFT::kUnknown ) {

      bool is_detvar = ntuple_type_is_detVar( temp_type );
      bool is_altCV = ntuple_type_is_altCV( temp_type );
      if ( !is_detvar && !is_altCV ) throw std::runtime_error( "Universe name "
        + univ_name + " matches a non-detVar and non-altCV file type."
        + " Handling of this situation is currently unimplemented." );

      if ( is_detvar && detvar_universes_.count(temp_type) ) {
        throw std::runtime_error( "detVar multisims are not currently"
          " supported" );
      }
      else if ( is_altCV && alt_cv_universes_.count(temp_type) ) {
        throw std::runtime_error( "altCV multisims are not currently"
          " supported" );
      }

      // Move the detector variation Universe object into the map
      if ( is_detvar ) {
        detvar_universes_[ temp_type ].reset( temp_univ.release() );
      }
      else { // is_altCV
        alt_cv_universes_[ temp_type ].reset( temp_univ.release() );
      }
    }
    // If we're working with fake data, then a single universe with a specific
    // name stores all of the MC information for the "data." Save it in the
    // dedicated fake data Universe object.
    else if ( univ_name == "FakeDataMC" ) {
      std::cout << "******* USING FAKE DATA *******\n";
      fake_data_universe_ = std::move( temp_univ );
    }
    else {
      // If we've made it here, then we're working with a universe
      // for a reweightable systematic variation

      // If we do not already have a map entry for this kind of universe,
      // then create one
      if ( !rw_universes_.count(univ_name) ) {
        rw_universes_[ univ_name ]
          = std::vector< std::unique_ptr<Universe> >();
      }

      // Move this universe into the map. Note that the automatic
      // sorting of keys in a ROOT TDirectoryFile ensures that the
      // universe ordering remains correct. We'll double-check that
      // below, though, just in case.
      auto& univ_vec = rw_universes_.at( univ_name );
      univ_vec.emplace_back( std::move(temp_univ) );

      // Verify that the new universe is placed in the expected
      // position in the vector. If there's a mismatch, something has
      // gone wrong and the universe ordering will not be preserved.
      int vec_index = univ_vec.size() - 1;
      if ( vec_index != univ_index ) {
        throw std::runtime_error( "Universe index mismatch encountered!" );
      }
    }

  } // TDirectoryFile keys (and 2D universe histograms)

  constexpr std::array< NFT, 2 > data_file_types = { NFT::kOnBNB,
    NFT::kExtBNB };

  for ( const auto& file_type : data_file_types ) {
    std::string data_name = fpm.ntuple_type_to_string( file_type );
    std::string hist_name = data_name + "_reco";

    std::string hist2d_name = data_name + "_reco2d";

    TH1D* hist = nullptr;
    total_subdir.GetObject( hist_name.c_str(), hist );

    TH2D* hist2d = nullptr;
    total_subdir.GetObject( hist2d_name.c_str(), hist2d );

    if ( !hist || !hist2d ) {
      throw std::runtime_error( "Missing data histogram for " + data_name );
    }

    hist->SetDirectory( nullptr );
    hist2d->SetDirectory( nullptr );

    if ( data_hists_.count(file_type) || data_hists2d_.count(file_type) ) {
      throw std::runtime_error( "Duplicate data histogram for "
        + data_name );
    }

    data_hists_[ file_type ].reset( hist );
    data_hists2d_[ file_type ].reset( hist2d );
  }

  TParameter< double >* temp_pot = nullptr;
  total_subdir.GetObject( "total_bnb_data_pot", temp_pot );
  if ( !temp_pot ) {
    throw std::runtime_error( "Missing BNB data POT value" );
  }

  total_bnb_data_pot_ = temp_pot->GetVal();

}

void SystematicsCalculator::build_universes( TDirectoryFile& root_tdir ) {

  // Set default values of flags used to signal the presence of fake data. If
  // fake data are detected, corresponding truth information will be stored and
  // a check will be performed to prevent mixing real and fake data together.
  bool using_fake_data = false;
  bool began_checking_for_fake_data = false;

  // Get some normalization factors here for easy use later.
  std::map< int, double > run_to_bnb_pot_map;
  std::map< int, double > run_to_bnb_trigs_map;
  std::map< int, double > run_to_ext_trigs_map;

  const auto& fpm = FilePropertiesManager::Instance();
  const auto& data_norm_map = fpm.data_norm_map();
  for ( const auto& run_and_type_pair : fpm.ntuple_file_map() ) {
    int run = run_and_type_pair.first;
    const auto& type_map = run_and_type_pair.second;

    const auto& bnb_file_set = type_map.at( NFT::kOnBNB );
    for ( const std::string& bnb_file : bnb_file_set ) {
      const auto& pot_and_trigs = data_norm_map.at( bnb_file );

      if ( !run_to_bnb_pot_map.count(run) ) {
        run_to_bnb_pot_map[ run ] = 0.;
        run_to_bnb_trigs_map[ run ] = 0.;
      }

      run_to_bnb_pot_map.at( run ) += pot_and_trigs.pot_;
      run_to_bnb_trigs_map.at( run ) += pot_and_trigs.trigger_count_;

    } // BNB data files

    const auto& ext_file_set = type_map.at( NFT::kExtBNB );
    for ( const std::string& ext_file : ext_file_set ) {
      const auto& pot_and_trigs = data_norm_map.at( ext_file );

      if ( !run_to_ext_trigs_map.count(run) ) {
        run_to_ext_trigs_map[ run ] = 0.;
      }

      run_to_ext_trigs_map.at( run ) += pot_and_trigs.trigger_count_;

    } // EXT files

  } // runs

  // Now that we have the accumulated POT over all BNB data runs, sum it
  // into a single number. This will be used to normalize the detVar MC
  // samples.
  total_bnb_data_pot_ = 0.;
  for ( const auto& pair : run_to_bnb_pot_map ) {
    total_bnb_data_pot_ += pair.second;
  }

  // Loop through the ntuple files for the various run / ntuple file type
  // pairs considered in the analysis. We will react differently in a run-
  // and type-dependent way.
  for ( const auto& run_and_type_pair : fpm.ntuple_file_map() ) {

    int run = run_and_type_pair.first;
    const auto& type_map = run_and_type_pair.second;

    for ( const auto& type_and_files_pair : type_map ) {
      const NFT& type = type_and_files_pair.first;
      const auto& file_set = type_and_files_pair.second;

      bool is_detVar = ntuple_type_is_detVar( type );
      bool is_altCV = ntuple_type_is_altCV( type );
      bool is_reweightable_mc = ntuple_type_is_reweightable_mc( type );
      bool is_mc = ntuple_type_is_mc( type );

      for ( const std::string& file_name : file_set ) {

        std::cout << "PROCESSING universes for " << file_name << '\n';

        // Default to assuming that the current ntuple file is not a fake data
        // sample. If it is a data sample (i.e., if is_mc == false), then the
        // value of this flag will be reconsidered below.
        bool is_fake_data = false;

        // Get the simulated or measured POT belonging to the current file.
        // This will be used to normalize the relevant histograms
        double file_pot = 0.;
        if ( is_mc ) {
          // MC files have the simulated POT stored alongside the ntuple
          // TODO: use the TDirectoryFile to handle this rather than
          // pulling it out of the original ntuple file
          TFile temp_mc_file( file_name.c_str(), "read" );
          TParameter<float>* temp_pot = nullptr;
          temp_mc_file.GetObject( "summed_pot", temp_pot );
          if ( !temp_pot ) throw std::runtime_error(
            "Missing POT in MC file!" );
          file_pot = temp_pot->GetVal();
        }
        else {
          // We can ask the FilePropertiesManager for the data POT values
          file_pot = fpm.data_norm_map().at( file_name ).pot_;
        }

        // Get the TDirectoryFile name used to store histograms for the
        // current ntuple file
        std::string subdir_name = ntuple_subfolder_from_file_name(
          file_name );

        TDirectoryFile* subdir = nullptr;
        root_tdir.GetObject( subdir_name.c_str(), subdir );
        if ( !subdir ) throw std::runtime_error(
          "Missing TDirectoryFile " + subdir_name );

        // For data, just add the reco-space event counts to the total,
        // scaling to the beam-on triggers in the case of EXT data
        if ( !is_mc ) {

          auto reco_hist = get_object_unique_ptr< TH1D >(
            "unweighted_0_reco", *subdir );

          auto reco_hist2d = get_object_unique_ptr< TH2D >(
            "unweighted_0_reco2d", *subdir );

          // If we're working with EXT data, scale it to the corresponding
          // number of triggers from the BNB data from the same run
          if ( type == NFT::kExtBNB ) {
            double bnb_trigs = run_to_bnb_trigs_map.at( run );
            double ext_trigs = run_to_ext_trigs_map.at( run );

            reco_hist->Scale( bnb_trigs / ext_trigs );
            reco_hist2d->Scale( bnb_trigs / ext_trigs );
          }

          // If we don't have a histogram in the map for this data type
          // yet, just clone the existing histogram.
          if ( !data_hists_.count(type) ) {
            TH1D* temp_clone = dynamic_cast<TH1D*>(
              reco_hist->Clone("temp_clone")
            );
            temp_clone->SetStats( false );
            temp_clone->SetDirectory( nullptr );
            // Note: here the map entry takes ownership of the histogram
            data_hists_[ type ].reset( temp_clone );
          }
          else {
            // Otherwise, just add its contribution to the existing
            // histogram
            data_hists_.at( type )->Add( reco_hist.get() );
          }

          if ( !data_hists2d_.count(type) ) {
            TH2D* temp_clone = dynamic_cast<TH2D*>(
              reco_hist2d->Clone("temp_clone")
            );
            temp_clone->SetStats( false );
            temp_clone->SetDirectory( nullptr );
            // Note: here the map entry takes ownership of the histogram
            data_hists2d_[ type ].reset( temp_clone );
          }
          else {
            // Otherwise, just add its contribution to the existing
            // histogram
            data_hists2d_.at( type )->Add( reco_hist2d.get() );
          }

          // For EXT data files (always assumed to be real data), no further
          // processing is needed, so just move to the next file
          if ( type == NFT::kExtBNB ) continue;

          // For real BNB data, we're also done. However, if we're working with
          // fake data, then we also want to store the truth information for
          // the current ntuple file. Check to see whether we have fake data or
          // not by trying to retrieve the unweighted "2d" histogram (which
          // stores event counts that simultaneously fall into a given true bin
          // and reco bin). It will only be present for fake data ntuples.

          auto fd_check_2d_hist = get_object_unique_ptr< TH2D >(
            "unweighted_0_2d", *subdir );

          is_fake_data = ( fd_check_2d_hist.get() != nullptr );

          // The BNB data files should either all be real data or all be fake
          // data. If any mixture occurs, then this suggests that something has
          // gone seriously wrong. Throw an exception to warn the user.
          if ( began_checking_for_fake_data ) {
            if ( is_fake_data != using_fake_data ) {
              throw std::runtime_error( "Mixed real and fake data ntuple files"
                " encountered in SystematicsCalculator::build_universes()" );
            }
          }
          else {
            using_fake_data = is_fake_data;
            began_checking_for_fake_data = true;
          }

          // If we are working with real BNB data, then we don't need to do the
          // other MC-ntuple-specific stuff below, so just move on to the next
          // file
          if ( !is_fake_data ) continue;

        } // data ntuple files

        // If we've made it here, then we're working with an MC ntuple
        // file. For these, all four histograms for the "unweighted"
        // universe are always evaluated. Use this to determine the number
        // of true and reco bins easily.
        auto temp_2d_hist = get_object_unique_ptr< TH2D >(
          "unweighted_0_2d", *subdir );

        // NOTE: the convention of the UniverseMaker class is to use
        // x as the true axis and y as the reco axis.
        int num_true_bins = temp_2d_hist->GetXaxis()->GetNbins();
        int num_reco_bins = temp_2d_hist->GetYaxis()->GetNbins();

        // Let's handle the fake BNB data samples first.
        if ( is_fake_data ) {

          // If this is our first fake BNB data ntuple file, then create
          // the Universe object that will store the full MC information
          // (truth and reco)
          if ( !fake_data_universe_ ) {

            fake_data_universe_ = std::make_unique< Universe >( "FakeDataMC",
              0, num_true_bins, num_reco_bins );

            set_stats_and_dir( *fake_data_universe_ );

          } // first fake BNB data ntuple file


          // Since all other histograms are scaled to the POT exposure
          // of the BNB data ones, no rescaling is needed here. Just add the
          // contributions of the current ntuple file's histograms to the
          // fake data Universe object.

          // Retrieve the histograms for the fake data universe (which
          // corresponds to the "unweighted" histogram name prefix) from
          // the current TDirectoryFile.
          // TODO: add the capability for fake data to use event weights
          // (both in the saved fake data Universe object and in the
          // corresponding "data" histogram of reco event counts
          std::string hist_name_prefix = "unweighted_0";

          auto h_reco = get_object_unique_ptr< TH1D >(
            (hist_name_prefix + "_reco"), *subdir );

          auto h_true = get_object_unique_ptr< TH1D >(
            (hist_name_prefix + "_true"), *subdir );

          auto h_2d = get_object_unique_ptr< TH2D >(
            (hist_name_prefix + "_2d"), *subdir );

          auto h_categ = get_object_unique_ptr< TH2D >(
            (hist_name_prefix + "_categ"), *subdir );

          auto h_reco2d = get_object_unique_ptr< TH2D >(
            (hist_name_prefix + "_reco2d"), *subdir );

          // Add their contributions to the owned histograms for the
          // current Universe object
          fake_data_universe_->hist_reco_->Add( h_reco.get() );
          fake_data_universe_->hist_true_->Add( h_true.get() );
          fake_data_universe_->hist_2d_->Add( h_2d.get() );
          fake_data_universe_->hist_categ_->Add( h_categ.get() );
          fake_data_universe_->hist_reco2d_->Add( h_reco2d.get() );

        } // fake data sample

        // Now we'll take care of the detVar and altCV samples.
        else if ( is_detVar || is_altCV ) {

          std::string dv_univ_name = fpm.ntuple_type_to_string( type );

          // Make a temporary new Universe object to store
          // (POT-scaled) detVar/altCV histograms (if needed)
          auto temp_univ = std::make_unique< Universe >( dv_univ_name,
            0, num_true_bins, num_reco_bins );

          // Temporary pointer that will allow us to treat the single-file
          // detector variation samples and the multi-file alternate CV samples
          // on the same footing below
          Universe* temp_univ_ptr = nullptr;

          // Check whether a prior altCV Universe exists in the map
          bool prior_altCV = alt_cv_universes_.count( type ) > 0;

          // Right now, we're assuming there's just one detVar ntuple file
          // per universe. If this assumption is violated, our scaling will
          // be screwed up. Prevent this from happening by throwing an
          // exception when a duplicate is encountered.
          // TODO: revisit this when you have detVar samples for all runs
          if ( is_detVar && detvar_universes_.count(type) ) {
            throw std::runtime_error( "Duplicate detVar ntuple file!" );
          }
          // For the alternate CV sample, if a previous universe already
          // exists in the map, then get access to it via a pointer
          else if ( is_altCV && prior_altCV ) {
            temp_univ_ptr = alt_cv_universes_.at( type ).get();
          }
          else {
            temp_univ_ptr = temp_univ.get();
          }

          // The detector variation and alternate CV universes are all labeled
          // "unweighted" in the response matrix files. Retrieve the
          // corresponding histograms.
          // TODO: revisit this assumption as appropriate
          auto hist_reco = get_object_unique_ptr< TH1D >(
            "unweighted_0_reco", *subdir );

          auto hist_true = get_object_unique_ptr< TH1D >(
            "unweighted_0_true", *subdir );

          auto hist_2d = get_object_unique_ptr< TH2D >(
            "unweighted_0_2d", *subdir );

          auto hist_categ = get_object_unique_ptr< TH2D >(
            "unweighted_0_categ", *subdir );

          auto hist_reco2d = get_object_unique_ptr< TH2D >(
            "unweighted_0_reco2d", *subdir );

          double temp_scale_factor = 1.;
          if ( is_altCV ) {
            // AltCV ntuple files are available for all runs, so scale
            // each individually to the BNB data POT for the current run
            double temp_run_pot = run_to_bnb_pot_map.at( run );
            temp_scale_factor = temp_run_pot / file_pot;
          }
          else {
            // Scale all detVar universe histograms from the simulated POT to
            // the *total* BNB data POT for all runs analyzed. Since we only
            // have detVar samples for Run 3b, we assume that they can be
            // applied globally in this step.
            // TODO: revisit this as appropriate
            temp_scale_factor = total_bnb_data_pot_ / file_pot;
          }

          // Apply the scaling factor defined above to all histograms that
          // will be owned by the new Universe
          hist_reco->Scale( temp_scale_factor );
          hist_true->Scale( temp_scale_factor );
          hist_2d->Scale( temp_scale_factor );
          hist_categ->Scale( temp_scale_factor );
          hist_reco2d->Scale( temp_scale_factor );

          // Add the scaled contents of these histograms to the
          // corresponding histograms in the new Universe object
          temp_univ_ptr->hist_reco_->Add( hist_reco.get() );
          temp_univ_ptr->hist_true_->Add( hist_true.get() );
          temp_univ_ptr->hist_2d_->Add( hist_2d.get() );
          temp_univ_ptr->hist_categ_->Add( hist_categ.get() );
          temp_univ_ptr->hist_reco2d_->Add( hist_reco2d.get() );

          // Adjust the owned histograms to avoid auto-deletion problems
          set_stats_and_dir( *temp_univ_ptr );

          // If one wasn't present before, then move the finished Universe
          // object into the map
          if ( is_detVar ) {
            detvar_universes_[ type ].reset( temp_univ.release() );
          }
          else if ( !prior_altCV ) { // is_altCV
            alt_cv_universes_[ type ].reset( temp_univ.release() );
          }

        } // detVar and altCV samples

        // Now handle the reweightable systematic universes
        else if ( is_reweightable_mc ) {

          // If this is our first reweightable MC ntuple file, then build
          // the map of reweighting universes from the 2D histogram keys in
          // its TDirectoryFile.
          // NOTE: I rely here on the reweighting universe definitions
          // being identical across all ntuples considered by the script.
          if ( rw_universes_.empty() ) {

            TList* universe_key_list = subdir->GetListOfKeys();
            int num_keys = universe_key_list->GetEntries();

            for ( int k = 0; k < num_keys; ++k ) {
              // To avoid double-counting universes, only create new
              // Universe objects for the 2D event count histograms
              std::string key = universe_key_list->At( k )->GetName();
              bool is_not_2d_hist = !has_ending( key, "_2d" );
              if ( is_not_2d_hist ) continue;

              // Get rid of the trailing "_2d" by deleting the last three
              // characters from the current key
              key.erase( key.length() - 3u );

              // The last underscore separates the universe name from its
              // index. Split the key into these two parts.
              size_t temp_idx = key.find_last_of( '_' );

              std::string univ_name = key.substr( 0, temp_idx );
              std::string univ_index_str = key.substr( temp_idx + 1u );

              int univ_index = std::stoi( univ_index_str );

              // We have what we need to create the new Universe object. Do
              // it! Note that its owned histograms are currently empty.
              auto temp_univ = std::make_unique<Universe>( univ_name,
                univ_index, num_true_bins, num_reco_bins );

              set_stats_and_dir( *temp_univ );

              // If we do not already have a map entry for this kind of
              // universe, then create one
              if ( !rw_universes_.count(univ_name) ) {
                rw_universes_[ univ_name ]
                  = std::vector< std::unique_ptr<Universe> >();
              }

              // Move this universe into the map. Note that the automatic
              // sorting of keys in a ROOT TDirectoryFile ensures that the
              // universe ordering remains correct.
              rw_universes_.at( univ_name ).emplace_back(
                std::move(temp_univ) );

            } // TDirectoryFile keys

          } // first reweightable MC ntuple file


          // For reweightable MC ntuple files, scale the histograms for
          // each universe to the BNB data POT for the current run
          double run_bnb_pot = run_to_bnb_pot_map.at( run );
          double rw_scale_factor = run_bnb_pot / file_pot;

          // Iterate over the reweighting universes, retrieve the
          // histograms for each, and add their POT-scaled contributions
          // from the current ntuple file to the total
          for ( auto& rw_pair : rw_universes_ ) {
            std::string univ_name = rw_pair.first;
            auto& univ_vec = rw_pair.second;

            for ( size_t u_idx = 0u; u_idx < univ_vec.size(); ++u_idx ) {
              // Get a reference to the current universe object
              auto& universe = *univ_vec.at( u_idx );

              // Double-check that the universe ordering is right. The
              // index in the map of universes should match the index
              // stored in the Universe object itself. If this check fails,
              // something went wrong with the key sorting imposed by the
              // input TDirectoryFile.
              if ( u_idx != universe.index_ ) throw std::runtime_error(
                "Universe sorting went wrong!" );

              // Retrieve the histograms for the current universe from the
              // current TDirectoryFile
              std::string hist_name_prefix = univ_name
                + '_' + std::to_string( u_idx );

              auto h_reco = get_object_unique_ptr< TH1D >(
                (hist_name_prefix + "_reco"), *subdir );

              auto h_true = get_object_unique_ptr< TH1D >(
                (hist_name_prefix + "_true"), *subdir );

              auto h_2d = get_object_unique_ptr< TH2D >(
                (hist_name_prefix + "_2d"), *subdir );

              auto h_categ = get_object_unique_ptr< TH2D >(
                (hist_name_prefix + "_categ"), *subdir );

              auto h_reco2d = get_object_unique_ptr< TH2D >(
                (hist_name_prefix + "_reco2d"), *subdir );

              // Scale these histograms to the appropriate BNB data POT for
              // the current run
              h_reco->Scale( rw_scale_factor );
              h_true->Scale( rw_scale_factor );
              h_2d->Scale( rw_scale_factor );
              h_categ->Scale( rw_scale_factor );
              h_reco2d->Scale( rw_scale_factor );

              // Add their contributions to the owned histograms for the
              // current Universe object
              universe.hist_reco_->Add( h_reco.get() );
              universe.hist_true_->Add( h_true.get() );
              universe.hist_2d_->Add( h_2d.get() );
              universe.hist_categ_->Add( h_categ.get() );
              universe.hist_reco2d_->Add( h_reco2d.get() );

            } // universes indices

          } // universe types

        } // reweightable MC samples

      } // ntuple file

    } // type

  } // run

  // Everything is ready to go with one possible exception: if we're working
  // with fake data ntuples, then the "data" histograms of reco-space event
  // counts will contain the "beam-on" (MC) contribution but not the "beam-off"
  // (EXT) contribution needed to form a full fake measurement. Correct this
  // situation if needed by adding the total EXT histograms to the BNB "data"
  // histograms.
  if ( using_fake_data ) {

    const TH1D* ext_hist = data_hists_.at( NFT::kExtBNB ).get(); // EXT data
    TH1D* bnb_hist = data_hists_.at( NFT::kOnBNB ).get();

    bnb_hist->Add( ext_hist );

    const TH2D* ext_hist2d = data_hists2d_.at( NFT::kExtBNB ).get(); // EXT data
    TH2D* bnb_hist2d = data_hists2d_.at( NFT::kOnBNB ).get();

    bnb_hist2d->Add( ext_hist2d );

    // In a real measurement, the (Poisson) statistical uncertainty on each BNB
    // data histogram bin would be simply the square root of the event count.
    // Imitate a real measurement by replacing the existing BNB "data"
    // histogram bin errors (based on MC and EXT data statistics) with the
    // square root of the bin contents. Note that no scaling subtleties come
    // into play because the BNB data histograms are never rescaled based upon
    // POT exposure. Include the under/overflow bins when reassigning the error
    // bars, just in case.
    // TODO: consider whether a different approach to assigning the fake data
    // statistical uncertainties is more appropriate
    int num_reco_bins = bnb_hist->GetNbinsX();
    for ( int rb = 0; rb <= num_reco_bins + 1; ++rb ) {
      double evts = std::max( 0., bnb_hist->GetBinContent(rb) );
      bnb_hist->SetBinError( rb, std::sqrt(evts) );

      for ( int rb2 = 0; rb2 <= num_reco_bins + 1; ++rb2 ) {
        double evts2d = std::max( 0., bnb_hist2d->GetBinContent(rb, rb2) );
        bnb_hist2d->SetBinError( rb, rb2, std::sqrt(evts2d) );
      }

    }

    std::cout << "******* USING FAKE DATA *******\n";
  }

  std::cout << "TOTAL BNB DATA POT = " << total_bnb_data_pot_ << '\n';
}

void SystematicsCalculator::save_universes( TDirectoryFile& out_tdf ) {

  // Make the requested output TDirectoryFile the active one
  out_tdf.cd();

  // Save the data histograms
  for ( const auto& pair : data_hists_ ) {

    NFT type = pair.first;
    auto& hist = pair.second;

    const auto& fpm = FilePropertiesManager::Instance();
    std::string univ_name = fpm.ntuple_type_to_string( type );
    hist->Write( (univ_name + "_reco").c_str() );
  }

  for ( const auto& pair : data_hists2d_ ) {

    NFT type = pair.first;
    auto& hist2d = pair.second;

    const auto& fpm = FilePropertiesManager::Instance();
    std::string univ_name = fpm.ntuple_type_to_string( type );
    hist2d->Write( (univ_name + "_reco2d").c_str() );
  }

  // Save the detector variation histograms
  for ( const auto& pair : detvar_universes_ ) {
    NFT type = pair.first;
    auto& universe = pair.second;

    universe->hist_reco_->Write();
    universe->hist_true_->Write();
    universe->hist_2d_->Write();
    universe->hist_categ_->Write();
    universe->hist_reco2d_->Write();
  }

  // Save the alternate CV MC histograms
  for ( const auto& pair : alt_cv_universes_ ) {
    NFT type = pair.first;
    auto& universe = pair.second;

    universe->hist_reco_->Write();
    universe->hist_true_->Write();
    universe->hist_2d_->Write();
    universe->hist_categ_->Write();
    universe->hist_reco2d_->Write();
  }

  // Save the reweightable systematic histograms
  for ( const auto& pair : rw_universes_ ) {

    const auto& univ_vec = pair.second;
    for ( const auto& universe : univ_vec ) {
      universe->hist_reco_->Write();
      universe->hist_true_->Write();
      universe->hist_2d_->Write();
      universe->hist_categ_->Write();
      universe->hist_reco2d_->Write();
    }

  }

  // Save the fake data universe histograms if they exist
  if ( fake_data_universe_ ) {
    fake_data_universe_->hist_reco_->Write();
    fake_data_universe_->hist_true_->Write();
    fake_data_universe_->hist_2d_->Write();
    fake_data_universe_->hist_categ_->Write();
    fake_data_universe_->hist_reco2d_->Write();
  }

  // Save the total BNB data POT for easy retrieval later
  TParameter< double > temp_pot( "total_bnb_data_pot",
    total_bnb_data_pot_ );

  temp_pot.Write();

}

CovMatrix SystematicsCalculator::make_covariance_matrix(
  const std::string& hist_name ) const
{
  int num_cm_bins = this->get_covariance_matrix_size();
  TH2D* hist = new TH2D( hist_name.c_str(),
    "covariance; reco bin; reco bin; covariance", num_cm_bins, 0.,
    num_cm_bins, num_cm_bins, 0., num_cm_bins );
  hist->SetDirectory( nullptr );
  hist->SetStats( false );

  CovMatrix result( hist );
  return result;
}

template < class UniversePointerContainer >
  void make_cov_mat( const SystematicsCalculator& sc, CovMatrix& cov_mat,
  const Universe& cv_univ, const UniversePointerContainer& universes,
  bool average_over_universes, bool is_flux_variation )
{
  // Get the total number of true bins and the covariance matrix dimension for
  // later reference
  size_t num_true_bins = sc.true_bins_.size();
  size_t num_cm_bins = sc.get_covariance_matrix_size();

  // Get the expected observable values in each reco bin in the CV universe
  std::vector< double > cv_reco_obs( num_cm_bins, 0. );

  for ( size_t rb = 0u; rb < num_cm_bins; ++rb ) {
    cv_reco_obs.at( rb ) = sc.evaluate_observable( cv_univ, rb );
  }

  // Loop over universes
  int num_universes = universes.size();
  for ( int u_idx = 0; u_idx < num_universes; ++u_idx ) {

    const auto& univ = universes.at( u_idx );

    // For flux variations, we also want to be able to adjust the integrated
    // flux used to compute the differential cross section. This is signaled by
    // passing a nonnegative index for a flux systematic universe to
    // SystematicsCalculator::evaluate_observable()
    int flux_u_idx = -1;
    if ( is_flux_variation ) flux_u_idx = u_idx;

    // Get the expected observable values in each reco bin in the
    // current universe.
    std::vector< double > univ_reco_obs( num_cm_bins, 0. );

    for ( size_t rb = 0u; rb < num_cm_bins; ++rb ) {
      univ_reco_obs.at( rb ) = sc.evaluate_observable( *univ,
        rb, flux_u_idx );
    }

    // We have all the needed ingredients to get the contribution of this
    // universe to the covariance matrix. Loop over each pair covariance matrix
    // elements and fill them.
    // TODO: the covariance matrix are symmetric by definition. You can
    // therefore make this more efficient by calculating only the subset of
    // elements that you need.
    for ( size_t a = 0u; a < num_cm_bins; ++a ) {

      double cv_a = cv_reco_obs.at( a );
      double univ_a = univ_reco_obs.at( a );

      for ( size_t b = 0u; b < num_cm_bins; ++b ) {

        double cv_b = cv_reco_obs.at( b );
        double univ_b = univ_reco_obs.at( b );

        double covariance  = ( cv_a - univ_a ) * ( cv_b - univ_b );

        // We cheat here by noting that the lower bound of each covariance
        // matrix TH2D bin is the bin index. Filling using the zero-based bin
        // indices and the covariance as the weight yields the desired behavior
        // (increment the existing element by the current covariance value) in
        // an easy-to-read (if slightly evil) way.
        cov_mat.cov_matrix_->Fill( a, b, covariance );
      } // reco bin index b
    } // reco bin index a

  } // universe

  // If requested, average the final covariance matrix elements over all
  // universes
  if ( average_over_universes ) {
    cov_mat.cov_matrix_->Scale( 1. / num_universes );
  }

}

// Overloaded version that takes a single alternate universe wrapped in a
// std::unique_ptr
void make_cov_mat( const SystematicsCalculator& sc, CovMatrix& cov_mat,
  const Universe& cv_univ,
  const Universe& alt_univ, bool average_over_universes = false,
  bool is_flux_variation = false )
{
  std::vector< const Universe* > temp_univ_vec;

  temp_univ_vec.emplace_back( &alt_univ );

  make_cov_mat( sc, cov_mat, cv_univ, temp_univ_vec, average_over_universes,
    is_flux_variation );
}

std::unique_ptr< CovMatrixMap > SystematicsCalculator::get_covariances() const
{
  // Make an empty map to store the covariance matrices
  auto matrix_map_ptr = std::make_unique< CovMatrixMap >();
  auto& matrix_map = *matrix_map_ptr;

  // Read in the definition of each covariance matrix and calculate it. Each
  // definition contains at least a name and a type specifier
  std::ifstream config_file( syst_config_file_name_ );
  std::string name, type;
  while ( config_file >> name >> type ) {

    CovMatrix temp_cov_mat = this->make_covariance_matrix( name );

    // If the current covariance matrix is defined as a sum of others, then
    // just add the existing ones together to compute it
    if ( type == "sum" ) {
      int count = 0;
      config_file >> count;
      std::string cm_name;
      for ( int cm = 0; cm < count; ++cm ) {
        config_file >> cm_name;
        if ( !matrix_map.count(cm_name) ) {
          throw std::runtime_error( "Undefined covariance matrix " + cm_name );
        }
        temp_cov_mat += matrix_map.at( cm_name );
      } // terms in the sum

    } // sum type

    else if ( type == "MCstat" ) {

      size_t num_cm_bins = this->get_covariance_matrix_size();
      const auto& cv_univ = this->cv_universe();

      for ( size_t rb1 = 0u; rb1 < num_cm_bins; ++rb1 ) {
        for ( size_t rb2 = 0u; rb2 < num_cm_bins; ++rb2 ) {
          // Calculate the MC statistical covariance for the current pair of
          // reco bins for the observable of interest. Use the CV universe.
          double mc_cov = this->evaluate_mc_stat_covariance( cv_univ,
            rb1, rb2 );

          // Note the one-based reco bin index used by ROOT histograms
          temp_cov_mat.cov_matrix_->SetBinContent( rb1 + 1, rb2 + 1, mc_cov );
        }
      }

    } // MCstat type

    else if ( type == "BNBstat" || type == "EXTstat" ) {

      bool use_ext = false;
      if ( type == "EXTstat" ) use_ext = true;

      size_t num_cm_bins = this->get_covariance_matrix_size();

      for ( size_t rb1 = 0u; rb1 < num_cm_bins; ++rb1 ) {
        for ( size_t rb2 = 0u; rb2 < num_cm_bins; ++rb2 ) {
          double stat_cov = this->evaluate_data_stat_covariance( rb1,
            rb2, use_ext );

          // Note the one-based reco bin index used by ROOT histograms
          temp_cov_mat.cov_matrix_->SetBinContent( rb1 + 1, rb2 + 1, stat_cov );
        }
      }

    } // BNBstat and EXTstat types

    else if ( type == "MCFullCorr" ) {
      // Read in the fractional uncertainty from the configuration file
      double frac_unc = 0.;
      config_file >> frac_unc;

      const double frac2 = std::pow( frac_unc, 2 );
      int num_cm_bins = this->get_covariance_matrix_size();

      const auto& cv_univ = this->cv_universe();
      for ( size_t a = 0u; a < num_cm_bins; ++a ) {

        double cv_a = this->evaluate_observable( cv_univ, a );

        for ( int b = 0u; b < num_cm_bins; ++b ) {

          double cv_b = this->evaluate_observable( cv_univ, b );

          double covariance = cv_a * cv_b * frac2;

          temp_cov_mat.cov_matrix_->SetBinContent( a + 1, b + 1, covariance );

        } // reco bin b

      } // reco bin a

    } // MCFullCorr type

    else if ( type == "DV" ) {
      // Get the detector variation type represented by the current universe
      std::string ntuple_type_str;
      config_file >> ntuple_type_str;

      const auto& fpm = FilePropertiesManager::Instance();
      auto ntuple_type = fpm.string_to_ntuple_type( ntuple_type_str );

      // Check that it's valid. If not, then complain.
      bool is_not_detVar = !ntuple_type_is_detVar( ntuple_type );
      if ( is_not_detVar ) {
        throw std::runtime_error( "Invalid NtupleFileType!" );
      }

      // Use a bare pointer for the CV universe so that we can reassign it
      // below if needed. References can't be reassigned after they are
      // initialized.
      const auto* detVar_cv_u = detvar_universes_.at( NFT::kDetVarMCCV ).get();
      const auto& detVar_alt_u = detvar_universes_.at( ntuple_type );

      // The Recomb2 and SCE variations use an alternate "extra CV" universe
      // since they were generated with smaller MC statistics.
      // TODO: revisit this if your detVar samples change in the future
      if ( ntuple_type == NFT::kDetVarMCSCE
        || ntuple_type == NFT::kDetVarMCRecomb2 )
      {
        detVar_cv_u = detvar_universes_.at( NFT::kDetVarMCCVExtra ).get();
      }

      make_cov_mat( *this, temp_cov_mat, *detVar_cv_u,
        *detVar_alt_u, false, false );
    } // DV type

    else if ( type == "RW" || type == "FluxRW" ) {

      // Treat flux variations in a special way by setting a flag
      bool is_flux_variation = false;
      if ( type == "FluxRW" ) is_flux_variation = true;

      // Get the key to use when looking up weights in the map of reweightable
      // systematic variation universes
      std::string weight_key;
      config_file >> weight_key;

      // Retrieve the vector of universes
      auto end = rw_universes_.cend();
      auto iter = rw_universes_.find( weight_key );
      if ( iter == end ) {
        throw std::runtime_error( "Missing weight key " + weight_key );
      }
      const auto& alt_univ_vec = iter->second;

      // Also read in the flag for whether we should average over universes
      // or not for the current covariance matrix
      bool avg_over_universes = false;
      config_file >> avg_over_universes;

      const auto& cv_univ = this->cv_universe();
      make_cov_mat( *this, temp_cov_mat, cv_univ, alt_univ_vec,
        avg_over_universes, is_flux_variation );

    } // RW and FluxRW types

    else if ( type == "AltUniv" ) {

      std::vector< const Universe* > alt_univ_vec;
      for ( const auto& univ_pair : alt_cv_universes_ ) {
        const auto* univ_ptr = univ_pair.second.get();
        alt_univ_vec.push_back( univ_ptr );
      }

      const auto& cv_univ = this->cv_universe();
      make_cov_mat( *this, temp_cov_mat, cv_univ, alt_univ_vec,
        true, false );
    }

    // Complain if we don't know how to calculate the requested covariance
    // matrix
    else throw std::runtime_error( "Unrecognized covariance matrix type \""
      + type + '\"' );

    // Add the finished covariance matrix to the map. If an entry already
    // exists with the same name, throw an exception.
    if ( matrix_map.count(name) ) {
      throw std::runtime_error( "Duplicate covariance matrix definition for "
        + name );
    }
    matrix_map[ name ] = std::move( temp_cov_mat );

  } // Covariance matrix definitions

  return matrix_map_ptr;
}

// Helper function that searches through the owned map of detector variation
// universes. If the input Universe shares the same memory address as one
// of the Universe objects stored in the map, then this function returns
// true. Otherwise, false is returned. This may be used to detect whether the
// detVarCV is needed in SystematicsCalculator::evaluate_observable().
//
// TODO: This is a bit of a hack. Revisit the technique used, particularly
// if it slows down calculations with this class significantly in the future.
bool SystematicsCalculator::is_detvar_universe( const Universe& univ ) const {

  using DetVarMapElement = std::pair< const NFT, std::unique_ptr<Universe> >;

  const auto begin = detvar_universes_.cbegin();
  const auto end = detvar_universes_.cend();

  const auto iter = std::find_if( begin, end,
    [ &univ ]( const DetVarMapElement& my_pair )
      -> bool { return my_pair.second.get() == &univ; }
  );

  bool found_detvar_universe = ( iter != end );
  return found_detvar_universe;
}

// Returns the smearceptance matrix in the requested Universe
// NOTE: This function assumes that all "ordinary" reco bins are listed before
// the sideband ones and that all signal true bins are listed before the
// background ones.
std::unique_ptr< TMatrixD >
  SystematicsCalculator::get_smearceptance_matrix( const Universe& univ ) const
{
  // The smearceptance matrix definition used here uses the reco bin as the row
  // index and the true bin as the column index. This ensures that multiplying
  // a column vector of true event counts by the matrix will yield a column
  // vector of reco event counts.
  auto smearcept = std::make_unique< TMatrixD >( num_ordinary_reco_bins_,
    num_signal_true_bins_ );

  for ( size_t r = 0u; r < num_ordinary_reco_bins_; ++r ) {
    for ( size_t t = 0u; t < num_signal_true_bins_; ++t ) {
      // Get the numerator and denominator of the smearceptance matrix element.
      // Note that we need to switch to one-based bin indices here to retrieve
      // the information from the ROOT histograms stored in the Universe object.
      double numer = univ.hist_2d_->GetBinContent( t + 1, r + 1 );
      double denom = univ.hist_true_->GetBinContent( t + 1 );
      // Store the result in the smearceptance matrix. Avoid dividing by zero
      // by setting any element with zero denominator to zero.
      double temp_element = 0.;
      if ( denom != 0. ) temp_element = numer / denom;
      smearcept->operator()( r, t ) = temp_element;
    }
  }

  return smearcept;
}

// Returns the event counts in each signal true bin in the central-value
// Universe.
// NOTE: This function assumes that all signal true bins are listed before
// the background ones.
std::unique_ptr< TMatrixD > SystematicsCalculator::get_cv_true_signal() const
{
  const auto& cv_univ = this->cv_universe();

  // The output TMatrixD is a column vector containing the event counts in each
  // signal true bin.
  auto result = std::make_unique< TMatrixD >( num_signal_true_bins_, 1 );

  for ( size_t t = 0u; t < num_signal_true_bins_; ++t ) {
    // Note that we need to switch to one-based bin indices here to retrieve
    // the event counts from the ROOT histogram stored in the Universe object.
    double temp_element = cv_univ.hist_true_->GetBinContent( t + 1 );
    result->operator()( t, 0 ) = temp_element;
  }

  return result;
}

// Returns the expected background in each ordinary reco bin (including both
// EXT and the central-value MC prediction for beam-correlated backgrounds)
// NOTE: This function assumes that all "ordinary" reco bins are listed before
// the sideband ones.
std::unique_ptr< TMatrixD >
  SystematicsCalculator::get_cv_ordinary_reco_helper( bool return_bkgd ) const
{
  int num_true_bins = true_bins_.size();
  const auto& cv_univ = this->cv_universe();
  const TH1D* ext_hist = data_hists_.at( NFT::kExtBNB ).get(); // EXT data

  auto result = std::make_unique< TMatrixD >( num_ordinary_reco_bins_, 1 );

  for ( int r = 0; r < num_ordinary_reco_bins_; ++r ) {

    // Start by tallying the EXT contribution in the current reco bin. Note
    // that the EXT data histogram bins have a one-based index.
    double bkgd_events = ext_hist->GetBinContent( r + 1 );

    // Also start out with zero signal events (the signal prediction comes
    // purely from MC)
    double signal_events = 0.;

    for ( int t = 0; t < num_true_bins; ++t ) {
      auto& tbin = true_bins_.at( t );

      // Make a temporary pointer that will be used to increment
      // either the signal or background event tally based on the
      // interpretation of the bin contents below
      double* temp_events_ptr = nullptr;

      // Set the pointer to the appropriate counter
      if ( tbin.type_ == kBackgroundTrueBin ) {
        temp_events_ptr = &bkgd_events;
      }
      else if ( tbin.type_ == kSignalTrueBin ) {
        temp_events_ptr = &signal_events;
      }
      else {
         throw std::runtime_error( "Bad true bin type in"
           " SystematicsCalculator::get_cv_ordinary_reco_helper()" );
      }

      // Tally the contribution from this true bin. Note that we need one-based
      // bin indices to retrieve this information from the TH2D owned by the
      // Universe object
      ( *temp_events_ptr ) += cv_univ.hist_2d_->GetBinContent( t + 1, r + 1 );
    }

    // We've looped through all of the true bins. Now assign the appropriate
    // event count to the vector of results
    result->operator()( r, 0 ) = return_bkgd ? bkgd_events : signal_events;
  }

  return result;
}

MeasuredEvents SystematicsCalculator::get_measured_events() const
{
  // First retrieve the central-value background prediction
  // (EXT + beam-correlated MC background) for all ordinary reco bins
  TMatrixD* ext_plus_mc_bkgd = this->get_cv_ordinary_reco_bkgd().release();

  // Also get the central-value signal prediction for all ordinary reco bins
  auto mc_signal = this->get_cv_ordinary_reco_signal();

  // Get the total covariance matrix on the reco-space EXT+MC prediction
  // (this will not change after subtraction of the central-value background)
  auto cov_map_ptr = this->get_covariances();
  auto temp_cov = cov_map_ptr->at( "total" ).get_matrix();

  // Extract just the covariance matrix block that describes the ordinary reco
  // bins
  TMatrixD* cov_mat = new TMatrixD(
    temp_cov->GetSub( 0, num_ordinary_reco_bins_ - 1,
      0, num_ordinary_reco_bins_ - 1 )
  );

  // Create the vector of measured event counts in the ordinary reco bins
  const TH1D* d_hist = data_hists_.at( NFT::kOnBNB ).get(); // BNB data
  TMatrixD ordinary_data( num_ordinary_reco_bins_, 1 );
  for ( int r = 0; r < num_ordinary_reco_bins_; ++r ) {
    // Switch to using the one-based TH1D index when retrieving these values
    double bnb_events = d_hist->GetBinContent( r + 1 );

    ordinary_data( r, 0 ) = bnb_events;
  }

  // Get the ordinary reco bin data column vector after subtracting the
  // central-value EXT+MC background prediction
  auto* reco_data_minus_bkgd = new TMatrixD( ordinary_data,
    TMatrixD::EMatrixCreatorsOp2::kMinus, *ext_plus_mc_bkgd );

  // Also get the total central-value EXT+MC prediction in each reco bin
  auto* ext_plus_mc_total = new TMatrixD( *mc_signal,
    TMatrixD::EMatrixCreatorsOp2::kPlus, *ext_plus_mc_bkgd );

  MeasuredEvents result( reco_data_minus_bkgd, ext_plus_mc_bkgd,
    ext_plus_mc_total, cov_mat );
  return result;
}
