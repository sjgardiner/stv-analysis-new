#pragma once

// ROOT includes
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

// STV analysis includes
#include "FilePropertiesManager.hh"
#include "IntegratedFluxUniverseManager.hh"
#include "ResponseMatrixMaker.hh"

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

  std::unique_ptr< TMatrixDSym > get_matrix() const {
    // Note that ROOT histogram bin indices are one-based to allow for
    // underflow. The TMatrixDSym element indices, on the other hand,
    // are zero-based.
    int num_reco_bins = cov_matrix_->GetNbinsX();
    auto result = std::make_unique< TMatrixDSym >( num_reco_bins );
    // TODO: consider doing something more efficient than setting each
    // element manually
    for ( int a = 0; a < num_reco_bins; ++a ) {
      for ( int b = 0; b < num_reco_bins; ++b ) {
        result->operator()( a, b ) = cov_matrix_->GetBinContent( a + 1, b + 1 );
      }
    }
    return result;
  }

};

using CovMatrixMap = std::map< std::string, CovMatrix >;

class SystematicsCalculator {

  public:

    using NFT = NtupleFileType;

    SystematicsCalculator( const std::string& input_respmat_file_name,
      const std::string& respmat_tdirectoryfile_name = "" )
    {
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

      // Check whether a set of POT-summed histograms for each universe
      // is already present in the input response matrix file. This is
      // signalled by a TDirectoryFile with a name matching the string
      // TOTAL_SUBFOLDER_NAME.
      TDirectoryFile* total_subdir = nullptr;
      root_tdir->GetObject( TOTAL_SUBFOLDER_NAME.c_str(), total_subdir );

      if ( !total_subdir ) {

        // We couldn't find the pre-computed POT-summed universe histograms,
        // so make them "on the fly" and store them in this object
        this->build_universes( *root_tdir );

        // Create a new TDirectoryFile as a subfolder to hold the POT-summed
        // universe histograms
        total_subdir = new TDirectoryFile( TOTAL_SUBFOLDER_NAME.c_str(),
          "response matrices", "", root_tdir );

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
      // response matrices
      std::string* true_bin_spec = nullptr;
      std::string* reco_bin_spec = nullptr;

      root_tdir->GetObject( TRUE_BIN_SPEC_NAME.c_str(), true_bin_spec );
      root_tdir->GetObject( RECO_BIN_SPEC_NAME.c_str(), reco_bin_spec );

      if ( !true_bin_spec || !reco_bin_spec ) {
        throw std::runtime_error( "Failed to load bin specifications" );
      }

      std::istringstream iss_true( *true_bin_spec );
      TrueBin temp_true_bin;
      while ( iss_true >> temp_true_bin ) {
        true_bins_.push_back( temp_true_bin );
      }

      std::istringstream iss_reco( *reco_bin_spec );
      RecoBin temp_reco_bin;
      while ( iss_reco >> temp_reco_bin ) {
        reco_bins_.push_back( temp_reco_bin );
      }

    }

    void load_universes( TDirectoryFile& total_subdir ) {

      const auto& fpm = FilePropertiesManager::Instance();

      // TODO: reduce code duplication between this function
      // and SystematicsCalculator::build_universes()
      TList* universe_key_list = total_subdir.GetListOfKeys();
      int num_keys = universe_key_list->GetEntries();

      // Loop over the keys in the TDirectoryFile. Build a universe object
      // for each key ending in "_2d" and store it in either the rw_universes_
      // map (for reweightable systematic universes) or the detvar_universes_
      // map (for detector systematic universes)
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

        total_subdir.GetObject( (key + "_true").c_str(), hist_true );
        total_subdir.GetObject( (key + "_reco").c_str(), hist_reco );
        total_subdir.GetObject( (key + "_2d").c_str(), hist_2d );

        if ( !hist_true || !hist_reco || !hist_2d ) {
          throw std::runtime_error( "Failed to retrieve histograms for the "
            + key + " universe" );
        }

        // Reconstruct the Universe object from the retrieved histograms
        auto temp_univ = std::make_unique< Universe >( univ_name, univ_index,
          hist_true, hist_reco, hist_2d );

        // Determine whether the current universe represents a detector
        // variation or a reweightable variation. We'll use this information to
        // decide where it should be stored.
        NFT temp_type = fpm.string_to_ntuple_type( univ_name );
        if ( temp_type != NFT::kUnknown ) {

          bool is_detvar = ntuple_type_is_detVar( temp_type );
          if ( !is_detvar ) throw std::runtime_error( "Universe name "
            + univ_name + " matches a non-detVar file type. Handling of"
            + " this situation is currently unimplemented." );

          if ( detvar_universes_.count(temp_type) ) {
            throw std::runtime_error( "detVar multisims are not currently"
              " supported" );
          }

          // Move the detector variation Universe object into the map
          detvar_universes_[ temp_type ].reset( temp_univ.release() );
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

        TH1D* hist = nullptr;
        total_subdir.GetObject( hist_name.c_str(), hist );
        hist->SetDirectory( nullptr );

        if ( !hist ) {
          throw std::runtime_error( "Missing data histogram for " + data_name );
        }

        if ( data_hists_.count(file_type) ) {
          throw std::runtime_error( "Duplicate data histogram for "
            + data_name );
        }

        data_hists_[ file_type ].reset( hist );
      }

      TParameter< double >* temp_pot = nullptr;
      total_subdir.GetObject( "total_bnb_data_pot", temp_pot );
      if ( !temp_pot ) {
        throw std::runtime_error( "Missing BNB data POT value" );
      }

      total_bnb_data_pot_ = temp_pot->GetVal();

    }

    void build_universes( TDirectoryFile& root_tdir ) {
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
          bool is_reweightable_mc = ntuple_type_is_reweightable_mc( type );
          bool is_mc = ntuple_type_is_mc( type );

          for ( const std::string& file_name : file_set ) {

            std::cout << "PROCESSING response matrices for "
              << file_name << '\n';

            // Get the simulated or measured POT belonging to the current file.
            // This will be used to normalize the relevant histograms
            double file_pot = 0.;
            if ( is_mc ) {
              // MC files have the simulated POT stored alongside the ntuple
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
              TH1D* reco_hist = nullptr;
              subdir->GetObject( "unweighted_0_reco", reco_hist );

              // If we're working with EXT data, scale it to the corresponding
              // number of triggers from the BNB data from the same run
              if ( type == NFT::kExtBNB ) {
                double bnb_trigs = run_to_bnb_trigs_map.at( run );
                double ext_trigs = run_to_ext_trigs_map.at( run );

                reco_hist->Scale( bnb_trigs / ext_trigs );
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
                data_hists_.at( type )->Add( reco_hist );
              }

              // We don't need to do the other MC-ntuple-specific stuff below,
              // so just move on to the next file
              continue;

            } // data ntuple files


            // If we've made it here, then we're working with an MC ntuple
            // file. For these, all three histograms for the "unweighted"
            // universe are always evaluated. Use this to determine the number
            // of true and reco bins easily.
            TH2D* temp_2d_hist = nullptr;
            subdir->GetObject( "unweighted_0_2d", temp_2d_hist );

            // NOTE: the convention of the ResponseMatrixMaker class is to use
            // x as the true axis and y as the reco axis.
            int num_true_bins = temp_2d_hist->GetXaxis()->GetNbins();
            int num_reco_bins = temp_2d_hist->GetYaxis()->GetNbins();

            // Let's handle the detVar samples first.
            if ( is_detVar ) {

              std::string dv_univ_name = fpm.ntuple_type_to_string( type );

              // Right now, we're assuming there's just one detVar ntuple file
              // per universe. If this assumption is violated, our scaling will
              // be screwed up. Prevent this from happening by throwing an
              // exception when a duplicate is encountered.
              if ( detvar_universes_.count(type) ) {
                throw std::runtime_error( "Duplicate detVar ntuple file!" );
              }

              // Otherwise, make a new Universe object to store the
              // (POT-scaled) detVar histograms.
              auto temp_univ_ptr = std::make_unique< Universe >( dv_univ_name,
                0, num_true_bins, num_reco_bins );

              // The detector variation universes are all labeled "unweighted"
              // in the response matrix files. Retrieve the corresponding
              // histograms.
              TH1D* hist_reco = nullptr;
              TH1D* hist_true = nullptr;
              TH2D* hist_2d = nullptr;

              subdir->GetObject( "unweighted_0_reco", hist_reco );
              subdir->GetObject( "unweighted_0_true", hist_true );
              subdir->GetObject( "unweighted_0_2d", hist_2d );

              // Scale all detVar universe histograms from the simulated POT to
              // the *total* BNB data POT for all runs analyzed. Since we only
              // have detVar samples for Run 3b, we assume that they can be
              // applied globally in this step.
              // TODO: revisit this as appropriate
              hist_reco->Scale( total_bnb_data_pot_ / file_pot );
              hist_true->Scale( total_bnb_data_pot_ / file_pot );
              hist_2d->Scale( total_bnb_data_pot_ / file_pot );

              // Add the scaled contents of these histograms to the
              // corresponding histograms in the new Universe object
              temp_univ_ptr->hist_reco_->Add( hist_reco );
              temp_univ_ptr->hist_true_->Add( hist_true );
              temp_univ_ptr->hist_2d_->Add( hist_2d );

              // Adjust the owned histograms to avoid auto-deletion problems
              set_stats_and_dir( *temp_univ_ptr );

              // Move the finished Universe object into the map
              detvar_universes_[ type ].reset( temp_univ_ptr.release() );

            } // detVar sample

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

                  TH1D* h_reco = nullptr;
                  TH1D* h_true = nullptr;
                  TH2D* h_2d = nullptr;

                  subdir->GetObject( (hist_name_prefix + "_reco").c_str(),
                    h_reco );

                  subdir->GetObject( (hist_name_prefix + "_true").c_str(),
                    h_true );

                  subdir->GetObject( (hist_name_prefix + "_2d").c_str(), h_2d );

                  // Scale these histograms to the appropriate BNB data POT for
                  // the current run
                  h_reco->Scale( rw_scale_factor );
                  h_true->Scale( rw_scale_factor );
                  h_2d->Scale( rw_scale_factor );

                  // Add their contributions to the owned histograms for the
                  // current Universe object
                  universe.hist_reco_->Add( h_reco );
                  universe.hist_true_->Add( h_true );
                  universe.hist_2d_->Add( h_2d );

                } // universes indices

              } // universe types

            } // reweightable MC samples

          } // ntuple file

        } // type

      } // run

      std::cout << "TOTAL BNB DATA POT = " << total_bnb_data_pot_ << '\n';
    }

    void save_universes( TDirectoryFile& out_tdf ) {

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

      // Save the detector variation histograms
      for ( const auto& pair : detvar_universes_ ) {
        NFT type = pair.first;
        auto& universe = pair.second;

        universe->hist_reco_->Write();
        universe->hist_true_->Write();
        universe->hist_2d_->Write();
      }

      // Save the reweightable systematic histograms
      for ( const auto& pair : rw_universes_ ) {

        const auto& univ_vec = pair.second;
        for ( const auto& universe : univ_vec ) {
          universe->hist_reco_->Write();
          universe->hist_true_->Write();
          universe->hist_2d_->Write();
        }

      }

      // Save the total BNB data POT for easy retrieval later
      TParameter< double > temp_pot( "total_bnb_data_pot",
        total_bnb_data_pot_ );

      temp_pot.Write();

    }

    const Universe& cv_universe() const {
      return *rw_universes_.at( CV_UNIV_NAME ).front();
    }

    std::unique_ptr< CovMatrixMap > get_covariances() const;

    // Helper functions for comparison to CCNp result
    double effective_efficiency( const Universe& univ, int reco_bin ) const;
    double scaling_factor( int reco_bin, int flux_universe_index = -1 ) const;
    double expected_mc_background( const Universe& univ, int reco_bin,
      bool stat_var = false ) const;
    double forward_folded_xsec( const Universe& univ, int reco_bin,
      int flux_universe_index = -1 ) const;

  //protected:

    CovMatrix make_covariance_matrix( const std::string& hist_name ) const;

    // Central value universe name
    const std::string CV_UNIV_NAME = "weight_TunedCentralValue_UBGenie";

    // Subdirectory name for the TDirectoryFile containing the POT-summed
    // histograms for the various universes across all analysis ntuples
    const std::string TOTAL_SUBFOLDER_NAME = "total";

    // Holds reco-space histograms for data (BNB and EXT) bin counts
    std::map< NFT, std::unique_ptr<TH1D> > data_hists_;

    // Holds universe objects for reweightable systematics
    std::map< std::string, std::vector< std::unique_ptr<Universe> > >
      rw_universes_;

    // Detector systematic universes (and the detVar CV) will be indexed using
    // ntuple file type values. We're currently scaling one set of ntuples
    // (from Run 3b) to the full dataset.
    // TODO: revisit this procedure if new detVar samples become available
    std::map< NFT, std::unique_ptr<Universe> > detvar_universes_;

    // True bin configuration that was used to compute the reponse matrices
    std::vector< TrueBin > true_bins_;

    // Reco bin configuration that was used to compute the reponse matrices
    std::vector< RecoBin > reco_bins_;

    // Total POT exposure for the analyzed BNB data
    double total_bnb_data_pot_ = 0.;
};

CovMatrix SystematicsCalculator::make_covariance_matrix(
  const std::string& hist_name ) const
{
  int num_reco_bins = reco_bins_.size();
  TH2D* hist = new TH2D( hist_name.c_str(),
    "covariance; reco bin; reco bin; covariance", num_reco_bins, 0.,
    num_reco_bins, num_reco_bins, 0., num_reco_bins );
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
  // Get the total number of true and reco bins for later reference
  size_t num_true_bins = sc.true_bins_.size();
  size_t num_reco_bins = sc.reco_bins_.size();

  // Get the expected differential cross-section in each reco bin in the CV
  // universe
  std::vector< double > cv_reco_xsecs( num_reco_bins, 0. );

  for ( size_t rb = 0u; rb < num_reco_bins; ++rb ) {
    cv_reco_xsecs.at( rb ) = sc.forward_folded_xsec( cv_univ, rb + 1 );
  }

  // Loop over universes
  int num_universes = universes.size();
  for ( int u_idx = 0; u_idx < num_universes; ++u_idx ) {

    const auto& univ = universes.at( u_idx );

    // For flux variations, we also want to adjust the integrated flux
    // used to compute the differential cross section. This is signaled
    // by passing a nonnegative index for a flux systematic universe
    // to SystematicsCalculator::forward_folded_xsec()
    int flux_u_idx = -1;
    if ( is_flux_variation ) flux_u_idx = u_idx;

    // Get the expected differential cross section in each reco bin in the
    // current universe.
    std::vector< double > univ_reco_xsecs( num_reco_bins, 0. );

    for ( size_t rb = 0u; rb < num_reco_bins; ++rb ) {
      univ_reco_xsecs.at( rb ) = sc.forward_folded_xsec( *univ,
        rb + 1, flux_u_idx );
    }

    // We have all the needed ingredients to get the contribution of this
    // universe to the covariance matrix. Loop over each pair of reco bins and
    // fill the corresponding covariance matrix elements.
    // TODO: the covariance matrix are symmetric by definition. You can
    // therefore make this more efficient by calculating only the subset of
    // elements that you need.
    for ( size_t a = 0u; a < num_reco_bins; ++a ) {

      double cv_a = cv_reco_xsecs.at( a );
      double univ_a = univ_reco_xsecs.at( a );

      for ( size_t b = 0u; b < num_reco_bins; ++b ) {

        double cv_b = cv_reco_xsecs.at( b );
        double univ_b = univ_reco_xsecs.at( b );

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

  // Look up the location of the configuration file that defines the various
  // covariance matrices that should be calculated
  const auto& fpm = FilePropertiesManager::Instance();
  std::string config_file_name = fpm.analysis_path() + "/systcalc.conf";

  // Read in the definition of each covariance matrix and calculate it. Each
  // definition contains at least a name and a type specifier
  std::ifstream config_file( config_file_name );
  std::string name, type;
  while ( config_file >> name >> type ) {

    CovMatrix temp_cov_mat = this->make_covariance_matrix( name );

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

      int num_reco_bins = reco_bins_.size();
      const auto& cv_univ = this->cv_universe();

      // TODO: optimize so that you can use the Sumw2 array entries instead.
      // This avoids needing to take many square roots only to square them
      // again.
      // NOTE: The reco bin index used here is one-based since ROOT histograms
      // always include an underflow bin
      for ( int rb = 1; rb <= num_reco_bins; ++rb ) {

        double err2 = this->expected_mc_background( cv_univ, rb, true );

        double eff = this->effective_efficiency( cv_univ, rb );
        double scaling = this->scaling_factor( rb );

        if ( eff <= 0. ) err2 = 0.;
        else {
          err2 *= std::pow( scaling / eff, 2 );
        }

        temp_cov_mat.cov_matrix_->SetBinContent( rb, rb, err2 );

      } // reco bins

    } // MCstat type

    else if ( type == "EXTstat" ) {

      const TH1D* ext_hist = data_hists_.at( NFT::kExtBNB ).get();
      int num_reco_bins = ext_hist->GetNbinsX();

      const auto& cv_univ = this->cv_universe();

      // Note the one-based bin numbering convention for TH1D
      for ( int rb = 1; rb <= num_reco_bins; ++rb ) {
        double err2 = ext_hist->GetBinError( rb );
        err2 *= err2;

        // TODO: reduce code duplication with the MCstat type
        double eff = this->effective_efficiency( cv_univ, rb );
        double scaling = this->scaling_factor( rb );

        if ( eff <= 0. ) err2 = 0.;
        else {
          err2 *= std::pow( scaling / eff, 2 );
        }

        temp_cov_mat.cov_matrix_->SetBinContent( rb, rb, err2 );
      } // reco bins

    } // EXTstat type

    else if ( type == "BNBstat" ) {

      // TODO: reduce code duplication with the EXTstat type
      const TH1D* bnb_hist = data_hists_.at( NFT::kOnBNB ).get();
      int num_reco_bins = bnb_hist->GetNbinsX();

      const auto& cv_univ = this->cv_universe();

      // Note the one-based bin numbering convention for TH1D
      for ( int rb = 1; rb <= num_reco_bins; ++rb ) {
        double err2 = bnb_hist->GetBinError( rb );
        err2 *= err2;

        // TODO: reduce code duplication with the MCstat type
        double eff = this->effective_efficiency( cv_univ, rb );
        double scaling = this->scaling_factor( rb );

        if ( eff <= 0. ) err2 = 0.;
        else {
          err2 *= std::pow( scaling / eff, 2 );
        }

        temp_cov_mat.cov_matrix_->SetBinContent( rb, rb, err2 );
      } // reco bins

    } // BNBstat type

    else if ( type == "MCFullCorr" ) {
      // Read in the fractional uncertainty from the configuration file
      double frac_unc = 0.;
      config_file >> frac_unc;

      const double frac2 = std::pow( frac_unc, 2 );
      int num_reco_bins = reco_bins_.size();

      const auto& cv_univ = this->cv_universe();
      for ( int a = 1; a <= num_reco_bins; ++a ) {

        double cv_a = this->forward_folded_xsec( cv_univ, a );

        for ( int b = 1; b <= num_reco_bins; ++b ) {

          double cv_b = this->forward_folded_xsec( cv_univ, b );

          double covariance = cv_a * cv_b * frac2;

          temp_cov_mat.cov_matrix_->SetBinContent( a, b, covariance );

        } // reco bin b

      } // reco bin a

    } // MCFullCorr type

    else if ( type == "DV" ) {
      // Get the detector variation type represented by the current universe
      std::string ntuple_type_str;
      config_file >> ntuple_type_str;
      auto ntuple_type = fpm.string_to_ntuple_type( ntuple_type_str );

      // Check that it's valid. If not, then complain.
      bool is_not_detVar = !ntuple_type_is_detVar( ntuple_type );
      if ( is_not_detVar ) {
        throw std::runtime_error( "Invalid NtupleFileType!" );
      }

      const auto& detVar_cv_u = detvar_universes_.at( NFT::kDetVarMCCV );
      const auto& detVar_alt_u = detvar_universes_.at( ntuple_type );

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


// NOTE: this uses a one-based reco bin index
double SystematicsCalculator::effective_efficiency(
  const Universe& univ, int reco_bin ) const
{
  int num_true_bins = true_bins_.size();
  int num_reco_bins = reco_bins_.size();

  if ( reco_bin > num_reco_bins ) {
    throw std::runtime_error( "Invalid reco bin" );
  }

  double numerator = 0.;
  double denominator = 0.;
  for ( int tb = 1; tb <= num_true_bins; ++tb ) {
    // Note that the true bin definitions have a zero-based index
    auto& tbin = true_bins_.at( tb - 1 );
    // If this isn't a signal true bin, just skip it
    if ( tbin.type_ != kSignalTrueBin ) continue;

    double num_j_gen = univ.hist_true_->GetBinContent( tb );
    double num_ij = univ.hist_2d_->GetBinContent( tb, reco_bin );

    double num_j_sel = 0.;
    for ( int rb = 1; rb <= num_reco_bins; ++rb ) {
      num_j_sel += univ.hist_2d_->GetBinContent( tb, rb );
    }

    double denom_term = num_ij * num_j_gen;
    if ( num_j_sel > 0. ) denom_term /= num_j_sel;
    else denom_term = 0.;

    numerator += num_ij;
    denominator += denom_term;
  }

  double eff = 0.;
  if ( denominator > 0. ) eff = numerator / denominator;

  return eff;
}

// NOTE: this uses a one-based reco bin index
double SystematicsCalculator::expected_mc_background(
  const Universe& univ, int reco_bin, bool stat_var ) const
{
  int num_true_bins = true_bins_.size();
  int num_reco_bins = reco_bins_.size();

  if ( reco_bin > num_reco_bins ) {
    throw std::runtime_error( "Invalid reco bin" );
  }

  double Bi = 0.;
  for ( int tb = 1; tb <= num_true_bins; ++tb ) {
    // Note that the true bin definitions have a zero-based index
    auto& tbin = true_bins_.at( tb - 1 );
    // If this isn't a background true bin, just skip it
    if ( tbin.type_ != kBackgroundTrueBin ) continue;

    // If this flag is set, then calculate the MC statistical variance
    // on the prediction instead of the event count
    if ( stat_var ) {
      double stat_err = univ.hist_2d_->GetBinError( tb, reco_bin );
      Bi += std::pow( stat_err, 2 );
    }
    else {
      Bi += univ.hist_2d_->GetBinContent( tb, reco_bin );
    }
  }

  return Bi;
}

// NOTE: this uses a one-based reco bin index
double SystematicsCalculator::scaling_factor( int reco_bin,
  int flux_universe_index ) const
{
  const auto& ifum = IntegratedFluxUniverseManager::Instance();
  double numu_flux = ifum.cv_numu_integrated_flux(); // numu / POT / cm^2

  if ( flux_universe_index >= 0 ) {
    numu_flux *= ifum.flux_factor( flux_universe_index );
  }

  // Convert to numu / cm^2 by multiplying by the total analyzed beam data POT
  numu_flux *= total_bnb_data_pot_;

  const auto& cv_univ = this->cv_universe();
  // TODO: fix this. The universe histograms are expressed in terms of
  // reco bin number, so the widths are always trivially one. You need
  // to make this code aware of the physics units on the x-axis.
  double reco_bin_width = cv_univ.hist_reco_->GetBinWidth( reco_bin );

  constexpr double NUM_TARGETS_ACTIVE_VOL = 1.0068e30;

  double factor = 1. / ( numu_flux * NUM_TARGETS_ACTIVE_VOL * reco_bin_width );
  return factor;
}

// NOTE: this uses a one-based reco bin index
double SystematicsCalculator::forward_folded_xsec( const Universe& univ,
  int reco_bin, int flux_universe_index ) const
{
  const TH1D* bnb_hist = data_hists_.at( NFT::kOnBNB ).get();
  double data_counts = bnb_hist->GetBinContent( reco_bin );

  const TH1D* ext_hist = data_hists_.at( NFT::kExtBNB ).get();
  double ext_counts = ext_hist->GetBinContent( reco_bin );

  double mc_bkgd_counts = this->expected_mc_background( univ, reco_bin );

  double eff = this->effective_efficiency( univ, reco_bin );
  double scaling = this->scaling_factor( reco_bin, flux_universe_index );

  double xsec = ( data_counts - ext_counts - mc_bkgd_counts ) * scaling / eff;
  return xsec;
}
