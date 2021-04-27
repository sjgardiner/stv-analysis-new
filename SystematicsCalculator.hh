#pragma once

// ROOT includes
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

// STV analysis includes
#include "FilePropertiesManager.hh"
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


class SystematicsCalculator {

  public:

    using NFT = NtupleFileType;

    SystematicsCalculator( const std::string& input_respmat_file_name,
      const std::string& respmat_tdirectoryfile_name = "" )
    {
      TFile in_tfile( input_respmat_file_name.c_str(), "read" );

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
      double total_bnb_data_pot = 0.;
      for ( const auto& pair : run_to_bnb_pot_map ) {
        total_bnb_data_pot += pair.second;
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
            root_tdir->GetObject( subdir_name.c_str(), subdir );
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
              hist_reco->Scale( total_bnb_data_pot / file_pot );
              hist_true->Scale( total_bnb_data_pot / file_pot );
              hist_2d->Scale( total_bnb_data_pot / file_pot );

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
                  Universe temp_univ( univ_name, univ_index,
                    num_true_bins, num_reco_bins );

                  set_stats_and_dir( temp_univ );

                  // If we do not already have a map entry for this kind of
                  // universe, then create one
                  if ( !rw_universes_.count(univ_name) ) {
                    rw_universes_[ univ_name ] = std::vector< Universe >();
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
                  auto& universe = univ_vec.at( u_idx );

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

      std::cout << "TOTAL BNB DATA POT = " << total_bnb_data_pot << '\n';
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
          universe.hist_reco_->Write();
          universe.hist_true_->Write();
          universe.hist_2d_->Write();
        }

      }

    }

  protected:

    // Holds reco-space histograms for data (BNB and EXT) bin counts
    std::map< NFT, std::unique_ptr<TH1D> > data_hists_;

    // Holds universe objects for reweightable systematics
    std::map< std::string, std::vector<Universe> > rw_universes_;

    // Detector systematic universes (and the detVar CV) will be indexed using
    // ntuple file type values. We're currently scaling one set of ntuples
    // (from Run 3b) to the full dataset.
    // TODO: revisit this procedure if new detVar samples become available
    std::map< NFT, std::unique_ptr<Universe> > detvar_universes_;

};
