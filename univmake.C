// Executable for generating systematic universe files for later analysis. It
// has been adapted from a similar ROOT macro.

// Standard library includes
#include <stdexcept>

// ROOT includes
#include "TBranch.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"

// STV analysis includes
#include "FilePropertiesManager.hh"
#include "MCC9SystematicsCalculator.hh"
#include "UniverseMaker.hh"

// Helper function that checks whether a given ROOT file represents an ntuple
// from a reweightable MC sample. This is done by checking for the presence of
// a branch whose name matches the TUNE_WEIGHT_NAME string defined in
// UniverseMaker.hh. Central-value GENIE MC samples are expected to have
// this branch. Real data, detector variation systematics samples and MC
// samples prepared using alternative generators (e.g., NuWro) are not expected
// to have this branch.
bool is_reweightable_mc_ntuple( const std::string& input_file_name ) {
  TFile temp_file( input_file_name.c_str(), "read" );
  TTree* stv_tree = nullptr;
  temp_file.GetObject( "stv_tree", stv_tree );
  if ( !stv_tree ) throw std::runtime_error( "Missing TTree \"stv_tree\" in"
    " the input ROOT file " + input_file_name );

  TBranch* cv_weight_br = stv_tree->GetBranch( TUNE_WEIGHT_NAME.c_str() );
  bool has_cv_weights = ( cv_weight_br != nullptr );
  return has_cv_weights;
}

int main( int argc, char* argv[] ) {

  if ( argc != 4 && argc != 5 ) {
    std::cout << "Usage: univmake LIST_FILE"
      << " UNIVMAKE_CONFIG_FILE OUTPUT_ROOT_FILE"
      << " [FILE_PROPERTIES_CONFIG_FILE]\n";
    return 1;
  }

  std::string list_file_name( argv[1] );
  std::string univmake_config_file_name( argv[2] );
  std::string output_file_name( argv[3] );

  // If the user specified an (optional) non-default configuration file for the
  // FilePropertiesManager on the command line, then load it here. Note that the
  // only place where the FilePropertiesManager configuration is relevant is in
  // the use of MCC9SystematicsCalculator to compute total event count
  // histograms (see below).
  auto& fpm = FilePropertiesManager::Instance();
  if ( argc == 5 ) {
    fpm.load_file_properties( argv[4] );
  }
  // Regardless of whether the default was used or not, retrieve the
  // name of the FilePropertiesManager configuration file that was
  // actually used
  std::string fp_config_file_name = fpm.config_file_name();

  std::cout << "Loaded FilePropertiesManager configuration from "
    << fp_config_file_name << '\n';

  // Read in the complete list of input ntuple files that should be processed
  std::ifstream in_file( list_file_name );
  std::vector< std::string > input_files;
  std::string temp_line;
  while ( std::getline(in_file, temp_line) ) {
    // Ignore lines that begin with the '#' character (this allows for
    // comments in the normalization table file
    if ( temp_line.front() == '#' ) continue;

    // Read in the ntuple file name from the beginning of the current line of
    // the list file. Any trailing line contents separated from the name by
    // whitespace will be ignored.
    std::string file_name;
    std::istringstream temp_ss( temp_line );
    temp_ss >> file_name;

    input_files.push_back( file_name );
  }

  std::cout << "Processing systematic universes for a total of "
    << input_files.size() << " input ntuple files\n";

  ROOT::EnableImplicitMT();

  // Store the name of the root TDirectoryFile created by the UniverseMaker
  // objects below. We will use it to ensure that the MCC9SystematicsCalculator
  // object used to calculate the total event counts will always be working with
  // the correct sets of universes.
  std::string tdirfile_name;
  bool set_tdirfile_name = false;

  for ( const auto& input_file_name : input_files ) {

    std::cout << "Calculating systematic universes for ntuple input file "
      << input_file_name << '\n';
    std::cout << "Loading UniverseMaker configuration from "
      << univmake_config_file_name << '\n';

    UniverseMaker univ_maker( univmake_config_file_name );

    univ_maker.add_input_file( input_file_name.c_str() );

    bool has_event_weights = is_reweightable_mc_ntuple( input_file_name );

    if ( has_event_weights ) {
      // If the check above was successful, then run all of the histogram
      // calculations in the usual way
      univ_maker.build_universes();
    }
    else {
      // Passing in the fake list of explicit branch names below instructs
      // the UniverseMaker class to ignore all event weights while
      // processing the current ntuple
      univ_maker.build_universes( { "FAKE_BRANCH_NAME" } );
    }

    univ_maker.save_histograms( output_file_name, input_file_name );

    // The root TDirectoryFile name is the same across all iterations of this
    // loop, so just set it once on the first iteration
    if ( !set_tdirfile_name ) {
      tdirfile_name = univ_maker.dir_name();
      set_tdirfile_name = true;
    }

  } // loop over input files

  // Use a temporary MCC9SystematicsCalculator object to automatically calculate the total
  // event counts in each universe across all input files. Since the
  // get_covariances() member function is never called, the specific
  // systematics configuration file used doesn't matter. The empty string
  // passed as the second argument to the constructor just instructs the
  // MCC9SystematicsCalculator class to use the default systematics configuration file.
  MCC9SystematicsCalculator unfolder( output_file_name, "", tdirfile_name );

  return 0;
}
