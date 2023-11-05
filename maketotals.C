// Executable that adds the totals directory to a file filled with histograms

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
  //fpm.load_file_properties(list_file_name);
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

  // Just use the first file in the config to get the name of the total dir
  const auto& input_file_name = input_files.front();
  UniverseMaker univ_maker( univmake_config_file_name );
  univ_maker.add_input_file( input_file_name.c_str() );  
  MCC9SystematicsCalculator unfolder( output_file_name, "", univ_maker.dir_name() );
  
  return 0;
}
