// Executable for generating response matrix files for later analysis. It
// has been adapted from a similar ROOT macro.

// ROOT includes
#include "TROOT.h"

// STV analysis includes
#include "FilePropertiesManager.hh"
#include "ResponseMatrixMaker.hh"


int main( int argc, char* argv[] ) {

  if ( argc != 3 ) {
    std::cout << "Usage: respmat CONFIG_FILE OUTPUT_ROOT_FILE\n";
    return 1;
  }

  std::string config_file_name( argv[1] );
  std::string output_file_name( argv[2] );

  // Build a vector of input ntuple file names using the singleton
  // FilePropertiesManager. For now, include every ntuple file that we
  // plan to use in the STV analysis.
  // TODO: consider adding command-line options to allow the user to specify a
  // subset of the full set of files managed by the FilePropertiesManager.
  std::vector< std::string > input_file_names;

  const auto& fpm = FilePropertiesManager::Instance();
  for ( const auto& run_and_type_pair : fpm.ntuple_file_map() ) {

    const auto& type_map = run_and_type_pair.second;

    for ( const auto& type_and_files_pair : type_map ) {
      const auto& file_set = type_and_files_pair.second;

      for ( const std::string& file_name : file_set ) {
        input_file_names.push_back( file_name );
      }
    }
  }

  ROOT::EnableImplicitMT();

  for ( const auto& input_file_name : input_file_names ) {

    std::cout << "Calculating response matrices for ntuple input file "
      << input_file_name << '\n';
    std::cout << "Loading ResponseMatrixMaker configuration from "
      << config_file_name << '\n';

    ResponseMatrixMaker resp_mat( config_file_name );

    resp_mat.add_input_file( input_file_name.c_str() );

    resp_mat.build_response_matrices();

    resp_mat.save_histograms( output_file_name, input_file_name );

  } // loop over input files

  return 0;
}
