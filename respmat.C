// Executable for generating response matrix files for later analysis. It
// has been adapted from a similar ROOT macro.

// ROOT includes
#include "TROOT.h"

// STV analysis includes
#include "FilePropertiesManager.hh"
#include "MCC9Unfolder.hh"
#include "ResponseMatrixMaker.hh"

int main( int argc, char* argv[] ) {

  if ( argc != 4 ) {
    std::cout << "Usage: respmat FILE_PROPERTIES_CONFIG_FILE"
      << " RESPMAT_CONFIG_FILE OUTPUT_ROOT_FILE\n";
    return 1;
  }

  std::string fp_config_file_name( argv[1] );
  std::string respmat_config_file_name( argv[2] );
  std::string output_file_name( argv[3] );

  // Build a vector of input ntuple file name/type pairs using the singleton
  // FilePropertiesManager. Use the file properties configuration file
  // specified by the user on the command line.
  std::vector< std::pair<std::string, NtupleFileType> > input_files;

  std::cout << "Loading FilePropertiesManager configuration from "
    << fp_config_file_name << '\n';

  auto& fpm = FilePropertiesManager::Instance();
  fpm.load_file_properties( fp_config_file_name );

  for ( const auto& run_and_type_pair : fpm.ntuple_file_map() ) {

    const auto& type_map = run_and_type_pair.second;

    for ( const auto& type_and_files_pair : type_map ) {
      const auto& type = type_and_files_pair.first;
      const auto& file_set = type_and_files_pair.second;

      for ( const std::string& file_name : file_set ) {
        input_files.emplace_back( file_name, type );
      }
    }
  }

  ROOT::EnableImplicitMT();

  // Store the name of the root TDirectoryFile created by the
  // ResponseMatrixMaker objects below. We will use it to ensure that
  // the MCC9Unfolder object used to calculate the total event counts
  // will always be working with the correct sets of universes.
  std::string tdirfile_name;
  bool set_tdirfile_name = false;

  for ( const auto& pair : input_files ) {

    const auto& input_file_name = pair.first;
    const auto& type = pair.second;

    std::cout << "Calculating response matrices for ntuple input file "
      << input_file_name << '\n';
    std::cout << "Loading ResponseMatrixMaker configuration from "
      << respmat_config_file_name << '\n';

    ResponseMatrixMaker resp_mat( respmat_config_file_name );

    resp_mat.add_input_file( input_file_name.c_str() );

    if ( ntuple_type_is_detVar(type) || !ntuple_type_is_mc(type) ) {
      // Ignore all event weights in the detVar post-processed ntuples
      // TODO: revisit this if you get new samples in which the CV correction
      // weights are calculated correctly
      resp_mat.build_response_matrices( { "FAKE_BRANCH_NAME" } );
    }
    else {
      resp_mat.build_response_matrices();
    }

    resp_mat.save_histograms( output_file_name, input_file_name );

    // The root TDirectoryFile name is the same across all iterations of this
    // loop, so just set it once on the first iteration
    if ( !set_tdirfile_name ) {
      tdirfile_name = resp_mat.dir_name();
      set_tdirfile_name = true;
    }

  } // loop over input files

  // Use a temporary MCC9Unfolder object to automatically calculate the total
  // event counts in each universe across all input files. Since the
  // get_covariances() member function is never called, the specific
  // systematics configuration file used doesn't matter. The empty string
  // passed as the second argument to the constructor just instructs the
  // MCC9Unfolder class to use the default systematics configuration file.
  MCC9Unfolder unfolder( output_file_name, "", tdirfile_name );

  return 0;
}
