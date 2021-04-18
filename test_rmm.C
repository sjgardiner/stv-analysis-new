// ROOT includes
#include "TROOT.h"

// STV analysis includes
#include "FilePropertiesManager.hh"
#include "ResponseMatrixMaker.hh"

void make_response_matrix( const std::vector<std::string>& input_file_names,
  const std::string& output_file_name )
{
  for ( const auto& input_file_name : input_file_names ) {

    ResponseMatrixMaker resp_mat( "myconfig.txt" );

    resp_mat.add_input_file( input_file_name.c_str() );

    resp_mat.build_response_matrices();

    //// Get the basename of the input file using the trick described here:
    //// https://stackoverflow.com/a/24386991
    //std::string input_basename = input_file_name.substr(
    //  input_file_name.find_last_of('/') + 1 );

    resp_mat.save_histograms( output_file_name, input_file_name );
  }
}

void test_rmm() {
  ROOT::EnableImplicitMT();

  make_response_matrix( { "/uboone/data/users/gardiner/new-ntuples-stv/"
    "stv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run2"
    "_reco2_D1D2_reco2.root" },

    "test_respmat_out2.root"
  );

  //make_response_matrix( "/uboone/data/users/gardiner/ntuples-stv/stv-run1_neutrinoselection_filt_numu_ALL.root" );

  //std::ifstream list_file( "foo2" );
  //std::string file_name;
  //while ( std::getline(list_file, file_name, '\n') ) {
  //  make_response_matrix( file_name );
  //}
}

int main() {
  test_rmm();
  return 0;
}
