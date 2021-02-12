#include "FilePropertiesManager.hh"
#include "MySelector.hh"
#include "ResponseMatrixMaker.hh"
#include "TreeUtils.hh"

void make_response_matrix( const std::string& input_file_name ) {

  ResponseMatrixMaker resp_mat( "myconfig.txt" );

  resp_mat.add_input_file( input_file_name.c_str() );

  // Get access to the ResponseMatrixMaker object's owned input TChain
  auto& input_chain = resp_mat.input_chain();

  std::vector< std::string > weight_names;
  weight_names.push_back( "weight_All_UBGenie" );
  //weight_names.push_back( "weight_AxFFCCQEshape_UBGenie" );

  // Find relevant entries in the input TChain for each of the reco space bins
  resp_mat.build_true_entry_lists();
  resp_mat.build_reco_entry_lists();

  size_t num_true_bins = resp_mat.true_bins().size();
  size_t num_reco_bins = resp_mat.reco_bins().size();

  for ( const auto& wn : weight_names ) {

    MySelector selector( wn );
    std::cout << wn << '\n';

    for ( size_t true_bin = 0u; true_bin < num_true_bins; ++true_bin ) {

      std::cout << "  true bin = " << true_bin << '\n';

      const auto& tel = resp_mat.true_entry_lists().at( true_bin );

      // Get the sum of the event weights in the current true bin
      input_chain.SetEntryList( tel.get() );
      input_chain.Process( &selector );
      input_chain.SetEntryList( nullptr );

      for ( size_t reco_bin = 0u; reco_bin < num_reco_bins; ++reco_bin ) {

        std::cout << "    reco bin = " << reco_bin << '\n';

        const auto& rel = resp_mat.reco_entry_lists().at( reco_bin );

        // Otherwise proceed normally with the summation of the event weights
        const auto& tel = resp_mat.true_entry_lists().at( true_bin );

        // Compute the intersection of the entry lists for the current true
        // and reco bins
        auto tr_el = tel * rel;

        // Get the sum of the event weights in the current 2D bin
        // in (true, reco) space
        input_chain.SetEntryList( tr_el.get() );
        input_chain.Process( &selector );
        input_chain.SetEntryList( nullptr );
      } // reco bin
    } // true bin

    for ( size_t reco_bin = 0u; reco_bin < num_reco_bins; ++reco_bin ) {

      std::cout << "  reco bin = " << reco_bin << '\n';

      const auto& rel = resp_mat.reco_entry_lists().at( reco_bin );

      // Get the sum of the event weights in the current reco bin
      input_chain.SetEntryList( rel.get() );
      input_chain.Process( &selector );
      input_chain.SetEntryList( nullptr );
    } // reco bin
  } // weight ID
}

void selector_test() {
  //ROOT::EnableImplicitMT();
  make_response_matrix( "/uboone/data/users/gardiner/ntuples-stv/stv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run2_reco2_D1D2_reco2.root" );
  //make_response_matrix( "/uboone/data/users/gardiner/ntuples-stv/stv-run1_neutrinoselection_filt_numu_ALL.root" );

  //std::ifstream list_file( "foo2" );
  //std::string file_name;
  //while ( std::getline(list_file, file_name, '\n') ) {
  //  make_response_matrix( file_name );
  //}

}
