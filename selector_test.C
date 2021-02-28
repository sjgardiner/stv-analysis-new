#include "FilePropertiesManager.hh"
#include "MySelector.hh"
#include "ResponseMatrixMaker.hh"
#include "TreeUtils.hh"

// Use this as the bin index for output TTree entries in which one of the bins
// is undefined
constexpr int UNDEFINED_BIN = -1;

void sum_event_weights( TTree& tree, MySelector& ms, TEntryList* el,
  NtupleFileType nf_type )
{
  // Actual summation of event weights only needs to be performed for CV MC
  // files. This function is a no-op otherwise.
  if ( nf_type != NtupleFileType::kNumuMC
    && nf_type != NtupleFileType::kIntrinsicNueMC
    && nf_type != NtupleFileType::kDirtMC ) return;

  //// Perform the summation
  tree.SetEntryList( el );
  tree.Process( &ms );
  tree.SetEntryList( nullptr );
}

void make_response_matrix( const std::string& input_file_name ) {

  // Check the type of file used as input. This will influence how the output
  // TTree is prepared.
  const auto& fpm = FilePropertiesManager::Instance();
  NtupleFileType nf_type = fpm.get_ntuple_file_type( input_file_name );

  ResponseMatrixMaker resp_mat( "myconfig.txt" );

  resp_mat.add_input_file( input_file_name.c_str() );

  // Get access to the ResponseMatrixMaker object's owned input TChain
  auto& input_chain = resp_mat.input_chain();

  //std::vector< std::string > weight_names;
  //weight_names.push_back( "weight_TunedCentralValue_UBGenie" );
  ////weight_names.push_back( "weight_All_UBGenie" );
  ////weight_names.push_back( "weight_AxFFCCQEshape_UBGenie" );

  // Find relevant entries in the input TChain for each of the reco space bins
  resp_mat.build_reco_entry_lists();
  size_t num_reco_bins = resp_mat.reco_bins().size();
  // Default to zero true bins. This will be correct as needed below
  size_t num_true_bins = 0u;

  // If we're working with MC events, then also find the relevant entries in
  // the input TChain for each of the true space bins
  if ( nf_type != NtupleFileType::kOnBNB
    && nf_type != NtupleFileType::kExtBNB )
  {
    resp_mat.build_true_entry_lists();
    num_true_bins = resp_mat.true_bins().size();
  }

  //MySelector selector( &weight_names );
  MySelector selector;

  // Configure the output TTree
  // TODO: set output file name in a smarter way
  // Get the basename of the input file using the trick described here:
  // https://stackoverflow.com/a/24386991
  std::string input_basename = input_file_name.substr(
    input_file_name.find_last_of('/') + 1 );

  std::string out_file_name = "/uboone/data/users/gardiner/ntuples-stv"
    "/resp/respmat-" + input_basename;

  TFile out_file( out_file_name.c_str(), "recreate" );
  TTree* out_tree = new TTree( "resp_mat_tree", "response matrix elements" );
  int tb, rb; // bin indices in true and reco space
  int count; // unweighted event count
  out_tree->Branch( "tbin", &tb, "tbin/I" );
  out_tree->Branch( "rbin", &rb, "rbin/I" );
  out_tree->Branch( "count", &count, "count/I" );
  bool need_to_create_other_branches = true;

  for ( size_t true_bin = 0u; true_bin < num_true_bins; ++true_bin ) {

    std::cout << "  true bin = " << true_bin << '\n';

    const auto& tel = resp_mat.true_entry_lists().at( true_bin );

    // Get the sum of the event weights in the current true bin
    sum_event_weights( input_chain, selector, tel.get(), nf_type );
    tb = true_bin;
    rb = UNDEFINED_BIN;
    count = tel->GetN();
    if ( need_to_create_other_branches ) {
      for ( auto& pair : selector.sums_of_weights() ) {
        std::string br_name = "sum_" + pair.first;
        set_object_output_branch_address( *out_tree,
          br_name, pair.second, true );
      }
      for ( auto& pair : selector.sums_of_weights2() ) {
        std::string br_name = "sum2_" + pair.first;
        set_object_output_branch_address( *out_tree,
          br_name, pair.second, true );
      }
      // Don't repeat branch creation now that we've already done it once
      need_to_create_other_branches = false;
    }
    out_tree->Fill();

    for ( size_t reco_bin = 0u; reco_bin < num_reco_bins; ++reco_bin ) {

      std::cout << "    reco bin = " << reco_bin << '\n';

      const auto& rel = resp_mat.reco_entry_lists().at( reco_bin );

      // Otherwise proceed normally with the summation of the event weights
      const auto& tel = resp_mat.true_entry_lists().at( true_bin );

      // Compute the intersection of the entry lists for the current true
      // and reco bins
      auto tr_el = resp_mat.get_intersection_entry_list( true_bin, reco_bin );
      //auto tr_el = tel * rel;

      std::cout << true_bin << ' ' << reco_bin << "    trel = " << tr_el->GetN() << '\n';

      // Get the sum of the event weights in the current 2D bin
      // in (true, reco) space
      sum_event_weights( input_chain, selector, tr_el.get(), nf_type );

      // Save the results to the output TTree
      tb = true_bin;
      rb = reco_bin;
      count = tr_el->GetN();
      out_tree->Fill();

    } // reco bin
  } // true bin

  for ( size_t reco_bin = 0u; reco_bin < num_reco_bins; ++reco_bin ) {

    std::cout << "  reco bin = " << reco_bin << '\n';

    const auto& rel = resp_mat.reco_entry_lists().at( reco_bin );

    // Get the sum of the event weights in the current reco bin.
    sum_event_weights( input_chain, selector, rel.get(), nf_type );

    // Save the results to the output TTree
    tb = UNDEFINED_BIN;
    rb = reco_bin;
    count = rel->GetN();
    out_tree->Fill();
  } // reco bin

  out_tree->Write();
}

void selector_test() {
  //ROOT::EnableImplicitMT();
  //make_response_matrix( "/uboone/data/users/gardiner/ntuples-stv/stv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run2_reco2_D1D2_reco2.root" );
  //make_response_matrix( "/uboone/data/users/gardiner/ntuples-stv/stv-run1_neutrinoselection_filt_numu_ALL.root" );

  std::ifstream list_file( "foo2" );
  std::string file_name;
  while ( std::getline(list_file, file_name, '\n') ) {
    make_response_matrix( file_name );
  }

}

int main() {
  selector_test();
  return 0;
}
