#include "WeightHandler.hh"

const std::string INPUT_DIR_NAME = "/uboone/data/users/gardiner/ntuples-stv";
const std::string OUTPUT_DIR_NAME = INPUT_DIR_NAME + "/weight_dumps";

void dump_them( const std::string& input_file_name ) {

  TChain events_ch( "stv_tree" );

  std::string full_input_file_name = INPUT_DIR_NAME + '/' + input_file_name;

  events_ch.Add( full_input_file_name.c_str() );

  // Create temporary storage for the systematic variation event weights
  WeightHandler wh;

  // Sync the branch addresses of the input TChain with the WeightHandler
  // object
  wh.set_branch_addresses( events_ch );

  std::string full_output_file_name = OUTPUT_DIR_NAME + '/'
    + "weights-" + input_file_name;

  TFile out_file( full_output_file_name.c_str(), "recreate" );

  TTree* out_tree = new TTree( "weight_tree", "weight_tree" );

  long events_entry = 0;
  bool need_to_create_branches = true;

  while ( true ) {

    if ( events_entry % 1000 == 0 ) {
      std::cout << "  Processing event #" << events_entry << '\n';
    }

    // TChain::LoadTree() returns the entry number that should be used with
    // the current TTree object, which (together with the TBranch objects
    // that it owns) doesn't know about the other TTrees in the TChain.
    // If the return value is negative, there was an I/O error, or we've
    // attempted to read past the end of the TChain.
    int local_entry = events_ch.LoadTree( events_entry );

    // If we've reached the end of the TChain (or encountered an I/O error),
    // then terminate the event loop
    if ( local_entry < 0 ) break;

    events_ch.GetEntry( events_entry );

    // Make separate branches for individual systematic variation
    // weights in the map. This is hacky and results in tons of branches,
    // but it helps to optimize response matrix evaluation downstream.
    // NOTE: I cheat here by relying on the branch addresses staying the
    // same throughout the entire event loop. Setting the branch addresses
    // during every iteration is very slow. This is fragile but appears to
    // work thanks to consistent weight vector sizes.
    if ( need_to_create_branches ) {
      for ( auto& pair : wh.weight_map() ) {

        const auto& weight_name = pair.first;
        auto& weight_vec = pair.second;

        size_t num_weights = weight_vec->size();
        for ( size_t w = 0u; w < num_weights; ++w ) {

          // Add the universe index to the branch name
          std::string weight_branch_name = weight_name + std::to_string( w );

          auto* weight_ptr = &weight_vec->at( w );

          std::string leaf_spec = weight_branch_name + "/D";

          // Set the branch address for this vector of weights
          set_output_branch_address( *out_tree, weight_branch_name,
            weight_ptr, need_to_create_branches, leaf_spec );
        }
      }

      need_to_create_branches = false;
    }

    out_tree->Fill();

    ++events_entry;
  }

  out_tree->Write();

}

void dump_weights() {
  std::ifstream in_file( "weight_dump_list.txt" );
  std::string input_file_name;
  while ( std::getline(in_file, input_file_name, '\n') ) {
    std::cout << "Processing event weights for " << input_file_name << '\n';
    dump_them( input_file_name );
  }
}
