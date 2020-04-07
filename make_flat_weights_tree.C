// Run with genie -l
#include <iostream>
#include <map>
#include <string>

void make_flat_weights_tree( const std::string& input_file_name,
  const std::string& output_file_name )
{

  genie::NtpMCEventRecord* gmcrec = nullptr;
  std::map< std::string, std::vector<double> >* weights_map = nullptr;

  TChain gtree_ch( "NuCCanalyzer/gtree" );
  gtree_ch.Add( input_file_name.c_str() );
  gtree_ch.SetBranchAddress( "weights_map", &weights_map );

  TFile out_file( output_file_name.c_str(), "recreate" );
  TTree* out_tree = new TTree( "flat_gtree", "Flat tree for GENIE + EventWeight" );

  // First make a branch for the GENIE events
  out_tree->Branch( "gmcrec", "genie::NtpMCEventRecord", &gmcrec );

  // Get the first map entry and scan through it. For each kind of weight in
  // the map, make a suitable output TTree branch.
  gtree_ch.GetEntry( 0 );
  for ( auto& pair : *weights_map ) {
    std::string label = pair.first;
    auto& weights = pair.second;

    // Here we assume that every entry has the same number of universes for a given
    // weight label. This allows us to use fixed-size arrays to create the branches.
    size_t num_weights = weights.size();
    std::string leaflist = label + '[' + std::to_string(num_weights) + "]/D";

    out_tree->Branch( label.c_str(), weights.data(), leaflist.c_str() );
  }

  // Okay, we're ready. Dump the weights into the new flat TTree structure
  for ( long e = 0; e < gtree_ch.GetEntries(); ++e ) {

    if ( e % 1000 == 0 ) std::cout << "Processing event #" << e << '\n';

    gtree_ch.GetEntry( e );

    // Update the branch addresses for this entry. The old ones are not guaranteed
    // to be valid after reading in a new std::map
    for ( auto& pair : *weights_map ) {
      std::string label = pair.first;
      out_tree->SetBranchAddress( label.c_str(), pair.second.data() );
    }

    out_tree->Fill();

    // Prevent memory leaks by deleting the owned genie::EventRecord object.
    // We need to do this for dumb ROOT reasons.
    delete gmcrec->event;
  }
}
