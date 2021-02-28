#include "FilePropertiesManager.hh"
#include "ResponseMatrixMaker.hh"

const std::string UNWEIGHTED_ID = "unweighted";
const std::string MAX_WEIGHT = "30.";

// Utility function that sums event weights (as computed according to a
// TTree::Draw cut expression) for TChain entries included in a TEntryList.
// This function saves both the sum itself and its squared statistical
// uncertainty (given by the sum of the squares of the weights).
inline void sum_event_weights( TChain& ch, const std::string& cut_expr,
  TEntryList* el, double& sum, double& sum_err2 )
{
  // Make a dummy histogram to hold the summed event weights. Using it
  // with TTree::Draw is much faster than summing manually.
  auto tmp_h = std::make_unique< TH1D >( "tmp_h", "tmp_h", 1, 0., 1. );

  // Tell the histogram to also store the sum of the squares of the event
  // weights
  tmp_h->Sumw2();

  // Configure the TChain to use the given entry list
  ch.SetEntryList( el );

  // Sum the event weights. Fill the center of the only histogram bin (0.5)
  // each time.
  ch.Draw( "0.5 >> tmp_h", cut_expr.c_str(), "e goff" );

  // Remove the entry list now that we're done with it
  ch.SetEntryList( nullptr );

  // Remember that, for a TH1, bin zero is the underflow bin. Retrieve the
  // summation results before returning.
  sum = tmp_h->GetBinContent( 1 );
  sum_err2 = tmp_h->GetSumw2()->At( 1 );
}

// Utility function used to check endings of (trimmed) weight labels based on
// branch names in the weights TTree
bool string_has_end( const std::string& str, const std::string& end ) {
  if ( str.length() >= end.length() ) {
    int comp = str.compare( str.length() - end.length(),
      end.length(), end );
    bool test_result = ( comp == 0 );
    return test_result;
  }
  return false;
}

// Utility function used to remove trailing digits (used as universe labels)
// from branch names in the weights TTree. Based on nice examples from
// https://stackoverflow.com/a/217605.
inline void trim_digits_from_end_in_place( std::string& s ) {
  s.erase( std::find_if(s.rbegin(), s.rend(),
    [](unsigned char ch) { return !std::isdigit(ch); }).base(), s.end() );
}

// TODO: include the rootino_fix weight as a correction to the central value
std::string get_weight_cut_expr( const std::string& wgt_name ) {
  std::string result;
  if ( string_has_end(wgt_name, "UBGenie") ) {
    result = "weight_splines_general_Spline0 * ";
  }
  else if ( wgt_name == "weight_flux_all"
    || wgt_name == "weight_reint_all"
    || wgt_name == "weight_xsr_scc_Fa3_SCC"
    || wgt_name == "weight_xsr_scc_Fv3_SCC" )
  {
    result = "weight_splines_general_Spline0 *"
      " weight_TunedCentralValue_UBGenie0 * ";
  }
  else if ( wgt_name == "weight_splines_general_Spline" ) {
    // No extra weight factors needed
    result = "";
  }
  else if ( wgt_name == UNWEIGHTED_ID ) {
    result = "1.0";
    return result;
  }
  else throw std::runtime_error( "Unrecognized weight name in "
    " get_weight_cut_expr()" );

  result += wgt_name;
  return result;
}

void make_response_matrix( const std::string& input_file_name ) {

  ResponseMatrixMaker resp_mat( "myconfig.txt" );

  resp_mat.add_input_file( input_file_name.c_str() );

  // Get access to the ResponseMatrixMaker object's owned input TChain
  auto& input_chain = resp_mat.input_chain();

  // Build a map of weight names (ignoring universe indices) by looping over
  // all of the branches and trimming off trailing digits. We'll use the map
  // values to keep track of universe counts for each unique weight name. For
  // data and detVar MC samples, this map will be empty except for a special
  // "unweighted" universe that we add below.
  std::map< std::string, unsigned > weight_names;

  // Check if we need to associate dumped event weights with the input TChain.
  // This is a requirement only for non-detVar MC samples.
  const auto& fpm = FilePropertiesManager::Instance();
  NtupleFileType nf_type = fpm.get_ntuple_file_type( input_file_name );

  // Set up a helper TChain to manage reading weights from the "weight dump"
  // input file (if needed)
  TChain weights_ch( "weight_tree" );

  if ( nf_type == NtupleFileType::kNumuMC
    || nf_type == NtupleFileType::kIntrinsicNueMC
    || nf_type == NtupleFileType::kDirtMC )
  {
    // Retrieve the name of the weight dump file from the FilePropertiesManager
    const auto& wgt_file_map = fpm.ntuple_to_weight_file_map();
    std::string wgt_file_name = wgt_file_map.at( input_file_name );

    // Add the weight file to the helper TChain
    weights_ch.Add( wgt_file_name.c_str() );

    auto* lob = weights_ch.GetListOfBranches();
    for ( int b = 0; b < lob->GetEntries(); ++b ) {

      // Get the name of the current branch
      auto* branch = dynamic_cast< TBranch* >( lob->At(b) );
      std::string br_name = branch->GetName();

      // Do the trimming
      trim_digits_from_end_in_place( br_name );

      // Check to see if the current weight name is already in the map
      auto iter = weight_names.find( br_name );
      // If it is, just increment the matching counter by one
      if ( iter != weight_names.end() ) ++iter->second;
      // If it isn't, add it and a new universe counter starting at 1
      else {
        weight_names[ br_name ] = 1u;
      }
    } // weight branch names

    // Add the weights TChain as a friend of the event chain for easy
    // manipulation of both
    input_chain.AddFriend( &weights_ch );

    //// DEBUG: test with CV and unweighted only
    //weight_names.clear();
    //weight_names[ "weight_TunedCentralValue_UBGenie" ] = 1u;

  } // MC event weights are needed

  // Unconditionally add a special "unweighted" universe to the map of weight
  // names
  weight_names[ UNWEIGHTED_ID ] = 1u;

  for ( const auto& pair : weight_names ) {
    std::cout << pair.first << ' ' << pair.second << '\n';
  }

  // Find relevant entries in the input TChain for each of the reco space bins
  resp_mat.build_reco_entry_lists();

  // If we're working with MC events, then also find the relevant entries in
  // the input TChain for each of the true space bins
  if ( nf_type != NtupleFileType::kOnBNB
    && nf_type != NtupleFileType::kExtBNB )
  {
    resp_mat.build_true_entry_lists();
  }

  // Configure the output TTree
  // TODO: set output file name in a smarter way
  // Get the basename of the input file using the trick described here:
  // https://stackoverflow.com/a/24386991
  std::string input_basename = input_file_name.substr(
    input_file_name.find_last_of('/') + 1 );

  std::string out_file_name = "/uboone/data/users/gardiner/ntuples-stv"
    "/resp/respmat-" + input_basename;

  TFile out_file( out_file_name.c_str(), "recreate" );

  double smear, smear_err2, swt, swt_err2, swr, swr_err2;
  unsigned true_bin, reco_bin, universe;
  std::string weight_id;

  TTree* out_tree = new TTree( "resp_mat_tree", "response matrix elements" );
  out_tree->Branch( "smear", &smear, "smear/D" );
  out_tree->Branch( "smear_err2", &smear_err2, "smear_err2/D" );
  out_tree->Branch( "swt", &swt, "swt/D" );
  out_tree->Branch( "swt_err2", &swt_err2, "swt_err2/D" );
  out_tree->Branch( "swr", &swr, "swr/D" );
  out_tree->Branch( "swr_err2", &swr_err2, "swr_err2/D" );
  out_tree->Branch( "true_bin", &true_bin, "true_bin/i" );
  out_tree->Branch( "reco_bin", &reco_bin, "reco_bin/i" );
  out_tree->Branch( "universe", &universe, "universe/i" );
  out_tree->Branch( "weight_id", &weight_id );

  size_t num_true_bins = resp_mat.true_bins().size();
  size_t num_reco_bins = resp_mat.reco_bins().size();

  for ( const auto& wn_pair : weight_names ) {

    weight_id = wn_pair.first;
    unsigned num_universes = wn_pair.second;

    std::string cut_base = get_weight_cut_expr( weight_id );

    for ( universe = 0u; universe < num_universes; ++universe ) {

      std::cout << weight_id << ' ' <<  universe << '\n';

      std::string cut_expr = cut_base + std::to_string( universe );

      // Add expressions to limit the systematic weights to sane values (e.g.,
      // finite, positive numbers)
      std::string safe_cut_expr = "std::isfinite(" + cut_expr + ") && "
        + cut_expr + " <= " + MAX_WEIGHT + " ? " + cut_expr + " : 1";

      std::cout << "CUT EXPR = \"" << cut_expr << "\"\n";

      for ( true_bin = 0u; true_bin < num_true_bins; ++true_bin ) {

        std::cout << "  true bin = " << true_bin << '\n';

        // If we're working with a data file, just zero out the
        // summed event weights for each true bin. We can't fill
        // those without truth information!
        if ( nf_type == NtupleFileType::kOnBNB
          || nf_type == NtupleFileType::kExtBNB )
        {
          swt = 0.;
          swt_err2 = 0.;
        }
        // Otherwise, process the event weights in true space normally
        else {
          const auto& tel = resp_mat.true_entry_lists().at( true_bin );

          // Get the sum of the event weights in the current true bin
          sum_event_weights( input_chain, safe_cut_expr,
            tel.get(), swt, swt_err2 );
        }

        for ( reco_bin = 0u; reco_bin < num_reco_bins; ++reco_bin ) {

          std::cout << "    reco bin = " << reco_bin << '\n';

          const auto& rel = resp_mat.reco_entry_lists().at( reco_bin );

          // Compute the sum of all of the event weights in the current reco
          // bin.
          sum_event_weights( input_chain, safe_cut_expr,
            rel.get(), swr, swr_err2 );

          // If we're working with a real data file, also zero out the smearing
          // matrix entry, i.e., the summed event weights for the 2D bin in
          // (true, reco) space.
          if ( nf_type == NtupleFileType::kOnBNB
            || nf_type == NtupleFileType::kExtBNB )
          {
            smear = 0.;
            smear_err2 = 0.;
          }
          // Otherwise proceed normally with the summation of the event weights
          else {
            const auto& tel = resp_mat.true_entry_lists().at( true_bin );

            // Compute the intersection of the entry lists for the current true
            // and reco bins
            auto tr_el = tel * rel;

            // Get the sum of the event weights in the current 2D bin
            // in (true, reco) space
            sum_event_weights( input_chain, safe_cut_expr,
              tr_el.get(), smear, smear_err2 );
          }

          std::cout << "      smear = " << smear << '\n';

          out_tree->Fill();

        } // reco bin
      } // true bin
    } // universe
  } // weight ID

  out_tree->Write();
}

void test_rmm() {
  ROOT::EnableImplicitMT();
  //make_response_matrix( "/uboone/data/users/gardiner/ntuples-stv/stv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run2_reco2_D1D2_reco2.root" );
  //make_response_matrix( "/uboone/data/users/gardiner/ntuples-stv/stv-run1_neutrinoselection_filt_numu_ALL.root" );

  std::ifstream list_file( "foo2" );
  std::string file_name;
  while ( std::getline(list_file, file_name, '\n') ) {
    make_response_matrix( file_name );
  }

}
