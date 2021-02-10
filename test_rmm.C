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

  std::cout << "POT = " << resp_mat.pot() << '\n';

  // TODO retrieve weight file programmatically
  TChain weights_ch( "weight_tree" );
  weights_ch.Add( "/uboone/data/users/gardiner/ntuples-stv/weight_dumps/weights-stv-prodgenie_dirt_overlay_v08_00_00_35_all_run2_reco2_reco2.root" );

  // Build a map of weight names (ignoring universe indices) by looping
  // over all of the branches and trimming off trailing digits. We'll use the
  // map values to keep track of universe counts for each unique weight name.
  std::map< std::string, unsigned > weight_names;

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
  }

  for ( const auto& pair : weight_names ) {
    std::cout << pair.first << ' ' << pair.second << '\n';
  }

  // Add the weights TChain as a friend of the event chain for easy
  // manipulation of both
  auto& input_chain = resp_mat.input_chain();
  input_chain.AddFriend( &weights_ch );

  // Find relevant entries in the input TChain for each of the true and reco
  // space bins
  resp_mat.build_entry_lists();

  // Configure the output TTree
  // TODO: set output file name
  TFile out_file( "/uboone/data/users/gardiner/resp_mat_test.root",
    "recreate" );

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

  // TESTING CODE
  size_t num_true_bins = resp_mat.true_bins().size();
  size_t num_reco_bins = resp_mat.reco_bins().size();

  // Add a special "unweighted" universe to the map of weight names
  weight_names[ UNWEIGHTED_ID ] = 1u;

  for ( const auto& wn_pair : weight_names ) {

    weight_id = wn_pair.first;
    unsigned num_universes = wn_pair.second;

    std::string cut_base = get_weight_cut_expr( weight_id );

    for ( universe = 0u; universe < num_universes; ++universe ) {

      // DEBUG
      if ( universe > 2 ) break;

      std::cout << weight_id << ' ' <<  universe << '\n';

      std::string cut_expr = cut_base + std::to_string( universe );

      // Add expressions to limit the systematic weights to sane values (e.g.,
      // finite, positive numbers)
      std::string safe_cut_expr = "std::isfinite(" + cut_expr + ") && "
        + cut_expr + " <= " + MAX_WEIGHT + " ? " + cut_expr + " : 1";

      std::cout << "CUT EXPR = \"" << cut_expr << "\"\n";

      for ( true_bin = 0u; true_bin < num_true_bins; ++true_bin ) {

        std::cout << "  true bin = " << true_bin << '\n';

        const auto& tel = resp_mat.true_entry_lists().at( true_bin );

        // Get the sum of the event weights in the current true bin
        sum_event_weights( input_chain, safe_cut_expr,
          tel.get(), swt, swt_err2 );

        for ( reco_bin = 0u; reco_bin < num_reco_bins; ++reco_bin ) {

          std::cout << "    reco bin = " << reco_bin << '\n';

          const auto& rel = resp_mat.reco_entry_lists().at( reco_bin );

          // For the first true bin only (to avoid redundant function calls in
          // later iterations of the reco bin loop), compute the sum of the
          // event weights in the current reco bin.
          if ( true_bin == 0u ) {
            sum_event_weights( input_chain, safe_cut_expr,
              rel.get(), swr, swr_err2 );
          }

          // Compute the intersection of the entry lists for the current true
          // and reco bins
          auto tr_el = tel * rel;

          // Get the sum of the event weights in the current 2D bin
          // in (true, reco) space
          sum_event_weights( input_chain, safe_cut_expr,
            tr_el.get(), smear, smear_err2 );

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
  make_response_matrix( "/uboone/data/users/gardiner/ntuples-stv/stv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run2_reco2_D1D2_reco2.root" );
}
