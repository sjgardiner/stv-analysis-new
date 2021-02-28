#include "ResponseMatrixMaker.hh"

void make_config() {

  //std::vector< double > dpT_edges = { 0., 0.1, 0.2, 0.3, 0.4, 0.5,
  //  0.6, 0.7, 0.8 };

  std::vector<double> dpT_edges = { 0., 0.0533333, 0.106667, 0.16, 0.213333, 0.266667, 0.32,
    0.373333, 0.426667, 0.48, 0.533333, 0.586667, 0.64, 0.693333, 0.746667, 0.8 };

  std::vector< int > signal_categories = { 1, 2, 3, 4 };

  std::string selection = "sel_CCNp0pi";
  std::string signal_def = "mc_is_signal";

  std::vector< std::string > background_defs = {
    "category == 5", "category == 6", "category == 7", "category == 8",
    "category == 9", "category == 10", "category == 11"
  };

  std::vector< TrueBin > true_bins;
  std::vector< RecoBin > reco_bins;

  // Configure kinematic limits for all of the signal bins
  for ( size_t pTbin = 0u; pTbin < dpT_edges.size() - 1u; ++pTbin ) {

    std::string pTlow_str = std::to_string( dpT_edges.at(pTbin) );
    std::string pThigh_str = std::to_string( dpT_edges.at(pTbin + 1u) );

    for ( const int& cat : signal_categories ) {

      std::string cat_str = "category == " + std::to_string( cat );

      std::string true_bin_def = signal_def + " && mc_delta_pT >= " + pTlow_str
        + " && mc_delta_pT < " + pThigh_str + " && " + cat_str;

      true_bins.emplace_back( true_bin_def, kSignalTrueBin );
    }

    std::string reco_bin_def = selection + " && delta_pT >= " + pTlow_str
      + " && delta_pT < " + pThigh_str;

    reco_bins.emplace_back( reco_bin_def );
  }

  // Add true bins for the background categories of interest
  for ( const auto& bdef : background_defs ) {
    true_bins.emplace_back( bdef, kBackgroundTrueBin );
  }

  // Dump this information to the output file
  std::ofstream out_file( "myconfig.txt" );
  out_file << "stv_tree\n";
  out_file << true_bins.size() << '\n';
  for ( const auto& tb : true_bins ) out_file << tb << '\n';

  out_file << reco_bins.size() << '\n';
  for ( const auto& rb : reco_bins ) out_file << rb << '\n';
}
