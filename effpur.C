// Efficiencies and purities for the MCC9 CCNp0pi STV analysis

//  const std::string& branch, const std::string& variable_name,
//  const std::string& hist_name, int num_bins, double x_min, double x_max)

// For a given analysis TTree (containing the events to use), compute the
// efficiency and purity given definitions for what events count as signal
// (based on MC truth information) and what events are selected (based on
// reconstructed information) in terms of our analysis TTree branch variables.
// Store the results in the output variables eff and pur.
void compute_eff_pur( TTree& stv_tree, const std::string& signal_cuts,
  const std::string& selection_cuts, double& eff, double& pur )
{
  // For computing efficiencies and purities, we need to only use MC events.
  // Unconditionally add this requirement to the cuts defined above.
  std::string signal = signal_cuts + " && is_mc && genie_ok ";
  std::string selection = selection_cuts + " && is_mc && genie_ok ";

  // These are actually integer counts, but since we will use them below to
  // compute ratios, intrinsically cast them to double-precision values for
  // convenience.
  double num_signal = stv_tree.Draw( "", signal.c_str(), "goff" );
  double num_selected = stv_tree.Draw( "", selection.c_str(), "goff" );
  double num_selected_signal = stv_tree.Draw( "",
    (signal + " && " + selection).c_str(), "goff" );

  eff = num_selected_signal / num_signal;
  pur = num_selected_signal / num_selected;

  //std::cout << "signal = " << num_signal << '\n';
  //std::cout << "selected = " << num_selected << '\n';
  //std::cout << "selected_signal = " << num_selected_signal << '\n';

}

const std::string signal_cuts = "mc_is_signal";
const std::string selection_cuts = "sel_CCNp0pi";


void effpur() {

  const std::vector< std::string > signal_defs = { "1",
    "mc_vertex_in_FV",
    "mc_vertex_in_FV && mc_neutrino_is_numu",
    "mc_vertex_in_FV && mc_neutrino_is_numu && mc_pmu_above_threshold",
    "mc_vertex_in_FV && mc_neutrino_is_numu && mc_pmu_above_threshold"
      " && mc_has_p_above_threshold",
    "mc_vertex_in_FV && mc_neutrino_is_numu && mc_pmu_above_threshold"
      " && mc_has_p_above_threshold && mc_no_fs_pi0",
    "mc_is_signal" };

  const std::vector< std::string > selection_defs = { "1",
  "sel_nu_mu_cc",
  "sel_nu_mu_cc && sel_no_reco_showers",
  "sel_nu_mu_cc && sel_no_reco_showers && sel_has_single_muon_candidate",
  "sel_nu_mu_cc && sel_no_reco_showers && sel_has_single_muon_candidate"
    " && sel_muon_above_threshold",
  "sel_nu_mu_cc && sel_no_reco_showers && sel_has_single_muon_candidate"
    " && sel_muon_above_threshold && sel_has_p_candidate",
  "sel_nu_mu_cc && sel_no_reco_showers && sel_has_single_muon_candidate"
    " && sel_muon_above_threshold && sel_has_p_candidate"
    " && sel_protons_contained ",
  "sel_nu_mu_cc && sel_no_reco_showers && sel_has_single_muon_candidate"
    " && sel_muon_above_threshold && sel_has_p_candidate"
    " && sel_protons_contained && sel_passed_proton_pid_cut",
  "sel_nu_mu_cc && sel_no_reco_showers && sel_has_single_muon_candidate"
    " && sel_muon_above_threshold && sel_has_p_candidate"
    " && sel_protons_contained && sel_passed_proton_pid_cut"
    " && sel_lead_p_passed_hits_cut",
  "sel_CCNp0pi" };

  TChain stv_ch( "stv_tree" );
  stv_ch.Add( "*stv.root" );

  size_t num_points = selection_defs.size();
  TGraph* eff_graph = new TGraph( num_points );
  TGraph* pur_graph = new TGraph( num_points );

  std::string signal = signal_defs.back();
  double eff, pur;
  for ( size_t k = 0u; k < num_points; ++k  ) {

    const auto& selection = selection_defs.at( k );
    compute_eff_pur( stv_ch, signal, selection, eff, pur );

    eff_graph->SetPoint( k, k + 1, eff );
    pur_graph->SetPoint( k, k + 1, pur );

    std::cout << "selection = " << selection << '\n';
    std::cout << "eff = " << eff << '\n';
    std::cout << "pur = " << pur << '\n';
    std::cout << "\n\n";
  }

  TCanvas* c1 = new TCanvas;
  eff_graph->SetLineColor(kBlue);
  eff_graph->SetMarkerColor(kBlue);
  eff_graph->SetLineWidth(3);
  eff_graph->SetMarkerStyle(20);
  eff_graph->GetYaxis()->SetRangeUser( 0., 1. );
  eff_graph->Draw( "alp" );

  pur_graph->SetLineColor(kRed);
  pur_graph->SetMarkerColor(kRed);
  pur_graph->SetLineWidth(3);
  pur_graph->SetMarkerStyle(20);
  pur_graph->Draw("same lp");

  TLegend* lg = new TLegend(0.65, 0.4, 0.85, 0.6);
  lg->AddEntry( eff_graph, "efficiency", "lp" );
  lg->AddEntry( pur_graph, "purity", "lp" );
  lg->Draw("same");
}
