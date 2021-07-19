// Efficiency histogram plotter
#include "../HistUtils.hh"

// Equivalent to mc_is_signal with phase space cuts removed for the muon
// and proton momenta
const std::string signal_cuts = "mc_vertex_in_FV && mc_neutrino_is_numu"
  " && mc_no_fs_mesons && category != 9";
// && mc_pmu_above_threshold
// && mc_lead_p_in_mom_range

// Equivalent to sel_CCNp0pi with phase space cuts removed for the muon
// and proton momenta
const std::string selection_cuts = "sel_nu_mu_cc && sel_no_reco_showers"
  " && sel_has_p_candidate && sel_passed_proton_pid_cut"
  " && sel_protons_contained";
//&& sel_muon_above_threshold
//&& sel_lead_p_passed_mom_cuts"

void plot_eff_hist( TTree& stv_tree, const std::string& branch,
  const std::string& variable_name, const std::string& hist_name,
  const std::string& unit_name, int num_bins, double x_min, double x_max,
  const std::string& mc_event_weight = DEFAULT_MC_EVENT_WEIGHT)
{
  // For computing efficiencies, we need to only use MC events. Unconditionally
  // add this requirement to the cuts defined above, just in case.
  std::string signal = signal_cuts + " && is_mc ";
  std::string selection = selection_cuts + " && is_mc ";

  std::string unit_part;
  if ( !unit_name.empty() ) unit_part = " [" + unit_name + ']';
  //std::string eff_title = "Efficiency in true " + variable_name + " bins; "
  //  + variable_name + unit_part + "; Efficiency";

  std::string eff_title = ";" + variable_name + unit_part + "; efficiency";

  // Create a histogram to store the efficiency values
  TH1D* eff_hist = new TH1D( hist_name.c_str(), eff_title.c_str(), num_bins,
    x_min, x_max );

  // Create another to use to temporarily store the denominator of the
  // efficiency
  std::string denom_hist_name = hist_name + "_denom";
  TH1D* eff_denom_hist = new TH1D( denom_hist_name.c_str(),
    eff_title.c_str(), num_bins, x_min, x_max );

  // Fill the histograms with event counts from the TTree containing MC
  // samples. Selected signal events go in the first histogram
  stv_tree.Draw( (branch + " >> " + hist_name).c_str(),
    (mc_event_weight + " * (" + signal + " && "
    + selection + ')').c_str(), "goff" );

  // All signal events go in the second histogram
  stv_tree.Draw( (branch + " >> " + denom_hist_name).c_str(),
    (mc_event_weight + " * (" + signal + ')').c_str(), "goff" );

  // Now we have both pieces. Divide the numerator by the denominator
  // to compute the efficiency
  eff_hist->Divide( eff_denom_hist );

  TCanvas* c1 = new TCanvas;
  eff_hist->SetStats( false );
  eff_hist->GetYaxis()->SetRangeUser( 0., 0.45 );

  c1->SetBottomMargin(0.15);

  eff_hist->SetLineWidth(4);
  eff_hist->SetLineColor(kBlack);
  eff_hist->GetYaxis()->SetRangeUser(0., 0.45);
  eff_hist->GetYaxis()->CenterTitle(true);
  eff_hist->GetYaxis()->SetTitleOffset(0.95);
  eff_hist->GetYaxis()->SetTitleSize(0.05);

  eff_hist->GetXaxis()->SetTitleOffset(1.2);
  eff_hist->GetXaxis()->CenterTitle(true);
  eff_hist->GetXaxis()->SetTitleSize(0.05);

  eff_hist->SetStats( false );
  eff_hist->Draw( );

  //c1->SaveAs( ( "eff_hist_" + hist_name + ".jpg").c_str() );
}

void eff_hist() {

  TChain stv_ch( "stv_tree" );
  stv_ch.Add( "/uboone/data/users/gardiner/ntuples-stv-MCC9InternalNote/stv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root" );
  stv_ch.Add( "/uboone/data/users/gardiner/ntuples-stv-MCC9InternalNote/stv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run2_reco2_D1D2_reco2.root" );
  stv_ch.Add( "/uboone/data/users/gardiner/ntuples-stv-MCC9InternalNote/stv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run3_reco2_G_reco2.root" );

  plot_eff_hist( stv_ch, "mc_p3_mu.Mag()", "p_{#mu}^{true}", "pmutrue",
    "GeV/c", 80, 0., 2.);

  plot_eff_hist( stv_ch, "mc_p3_lead_p.Mag()", "p_{p}^{true}",
    "pptrue", "GeV/c", 120, 0., 1.5 );
}
