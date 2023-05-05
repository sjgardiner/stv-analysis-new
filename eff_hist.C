// Efficiency histogram plotter

//const std::string signal_cuts = "mc_is_signal";
 // " || ( category == 6 && @mc_p3_p_vec.size() > 0 )";
//const std::string selection_cuts = "sel_CCNp0pi";
//const std::string selection_cuts = "sel_nu_mu_cc"
//" && sel_no_reco_showers && sel_muon_above_threshold"
//" && sel_has_p_candidate && sel_passed_proton_pid_cut"
//" && sel_protons_contained ";

const std::string signal_cuts = "mc_vertex_in_FV && mc_neutrino_is_numu"
  " && mc_no_fs_mesons && Sum$( mc_pdg == 13 ) == 1"
  " && Sum$( mc_pdg == 2212 ) > 0";
  //" && mc_muon_in_mom_range && mc_lead_p_in_mom_range"

const std::string selection_cuts = "sel_nu_mu_cc && sel_muon_contained"
  " && sel_muon_quality_ok"
  //" && sel_muon_passed_mom_cuts && sel_lead_p_passed_mom_cuts"
  " && sel_no_reco_showers && sel_has_p_candidate && sel_protons_contained "
  " && sel_passed_proton_pid_cut";

void plot_eff_hist( TTree& stv_tree, const std::string& branch,
  const std::string& variable_name, const std::string& hist_name,
  const std::string& unit_name, int num_bins, double x_min, double x_max)
{
  // For computing efficiencies and purities, we need to only use MC events.
  // Unconditionally add this requirement to the cuts defined above.
  std::string signal = signal_cuts + " && is_mc ";
  std::string selection = selection_cuts + " && is_mc ";

  std::string unit_part;
  if ( !unit_name.empty() ) unit_part = " [" + unit_name + ']';
  std::string eff_title = "Efficiency in true " + variable_name + " bins; "
    + variable_name + unit_part + "; Efficiency";
  TH1D* eff_hist = new TH1D( hist_name.c_str(), eff_title.c_str(), num_bins,
    x_min, x_max );

  for ( int b = 1; b <= num_bins; ++b ) {
    double xlow = eff_hist->GetBinLowEdge( b );
    double xhigh = eff_hist->GetBinLowEdge( b + 1 );

    // Define a new set of cuts requiring the event to fall inside the current
    // bin
    std::stringstream ss;
    ss << branch << " >= " << xlow << " && "
      << branch << " < " << xhigh;
    std::string bin_cut = ss.str();

    // Combine the bin cuts with the signal and selection cuts
    std::string bin_signal = signal + " && " + bin_cut;
    std::string bin_selection = selection + " && " + bin_cut;

    // These are actually integer counts, but since we will use them below to
    // compute ratios, intrinsically cast them to double-precision values for
    // convenience.
    double num_signal = stv_tree.Draw( "", bin_signal.c_str(), "goff" );
    double num_selected_signal = stv_tree.Draw( "",
      (bin_signal + " && " + selection).c_str(), "goff" );

    // Compute the efficiency for the current bin
    double eff = 0.;
    double eff_stat_err = 0.;
    if ( num_signal > 0. && num_selected_signal > 0. ) {
      eff = num_selected_signal / num_signal;
      eff_stat_err = eff * std::sqrt( (1. / num_selected_signal)
        + (1. / num_signal) );
    }

    eff_hist->SetBinContent( b, eff );
    eff_hist->SetBinError( b, eff_stat_err );
  }

  TCanvas* c1 = new TCanvas;
  eff_hist->SetStats( false );
  eff_hist->GetYaxis()->SetRangeUser( 0., 1. );

  c1->SetBottomMargin(0.15);

  eff_hist->SetLineWidth(4);
  eff_hist->SetLineColor(kBlack);
  eff_hist->GetYaxis()->SetRangeUser(0., 1.);
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
  stv_ch.Add( "/uboone/data/users/gardiner/ntuples-stv/stv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root" );
  stv_ch.Add( "/uboone/data/users/gardiner/ntuples-stv/stv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run2_reco2_D1D2_reco2.root" );
  stv_ch.Add( "/uboone/data/users/gardiner/ntuples-stv/stv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run3_reco2_G_reco2.root" );

  //// Fake data samples (NuWro)
  //stv_ch.Add( "/uboone/data/users/gardiner/ntuples-stv/stv-high_stat_prodgenie_bnb_nu_overlay_DetVar_Run1_NuWro_reco2_reco2.root" );
  //stv_ch.Add( "/uboone/data/users/gardiner/ntuples-stv/stv-high_stat_prodgenie_bnb_nu_overlay_DetVar_Run2_NuWro_reco2_reco2.root" );
  //stv_ch.Add( "/uboone/data/users/gardiner/ntuples-stv/stv-high_stat_prodgenie_bnb_nu_overlay_DetVar_Run3_NuWro_reco2_reco2.root" );

 plot_eff_hist( stv_ch, "mc_p3_mu.Mag()", "p_{#mu}^{true}", "pmutrue", "GeV", 80, 0., 2.);

 //plot_eff_hist( stv_ch, "mc_p3_lead_p.Mag()", "p_{p}^{true}", "pptrue", "GeV", 80, 0., 1.2 );

  //plot_eff_hist( stv_ch, "mc_p3_mu.CosTheta()", "cos#theta_{#mu}^{true}", "cthmutrue", "", 25, -1., 1.);
  //plot_eff_hist( stv_ch, "mc_p3_mu.Phi()", "#phi_{#mu}^{true}", "phimutrue", "rad", 25, 0., M_PI);

  //plot_eff_hist( stv_ch, "mc_p3_lead_p.Mag()", "p_{lead p}^{true}", "pptrue", "GeV", 25, 0., 2.);
  //plot_eff_hist( stv_ch, "mc_p3_lead_p.CosTheta()", "cos#theta_{lead p}^{true}", "cthptrue", "", 25, -1., 1.);
  //plot_eff_hist( stv_ch, "mc_p3_lead_p.Phi()", "#phi_{lead p}^{true}", "phiptrue", "rad", 25, 0., M_PI);

  //plot_eff_hist( stv_ch, "mc_p3_lead_p.X() / mc_p3_lead_p.Mag()", "px_{lead p}^{true} / p_{lead p}^{true}", "pxp", "", 25, 0., 1.);
  //plot_eff_hist( stv_ch, "mc_p3_lead_p.Y() / mc_p3_lead_p.Mag()", "py_{lead p}^{true} / p_{lead p}^{true}", "pyp", "", 25, 0., 1.);
  //plot_eff_hist( stv_ch, "mc_p3_lead_p.Z() / mc_p3_lead_p.Mag()", "pz_{lead p}^{true} / p_{lead p}^{true}", "pzp", "", 25, 0., 1.);

  //plot_eff_hist( stv_ch, "mc_delta_pT", "#deltap_{T}", "deltapT", "GeV", 25, 0., 2.);
  //plot_eff_hist( stv_ch, "mc_delta_phiT", "#delta#phi_{T}", "deltaphiT", "rad", 25, 0., M_PI);
  //plot_eff_hist( stv_ch, "mc_delta_alphaT", "#delta#alpha_{T}", "deltaalphaT", "rad", 25, 0., M_PI);
  //plot_eff_hist( stv_ch, "mc_delta_pL", "#deltap_{L}", "deltapL", "rad", 25, 0., 2.);
  //plot_eff_hist( stv_ch, "mc_pn", "p_{n}", "pn", "rad", 25, 0., 2.);

}
