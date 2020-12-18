// STV analysis includes
#include "EventCategory.hh"

const std::string BNB_DATA_FILE = "/uboone/data/users/gardiner/ntuples-stv/stv-data_bnb_mcc9.1_v08_00_00_25_reco2_C1_beam_good_reco2_5e19.root";
const std::string EXT_BNB_DATA_FILE = "/uboone/data/users/gardiner/ntuples-stv/stv-data_extbnb_mcc9.1_v08_00_00_25_reco2_C_all_reco2.root";

constexpr double EXT_OFF_DATA = 65498807.0;
constexpr double E1DCNT_ON_DATA = 10080350.0;
constexpr double POT_ON_DATA = 4.54e+19;
constexpr double POT_OFF_DATA = (EXT_OFF_DATA / E1DCNT_ON_DATA) * POT_ON_DATA;

// Cuts to use when filling histograms with selected events
const std::string selection = "sel_CCNp0pi"; //" && @p3_p_vec.size() > 2";
//const std::string selection = "sel_nu_mu_cc && sel_has_muon_candidate"; // && reco_muon_contained";

// Weight to apply to MC events when filling histograms
//const std::string mc_event_weight = "1.";
//const std::string mc_event_weight = "spline_weight";
const std::string mc_event_weight = "spline_weight * (std::isfinite("
  "tuned_cv_weight) && tuned_cv_weight <= 100. ? tuned_cv_weight : 1)";

void make_plots(const std::string& hist_name_prefix, const std::string& branch,
  const std::string& var_name, double xmin, double xmax, int Nbins,
  const std::vector<std::string>& mc_file_names)
{
  const EventCategoryInterpreter& eci = EventCategoryInterpreter::Instance();

  auto* c1 = new TCanvas; //( "c1", "", 800, 600 );
  //c1->SetLeftMargin( 0.12 );
  //c1->SetBottomMargin( 1.49 );

  std::stringstream temp_ss;
  temp_ss << POT_ON_DATA;

  std::string plot_title = var_name + ", MCC9, Run 1; " + var_name
    + "; #frac{Selected Events}{" + temp_ss.str() + " POT}";

  TFile off_data_file( EXT_BNB_DATA_FILE.c_str(), "read" );
  TTree* off_data_tree = nullptr;
  off_data_file.GetObject( "stv_tree", off_data_tree );

  std::string off_data_hist_name = hist_name_prefix + "-ext";
  TH1D* off_data_hist = new TH1D( off_data_hist_name.c_str(),
    plot_title.c_str(), Nbins, xmin, xmax );
  off_data_tree->Draw( (branch + " >> " + off_data_hist_name).c_str(),
    selection.c_str() );
  off_data_hist->Scale( POT_ON_DATA / POT_OFF_DATA );
  off_data_hist->SetDirectory( nullptr );

  eci.set_ext_histogram_style( off_data_hist );

  TFile on_data_file( BNB_DATA_FILE.c_str(), "read" );
  TTree* on_data_tree = nullptr;
  on_data_file.GetObject( "stv_tree", on_data_tree );

  std::string on_data_hist_name = hist_name_prefix + "-on";
  TH1D* on_data_hist = new TH1D( on_data_hist_name.c_str(),
    plot_title.c_str(), Nbins, xmin, xmax);
  on_data_tree->Draw( (branch + " >> " + on_data_hist_name).c_str(),
    selection.c_str() );

  on_data_hist->Scale( 1. );
  on_data_hist->SetDirectory( nullptr );

  eci.set_bnb_data_histogram_style( on_data_hist );

  // Initialize empty stacked histograms by MC event category
  std::map< EventCategory, TH1D* > mc_hists;
  // Loop over all MC event categories
  for ( const auto& pair : eci.label_map() ) {
    EventCategory cat = pair.first;
    std::string cat_label = pair.second;

    std::string temp_mc_hist_name = hist_name_prefix + "-temp_mc"
      + std::to_string( cat );

    TH1D* temp_mc_hist = new TH1D( temp_mc_hist_name.c_str(),
      plot_title.c_str(), Nbins, xmin, xmax );

    mc_hists[ cat ] = temp_mc_hist;

    temp_mc_hist->SetDirectory( nullptr );
    eci.set_mc_histogram_style( cat, temp_mc_hist );
  }

  // Loop over the different MC files and collect their contributions.
  // We have to handle them separately in order to get the POT normalization
  // correct.
  int dummy_counter = 0;
  for ( const auto& mc_file_name : mc_file_names ) {

    // Get the POT values from the current input MC file
    TFile temp_mc_file( mc_file_name.c_str(), "read" );
    TParameter<float>* temp_pot = nullptr;
    temp_mc_file.GetObject( "summed_pot", temp_pot );

    double mc_pot = temp_pot->GetVal();

    // Use a temporary TChain to analyze the MC events
    TChain mc_ch( "stv_tree" );
    mc_ch.Add( mc_file_name.c_str() );

    // Add this file's contribution to the stacked histograms by MC event
    // category
    for ( const auto& pair : eci.label_map() ) {

      EventCategory ec = pair.first;

      std::string temp_mc_hist_name = hist_name_prefix + "-temp_mc"
        + std::to_string(ec) + "_number" + std::to_string(dummy_counter);

      ++dummy_counter;

      TH1D* temp_mc_hist = new TH1D( temp_mc_hist_name.c_str(),
        plot_title.c_str(), Nbins, xmin, xmax );

      mc_ch.Draw( (branch + " >> " + temp_mc_hist_name).c_str(),
        (mc_event_weight + "*(" + selection + " && category == "
        + std::to_string(ec) + ')').c_str() );

      // Scale to the same exposure as the beam on data
      temp_mc_hist->Scale( POT_ON_DATA / mc_pot );

      // Add this histogram's contribution (now properly scaled) to the total
      mc_hists.at( ec )->Add( temp_mc_hist );

      // We don't need the temporary histogram anymore, so just get rid of it
      delete temp_mc_hist;
    }
  } // loop over MC files

  TPad* pad1 = new TPad( "pad1", "", 0.0, 0.23, 1.0, 1.0 );
  pad1->SetBottomMargin( 0 );
  pad1->SetRightMargin( 0.06 );
  pad1->SetLeftMargin( 0.13 );
  pad1->SetGridx();
  pad1->Draw();
  pad1->cd();

  //on_data_hist->GetYaxis()->SetRangeUser(0., 300.);
  on_data_hist->Draw( "E1" );

  // Stack of categorized MC predictions plus extBNB contribution
  THStack* stacked_hist = new THStack( "mc", "" );

  // Sum all contributions into this TH1D so that we can get the overall
  // statistical uncertainty easily
  TH1D* stat_err_hist = new TH1D(
    ("stat_err_hist_" + hist_name_prefix).c_str(), "", Nbins, xmin, xmax );

  stacked_hist->Add( off_data_hist );
  stat_err_hist->Add( off_data_hist );

  for ( auto citer = mc_hists.crbegin(); citer != mc_hists.crend();
    ++citer )
  {
    TH1D* hist = citer->second;
    stacked_hist->Add( hist );
    stat_err_hist->Add( hist );
  }

  stacked_hist->Draw( "hist same" );

  on_data_hist->Draw( "E1 same" );

  eci.set_stat_err_histogram_style( stat_err_hist );
  stat_err_hist->Draw( "E2 same" );

  TLegend* lg = new TLegend( 0.64, 0.42, 0.94, 0.85 );
  lg->AddEntry( on_data_hist, "Data (beam on)", "lp" );
  lg->AddEntry( stat_err_hist, "Statistical uncertainty", "f" );


  double total_events = stat_err_hist->Integral();
  for ( const auto& pair : eci.label_map() ) {
    EventCategory ec = pair.first;
    std::string label = pair.second;

    TH1* category_hist = mc_hists.at( ec );

    // Use TH1::Integral() to account for CV reweighting correctly
    double events_in_category = category_hist->Integral();
    double category_percentage = events_in_category * 100. / total_events;

    std::string cat_pct_label = Form( "%.2f%#%", category_percentage );

    lg->AddEntry( category_hist, (label + ", " + cat_pct_label).c_str(), "f" );
  }

  double beam_off_events = off_data_hist->Integral();
  double beam_off_percentage = beam_off_events * 100. / total_events;

  std::string off_pct_label = Form( "%.2f%#%", beam_off_percentage );

  lg->AddEntry( off_data_hist, ("Data (beam off), "
    + off_pct_label).c_str(), "f" );

  lg->SetBorderSize( 0 );

  lg->Draw( "same" );

  // Ratio plot
  c1->cd(); // Go back from pad1 to main canvas c1

  TPad* pad2 = new TPad( "pad2", "", 0, 0.01, 1.0, 0.23 );
  pad2->SetTopMargin( 0 );
  pad2->SetFrameFillStyle( 4000 );
  pad2->SetBottomMargin( 0.38 );
  pad2->SetRightMargin( 0.06 );
  pad2->SetLeftMargin( 0.13 );

  pad2->SetGridx();
  pad2->Draw();
  pad2->cd(); // change current pad to pad2

  // Ratio plot
  TH1D* h_ratio = dynamic_cast<TH1D*>( on_data_hist->Clone("h_ratio") );
  h_ratio->SetStats( false );
  h_ratio->Divide( stat_err_hist );
  h_ratio->SetLineWidth( 2 );
  h_ratio->SetLineColor( kBlack );
  h_ratio->SetMarkerStyle( kFullCircle );
  h_ratio->SetMarkerSize( 0.8 );
  h_ratio->SetTitle( "" );

  // x-axis
  h_ratio->GetXaxis()->SetTitle( on_data_hist->GetXaxis()->GetTitle() );
  h_ratio->GetXaxis()->CenterTitle( true );
  h_ratio->GetXaxis()->SetLabelSize( 0.12 );
  h_ratio->GetXaxis()->SetTitleSize( 0.18 );
  h_ratio->GetXaxis()->SetTickLength( 0.05 );
  h_ratio->GetXaxis()->SetTitleOffset( 0.9 );

  // y-axis
  h_ratio->GetYaxis()->SetTitle( "#frac{Beam ON}{Beam OFF + MC}" );
  h_ratio->GetYaxis()->CenterTitle( true );
  h_ratio->GetYaxis()->SetLabelSize( 0.08 );
  h_ratio->GetYaxis()->SetTitleSize( 0.085 );
  h_ratio->GetYaxis()->SetTitleOffset( 0.5 );

  h_ratio->Draw( "E1" );

  gStyle->SetGridColor( 17 );

  // Adjust y-axis
  double ratio_max = h_ratio->GetBinContent( h_ratio->GetMaximumBin() );
  double ratio_min = h_ratio->GetBinContent( h_ratio->GetMinimumBin() );

  //h_ratio->GetYaxis()->SetRangeUser( ratio_min - ratio_min*0.2, ratio_max - ratio_max*0.8 );
  h_ratio->SetMaximum( ratio_max + ratio_max*0.15 );
  h_ratio->SetMinimum( ratio_min - ratio_min*0.2 );

  gPad->Update();

  // Draw a horizontal dashed line at ratio == 1
  TLine* line = new TLine( h_ratio->GetXaxis()->GetXmin(), 1.0,
    h_ratio->GetXaxis()->GetXmax(), 1.0 );
  line->SetLineColor( kBlack );
  line->SetLineStyle( 9 ); // dashed
  line->Draw();

  c1->Update();
}

void plots() {

  std::vector<std::string> mc_file_names = {
"/uboone/data/users/gardiner/ntuples-stv/stv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root",
"/uboone/data/users/gardiner/ntuples-stv/stv-prodgenie_bnb_intrinsice_nue_uboone_overlay_mcc9.1_v08_00_00_26_run1_reco2_reco2.root",
"/uboone/data/users/gardiner/ntuples-stv/stv-prodgenie_bnb_dirt_overlay_mcc9.1_v08_00_00_26_run1_reco2_reco2.root"
  };

  //make_plots("cthmu", "p3_mu.CosTheta()", "cos#theta_{#mu}", -1., 1., 20, mc_file_names);

  make_plots("phimu", "p3_mu.Phi()", "#phi_{#mu}", 0., M_PI, 40, mc_file_names);

  //make_plots("cthp", "p3_lead_p.CosTheta()", "cos#theta_{p}", -1., 1., 20, mc_file_names);

  //make_plots( "sumTp",
  //  "Sum$(TMath::Sqrt(p3_p_vec.Mag2() + 0.93827208*0.93827208) - 0.93827208)",
  //  "#Sigma Tp", 0., 0.8, 15, mc_file_names );

  //make_plots("delta_alphaT", "delta_alphaT", "reco #delta#alpha_{T}", 0., M_PI, 15, mc_file_names);

  //make_plots( "pn", "pn", "reco p_{n} (GeV)", 0., 1., 25, mc_file_names );
  //make_plots("delta_phiT", "delta_phiT", "reco #delta#phi_{T}", 0., M_PI, 20, mc_file_names);

  //make_plots("trk_len_mu", "trk_len_v[muon_candidate_idx]", "muon candidate track length (uncontained only)", 0., 500., 15, mc_file_names);

  //make_plots("trk_mom_mu", "trk_range_muon_mom_v[muon_candidate_idx]", "muon candidate momentum (uncontained only)", 0., 1.5, 15, mc_file_names);

  //make_plots("trk_mom_mu", "trk_mcs_muon_mom_v[muon_candidate_idx]", "muon candidate momentum (contained only)", 0., 1.5, 15, mc_file_names);

  make_plots("pp", "p3_lead_p.Mag()", "reco p_{lead p} (GeV)", 0.2, 1.2, 20, mc_file_names);

  make_plots("delta_pT", "delta_pT", "reco #deltap_{T} (GeV)", 0., 1., 15, mc_file_names);

}
