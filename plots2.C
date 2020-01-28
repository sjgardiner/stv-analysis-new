constexpr double EXT_OFF_DATA = 13625641.0;
constexpr double E1DCNT_ON_DATA = 10384652.0;
constexpr double POT_ON_DATA = 4.449e+19;
constexpr double POT_OFF_DATA = (EXT_OFF_DATA / E1DCNT_ON_DATA) * POT_ON_DATA;

const std::string selection = "sel_CCNp0pi";
//const std::string selection = "sel_CCNp0pi * spline_weight"
// " * (std::isfinite(tuned_cv_weight) && tuned_cv_weight <= 100. ? tuned_cv_weight : 1)";

std::string cat_to_label( int category ) {
  if ( category == 1 ) return "signal CCQE";
  if ( category == 2 ) return "signal CCMEC";
  if ( category == 3 ) return "signal CCRES";
  if ( category == 4 ) return "signal other";
  if ( category == 5 ) return "out FV";
  if ( category == 6 ) return "NC";
  if ( category == 7 ) return "#nu_{#mu} CCN#pi";
  if ( category == 8 ) return "#nu_{e} CC";
  if ( category == 9 ) return "other";
  return "other";
}

void make_plots(const std::string& hist_name_prefix, const std::string& branch,
  const std::string& var_name, double xmin, double xmax, int Nbins,
  const std::vector<std::string>& mc_file_names)
{
  auto* c1 = new TCanvas;
  c1->SetLeftMargin(0.12);
  c1->SetBottomMargin(0.13);

  std::stringstream temp_ss;
  temp_ss << POT_ON_DATA;

  std::string plot_title = var_name + ", MCC9, Run 1; " + var_name + "; events / "
    + temp_ss.str() + " POT";

  TFile off_data_file( "off_data_stv.root", "read" );
  TTree* off_data_tree = nullptr;
  off_data_file.GetObject( "stv_tree", off_data_tree );

  std::string off_data_hist_name = hist_name_prefix + "-ext";
  TH1D* off_data_hist = new TH1D( off_data_hist_name.c_str(), plot_title.c_str(),
    Nbins, xmin, xmax );
  off_data_tree->Draw((branch + " >> " + off_data_hist_name).c_str(), selection.c_str());
  off_data_hist->Scale(POT_ON_DATA / POT_OFF_DATA);
  off_data_hist->SetFillColor( 44 );
  off_data_hist->SetStats(false);
  off_data_hist->SetDirectory( nullptr );

  TFile on_data_file( "on_data_stv.root", "read" );
  TTree* on_data_tree = nullptr;
  on_data_file.GetObject( "stv_tree", on_data_tree );

  std::string on_data_hist_name = hist_name_prefix + "-on";
  TH1D* on_data_hist = new TH1D(on_data_hist_name.c_str(), plot_title.c_str(), Nbins, xmin, xmax);
  on_data_tree->Draw((branch + " >> " + on_data_hist_name).c_str(), selection.c_str());
  on_data_hist->Scale(1.);
  on_data_hist->SetLineColor(kBlack);
  on_data_hist->SetLineWidth(3);
  on_data_hist->SetStats(false);
  on_data_hist->GetXaxis()->SetTitleOffset(1.);
  on_data_hist->GetXaxis()->SetTitleSize(0.05);
  on_data_hist->GetYaxis()->SetTitleSize(0.05);
  on_data_hist->SetDirectory( nullptr );

  // Initialize empty stacked histograms by MC event category
  // TODO: redo this differently (with a std::map, perhaps?)
  std::vector<TH1D*> mc_hists;
  for ( int cat = 1; cat <= 9; ++cat ) {
    std::string temp_mc_hist_name = hist_name_prefix + "-temp_mc" + std::to_string(cat);
    TH1D* temp_mc_hist = new TH1D(temp_mc_hist_name.c_str(), ("; " + var_name + "; events / POT").c_str(),
      Nbins, xmin, xmax);
    mc_hists.push_back( temp_mc_hist );
    temp_mc_hist->SetFillColor( cat == 9 ? 11 : cat + 1 );
    temp_mc_hist->SetStats(false);
    temp_mc_hist->SetDirectory( nullptr );
  }

  // Loop over the different MC files and collect their contributions.
  // We have to handle them separately in order to get the POT normalization correct.
  for ( const auto& mc_file_name : mc_file_names ) {
    // Get the POT values from the current input MC file
    TFile temp_mc_file( mc_file_name.c_str(), "read" );
    TParameter<float>* temp_pot = nullptr;
    temp_mc_file.GetObject( "summed_pot", temp_pot );

    double mc_pot = temp_pot->GetVal();

    // Use a TChain to analyze the MC events
    TChain mc_ch( "stv_tree" );
    mc_ch.Add( mc_file_name.c_str() );

    // Add this file's contribution to the stacked histograms by MC event
    // category
    for ( int cat = 1; cat <= 9; ++cat ) {
      std::string temp_mc_hist_name = hist_name_prefix + "-temp_mc" + std::to_string(cat) + mc_file_name;
      TH1D* temp_mc_hist = new TH1D(temp_mc_hist_name.c_str(), ("; " + var_name + "; events / POT").c_str(),
        Nbins, xmin, xmax);
      mc_ch.Draw( (branch + " >> " + temp_mc_hist_name).c_str(), (selection + " && category == " + std::to_string(cat)).c_str() );

      // Scale to the same exposure as the beam on data
      temp_mc_hist->Scale( POT_ON_DATA / mc_pot );

      // Add this histogram's contribution (now properly scaled) to the total
      mc_hists.at(cat - 1)->Add( temp_mc_hist );
    }
  } // loop over MC files

   TPad* pad1 = new TPad( "pad1", "", 0.028, 0.3, 1.0, 1.0 );
   pad1->SetBottomMargin(0);
//   pad1->SetRightMargin(0.06);
//   pad1->SetLeftMargin(0.13);
//   pad1->SetGridx();
   pad1->Draw();
   pad1->cd();

  on_data_hist->Draw();

  THStack* stacked_hist = new THStack("mc", "");
  TH1D* stacked_histo = new TH1D( "stacked_histo", "", Nbins, xmin, xmax );

  stacked_hist->Add( off_data_hist );
  stacked_histo->Add( off_data_hist );

  for ( const auto& hist : mc_hists ) {
    stacked_hist->Add( hist );
    stacked_histo->Add( hist );
  }

   

  stacked_hist->Draw("hist same");

  on_data_hist->Draw("same");


  TLegend* lg = new TLegend(0.7, 0.5, 0.89, 0.88);
  lg->AddEntry(on_data_hist, "on data", "l");
  lg->AddEntry(off_data_hist, "EXT", "f");
  for ( int cat = 1; cat <= 9; ++cat ) {
    lg->AddEntry( mc_hists.at(cat - 1), cat_to_label(cat).c_str(), "f" );
  }
  lg->Draw("same");

  // Ratio plot
   c1->cd();

   TPad* pad2 = new TPad( "pad2", "", 0, 0.02, 0.96, 0.28);
   pad2->SetTopMargin(0);
   pad2->SetFrameFillStyle(4000);
   pad2->SetBottomMargin(0.45);
   pad2->SetRightMargin(0.06);
   pad2->SetLeftMargin(0.13);

   pad2->SetGridx();
   pad2->Draw();
   pad2->cd();

   // Ratio plot
   TH1D *h_ratio = (TH1D*)on_data_hist->Clone( "h_ratio" );
   h_ratio->SetStats(0);

   h_ratio->Divide( stacked_histo );
   h_ratio->SetLineWidth(1);
   h_ratio->SetLineColor(kBlack);
   h_ratio->SetMarkerStyle(kFullCircle);
   h_ratio->SetMarkerSize(0.6);

   h_ratio->SetTitle("");

   h_ratio->GetYaxis()->SetTitle("Ratio");
//   h_ratio->GetXaxis()->SetTitle( plot_title.c_str() );

//  h_ratio->GetXaxis()->CenterTitle(true);
//   h_ratio->GetXaxis()->SetLabelFont(42);
   h_ratio->GetXaxis()->SetLabelSize(0.12);
   h_ratio->GetXaxis()->SetTitleSize(0.18);
   h_ratio->GetXaxis()->SetTickLength(0.05);
   h_ratio->GetXaxis()->SetTitleOffset(0.9);
//   h_ratio->GetXaxis()->SetTitleFont(42);

   h_ratio->GetYaxis()->CenterTitle(true);
//   h_ratio->GetYaxis()->SetLabelFont(42);
   h_ratio->GetYaxis()->SetLabelSize(0.08);
   h_ratio->GetYaxis()->SetTitleSize(0.16);
   h_ratio->GetYaxis()->SetTitleOffset(0.27);
   h_ratio->GetYaxis()->SetRangeUser(0.6, 1.4);
//   h_ratio->GetYaxis()->SetTitleFont(42);

   h_ratio->Draw("E1");

    TLine *line = new TLine(h_ratio->GetXaxis()->GetXmin(), 1, h_ratio->GetXaxis()->GetXmax(), 1);
    line->SetLineColor(kBlack);
    line->SetLineStyle(9); // dashed
    line->Draw();

   c1->Update();

    c1->SaveAs( ( hist_name_prefix + ".jpg").c_str() );

}

void plots2() {
  //make_plots("delta_alphaT", "reco #delta#alpha_{T} (GeV)", 0., M_PI, 10);
  //make_plots("delta_phiT", "reco #delta#phi_{T} (GeV)", 0., M_PI, 10);

  std::vector<std::string> mc_file_names = { "dirt_stv.root",
    "ncdelta_stv.root", "nue_overlay_stv.root",
    "numu_overlay_stv.root" };

  //make_plots( "pn", "reco p_{n} (GeV)", 0., 2., 10, mc_file_names );
 // make_plots("delta_pT", "reco #deltap_{T} (GeV)", 0., 2., 10, mc_file_names);
//  make_plots("delta_phiT", "reco #delta_{#alphaT}", 0., M_PI, 10, mc_file_names);
//  make_plots("delta_phiT", "reco #delta_{#alphaT}", 0., M_PI, 10, mc_file_names);
  make_plots("pmu", "p3_mu.Mag()", "reco p_{#mu} (GeV)", 0., 2., 10, mc_file_names);
//  make_plots("pp", "p3_lead_p.Mag()", "reco p_{lead p} (GeV)", 0., 2., 10, mc_file_names);
}
