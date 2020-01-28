// Adapted from MigrationMatrices.cpp from Afro

std::string variable;
std::string variable_addition;
std::string variable_title;
std::string xmax;
std::string ymax;
int FONTSTYLE = 62; // Arial


void resolution() {

  // Give names to resolution plots according to variables
  vector<TString> PlotNames;
  PlotNames.push_back( "delta_pT_plot2D" );
        PlotNames.push_back( "delta_phiT_plot2D" );
        PlotNames.push_back( "delta_alphaT_plot2D" );
        PlotNames.push_back( "p_mu_truth_vs_mcs_plot2D" );
  PlotNames.push_back( "p_p_truth_vs_mcs_plot2D" );


  // Name of samples for each plot
  vector<TString> NameOfSamples;;
  NameOfSamples.push_back( "SampleName" );

  const int num_N2DPlots = PlotNames.size();
  const int num_NSamples = NameOfSamples.size();

  cout << "Number of 2D-plots: " << num_N2DPlots << endl;
        cout << "Number of samples: " << num_NSamples << endl;

  TCanvas *PlotCanvas[1][5] = {};
  TFile * FileSample[1] = {};
  TH1D* h_1DPlots[1][5] = {};

  std::string output_file_name = "plots_resolution.root";
  TFile *output_file = new TFile( output_file_name.c_str(), "recreate" );

  // Open files that contain samples
  for ( int which_sample = 0; which_sample < num_NSamples; which_sample++ ) {
    FileSample[which_sample] = TFile::Open( "numu_overlay_stv.root" );
  }

  for ( int which_sample = 0; which_sample < num_NSamples; which_sample++ ) {

    for ( int which_plot = 0; which_plot < num_N2DPlots; which_plot++ ) {


      if ( which_plot == 0 ) {
        variable = "delta_pT";
                                variable_addition = "";
        variable_title = "#delta_{p_{T}}";
        xmax = "2";
        ymax = "2";
      }

                     if ( which_plot == 1 ) {
                                variable = "delta_phiT";
                                variable_addition = "";
                                variable_title = "#delta_{#phi_{T}}";
        xmax = "3.141592654";
                                ymax = xmax;
                        }

                     if ( which_plot == 2 ) {
                                variable = "delta_alphaT";
                                variable_addition = "";
                                variable_title = "#delta_{#alpha_{T}}";
                                xmax = "3.141592654";
                                ymax = xmax;
                        }

                     if ( which_plot == 3 ) {
                                variable = "p3_mu.Mag()";
                                variable_addition = "";
                                variable_title = "p_{#mu_{mcs}}";
                                xmax = "2";
                                ymax = "2";
                        }

                     if ( which_plot == 4 ) {
                                variable = "p3_lead_p.Mag()";
                                variable_addition = "";
                                variable_title = "p_{p_{range}}";
                                xmax = "2";
                                ymax = "2";
                        }

      // Create canvas
      PlotCanvas[which_sample][which_plot] = new TCanvas( PlotNames[which_plot] + NameOfSamples[which_sample], PlotNames[which_plot] + NameOfSamples[which_sample] );

      PlotCanvas[which_sample][which_plot]->SetBottomMargin(0.15);
                        PlotCanvas[which_sample][which_plot]->SetLeftMargin(0.13);

      PlotCanvas[which_sample][which_plot]->cd(); // change directory to that canvas

      // Create TTree
                  TTree *stv_tree = (TTree*)FileSample[which_sample]->Get("stv_tree");
                  gROOT->cd(); // change directory to that TTree

      // Decide what to put into TH1D
                  stv_tree->Draw( ( std::string("mc_") + variable.c_str() + std::string(":") + variable.c_str() + variable_addition.c_str() + std::string(">>h(10,0,") + xmax.c_str() + std::string(",10,0,") + ymax.c_str() + std::string(")") ).c_str(), "is_mc && mc_is_signal && sel_CCNp0pi && genie_ok", "colz" );

                        gStyle->SetPaintTextFormat("4.2f"); // refers to "text" drawing option (numbers are rounded)


      // Save it in two dimensional histogram array (matrix)
      h_1DPlots[which_sample][which_plot] = (TH1D*)(gDirectory->Get("h") );

      // Titles
                        h_1DPlots[which_sample][which_plot]->SetTitle( ( std::string("Resolution;") + std::string("MC Truth ") + variable_title.c_str() + std::string("^{truth};") + std::string("Reconstructed ") + variable_title.c_str() + std::string("^{reco};") ).c_str() );
//                        h_1DPlots[which_sample][which_plot]->GetXaxis()->SetTitle( ( std::string("mc_") + variable.c_str() ).c_str() );


      // Rebin histogram
//      h_1DPlots[which_sample][which_plot]->Rebin(20, "d");


      // Aesthetics
      h_1DPlots[which_sample][which_plot]->GetXaxis()->SetTitleFont( FONTSTYLE);
      h_1DPlots[which_sample][which_plot]->GetYaxis()->SetTitleFont( FONTSTYLE );
      h_1DPlots[which_sample][which_plot]->GetXaxis()->SetTitleSize( 0.05 );
      h_1DPlots[which_sample][which_plot]->GetYaxis()->SetTitleSize( 0.05 );
      h_1DPlots[which_sample][which_plot]->GetXaxis()->SetLabelFont( FONTSTYLE );
      h_1DPlots[which_sample][which_plot]->GetYaxis()->SetLabelFont( FONTSTYLE );
      h_1DPlots[which_sample][which_plot]->GetZaxis()->SetLabelFont( FONTSTYLE );
      h_1DPlots[which_sample][which_plot]->GetZaxis()->SetLabelSize( 0.03 );
                        h_1DPlots[which_sample][which_plot]->GetXaxis()->CenterTitle();
                        h_1DPlots[which_sample][which_plot]->GetYaxis()->CenterTitle();

      gStyle->SetOptStat(0);

      h_1DPlots[which_sample][which_plot]->SetMarkerSize(1.8); // Size of text
                        h_1DPlots[which_sample][which_plot]->SetMarkerColor(kWhite); // Color of text

      // Adjust bins
      int num_bin_x = h_1DPlots[which_sample][which_plot]->GetXaxis()->GetNbins(); // Get bins in x-dimension
//      cout << "Number of x-axis bins: " << num_bin_x << endl;
      int num_bin_y = h_1DPlots[which_sample][which_plot]->GetYaxis()->GetNbins(); // Get bins in y-dimension

      for ( int which_bin_x = 0; which_bin_x < num_bin_x; which_bin_x++ ) {

        int num_events_current_xy_bin = 0;

        for ( int which_bin_y = 0; which_bin_y < num_bin_y; which_bin_y ++ ) {

          num_events_current_xy_bin += h_1DPlots[which_sample][which_plot]->GetBinContent( which_bin_x + 1, which_bin_y + 1 );
        }

        for ( int which_bin_y = 0; which_bin_y < num_bin_y; which_bin_y++ ) {

          if ( num_events_current_xy_bin > 0 ) {

            // Volume normalise
            double xy_normalisation = h_1DPlots[which_sample][which_plot]->GetBinContent( which_bin_x + 1, which_bin_y + 1 )/double( num_events_current_xy_bin );

            // Set bin content to volume-normalised entry
            h_1DPlots[which_sample][which_plot]->SetBinContent( which_bin_x + 1, which_bin_y + 1, xy_normalisation );

            // Calculate error
            double error = std::sqrt( TMath::Power(h_1DPlots[which_sample][which_plot]->GetBinError( which_bin_x + 1, which_bin_y + 1 )/double( num_events_current_xy_bin ),2.) + TMath::Power( h_1DPlots[which_sample][which_plot]->GetBinContent( which_bin_x + 1, which_bin_y + 1 ) * sqrt( num_events_current_xy_bin )/double( num_events_current_xy_bin*num_events_current_xy_bin),2.));
            h_1DPlots[which_sample][which_plot]->SetBinError( which_bin_x + 1, which_bin_y + 1, error );
          }

          else {

            h_1DPlots[which_sample][which_plot]->SetBinContent( which_bin_x + 1, which_bin_y + 1, 0. );
            h_1DPlots[which_sample][which_plot]->SetBinError( which_bin_x + 1, which_bin_y + 1, 0. );

          }

        }

      }


      //h_1DPlots[which_sample][which_plot]->GetZaxis()->SetRangeUser( 0,2. );

      // Write to output file resolution.root
            output_file->cd();
            h_1DPlots[which_sample][which_plot]->Draw("text colz"); // Option E in addition to that also writes error in box
                        PlotCanvas[which_sample][which_plot]->Write();

      //PlotCanvas[which_sample][which_plot]->SaveAs( "Plots/resolution/" + PlotNames[which_plot] + NameOfSamples[which_sample], PlotNames[which_plot] + NameOfSamples[which_sample] + ".pdf" );

      PlotCanvas[which_sample][which_plot]->Close();

    } // End of loop over plots

  } // End of loop over samples

} // End of void

