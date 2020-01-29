// Adapted from MigrationMatrices.cpp from Afro

// ROOT integer code for Arial font
constexpr int FONT_STYLE = 62; // Arial

// When bins have zero content, set them to this very small value so that the
// colz style will still paint them
constexpr double REALLY_SMALL = 1e-11;

void make_resolution_plot( const std::string& stv_file_name,
  const std::string& branch, const std::string& variable_name,
  const std::string& hist_name, int num_bins, double x_min, double x_max)
{
  // Assume that the user supplied the name of the reco observable
  // in the branch variable. Prepend "mc_" to it and assume that
  // the corresponding MC truth quantity is stored there.
  // TODO: revisit this as needed
  std::string mc_branch = "mc_" + branch;

  std::string plot_title = "Resolution; MC Truth " + variable_name
    + "^{truth}; Reconstructed " + variable_name + "^{reco}";

  TFile stv_file( stv_file_name.c_str(), "read" );
  TTree* stv_tree = nullptr;
  stv_file.GetObject( "stv_tree", stv_tree );

  // Round all numbers to this precision when rendering them on plots
  gStyle->SetPaintTextFormat( "4.2f" );

  TCanvas* c1 = new TCanvas;
  c1->SetBottomMargin( 0.15 );
  c1->SetLeftMargin( 0.13 );
  c1->cd();

  // Define a 2D histogram to hold the information for the resolution plot. Use
  // the same binning for reco and truth variables (the y and x axes,
  // respectively)
  TH2D* hist = new TH2D( hist_name.c_str(), plot_title.c_str(), num_bins,
    x_min, x_max, num_bins, x_min, x_max );

  std::string plot_invocation = mc_branch + ':' + branch + ">> " + hist_name;
  stv_tree->Draw( plot_invocation.c_str(), "is_mc && mc_is_signal"
    " && sel_CCNp0pi && genie_ok", "colz" );

  hist->SetDirectory( nullptr );

  // Histogram style options
  hist->GetXaxis()->SetTitleFont( FONT_STYLE);
  hist->GetYaxis()->SetTitleFont( FONT_STYLE );
  hist->GetXaxis()->SetTitleSize( 0.05 );
  hist->GetYaxis()->SetTitleSize( 0.05 );
  hist->GetXaxis()->SetLabelFont( FONT_STYLE );
  hist->GetYaxis()->SetLabelFont( FONT_STYLE );
  hist->GetZaxis()->SetLabelFont( FONT_STYLE );
  hist->GetZaxis()->SetLabelSize( 0.03 );
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
  hist->GetXaxis()->SetTitleOffset( 1.2 );
  hist->GetYaxis()->SetTitleOffset( 1.1 );
  hist->SetStats( false );
  hist->SetMarkerSize( 1.8 ); // text size
  hist->SetMarkerColor( kWhite ); // text color

  // Normalize bins
  int num_bins_x = hist->GetXaxis()->GetNbins();
  int num_bins_y = hist->GetYaxis()->GetNbins();

  // Loop over the x bins. Include the underflow (index zero) and overflow
  // (index num_bins_x + 1) bins.
  for ( int bx = 0; bx <= num_bins_x + 1; ++bx ) {

    // For the current x bin, compute the sum of all y bins.
    double y_sum = 0.;
    for ( int by = 0; by <= num_bins_y + 1; ++by ) {
      y_sum += hist->GetBinContent( bx, by );
    }

    // Normalize each of the y bins so that the sum over y is unity.
    for ( int by = 0; by <= num_bins_y + 1; ++by ) {

      // To avoid dividing by zero, set the bins contents and
      // errors to zero if the sum of the y bins is not positive.
      if ( y_sum <= 0. ) {
        hist->SetBinContent( bx, by, REALLY_SMALL );
        hist->SetBinError( bx, by, REALLY_SMALL );
      }
      else {
        // Otherwise, normalize in the usual way
        double bc = hist->GetBinContent( bx, by );
        double berr = hist->GetBinError( bx, by );

        double content = std::max( bc / y_sum, REALLY_SMALL );

        // Apply the correct normalization to the bin error
        double error = std::sqrt( std::pow(berr / y_sum, 2)
          + std::pow(bc * std::sqrt(y_sum) / std::pow(y_sum, 2), 2) );

        // Update the content and error of the current bin appropriately
        hist->SetBinContent( bx, by, content );
        hist->SetBinError( bx, by, error );
      }
    }
  }

  // TODO: revisit this
  // Set the overflow bin in x and y to a negative value. This will force
  // drawing of all bins with zero content in the "colz" style.
  //c1->SetFrameFillColor(TColor::GetColorPalette(0));
  //hist->SetBinContent( num_bins_x, num_bins_y, -1e-5 );

  // Note: adding option "E" to this also displays error in each bin
  hist->Draw("text colz");

}

void resolution() {

  make_resolution_plot( "numu_overlay_stv.root", "p3_mu.Mag()", "p_{#mu}", "pmu", 10, 0., 2. );
  make_resolution_plot( "numu_overlay_stv.root", "p3_mu.CosTheta()", "cos#theta_{#mu}", "cthmu", 10, -1., 1. );
  make_resolution_plot( "numu_overlay_stv.root", "p3_mu.Phi()", "#phi_{#mu}", "phimu", 10, 0., M_PI );

  make_resolution_plot( "numu_overlay_stv.root", "p3_lead_p.Mag()", "p_{lead p}", "pp", 10, 0., 2. );
  make_resolution_plot( "numu_overlay_stv.root", "p3_lead_p.CosTheta()", "cos#theta_{lead p}", "cthp", 10, -1., 1. );
  make_resolution_plot( "numu_overlay_stv.root", "p3_lead_p.Phi()", "#phi_{lead_p}", "phip", 10, 0., M_PI );

  make_resolution_plot( "numu_overlay_stv.root", "delta_pT", "#deltap_{T}", "pT", 10, 0., 2. );
  make_resolution_plot( "numu_overlay_stv.root", "delta_phiT", "#delta#phi_{T}", "phiT", 10, 0., M_PI );
  make_resolution_plot( "numu_overlay_stv.root", "delta_alphaT", "#delta#alpha_{T}", "alphaT", 10, 0., M_PI );
  make_resolution_plot( "numu_overlay_stv.root", "delta_pL", "#deltap_{L}", "pL", 10, 0., 2. );
  make_resolution_plot( "numu_overlay_stv.root", "pn", "p_{n}", "pn", 10, 0., 2. );

}
