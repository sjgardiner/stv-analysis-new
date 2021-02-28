// KEY: TParameter<float> summed_pot;1
// KEY: TTree    stv_tree;2      STV analysis tree

constexpr int NUM_BINS = 15;
constexpr double X_MIN = 0.;
constexpr double X_MAX = 0.8;

struct HistInfo {
  HistInfo() {}

  double pot_ = 0.;
  std::string file_name_;
  TH1D* hist_ = nullptr;
};

void detSyst() {

  // Keys are variation names, values are HistInfo objects containing
  // histograms and metadata built using the STV analysis TTree from the
  // corresponding file
  std::map< std::string, HistInfo > var_map;

  std::ifstream list_file( "detSystList.txt" );
  std::string in_file_name;
  while ( std::getline(list_file, in_file_name) ) {

    // Skip lines beginning with '#'
    if ( in_file_name.front() == '#' ) continue;

    std::string variation_name;

    // Get the name of the detector variation sample using the file name
    std::string first_delimiter( "DetVar_" );
    size_t first = in_file_name.find( first_delimiter );

    if ( first == std::string::npos ) {
      variation_name = "overlay";
    }
    else {
      size_t last = std::min( in_file_name.find("_v08"),
        in_file_name.find("_reco2") );
      size_t end_of_first = first + first_delimiter.length();

      variation_name = in_file_name.substr( end_of_first,
        last - end_of_first );
    }

    TFile in_file( in_file_name.c_str(), "read" );

    TParameter<float>* pot = nullptr;
    TTree* in_tree = nullptr;

    in_file.GetObject( "summed_pot", pot );
    in_file.GetObject( "stv_tree", in_tree );

    std::string temp_hist_name = variation_name + "_h";
    TH1D* temp_hist = new TH1D( temp_hist_name.c_str(),
      "", NUM_BINS, X_MIN, X_MAX );

    //in_tree->Draw( ("mc_delta_pT >> " + temp_hist_name).c_str(),
    //  "mc_is_signal * spline_weight", "goff" );

    in_tree->Draw( ("delta_pT >> " + temp_hist_name).c_str(),
      "sel_CCNp0pi * spline_weight", "goff" );

    // Prevent auto-deletion of the histogram when the TFile object it
    // initially belongs to goes out of scope
    temp_hist->SetDirectory( nullptr );

    HistInfo temp_info;
    temp_info.pot_ = pot->GetVal();
    temp_info.hist_ = temp_hist;
    temp_info.file_name_ = in_file_name;

    var_map[ variation_name ] = temp_info;
  }

  // Scale all histograms to the same POT as the detVar CV sample
  int color = 1;
  double CV_pot = var_map.at( "CV" ).pot_;
  for ( auto& pair : var_map ) {
    auto& info = pair.second;
    double pot = info.pot_;
    info.hist_->Scale( CV_pot / pot );
    info.hist_->SetLineWidth( 3 );
    info.hist_->SetLineColor( color );
    info.hist_->SetStats( false );
    if ( color < 9 ) ++color;
    else if ( color == 9 ) color += 11;
    else color += 10;
  }

  auto* cv_hist = var_map.at( "CV" ).hist_;
  cv_hist->SetTitle( "; reco #deltap_{T} [GeV]; detVar CV events (run 3)" );

  // Draw all of the histograms on the same plot
  TLegend* lg = new TLegend( 0.7, 0.7, 0.9, 0.9 );
  cv_hist->Draw( "hist e" );
  cv_hist->GetYaxis()->SetRangeUser( 0., 750. );
  cv_hist->GetYaxis()->SetTitleOffset( 1.05 );
  cv_hist->GetYaxis()->SetTitleSize( 0.045 );
  cv_hist->GetXaxis()->SetTitleOffset( 1.0 );
  cv_hist->GetXaxis()->SetTitleSize( 0.045 );

  for ( auto& pair : var_map ) {
    auto* hist = pair.second.hist_;
    hist->Draw( "same e" );
    lg->AddEntry( hist, pair.first.c_str(), "le" );
  }

  lg->Draw( "same" );

  //// Add fractional uncertainties in quadrature for each bin
  //TH1D* frac_sigmas = new TH1D( "frac_sigmas", "frac_sigmas",
  //  NUM_BINS, X_MIN, X_MAX );

  //auto& cv_hist_info = var_map.at( "CV" );
  //for ( const auto& pair : var_map ) {
  //  std::string label = pair.first;
  //  auto& info = pair.second;

  //  for ( int b = 1; b <= NUM_BINS; ++b ) {
  //    double cv = cv_hist_info.hist_->GetBinContent( b );
  //    double var = info.hist_->GetBinContent( b );
  //    double sigma2 = std::pow( var - cv, 2 );

  //    frac_sigmas->Fill( frac_sigmas->GetBinLowEdge(b), sigma2 );
  //  }
  //}

  //for ( int b = 1; b <= NUM_BINS; ++b ) {
  //  double cv = cv_hist_info.hist_->GetBinContent( b );
  //  double sum_sigma2 = frac_sigmas->GetBinContent( b );
  //  std::cout << sum_sigma2 << '\n';

  //  double frac_err = std::sqrt( sum_sigma2 ) / cv;
  //  frac_sigmas->SetBinContent( b, frac_err );
  //}

  //TCanvas* cfrac = new TCanvas;
  //frac_sigmas->Draw( "hist" );

}
