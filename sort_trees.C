struct EventLabel {

  EventLabel() {}

  EventLabel(unsigned int ev, unsigned int r, unsigned int sr, unsigned int sec,
    unsigned int nsec, int mct_idx) : event( ev ), run( r ), subrun( sr ),
    evt_time_sec( sec ), evt_time_ns( nsec ), mctruth_index( mct_idx ) {}

  // Needed to use this struct as the key type for a std::map
  bool operator<( const EventLabel& other ) const {
    if ( run != other.run ) return run < other.run;
    if ( subrun != other.subrun ) return subrun < other.subrun;
    if ( event != other.event ) return event < other.event;
    if ( evt_time_sec != other.evt_time_sec )
      return evt_time_sec < other.evt_time_sec;
    if ( evt_time_ns != other.evt_time_ns )
      return evt_time_ns < other.evt_time_ns;
    if ( mctruth_index != other.mctruth_index )
      return mctruth_index < other.mctruth_index;
    return false;
  }

  bool operator==( const EventLabel& other ) const {
    if ( run != other.run ) return false;
    if ( subrun != other.subrun ) return false;
    if ( event != other.event ) return false;
    if ( evt_time_sec != other.evt_time_sec ) return false;
    if ( evt_time_ns != other.evt_time_ns ) return false;
    if ( mctruth_index != other.mctruth_index ) return false;
    return true;
  }

  unsigned int event = 0u;
  unsigned int run = 0u;
  unsigned int subrun = 0u;
  unsigned int evt_time_sec = 0u;
  unsigned int evt_time_ns = 0u;
  int mctruth_index = 0;
};

std::ostream& operator<<( std::ostream& os, const EventLabel& el ) {
  os << el.run << ' ' << el.subrun << ' ' << el.event << ' '
    << el.evt_time_sec << ' ' << el.evt_time_ns << ' ' << el.mctruth_index;
  return os;
}

void sort_trees( /* const std::string& ghep_file_name, */
  const std::string& nucc_analyzer_file_name,
  const std::string& output_file_name)
{
  // Use all saved GENIE MC events from the SAM definitions employed by Wouter
  TChain* gtree = new TChain( "gtree" );
  gtree->Add("/uboone/app/users/gardiner/myroot2/merged_files/*");

  TChain* ubtree = new TChain( "ub_tune_cv" );
  ubtree->Add("/uboone/app/users/gardiner/myroot2/merged_files/*");

  //TFile* ghep_file = new TFile( ghep_file_name.c_str(), "read" );
  //TTree* gtree = nullptr;
  //ghep_file->GetObject( "gtree", gtree );
  //if ( !gtree ) return;

  // Allow access to the MicroBooNE CV tune event weights using a friend TTree
  //gtree->AddFriend( "ub_tune_cv" );
  gtree->AddFriend( ubtree );

  EventLabel event_label;

  // Disable retrieval (for now) of the branches that we don't need for
  // indexing
  gtree->SetBranchStatus( "gmcrec", false );
  gtree->SetBranchStatus( "orig_filename", false );
  gtree->SetBranchStatus( "orig_evtnum", false );

  gtree->SetBranchAddress( "art_event", &event_label.event );
  gtree->SetBranchAddress( "art_run", &event_label.run );
  gtree->SetBranchAddress( "art_subrun", &event_label.subrun );
  gtree->SetBranchAddress( "evt_time_sec", &event_label.evt_time_sec );
  gtree->SetBranchAddress( "evt_time_nsec", &event_label.evt_time_ns );
  gtree->SetBranchAddress( "mctruth_index", &event_label.mctruth_index );

  // For speed, build a map-based index of each of the GHEP event labels.
  // Keys are event labels, values are entry indices in the GHEP TTree
  std::map< EventLabel, long > ghep_entry_map;
  for ( long e = 0; e < gtree->GetEntries(); ++e ) {
    if ( e % 10000 == 0 ) std::cout << "Indexing GHEP entry #" << e << '\n';
    gtree->GetEntry( e );

    ghep_entry_map[ event_label ] = e;
  }

  // We're now ready to match the GHEP events to the
  // ones processed by the NuCCAnalyzer
  TFile* nucc_file = new TFile( nucc_analyzer_file_name.c_str(), "read" );
  TTree* nucc_tree = nullptr;
  nucc_file->GetObject( "NuCCanalyzer/Event", nucc_tree );
  if ( !nucc_tree ) return;

  EventLabel ev_label;
  nucc_tree->SetBranchAddress( "event", &ev_label.event );
  nucc_tree->SetBranchAddress( "run", &ev_label.run );
  nucc_tree->SetBranchAddress( "subrun", &ev_label.subrun );
  nucc_tree->SetBranchAddress( "evt_time_sec", &ev_label.evt_time_sec );
  nucc_tree->SetBranchAddress( "evt_time_nsec", &ev_label.evt_time_ns );

  // TODO: revisit this. I believe that Wouter's analyzer always assumes that
  // the first MCTruth object in the event is the only one present.
  ev_label.mctruth_index = 0;

  // Before continuing, enable retrieval of the GENIE event record
  // from the GHEP tree
  genie::NtpMCEventRecord* gmcrec = nullptr;
  gtree->SetBranchStatus( "gmcrec", true );
  gtree->SetBranchAddress( "gmcrec", &gmcrec );

  // Also get the MicroBooNE CV tune event weights using the
  // friend TTree
  float cv_weight;
  gtree->SetBranchAddress( "cv_weight", &cv_weight );

  // Now prepare an output TFile and TTree to store the re-ordered
  // event records and weights

  // Flag that indicates whether the GENIE event was successfully matched to an
  // event from Wouter's Event TTree
  bool genie_ok = false;
  genie::EventRecord* evrec = nullptr;
  TFile* out_file = new TFile( output_file_name.c_str(), "recreate" );
  TTree* out_tree = new TTree( "stv_gtree", "GENIE truth info for STV analysis");
  out_tree->Branch( "evrec", &evrec );
  out_tree->Branch( "cv_weight", &cv_weight );
  out_tree->Branch( "genie_ok", &genie_ok );

  std::map<long, EventLabel> bad_nucc_entries;
  for ( long e = 0; e < nucc_tree->GetEntries(); ++e ) {
    nucc_tree->GetEntry( e );

    bool found_match = static_cast<bool>( ghep_entry_map.count(ev_label) );
    if ( found_match ) {
      // Find the entry of the matching GENIE event record
      long ghep_entry = ghep_entry_map.at( ev_label );
      std::cout << "Matched nucc entry #" << e << " to GHEP #" <<
        ghep_entry << '\n';

      // Retrieve it and set the "ok" flag
      gtree->GetEntry( ghep_entry );
      evrec = gmcrec->event;
      genie_ok = true;
    }
    else {
      bad_nucc_entries[ e ] = ev_label;
      std::cout << "Bad NuCCAnalyzer match!";

      evrec = nullptr;
      cv_weight = 0.; // TODO: Revisit. Should this be 1 instead?
      genie_ok = false;
    }

    // Store the results in the output TTree
    out_tree->Fill();

    // Clean up by freeing memory as needed
    if ( evrec ) delete evrec;
  }

  out_tree->Write();

  std::cout << "\n\nBad entries found: " << bad_nucc_entries.size() << '\n';
  for ( const auto& pair : bad_nucc_entries ) {
    std::cout << "  Entry " << pair.first << " had " << pair.second << '\n';
  }

}
