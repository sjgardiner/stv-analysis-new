#include "ResponseMatrix.hh"

// Define the missing intersection operator for TEntryList objects. Work
// with smart pointers because we're using them elsewhere.
// See tinyurl.com/79q6ztbt for details.
std::unique_ptr< TEntryList > operator*( const std::unique_ptr<TEntryList>& t1,
  const std::unique_ptr<TEntryList>& t2 )
{
  auto intersection = std::make_unique< TEntryList >( *t1 );
  auto temp_copy_t1 = std::make_unique< TEntryList >( *t1 );
  temp_copy_t1->Subtract( t2.get() );
  intersection->Subtract( temp_copy_t1.get() );

  std::string merged_name = t1->GetName();
  merged_name += "_and_";
  merged_name += t2->GetName();

  intersection->SetName( merged_name.c_str() );
  return intersection;
}

// By default, weight the MC events using the MicroBooNE CV tune
const std::string DEFAULT_MC_EVENT_WEIGHT = "spline_weight * (std::isfinite("
  "tuned_cv_weight) && tuned_cv_weight <= 100. ? tuned_cv_weight : 1)";

ResponseMatrix* resp_mat = nullptr;

void test_response_matrix2() {

  resp_mat = new ResponseMatrix( "myconfig.txt" );
  resp_mat->set_overall_mc_weight( DEFAULT_MC_EVENT_WEIGHT );

  TFile* in_tfile = new TFile( "/uboone/data/users/gardiner/ntuples-stv/stv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root", "read" );

  TTree* stv_tree;
  in_tfile->GetObject( "stv_tree", stv_tree );

  TParameter<float>* summed_pot;
  in_tfile->GetObject( "summed_pot", summed_pot );

  //resp_mat->fill( *stv_tree, summed_pot->GetVal() );

  const auto& true_bins = resp_mat->true_bins();
  const auto& reco_bins = resp_mat->reco_bins();

  // Ensures that every TEntryList that is automatically generated
  // has a unique ROOT name
  static int dummy_tb_counter = 0;
  static int dummy_rb_counter = 0;

  std::vector< std::unique_ptr<TEntryList> > true_entry_lists;
  std::vector< std::unique_ptr<TEntryList> > reco_entry_lists;

  for ( const auto& tb : true_bins ) {

    std::string entry_list_name = "tb_entry_list"
      + std::to_string( dummy_tb_counter );

    std::string temp_var_expr = ">> " + entry_list_name;

    stv_tree->Draw( temp_var_expr.c_str(),
      tb.signal_cuts_.c_str(), "entrylist" );

    auto* temp_el = dynamic_cast<TEntryList*>(
      gDirectory->Get(entry_list_name.c_str()) );

    true_entry_lists.emplace_back( temp_el );

    ++dummy_tb_counter;
  }

  for ( const auto& rb : reco_bins ) {

    std::string entry_list_name = "rb_entry_list"
      + std::to_string( dummy_rb_counter );

    std::string temp_var_expr = ">> " + entry_list_name;

    stv_tree->Draw( temp_var_expr.c_str(),
      rb.selection_cuts_.c_str(), "entrylist" );

    auto* temp_el = dynamic_cast<TEntryList*>(
      gDirectory->Get(entry_list_name.c_str()) );

    reco_entry_lists.emplace_back( temp_el );

    ++dummy_rb_counter;
  }

  for ( const auto& tel : true_entry_lists ) {
    std::cout << tel->GetName() << '\n';
    tel->Print();
    std::cout << '\n';
  }

  for ( const auto& rel : reco_entry_lists ) {
    std::cout << rel->GetName() << '\n';
    rel->Print();
    std::cout << '\n';
  }

  for ( const auto& tel : true_entry_lists ) {
    for ( const auto& rel : reco_entry_lists ) {
      // Compute the intersection of the entry lists
      auto inter_tr_el = tel * rel;
      std::cout << '\"' << inter_tr_el->GetName() << "\"\n";
      inter_tr_el->Print();
      std::cout << '\n';
    }
  }

}
