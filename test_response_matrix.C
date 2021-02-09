#include "ResponseMatrix.hh"

// By default, weight the MC events using the MicroBooNE CV tune
const std::string DEFAULT_MC_EVENT_WEIGHT = "spline_weight * (std::isfinite("
  "tuned_cv_weight) && tuned_cv_weight <= 100. ? tuned_cv_weight : 1)";

ResponseMatrix* resp_mat = nullptr;

void test_response_matrix() {

  resp_mat = new ResponseMatrix( "myconfig.txt" );
  resp_mat->set_overall_mc_weight( DEFAULT_MC_EVENT_WEIGHT );

  TFile* in_tfile = new TFile( "/uboone/data/users/gardiner/ntuples-stv/stv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root", "read" );

  TTree* stv_tree;
  in_tfile->GetObject( "stv_tree", stv_tree );

  TParameter<float>* summed_pot;
  in_tfile->GetObject( "summed_pot", summed_pot );

  resp_mat->fill( *stv_tree, summed_pot->GetVal() );

}
