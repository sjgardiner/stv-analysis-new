#include "ResponseMatrixMaker.hh"

ResponseMatrixMaker* resp_mat = nullptr;

void test_rmm() {

  resp_mat = new ResponseMatrixMaker( "myconfig.txt" );

  resp_mat->add_input_file( "/uboone/data/users/gardiner/ntuples-stv/stv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root" );

  //resp_mat->add_input_file( "/uboone/data/users/gardiner/ntuples-stv/stv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run2_reco2_D1D2_reco2.root" );

  std::cout << "POT = " << resp_mat->pot() << '\n';

  //const auto& true_bins = resp_mat->true_bins();
  //const auto& reco_bins = resp_mat->reco_bins();

  resp_mat->build_entry_lists();
  resp_mat->build_response_matrices();
}
