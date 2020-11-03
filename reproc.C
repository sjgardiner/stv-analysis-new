#include "analyzer.C"

void reproc() {

  // Process the on beam data
  analyzer("/uboone/data/users/davidc/searchingfornues/v08_00_00_43/0702/run1/data_bnb_mcc9.1_v08_00_00_25_reco2_C1_beam_good_reco2_5e19.root",
    "on_data_stv.root");

  // Process the off beam data (ext C1 + C2)
  analyzer("/uboone/data/users/davidc/searchingfornues/v08_00_00_43/0702/run1/data_extbnb_mcc9.1_v08_00_00_25_reco2_C_all_reco2.root",
    "off_data_stv.root");

  // Process the numu overlay MC
  analyzer("/pnfs/uboone/persistent/users/davidc/searchingfornues/v08_00_00_43/0928/prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root",
    "numu_overlay_stv.root");

  // Process the nue overlay MC
  analyzer("/pnfs/uboone/persistent/users/davidc/searchingfornues/v08_00_00_43/0928/prodgenie_bnb_intrinsice_nue_uboone_overlay_mcc9.1_v08_00_00_26_run1_reco2_reco2.root",
    "nue_overlay_stv.root");

  // Process the dirt MC
  analyzer("/pnfs/uboone/persistent/users/davidc/searchingfornues/v08_00_00_43/0928/prodgenie_bnb_dirt_overlay_mcc9.1_v08_00_00_26_run1_reco2_reco2.root",
    "dirt_stv.root");

}
