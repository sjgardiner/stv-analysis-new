#include "analyzer.C"

void reproc() {

  // Process the on beam data (use a dummy GENIE file)
  analyzer("/uboone/app/users/gardiner/myroot/v20/run1/nucc_on_data_run1_mcc9.root",
    "on_data_stv.root", "/uboone/app/users/gardiner/myroot2/sorted_ghep_files/nucc_nu_overlay_run1_big_mcc9_genie.root");

  // Process the off beam data (use a dummy GENIE file)
  analyzer("/uboone/app/users/gardiner/myroot/v20/run1/nucc_off_data_run1_mcc9.root",
    "off_data_stv.root", "/uboone/app/users/gardiner/myroot2/sorted_ghep_files/nucc_nu_overlay_run1_big_mcc9_genie.root");

  // Process the numu overlay MC
  analyzer("/uboone/app/users/gardiner/myroot/v20/run1/nucc_nu_overlay_run1_big_mcc9.root",
    "numu_overlay_stv.root", "/uboone/app/users/gardiner/myroot2/sorted_ghep_files/nucc_nu_overlay_run1_big_mcc9_genie.root");

  // Process the nue overlay MC
  analyzer("/uboone/app/users/gardiner/myroot/v20/run1/nucc_nue_overlay_run1_mcc9.root",
    "nue_overlay_stv.root", "/uboone/app/users/gardiner/myroot2/sorted_ghep_files/nucc_nue_overlay_run1_mcc9_genie.root");

  // Process the dirt MC
  analyzer("/uboone/app/users/gardiner/myroot/v20/run1/nucc_dirt_overlay_run1_mcc9.root",
    "dirt_stv.root", "/uboone/app/users/gardiner/myroot2/sorted_ghep_files/nucc_dirt_overlay_run1_mcc9_genie.root");

  // Process the NC Delta MC
  analyzer("/uboone/app/users/gardiner/myroot/v20/run1/nucc_ncdelta_overlay_run1_mcc9.root",
    "ncdelta_stv.root", "/uboone/app/users/gardiner/myroot2/sorted_ghep_files/nucc_ncdelta_overlay_run1_mcc9_genie.root");
}
