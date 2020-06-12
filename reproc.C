#include "analyzer.C"

void reproc() {

  // Process the on beam data
  analyzer("/pnfs/uboone/persistent/users/gardiner/stv-ntuples-new/bnb_5e19_run1.root",
    "on_data_stv.root");

  // Process the off beam data
  // TODO: add extC2 for more statistics
  analyzer("/pnfs/uboone/persistent/users/gardiner/stv-ntuples-new/extC1_run1.root",
    "off_data_stv.root");

  // Process the numu overlay MC
  analyzer("/pnfs/uboone/persistent/users/gardiner/stv-ntuples-new/numu_run1.root",
    "numu_overlay_stv.root");

  // Process the nue overlay MC
  analyzer("/pnfs/uboone/persistent/users/gardiner/stv-ntuples-new/nue_run1.root",
    "nue_overlay_stv.root");

  // Process the dirt MC
  analyzer("/pnfs/uboone/persistent/users/gardiner/stv-ntuples-new/dirt_run1.root",
    "dirt_stv.root");

}
