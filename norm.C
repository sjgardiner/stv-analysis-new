// ROOT includes
#include "TFile.h"

// STV analysis includes
#include "SystematicsCalculator.hh"

void norm() {

  SystematicsCalculator sc(
    "/uboone/data/users/gardiner/respmat-myconfig_delta_pTx.root",
    "mcc9_delta_pTx"
  );

  // Write the final histograms to the output TFile
  TFile out_tfile( "myout.root", "recreate" );

  sc.save_universes( out_tfile );

}
