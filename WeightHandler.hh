#pragma once

// Standard library includes
#include <array>
#include <map>
#include <string>

// ROOT includes
#include "TTree.h"

// STV analysis includes
#include "TreeUtils.hh"

constexpr std::array< const char*, 15 > UNIVERSE_WEIGHT_NAMES = {
  "weight_All_UBGenie",
  "weight_AxFFCCQEshape_UBGenie",
  "weight_DecayAngMEC_UBGenie",
  "weight_NormCCCOH_UBGenie",
  "weight_NormNCCOH_UBGenie",
  "weight_RPA_CCQE_UBGenie",
  "weight_ThetaDelta2NRad_UBGenie",
  "weight_Theta_Delta2Npi_UBGenie",
  "weight_TunedCentralValue_UBGenie",
  "weight_VecFFCCQEshape_UBGenie",
  "weight_XSecShape_CCMEC_UBGenie",
  "weight_flux_all",
  "weight_reint_all",
  "weight_xsr_scc_Fa3_SCC",
  "weight_xsr_scc_Fv3_SCC"
};

constexpr std::array< const char*, 1 > BUG_FIX_WEIGHT_NAMES = {
  "weight_splines_general_Spline"
  // "weight_RootinoFix_UBGenie"
};

// Class that provides temporary storage for event weights being processed by a
// ResponseMatrixMaker object
class WeightHandler {
  public:

    WeightHandler();

    void set_branch_addresses( TTree& in_tree );

    // Systematic variation weights for each universe
    std::map< std::string, MyPointer< std::vector<double> > >
      universe_weight_map_;

    // Bug-fix weights that should be applied unconditionally in every
    // universe (to correct a problem with the central value)
    std::map< std::string, MyPointer< std::vector<double> > >
      bug_fix_weight_map_;

};

WeightHandler::WeightHandler() {
  for ( const std::string& name : UNIVERSE_WEIGHT_NAMES ) {
    universe_weight_map_[ name ] = MyPointer< std::vector<double> >();
  }

  for ( const std::string& name : BUG_FIX_WEIGHT_NAMES ) {
    bug_fix_weight_map_[ name ] = MyPointer< std::vector<double> >();
  }
}

void WeightHandler::set_branch_addresses( TTree& in_tree )
{
  // The maps owned by the WeightHandler object use branch names as the keys
  // and the objects to be filled as the values. We can therefore just iterate
  // over them to set the branch addresses for each entry.
  for ( auto& pair : universe_weight_map_ ) {
    set_object_input_branch_address( in_tree, pair.first, pair.second );
  }

  for ( auto& pair : bug_fix_weight_map_ ) {
    set_object_input_branch_address( in_tree, pair.first, pair.second );
  }
}
