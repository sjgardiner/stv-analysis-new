#pragma once

// Standard library includes
#include <array>
#include <map>
#include <string>

// ROOT includes
#include "TTree.h"

// STV analysis includes
#include "TreeUtils.hh"

// Class that provides temporary storage for event weights being processed by a
// ResponseMatrixMaker object
class WeightHandler {
  public:

    WeightHandler() {};

    void set_branch_addresses( TTree& in_tree );

    // Access the owned map
    inline const auto& weight_map() const { return weight_map_; }
    inline auto& weight_map() { return weight_map_; }

  protected:

    // Keys are branch names in the input TTree, values point to vectors of
    // event weights
    std::map< std::string, MyPointer< std::vector<double> > > weight_map_;
};

void WeightHandler::set_branch_addresses( TTree& in_tree ) {

  // Delete any pre-existing contents of the weight map
  weight_map_.clear();

  // Loop over each of the branches of the input TTree
  auto* lob = in_tree.GetListOfBranches();
  for ( int b = 0; b < lob->GetEntries(); ++b ) {

    auto* branch = dynamic_cast< TBranch* >( lob->At(b) );

    // Add a new vector in the weight map for every branch whose name begins
    // with "weight_". Assume that all such branches store a
    // std::vector<double> object. Set the branch address so that the
    // vector can accept input from the TTree.
    const std::string wgt_br_prefix = "weight_";
    std::string br_name = branch->GetName();
    if ( br_name.compare(0, wgt_br_prefix.size(), wgt_br_prefix) == 0 ) {

      weight_map_[ br_name ] = MyPointer< std::vector<double> >();

      auto& wgt_vec = weight_map_.at( br_name );
      set_object_input_branch_address( in_tree, br_name, wgt_vec );

    }
  }
}
