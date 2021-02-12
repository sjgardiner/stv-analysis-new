#pragma once

// Standard library includes
#include <cmath>

// Needed for the MyPointer class template
#include "TreeUtils.hh"

constexpr double MIN_WEIGHT = 0.;
constexpr double MAX_WEIGHT = 30.;

// Event weights that are below MIN_WEIGHT, above MAX_WEIGHT, infinite, or NaN
// are reset to unity by this function. Other weights are returned unaltered.
inline double safe_weight( double w ) {
  if ( std::isfinite(w) && w >= MIN_WEIGHT && w <= MAX_WEIGHT ) return w;
  else return 1.0;
}

class MySelector : public TSelector {

  public:

    MySelector( const std::string& br_name ) : TSelector(),
      branch_name_( br_name ) {}

    // TSelector interface
    virtual void Begin( TTree* );
    virtual void Init( TTree* tree );
    virtual Bool_t Notify();
    virtual void SlaveBegin( TTree* tree );
    virtual void SlaveTerminate();
    virtual Bool_t Process( Long64_t entry );
    virtual void Terminate();
    inline virtual int Version() const { return 2; }

  protected:

    // Input TTree or TChain (not owned by this class)
    TTree* chain_ = nullptr;

    // Branches of interest for the analysis
    TBranch* br_weight_vec_ = nullptr;

    // Temporary storage to use for reading from the TChain
    MyPointer< std::vector<double> > weight_vec_;

    // Vector that holds summed event weights
    std::vector< double > sums_of_weights_;

    // Vector that holds summed squares of event weights
    std::vector< double > sums_of_weights2_;

    // Name of the branch to be summed
    std::string branch_name_;

};

void MySelector::Init( TTree* tree ) {

  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses of the tree
  // will be set. It is normaly not necessary to make changes to the
  // generated code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running with PROOF.

  if ( !tree ) return;

  chain_ = tree;

  set_object_input_branch_address( *chain_, branch_name_, weight_vec_ );

}

Bool_t MySelector::Notify() {

  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. Typically here the branch pointers
  // will be retrieved. It is normaly not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed.

  if ( chain_ ) {
    // Get a pointer to the branch of interest
    br_weight_vec_ = chain_->GetBranch( branch_name_.c_str() );
  }

  return kTRUE;
}

Bool_t MySelector::Process( Long64_t entry ) {
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either TTree::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.

  // WARNING when a selector is used with a TChain, you must use
  //  the pointer to the current TTree to call GetEntry(entry).
  //  The entry is always the local entry number in the current tree.
  //  Assuming that fChain is the pointer to the TChain being processed,
  //  use fChain->GetTree()->GetEntry(entry).

  // All we need is the single branch of interest. Get the current entry.
  br_weight_vec_->GetEntry( entry );

  // If we haven't initialized the sum vectors yet, we can do so now because
  // we have the first entry (and therefore the needed size)
  if ( sums_of_weights_.empty() ) {
    // Initialize the sum vectors with the appropriate number of elements (all
    // zero to start)
    size_t num_elements = weight_vec_->size();
    sums_of_weights_.assign( num_elements, 0. );
    sums_of_weights2_.assign( num_elements, 0. );
  }

  // Add the event weights in the current entry to the sum vectors
  for ( size_t w = 0u; w < weight_vec_->size(); ++w ) {
    // No need to use the slightly slower "at" here since we're directly
    // looping over the weight vector
    double wgt = weight_vec_->operator[]( w );
    double safe_wgt = safe_weight( wgt );

    // Use "at" here, just in case
    sums_of_weights_.at( w ) += safe_wgt;
    sums_of_weights2_.at( w ) += safe_wgt;
  }

  return kTRUE;
}

void MySelector::Begin( TTree* ) {
  // Executes at start of query only on the client
}

void MySelector::SlaveBegin( TTree* tree ) {
  // Executes at start of query (on all worker nodes if using PROOF)
  // Make histograms here, etc. as needed

  // Reset the owned vectors of summed event weights
  sums_of_weights_.clear();
  sums_of_weights2_.clear();
}

void MySelector::SlaveTerminate() {
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
}

void MySelector::Terminate() {
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  std::cout << sums_of_weights_.size() << ' '
    << sums_of_weights_.front() << '\n';
}
