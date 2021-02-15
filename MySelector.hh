#pragma once

// Standard library includes
#include <cmath>

// STV analysis includes
#include "TreeUtils.hh" // Needed for the MyPointer class template
#include "WeightHandler.hh"

constexpr double MIN_WEIGHT = 0.;
constexpr double MAX_WEIGHT = 30.;

// Branch names for special event weights
const std::string SPLINE_WEIGHT_NAME = "weight_splines_general_Spline";
const std::string TUNE_WEIGHT_NAME = "weight_TunedCentralValue_UBGenie";

// Event weights that are below MIN_WEIGHT, above MAX_WEIGHT, infinite, or NaN
// are reset to unity by this function. Other weights are returned unaltered.
inline double safe_weight( double w ) {
  if ( std::isfinite(w) && w >= MIN_WEIGHT && w <= MAX_WEIGHT ) return w;
  else return 1.0;
}

// Utility function used to check endings of (trimmed) weight labels based on
// branch names in the weights TTree
bool string_has_end( const std::string& str, const std::string& end ) {
  if ( str.length() >= end.length() ) {
    int comp = str.compare( str.length() - end.length(),
      end.length(), end );
    bool test_result = ( comp == 0 );
    return test_result;
  }
  return false;
}

class MySelector : public TSelector {

  public:

    MySelector( const std::vector<std::string>* branch_name_vec = nullptr )
      : TSelector(), branch_names_( branch_name_vec ) {}

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

    inline void GetBranchEntries( Long64_t entry ) {
      for ( auto* br : input_branches_ ) {
        br->GetEntry( entry );
      }
    }

    void apply_cv_correction_weights( const std::string& wgt_name, double& wgt,
      double spline_weight, double tune_weight );

    void InitializeSumVectors();
    void SetBranchPointers();

    // Input TTree or TChain (not owned by this class)
    TTree* chain_ = nullptr;

    // Optional vector of branch names to include in the output
    // (also not owned by this class)
    const std::vector< std::string >* branch_names_;

    // Pointers to branches of interest for the analysis
    std::vector< TBranch* > input_branches_;

    // Automatically manages temporary storage to use for reading event weights
    // from the TChain
    WeightHandler wh_;

    // Map that holds summed event weights organized by TTree branch
    std::map< std::string, std::vector<double> > sums_of_weights_;

    // Map that holds summed squares of event weights organized by TTree branch
    std::map< std::string, std::vector<double> > sums_of_weights2_;

};

void MySelector::Init( TTree* tree ) {

  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses of the tree
  // will be set. Init() will be called many times when running with PROOF.

  if ( !tree ) return;

  // Store the input TTree for later use
  chain_ = tree;

  // Set up the branch addresses
  wh_.set_branch_addresses( *tree, branch_names_ );

  // Make sure that we always have branches set up for the CV correction
  // weights, i.e., the spline and tune weights
  wh_.add_branch( *tree, SPLINE_WEIGHT_NAME );
  wh_.add_branch( *tree, TUNE_WEIGHT_NAME );

  // Get pointers to each branch of interest
  this->SetBranchPointers();

  // Get the vectors ready to store summed event weights
  InitializeSumVectors();

}

Bool_t MySelector::Notify() {

  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. Typically here the branch pointers
  // will be retrieved.
  this->SetBranchPointers();

  return kTRUE;
}

// Multiplies a given weight value by extra correction factors as appropriate.
// TODO: include the rootino_fix weight as a correction to the central value
void MySelector::apply_cv_correction_weights( const std::string& wgt_name,
  double& wgt, double spline_weight, double tune_weight )
{
  if ( string_has_end(wgt_name, "UBGenie") ) {
    wgt *= spline_weight;
  }
  else if ( wgt_name == "weight_flux_all"
    || wgt_name == "weight_reint_all"
    || wgt_name == "weight_xsr_scc_Fa3_SCC"
    || wgt_name == "weight_xsr_scc_Fv3_SCC" )
  {
    wgt *= spline_weight * tune_weight;
  }
  else if ( wgt_name == SPLINE_WEIGHT_NAME ) {
    // No extra weight factors needed
    return;
  }
  else throw std::runtime_error( "Unrecognized weight name" );
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

  // Get the current entry for all input TTree branches of interest
  this->GetBranchEntries( entry );

  // Get the CV correction weights here for potentially frequent re-use below
  double spline_weight = wh_.weight_map().at(
    SPLINE_WEIGHT_NAME )->front();

  double tune_weight = wh_.weight_map().at(
    TUNE_WEIGHT_NAME )->front();

  for ( const auto& pair : wh_.weight_map() ) {

    const std::string& br_name = pair.first;
    const auto& weight_vec = pair.second;

    auto& sum_vec = sums_of_weights_.at( br_name );
    auto& sum2_vec = sums_of_weights2_.at( br_name );

    // Add the event weights in the current entry to the sum vectors for
    // the appropriate branch
    for ( size_t w = 0u; w < weight_vec->size(); ++w ) {
      // No need to use the slightly slower "at" here since we're directly
      // looping over the weight vector
      double wgt = weight_vec->operator[]( w );

      // Multiply by any needed CV correction weights
      this->apply_cv_correction_weights( br_name, wgt,
        spline_weight, tune_weight );

      // Deal with NaNs, etc. to make a "safe weight" in all cases
      double safe_wgt = safe_weight( wgt );

      // Use "at" here, just in case
      sum_vec.at( w ) += safe_wgt;
      sum2_vec.at( w ) += safe_wgt;
    } // universe loop
  } // weight branch loop

  return kTRUE;
}

void MySelector::Begin( TTree* ) {
  // Executes at start of query only on the client
}

void MySelector::SlaveBegin( TTree* tree ) {
  // Executes at start of query (on all worker nodes if using PROOF)
  // Make histograms here, etc. as needed

  // Zero out the elements of all vectors of summed event weights
  for ( auto& pair : sums_of_weights_ ) {
    auto& sum_vec = pair.second;
    for ( double& sum : sum_vec ) sum = 0.;
  }
  for ( auto& pair : sums_of_weights2_ ) {
    auto& sum2_vec = pair.second;
    for ( double& sum2 : sum2_vec ) sum2 = 0.;
  }
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
  for ( const auto& pair : sums_of_weights_ ) {
    std::cout << pair.first << ' ' << pair.second.size()
      << ' ' << pair.second.front() << '\n';
  }
}

void MySelector::InitializeSumVectors() {
  // Initialize the sum vectors with the appropriate number of elements (all
  // zero to start). To get the sizes, load the first entry in each branch so
  // that we have example vectors available. Here and elsewhere in this class,
  // we assume that the vectors' sizes do not change throughout the TTree
  // entries.
  GetBranchEntries( 0 );

  for ( auto& pair : wh_.weight_map() ) {
    const std::string& br_name = pair.first;
    const auto& weight_vec = pair.second;
    size_t num_elements = weight_vec->size();

    sums_of_weights_[ br_name ] = std::vector<double>( num_elements, 0. );
    sums_of_weights2_[ br_name ] = std::vector<double>( num_elements, 0. );
  }
}

void MySelector::SetBranchPointers() {
  input_branches_.clear();

  if ( chain_ ) {
    // Get pointers to the branch(es) of interest
    for ( const auto& pair : wh_.weight_map() ) {
      const std::string& br_name = pair.first;
      TBranch* br = chain_->GetBranch( br_name.c_str() );
      input_branches_.push_back( br );
    }
  }
}
