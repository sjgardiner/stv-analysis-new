#pragma once

// Needed for the weight range limits used in DEFAULT_MC_EVENT_WEIGHT. We pull
// them from the UniverseMaker header file to ensure consistency
// with the safe_weight() function.
#include "UniverseMaker.hh"

// **** Helper code to facilitate making histograms ****

// By default, weight the MC events using the MicroBooNE CV tune. The
// TTree::Draw() expression below includes some workarounds for problematic
// weights. They are equivalent to the safe_weight() function defined in
// UniverseMaker.hh.
const std::string DEFAULT_MC_EVENT_WEIGHT = "(std::isfinite(spline_weight*"
  "tuned_cv_weight) && spline_weight*tuned_cv_weight >= "
  + std::to_string( MIN_WEIGHT ) + " && spline_weight*"
  "tuned_cv_weight <= " + std::to_string( MAX_WEIGHT )
  + " ? spline_weight*tuned_cv_weight : 1)";

// Generates a vector of bin low edges equivalent to the approach used by the
// TH1 constructor that takes xmin and xmax in addition to the number of bins
std::vector<double> get_bin_low_edges( double xmin, double xmax, int Nbins )
{
  std::vector<double> bin_low_edges;
  double bin_step = ( xmax - xmin ) / Nbins;
  for ( int b = 0; b <= Nbins; ++b ) {
    double low_edge = xmin + b*bin_step;
    bin_low_edges.push_back( low_edge );
  }

  return bin_low_edges;
}
