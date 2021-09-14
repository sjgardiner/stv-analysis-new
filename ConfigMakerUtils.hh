#pragma once
// Provides variable definitions for possible use in multiple configuration
// file writer scripts

#include <map>
#include <string>
#include <vector>

// These strings provide non-default selections (in TTree::Draw format) for
// sideband control samples.
const std::string DIRT_SIDEBAND_SELECTION =
  "!sel_reco_vertex_in_FV && sel_pfp_starts_in_PCV"
  " && sel_has_muon_candidate && sel_topo_cut_passed"
  " && sel_no_reco_showers && sel_muon_above_threshold"
  " && sel_has_p_candidate && sel_passed_proton_pid_cut"
  " && sel_protons_contained && sel_lead_p_passed_mom_cuts";

const std::string NC_SIDEBAND_SELECTION =
  "sel_reco_vertex_in_FV && sel_pfp_starts_in_PCV"
  " && !sel_has_muon_candidate && sel_topo_cut_passed"
  " && sel_no_reco_showers && sel_has_p_candidate"
  " && sel_passed_proton_pid_cut && sel_protons_contained"
  " && sel_lead_p_passed_mom_cuts";

const std::string CCNPI_SIDEBAND_SELECTION =
  "sel_reco_vertex_in_FV && sel_pfp_starts_in_PCV"
  " && sel_has_muon_candidate && sel_topo_cut_passed"
  " && sel_no_reco_showers && sel_muon_above_threshold"
  " && sel_has_p_candidate && sel_protons_contained"
  " && sel_lead_p_passed_mom_cuts"
  " && trk_llr_pid_score_v[ lead_p_candidate_idx ] > 0.2 ";

// Bin definitions for the 2D cross-section measurements

  // Using floating-point numbers as std::map keys is admittedly evil, but
  // it's safe in this case: all we'll do with this map is iterate over the
  // elements. Keys are muon momentum bin edges, values are muon scattering
  // cosine bin edges.
  std::map< double, std::vector<double> > MUON_2D_BIN_EDGES = {

    // No need for an underflow bin: due to the signal definition, all muons
    // with reco momentum below 0.1 GeV/c will be lost

    { 0.1,  { -1., 0., 1. } },

    { 0.17, { -1., -0.2, 0.4, 1. } },

    { 0.21, { -1, -0.2, 0.4, 1. } },

    { 0.24, { -1, -0.1, 0.5, 1. } },

    { 0.27, { -1, -0.1, 0.35, 0.6, 1. } },

    { 0.3,  { -1, -0.4, -0.1, 0.1, 0.35, 0.5, 0.7, 0.85, 1. } },

    { 0.38, { -1, 0, 0.5, 0.65, 0.8, 0.92, 1.00 } },

    { 0.48, { -1, 0.2, 0.5, 0.65, 0.8, 0.875, 0.950, 1.00 } },

    { 0.75, { -1, 0.5, 0.8, 0.875, 0.950, 1.00 } },

    { 1.14, { -1, 0.85, 0.9, 0.950, 1.00 } },

    // Upper edge of the last bin. The script will create an overflow bin above
    // this.
    { 2.5, {} }

  };

  // Using floating-point numbers as std::map keys is admittedly evil, but it's
  // safe in this case: all we'll do with this map is iterate over the
  // elements. Keys are proton momentum bin edges, values are proton cosine bin
  // edges.
  std::map< double, std::vector<double> > PROTON_2D_BIN_EDGES = {

    // No need for an underflow bin: due to the signal definition, all leading
    // protons with reco momentum below 0.25 GeV/c will be lost
    { 0.250, { -1, -0.5, 0.1, 0.6, 1.0 } },
    { 0.325, { -1, -0.7, -0.4, 0, 0.4, 0.6, 0.8, 1.0 } },
    { 0.4,   { -1, -0.6, -0.2, 0.2, 0.5, 0.65, 0.85, 1.0 } },
    { 0.45,  { -1, -0.4, 0, 0.2, 0.4, 0.55, 0.65, 0.8, 0.92, 1.0 } },
    { 0.5,   { -1, -0.2, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 } },
    { 0.550, { -1, 0, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 } },
    { 0.6,   { -1, 0.1, 0.37, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 } },
    { 0.65,  { -1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 } },
    { 0.7,   { -1, 0.45, 0.65, 0.75, 0.82, 0.9, 1.0 } },
    { 0.75,  { -1, 0.55, 0.7, 0.8, 0.87, 1.0 } },
    { 0.8,   { -1, 0.65, 0.78, 0.89, 1.0 } },
    { 0.85,  { -1, 0.73, 0.86, 1.0 } },
    { 0.9,   { -1, 0.77, 0.88, 1.0 } },
    { 0.975, { -1, 0.84, 1.0 } },

    // Upper edge of the last bin. We don't need an overflow bin because the
    // signal definition excludes any leading protons with momenta above 1.2
    // GeV/c
    { 1.2, {} }

  };
