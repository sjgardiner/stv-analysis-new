#include "TVector3.h"
//// Only load the library in this way if we're using this code from inside the
//// ROOT C++ interpreter. We could check for __CINT__ as well, but the specific
//// R__LOAD_LIBRARY approach used here only works with ROOT 6.
//#ifdef __CLING__
//// Pre-load the definition of the TTreeFormula::EvalInstance function from the
//// TreePlayer shared library. This approach is based on the trick mentioned
//// here: https://tinyurl.com/2s4yuzxm
//R__LOAD_LIBRARY(libTreePlayer.so)
//#endif

// STV analysis includes
#include "FiducialVolume.hh"
#include "TreeUtils.hh"

constexpr int MUON = 13;
constexpr double TRACK_SCORE_CUT = 0.5;

// Boundaries of the containment volume in cm. Chosen to leave a 10-cm border
// between all edges and the active volume boundary
double VOL_X_MIN =   00.;
double VOL_X_MAX =  256.4;

double VOL_Y_MIN = -116.5;
double VOL_Y_MAX =  116.5;

double VOL_Z_MIN =   0.;
double VOL_Z_MAX = 1036.8;

bool in_active_volume( double x, double y, double z ) {
  bool x_inside_V = ( VOL_X_MIN < x ) && ( x < VOL_X_MAX );
  bool y_inside_V = ( VOL_Y_MIN < y ) && ( y < VOL_Y_MAX );
  bool z_inside_V = ( VOL_Z_MIN < z ) && ( z < VOL_Z_MAX );
  return ( x_inside_V && y_inside_V && z_inside_V );
}

constexpr double PCV_OFFSET = 10.; // cm

bool in_pcv( double x, double y, double z ) {
  bool x_inside_V = ( VOL_X_MIN + PCV_OFFSET < x )
    && ( x < VOL_X_MAX - PCV_OFFSET );
  bool y_inside_V = ( VOL_Y_MIN + PCV_OFFSET < y )
    && ( y < VOL_Y_MAX - PCV_OFFSET );
  bool z_inside_V = ( VOL_Z_MIN + PCV_OFFSET < z )
    && ( z < VOL_Z_MAX - PCV_OFFSET );
  return ( x_inside_V && y_inside_V && z_inside_V );
}

struct TempEvent {

  TempEvent() {}

  void set_branch_addresses( TTree& etree );

  MyPointer< std::vector<int> > pfp_true_pdg_;
  MyPointer< std::vector<float> > pfp_true_px_;
  MyPointer< std::vector<float> > pfp_true_py_;
  MyPointer< std::vector<float> > pfp_true_pz_;
  MyPointer< std::vector<float> > pfp_track_score_;
  MyPointer< std::vector<float> > track_mcs_mom_mu_;
  MyPointer< std::vector<float> > track_startx_;
  MyPointer< std::vector<float> > track_starty_;
  MyPointer< std::vector<float> > track_startz_;
  MyPointer< std::vector<float> > track_endx_;
  MyPointer< std::vector<float> > track_endy_;
  MyPointer< std::vector<float> > track_endz_;
  MyPointer< std::vector<float> > track_dirx_;
  MyPointer< std::vector<float> > track_diry_;
  MyPointer< std::vector<float> > track_dirz_;
  MyPointer< std::vector<float> > track_length_;
};

void TempEvent::set_branch_addresses( TTree& etree ) {

  set_object_input_branch_address( etree, "backtracked_pdg", pfp_true_pdg_ );
  set_object_input_branch_address( etree, "backtracked_px", pfp_true_px_ );
  set_object_input_branch_address( etree, "backtracked_py", pfp_true_py_ );
  set_object_input_branch_address( etree, "backtracked_pz", pfp_true_pz_ );
  set_object_input_branch_address( etree, "trk_mcs_muon_mom_v",
    track_mcs_mom_mu_ );
  set_object_input_branch_address( etree, "trk_score_v", pfp_track_score_ );

  set_object_input_branch_address( etree, "trk_sce_start_x_v", track_startx_ );
  set_object_input_branch_address( etree, "trk_sce_start_y_v", track_starty_ );
  set_object_input_branch_address( etree, "trk_sce_start_z_v", track_startz_ );

  set_object_input_branch_address( etree, "trk_sce_end_x_v", track_endx_ );
  set_object_input_branch_address( etree, "trk_sce_end_y_v", track_endy_ );
  set_object_input_branch_address( etree, "trk_sce_end_z_v", track_endz_ );

  set_object_input_branch_address( etree, "trk_dir_x_v", track_dirx_ );
  set_object_input_branch_address( etree, "trk_dir_y_v", track_diry_ );
  set_object_input_branch_address( etree, "trk_dir_z_v", track_dirz_ );

  set_object_input_branch_address( etree, "trk_len_v", track_length_ );
}

struct Rectangle {
  Rectangle( const TVector3& norm,
    const TVector3& min, const TVector3& max ) : normal_( norm ),
    lower_( min ), upper_( max ) {}

  // Vector normal to the plane formed by this face
  TVector3 normal_;

  // Minimum (x,y,z) values for the face
  TVector3 lower_;

  // Maximum (x,y,z) values for the face
  TVector3 upper_;
};

const std::array< Rectangle, 6 > VOL_FACES = {
  // front
  Rectangle( TVector3(0., 0., -1.),
    TVector3(VOL_X_MIN, VOL_Y_MIN, VOL_Z_MIN),
    TVector3(VOL_X_MAX, VOL_Y_MAX, VOL_Z_MIN) ),
  // back
  Rectangle( TVector3(0., 0., 1.),
    TVector3(VOL_X_MIN, VOL_Y_MIN, VOL_Z_MAX),
    TVector3(VOL_X_MAX, VOL_Y_MAX, VOL_Z_MAX) ),
  // left
  Rectangle( TVector3(-1., 0., 0.),
    TVector3(VOL_X_MIN, VOL_Y_MIN, VOL_Z_MIN),
    TVector3(VOL_X_MIN, VOL_Y_MAX, VOL_Z_MAX) ),
  // right
  Rectangle( TVector3(1., 0., 0.),
    TVector3(VOL_X_MAX, VOL_Y_MIN, VOL_Z_MIN),
    TVector3(VOL_X_MAX, VOL_Y_MAX, VOL_Z_MAX) ),
  // top
  Rectangle( TVector3(0., 1., 0.),
    TVector3(VOL_X_MIN, VOL_Y_MAX, VOL_Z_MIN),
    TVector3(VOL_X_MAX, VOL_Y_MAX, VOL_Z_MAX) ),
  // bottom
  Rectangle( TVector3(0., -1., 0.),
    TVector3(VOL_X_MIN, VOL_Y_MIN, VOL_Z_MIN),
    TVector3(VOL_X_MAX, VOL_Y_MIN, VOL_Z_MAX) )
};

// Ray position as a function of a length parameter s (cm)
TVector3 ray_pos( const TVector3& start_pos, const TVector3& dir, double s ) {
  TVector3 result = start_pos;
  result += s * dir;
  return result;
}

double track_length_in_vol( const TVector3& start_pos, const TVector3& mom3,
  int& num_intersections )
{
  // Create a unit vector in the indicated direction
  TVector3 unit_mom3 = mom3.Unit();

  // Lengths from the starting position to the point(s) of intersection with
  // the volume of interest
  std::vector<double> intersection_lengths;

  // Approach based on the section "Line-Plane Intersection" from
  // http://geomalgorithms.com/a05-_intersect-1.html
  for ( const auto& face : VOL_FACES ) {

    // Check that the face normal vector and the ray direction
    // are not perpendicular
    double n_dot_u = face.normal_.Dot( unit_mom3 );
    // If they are, then they are either coincident or disjoint, so we can
    // skip searching for the point of intersection.
    if ( n_dot_u == 0. ) continue;

    // Use the lower corner of the current face as a reference point on the
    // plane. Solve for the distance s from the starting point to the
    // point of intersection between the ray and the plane.
    TVector3 minus_w = face.lower_ - start_pos;
    double s = face.normal_.Dot( minus_w ) / n_dot_u;

    // Negative s values are *behind* the starting point, so skip them.
    if ( s < 0. ) continue;

    // We have found the distance from the starting point to the intersection
    // with the plane. Use it to evaluate the intersection position.
    TVector3 intersect_pos = ray_pos( start_pos, unit_mom3, s );

    // If the intersection point lies outside of the face boundaries, then skip
    // to the next face. Don't bother to check the dimension perpendicular to
    // the face.
    if ( face.lower_.X() != face.upper_.X() ) {
      if ( intersect_pos.X() < face.lower_.X() ) continue;
      if ( intersect_pos.X() > face.upper_.X() ) continue;
    }
    if ( face.lower_.Y() != face.upper_.Y() ) {
      if ( intersect_pos.Y() < face.lower_.Y() ) continue;
      if ( intersect_pos.Y() > face.upper_.Y() ) continue;
    }
    if ( face.lower_.Z() != face.upper_.Z() ) {
      if ( intersect_pos.Z() < face.lower_.Z() ) continue;
      if ( intersect_pos.Z() > face.upper_.Z() ) continue;
    }

    // We've passed all the criteria for a good intersection point, so add
    // it to the vector of distances.
    intersection_lengths.push_back( s );

    // We're working with a volume that is a rectangular prism, so if we've
    // already found two intersections, then we're done.
    if ( intersection_lengths.size() == 2u ) break;
  }

  // If there are no intersections with the faces, then the path length within
  // the volume is trivially zero.
  num_intersections = intersection_lengths.size();
  if ( num_intersections <= 0 ) return 0.;

  // If there is only one intersection, then the starting point is within
  // the volume.
  else if ( num_intersections == 1 ) {
    return intersection_lengths.front();
  }

  // We shouldn't have more than two intersections for a volume that is a
  // rectangular prism, except at the corners (where we can double-count). We
  // prevent this via the break statement above.
  else if ( num_intersections > 2 ) {
    std::cout << "ERROR: too many intersections with volume boundaries\n";
    return 0.;
  }

  // If we've made it here, then there are exactly two intersection points. The
  // difference in ray lengths is the path length for the ray within the
  // volume.
  double s1 = intersection_lengths.front();
  double s2 = intersection_lengths.back();
  double length_inside = std::abs( s2 - s1 );

  return length_inside;
}

void mcs_study() {

  TChain stv_tree( "stv_tree" );
  stv_tree.Add( "/uboone/data/users/gardiner/ntuples-stv-MCC9InternalNote/"
    "stv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_"
    "run1_reco2_reco2.root" );

  TFile out_file( "my_out_av.root", "recreate" );
  TTree* out_tree = new TTree( "out_tree", "out tree" );
  double trkl, pmu, mc_pmu, full_trkl;
  int num_intersections;
  out_tree->Branch( "trkl", &trkl, "trkl/D" );
  out_tree->Branch( "full_trkl", &full_trkl, "full_trkl/D" );
  out_tree->Branch( "num_inter", &num_intersections, "num_inter/I" );
  out_tree->Branch( "pmu", &pmu, "pmu/D" );
  out_tree->Branch( "mc_pmu", &mc_pmu, "mc_pmu/D" );

  // Loop in a weird way to avoid calling the potentially slow function
  // TChain::GetEntries()
  long entry = 0;
  while ( true ) {

    // Create a new TempEvent object. This will reset all analysis variables
    // for the current event.
    TempEvent cur_event;

    // Set branch addresses for the member variables that will be read
    // directly from the Event TTree.
    cur_event.set_branch_addresses( stv_tree );

    int local_entry = stv_tree.LoadTree( entry );
    if ( local_entry < 0 ) break;

    // Load all of the branches for which we've called
    // TChain::SetBranchAddress() above
    stv_tree.GetEntry( entry );
    ++entry;

    size_t num_pfparticles = cur_event.pfp_true_pdg_->size();
    for ( size_t p = 0u; p < num_pfparticles; ++p ) {
      // Skip all PFParticles except those corresponding to true simulated
      // muons
      int pdg = cur_event.pfp_true_pdg_->at( p );
      if ( pdg != MUON ) continue;

      // Skip muons reconstructed as showers
      double track_score = cur_event.pfp_track_score_->at( p );
      if ( track_score < TRACK_SCORE_CUT ) continue;

      // Skip muons whose reco start points aren't contained
      double start_x = cur_event.track_startx_->at( p );
      double start_y = cur_event.track_starty_->at( p );
      double start_z = cur_event.track_startz_->at( p );
      bool start_in_vol = in_pcv( start_x, start_y, start_z );
      if ( !start_in_vol ) continue;

      // Skip muons whose reco end points are contained
      double end_x = cur_event.track_endx_->at( p );
      double end_y = cur_event.track_endy_->at( p );
      double end_z = cur_event.track_endz_->at( p );
      bool end_in_vol = in_pcv( end_x, end_y, end_z );
      if ( end_in_vol ) continue;

      // Get the true muon momentum
      double true_px = cur_event.pfp_true_px_->at( p );
      double true_py = cur_event.pfp_true_py_->at( p );
      double true_pz = cur_event.pfp_true_pz_->at( p );

      TVector3 mc_p3mu( true_px, true_py, true_pz );
      mc_pmu = mc_p3mu.Mag();

      // Get the reconstructed muon momentum (via MCS)
      pmu = cur_event.track_mcs_mom_mu_->at( p );

      //// Get the reconstructed direction of the muon track
      //double dir_x = cur_event.track_dirx_->at( p );
      //double dir_y = cur_event.track_diry_->at( p );
      //double dir_z = cur_event.track_dirz_->at( p );

      double dir_x = end_x - start_x;
      double dir_y = end_y - start_y;
      double dir_z = end_z - start_z;

      // Get the reco track length that lies within the volume of interest
      TVector3 start_point3( start_x, start_y, start_z );
      TVector3 dir3( dir_x, dir_y, dir_z );
      trkl = track_length_in_vol( start_point3, dir3, num_intersections );

      // Get the full length of the reco track
      full_trkl = cur_event.track_length_->at( p );

      // Save the result to the output TTree
      out_tree->Fill();
      std::cout << "entry = " << entry << ", trkl = " << trkl
        << ", intersections = " << num_intersections << '\n';
    }
  }

  out_tree->Write();
}
