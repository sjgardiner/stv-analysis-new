// Run with genie -l to load the class dictionaries needed to use
// genie::flux::GSimpleNtpEntry and genie::flux::GSimpleNtpMeta

// ROOT craps out for dumb reasons when using the real header file, so just
// reproduce the class declarations here as a stopgap solution.


namespace genie {
  namespace flux {
    /// Small persistable C-struct -like classes that makes up the SimpleNtpFlux
    /// ntuple.  This is only valid for a particular flux window (no reweighting,
    /// no coordinate transformation available).
    ///
    /// Order elements from largest to smallest for ROOT alignment purposes

    /// GSimpleNtpEntry
    /// =========================
    /// This is the only required branch ("entry") of the "flux" tree
    class GSimpleNtpEntry {
    public:
      GSimpleNtpEntry();
      /* allow default copy constructor ... for now nothing special
         GSimpleNtpEntry(const GSimpleNtpEntry & info);
      */
      virtual ~GSimpleNtpEntry() { };
      void Reset();
      void Print(const Option_t* opt = "") const;
      friend ostream & operator << (ostream & stream, const GSimpleNtpEntry & info);

      Double_t   wgt;      ///< nu weight

      Double_t   vtxx;     ///< x position in lab frame
      Double_t   vtxy;     ///< y position in lab frame
      Double_t   vtxz;     ///< z position in lab frame
      Double_t   dist;     ///< distance from hadron decay

      Double_t   px;       ///< x momentum in lab frame
      Double_t   py;       ///< y momentum in lab frame
      Double_t   pz;       ///< z momentum in lab frame
      Double_t   E;        ///< energy in lab frame

      Int_t      pdg;      ///< nu pdg-code
      UInt_t     metakey;  ///< key to meta data

      ClassDef(GSimpleNtpEntry,1)
    };

    /// GSimpleNtpMeta
    /// =========================
    /// A small persistable C-struct -like class that holds metadata
    /// about the the SimpleNtpFlux ntple.
    ///
    class GSimpleNtpMeta: public TObject {
    public:
      GSimpleNtpMeta();
      /* allow default copy constructor ... for now nothing special
         GSimpleNtpMeta(const GSimpleNtpMeta & info);
      */
      virtual ~GSimpleNtpMeta();

      void Reset();
      void AddFlavor(Int_t nupdg);
      void Print(const Option_t* opt = "") const;
      friend ostream & operator << (ostream & stream, const GSimpleNtpMeta & info);

      std::vector<Int_t>  pdglist; ///< list of neutrino flavors

      Double_t maxEnergy;   ///< maximum energy
      Double_t minWgt;      ///< minimum weight
      Double_t maxWgt;      ///< maximum weight
      Double_t protons;     ///< represented number of protons-on-target

      Double_t windowBase[3]; ///< x,y,z position of window base point
      Double_t windowDir1[3]; ///< dx,dy,dz of window direction 1
      Double_t windowDir2[3]; ///< dx,dy,dz of window direction 2

      std::vector<std::string>    auxintname;  ///< tagname of aux ints associated w/ entry
      std::vector<std::string>    auxdblname;  ///< tagname of aux doubles associated w/ entry
      std::vector<std::string>    infiles; ///< list of input files

      Int_t    seed;     ///< random seed used in generation
      UInt_t   metakey;  ///< index key to tie to individual entries

      static UInt_t mxfileprint;  ///< allow user to limit # of files to print

      ClassDef(GSimpleNtpMeta,1)
    };
  }
}

// STV analysis includes
#include "FiducialVolume.hh"
//#include "ActiveVolume.hh"

constexpr double METER_TO_CENTIMETER = 1e2;

struct FiducialVolumeFace {
  FiducialVolumeFace( const TVector3& norm,
    const TVector3& min, const TVector3& max ) : normal_( norm ),
    lower_( min ), upper_( max ) {}

  // Vector normal to the plane formed by this face
  TVector3 normal_;

  // Minimum (x,y,z) values for the face
  TVector3 lower_;

  // Maximum (x,y,z) values for the face
  TVector3 upper_;
};

const std::array< FiducialVolumeFace, 6 > FV_FACES = {
  // front
  FiducialVolumeFace( TVector3(0., 0., -1.),
    TVector3(FV_X_MIN, FV_Y_MIN, FV_Z_MIN),
    TVector3(FV_X_MAX, FV_Y_MAX, FV_Z_MIN) ),
  // back
  FiducialVolumeFace( TVector3(0., 0., 1.),
    TVector3(FV_X_MIN, FV_Y_MIN, FV_Z_MAX),
    TVector3(FV_X_MAX, FV_Y_MAX, FV_Z_MAX) ),
  // left
  FiducialVolumeFace( TVector3(-1., 0., 0.),
    TVector3(FV_X_MIN, FV_Y_MIN, FV_Z_MIN),
    TVector3(FV_X_MIN, FV_Y_MAX, FV_Z_MAX) ),
  // right
  FiducialVolumeFace( TVector3(1., 0., 0.),
    TVector3(FV_X_MAX, FV_Y_MIN, FV_Z_MIN),
    TVector3(FV_X_MAX, FV_Y_MAX, FV_Z_MAX) ),
  // top
  FiducialVolumeFace( TVector3(0., 1., 0.),
    TVector3(FV_X_MIN, FV_Y_MAX, FV_Z_MIN),
    TVector3(FV_X_MAX, FV_Y_MAX, FV_Z_MAX) ),
  // bottom
  FiducialVolumeFace( TVector3(0., -1., 0.),
    TVector3(FV_X_MIN, FV_Y_MIN, FV_Z_MIN),
    TVector3(FV_X_MAX, FV_Y_MIN, FV_Z_MAX) )
};

// Neutrino ray position as a function of a length parameter s (cm)
TVector3 ray_pos( const TVector3& vtx_pos, const TVector3& dir, double s ) {
  TVector3 result = vtx_pos;
  result += s * dir;
  return result;
}

double track_length_in_FV( const TVector3& vtx_pos, const TVector3& mom3 ) {
  // Create a unit vector in the neutrino's direction of flight
  TVector3 unit_mom3 = mom3.Unit();

  // Lengths from the neutrino ray vertex (on the flux window) to
  // the point(s) of intersection with the fiducial volume
  std::vector<double> intersection_lengths;

  // Approach based on the section "Line-Plane Intersection" from
  // http://geomalgorithms.com/a05-_intersect-1.html
  for ( const auto& face : FV_FACES ) {

    // Check that the face normal vector and the neutrino ray direction
    // are not perpendicular
    double n_dot_u = face.normal_.Dot( unit_mom3 );
    // If they are, then they are either coincident or disjoint, so we can
    // skip searching for the point of intersection.
    if ( n_dot_u == 0. ) continue;

    // Use the lower corner of the current face as a reference point on the
    // plane. Solve for the distance s from the neutrino ray vertex to the
    // point of intersection between the ray and the plane.
    TVector3 minus_w = face.lower_ - vtx_pos;
    double s = face.normal_.Dot( minus_w ) / n_dot_u;

    // Negative s values are *behind* the flux window, so skip them.
    if ( s < 0. ) continue;

    // We have found the distance from the neutrino ray origin to the
    // intersection with the plane. Use it to evaluate the intersection
    // position.
    TVector3 intersect_pos = ray_pos( vtx_pos, unit_mom3, s );

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

    // We're working with a rectangular volume, so if we've already found two
    // intersections, then we're done.
    if ( intersection_lengths.size() == 2u ) break;
  }

  // If there are no intersections with the faces, then the path length within
  // the fiducial volume is trivially zero.
  size_t num_intersections = intersection_lengths.size();
  if ( num_intersections == 0u ) return 0.;

  // If there is only one intersection, then the neutrino ray origin is within
  // the fiducial volume. This shouldn't ever happen. Complain and return zero.
  else if ( num_intersections == 1u ) {
    std::cout << "ERROR: neutrino ray origin is within fiducial volume\n";
    return 0.;
  }

  // We shouldn't have more than two intersections for a fiducial volume that is
  // a rectangular prism, except at the corners (where we can double-count). We
  // prevent this via the break statement above.
  else if ( num_intersections > 2u ) {
    std::cout << "ERROR: too many intersections with fiducial volume\n";
    return 0.;
  }

  // If we've made it here, then there are exactly two intersection points. The
  // difference in ray lengths is the path length for the ray within the
  // fiducial volume.
  double s1 = intersection_lengths.front();
  double s2 = intersection_lengths.back();
  double length_inside = std::abs( s2 - s1 );

  return length_inside;
}

void testme() {
  TVector3 vtx3( 50., 0., -10. );
  TVector3 p3( 0., 0., 2. );
  double trkl = track_length_in_FV( vtx3, p3 );
  std::cout << "trkl = " << trkl << '\n';
}

void fv_flux_gsimple() {

  TChain meta( "meta" );
  TChain flux( "flux" );
  flux.Add( "/pnfs/uboone/persistent/uboonebeam/bnb_gsimple/bnb_gsimple_"
    "fluxes_01.09.2019_463/converted_beammc_wincorr_*.root" );
  meta.Add( "/pnfs/uboone/persistent/uboonebeam/bnb_gsimple/bnb_gsimple_"
    "fluxes_01.09.2019_463/converted_beammc_wincorr_*.root" );

  // Compute the total simulated beam exposure in protons-on-target (POT)
  // by summing the entries on the "protons" branch of the "meta" TTree
  // (one entry per file)
  genie::flux::GSimpleNtpMeta* meta_info = nullptr;
  meta.SetBranchAddress( "meta", &meta_info );
  double total_POT = 0.;

  // Loop over the "meta" TChain in a weird way since TChain::GetEntries() can
  // be slow when working with many files
  long meta_entry = 0;
  while ( true ) {

    // TChain::LoadTree() returns the entry number that should be used with
    // the current TTree object, which (together with the TBranch objects
    // that it owns) doesn't know about the other TTrees in the TChain.
    // If the return value is negative, there was an I/O error, or we've
    // attempted to read past the end of the TChain.
    int local_entry = meta.LoadTree( meta_entry );

    // If we've reached the end of the TChain (or encountered an I/O error),
    // then terminate the event loop
    if ( local_entry < 0 ) break;

    // Load all of the branches for which we've called
    // TChain::SetBranchAddress() above
    meta.GetEntry( meta_entry );

    total_POT += meta_info->protons;

    ++meta_entry;
  }

  std::cout << "A total beam exposure of " << total_POT << " was simulated.\n";

  // Set up the object needed to read the neutrino ray information
  // from each "flux" TTree entry
  genie::flux::GSimpleNtpEntry flux_info;
  genie::flux::GSimpleNtpEntry* flux_info_ptr = &flux_info;
  flux.SetBranchAddress( "entry", &flux_info_ptr );

  // Set up the output TFile and TTree
  double trkl = 0.;
  TFile* out_file = new TFile( "/uboone/data/users/gardiner/nu_rays.root", "recreate" );
  //TFile* out_file = new TFile( "/uboone/data/users/gardiner/nu_rays_activeVol.root", "recreate" );
  TTree* out_tree = new TTree( "nu_ray_tree", "nu_ray_tree" );
  out_tree->Branch( "E", &flux_info.E, "E/D" );
  out_tree->Branch( "px", &flux_info.px, "px/D" );
  out_tree->Branch( "py", &flux_info.py, "py/D" );
  out_tree->Branch( "pz", &flux_info.pz, "pz/D" );
  out_tree->Branch( "pdg", &flux_info.pdg, "pdg/I" );
  out_tree->Branch( "vtxx", &flux_info.vtxx, "vtxx/D" );
  out_tree->Branch( "vtxy", &flux_info.vtxy, "vtxy/D" );
  out_tree->Branch( "vtxz", &flux_info.vtxz, "vtxz/D" );
  out_tree->Branch( "wgt", &flux_info.wgt, "wgt/D" );
  out_tree->Branch( "trkl", &trkl, "trkl/D" );

  TDatime start_time;
  std::cout << start_time.AsString() << '\n';

  // Use the same looping technique for the "flux" TChain
  long flux_entry = 0;
  while ( true ) {

    int local_entry = flux.LoadTree( flux_entry );
    if ( local_entry < 0 ) break;

    // Load all of the branches for which we've called
    // TChain::SetBranchAddress() above
    flux.GetEntry( flux_entry );
    ++flux_entry;

    // Convert the coordinates of the neutrino ray vertex on the flux window
    // from meters (as given in the gsimple files) to centimeters (as used
    // in this script)
    flux_info.vtxx *= METER_TO_CENTIMETER;
    flux_info.vtxy *= METER_TO_CENTIMETER;
    flux_info.vtxz *= METER_TO_CENTIMETER;

    TVector3 vtx3( flux_info.vtxx, flux_info.vtxy, flux_info.vtxz );
    TVector3 p3( flux_info.px, flux_info.py, flux_info.pz );

    trkl = track_length_in_FV( vtx3, p3 );

    if ( flux_entry % 100000 == 0 ) {
      std::cout << "flux entry = " << flux_entry << ", trkl = " << trkl << '\n';
    }
    out_tree->Fill();
  }

  out_tree->Write();

  TDatime end_time;
  std::cout << end_time.AsString() << '\n';

}
