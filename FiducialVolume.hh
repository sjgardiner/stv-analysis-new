#pragma once
// Definition of the fiducial volume for the analysis, together with a utility
// function template to check if a given point lies within it or not

// Boundaries of the neutrino vertex fiducial volume (cm)
// This is handled the same way for reco and in MC
double FV_X_MIN =   10.;
double FV_X_MAX =  246.35;

double FV_Y_MIN = -106.5;
double FV_Y_MAX =  106.5;

double FV_Z_MIN =   10.;
double FV_Z_MAX =  968.8;

// Use a template here so that this function can take float or double values as
// input
template <typename Number> bool point_inside_FV( Number x, Number y, Number z )
{
  bool x_inside_FV = ( FV_X_MIN < x ) && ( x < FV_X_MAX );
  bool y_inside_FV = ( FV_Y_MIN < y ) && ( y < FV_Y_MAX );
  bool z_inside_FV = ( FV_Z_MIN < z ) && ( z < FV_Z_MAX );
  return ( x_inside_FV && y_inside_FV && z_inside_FV );
}

inline bool point_inside_FV( const TVector3& pos ) {
  return point_inside_FV( pos.X(), pos.Y(), pos.Z() );
}
