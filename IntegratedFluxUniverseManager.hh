#pragma once

// Standard library includes
#include <memory>
#include <string>
#include <vector>

// ROOT includes
#include "TFile.h"
#include "TH1D.h"

constexpr size_t NUM_FLUX_UNIVERSES = 1000u;

class IntegratedFluxUniverseManager {

  public:

    // Deleted copy constructor
    IntegratedFluxUniverseManager(
      const IntegratedFluxUniverseManager& ) = delete;

    // Deleted move constructor
    IntegratedFluxUniverseManager( IntegratedFluxUniverseManager&& ) = delete;

    // Deleted copy assignment operator
    IntegratedFluxUniverseManager& operator=(
      const IntegratedFluxUniverseManager& ) = delete;

    // Deleted move assignment operator
    IntegratedFluxUniverseManager& operator=(
      IntegratedFluxUniverseManager&& ) = delete;

    // Get a const reference to the singleton instance of the
    // IntegratedFluxUniverseManager
    static const IntegratedFluxUniverseManager& Instance() {
      // Create the singleton object using a static variable. This ensures
      // that it is only created once.
      static std::unique_ptr< IntegratedFluxUniverseManager >
        the_instance( new IntegratedFluxUniverseManager() );

      // Return a reference to the singleton instance
      return *the_instance;
    }

    inline double flux_factor( size_t universe ) const {
      return integrated_flux_factors_.at( universe );
    }

    // Return the CV integrated numu flux (numu / POT / cm^2)
    inline double cv_numu_integrated_flux() const { return numu_CV_flux_; }

  protected:

    // Scaling factors used to vary the integrated flux across flux universes.
    // Upon construction of the singleton object, this vector is filled with
    // the scaling factors. They are indexed by universe number (zero-based).
    std::vector< double > integrated_flux_factors_;

    // Integrated central-value BNB numu flux for the MicroBooNE active volume
    // (numu / POT / cm^2)
    double numu_CV_flux_;

    IntegratedFluxUniverseManager() {
      // Get the latest MCC9 flux histograms from the standard location
      //TFile in_file( "/pnfs/uboone/persistent/uboonebeam/bnb_gsimple/"
      //  "bnb_gsimple_fluxes_01.09.2019_463_hist/MCC9_FluxHist_volTPC"
      //  "Active.root" );

      // Use a local copy in /uboone/app/ to deal with /pnfs outages
      TFile in_file( "/uboone/app/users/gardiner/stv/analysis/"
        "MCC9_FluxHist_volTPCActive.root" );

      // Get the integrated CV numu flux. Note that, since factors like the
      // bin width will cancel out in the ratios that we'll take below, we
      // don't care about units here, nor do we care about the distinction
      // between TH1::Integral() and TH1::Integral( "width" ). All we need
      // to do is compute the integrated flux using a consistent procedure
      // between the CV and each of the systematic variation universes.
      TH1D* cv_hist = dynamic_cast< TH1D* >( in_file.Get("hEnumu_cv") );
      double cv_integral = cv_hist->Integral();

      // See /pnfs/uboone/persistent/uboonebeam/bnb_gsimple/
      // bnb_gsimple_fluxes_01.09.2019_463_hist/readme.txt for an explanation
      // of this calculation.
      numu_CV_flux_ = cv_integral / ( 4997.*5e8 ) / ( 256.35*233. );

      for ( size_t u = 0u; u < NUM_FLUX_UNIVERSES; ++u ) {
        std::string hist_name = "numu_ms_total/hEnumu_ms_"
          + std::to_string( u );

        TH1D* univ_hist = dynamic_cast< TH1D* >(
          in_file.Get( hist_name.c_str() )
        );
        double univ_integral = univ_hist->Integral();

        integrated_flux_factors_.push_back( univ_integral / cv_integral );
      } // flux universes
    }
};
