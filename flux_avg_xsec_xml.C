// Run with genie -l
#include <fstream>
#include <map>
#include <sstream>
#include <string>

#include "Framework/Conventions/Units.h"
#include "Framework/EventGen/GEVGDriver.h"
#include "Framework/EventGen/InteractionList.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Interaction/InitialState.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/XSecSplineList.h"
#include "Framework/Utils/RunOpt.h"

int A_TARGET = 40;
const double FLUX_INTEGRAL_RELATIVE_TOLERANCE = 1e-8;

double total_xsec( const genie::GEVGDriver& evg_driver, double Ev )
{
  // Get the list of available interactions from the event generation driver
  const genie::InteractionList* ilist = evg_driver.Interactions();

  double xsec_sum = 0.;

  for ( const auto* interaction : *ilist ) {
    const genie::Spline* spl = evg_driver.XSecSpline( interaction );
    xsec_sum += spl->Evaluate( Ev ) / genie::units::cm2;
  }

  return xsec_sum;
}

genie::GEVGDriver configure_evg_driver( int probe_pdg, int target_pdg )
{
  // Create and configure an event generation driver for the requested probe +
  // target combination
  genie::InitialState init_state( target_pdg, probe_pdg );

  auto* ro = genie::RunOpt::Instance();

  genie::GEVGDriver evg_driver;
  evg_driver.SetEventGeneratorList( ro->EventGeneratorList() );
  evg_driver.Configure( init_state );
  evg_driver.CreateSplines();
  //evg_driver.CreateXSecSumSpline (100, gEmin, gEmax);

  return evg_driver;
}

void load_splines(const std::string& spline_xml_filename,
  const std::string& tune_name, const std::string& event_generator_list )
{
  genie::utils::app_init::MesgThresholds( "Messenger_whisper.xml" );

  auto* ro = genie::RunOpt::Instance();
  ro->SetTuneName( tune_name );
  ro->SetEventGeneratorList( event_generator_list );
  ro->BuildTune();

  // Load the cross section splines from the input XML file
  auto* xspl_list = genie::XSecSplineList::Instance();
  xspl_list->SetCurrentTune( tune_name );
  xspl_list->LoadFromXml( spline_xml_filename );

  if ( !xspl_list->HasSplineFromTune(tune_name) ) {
    std::cout << "Couldn't find splines for tune " << tune_name
      << " in the XML file " << spline_xml_filename;
    std::exit( 1 );
  }

}

void flux_avg_xsec_xml()
{

  // Default GENIE model for MicroBooNE
 /* load_splines( "/cvmfs/uboone.opensciencegrid.org/products/genie_xsec/v3_00"
    "_04_ub2/NULL/G1810a0211a-k250-e1000/data/gxspl-FNALsmall.xml",
    "G18_10a_02_11a", "CCinclMEC" );*/

  // Alternate GENIE model for MicroBooNE
<<<<<<< HEAD
  load_splines("/cvmfs/larsoft.opensciencegrid.org/products/genie_xsec/v3_00_04a/NULL/G0000b00000-k250-e1000/data/gxspl-FNALsmall.xml");
=======
  load_splines("/cvmfs/larsoft.opensciencegrid.org/products/genie_xsec/v3_00_04a/NULL/G0000b00000-k250-e1000/data/gxspl-FNALsmall.xml","G00_00b_00_000","CCinclMEC");
>>>>>>> 4ec0b6fc40041b4e79f334b4e106611470e1840d
  // TODO: add me, then comment out the default model

  auto evgd_numu = configure_evg_driver( 14, 1000180400 );

  TFile flux_hist_file( "/uboone/data/users/gardiner/linh/BNB_uBooNE_numu_flux_2019.root", "read" );
  TH1* flux_hist_numu = nullptr;
  flux_hist_file.GetObject( "numu", flux_hist_numu );

  // Normalize the flux histograms to unity to use for cross section weighting
  double flux_integral_numu = flux_hist_numu->Integral( "width" );
  //std::cout << "Old flux norm = " << flux_integral_numu << '\n';
  flux_hist_numu->Scale( 1. / flux_integral_numu );
  //std::cout << "New flux norm = " << flux_hist_numu->Integral("width") << '\n';

  double max_energy_numu = flux_hist_numu->GetBinLowEdge( flux_hist_numu->GetNbinsX() + 1 );

  // Functions to use for numerical integration. Element x[0] is the
  // neutrino energy
  TF1 xsec_weighted_flux_func_numu("temp_func", [&](double* x, double*)
    {
      int flux_bin = flux_hist_numu->FindBin( x[0] );
      double flux = flux_hist_numu->GetBinContent( flux_bin );

      double xsec = total_xsec( evgd_numu, x[0] );

      return flux * xsec;
    }, 0., max_energy_numu, 0);

    //const genie::InteractionList* ilist = evgd_numu.Interactions();
    //for ( const auto* interaction : *ilist ) {
    //  std::cout << interaction->AsString() << '\n';
    //}

    double flux_avg_xsec_numu = xsec_weighted_flux_func_numu.Integral(0., max_energy_numu, FLUX_INTEGRAL_RELATIVE_TOLERANCE);

    std::cout << "flux_avg_xsec_numu = " << flux_avg_xsec_numu << " cm^2 / nucleus\n";
    double flux_avg_xsec_per_nucleon_numu = flux_avg_xsec_numu / A_TARGET;
    std::cout << "flux_avg_xsec_numu = " << flux_avg_xsec_per_nucleon_numu << " cm^2 / nucleon\n";

    TFile out_file( "flux_avg_xsecs.root", "recreate" );

    TParameter<double> flux_avg_numu_ccincl_xsec( "flux_avg_ccincl_xsec_per_nucleon_numu",
      flux_avg_xsec_per_nucleon_numu );
    flux_avg_numu_ccincl_xsec.Write();
}

