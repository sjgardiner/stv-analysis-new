#pragma once

#include "TH1.h"

// Enum used to label event categories of interest for analysis plots
enum EventCategory {

  // Unable to categorize (e.g., because the event is real data and thus
  // has no MC truth information)
  kUnknown = 0,

  // Signal events broken down by underlying reaction mode
  kSignalCCQE = 1,
  kSignalCCMEC = 2,
  kSignalCCRES = 3,
  kSignalOther = 4,

  // True numu CC event with at least one final-state pion above threshold
  kNuMuCCNpi = 5,

  // True numu CC event with zero final-state pions above threshold and
  // zero final-state protons above threshold
  kNuMuCC0pi0p = 6,

  // Any true numu CC event which does not satisfy the criteria for inclusion
  // in one of the other categories above
  kNuMuCCOther = 7,

  // True nue CC event
  kNuECC = 8,

  // True neutral current event for any neutrino flavor
  kNC = 9,

  // True neutrino vertex (any reaction mode and flavor combination) is outside
  // of the fiducial volume
  kOOFV = 10,

  // All events that do not fall within any of the other categories (e.g.,
  // numubar CC)
  kOther = 11

};

// Singleton class that helps manipulate EventCategory enum values
class EventCategoryInterpreter {

  public:

    // This is a singleton class, so we'll delete the copy constructor
    // the move constructor, and the assignment operators
    EventCategoryInterpreter( const EventCategoryInterpreter& ) = delete;
    EventCategoryInterpreter( EventCategoryInterpreter&& ) = delete;
    EventCategoryInterpreter& operator=( const EventCategoryInterpreter& )
      = delete;
    EventCategoryInterpreter& operator=( EventCategoryInterpreter&& )
      = delete;

    // Get a const reference to the singleton instance of the
    // EventCategoryInterpreter
    inline static const EventCategoryInterpreter& Instance() {

      // Create the EventCategoryInterpreter object using a static variable.
      // This ensures that the singleton instance is only created once.
      static std::unique_ptr<EventCategoryInterpreter>
        the_instance( new EventCategoryInterpreter() );

      // Return a reference to the singleton instance
      return *the_instance;
    }

    inline const std::map< EventCategory, std::string >& label_map() const
      { return event_category_to_label_map_; }

    inline std::string label( EventCategory ec ) const
      { return event_category_to_label_map_.at( ec ); }

    inline int color_code( EventCategory ec ) const
      { return event_category_to_color_map_.at( ec ); }

    inline void set_mc_histogram_style( EventCategory ec, TH1* mc_hist ) const
    {
      int color = color_code( ec );
      mc_hist->SetFillColor( color );
      mc_hist->SetLineColor( color );
      mc_hist->SetStats( false );
    }

    inline void set_ext_histogram_style( TH1* ext_hist ) const {
      ext_hist->SetFillColor( 28 );
      ext_hist->SetLineColor( 28 );
      ext_hist->SetLineWidth( 2 );
      ext_hist->SetFillStyle( 3005 );
      ext_hist->SetStats( false );
    }

    inline void set_bnb_data_histogram_style( TH1* bnb_hist ) const {

      bnb_hist->SetLineColor( kBlack );
      bnb_hist->SetLineWidth( 3 );
      bnb_hist->SetMarkerStyle( kFullCircle );
      bnb_hist->SetMarkerSize( 0.8 );
      bnb_hist->SetStats( false );

      bnb_hist->GetXaxis()->SetTitleOffset( 0.0 );
      bnb_hist->GetXaxis()->SetTitleSize( 0.0 );
      bnb_hist->GetYaxis()->SetTitleSize( 0.05 );
      bnb_hist->GetYaxis()->CenterTitle( true );
      bnb_hist->GetXaxis()->SetLabelSize( 0.0 );

      // This prevents the first y-axis label label (0) to be clipped by the
      // ratio plot
      bnb_hist->SetMinimum( 1e-3 );
    }

    inline void set_stat_err_histogram_style( TH1* stat_err_hist ) const {
      stat_err_hist->SetFillColor( kBlack );
      stat_err_hist->SetLineColor( kBlack );
      stat_err_hist->SetLineWidth( 2 );
      stat_err_hist->SetFillStyle( 3004 );
    }

  private:

    EventCategoryInterpreter() {}

    std::map< EventCategory, std::string > event_category_to_label_map_ = {
      { kUnknown, "Unknown" },
      { kSignalCCQE, "Signal (CCQE)" },
      { kSignalCCMEC, "Signal (CCMEC)" },
      { kSignalCCRES, "Signal (CCRES)" },
      { kSignalOther, "Signal (Other)" },
      { kNuMuCCNpi, "#nu_{#mu} CCN#pi" },
      { kNuMuCC0pi0p, "#nu_{#mu} CC0#pi0p" },
      { kNuMuCCOther, "Other #nu_{#mu} CC" },
      { kNuECC, "#nu_{e} CC" },
      { kNC, "NC" },
      { kOOFV, "Out FV" },
      { kOther, "Other" }
    };

    std::map< EventCategory, int > event_category_to_color_map_ = {
      { kUnknown, kGray },
      { kSignalCCQE, kGreen },
      { kSignalCCMEC, kGreen + 1 },
      { kSignalCCRES, kGreen + 2 },
      { kSignalOther, kGreen + 3 },
      { kNuMuCCNpi, kAzure - 2 },
      { kNuMuCC0pi0p, kAzure - 1 },
      { kNuMuCCOther, kAzure },
      { kNuECC, kViolet },
      { kNC, kOrange },
      { kOOFV, kRed + 3 },
      { kOther, kRed + 1 }
    };
};
