#ifndef FpmcInterface_Fpmc_h
#define FpmcInterface_Fpmc_h

/** \class Fpmc
 *
 * Generates Fpmc HepMC events
 *
 ***************************************/

#include "HepMC/GenEvent.h"
#include "HepMC/HEPEVT_Wrapper.h"
#include "HepMC/Version.h"

#ifndef HEPMC_VERSION_CODE // HepMC v2
#define HEPMC_VERSION2

#include "HepMC/IO_HERWIG.h"
#include "HepMC/PdfInfo.h"

#else // HepMC v>=3

#include "HepMC/WriterAscii.h"
#include "HepMC/GenPdfInfo.h"

#endif

#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"

#include "FpmcParameters.h"

#include <vector>
#include <string>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <iostream>

namespace fpmc
{
  class Fpmc
  {
    public:
      Fpmc( double, const char* );
      ~Fpmc();

      void begin();
      bool run();
      void end();
      void write( const char* );

      /// Retrieve the last event generated
      const HepMC::GenEvent* event() const { return hepMCEvt_.get(); } 

    private:
      void initialiseParams();
#ifdef HEPMC_VERSION2
      HepMC::IO_HERWIG conv_;
#endif
      /// Last event generated
      std::shared_ptr<HepMC::GenEvent> hepMCEvt_;
 
      /// HERWIG verbosity
      unsigned int herwigVerbosity_;
      /// HepMC verbosity
      bool hepMCVerbosity_;
      /// Events to print if verbosity
      unsigned int maxEventsToPrint_;

      /// Number of events already generated
      unsigned int event_;
      /// Centre of mass energy for the initial system
      double comEnergy_;

      /// List of parameters obtained from the steering card
      FpmcParameters params_;

      /// Enable/disable the hadronisation
      bool hadronize_;
      /// Enable/disable the extra printout
      bool debug_;
      std::ostream& dbg_;
  };
} 

#endif
