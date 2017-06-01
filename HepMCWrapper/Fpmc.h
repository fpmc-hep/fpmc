#ifndef FpmcInterface_Fpmc_h
#define FpmcInterface_Fpmc_h

/** \class Fpmc
 *
 * Generates Fpmc HepMC events
 *
 ***************************************/

#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenPdfInfo.h"
//#include "HepMC/HerwigWrapper.h"
#include "HepMC/HEPEVT_Wrapper.h"
#include "HepMC/Version.h"

#ifdef HEPMC_VERSION_CODE // HepMC v>=3
#include "HepMC/WriterAscii.h"
#else // HepMC v2
#define HEPMC_VERSION2
#include "HepMC/IO_HERWIG.h"
#endif

#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"

#include "herwig.h"
#include "fpmc.h"

#include "FpmcParameters.h"

#include <vector>
#include <string>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <cstring>

namespace fpmc
{
  class Fpmc
  {
    public:
      Fpmc( double, long int, const char* );
      ~Fpmc();

      void begin();
      bool run();
      void end();
      void write( const char* );

      const HepMC::GenEvent* event() const { return hepMCEvt_.get(); } 

    private:
      void initialiseParams();
#ifdef HEPMC_VERSION2
      HepMC::IO_HERWIG conv_;
#endif
      std::shared_ptr<HepMC::GenEvent> hepMCEvt_;
 
      /// HERWIG verbosity
      unsigned int herwigVerbosity_;
      /// HepMC verbosity
      bool hepMCVerbosity_;
      /// Events to print if verbosity
      unsigned int maxEventsToPrint_;

      unsigned int event_;
      double comEnergy_;

      FpmcParameters params_;

      bool hadronize_;
      bool debug_;
  };
} 

#endif
