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

#include "Fpmc.h"
#include "FpmcParameters.h"

#include <string>
#include <sstream>
#include <stdexcept>
#include <iostream>

namespace fpmc
{
  class Wrapper : public Fpmc
  {
    public:
      Wrapper( double, const char* );
      ~Wrapper();

      /// Retrieve the last event generated
      const HepMC::GenEvent* event();
      /// Write the last event generated onto a file
      void write( const char* );

    private:
#ifdef HEPMC_VERSION2
      HepMC::IO_HERWIG conv_;
#endif
      /// Last event generated
      std::shared_ptr<HepMC::GenEvent> hepMCEvt_;
 
      /// HepMC verbosity
      bool hepMCVerbosity_;
  };
} 

#endif
