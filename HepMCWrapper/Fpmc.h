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
#include "HepMC/HerwigWrapper.h"
#include "HepMC/HEPEVT_Wrapper.h"
#include "HepMC/Version.h"

#ifdef HEPMC_VERSION_CODE // HepMC v>=3
#define HEPMC_VERSION3
#include "HepMC/WriterAscii.h"
#else // HepMC v2
#include "HepMC/IO_HERWIG.h"
#endif

#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"

#include "fostream.h"
#include "herwig.h"
#include "fpmc.h"

#include <vector>
#include <string>
#include <ostream>
#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <cstring>

namespace fpmc
{
  class Fpmc{
  public:
    Fpmc(double, long int, std::vector<std::string> const&);
    ~Fpmc();

    void begin();
    bool run();
    void end();
    void write( const char* );

    const HepMC::GenEvent* event() const { return hepMCEvt_.get(); } 

  private:
#ifndef HEPMC_VERSION3
    HepMC::IO_HERWIG conv_;
#else
    HepMC::HEPEVT_Wrapper conv_;
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
    // Not used temporarily (take from datacard)
    long int seed_;

    std::vector<std::string> params_;

    bool hadronize_;
    bool debug_;

    CLHEP::HepRandomEngine* randomEngine_;
    CLHEP::RandFlat*        randomGenerator_; 
  };
} 

#endif
