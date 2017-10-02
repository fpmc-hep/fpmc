#ifndef Fpmc_Fpmc_h
#define Fpmc_Fpmc_h

/** \class Fpmc
 *
 * Main FPMC object
 *
 ***************************************/

#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"

#include "herwig.h"
#include "Fpmc/interface/fpmc.h"

#include "FpmcParameters.h"

#include <string>
#include <sstream>
#include <iostream>

namespace fpmc
{
  class Fpmc
  {
    public:
      Fpmc();
      Fpmc( const char* card );
      Fpmc( const FpmcParameters& params );
      ~Fpmc();

      FpmcParameters& parameters() { return params_; }

      void initialise();
      double crossSection() const;

      bool next();

    protected:
      void initHerwig();
 
      /// HERWIG verbosity
      unsigned int herwigVerbosity_;
      /// Events to print if verbosity
      unsigned int maxEventsToPrint_;

      bool initialised_;
      /// Number of events already generated
      unsigned int event_;

      /// List of parameters obtained from the steering card
      FpmcParameters params_;

      /// Enable/disable the extra printout
      bool debug_;
      std::ostream& dbg_;
  };
} 

#endif
