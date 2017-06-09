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
#include "fpmc.h"

#include "FpmcParameters.h"

#include <string>
#include <sstream>
#include <iostream>

namespace fpmc
{
  class Fpmc
  {
    public:
      Fpmc( const char* card );
      Fpmc( const FpmcParameters& params );
      ~Fpmc();

      bool next( hwevnt_t& evt );

    protected:
      void init() const;
      void initHerwig() const;
 
      /// HERWIG verbosity
      unsigned int herwigVerbosity_;
      /// Events to print if verbosity
      unsigned int maxEventsToPrint_;

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
