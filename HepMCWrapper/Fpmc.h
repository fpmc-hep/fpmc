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
#include <map>
#include <regex>

namespace fpmc
{
  class Fpmc
  {
    public:
      class Parameters : private std::map<std::string,std::string>
      {
        public:
          Parameters() : std::map<std::string,std::string>() {}
          Parameters( const std::map<std::string,std::string>& map ) : std::map<std::string,std::string>( map ) {}
          ~Parameters() {}

	  std::map<std::string,std::string> map() const { return *this; }

          void add( const std::string& key, const std::string& value ) {
	    auto pair = find( key );
	    if ( pair==end() ) insert( std::make_pair( key, value ) );
	    else pair->second = value;
	  }
	  void add( const std::string& key, float value ) { std::ostringstream os; os << value; add( key, os.str() ); }
	  //void add( const std::string& key, int value ) { std::ostringstream os; os << value; add( key, os.str() ); }
          bool has( const std::string& key ) const { return find( key )!=end(); }

	  void dump( std::ostream& os=std::cout ) const {
	    os << size() << " key(s) in the parameters list" << std::endl;
	    for ( const auto& pair : map() ) { os << "[" << pair.first << "] \"" << pair.second << "\"" << std::endl; }
	  }

          const std::string getString( const std::string& key ) const {
            const auto it = find( key );
            if ( it==end() ) return "";
            return it->second;
          }
          int getInt( const std::string& key ) const {
            if ( !has( key ) ) return 0;
            return atoi( getString( key ).c_str() );
          }
          long getLong( const std::string& key ) const {
            if ( !has( key ) ) return 0;
            return atol( getString( key ).c_str() );
          }
          unsigned int getFloat( const std::string& key ) const {
            if ( !has( key ) ) return 0.;
            return atof( getString( key ).c_str() );
          }
      };

      Fpmc( double, long int, const char* );
      ~Fpmc();

      static void parseInputCard( const char*, Parameters& );
      static void writeInputCard( const char*, const Parameters& );

      void begin();
      bool run();
      void end();
      void write( const char* );

      const HepMC::GenEvent* event() const { return hepMCEvt_.get(); } 

    private:
      void initialiseParams();
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

      Parameters params_;

      bool hadronize_;
      bool debug_;
  };
} 

#endif
