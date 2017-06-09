#ifndef FpmcInterface_FpmcParameters_h
#define FpmcInterface_FpmcParameters_h

#include "herwig.h"
#include "fpmc.h"

#include <string>
#include <ostream>
#include <sstream>
#include <algorithm>
#include <map>
#include <regex>

#include <fstream>
#include <iostream>
#include <iomanip>

#include <type_traits>

namespace fpmc
{
  /** \class Parameters
   *
   */
  class FpmcParameters : private std::map<std::string,std::string>
  {
    public:
      FpmcParameters();
      FpmcParameters( const std::map<std::string,std::string>& params ) : std::map<std::string,std::string>( params ) {}
      ~FpmcParameters() {}

      static FpmcParameters parseCard( const char* filename );
      void writeCard( const char* filename ) const;

      std::map<std::string,std::string> map() const { return *this; }

      bool hadronise() const { return strcmp( getString( "hadr" ).c_str(), "Y" )==0; }

      void setProcessId( unsigned int iproc ) { add( "iproc", iproc ); }
      unsigned int processId() const { return getInt( "iproc" ); }

      void setInitialFlux( unsigned int nflux ) { add( "nflux", nflux ); }
      unsigned int initialFlux() const { return getInt( "nflux" ); }

      void setSqrtS( double sqrts ) { add( "ecms", sqrts ); }
      double sqrtS() const { return getFloat( "ecms" ); }
      
      void fetchHWPRAM( hwpram_t& ) const;
      void fetchHWPROP( hwprop_t& ) const;
      void fetchHWHARD( hwhard_t& ) const;
      void fetchXSECT( xsect_t& ) const;
      void fetchPDFS( pdfs_t& ) const;
      void fetchAAANOMAL( aaanomal_t& ) const;
      void fetchAAEXOTICAL( aaexotical_t& ) const;
      void fetchCHIDEFPMC( chidefpmc_t& ) const;
      void fetchKMR2FPMC( kmr2fpmc_t& ) const;
      void fetchION( ion_t& ) const;
      void fetchCYFFLONG1( cyfflong1_t& ) const;
      void fetchCHIDEPATH( chidepath_t& ) const;
      void fetchCHIDEENV( chideenv_t& ) const;
      
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
        for ( const auto& pair : map() ) { os << "[" << std::setw( 12 ) << pair.first << "]\t\"" << pair.second << "\"" << std::endl; }
      }

      const std::string getString( const std::string& key ) const {
        const auto it = find( key );
        if ( it==end() ) return "";
        return it->second;
      }
      int getInt( const std::string& key ) const {
        if ( !has( key ) ) return 0;
        return std::stoi( getString( key ) );
      }
      long getLong( const std::string& key ) const {
        if ( !has( key ) ) return 0;
        return std::stol( getString( key ) );
      }
      double getFloat( const std::string& key ) const {
        if ( !has( key ) ) return 0.;
        return std::stod( getString( key ) );
      }
  };
} 

#endif
