#ifndef FpmcInterface_FpmcParameters_h
#define FpmcInterface_FpmcParameters_h

#include "FpmcTypes.h"

#include "Fpmc/interface/herwig.h"
#include "Fpmc/interface/fpmc.h"

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

      void validate();

      static FpmcParameters parseCard( const char* filename );
      void writeCard( const char* filename ) const;

      std::map<std::string,std::string> map() const { return *this; }

      bool hadronise() const { return strcmp( getString( "hadr" ).c_str(), "Y" )==0; }

      //----- unique ID of the process to generate

      void setProcessId( unsigned int iproc ) { add( "iproc", iproc ); }
      unsigned int processId() const { return getInt( "iproc" ); }

      //----- type of process involved

      void setProcessType( const char* typepr ) { add( "typepr", typepr ); }
      void setProcessType( const ProcessType& typepr ) {
        switch ( typepr ) {
          case ExclusiveProcess: add( "typepr", "EXC" ); break;
          case InclusiveProcess: add( "typepr", "INC" ); break;
          default: break;
        }
      }
      ProcessType processType() const {
        const std::string typepr = getString( "typepr" );
        if ( typepr.compare( "EXC" )==0 ) return ExclusiveProcess;
        if ( typepr.compare( "INC" )==0 ) return InclusiveProcess;
        return InvalidProcess;
      }

      //----- type of interaction involved

      void setInteractionType( const char* typint ) { add( "typint", typint ); }
      void setInteractionType( const InteractionType& typint ) {
        switch ( typint ) {
          case QED: add( "typint", "QED" ); break;
          case QCD: add( "typint", "QCD" ); break;
          default: break;
        }
      }
      InteractionType interactionType() const {
        const std::string typint = getString( "typint" );
        if ( typint.compare( "QED" )==0 ) return QED;
        if ( typint.compare( "QCD" )==0 ) return QCD;
        return InvalidInteraction;
      }

      //----- type of intermediate particles flux

      void setIntermediateFlux( const Flux& nflux ) { add( "nflux", nflux ); }
      Flux intermediateFlux() const { return static_cast<Flux>( getInt( "nflux" ) ); }

      //----- centre of mass energy of the initial system

      void setSqrtS( double sqrts ) { add( "ecms", sqrts ); }
      double sqrtS() const { return getFloat( "ecms" ); }

      //----- events kinematics

      void setPtRange( double ptmin, double ptmax=0. ) {
        add( "ptmin", ptmin );
        if ( ptmax>0. ) add( "ptmax", ptmax );
      }
      double ptMin() const { return getFloat( "ptmin" ); }
      double ptMax() const { return getFloat( "ptmax" ); }

      void setEtaRange( double etamin, double etamax ) {
        add( "yjmin", etamin );
        add( "yjmax", etamax );
      }
      double etaMin() const { return getFloat( "yjmin" ); }
      double etaMax() const { return getFloat( "yjmax" ); }

      //----- PDF fits

      void setPDFfits( const PDFfits& fits ) { add( "ifit", fits ); }
      PDFfits pdfFits() const { return static_cast<PDFfits>( getInt( "ifit" ) ); }

      //----- full common blocks population

      void fetchHWPROC( hwproc_t& ) const;
      void fetchHWBMCH( hwbmch_t& ) const;
      void fetchHWPRAM( hwpram_t& ) const;
      void fetchHWPROP( hwprop_t& ) const;
      void fetchHWHARD( hwhard_t& ) const;
      void fetchHWPRCH( hwprch_t& ) const;
      void fetchPRTYPE( prtype_t& ) const;
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
