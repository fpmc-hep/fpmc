#include "Fpmc/HepMCWrapper.h"
#include "Fpmc/ArgsParser.h"

#include <iostream>
#include <iomanip>

#ifdef HEPMC_VERSION2
#include "HepMC/IO_GenEvent.h"
#else
#include "HepMC/WriterAscii.h"
#endif

using namespace std;

int main( int argc, char* argv[] )
{
  fpmc::ArgsParser args( argc, argv, { "cfg", "nevents", "comenergy" }, { "fileout" } );
  const auto required_params = args.required_parameters(), optional_params = args.optional_parameters();

  //----------------------------
  // Required parameters 
  string datacard_ = required_params.at( "cfg" );
  unsigned int maxEvents_ = stoi( required_params.at( "nevents" ) );
  double comEnergy_ = stof( required_params.at( "comenergy" ) );

  // Optional parameters 
  string outputFileName_ = "fpmc.hepmc";
  if ( optional_params.count("fileout") != 0 ) outputFileName_ = optional_params.at( "fileout" );

  stringstream oss;
  oss  << "=========================================================" << endl
       << "FPMC (Wrapper) will initialize with parameters: " << endl
       << "  Datacard:    " << datacard_ << endl
       << "  N events:    " << maxEvents_ << endl
       << "  COM energy:  " << comEnergy_ << endl
       << "  Output file: " << outputFileName_ << endl
       << "=========================================================" << endl;
  cout << oss.str();

  fpmc::HepMCWrapper generator( comEnergy_, datacard_.c_str() );


#ifdef HEPMC_VERSION2
  HepMC::IO_GenEvent output( outputFileName_, ios::out );
#else
  HepMC::WriterAscii output( outputFileName_ );
#endif
  for( unsigned int evt = 0; evt < maxEvents_; ++evt ) {
    cout << "[FPMC Wrapper] Processing event " << ( evt+1 ) << endl;
#ifdef HEPMC_VERSION2
    output.write_event( generator.event() );
#else
    output.write_event( *generator.event() );
#endif
  }
  return 0;
}
