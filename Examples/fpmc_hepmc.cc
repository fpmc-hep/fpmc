#include "Fpmc/HepMCWrapper.h"
#include "Fpmc/ArgsParser.h"

#include <iostream>

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
  string datacard = required_params.at( "cfg" );
  unsigned int maxEvents = stoi( required_params.at( "nevents" ) );
  double comEnergy = stof( required_params.at( "comenergy" ) );

  // Optional parameters
  string outputFileName = "fpmc.hepmc";
  if ( optional_params.count( "fileout" ) > 0 )
    outputFileName = optional_params.at( "fileout" );

  cout << "=========================================================" << endl
       << "FPMC (Wrapper) will initialize with parameters: " << endl
       << "  Datacard:    " << datacard << endl
       << "  N events:    " << maxEvents << endl
       << "  COM energy:  " << comEnergy << endl
       << "  Output file: " << outputFileName << endl
       << "=========================================================" << endl;

  fpmc::HepMCWrapper generator( comEnergy, datacard.c_str() );

#ifdef HEPMC_VERSION2
  std::ofstream output_file( outputFileName, ios::out );
  HepMC::IO_GenEvent output( output_file );
#else
  HepMC::WriterAscii output( outputFileName );
#endif
  for ( unsigned int evt = 0; evt < maxEvents; ++evt ) {
    cout << "[FPMC Wrapper] Processing event " << ( evt+1 ) << endl;
#ifdef HEPMC_VERSION2
    output.write_event( generator.event() );
#else
    output.write_event( *generator.event() );
#endif
  }
  //cout << "Cross section: " << generator.crossSection() << " +/- " << generator.crossSection() / sqrt( maxEvents ) << " pb"<< endl;
  return 0;
}
