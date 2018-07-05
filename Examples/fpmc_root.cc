#include "Fpmc/Fpmc.h"
#include "Fpmc/ArgsParser.h"

#include <TFile.h>
#include <TTree.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>

#include <iostream>

using namespace std;

extern "C"
{
  //COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
  //   & JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
  extern struct {
    int nevhep, nhep, isthep[nmxhep], idhep[nmxhep], jmohep[nmxhep][2], jdahep[nmxhep][2];
    double phep[nmxhep][5], vhep[nmxhep][4];
  } hepevt_;
}

int main( int argc, char* argv[] )
{
  const fpmc::ArgsParser args( argc, argv, { "cfg", "nevents" }, { "fileout", "comenergy" } );
  const auto required_params = args.required_parameters(),
             optional_params = args.optional_parameters();

  fpmc::Fpmc gen( required_params.at( "cfg" ).c_str() );

  // Optional parameters
  const string outputFileName = ( optional_params.count( "fileout" ) > 0 )
    ? optional_params.at( "fileout" )
    : "fpmc.root";

  cout << "=========================================================" << endl
       << " FPMC (ROOT) successfully initialised" << endl
       << "=========================================================" << endl;

  if ( optional_params.count( "comenergy" ) > 0 )
    gen.parameters().setSqrtS( stof( optional_params.at( "comenergy" ) ) ); // in GeV

  gen.parameters().dump();

  unique_ptr<TFile> file( TFile::Open( outputFileName.c_str(), "recreate" ) );
  //file->SetCompressionLevel( 1 );
  TTree tr( "events", "A tree containing FPMC events" );
  vector<unsigned short> status;
  tr.Branch( "status", &status );
  vector<short> pdgid;
  tr.Branch( "pdgid", &pdgid );
  vector<ROOT::Math::XYZTVector> momentum;
  vector<ROOT::Math::XYZTVector>* pMom = &momentum;
  tr.Branch( "momentum", "std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >", &pMom );
  vector<unsigned short> mother1, mother2;
  tr.Branch( "mother1", &mother1 );
  tr.Branch( "mother2", &mother2 );

  for ( unsigned long evt = 0; evt < stoul( required_params.at( "nevents" ) ); ++evt ) {
    if ( ( evt+1 ) % 10000 == 0 )
      cout << "[FPMC] Processing event " << ( evt+1 ) << endl;
    gen.next();
    if ( hwevnt_.IERROR )
      continue;
    pMom->clear();
    pMom->reserve( hepevt_.nhep );
    status.clear(); status.reserve( hepevt_.nhep );
    pdgid.clear(); pdgid.reserve( hepevt_.nhep );
    mother1.clear(); mother1.reserve( hepevt_.nhep );
    mother2.clear(); mother2.reserve( hepevt_.nhep );
    for ( int i = 0; i < hepevt_.nhep; ++i ) {
      status.emplace_back( hepevt_.isthep[i] );
      pdgid.emplace_back( hepevt_.idhep[i] );
      pMom->emplace_back( hepevt_.phep[i][0], hepevt_.phep[i][1], hepevt_.phep[i][2], hepevt_.phep[i][3] );
      mother1.emplace_back( hepevt_.jmohep[i][0] );
      mother2.emplace_back( hepevt_.jmohep[i][1] );
    }
    tr.Fill();
  }
  file->Write();

  return 0;
}

