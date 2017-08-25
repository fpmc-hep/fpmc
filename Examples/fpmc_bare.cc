#include "Fpmc/Fpmc.h"

#include <iostream>

using namespace std;

int main( int argc, char* argv[] )
{
  //fpmc::Fpmc gen( argv[1] );
  fpmc::Fpmc gen;
  gen.parameters().setSqrtS( 13.e3 ); // in GeV
  gen.parameters().setProcessId( 16008 ); // QED yy->mu,mu
  gen.parameters().setProcessType( fpmc::ExclusiveProcess );
  gen.parameters().setInteractionType( fpmc::QED );
  gen.parameters().setIntermediateFlux( fpmc::PhotonPhotonBudnevCoherent );
  gen.parameters().setPtRange( 10. ); // no upper limit

  gen.parameters().dump();

  for( unsigned int evt = 0; evt < 10e5; ++evt ) {
    cout << "[FPMC] Processing event " << ( evt+1 ) << endl;
    gen.next();
  }
  return 0;
}
