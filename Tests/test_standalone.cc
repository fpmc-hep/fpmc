#include "Fpmc.h"
#include "dummy_hwaend.h"

int main()
{
  fpmc::Fpmc gen;

  //--- parameterise the run

  gen.parameters().setSqrtS( 13.e3 ); // in GeV
  gen.parameters().setProcessId( 16008 ); // QED yy->mu,mu
  gen.parameters().setProcessType( fpmc::ExclusiveProcess );
  gen.parameters().setInteractionType( fpmc::QED );
  gen.parameters().setIntermediateFlux( fpmc::PhotonPhotonBudnevCoherent );
  gen.parameters().setPtRange( 10. ); // no upper limit

  //--- initialise the run

  gen.initialise();

  //--- generate some events

  hwevnt_t event;
  unsigned int i=0, num_events = 100;

  do {
    if ( gen.next( event ) ) i++;
  } while ( i<num_events );
  std::cout << "cross section after " << num_events << " events generated: " << gen.crossSection() << " pb" << std::endl;
}
