#include "Fpmc.h"
#include "dummy_hwaend.h"

int main()
{
  fpmc::Fpmc gen;
  gen.parameters().setProcessId( 16008 ); // QED yy->mu,mu
  gen.parameters().setInitialFlux( 15 );
  gen.initialise();

  hwevnt_t event;
  unsigned int i=0, num_events = 100;

  do {
    if ( gen.next( event ) ) i++;
  } while ( i<num_events );
  std::cout << "cross section after " << num_events << " events generated: " << gen.crossSection() << " pb" << std::endl;
}
