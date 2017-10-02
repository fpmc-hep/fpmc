#include "Fpmc/Fpmc.h"

#include <iostream>

using namespace std;

int main( int argc, char* argv[] )
{
  fpmc::Fpmc generator( argv[1] );

  generator.parameters().dump();

  hwevnt_t ev;
  for( unsigned int evt = 0; evt < 10; ++evt ) {
    cout << "[FPMC] Processing event " << ( evt+1 ) << endl;
    generator.next( ev );
  }
  return 0;
}
