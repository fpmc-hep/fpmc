#include "Fpmc/ArgsParser.h"

int main( int argc, char* argv[] )
{
  fpmc::ArgsParser args( argc, argv, { "cfg", "nevents", "comenergy" }, { "fileout" } );

  return 0;
}
