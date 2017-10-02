#!/bin/sh

g++ test_standalone.cc -o test_standalone \
  -L../build/ -lFPMC \
  `gsl-config --libs` \
  `clhep-config --libs` \
  -I../Fpmc -I../Herwig
