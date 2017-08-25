First search for the process of interest in the `Datacards` directory.
If the process you are looking for is not included there, look at the FPMC documentation or contact the developers.

To run FPMC, first:
```sh
mkdir build
cd build
cmake {PATH_TO_FPMC_SOURCES}
make
```

For example, for the process &gamma;&gamma; &rarr; Spin 0 neutral resonance &varphi; &rarr; &gamma;&gamma; (With intact protons), use the command:
```sh
./fpmc-hepmc --cfg Datacards/dataQED_AASpin0EvenResonances --comenergy 13000 --nevents 1000
```

You can check some basic parameters of the collision and the coupling by editing the datacard, for example:

- the acceptance for the protons,
- the coupling value(s),
- the centre of mass energy,
- the output filename (default: `tmpntuple.ntp`)`
- etc.

You will find `fpmc.hepmc` file as an output.

Note:

* If you are working on the lxplus nodes, use `source /cvmfs/sft.cern.ch/lcg/external/gcc/6.1.0/x86_64-slc6/setup.sh` to set up the proper `gcc` environment.

