First search for the process of interest in the `Datacards` directory.
If the process you are looking for is not included there, look at the FPMC documentation or contact the developers.

To run FPMC, first:
```sh
make clean
make
```

For example, for the process &gamma;&gamma; &rarr; Spin 0 neutral resonance &varphi; &rarr; &gamma;&gamma; (With intact protons), use the command:
```sh
./fpmc <Datacards/dataQED_AASpin0EvenResonances
```

You can check some basic parameters of the collision and the coupling by editing the datacard, for example:

- the acceptance for the protons,
- the coupling value(s),
- the centre of mass energy,
- the output filename (default: `tmpntuple.ntp`)`
- etc.

You will find `tmpntuple.ntp` file as an output.
You can convert this file to a ROOT file with the `h2root` command:
```sh
h2root tmpntuple.ntp yourfilename.root
```
This ROOT file has all the basic kinematical quantities of the individual photons and protons.
You can then run an analysis code for this ROOT file.
In case you want to do a quick check instead (Look at the hard subprocess, parton showering, hadronization, etc. event by event) you can write the output to a `.dat` file with the command
```sh
./fpmc <Datacards/dataQED_AASpin0EvenResonances>yourfilename.dat
```

Note:

* If you are working on the lxplus nodes, use `setup_lxplus.sh` to set up the proper CERNLIB environment.

