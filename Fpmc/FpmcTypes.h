#ifndef Fpmc_FpmcTypes_h
#define Fpmc_FpmcTypes_h

namespace fpmc
{
  enum Flux
  {
    // Exclusive QCD fluxes
    PomeronFlux = 9,
    ReggeonFlux = 10,
    KMRExclusive = 16,
    CHIDeExclusive = 17,
    PomeronReggeon = 19,
    PomeronReggeonPP = 21,
    // QED fluxes
    PhotonPhotonCahnJackson = 12,
    PhotonPhotonAA = 13,
    PhotonPhotonPPPapageorgiou = 14,
    PhotonPhotonBudnevCoherent = 15,
    PhotonPhotonPA = 23,
    // Combined fluxes
    PhotonPomeronPP = 20,
    PomeronPhotonPP = 22,
    PhotonPomeronPA = 26,
    PomeronPhotonPA = 25
  };

  enum ProcessType { InvalidProcess = -1, ExclusiveProcess = 0, InclusiveProcess = 1 };
  enum InteractionType { InvalidInteraction = -1, QED = 0, QCD = 1 };
}

#endif
