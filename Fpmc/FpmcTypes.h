#ifndef Fpmc_FpmcTypes_h
#define Fpmc_FpmcTypes_h

#include <iostream>

namespace fpmc
{
  enum Flux
  {
    // Exclusive QCD fluxes
    PomeronFlux = 9, ///< Cox-Forshaw (DPE)
    ReggeonFlux = 10, ///< Bialas-Landshoff (exc. DPE) or Boonekamp et al (inc. DPE)
    KMRExclusive = 16, ///<KMR flux (tables from L.Lonnblad or ExHume)
    CHIDeExclusive = 17,
    PomeronReggeon = 19,
    PomeronReggeonPP = 21,
    // QED fluxes
    PhotonPhotonCahnJackson = 12, ///< Cahn, Jackson (heavy ions, QED)
    PhotonPhotonAA = 13, ///< Drees, Ellis, Zeppenfled (heavy ions, QED)
    PhotonPhotonPPPapageorgiou = 14, ///< Papageorgiu (proton, QED)
    PhotonPhotonBudnevCoherent = 15, ///< Budnev flux (RECOMMENDED for proton, QED)
    PhotonPhotonPA = 23,
    // Combined fluxes
    PhotonPomeronPP = 20,
    PomeronPhotonPP = 22,
    PhotonPomeronPA = 26,
    PomeronPhotonPA = 25
  };
  /// PDF sets to use
  /// (according to hep/ph 0609291)
  enum PDFfits
  {
    //NLOfitH1 = 2, LOfitH1 = 5, ExtendedH1 = 8, //legacy
    H1 = 10, Zeus = 20, H1ZeusCombined = 30
  };

  enum ProcessType { InvalidProcess = -1, ExclusiveProcess = 0, InclusiveProcess = 1 };
  enum InteractionType { InvalidInteraction = -1, QED = 0, QCD = 1 };

  inline std::ostream& operator<<( std::ostream& os, const ProcessType& type ) {
    switch ( type ) {
      case InvalidProcess: return os << "INV";
      case ExclusiveProcess: return os << "EXC";
      case InclusiveProcess: return os << "INC";
    }
    return os;
  }
  inline std::ostream& operator<<( std::ostream& os, const InteractionType& type ) {
    switch ( type ) {
      case InvalidInteraction: return os << "INV";
      case QED: return os << "QED";
      case QCD: return os << "QCD";
    }
    return os;
  }
}

#endif
