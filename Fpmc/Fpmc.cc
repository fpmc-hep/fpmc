#include "Fpmc.h"

namespace fpmc
{
  Fpmc::Fpmc() :
    herwigVerbosity_( 1 ), maxEventsToPrint_( 2 ),
    initialised_( false ),
    debug_( false ), dbg_( std::cout )
  {}

  Fpmc::Fpmc( const char* card ) :
    herwigVerbosity_( 1 ), maxEventsToPrint_( 2 ),
    initialised_( false ),
    params_( FpmcParameters::parseCard( card ) ),
    debug_( false ), dbg_( std::cout )
  {}

  Fpmc::Fpmc( const FpmcParameters& params ) :
    herwigVerbosity_( 1 ), maxEventsToPrint_( 2 ),
    initialised_( false ),
    params_( params ),
    debug_( false ), dbg_( std::cout )
  {}

  Fpmc::~Fpmc()
  {}

  void
  Fpmc::initialise()
  {
    //--- print the FPMC greetings
    fpmc_welcome();

    if ( debug_ ) {
      dbg_ << "UTYPEPR = " << params_.getString( "typepr" ) << std::endl
           << "UTYPINT = " << params_.getString( "typint" ) << std::endl
           << "UTMASS  = " << params_.getFloat( "tmass" ) << std::endl;
    }
    initHerwig();
    params_.writeCard( "lastrun.card" );
  }

  void
  Fpmc::initHerwig()
  {
    if ( initialised_ ) dbg_ << "WARNING: Herwig already initialised for this run!" << std::endl;

    dbg_ << "Initializing HERWIG/FPMC" << std::endl;

    //--- call hwudat to set up HERWIG block data
    //hwudat();

    if ( params_.has( "part1" ) ) params_.getString( "part1" ).copy( hwbmch_.PART1, 8 );
    if ( params_.has( "part2" ) ) params_.getString( "part2" ).copy( hwbmch_.PART2, 8 );

    hwproc_.PBEAM1 = hwproc_.PBEAM2 = params_.sqrtS()/2.;

    if ( params_.has( "typepr" ) ) params_.getString( "typepr" ).copy( prtype_.TYPEPR, 3 );
    if ( params_.has( "typint" ) ) params_.getString( "typint" ).copy( prtype_.TYPINT, 3 );

    if ( params_.has( "iproc" ) ) hwproc_.IPROC = params_.processId();

    if ( debug_ ) {
      dbg_ << "PART1  = '" << hwbmch_.PART1 << "'" << std::endl
           << "PART2  = '" << hwbmch_.PART2 << "'" << std::endl
           << "TYPEPR = " << prtype_.TYPEPR << std::endl
           << "TYPINT = " << prtype_.TYPINT << std::endl
           << "IPROC  = " << hwproc_.IPROC << std::endl;
    }
    hwigin_();

    params_.fetchHWPRAM( hwpram_ );
    hwpram_.IPRINT = herwigVerbosity_;
    hwpram_.EFFMIN = 1.e-6;
    hwpram_.MODPDF[1] = -111;
    hwpram_.MODPDF[0] = -111;
    hwpram_.LWSUD = 0; // don't write Sudakov form factors

    //--- use random seeds from datacard
    long seed0 = ( params_.has( "nrn1" ) ) ? params_.getLong( "nrn1" ) : -1L;
    long seed1 = ( params_.has( "nrn2" ) ) ? params_.getLong( "nrn2" ) : -1L;
    if ( seed0<0 || seed1<0 ) {
      auto randomEngine = std::make_shared<CLHEP::HepJamesRandom>();

      if ( seed0<0 ) seed0 = CLHEP::RandFlat::shoot( randomEngine.get(), 1L, 10000L );
      if ( seed1<0 ) seed1 = CLHEP::RandFlat::shoot( randomEngine.get(), 1L, 10000L );

      if ( debug_ ) dbg_ << "SEEDS: " << seed0 << ", " << seed1 << std::endl;
    }

    hwevnt_.NRN[0] = seed0;
    hwevnt_.NRN[1] = seed1;

    hwevnt_.MAXER = 100000000; // O(inf)
    hwdspn_.LWDEC = 0; // don't write three/four body decays
    // (no fort.77 and fort.88 ...)

    // Init LHAPDF glue
    const std::string pdfset( "HWLHAPDF" );
    for ( unsigned int i=0; i<2; i++ ) {
      pdfset.copy( hwprch_.AUTPDF[i], 8 );
    }

    hwevnt_.MAXPR = maxEventsToPrint_;

    //--- feed the parameters from the steering card to the generator core

    params_.fetchHWPROP( hwprop_ );
    params_.fetchHWHARD( hwhard_ );

    params_.fetchXSECT( xsect_ );
    params_.fetchPDFS( pdfs_ );
    params_.fetchAAANOMAL( aaanomal_ );
    params_.fetchAAEXOTICAL( aaexotical_ );

    params_.fetchCHIDEFPMC( chidefpmc_ );
    chidefpmc_.CHIDeS = params_.sqrtS()*params_.sqrtS();
    params_.fetchKMR2FPMC( kmr2fpmc_ );

    params_.fetchION( ion_ );
    params_.fetchCHIDEENV( chideenv_ );

    //--- initialize model/pdf dependant parameters
    hwmodini();

    //--- compute parameter dependent constants
    hwuinc();

    //--- check POMWIG Settings + Initialisations for consistency
    hwchek();

    //--- call HWUSTA to make any particle stable
    int iopt = 1;
    int iwig = 0;
    char nwig[9] = "        ";

    int ipdg = 111;
    hwuidt( &iopt, &ipdg, &iwig, nwig );
    if ( ipdg ) hwusta(nwig, 1);
 
    //--- initialize elementary process
    hweini();

    //--- initialize event record fixing
    // this will replace the beam electrons by protons, radiated photons by pomerons/reggeons, etc.
    int init = 1;
    hwfxer( &init );

    initialised_ = true;
  }

  double
  Fpmc::crossSection() const
  {
    if ( !initialised_ ) {
      dbg_ << "WARNING: Herwig/FPMC not yet initialised! Invalid cross-section retrieved." << std::endl;
      return 0.;
    }
    return hwevnt_.EVWGT * 1.e3; // return in pb
  }

  bool
  Fpmc::next( hwevnt_t& event )
  {
    //--- call herwig routines to create HEPEVT

    hwuine(); // initialize event

    hwepro(); // generate hard subprocess

    if ( params_.hadronise() ) {
      hwbgen();	// parton cascades

      hwdhob();	// heavy quark decays
      hwcfor();	// cluster formation
      hwcdec();	// cluster decays

      hwdhad();	// unstable particle decays
      hwdhvy();	// heavy flavour decays
      hwmevt();	// soft underlying event		
    }

    hwufne(); // finalize event

    int init = 0;
    hwfxer( &init ); // fix event record (i.e. restore correct intermediate states)

    event = hwevnt_;

    if ( hwevnt_.IERROR ) return false;

    //--- increment the events counter
    ++event_;

    return true;
  }
}
