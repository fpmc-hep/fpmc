#include "Fpmc.h"

#ifndef HEPMC_VERSION2
extern HEPEVT hepevt_;
#endif

namespace fpmc
{
  Fpmc::Fpmc( double comEnergy, const char* card ) :
    herwigVerbosity_( 1 ), hepMCVerbosity_( true ), maxEventsToPrint_( 2 ),
    comEnergy_( comEnergy ),
    params_( FpmcParameters::parseCard( card ) ),
    hadronize_( true ), debug_( false ), dbg_( std::cout )
  {
    params_.dump();
  }

  Fpmc::~Fpmc()
  {}

  void
  Fpmc::write( const char* out )
  {
#ifdef HEPMC_VERSION2
    std::ofstream output( out );
    hepMCEvt_->write( output );
    output.close();
#else // HepMC v>=3
    HepMC::WriterAscii output( out );
    output.write_event( *hepMCEvt_ );
    output.close();
#endif
  }

  void
  Fpmc::begin()
  {
    //--- use random seeds from datacard
    long seed0 = ( params_.has( "nrn1" ) ) ? params_.getLong( "nrn1" ) : -1L;
    long seed1 = ( params_.has( "nrn2" ) ) ? params_.getLong( "nrn2" ) : -1L;
    if ( seed0<0 || seed1<0 ) {
      auto randomEngine = std::make_shared<CLHEP::HepJamesRandom>();

      if ( seed0<0 ) seed0 = CLHEP::RandFlat::shoot( randomEngine.get(), 1L, 10000L );
      if ( seed1<0 ) seed1 = CLHEP::RandFlat::shoot( randomEngine.get(), 1L, 10000L );

      if ( debug_ ) dbg_ << "[FPMC Wrapper] SEEDS: " << seed0 << ", " << seed1 << std::endl;
    }

    //--- print the FPMC greetings
    fpmc_welcome();

    if ( debug_ ) {
      dbg_ << "[FPMC Wrapper] UTYPEPR = " << params_.getString( "typepr" ) << std::endl
           << "               UTYPINT = " << params_.getString( "typint" ) << std::endl
           << "               UTMASS  = " << params_.getFloat( "tmass" ) << std::endl;
    }

    dbg_ << "[FPMC Wrapper] Initializing HERWIG/FPMC" << std::endl;

    //--- call hwudat to set up HERWIG block data
    //hwudat();

    if ( params_.has( "part1" ) ) params_.getString( "part1" ).copy( hwbmch_.PART1, 8 );
    if ( params_.has( "part2" ) ) params_.getString( "part2" ).copy( hwbmch_.PART2, 8 );

    hwproc_.PBEAM1 = comEnergy_/2.;
    hwproc_.PBEAM2 = comEnergy_/2.;

    if ( params_.has( "typepr" ) ) params_.getString( "typepr" ).copy( prtype_.TYPEPR, 3 );
    if ( params_.has( "typint" ) ) params_.getString( "typint" ).copy( prtype_.TYPINT, 3 );

    if ( params_.has( "iproc" ) ) hwproc_.IPROC = params_.getInt( "iproc" );

    hadronize_ = strcmp( params_.getString( "hadr" ).c_str(), "Y" )==0;  
    dbg_ << "[FPMC Wrapper] Run hadronization/showering: " << hadronize_ << std::endl;
  
    if ( debug_ ) {
      dbg_ << "[FPMC Wrapper] PART1  = '" << hwbmch_.PART1 << "'" << std::endl
           << "[FPMC Wrapper] PART2  = '" << hwbmch_.PART2 << "'" << std::endl
           << "[FPMC Wrapper] TYPEPR = " << prtype_.TYPEPR << std::endl
           << "[FPMC Wrapper] TYPINT = " << prtype_.TYPINT << std::endl
           << "[FPMC Wrapper] IPROC  = " << hwproc_.IPROC << std::endl;
    }
    hwigin_();

    params_.fetchHWPRAM( hwpram_ );
    hwpram_.IPRINT = herwigVerbosity_;
    hwpram_.EFFMIN = 1.e-6;
    hwpram_.MODPDF[1] = -111;
    hwpram_.MODPDF[0] = -111;
    hwpram_.LWSUD = 0; // don't write Sudakov form factors

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
    chidefpmc_.CHIDeS = comEnergy_*comEnergy_;
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

    params_.writeCard( "lastrun.card" );
  }

  void
  Fpmc::end()
  {}

  bool
  Fpmc::run()
  {
#ifndef HEPMC_VERSION2
    HepMC::HEPEVT_Wrapper::set_hepevt_address( ( char* )&hepevt_ );
    HepMC::HEPEVT_Wrapper::zero_everything();
#endif
    //--- call herwig routines to create HEPEVT

    hwuine(); // initialize event

    hwepro(); // generate hard subprocess

    if ( hadronize_ ) {
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

    if ( hwevnt_.IERROR ) return false;
#ifdef HEPMC_VERSION2
    hepMCEvt_ = std::make_shared<HepMC::GenEvent>( *conv_.read_next_event() );
#else
    std::cout << "---> " << HepMC::HEPEVT_Wrapper::number_entries() << std::endl;
    //HepMC::HEPEVT_Wrapper::set_number_entries( 4000 );
    //HepMC::HEPEVT_Wrapper::set_sizeof_real( 8 );
    HepMC::HEPEVT_Wrapper::print_hepevt();
    HepMC::HEPEVT_Wrapper::HEPEVT_to_GenEvent( hepMCEvt_.get() );
#endif

    hepMCEvt_->set_event_number( event_-1 );
#ifdef HEPMC_VERSION2
    hepMCEvt_->set_signal_process_id( hwproc_.IPROC );
    hepMCEvt_->set_event_scale( -1. );
#endif
    hepMCEvt_->weights().push_back( hwevnt_.EVWGT );

#ifdef HEPMC_VERSION2
    HepMC::PdfInfo pdfInfo;
    pdfInfo.set_x1( hwhard_.XX[0] );
    pdfInfo.set_x2( hwhard_.XX[1] );
    pdfInfo.set_scalePDF( hwhard_.EMSCA );
    hepMCEvt_->set_pdf_info( pdfInfo );
#else
    auto pdf_info = std::make_shared<HepMC::GenPdfInfo>();
    pdf_info->set( 0, 0, hwhard_.XX[0], hwhard_.XX[1], hwhard_.EMSCA, 0, 0 );
    hepMCEvt_->set_pdf_info( pdf_info );
#endif

    /*HepMC::GenParticle* incomingParton = NULL;
    HepMC::GenParticle* targetParton = NULL;
    // find incoming parton (first entry with IST=121)
    for(HepMC::GenEvent::particle_const_iterator it = event()->particles_begin(); (it != event()->particles_end() && incomingParton==NULL); it++)
      if((*it)->status()==121) incomingParton = (*it);
  
    // find target parton (first entry with IST=122)
    for(HepMC::GenEvent::particle_const_iterator it = event()->particles_begin(); (it != event()->particles_end() && targetParton==NULL); it++)
      if((*it)->status()==122) targetParton = (*it);*/

#ifdef HEPMC_VERSION2
    //******** Verbosity ********
    if ( event_<=maxEventsToPrint_ && hepMCVerbosity_ ) {
      // Prints HepMC event
      dbg_ << "\n----------------------" << std::endl
           << "Event process id = " << hepMCEvt_->signal_process_id() << std::endl;
      hepMCEvt_->print();
    }
#endif

    //--- increment the events counter
    ++event_;

    return true;
  }
}
