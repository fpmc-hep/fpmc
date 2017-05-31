#include "Fpmc.h"

namespace fpmc
{
  Fpmc::Fpmc( double comEnergy, long int seed, const char* card ) :
    herwigVerbosity_( 1 ), hepMCVerbosity_( true ), maxEventsToPrint_( 2 ),
    comEnergy_( comEnergy ),
    hadronize_( true ), debug_( false )
  {
    initialiseParams();
    parseInputCard( card, params_ );
  }

  Fpmc::~Fpmc()
  {}

  void
  Fpmc::write( const char* out )
  {
#ifndef HEPMC_VERSION3
    hepMCEvt_->write( out );
#else // HepMC v>=3
    HepMC::WriterAscii output( out );
    output.write_event( *hepMCEvt_ );
    output.close();
#endif
  }

  void
  Fpmc::begin()
  {
    // Use random seeds from datacard
    // ===============================
    // Using CLHEP engines
    long seed0 = ( params_.has( "nrn1" ) ) ? params_.getLong( "nrn1" ) : -1L;
    long seed1 = ( params_.has( "nrn2" ) ) ? params_.getLong( "nrn2" ) : -1L;
    if ( seed0<0 || seed1<0 ) {
      auto randomEngine = std::make_shared<CLHEP::HepJamesRandom>();

      if ( seed0<0 ) seed0 = CLHEP::RandFlat::shoot( randomEngine.get(), 1L, 10000L );
      if ( seed1<0 ) seed1 = CLHEP::RandFlat::shoot( randomEngine.get(), 1L, 10000L );

      std::cout << "[FPMC Wrapper] SEEDS: " << seed0 << ", " << seed1 << std::endl;
    }
    // ===============================

    fpmc_welcome();

    if ( debug_ ) {
      std::cout << "[FPMC Wrapper] UTYPEPR = " << params_.getString( "typepr" ) << std::endl
		<< "               UTYPINT = " << params_.getString( "typint" ) << std::endl
		<< "               UTMASS  = " << params_.getFloat( "tmass" ) << std::endl;
    }

    std::cout << "[FPMC Wrapper] Initializing HERWIG/FPMC" << std::endl;

    // Call hwudat to set up HERWIG block data
    //hwudat();

    if ( params_.has( "part1" ) ) params_.getString( "part1" ).copy( hwbmch.PART1, 8 );
    if ( params_.has( "part2" ) ) params_.getString( "part2" ).copy( hwbmch.PART2, 8 );

    hwproc.PBEAM1 = comEnergy_/2.;
    hwproc.PBEAM2 = comEnergy_/2.;

    if ( params_.has( "typepr" ) ) params_.getString( "typepr" ).copy( prtype.TYPEPR, 3 );
    if ( params_.has( "typint" ) ) params_.getString( "typint" ).copy( prtype.TYPINT, 3 );

    if ( params_.has( "iproc" ) ) hwproc.IPROC = params_.getInt( "iproc" );

    //ANSWER=UHADR 
    //
    hadronize_ = strcmp( params_.getString( "hadr" ).c_str(), "Y" )==0;  
    std::cout << "[FPMC Wrapper] Run hadronization/showering: " << hadronize_ << std::endl; 
  
    if ( debug_ ) {
      std::cout << "[FPMC Wrapper] PART1  = '" << hwbmch.PART1 << "'" << std::endl
                << "[FPMC Wrapper] PART2  = '" << hwbmch.PART2 << "'" << std::endl
                << "[FPMC Wrapper] TYPEPR = " << prtype.TYPEPR << std::endl
                << "[FPMC Wrapper] TYPINT = " << prtype.TYPINT << std::endl
                << "[FPMC Wrapper] IPROC  = " << hwproc.IPROC << std::endl;
    }
    //CALL HWIGIN
    //
    hwigin();

    hwevnt.NRN[0] = seed0;
    hwevnt.NRN[1] = seed1;
    hwpram.EFFMIN = 1.e-6;
    hwevnt.MAXER = 100000000; // O(inf)
    hwpram.LWSUD = 0;         // don't write Sudakov form factors
    hwdspn.LWDEC = 0;         // don't write three/four body decays
    // (no fort.77 and fort.88 ...)

    // Init LHAPDF glue
    std::memset(hwprch.AUTPDF, ' ', sizeof(hwprch.AUTPDF));
    for ( unsigned int i = 0; i < 2; i++ ) {
      hwpram.MODPDF[i] = -111;
      std::memcpy(hwprch.AUTPDF[i], "HWLHAPDF", 8);
    }

    hwevnt.MAXPR = maxEventsToPrint_;
    hwpram.IPRINT = herwigVerbosity_;

    if ( params_.has( "modpdf1" ) ) hwpram.MODPDF[0] = params_.getInt( "modpdf1" );
    if ( params_.has( "modpdf2" ) ) hwpram.MODPDF[1] = params_.getInt( "modpdf2" );

    if ( params_.has( "hmass" ) ) hwprop.RMASS[201] = params_.getFloat( "hmass" ); // higgs mass
    if ( params_.has( "tmass" ) ) hwprop.RMASS[6]   = params_.getFloat( "tmass" ); // top mass
    if ( params_.has( "wmass" ) ) hwprop.RMASS[198] = params_.getFloat( "wmass" ); // W mass
    if ( params_.has( "mst1" ) ) hwprop.RMASS[406] = params_.getFloat( "mst1" ); // stop1 mass
    if ( params_.has( "msb1" ) ) hwprop.RMASS[405] = params_.getFloat( "msb1" ); // stop2 mass

    if ( params_.has( "q2wwmn" ) ) hwhard.Q2WWMN = params_.getFloat( "q2wwmn" );
    if ( params_.has( "q2wwmx" ) ) hwhard.Q2WWMX = params_.getFloat( "q2wwmx" );
    if ( params_.has( "ywwmin" ) ) hwhard.YWWMIN = params_.getFloat( "ywwmin" );
    if ( params_.has( "ywwmax" ) ) hwhard.YWWMAX = params_.getFloat( "ywwmax" );
    if ( params_.has( "yjmin" ) ) hwhard.YJMIN = params_.getFloat( "yjmin" );
    if ( params_.has( "yjmax" ) ) hwhard.YJMAX = params_.getFloat( "yjmax" );
    if ( params_.has( "ptmin" ) ) hwhard.PTMIN = params_.getFloat( "ptmin" );
    if ( params_.has( "ptmax" ) ) hwhard.PTMAX = params_.getFloat( "ptmax" );
    if ( params_.has( "emmin" ) ) hwhard.EMMIN = params_.getFloat( "emmin" );

    if ( params_.has( "nflux" ) ) xsect.NFLUX = params_.getInt( "nflux" );
    if ( params_.has( "ifit" ) ) pdfs.IFITPDF = params_.getInt( "ifit" );

    if ( params_.has( "aaanom" ) ) aaanomal.AAANOM = params_.getInt( "aaanom" );
    if ( params_.has( "dkappa" ) ) aaanomal.D_KAPPA = params_.getFloat( "dkappa" );
    if ( params_.has( "dlambda" ) ) aaanomal.LAMBDA = params_.getFloat( "dlambda" );
    if ( params_.has( "a0w" ) ) aaanomal.A0W = params_.getFloat( "a0w" );
    if ( params_.has( "acw" ) ) aaanomal.ACW = params_.getFloat( "acw" );
    if ( params_.has( "a0z" ) ) aaanomal.A0Z = params_.getFloat( "a0z" );
    if ( params_.has( "acz" ) ) aaanomal.ACZ = params_.getFloat( "acz" );
    if ( params_.has( "a1a" ) ) aaanomal.A1A = params_.getFloat( "a1a" );
    if ( params_.has( "a2a" ) ) aaanomal.A2A = params_.getFloat( "a2a" );
    if ( params_.has( "anomcutoff" ) ) aaanomal.ANOMCUTOFF = params_.getInt( "anomcutoff" );

    if ( params_.has( "aaexotic" ) ) aaexotical.AAEXOTIC = params_.getInt( "aaexotic" );
    if ( params_.has( "aam" ) ) aaexotical.AAM = params_.getFloat( "aam" );
    if ( params_.has( "aaq" ) ) aaexotical.AAQ = params_.getFloat( "aaq" );
    if ( params_.has( "aan" ) ) aaexotical.AAN = params_.getFloat( "aan" );
    if ( params_.has( "aaf0" ) ) aaexotical.AAF0 = params_.getFloat( "aaf0" );
    if ( params_.has( "aaw" ) ) aaexotical.AAW  = params_.getFloat( "aaw" );
    if ( params_.has( "aaa2" ) ) aaexotical.AAA2 = params_.getFloat( "aaa2" );

    if ( params_.has( "chideiglu" ) ) chidefpmc.CHIDeIGLU = params_.getInt( "chideiglu" );
    if ( params_.has( "chidex" ) ) chidefpmc.CHIDeX = params_.getFloat( "chidex" );
    if ( params_.has( "chidexp" ) ) chidefpmc.CHIDeXP = params_.getFloat( "chidexp" );
    if ( params_.has( "chides2" ) ) chidefpmc.CHIDeS2 = params_.getFloat( "chides2" );
    chidefpmc.CHIDeS = comEnergy_*comEnergy_;
  
    if ( params_.has( "kmr2delta" ) ) kmr2fpmc.KMR2DELTA = params_.getFloat( "kmr2delta" );
    if ( params_.has( "kmr2q2cut" ) ) kmr2fpmc.KMR2Q2CUT = params_.getFloat( "kmr2q2cut" );
    if ( params_.has( "kmr2surv" ) ) kmr2fpmc.KMR2SURV = params_.getFloat( "kmr2surv" );
    if ( params_.has( "kmr2scale" ) ) kmr2fpmc.KMR2SCALE = params_.getFloat( "kmr2scale" );

    if ( params_.has( "isoftm" ) ) xsect.ISOFTM = params_.getInt( "isoftm" );

    if ( params_.has( "zion" ) ) ion.ZION = params_.getInt( "zion" );
    if ( params_.has( "aion" ) ) ion.AION = params_.getInt( "aion" );
    if ( params_.has( "ubmin" ) ) ion.RBMIN = params_.getFloat( "ubmin" );

    //--- Initialize model/pdf dependant parameters
    hwmodini();

    //--- Compute parameter dependent constants
    hwuinc();

    //--- Check POMWIG Settings + Initialisations for consistency
    //
    hwchek();

    //--- Call HWUSTA to make any particle stable
    int iopt = 1;
    int iwig = 0;
    char nwig[9] = "        ";

    int ipdg = 111;
    hwuidt( &iopt, &ipdg, &iwig, nwig );
    if ( ipdg ) hwusta(nwig, 1);
 
    //--- Initialize elementary process
    hweini();

    //---Initialize event record fixing : this will replace the beam 
    //   electrons by protons, radiated photons by pomerons/reggeons etc
    int init = 1;
    hwfxer( &init );
  }

  void
  Fpmc::end()
  {}

  bool
  Fpmc::run()
  {
    // Call herwig routines to create HEPEVT

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

    if ( hwevnt.IERROR ) return false;
#ifndef HEPMC_VERSION3
    hepMCEvt_ = std::make_shared<HepMC::GenEvent>( *conv_.read_next_event() );
#else
    conv_.print_hepevt();
    conv_.HEPEVT_to_GenEvent( hepMCEvt_.get() );
#endif
    ++event_;

    hepMCEvt_->set_event_number(event_ - 1);
#ifndef HEPMC_VERSION3
    hepMCEvt_->set_signal_process_id( hwproc.IPROC );
    hepMCEvt_->set_event_scale( -1. );
#endif
    hepMCEvt_->weights().push_back( hwevnt.EVWGT );

#ifndef HEPMC_VERSION3
    HepMC::PdfInfo pdfInfo;
    pdfInfo.set_x1( hwhard.XX[0] );
    pdfInfo.set_x2( hwhard.XX[1] );
    pdfInfo.set_scalePDF( hwhard.EMSCA );
    hepMCEvt_->set_pdf_info( pdfInfo );
#else
    auto pdf_info = std::make_shared<HepMC::GenPdfInfo>();
    pdf_info->set( 0, 0, hwhard.XX[0], hwhard.XX[1], hwhard.EMSCA, 0, 0 );
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


    //******** Verbosity ********
    if ( event_<=maxEventsToPrint_ && hepMCVerbosity_ ) {
#ifndef HEPMC_VERSION3
      // Prints HepMC event
      std::ostringstream oss;
      oss << "\n----------------------" << endl	
          << "Event process id = " << hepMCEvt_->signal_process_id() << endl; 
      std::cout << oss.str();
      hepMCEvt_->print();
#endif
    }

    return true;
  }

  void
  Fpmc::parseInputCard( const char* filename, Fpmc::Parameters& params )
  {
    std::ifstream card( filename );
    std::string buf;
    //std::regex rgx_parse( "(\\w+)\\s+(\\S+)" );
    std::regex rgx_parse( "(\\w+)[ ']+([^'\\n\\t]+)" );
    std::smatch match;
    while ( !card.eof() ) {
      std::getline( card, buf );
      if ( !std::regex_match( buf, match, rgx_parse ) ) continue;
      std::string key = match[1];
      std::transform( key.begin(), key.end(), key.begin(), ::tolower );
      params.add( key, match[2] );
    }
  }

  void
  Fpmc::writeInputCard( const char* output, const Fpmc::Parameters& params )
  {
    std::ofstream out( output );
    for ( const auto& pair : params.map() ) {
      std::ostringstream os;
      std::string key;
      std::transform( pair.first.begin(), pair.first.end(), key.begin(), ::toupper );
      os << std::setw( 10 ) << key
         << std::setw( 50 ) << pair.second << std::endl;
      out << os.str();
    }
    out.close();
  }

  void
  Fpmc::initialiseParams()
  {
    params_.add( "rmass", 0. );
    params_.add( "wmass", 80.425 );
    params_.add( "rmass", 0. );
    params_.add( "wmass", 80.425 );
    params_.add( "hmass", 125.0 );
    params_.add( "tmass", 174.3 );
    params_.add( "mst1", 250. );
    params_.add( "msb1", 250. );
    params_.add( "ecms", 14.e3 );
    params_.add( "yjmin", -6. );
    params_.add( "yjmax", 6. );
    params_.add( "ptmin", 0. );
    params_.add( "ptmax", 1.e8 );
    params_.add( "emmin", 10. );
    params_.add( "emmax", 1.e8 );
    params_.add( "dkappa", 0. );
    params_.add( "acw", 0. );
    params_.add( "a0w", 0. );
    params_.add( "a0z", 0. );
    params_.add( "acz", 0. );
    params_.add( "a1a", 0. );
    params_.add( "a2a", 0. );
    params_.add( "aam", 0. );
    params_.add( "aaq", 0. );
    params_.add( "aan", 0. );
    params_.add( "aaf0", 0. );
    params_.add( "aaf0z", 0. );
    params_.add( "aaf0w", 0. );
    params_.add( "aaf0zg", 0. );
    params_.add( "aaw", 0. );
    params_.add( "aaa2", 0. );
    params_.add( "chidex", -1. );
    params_.add( "chidexp", -1. );
    params_.add( "chides2", -1. );
    params_.add( "xi1min", -1. );
    params_.add( "xi1max", -1. );
    params_.add( "xi2min", -1. );
    params_.add( "xi2max", -1. );
    params_.add( "chidegapmin", 0. );
    params_.add( "chidegapmax", 0. );
    params_.add( "kmr2q2cut", 2. );
    params_.add( "kmr2surv", 0.3 );
    params_.add( "kmr2scale", 0.618 );
    //
    params_.add( "dlambda", 0. );
    params_.add( "anomcutoff", -1 );
    params_.add( "ywwmin", 0. );
    params_.add( "ywwmax", 0.1 );
    params_.add( "q2wwmn", 0. );
    params_.add( "q2wwmx", 4. );
    //
    params_.add( "output", 1 );
    params_.add( "outputlhe", 0 );
    params_.add( "maxev", 1000 );
    params_.add( "iproc", 16010 );
    params_.add( "nflux", 15 );
    params_.add( "nrn1", 33799 );
    params_.add( "nrn2", 11799 );
    params_.add( "ifit", 10 );
    params_.add( "isoftm", 1 );
    params_.add( "zion", 1 );
    params_.add( "aion", 1 );
    params_.add( "bmin", 1. );
    params_.add( "aaanom", 0 );
    params_.add( "aaexotic", 0 );
    params_.add( "chide_iglu", -1 );
    params_.add( "kmr2_delta", 1 );
    //
    params_.add( "hadr", "Y" );
    //
    params_.add( "typepr", "EXC" );
    //
    params_.add( "typint", "QED" );
    //
    params_.add( "part1", "E+" );
    //
    params_.add( "part2", "E+" );
    //
    params_.add( "modpdf1", -1 );
    params_.add( "modpdf2", -1 );
    //
    params_.add( "ntname", "tmpntuple.ntp" );
    params_.add( "chidepath", "External/CHIDe/Data/" );
    //
    params_.add( "lhefile", "FPMC.lhe" );
  }
}
