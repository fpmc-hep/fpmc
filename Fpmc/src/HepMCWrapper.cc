#include "Fpmc/HepMCWrapper.h"

#ifndef HEPMC_VERSION2
extern HEPEVT hepevt_;
#endif

namespace fpmc
{
  HepMCWrapper::HepMCWrapper( double comEnergy, const char* card ) :
    Fpmc( card ), hepMCVerbosity_( true )
  {
    params_.setSqrtS( comEnergy );
    params_.dump();
  }

  HepMCWrapper::~HepMCWrapper()
  {}

  const HepMC::GenEvent*
  HepMCWrapper::event()
  {
    //----- start by generating the next event with FPMC

    hwevnt_t last_evt;
    if ( !Fpmc::next( last_evt ) ) return 0;
    if ( last_evt.IERROR ) return 0;

#ifdef HEPMC_VERSION2
    hepMCEvt_ = std::make_shared<HepMC::GenEvent>( *conv_.read_next_event() );
#else
    HepMC::HEPEVT_Wrapper::print_hepevt();
    HepMC::HEPEVT_Wrapper::HEPEVT_to_GenEvent( hepMCEvt_.get() );
#endif

    hepMCEvt_->set_event_number( event_-1 );
#ifdef HEPMC_VERSION2
    hepMCEvt_->set_signal_process_id( params_.processId() );
    hepMCEvt_->set_event_scale( -1. );
#endif
    hepMCEvt_->weights().push_back( last_evt.EVWGT );

#ifdef HEPMC_VERSION2
    HepMC::PdfInfo pdfInfo;
    pdfInfo.set_x1( hwhard_.XX[0] );
    pdfInfo.set_x2( hwhard_.XX[1] );
    pdfInfo.set_scalePDF( hwhard_.EMSCA ); //FIXME
    hepMCEvt_->set_pdf_info( pdfInfo );
#else
    auto pdf_info = std::make_shared<HepMC::GenPdfInfo>();
    pdf_info->set( 0, 0, hwhard_.XX[0], hwhard_.XX[1], hwhard_.EMSCA, 0, 0 ); //FIXME
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

    return hepMCEvt_.get();
  }

  void
  HepMCWrapper::write( const char* out )
  {
    if ( !hepMCEvt_.get() ) return;

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
}
