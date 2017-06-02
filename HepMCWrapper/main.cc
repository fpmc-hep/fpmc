#include "Fpmc.h"

#include <vector>
#include <string>
#include <map>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <iomanip>

#ifdef HEPMC_VERSION2
#include "HepMC/IO_GenEvent.h"
#else
#include "HepMC/WriterAscii.h"
#endif

using namespace std;

void print_help(vector<string> const& required_parameters, vector<string> const& optional_parameters)
{
   stringstream oss;
   oss << "Usage: fpmc-hepmc ";
   vector<string>::const_iterator it_par = required_parameters.begin();
   vector<string>::const_iterator pars_end = required_parameters.end();
   for(; it_par != pars_end; ++it_par){
      oss << "--" <<  *it_par << " <" << *it_par << "> ";
   }
   it_par = optional_parameters.begin();
   pars_end = optional_parameters.end();
   for(; it_par != pars_end; ++it_par){
      oss << "--" <<  *it_par << " <" << *it_par << "> ";
   }
   oss << endl;
   cout << oss.str(); 
}

int main( int argc, char* argv[] )
{
   vector<string> required_parameters_;
   required_parameters_.push_back("cfg");
   required_parameters_.push_back("nevents");
   required_parameters_.push_back("comenergy");
   
   vector<string> optional_parameters_;
   optional_parameters_.push_back("fileout");

   // Read command line parameters 
   vector<string> command_line_parameters_;
   if(argc > 1){
      command_line_parameters_.resize(argc - 1);
      copy(argv + 1, argv + argc, command_line_parameters_.begin());
   }

   // Help option
   vector<string> help_str; help_str.push_back("-h"); help_str.push_back("--help");
   for(vector<string>::const_iterator it_help = help_str.begin();
                                      it_help != help_str.end(); ++it_help){
      if( find(command_line_parameters_.begin(), command_line_parameters_.end(), *it_help) != command_line_parameters_.end() ){ 
         print_help(required_parameters_,optional_parameters_);
         return 0;
      }
   }
 
   // Read required parameters
   map<string,string> required_parameters_map_;
   for(vector<string>::const_iterator it_par = required_parameters_.begin();
                                      it_par != required_parameters_.end(); ++it_par){
      stringstream par_ss; par_ss << "--"; par_ss << *it_par;
      vector<string>::const_iterator it_par_key = find(command_line_parameters_.begin(), command_line_parameters_.end(), par_ss.str());

      if( it_par_key == command_line_parameters_.end() ){
	 stringstream oss;
	 oss << "ERROR: The following parameter was not set: " << *it_par << endl;      
	 throw runtime_error( oss.str() );
      }

      vector<string>::const_iterator it_par_value = it_par_key + 1;

      if(  it_par_value == command_line_parameters_.end() ||
           find(required_parameters_.begin(), required_parameters_.end(), *it_par_value) != required_parameters_.end() ){
	 stringstream oss;
	 oss << "ERROR: Invalid value for parameter: " << *it_par << endl;      
	 throw runtime_error( oss.str() );

      }
      required_parameters_map_[*it_par] = *it_par_value;
   }
   
   // Read optional parameters
   map<string,string> optional_parameters_map_;
   for(vector<string>::const_iterator it_par = optional_parameters_.begin();
                                      it_par != optional_parameters_.end(); ++it_par){
      stringstream par_ss; par_ss << "--"; par_ss << *it_par;
      vector<string>::const_iterator it_par_key = find(command_line_parameters_.begin(), command_line_parameters_.end(), par_ss.str());

      if( it_par_key == command_line_parameters_.end() ) continue; // Parameter not set

      vector<string>::const_iterator it_par_value = it_par_key + 1;

      if(  it_par_value == command_line_parameters_.end() ||
           find(optional_parameters_.begin(), optional_parameters_.end(), *it_par_value) != optional_parameters_.end() ){
	 stringstream oss;
	 oss << "ERROR: Invalid value for parameter: " << *it_par << endl;      
	 throw runtime_error( oss.str() );
      }

      optional_parameters_map_[*it_par] = *it_par_value;
   }

   //----------------------------
   // Required parameters 
   string datacard_ = required_parameters_map_["cfg"];
   unsigned int maxEvents_ = atoi( required_parameters_map_["nevents"].c_str() );
   double comEnergy_ = atof( required_parameters_map_["comenergy"].c_str() );
   
   // Optional parameters 
   string outputFileName_ = "fpmc.hepmc";
   if( optional_parameters_map_.find("fileout") != optional_parameters_map_.end() ) outputFileName_ = optional_parameters_map_["fileout"];

   stringstream oss;
   oss  << "=========================================================" << endl
        << "FPMC (Wrapper) will initialize with parameters: " << endl
        << "  Datacard:    " << datacard_ << endl
        << "  N events:    " << maxEvents_ << endl
        << "  COM energy:  " << comEnergy_ << endl
        << "  Output file: " << outputFileName_ << endl
        << "=========================================================" << endl;
   cout << oss.str();

   fpmc::Fpmc* generator = new fpmc::Fpmc( comEnergy_, datacard_.c_str() );
   generator->begin();

#ifdef HEPMC_VERSION2
   HepMC::IO_GenEvent output( outputFileName_, ios::out );
#else
   HepMC::WriterAscii output( outputFileName_ );
#endif
   for(unsigned int evt = 0; evt < maxEvents_; ++evt){
      cout << "[FPMC Wrapper] Processing event " << (evt + 1) << endl;
      bool success = generator->run();
      if(!success){
         cout << "[FPMC Wrapper] WARNING: Event " << (evt + 1) << " failed." << endl;
         continue;
      }
#ifdef HEPMC_VERSION2
      output.write_event( generator->event() );
#else
      output.write_event( *generator->event() );
#endif
   }    
   generator->end();

   return 0;
}
