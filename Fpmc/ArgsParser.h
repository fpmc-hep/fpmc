#ifndef Fpmc_ArgsParser_h
#define Fpmc_ArgsParser_h

#include <vector>
#include <string>
#include <map>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <iostream>

namespace fpmc
{
  class ArgsParser
  {
    public:
      ArgsParser( int argc, char* argv[], const std::vector<std::string>& required_parameters, const std::vector<std::string>& optional_parameters ) :
        help_str_( { "-h", "--help" } ),
        required_params_( required_parameters ), optional_params_( optional_parameters ) {
        command_name_ = argv[0];
        if ( argc >1  ) {
          args_.resize( argc-1 );
          std::copy( argv+1, argv+argc, args_.begin() );
        }
        for ( const auto& str : help_str_ ) {
          if ( find( args_.begin(), args_.end(), str ) != args_.end() ) {
            print_help();
            exit( 0 );
          }
        }
      }
      /// Read required parameters
      std::map<std::string,std::string> required_parameters() const {
        std::map<std::string,std::string> out;
        for ( const auto& par : required_params_ ) {
          std::ostringstream par_ss; par_ss << "--" << par;
          const auto key = find( args_.begin(), args_.end(), par_ss.str() );
          if ( key == args_.end() ) {
            std::ostringstream oss; oss << "ERROR: The following parameter was not set: " << par << std::endl;
            throw std::runtime_error( oss.str() );
          }

          const auto value = key + 1;
          if ( value == args_.end() || find( required_params_.begin(), required_params_.end(), *value ) != required_params_.end() ) {
            std::ostringstream oss; oss << "ERROR: Invalid value for parameter: " << par << std::endl;
            throw std::runtime_error( oss.str() );
          }
          out[par] = *value;
        }
        return out;
      }
      /// Read optional parameters
      std::map<std::string,std::string> optional_parameters() const {
        std::map<std::string,std::string> out;
        for ( const auto& par : optional_params_ ) {
          std::ostringstream par_ss; par_ss << "--" << par;
          const auto key = find( args_.begin(), args_.end(), par_ss.str() );
          if ( key == args_.end() ) continue; // Parameter not set

          const auto value = key + 1;
          if ( value == args_.end() || find( optional_params_.begin(), optional_params_.end(), *value ) != optional_params_.end() ) {
            std::ostringstream oss; oss << "ERROR: Invalid value for parameter: " << par << std::endl;
            throw std::runtime_error( oss.str() );
          }
          out[par] = *value;
        }
        return out;
      }
      /// Show usage
      void print_help() const {
        std::ostringstream oss;
        oss << "Usage: " << command_name_ << " ";
        for ( const auto& par : required_params_ ) {
          oss << "--" <<  par << " <" << par << "> ";
        }
        for ( const auto& par : optional_params_ ) {
          oss << "--" <<  par << " [" << par << "] ";
        }
        oss << std::endl;
        std::cout << oss.str();
      }
    private:
      std::string command_name_;
      const std::vector<std::string> help_str_;
      std::vector<std::string> required_params_, optional_params_, args_;
  };
}

#endif

