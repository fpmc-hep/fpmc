#include "FpmcParameters.h"

#include <string>

namespace fpmc
{
  FpmcParameters::FpmcParameters() :
    std::map<std::string,std::string>( {
      // particles masses
      { "rmass", "0.0" }, { "wmass", "80.425" }, { "hmass", "125.0" }, { "tmass", "174.3" }, { "mst1", "250.0" }, { "msb1", "250.0" },
      // kinematics cuts
      { "ecms", "14000.0" }, { "yjmin", "-6.0" }, { "yjmax", "6.0" }, { "ptmin", "0.0" }, { "ptmax", "1.e8" }, { "emmin", "10.0" }, { "emmax", "1.e8" },
      { "ywwmin", "0.0" }, { "ywwmax", "0.1" }, { "q2wwmn", "0.0" }, { "q2wwmx", "4.0" },
      // events production
      { "maxev", "1000" },
      { "output", "1" }, { "outputlhe", "0" }, { "ntname", "tmpntuple.ntp" }, { "lhefile", "FPMC.lhe" },
      { "hadr", "Y" },
      // process definition
      { "part1", "E+" }, { "part2", "E+" },
      { "zion", "1" }, { "aion", "1" }, { "bmin", "1.0" },
      { "iproc", "16010" }, { "nflux", "15" }, { "ifit", "10" }, { "isoftm", "1" }, { "typepr", "EXC" }, { "typint", "QED" },
      { "nrn1", "33799" }, { "nrn2", "11799" },
      // anomalous couplings
      { "dkappa", "0.0" }, { "acw", "0.0" }, { "a0w", "0.0" }, { "a0z", "0.0" }, { "acz", "0.0" }, { "a1a", "0.0" }, { "a2a", "0.0" },
      { "aaanom", "0" }, { "aaexotic", "0" },
      { "aam", "0.0" }, { "aaq", "0.0" }, { "aan", "0.0" },
      { "aaf0", "0.0" }, { "aaf0z", "0.0" }, { "aaf0w", "0.0" }, { "aaf0zg", "0.0" }, { "aaw", "0.0" }, { "aaa2", "0.0" },
      { "dlambda", "0.0" }, { "anomcutoff", "-1.0" }, {"xi1ttbar", "0.0"}, {"xi2ttbar", "0.0"}, {"xi3ttbar", "0.0"}, 
      {"xi4ttbar", "0.0"}, {"xi5ttbar", "0.0"}, {"xi6ttbar", "0.0"}, 
      // CHIDe parameters
      { "chidex", "-1.0" }, { "chidexp", "-1.0" }, { "chides2", "-1.0" }, { "chidegapmin", "0.0" }, { "chidegapmax", "0.0" }, { "chide_iglu", "-1" }, { "chidepath", "External/CHIDe/Data/" },
      { "xi1min", "-1.0" }, { "xi1max", "-1.0" }, { "xi2min", "-1.0" }, { "xi2max", "-1.0" },
      // KMR parameters
      { "kmr2q2cut", "2.0" }, { "kmr2surv", "0.3" }, { "kmr2scale", "0.618" }, { "kmr2_delta", "1" },
      // PDF sets
      { "autpdf1", "HWLHAPDF" }, { "autpdf2", "HWLHAPDF" }, { "modpdf1", "-1" }, { "modpdf2", "-1" }
    } ) { validate(); }

  FpmcParameters
  FpmcParameters::parseCard( const char* filename )
  {
    FpmcParameters param;
    std::ifstream card( filename );
    std::string buf;
    std::regex rgx_parse( "(\\w+)\\s+(\\S+)" );
    //std::regex rgx_parse( "(\\w+)[ ']+([^'\\n\\t]+)" );
    std::smatch match;
    while ( !card.eof() ) {
      std::getline( card, buf );
      if ( !std::regex_match( buf, match, rgx_parse ) ) continue;
      std::string key = match[1], value = match[2];
      value.erase( std::find_if( value.rbegin(), value.rend(), std::bind1st( std::not_equal_to<char>(), ' ' ) ).base(), value.end() ); // remove trailing whitespaces
      std::transform( key.begin(), key.end(), key.begin(), ::tolower ); // transform the key to lowercase
      value.erase( std::remove( value.begin(), value.end(), '\''), value.end() );
      param.add( key, value );
    }
    param.validate();
    return param;
  }

  void
  FpmcParameters::writeCard( const char* output ) const
  {
    std::ofstream out( output );
    for ( const auto& pair : map() ) {
      std::ostringstream os;
      std::string key = pair.first, value = pair.second;
      std::transform( key.begin(), key.end(), key.begin(), ::toupper );
      char* buf; strtod( value.c_str(), &buf );
      if ( *buf!=0 ) value = "'" + value + "'";
      os << std::left
         << std::setw( 15 ) << key
         << std::setw( 50 ) << value << std::endl;
      out << os.str();
    }
    out.close();
  }

  void
  FpmcParameters::validate()
  {
    if ( has( "xi1min" ) && getFloat( "xi1min" ) < 0.0 ) add( "xi1min", getFloat( "ywwmin" ) );
    if ( has( "xi2min" ) && getFloat( "xi2min" ) < 0.0 ) add( "xi2min", getFloat( "ywwmin" ) );
    if ( has( "xi1max" ) && getFloat( "xi1max" ) < 0.0 ) add( "xi1max", getFloat( "ywwmax" ) );
    if ( has( "xi2max" ) && getFloat( "xi2max" ) < 0.0 ) add( "xi2max", getFloat( "ywwmax" ) );
    if ( has( "nflux" ) ) {
      const unsigned int nflux = getInt( "nflux" );
      if ( nflux == 18 ) {
        if ( has( "chide_iglu" ) && getFloat( "chide_iglu" ) < 0 ) add( "chide_iglu", 4 ); // impact factor parameterisation
        // scaling lower limit of Sudakov factor integration
        if ( has( "chidexp" ) && getFloat( "chidexp" ) < 0.0 ) {
          if ( has( "iproc" ) && getInt( "iproc" ) == 16059 ) add( "chidexp", 1.0 );
          else add( "chidexp", 0.5 );
        }
        // rapidity gap survival probability
        if ( has( "chides2" ) && getFloat( "chides2" ) < 0.0 ) {
          if ( has( "ecms" ) ) {
            const float ecms = getFloat( "ecms" );
            if ( labs( ecms-14.e3 ) < 1.0 ) add( "chides2", 0.03 ); // ~LHC
            else if ( labs( ecms-1960. ) < 1.0 ) add( "chides2", 0.10 ); // ~Tevatron
            else throw std::runtime_error( "ERROR: CHIDeS2 not specified!\n"
                                           "  For nonstandard ECMS (different than 1960 and 14000)\n"
                                           "  CHIDeS2 must be specified in data card." );
          }
        }
        if ( has( "iproc" ) ) {
          const unsigned int iproc = getInt( "iproc" );
          if ( iproc == 16012 && has( "chidex" ) && getFloat( "chidex" ) < 0.0 ) add( "chidex", 0.5 );
          if ( iproc == 16059 && has( "chidex" ) && getFloat( "chidex" ) < 0.0 ) add( "chidex", 0.5 );
          if ( iproc == 19999 && has( "chidex" ) && getFloat( "chidex" ) < 0.0 ) add( "chidex", 1.0 );
        }
      }
      else if ( nflux == 16 ) {
        if ( has( "chides2" ) && getFloat( "chides2" ) > 0.0 ) add( "kmr2surv", getFloat( "chides2" ) );
        else add( "kmr2surv", 0.03 );
      }
    }
  }

  void
  FpmcParameters::fetchHWBMCH( hwbmch_t& hwbmch ) const
  {
    std::string def_str( "        " );
    def_str.copy( hwbmch.PART1, 8 );
    def_str.copy( hwbmch.PART2, 8 );
    if ( has( "part1" ) ) getString( "part1" ).copy( hwbmch.PART1, 8 );
    if ( has( "part2" ) ) getString( "part2" ).copy( hwbmch.PART2, 8 );
  }

  void
  FpmcParameters::fetchHWPROC( hwproc_t& hwproc ) const
  {
    hwproc.PBEAM1 = hwproc.PBEAM2 = sqrtS()/2.;
    if ( has( "iproc" ) ) hwproc.IPROC = processId();
    if ( has( "maxev" ) ) hwproc.MAXEV = getInt( "maxev" );
  }

  void
  FpmcParameters::fetchHWPRAM( hwpram_t& hwpram ) const
  {
    if ( has( "modpdf1" ) ) hwpram.MODPDF[0] = getInt( "modpdf1" );
    if ( has( "modpdf2" ) ) hwpram.MODPDF[1] = getInt( "modpdf2" );
  }

  void
  FpmcParameters::fetchHWPRCH( hwprch_t& hwprch ) const
  {
    if ( has( "autpdf1" ) ) getString( "autpdf1" ).copy( hwprch.AUTPDF[0], 8 );
    if ( has( "autpdf2" ) ) getString( "autpdf2" ).copy( hwprch.AUTPDF[1], 8 );
  }

  void
  FpmcParameters::fetchHWPROP( hwprop_t& hwprop ) const
  {
    if ( has( "tmass" ) ) hwprop.RMASS[6] = hwprop.RMASS[12] = getFloat( "tmass" ); // top mass
    if ( has( "wmass" ) ) hwprop.RMASS[198] = hwprop.RMASS[199] = getFloat( "wmass" ); // W+/- mass
    if ( has( "hmass" ) ) hwprop.RMASS[201] = getFloat( "hmass" ); // higgs mass
    if ( has( "msb1" ) ) hwprop.RMASS[405] = getFloat( "msb1" ); // stop2 mass
    if ( has( "mst1" ) ) hwprop.RMASS[406] = getFloat( "mst1" ); // stop1 mass
  }

  void
  FpmcParameters::fetchHWHARD( hwhard_t& hwhard ) const
  {
    if ( has( "q2wwmn" ) ) hwhard.Q2WWMN = getFloat( "q2wwmn" );
    if ( has( "q2wwmx" ) ) hwhard.Q2WWMX = getFloat( "q2wwmx" );
    // beam momentum loss range
    if ( has( "ywwmin" ) ) hwhard.YWWMIN = getFloat( "ywwmin" );
    if ( has( "ywwmax" ) ) hwhard.YWWMAX = getFloat( "ywwmax" );
    // central products rapidity, pt, and mass cuts
    if ( has( "yjmin" ) ) hwhard.YJMIN = getFloat( "yjmin" );
    if ( has( "yjmax" ) ) hwhard.YJMAX = getFloat( "yjmax" );
    if ( has( "ptmin" ) ) hwhard.PTMIN = getFloat( "ptmin" );
    if ( has( "ptmax" ) ) hwhard.PTMAX = getFloat( "ptmax" );
    if ( has( "emmin" ) ) hwhard.EMMIN = getFloat( "emmin" );
  }

  void
  FpmcParameters::fetchXSECT( xsect_t& xsect ) const
  {
    //if ( has( "kmr2surv" ) ) xsect.GAPSPR = getFloat( "kmr2surv" );
    if ( has( "nflux" ) ) xsect.NFLUX = getInt( "nflux" );
    if ( has( "isoftm" ) ) xsect.ISOFTM = getInt( "isoftm" );
  }

  void
  FpmcParameters::fetchPDFS( pdfs_t& pdfs ) const
  {
    if ( has( "ifit" ) ) pdfs.IFITPDF = getInt( "ifit" );
  }

  void
  FpmcParameters::fetchPRTYPE( prtype_t& prtype ) const
  {
    if ( has( "typepr" ) ) getString( "typepr" ).copy( prtype.TYPEPR, 3 );
    if ( has( "typint" ) ) getString( "typint" ).copy( prtype.TYPINT, 3 );
  }

  void
  FpmcParameters::fetchAAANOMAL( aaanomal_t& aaanomal ) const
  {
    if ( has( "aaanom" ) ) aaanomal.AAANOM = getInt( "aaanom" );
    if ( has( "dkappa" ) ) aaanomal.D_KAPPA = getFloat( "dkappa" );
    if ( has( "lambda" ) ) aaanomal.LAMBDA = getFloat( "lambda" );
    if ( has( "a0w" ) ) aaanomal.A0W = getFloat( "a0w" );
    if ( has( "acw" ) ) aaanomal.ACW = getFloat( "acw" );
    if ( has( "a0z" ) ) aaanomal.A0Z = getFloat( "a0z" );
    if ( has( "acz" ) ) aaanomal.ACZ = getFloat( "acz" );
    if ( has( "a1a" ) ) aaanomal.A1A = getFloat( "a1a" );
    if ( has( "a2a" ) ) aaanomal.A2A = getFloat( "a2a" );
    if ( has( "xi1ttbar" ) ) aaanomal.XI1TTBAR = getFloat( "xi1ttbar" );
    if ( has( "xi2ttbar" ) ) aaanomal.XI2TTBAR = getFloat( "xi2ttbar" );
    if ( has( "xi3ttbar" ) ) aaanomal.XI3TTBAR = getFloat( "xi3ttbar" );
    if ( has( "xi4ttbar" ) ) aaanomal.XI4TTBAR = getFloat( "xi4ttbar" );
    if ( has( "xi5ttbar" ) ) aaanomal.XI5TTBAR = getFloat( "xi5ttbar" );
    if ( has( "xi6ttbar" ) ) aaanomal.XI6TTBAR = getFloat( "xi6ttbar" );
    if ( has( "anomcutoff" ) ) aaanomal.ANOMCUTOFF = getFloat( "anomcutoff" );
  }

  void
  FpmcParameters::fetchAAEXOTICAL( aaexotical_t& aaexotical ) const
  {
    if ( has( "aaexotic" ) ) aaexotical.AAEXOTIC = getInt( "aaexotic" );
    if ( has( "aam" ) ) aaexotical.AAM = getFloat( "aam" );
    if ( has( "aaq" ) ) aaexotical.AAQ = getFloat( "aaq" );
    if ( has( "aan" ) ) aaexotical.AAN = getFloat( "aan" );
    if ( has( "aaf0" ) ) aaexotical.AAF0 = getFloat( "aaf0" );
    if ( has( "aaf0w" ) ) aaexotical.AAF0W = getFloat( "aaf0w" );
    if ( has( "aaf0z" ) ) aaexotical.AAF0Z = getFloat( "aaf0z" );
    if ( has( "aaf0zg" ) ) aaexotical.AAF0ZG = getFloat( "aaf0zg" );
    if ( has( "aaw" ) ) aaexotical.AAW  = getFloat( "aaw" );
    if ( has( "aaa2" ) ) aaexotical.AAA2 = getFloat( "aaa2" );
  }

  void
  FpmcParameters::fetchCHIDEFPMC( chidefpmc_t& chidefpmc ) const
  {
    if ( has( "chideiglu" ) ) chidefpmc.CHIDeIGLU = getInt( "chideiglu" );
    if ( has( "chidex" ) ) chidefpmc.CHIDeX = getFloat( "chidex" );
    if ( has( "chidexp" ) ) chidefpmc.CHIDeXP = getFloat( "chidexp" );
    if ( has( "chides2" ) ) chidefpmc.CHIDeS2 = getFloat( "chides2" );
  }

  void
  FpmcParameters::fetchKMR2FPMC( kmr2fpmc_t& kmr2fpmc ) const
  {
    if ( has( "kmr2delta" ) ) kmr2fpmc.KMR2DELTA = getFloat( "kmr2delta" );
    if ( has( "kmr2q2cut" ) ) kmr2fpmc.KMR2Q2CUT = getFloat( "kmr2q2cut" );
    if ( has( "kmr2surv" ) ) kmr2fpmc.KMR2SURV = getFloat( "kmr2surv" );
    if ( has( "kmr2scale" ) ) kmr2fpmc.KMR2SCALE = getFloat( "kmr2scale" );
  }

  void
  FpmcParameters::fetchION( ion_t& ion ) const
  {
    if ( has( "zion" ) ) ion.ZION = getInt( "zion" );
    if ( has( "aion" ) ) ion.AION = getInt( "aion" );
    if ( has( "ubmin" ) ) ion.RBMIN = getFloat( "ubmin" );
  }

  void
  FpmcParameters::fetchCYFFLONG1( cyfflong1_t& cyfflong1 ) const
  {
    if ( has( "chidepath" ) ) getString( "chidepath" ).copy( cyfflong1.CHIDePATH, 32 );
  }

  void
  FpmcParameters::fetchCHIDEPATH( chidepath_t& chidepath ) const
  {
    if ( has( "chidepath" ) ) getString( "chidepath" ).copy( chidepath.sudatab, 50 );
  }

  void
  FpmcParameters::fetchCHIDEENV( chideenv_t& chideenv ) const
  {
    if ( has( "chidepath" ) ) getString( "chidepath" ).copy( chideenv.CHIDe_PATH, 500 );
  }
}
