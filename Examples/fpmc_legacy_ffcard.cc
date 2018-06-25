#include <iostream>
#include <string>
#include <map>
#include <unistd.h>

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif
  map<string,string> m_data;

  void ffinit_( int& card )
  {
    string data, key;
    unsigned short i = 0;
    while ( cin >> data ) {
      if ( i % 2 == 0 )
        key = data;
      else
        m_data[key] = data;
      ++i;
    }
  }

  void ffset_( const char* skey, int& idum ) {}

  void ffkey_( char* key, int& address, int& length, char* type, int key_len, int type_len )
  {
    string skey( key, key_len ), stype( type, 1 );
    if ( m_data.count( skey ) == 0 )
      return;
    const auto value = m_data.at( skey );
    if ( stype == "i" )
      address = stoi( value );
    else if ( stype == "r" ) {
      auto val = stof( value );
      address = reinterpret_cast<int&>( val );
    }
    else if ( stype == "m" ) {
      const size_t size = value.size()-2; // strip the trailing "'" + '\0'
      char val[size];
      value.copy( val, size, 1 ); // strip the leading "'"
      address = reinterpret_cast<int&>( val );
    }
  }

  void ffgo_() {}
#ifdef __cplusplus
}
#endif

