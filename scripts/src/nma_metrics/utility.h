# ifndef __UTILITY_H
# define __UTILITY_H

# include <map>
# include <locale>
# include <cmath>
# include <cctype>
# include <vector>
# include <string>
# include <cstdio>
# include <cstdlib>
# include <sstream>
# include <functional>
# include <algorithm>
# include <sys/stat.h>
# include <sys/types.h>
using namespace std;

string base_directory();
bool is_file( char const* );
bool is_file( string const& );
bool is_integer( string const& );
bool is_float( string const& );
bool is_dir( char const* );
bool is_dir( string const& );
vector<string> string_split(string const& , string const& );


string base_directory()
{
   return "/home/sumanta/Work/codelite_base/nma_metrics/volume_fluctuation";
}


bool is_file( char const* filename ){
  FILE *fp;
  if( (fp = fopen(filename, "r")) != 0 ){
     fclose(fp);
     return true;
  }
  return false;
}


bool is_file( string const& filename ){
   return is_file(filename.c_str());
}


bool is_dir( char const* dir ){
    struct stat info;
    if( stat(dir, &info) != 0)
        return false;
    else if( info.st_mode & S_IFDIR )
        return true;
    return false;
}


bool is_dir( string const& dir ){
    return is_dir(dir.c_str());
}


vector<string> string_split(string const& sentence, string const& delimiter )
{
    vector<string> words;
    string  s(sentence);
    size_t pos = 0;
    std::string token;
    while ((pos = s.find(delimiter)) != std::string::npos) {
      token = s.substr(0, pos);
      words.push_back(token);
      s.erase(0, pos + delimiter.length());
    }
    if( s.size() > 0 ) words.push_back(s);
    return words;
}


bool is_integer( string const& str )
{
  istringstream iss(str);
  int f;
  iss >> noskipws >> f;
  return iss.eof() && ! iss.fail();
}


bool is_float( string const& str )
{
	istringstream iss(str);
	float f;
	iss >> noskipws >> f;
	return iss.eof() && ! iss.fail();
}


template<typename Iterator>
float mean( Iterator start, Iterator end ){
  float s = 0.f;
  int   n = 0;
  for( ; start != end; start++ ){
    s += static_cast<float>( *start );
    n++;
  }
  return (n > 0)?(s/n):s;
}


template<typename Iterator>
float sd( Iterator start, Iterator end ){
   float s  = 0.f;
   float s2 = 0.f;
   int   n  = 0.f;
   for( ; start != end; start++ ){
     float x = static_cast<float>(*start);
     s += x;
     s2 += x*x;
     n++;
   } 
   return (n>0)?static_cast<float>( sqrt( (s2/n) - (s/n)*(s/n) ) ):0.f; 
}

template<typename Key, typename Value>
bool is_in( map<Key,Value> const& m, Key const& key){
    return (m.find(key) != m.end());
}

template<typename Value>
bool is_in( vector<Value> const& v, Value const& d){
	return (find(v.begin(), v.end(), d) != v.end());
}

template<typename Key, typename Value>
vector<Key> keys( map<Key,Value> const& m )
{
	vector<Key> names;
	for(typename map<Key,Value>::const_iterator i = m.begin(); i != m.end(); ++i )
		names.push_back(i->first);
	return names;
}

static inline std::string &ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
            std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

static inline std::string &rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
            std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}

inline std::string &trim(std::string &s) {
    return ltrim(rtrim(s));
}

# endif
