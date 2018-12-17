# ifndef __CONFIG_READER_H
# define __CONFIG_READER_H

# include <set>
# include <map>
# include <cstdlib>
# include <cstring>
# include <cassert>
# include <iomanip>
# include <iostream>
# include "utility.h"
# include "geometry.h"
# include "trajectory.h"
# include "trajectory_reader.h"
using namespace std;

# ifndef CONFIG_LINEWIDTH
# define CONFIG_LINEWIDTH 200
# endif

struct BoundingBoxConfig{
  Coordinate box_max;
  Coordinate box_min;
};

BoundingBoxConfig read_boundingbox_configfile( char const* filename )
{
  ifstream f(filename);
  assert( f.good() );
  
  map<string, bool> map_check;
  
  char   line[CONFIG_LINEWIDTH + 1];
  
  set<string>  filenames;
  set<int>     residue_ids;
  vector<Trajectory> trajectories;
  bzero( line, CONFIG_LINEWIDTH + 1);
  BoundingBoxConfig config;
  
  map_check["min_coord"] = false;
  map_check["max_coord"] = false;
  map_check["trajectory"] = false;
  map_check["residues"] = false;
  
  while( f.getline(line, CONFIG_LINEWIDTH) ){
      vector<string> words = string_split(string(line),",");
      if( words.size() > 0 ){
         if(words[0] == "m" ){ //maximum coordinate
           assert( words.size() == 4 );
           assert( is_float(words[1]) && is_float(words[2]) && is_float(words[3]) );
		   config.box_min = Coordinate( atof(words[1].c_str()), 
										atof(words[2].c_str()), 
										atof(words[3].c_str()) );
		   map_check["min_coord"] = true;
         }else if( words[0] ==  "n" ){ //minimum coordinate
           assert( words.size() == 4 );
           assert( is_float(words[1]) && is_float(words[2]) && is_float(words[3]) );
		   config.box_max = Coordinate( atof(words[1].c_str()), 
										atof(words[2].c_str()), 
										atof(words[3].c_str()) );
		   map_check["max_coord"] = true;
         }else if( words[0] == "t" ){ //trajectory file
		   assert( words.size() > 1 );
		   for( int i=1; i < words.size(); ++i )
			   if( is_file( words[i]) )
				   filenames.insert(words[i]);
				   
		   map_check["trajectory"] = filenames.size() > 0;
		 }else if( words[0] == "r" ){ // residue ids
		   assert( words.size() > 1);
		   for( int i=1; i < words.size(); ++i )
			   if( is_integer(words[i]) )
				   residue_ids.insert(atoi(words[i].c_str()));
				   
		   map_check["residues"] = residue_ids.size() > 2;
		 } 
      }
      bzero( line, CONFIG_LINEWIDTH);
  }
  f.close();
  if( (map_check["min_coord"] && map_check["max_coord"]) && 
      !(map_check["residues"] && map_check["trajectory"]) )
	  return config;

  assert( map_check["residues"] && map_check["trajectory"] );

  for( typename set<string>::const_iterator i = filenames.begin(); i != filenames.end(); ++i)
     trajectories.push_back( ::read_trajectory(*i) );

  int counter = 0;
  for( typename vector<Trajectory>::const_iterator j=trajectories.begin(); j != trajectories.end(); ++j)
  {
	 for( typename Trajectory::const_iterator k = j->begin(); k != j->end(); ++k )
	 {
       for( typename set<int>::const_iterator i = residue_ids.begin(); i != residue_ids.end(); ++i )
	   {
		 assert( find( k->residueIds.begin(),k->residueIds.end(),*i ) != k->residueIds.end() );
		 Coordinate crd = ::fetch_coordinate(*k, *i);
		 if(counter == 0)
		 {
			 if( ! map_check["max_coord"] ) 
				 config.box_max = crd;
			 if( ! map_check["min_coord"] ) 
				 config.box_min = crd;
		 }else{
			 if( ! map_check["max_coord"] )
			 {
				if(config.box_max.x < crd.x) config.box_max.x = crd.x;
				if(config.box_max.y < crd.y) config.box_max.y = crd.y;
				if(config.box_max.z < crd.z) config.box_max.z = crd.z;
			 }
			 if( ! map_check["min_coord"] )
			 {
				if(config.box_min.x > crd.x) config.box_min.x = crd.x;
				if(config.box_min.y > crd.y) config.box_min.y = crd.y;
				if(config.box_min.z > crd.z) config.box_min.z = crd.z;
			 }
		 }
		 counter++;
	   }
     }
  }
  return config;
}


ostream& operator << (ostream& os, BoundingBoxConfig const& config )
{
	os << "m" << "," << fixed << setprecision(3) << config.box_min.x 
			  << "," << fixed << setprecision(3) << config.box_min.y 
			  << "," << fixed << setprecision(3) << config.box_min.z << endl 
	   << "n" << "," << fixed << setprecision(3) << config.box_max.x 
	          << "," << fixed << setprecision(3) << config.box_max.y 
			  << "," << fixed << setprecision(3) << config.box_max.z;
	return os;
}


# endif
