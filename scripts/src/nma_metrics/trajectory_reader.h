# ifndef __TRAJECTORY_READER_H
# define __TRAJECTORY_READER_H

# include <set>
# include <string>
# include <cstdio>
# include <vector>
# include <cstring>
# include <cstdlib>
# include <cassert>
# include <fstream>
# include "trajectory.h"
# include "coordinate.h"
using namespace std;

# define PDB_LINEWIDTH 100


struct Config{
  Trajectory trajectory;
  ResidueIds bounding_residues;
};

bool       is_file( char const* );
Trajectory read_trajectory( const string& );
Config     read_config( char const* );

bool to_append( vector<Snapshot>& trj, Snapshot& snapshot) {
   return ( valid_snapshot(snapshot) && (trj.size() == 0 || snapshot_size(trj[0]) == snapshot_size(snapshot)) );
}

Trajectory read_trajectory( const string& trjfile ){
   char line[PDB_LINEWIDTH + 1];
   Trajectory trajectory;
   bzero( (void*)line, PDB_LINEWIDTH + 1);
   assert( is_file(trjfile.c_str()) );
   ifstream f(trjfile.c_str());

   string residueType;
   int    resid;
   float  x,y,z;
   Snapshot snapshot;
   bool   append = false;

   while( f.getline(line, PDB_LINEWIDTH) ){
     if( strncmp(line, "MODEL", 5) == 0 ){
         if( to_append(trajectory, snapshot) ) 
           trajectory.push_back(snapshot);
         clear_snapshot(snapshot);
     }
     if( (strncmp(line, "ATOM  ",6 ) == 0) && (string(line).substr(13,2) == "CA" ) ){
        residueType = string(line).substr(17,3);
        resid = atoi(string(line).substr(22,4).c_str());
        x = atof(string(line).substr(30,8).c_str());
        y = atof(string(line).substr(38,8).c_str());
        z = atof(string(line).substr(46,8).c_str());
        snapshot.residueType.push_back(residueType);
        snapshot.residueIds.push_back(resid);
        snapshot.coordinates.push_back( Coordinate(x,y,z));
     }
     bzero( (void*)line, PDB_LINEWIDTH + 1);
   }
   if( snapshot_size(snapshot) > 0 && to_append(trajectory, snapshot) ){
      trajectory.push_back(snapshot);
      clear_snapshot(snapshot);
   }
   return trajectory;
}

bool is_file( char const* filename ){
  FILE *fp;
  if( (fp = fopen(filename, "r")) != 0 ){
     fclose(fp);
     return true;
  }
  return false;
}

Config read_configfile( char const* filename ){
  ifstream f(filename);
  assert( f.good() );
  string   trjfile;
  int      resId;
  set<int> residues;
  f >> trjfile;
  while( f >> resId ){
    residues.insert(resId);
  }
  assert( residues.size() > 3 );
  Config config;
  config.trajectory = read_trajectory( trjfile );
  config.bounding_residues = ResidueIds(residues.begin(), residues.end());
  assert( trajectory_size( config.trajectory ) > 0 );

  for(int i=0; i < config.bounding_residues.size(); ++i){
     assert( find( config.trajectory[0].residueIds.begin(), 
                   config.trajectory[0].residueIds.end(), 
                   config.bounding_residues[i]) != config.trajectory[0].residueIds.end());
  }
  return config;
}

# endif
