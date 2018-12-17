# ifndef __PDB_H
# define __PDB_H

# include <map>
# include <cassert>
# include <cstring>
# include <string>
# include <vector>
# include <fstream>
# include "utility.h"
# include "geometry.h"
using namespace std;

# ifndef PDB_LINEWIDTH
# define PDB_LINEWIDTH 100
# endif

# ifndef INVALID_RESIDUE
# define INVALID_RESIDUE -99999
# endif

# ifndef INVALID_INDEX
# define INVALID_INDEX -1
# endif

struct Residue{
    int                     residue_id;
    string                  residue_name;
    std::vector<string>     atom_names;
    std::vector<Coordinate> coordinates;
};

struct PDBChain{
    char  chain;
    std::vector<Residue> residues;
};

inline bool check_residue(Residue const& residue ){
    return residue.atom_names.size() == residue.coordinates.size();
}

inline int residue_size( Residue const& residue ){
    return static_cast<int>(residue.atom_names.size());
}

inline int chain_size( PDBChain const& chain ){
    return static_cast<int>(chain.residues.size());
}

int search_atom( Residue const& residue, const string& atom_name ){
    int n = residue_size(residue);
    for(int i=0; i < n ; ++i )
        if( residue.atom_names[i] == atom_name )
            return i;
    return INVALID_INDEX;
}

Coordinate get_coordinate( Residue const& residue, const string& atom_name ){
    int i;
    assert( (i = search_atom(residue, atom_name)) != INVALID_INDEX );
    return residue.coordinates[i];
}


static inline void clean_residue( Residue& residue ){
    residue.residue_id = INVALID_RESIDUE;
    residue.residue_name = "UNK";
    residue.atom_names.clear();
    residue.coordinates.clear();
}


PDBChain read_pdb( string const& pdbfile, 
                   const char chain )
{
   assert( is_file(pdbfile) );
   
   PDBChain molecule;
   Residue  residue;
   
   char line[PDB_LINEWIDTH + 1];
   bzero(line, PDB_LINEWIDTH + 1);
   ifstream f(pdbfile.c_str());

   molecule.chain = chain;

   string atom_name, residue_name;
   int resid, last_resid = INVALID_RESIDUE;
   float x,y,z;
   
   while( f.getline(line, PDB_LINEWIDTH) )
   {
     if( (strncmp(line, "ATOM  ",6 ) == 0) && 
         (line[21] == chain) && 
         (line[16] == ' ' || line[16] == 'A') )
    {
        atom_name = string(line).substr(12,4);
        trim(atom_name);
        residue_name = string(line).substr(17,3);
        trim(residue_name);
        resid = atoi(string(line).substr(22,4).c_str());
        x = atof(string(line).substr(30,8).c_str());
        y = atof(string(line).substr(38,8).c_str());
        z = atof(string(line).substr(46,8).c_str());
        if( resid != last_resid ){
            if(residue_size(residue) > 0){
                molecule.residues.push_back(residue);
                clean_residue(residue);
                residue.residue_id = resid;
                residue.residue_name = residue_name;
            }
            last_resid = resid;
        }
        residue.coordinates.push_back( Coordinate(x,y,z) );
        residue.atom_names.push_back(atom_name);
     }
   }
   if( residue_size(residue) > 0 ){
       molecule.residues.push_back(residue);
       clean_residue(residue);
   }
   return molecule;  
}

# endif