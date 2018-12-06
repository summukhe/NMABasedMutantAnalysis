# include <cassert>
# include <cstdlib>
# include <cstdio>
# include <vector>
# include <string>
# include <cstring>
# include <fstream>
# include <iostream>
# include "utility.h"
# include "partial_charge.h"
using namespace std;

bool read_target_files( string const& pdb_dir, 
                        string const& list_file, 
                        vector<string>& pdbfiles,
                        vector<string>& pdbids,
                        vector<char>& chains )
{
    assert( is_dir(pdb_dir) );
    assert( is_file(list_file)); 
    pdbfiles.clear();
    pdbids.clear();
    chains.clear();
    int hit = 0;
    char line[100];
    char pdbfile[100];
    bzero(line, 100);
    bzero(pdbfile,100);
    ifstream f(list_file.c_str());
    while( f.getline(line, 100) ){
        vector<string> flds = string_split(string(line), "," );
        if( flds.size() == 2 && flds[0].size() == 4 && flds[1].size() == 1 ){
            sprintf(pdbfile, "%s/%s.pdb", pdb_dir.c_str(), flds[0].c_str());
            if( is_file(pdbfile) ){
                pdbfiles.push_back( string(pdbfile) );
                pdbids.push_back( flds[0] );
                chains.push_back( flds[1][0] );
                hit++;
            }
        }
        bzero(line, 100);
        bzero(pdbfile,100);
    }
    return hit > 0;
}


int main(int argc, char** argv)
{
    if( argc != 3 || ! is_dir(argv[1]) || ! is_file(argv[2]) ) {
        fprintf(stderr, "Usage: %s <pdb folder> <pdb file and chain list>\n", argv[0]);
        exit(1);
    }
    string pdb_dir(argv[1]);
    string list_file(argv[2]);
    vector<string> pdbs, pdbfiles;
    vector<char> chains;
    assert( read_target_files(pdb_dir, list_file, pdbfiles, pdbs, chains) );
    int n = pdbs.size();
    assert( pdbs.size() == chains.size() );
    assert( pdbfiles.size() == pdbs.size() );
    
    for( int i=0; i < n; ++i ){
        vector<ResidueAtomPosition> v = residue_positions_from_pdb( pdbfiles[i], chains[i]);
        for(int j=0; j < v.size(); ++j)
          cout << pdbs[i] << " " << chains[i] << " " << v[j] << endl;
    }
    return 0;
}