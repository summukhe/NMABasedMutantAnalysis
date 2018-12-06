# ifndef __PARTIAL_CHARGE_H
# define __PARTIAL_CHARGE_H

# include <fstream>
# include <vector>
# include <string>
# include <cstdlib>
# include <iomanip>
# include <iostream>
# include "pdb.h"
# include "utility.h"
# include "geometry.h"
# include "coordinate.h"
using namespace std;

# ifndef LINE_WIDTH
# define LINE_WIDTH 100
# endif

# ifndef M_PI
# define M_PI 3.14159265
# endif

# ifndef INVALID_INDEX
# define INVALID_INDEX -1
# endif

typedef map<string,float>              AtomChargeLookup;
typedef map<string,AtomChargeLookup>   ResidueChargeLookup;

struct ResidueAtomPosition{
    string residue_pre;
    string residue_name;
    string residue_post;
    string atom_name;
    float  atom_distance;
    float  calpha_angle;
    float  pseudo_angle;
    float  pseudo_dihed;
public:
    ResidueAtomPosition( string const& pre,
                         string const& rname, 
                         string const& post,
                         string const& aname,
                         const float dist,
                         const float cangle,
                         const float pangle,
                         const float pdihed):residue_pre(pre),
                                             residue_name(rname),
                                             residue_post(post),
                                             atom_name(aname),
                                             atom_distance(dist),
                                             calpha_angle(cangle),
                                             pseudo_angle(pangle),
                                             pseudo_dihed(pdihed)
    {}
};


ResidueChargeLookup read_charge_lookup( string const& );

Coordinate place_pseudo_sidechain( Coordinate const& , 
                                   Coordinate const& , 
                                   Coordinate const& ,
                                   const float = 1.f );


ostream& operator << (ostream& , 
                      ResidueAtomPosition const& );

vector<ResidueAtomPosition> residue_positions_from_pdb( string const& , 
                                                        const char ,
                                                        const AngleFormat = __degree);


ResidueChargeLookup  read_charge_lookup( string const& filename )
{
    assert( is_file( filename ) );
    ifstream f(filename.c_str());
    string sep=",";
    char line[LINE_WIDTH];
    bzero(line, LINE_WIDTH);
    ResidueChargeLookup lookup;
    
    while( f.getline(line, LINE_WIDTH)){
        vector<string> flds = string_split(string(line), sep);
        if( flds.size() == 4 ){
            string resname = flds[0]; 
            string atomname = flds[1];
            float  charge = atof(flds[3].c_str());
            if( ! is_in(lookup, resname) )
                lookup[resname] = AtomChargeLookup();
            lookup[resname][atomname] = charge;
        }
    }
    f.close();
    return lookup;
}

Coordinate place_pseudo_sidechain( Coordinate const& ca_1m, 
                                   Coordinate const& ca_0, 
                                   Coordinate const& ca_1p,
                                   const float dist )
{
    assert( dist > 0. );
    Vector3D v1 = connecting_vector(ca_1m, ca_0).unit_vector();
    Vector3D v2 = connecting_vector(ca_1p, ca_0).unit_vector();
    Vector3D r  = (v1 + v2).unit_vector();
    return Coordinate( ca_0[0] + dist * r[0], ca_0[1] + dist * r[1], ca_0[2] + dist * r[2]);
}
    


ostream& operator << (ostream& os, ResidueAtomPosition const& pos){
    os << pos.residue_pre << " "
       << pos.residue_name << " "
       << pos.residue_post << " "
       << pos.atom_name <<" "
       << setw(6) << setprecision(2)<< fixed << pos.calpha_angle << " "
       << setw(6) << setprecision(2)<< fixed << pos.atom_distance << " "
       << setw(6) << setprecision(2)<< fixed << pos.pseudo_angle << " "
       << setw(6) << setprecision(2)<< fixed << pos.pseudo_dihed;
    return os;
}
    

vector<ResidueAtomPosition> residue_positions_from_pdb( string const& filename, 
                                                        const char chain,
                                                        const AngleFormat angtype ){
    PDBChain molecule = read_pdb(filename, chain);
    vector<ResidueAtomPosition> results;
    int n = chain_size(molecule);
    for(int i=1; i < n-1; ++i ){
        if( (molecule.residues[i].residue_id - molecule.residues[i-1].residue_id == 1 ) &&
            (molecule.residues[i+1].residue_id - molecule.residues[i].residue_id == 1) && 
            (search_atom(molecule.residues[i-1],"CA") != INVALID_INDEX) && 
            (search_atom(molecule.residues[i], "CA") != INVALID_INDEX) && 
            (search_atom(molecule.residues[i+1], "CA") != INVALID_INDEX) )
        {
            Coordinate c1 = get_coordinate(molecule.residues[i-1], "CA");
            Coordinate c2 = get_coordinate(molecule.residues[i], "CA");
            Coordinate c3 = get_coordinate(molecule.residues[i+1], "CA");
            
            Coordinate pseudo_atom = place_pseudo_sidechain(c1,c2,c3);
            float ang = calculate_angle(c1, c2, c3, angtype);
            int m = residue_size( molecule.residues[i] );
            for(int j=0; j < m; ++j ){
                if( molecule.residues[i].atom_names[j] == "CA")
                    continue;
                float dist = euclidean_distance(pseudo_atom, molecule.residues[i].coordinates[j]);
                float dihed = dihedral_angle(c1, c2, pseudo_atom, 
                                             molecule.residues[i].coordinates[j], angtype);
                float pangle = calculate_angle(c2, pseudo_atom, molecule.residues[i].coordinates[j], angtype);
                results.push_back(ResidueAtomPosition(molecule.residues[i-1].residue_name,
                                                      molecule.residues[i].residue_name,
                                                      molecule.residues[i+1].residue_name,
                                                      molecule.residues[i].atom_names[j],
                                                      dist,
                                                      ang,
                                                      pangle,
                                                      dihed));
            }
        }
    }
    return results;
}

# endif