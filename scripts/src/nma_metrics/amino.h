# ifndef __AMINO_H
# define __AMINO_H

# include <cmath>
# include <string>
# include <cstring>
# include <cstdlib>
# include <stdexcept>
# include <iomanip>
# include <iostream>
# include "residue_volume.h"
using namespace std;

# ifndef NAMINOS
# define NAMINOS  20
# endif

# ifndef EPS
# define EPS 1e-5
# endif

# ifndef CACA_MAX_DISTANCE
# define CACA_MAX_DISTANCE 4.0f
# endif

# ifndef CACA_MIN_DISTANCE
# define CACA_MIN_DISTANCE 3.6f
# endif

struct AminoAcid{
    string name;
    vector<string> atoms;
    vector<float>  partial_charges;
    float          sidechain_radius;
};


static const string Amino3Letter[   ] = {"ALA", "CYS", "ASP", "GLU", 
                                         "PHE", "GLY", "HIS", "ILE", 
                                         "LYS", "LEU", "MET", "ASN", 
                                         "PRO", "GLN", "ARG", "SER", 
                                         "THR", "VAL", "TRP", "TYR" };
                                       
static const char Amino1Letter[ ] =  {'A', 'C', 'D', 'E', 
                                      'F', 'G', 'H', 'I',
                                      'K', 'L', 'M', 'N', 
                                      'P', 'Q', 'R', 'S', 
                                      'T', 'V', 'W', 'Y'};

static const float HydroPathyIndex[ ] = { 1.8,  2.5, -3.5, -3.5, 
                                          2.8, -0.4, -3.2,  4.5, 
									     -3.9,  3.8,  1.9, -3.5, 
										 -1.6, -3.5, -4.5, -0.8, 
										 -0.7,  4.2, -0.9, -1.3};  
										  
static const float ResidueVolume[ ] = { 88.6, 108.5, 111.1, 138.4, 
                                       189.9,  60.1, 153.2, 166.7, 
									   168.6, 166.7, 162.9, 114.1, 
									   112.7, 143.8, 173.4,  89.0, 
									   116.1, 140.0, 227.8, 193.6};
                                     
static const string ALA_AMINO[] = {"N", "H", "CA", "C", "O", "CB"};
static const string CYS_AMINO[] = {"N", "H", "CA", "C", "O", "CB", "SG"};
static const string ASP_AMINO[] = {"N", "H", "CA", "C", "O", "CB", "CG", "OD1", "OD2"};
static const string GLU_AMINO[] = {"N", "H", "CA", "C", "O", "CB", "CG", "CD",  "OE1", "OE2"};
static const string PHE_AMINO[] = {"N", "H", "CA", "C", "O", "CB", "CG", "CD1", "HD1", "CD2", "CE1", "HE1", "CE2", "HE2", "CZ"};
static const string GLY_AMINO[] = {"N", "H", "CA", "C", "O"};
static const string HIS_AMINO[] = {"N", "H", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"};
static const string ILE_AMINO[] = {"N", "H", "CA", "C", "O", "CB", "CG1", "CG2"};
static const string LYS_AMINO[] = {"N", "H", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ", "HZ1", "HZ2"};
static const string LEU_AMINO[] = {"N", "H", "CA", "C", "O", "CB", "CG", "CD1", "CD2"};
static const string MET_AMINO[] = {"N", "H", "CA", "C", "O", "CB", "CG", "SD", "CE"};
static const string ASN_AMINO[] = {"N", "H", "CA", "C", "O", "CB", "CG", "OD1", "ND2", "HD21", "HD22"};
static const string PRO_AMINO[] = {"N", "CA", "C", "O", "CB", "CG", "CD"};
static const string GLN_AMINO[] = {"N", "H", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2", "HE21", "HE22"};
static const string ARG_AMINO[] = {"N", "H", "CA", "C", "O", "CB", "CG", "CD", "CG", "CZ", "NE", "HE", "NH1", "HH11", "HH12", "NH2", "HH21", "HH22"};
static const string SER_AMINO[] = {"N", "H", "CA", "C", "O", "CB", "OG", "HG"};
static const string THR_AMINO[] = {"N", "H", "CA", "C", "O", "CB", "OG1", "HG1", "CG2"};
static const string VAL_AMINO[] = {"N", "H", "CA", "C", "O", "CB", "CG1", "CG2"};
static const string TRP_AMINO[] = {"N", "H", "CA", "C", "O", "CB", "CD1", "CD2", "CG", "NE1", "HE1", "CE2", "CE3", "CZ2", "CZ3", "CH2" };
static const string TYR_AMINO[] = {"N", "H", "CA", "C", "O", "CB", "CD1", "CD2", "CG", "CE1", "CE2", "CZ", "OH", "HH"};


static const float ALA_PARTIAL[] = {-0.280,  0.280,  0.000,  0.380, -0.380,  0.000};
static const float CYS_PARTIAL[] = {-0.280,  0.280,  0.000,  0.380, -0.380, -0.100, -0.400};
static const float ASP_PARTIAL[] = {-0.280,  0.280,  0.000,  0.380, -0.380,  0.000,  0.270, -0.635, -0.635};
static const float GLU_PARTIAL[] = {-0.280,  0.280,  0.000,  0.380, -0.380,  0.000,  0.000,  0.270, -0.635, -0.635};
static const float PHE_PARTIAL[] = {-0.280,  0.280,  0.000,  0.380, -0.380,  0.000,  0.000, -0.100,  0.100, -0.100, -0.100,  0.100, -0.100,  0.100, -0.100};
static const float GLY_PARTIAL[] = {-0.280,  0.280,  0.000,  0.380, -0.380};
static const float HIS_PARTIAL[] = {-0.280,  0.280,  0.000,  0.380, -0.380,  0.000,  0.000,  0.000,  0.130,  0.260, -0.580};
static const float ILE_PARTIAL[] = {-0.280,  0.280,  0.000,  0.380, -0.380,  0.000,  0.000,  0.000};
static const float LYS_PARTIAL[] = {-0.280,  0.280,  0.000,  0.380, -0.380,  0.000,  0.000,  0.000,  0.000, -0.830,  0.415,  0.415};
static const float LEU_PARTIAL[] = {-0.280,  0.280,  0.000,  0.380, -0.380,  0.000,  0.000,  0.000,  0.000};
static const float MET_PARTIAL[] = {-0.280,  0.280,  0.000,  0.380, -0.380,  0.000,  0.000,  0.000,  0.000};
static const float ASN_PARTIAL[] = {-0.280,  0.280,  0.000,  0.380, -0.380,  0.000,  0.380, -0.380, -0.830,  0.415,  0.415};
static const float PRO_PARTIAL[] = { 0.000,  0.000,  0.380, -0.380,  0.000,  0.000,  0.000};
static const float GLN_PARTIAL[] = {-0.280,  0.280,  0.000,  0.380, -0.380,  0.000,  0.000,  0.380, -0.380, -0.830,  0.415,  0.415};
static const float ARG_PARTIAL[] = {-0.280,  0.280,  0.000,  0.380, -0.380,  0.000,  0.000,  0.090,  0.000,  0.340, -0.110,  0.240, -0.260,  0.240,  0.240, -0.260,  0.240,  0.240};
static const float SER_PARTIAL[] = {-0.280,  0.280,  0.000,  0.380, -0.380,  0.150, -0.548,  0.398};
static const float THR_PARTIAL[] = {-0.280,  0.280,  0.000,  0.380, -0.380,  0.150, -0.548,  0.398,  0.000};
static const float VAL_PARTIAL[] = {-0.280,  0.280,  0.000,  0.380, -0.380,  0.000,  0.000,  0.000};
static const float TRP_PARTIAL[] = {-0.280,  0.280,  0.000,  0.380, -0.380,  0.000, -0.100,  0.000, -0.140, -0.050,  0.190,  0.000, -0.100, -0.100, -0.100, -0.100};
static const float TYR_PARTIAL[] = {-0.280,  0.280,  0.000,  0.380, -0.380,  0.000, -0.100, -0.100,  0.000, -0.100, -0.100,  0.150, -0.548,  0.398};


string name(AminoAcid const& a){
    return a.name;
}

int natoms(AminoAcid const& a){
    return static_cast<int>(a.atoms.size());
}

ostream& operator << (ostream& os, AminoAcid const& a ){
    os << a.name << ": [" << natoms(a) << " , " << fixed << setprecision(3) << a.sidechain_radius << "]";
    return os;
}

bool is_amino(const string& a){
    if( a.size() != 3 ) return false;
    for( int i=0; i < NAMINOS; ++i )
        if( Amino3Letter[i] == a )
            return true;
    return false;
}

bool is_amino( const char* a) {
    string s(a);
    return is_amino(s);
}

bool is_amino( const char a ){
    for(int i=0; i < NAMINOS; ++i)
        if( Amino1Letter[i] == a )
            return true;
    return false;
}

string amino_abbr( const char a){
    for(int i=0; i < NAMINOS; ++i)
        if( Amino1Letter[i] == a )
            return Amino3Letter[i];
    return string("UNK");
}

char amino_abbr( const string& a){
    for(int i=0; i < NAMINOS; ++i)
        if( Amino3Letter[i] == a )
            return Amino1Letter[i];
    return 'X';
}

AminoAcid get_amino(const char c){
    assert( is_amino(c) );
    AminoAcid amino;
    amino.name = amino_abbr(c);
    int m,n;
    switch( c ){
        case 'A':
           n = sizeof(ALA_AMINO)/sizeof(string);
           m = sizeof(ALA_PARTIAL)/sizeof(float);
           assert( n == m );
           for(int i=0; i < n; ++i ){
               amino.atoms.push_back(ALA_AMINO[i]);
               amino.partial_charges.push_back(ALA_PARTIAL[i]);
           }
           break;
        case 'C':
           n = sizeof(CYS_AMINO)/sizeof(string);
           m = sizeof(CYS_PARTIAL)/sizeof(float);
           assert( n == m );
           for(int i=0; i < n; ++i ){
               amino.atoms.push_back(CYS_AMINO[i]);
               amino.partial_charges.push_back(CYS_PARTIAL[i]);
           }
           break;
        case 'D':
           n = sizeof(ASP_AMINO)/sizeof(string);
           m = sizeof(ASP_PARTIAL)/sizeof(float);
           assert( n == m );
           for(int i=0; i < n; ++i ){
               amino.atoms.push_back(ASP_AMINO[i]);
               amino.partial_charges.push_back(ASP_PARTIAL[i]);
           }
           break;
        case 'E':
           n = sizeof(GLU_AMINO)/sizeof(string);
           m = sizeof(GLU_PARTIAL)/sizeof(float);
           assert( n == m );
           for(int i=0; i < n; ++i ){
               amino.atoms.push_back(GLU_AMINO[i]);
               amino.partial_charges.push_back(GLU_PARTIAL[i]);
           }
           break;
        case 'F':
           n = sizeof(PHE_AMINO)/sizeof(string);
           m = sizeof(PHE_PARTIAL)/sizeof(float);
           assert( n == m );
           for(int i=0; i < n; ++i ){
               amino.atoms.push_back(PHE_AMINO[i]);
               amino.partial_charges.push_back(PHE_PARTIAL[i]);
           }
           break;
        case 'G':
           n = sizeof(GLY_AMINO)/sizeof(string);
           m = sizeof(GLY_PARTIAL)/sizeof(float);
           assert( n == m );
           for(int i=0; i < n; ++i ){
               amino.atoms.push_back(GLY_AMINO[i]);
               amino.partial_charges.push_back(GLY_PARTIAL[i]);
           }
           break;
        case 'H':
           n = sizeof(HIS_AMINO)/sizeof(string);
           m = sizeof(HIS_PARTIAL)/sizeof(float);
           assert( n == m );
           for(int i=0; i < n; ++i ){
               amino.atoms.push_back(HIS_AMINO[i]);
               amino.partial_charges.push_back(HIS_PARTIAL[i]);
           }
           break;
        case 'I':
           n = sizeof(ILE_AMINO)/sizeof(string);
           m = sizeof(ILE_PARTIAL)/sizeof(float);
           assert( n == m );
           for(int i=0; i < n; ++i ){
               amino.atoms.push_back(ILE_AMINO[i]);
               amino.partial_charges.push_back(ILE_PARTIAL[i]);
           }
           break;
        case 'K':
           n = sizeof(LYS_AMINO)/sizeof(string);
           m = sizeof(LYS_PARTIAL)/sizeof(float);
           assert( n == m );
           for(int i=0; i < n; ++i ){
               amino.atoms.push_back(LYS_AMINO[i]);
               amino.partial_charges.push_back(LYS_PARTIAL[i]);
           }
           break;
        case 'L':
           n = sizeof(LEU_AMINO)/sizeof(string);
           m = sizeof(LEU_PARTIAL)/sizeof(float);
           assert( n == m );
           for(int i=0; i < n; ++i ){
               amino.atoms.push_back(LEU_AMINO[i]);
               amino.partial_charges.push_back(LEU_PARTIAL[i]);
           }
           break;
        case 'M':
           n = sizeof(MET_AMINO)/sizeof(string);
           m = sizeof(MET_PARTIAL)/sizeof(float);
           assert( n == m );
           for(int i=0; i < n; ++i ){
               amino.atoms.push_back(MET_AMINO[i]);
               amino.partial_charges.push_back(MET_PARTIAL[i]);
           }
           break;
        case 'N':
           n = sizeof(ASN_AMINO)/sizeof(string);
           m = sizeof(ASN_PARTIAL)/sizeof(float);
           assert( n == m );
           for(int i=0; i < n; ++i ){
               amino.atoms.push_back(ASN_AMINO[i]);
               amino.partial_charges.push_back(ASN_PARTIAL[i]);
           }
           break;
        case 'P':
           n = sizeof(PRO_AMINO)/sizeof(string);
           m = sizeof(PRO_PARTIAL)/sizeof(float);
           assert( n == m );
           for(int i=0; i < n; ++i ){
               amino.atoms.push_back(PRO_AMINO[i]);
               amino.partial_charges.push_back(PRO_PARTIAL[i]);
           }
           break;
        case 'Q':
           n = sizeof(GLN_AMINO)/sizeof(string);
           m = sizeof(GLN_PARTIAL)/sizeof(float);
           assert( n == m );
           for(int i=0; i < n; ++i ){
               amino.atoms.push_back(GLN_AMINO[i]);
               amino.partial_charges.push_back(GLN_PARTIAL[i]);
           }
           break;
        case 'R':
           n = sizeof(ARG_AMINO)/sizeof(string);
           m = sizeof(ARG_PARTIAL)/sizeof(float);
           assert( n == m );
           for(int i=0; i < n; ++i ){
               amino.atoms.push_back(ARG_AMINO[i]);
               amino.partial_charges.push_back(ARG_PARTIAL[i]);
           }
           break;
        case 'S':
           n = sizeof(SER_AMINO)/sizeof(string);
           m = sizeof(SER_PARTIAL)/sizeof(float);
           assert( n == m );
           for(int i=0; i < n; ++i ){
               amino.atoms.push_back(SER_AMINO[i]);
               amino.partial_charges.push_back(SER_PARTIAL[i]);
           }
           break;
        case 'T':
           n = sizeof(THR_AMINO)/sizeof(string);
           m = sizeof(THR_PARTIAL)/sizeof(float);
           assert( n == m );
           for(int i=0; i < n; ++i ){
               amino.atoms.push_back(THR_AMINO[i]);
               amino.partial_charges.push_back(THR_PARTIAL[i]);
           }
           break;
        case 'V':
           n = sizeof(VAL_AMINO)/sizeof(string);
           m = sizeof(VAL_PARTIAL)/sizeof(float);
           assert( n == m );
           for(int i=0; i < n; ++i ){
               amino.atoms.push_back(VAL_AMINO[i]);
               amino.partial_charges.push_back(VAL_PARTIAL[i]);
           }
           break;
        case 'W':
           n = sizeof(TRP_AMINO)/sizeof(string);
           m = sizeof(TRP_PARTIAL)/sizeof(float);
           assert( n == m );
           for(int i=0; i < n; ++i ){
               amino.atoms.push_back(TRP_AMINO[i]);
               amino.partial_charges.push_back(TRP_PARTIAL[i]);
           }
           break;
        case 'Y':
           n = sizeof(TYR_AMINO)/sizeof(string);
           m = sizeof(TYR_PARTIAL)/sizeof(float);
           assert( n == m );
           for(int i=0; i < n; ++i ){
               amino.atoms.push_back(TYR_AMINO[i]);
               amino.partial_charges.push_back(TYR_PARTIAL[i]);
           }
           break;
        default:
           throw "Error: Invalid Amino Requested!!";
    }
    amino.sidechain_radius = get_radius(amino.name);
    return amino;
}


AminoAcid get_amino(string const& a)
{
    return get_amino(amino_abbr(a));
}


vector<string> charged_atoms( AminoAcid const& a )
{
    int n = natoms(a);
    vector<string> catoms;
    for(int i=0; i < n; ++i )
        if( abs(a.partial_charges[i]) > EPS )
            catoms.push_back(a.atoms[i]);
    return catoms;
}

float get_volume( string const& res )
{
	int n = sizeof(Amino3Letter)/sizeof(string);
	for( int i=0; i < n; ++i )
		if( res == Amino3Letter[i] )
			return ResidueVolume[i];
	return 0.f;
}


float get_volume( AminoAcid const& aa )
{
	return ::get_volume(::name(aa));
}


float get_hydropathy( string const& res )
{
	int n = sizeof(Amino3Letter)/sizeof(string);
	for( int i=0; i < n; ++i )
		if( res == Amino3Letter[i] )
			return HydroPathyIndex[i];
	return 0.f;
}

float get_hydropathy( AminoAcid const& aa )
{
	return ::get_hydropathy(::name(aa));
}


# endif