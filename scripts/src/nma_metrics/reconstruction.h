# ifndef __RECONSTRUCTION_H
# define __RECONSTRUCTION_H

# include <map>
# include <vector>
# include <string>
# include <cstring>
# include <cassert>
# include <fstream>
# include <iomanip>
# include <iostream>
# include "amino.h"
# include "utility.h"
# include "geometry.h"
# include "coordinate.h"
# include "partial_charge.h"
using namespace std;

# ifndef LINE_WIDTH
# define LINE_WIDTH 100
# endif

# ifndef BIN_WIDTH
# define BIN_WIDTH 30
# endif


struct Configuration{
   Coordinate coordinate;
   float      prob;
   float      charge;
};

string lookup_conversion( const float angle){
    assert( angle >= 0 && angle <= 180 );
    int n = static_cast<int>(angle / BIN_WIDTH);
    char id[10];
    bzero(id, 10);
    sprintf(id, "N%d", n);
    return string(id);
}


struct ReconstructionUnit{
public:
    string lookup;
    float  distance;
    float  angle;
    float  dihedral;
    float  prob;
public:
    ReconstructionUnit(string const& n, 
                       const float d, 
                       const float a, 
                       const float t, 
                       const float p):lookup(n),distance(d),angle(a),dihedral(t),prob(p)
    {}
};

ostream& operator << (ostream& os, ReconstructionUnit const& ru ){
    os << "[" << fixed << setprecision(3) << ru.distance <<" , "
       << fixed << setprecision(3) << ru.angle  << " , "
       << fixed << setprecision(3) << ru.dihedral << "]";
    return os;
}

typedef vector<ReconstructionUnit>  BuildingBlocks;
typedef map<string, BuildingBlocks> AminoBlocks;


string base_directory()
{
    return "/home/sumanta/Work/codelite_base/nma_metrics/volume_fluctuation";
}

map<string, AminoBlocks> read_reconstruction_dictionary(){
    string dict = base_directory() + "/data/AtomModels.csv";
    assert( is_file(dict) );
    ifstream f(dict.c_str());
    char line[LINE_WIDTH];
    bzero(line, LINE_WIDTH);
    int counter = 0;
    map<string, AminoBlocks> library;
    AminoBlocks blocks;
    string lastName, lastAtom;
    while( f.getline(line, LINE_WIDTH) ){
        if(counter > 0){
            vector<string> flds = string_split(string(line), ",");
            if( flds.size() == 7 ){
                string name = flds[0];
                string atom = flds[1];
                string cls = flds[2];
                float dist = atof(flds[3].c_str());
                float angle = atof(flds[4].c_str());
                float dihed = atof(flds[5].c_str());
                float prob = atof(flds[6].c_str());
        
                if( ! is_in(library, name) )
                    library[name] = AminoBlocks();
                if( ! is_in(library[name], atom) )
                    library[name][atom] = BuildingBlocks();
                library[name][atom].push_back(ReconstructionUnit(cls, dist, angle, dihed, prob) );
                
            }
        }
        counter++;
    }
    return library;
}

class ReconstructionLibrary{
protected:
    map<string, AminoBlocks> library;
    static ReconstructionLibrary* m_inst;
private:
    ReconstructionLibrary()
    {
        this->library = read_reconstruction_dictionary();
    }
    
public:
    static ReconstructionLibrary const& get_instance()
    {
        if( ReconstructionLibrary::m_inst == 0 )
            ReconstructionLibrary::m_inst = new ReconstructionLibrary();
        return *ReconstructionLibrary::m_inst;
    }
    
    map<string, AminoBlocks> get_library()const
    {
        return this->library;
    }
};

Coordinate reconstruct_coordinate( Coordinate const& ak, 
                                   Coordinate const& aj, 
                                   Coordinate const& ai,
                                   const float dist,
                                   const float angle,
                                   const float dihed)
{
    Vector3D bij = connecting_vector(aj, ai).unit_vector();
	Vector3D nijk = cross( connecting_vector(aj, ai),
						   connecting_vector(ak, aj));
	Matrix3D dihedRot = rotation_matrix(bij, dihed);
	Matrix3D bondRot = rotation_matrix(nijk, -angle);
	Matrix3D m =  dihedRot * bondRot;
    Vector3D v = m * connecting_vector(aj, ai).unit_vector();
	return Coordinate( ai[0] + v[0]*dist, 
					   ai[1] + v[1]*dist, 
					   ai[2] + v[2]*dist);
}


class CaReconstructionUnit{

    typedef vector<Configuration>   Configurations;

    protected:
       AminoAcid                    amino;
       map<string, Configurations>  configs;
    
    private:
       CaReconstructionUnit();
     
    public:
        CaReconstructionUnit( Coordinate const& ca0, 
                              Coordinate const& ca1, 
                              Coordinate const& ca2,
                              string const& aa):amino( get_amino(aa) ),
							                    configs()
        {
            ReconstructionLibrary const& rl = ReconstructionLibrary::get_instance();
            float angle = calculate_angle(ca0, ca1, ca2, __radian);
            Coordinate c = place_pseudo_sidechain(ca0, ca1, ca2);
            string lp = lookup_conversion(angle);
            int n = ::natoms(this->amino);
            for(int i=0; i < n; ++i){
                string aname = this->amino.atoms[i];
            }
        }
        
        CaReconstructionUnit( CaReconstructionUnit const& ru):amino(ru.amino), configs(ru.configs)
        {
        }
     
        string  name()const 
        { 
            return ::name(this->amino); 
        }
        
        int  natoms()const
        {
            return ::natoms(this->amino);
        }
        
        vector<string>  atom_names()const
        {
            return this->amino.atoms;
        }
        
        vector<float> atom_charges()const
        {
            return this->amino.partial_charges;
        }
        
        vector<string>  charged_atoms()const
        {
            return ::charged_atoms(this->amino);
        }
        
        bool  is_atom( string const& a )const
        {
            return is_in(this->configs, a);
        }
        
        int   nconfigurations( string const& a )const
        {
            if( is_in(this->configs, a) ){
				return this->configs.find(a)->second.size();
			}
            return 0;
        }
        
        Configuration configuration( string const& a, const int i )const
        {
            assert( this->is_atom(a) );
            assert( i >= 0 &&  i < this->nconfigurations(a));
            return this->configs.find(a)->second[i] ;
        }
};


# endif