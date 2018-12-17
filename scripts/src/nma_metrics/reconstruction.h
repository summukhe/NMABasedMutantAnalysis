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
# include "trajectory.h"
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
public:
   Configuration( Coordinate const& c, const float q, const float p=1.f):coordinate(c),prob(p),charge(q)
   {
	   assert( prob >= 0.f && prob <= 1.f);
   }
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

typedef std::vector<ReconstructionUnit>  BuildingBlocks;
typedef std::map<string, BuildingBlocks> AminoBlocks;


std::map<string, AminoBlocks> read_reconstruction_dictionary(){
    string dict = base_directory() + "/data/AtomModels.csv";
    assert( is_file(dict) );
    ifstream f(dict.c_str());
    char line[LINE_WIDTH];
    bzero(line, LINE_WIDTH);
    int counter = 0;
    std::map<string, AminoBlocks> library;
    AminoBlocks blocks;
    string lastName, lastAtom;
    while( f.getline(line, LINE_WIDTH) ){
        if(counter > 0){
            std::vector<string> flds = string_split(string(line), ",");
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
    std::map<string, AminoBlocks> m_lib;
    static ReconstructionLibrary* m_inst;
private:
    ReconstructionLibrary()
    {
        this->m_lib = read_reconstruction_dictionary();
    }
    
public:
    static ReconstructionLibrary const& get_instance()
    {
        if( ReconstructionLibrary::m_inst == 0 )
            ReconstructionLibrary::m_inst = new ReconstructionLibrary();
        return *ReconstructionLibrary::m_inst;
    }
    
    int size()const
    {
        return static_cast<int>(this->m_lib.size());
    }
	
	std::vector<string> residues()const
	{
		return keys(this->m_lib);
	}
	
	int natoms(string const& n)const
	{
		if( is_in(this->m_lib, n))
			return static_cast<int>( this->m_lib.find(n)->second.size() );
		return 0;
	}
	
	std::vector<string> atom_list(string const& n)const
	{
		std::vector<string> v;
		if(is_in(this->m_lib, n))
			v = keys( this->m_lib.find(n)->second );
		return v;
	}
	
	int nconfigurations(string const& rname, string const& aname )const
	{
		if( is_in(this->m_lib, rname) )
			if( is_in( this->m_lib.find(rname)->second, aname) )
				return static_cast<int>(this->m_lib.find(rname)->second.find(aname)->second.size());
		return 0;
	}
	
	std::vector<ReconstructionUnit> configurations(string const& rname, string const& aname)const
	{
		std::vector<ReconstructionUnit> v;
		if( is_in(this->m_lib, rname) )
			if( is_in( this->m_lib.find(rname)->second, aname) )
				v = this->m_lib.find(rname)->second.find(aname)->second;
		return v;
	}
};

ReconstructionLibrary* ReconstructionLibrary::m_inst = 0;


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


class CaReconstructionUnit
{
    typedef std::vector<Configuration>   Configurations;

    protected:
       AminoAcid                    amino;
       std::map<string, Configurations>  configs;
    
    public:
        CaReconstructionUnit():amino(::get_amino("ALA")),configs()
		{}
		
		CaReconstructionUnit( CaReconstructionUnit const& ru):amino(ru.amino),
															  configs(ru.configs)
        {}
		
		CaReconstructionUnit(Coordinate const& ca0, 
                             Coordinate const& ca1, 
                             Coordinate const& ca2,
                             string const& aa):amino(::get_amino(aa)),configs()
		{
			this->build(ca0, ca1, ca2, aa);
		}

		bool build(Coordinate const& ca0, 
                   Coordinate const& ca1, 
                   Coordinate const& ca2,
                   string const& aa)
        {
			this->amino = ::get_amino(aa);
			this->configs.clear();
            ReconstructionLibrary const& rl = ReconstructionLibrary::get_instance();
            float angle = calculate_angle(ca0, ca1, ca2, __degree);
            Coordinate c = place_pseudo_sidechain(ca0, ca1, ca2);
			std::vector<string> atoms = rl.atom_list(aa);
            string lp = lookup_conversion(angle);
			PartialChargeLibrary const& ql = PartialChargeLibrary::get_instance();
            int n = ::natoms(this->amino);
            for(int i=0; i < n; ++i)
			{
                string aname = this->amino.atoms[i];
				Configurations pos;
				if( aname == "CA" )
				{
					pos.push_back(Configuration(ca1,ql.get_charge(aa,aname)));
				}else
				{
					std::vector<ReconstructionUnit> ru = rl.configurations(aa, aname); 
					for( int j = 0; j < ru.size(); ++j )
					{
						if( ru[j].lookup == lp )
						{
						  Coordinate crd = reconstruct_coordinate(ca0, ca1, c, 
						                                          ru[j].distance, 
																  DEG2RAD(ru[j].angle), 
																  DEG2RAD(ru[j].dihedral) );
						  pos.push_back(Configuration(crd, 
						                              ql.get_charge(aa,aname), 
													  ru[j].prob));
						}
					 
					}
				}
				configs[aname] = pos;
            }
        }
     
        string  name()const 
        { 
            return ::name(this->amino); 
        }
        
        int  natoms()const
        {
            return ::natoms(this->amino);
        }
        
        std::vector<string>  atom_names()const
        {
            return this->amino.atoms;
        }
        
        std::vector<float> atom_charges()const
        {
            return this->amino.partial_charges;
        }
        
        std::vector<string>  charged_atoms()const
        {
            return ::charged_atoms(this->amino);
        }
        
        bool is_atom( string const& a )const
        {
            return is_in(this->configs, a);
        }
        
        int nconfigurations( string const& a )const
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


class CASequence{
    protected:
	   std::vector<int>                    m_resids;
	   std::map<int,AminoAcid>             m_residues;
	   std::map<int,Coordinate>            m_ca_pos;
	   std::map<int,CaReconstructionUnit>  m_reconstruct;
    public:
	   CASequence():m_resids(),m_residues(),m_ca_pos(),m_reconstruct()
	   {}
	   
	   int size()const
	   {
		   return static_cast<int>(m_resids.size());
	   }
	   
	   bool add_residue(const int resid, string const& resname, Coordinate const& crd)
	   {
		   int n = this->size() - 1;
		   if( is_in(this->m_residues, resid) ) 
		       return false;
			   
		   if( n > 0 && resid - this->m_resids.back() != 1 )
			   return false;
			   
		   if( ! ::is_amino(resname) )
			   return false;
		   
		   if( n > 1 ){
			   float d = ::euclidean_distance(this->m_ca_pos[this->m_resids[n]], 
			                                this->m_ca_pos[this->m_resids[n-1]]);
											
			/*								
			   if( d < CACA_MIN_DISTANCE ||  d > CACA_MAX_DISTANCE )
				   cerr << "Error: CA("<< this->m_resids[n-1]<< ") - CA("<< this->m_resids[n] <<") distance : " << d << endl;
			 */  
		   }
			   
		   this->m_resids.push_back(resid);
		   this->m_residues[resid] = ::get_amino(resname);
		   this->m_ca_pos[resid] = crd;
		   if( this->size() > 2 )
		   {
			   n = this->size() - 1;
			   this->m_reconstruct[this->m_resids[n-2]] = CaReconstructionUnit();
			   
			   this->m_reconstruct[this->m_resids[n-2]].build(this->m_ca_pos[this->m_resids[n-3]],
			                                                  this->m_ca_pos[this->m_resids[n-2]],
														      this->m_ca_pos[this->m_resids[n-1]],   
															  ::name(this->m_residues[this->m_resids[n-2]]));
		   }
	   }
	   
	   std::vector<int> residue_ids( bool reconstruct_only = false )const 
	   {
		   std::vector<int> vResidues;
		   
		   if( reconstruct_only )
			   vResidues = ::keys( this->m_reconstruct );
		   else
			   vResidues = this->m_resids;
		   return vResidues;
	   }
		
	   AminoAcid const& get_amino( const int resid )const
	   {
		   assert( ::is_in(this->m_residues, resid) );
		   return this->m_residues.find(resid)->second;
	   }
		
	   Coordinate get_ca_coordinate( const int resid )const
	   {
		   assert( ::is_in(this->m_ca_pos, resid) );
		   return this->m_ca_pos.find(resid)->second;
	   }
	   
	   std::vector<int> amino_positions( const string& aa )const
	   {
		   std::vector<int> hits;
		   for( auto it = this->m_residues.begin(); it != this->m_residues.end(); ++it )
			 if( aa == ::name(it->second) )
			    hits.push_back(it->first);
		   return hits; 
	   }
	   
	   int count_aminos( const string& aa )const
	   {
		   return static_cast<int>( this->amino_positions(aa).size() );
	   }
	   

	   CaReconstructionUnit const& residue_config(const int resid)const
	   {
		   assert( ::is_in( this->m_reconstruct , resid) );
		   return this->m_reconstruct.find(resid)->second;
	   }
};


CASequence convert_snapshot2caseq( Snapshot const& s )
{
	int n = snapshot_size(s);
	CASequence seq;
	for( int i=0; i < n; ++i )
	  seq.add_residue(s.residueIds[i], 
	                  s.residueType[i], 
					  s.coordinates[i]);
	return seq;
}

# endif