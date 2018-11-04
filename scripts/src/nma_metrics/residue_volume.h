# ifndef __RESIDUE_VOLUME_H
# define __RESIDUE_VOLUME_H

# include <cmath>
# include <string>
# include <cstring>
# include <cassert>
# include <cstdlib>
# include "coordinate.h"
using namespace std;

# ifndef EPS
# define EPS 1e-3
# endif

const string Aminos[] = {"ALA", "ARG", "ASN", "ASP", 
                         "CYS", "GLN", "GLU", "GLY", 
                         "HIS", "ILE", "LEU", "LYS", 
                         "MET", "PHE", "PRO", "SER", 
                         "THR", "TRP", "TYR", "VAL"};

const float  Radius[] = { 1.027, 3.794, 2.178, 1.865,
                          1.640, 2.626, 2.503, 0.485,
                          2.623, 2.283, 2.318, 3.477, 
                          2.885, 3.004, 1.605, 1.579, 
                          1.686, 3.491, 3.382, 1.685};


const float  CARadius = 1.7;

int naminos() { return static_cast<int>(sizeof(Aminos)/sizeof(string)); }

int amino_index( string const& resname ){
   int n = naminos();
   for(int i=0; i < n; ++i)
      if(Aminos[i] == resname )
         return i;
   return -1;
}

float get_radius( string const& resname ){
   int n = amino_index(resname);
   return (n != -1)?Radius[n]:0;
}

struct Vector3D{
  float xi;
  float yi;
  float zi;
  public:
      Vector3D( const float x = 0.f, const float y = 0.f, const float z = 0.f):xi(x),yi(y),zi(z)
      {}
};

struct SideChainModel{
   Coordinate center;
   float R;   
   public: 
       SideChainModel(float x=0.f, float y=0.f, float z=0.f, float r=0.f):center(x,y,z),R(r)
       {}
};


float  vector_norm(Vector3D const& v ){
  return sqrt( (v.xi * v.xi) + (v.yi * v.yi) + (v.zi * v.zi) );
}

Vector3D unit_vector( Vector3D const& v ){
   float n = vector_norm(v);
   assert( n > EPS );
   Vector3D u(v.xi, v.yi, v.zi);
   u.xi /= n;
   u.yi /= n;
   u.zi /= n;
   return u;
}

Vector3D operator + ( Vector3D const& u, Vector3D const& v ){
   Vector3D w(u.xi + v.xi, u.yi + v.yi, u.zi + v.zi);
   return w;
}

Vector3D operator - ( Vector3D const& u, Vector3D const& v){
   Vector3D w(u.xi - v.xi, u.yi - v.yi, u.zi - v.zi);
   return w;
}

float dot( Vector3D const& u, Vector3D const& v){
   return u.xi * v.xi + u.yi * v.yi + u.zi * v.zi;
}

Vector3D operator * ( const float x, Vector3D const& v){
   Vector3D w(x*v.xi, x*v.yi, x*v.zi);
   return w;
}

Vector3D operator * (Vector3D const& v, const float x){
   Vector3D w(x*v.xi, x*v.yi, x*v.zi);
   return w;
}

SideChainModel fix_sidechain( Coordinate const& ca_0, Coordinate const& ca_1, Coordinate const& ca_2 , string const& resname ){
   float r = get_radius(resname);
   Vector3D u(ca_0.x, ca_0.y, ca_0.z), v(ca_1.x, ca_1.y, ca_1.z), w(ca_2.x, ca_2.y, ca_2.z);
   Vector3D uv = unit_vector(v - u);
   Vector3D wv = unit_vector(v - w);
   Vector3D ov = unit_vector(uv + wv);
   float x = ca_1.x + ov.xi * r;
   float y = ca_1.y + ov.yi * r;
   float z = ca_1.z + ov.zi * r;
   return SideChainModel(x,y,z,r); 
}



# endif
