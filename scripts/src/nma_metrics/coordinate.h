# ifndef __COORDINATE_H
# define __COORDINATE_H

# include <cmath>
# include <iomanip>
# include <iostream>
using namespace std;

struct Coordinate{
  float x, y, z;
  public:
    Coordinate(const float _x = 0.f,
               const float _y = 0.f,
               const float _z = 0.f ):x(_x),
                                      y(_y),
                                      z(_z)
    {}
};

float euclidean_distance( Coordinate const& crd1, Coordinate const& crd2 ){
  return static_cast<float>( sqrt( (crd1.x - crd2.x) * (crd1.x - crd2.x) +
                                   (crd1.y - crd2.y) * (crd1.y - crd2.y) +
                                   (crd1.z - crd2.z) * (crd1.z - crd2.z) ) );
}


ostream& operator << (ostream& os, Coordinate const& crd ){
   os << fixed << setprecision(3) << setw(8) << crd.x << setw(8) << crd.y << setw(8) << crd.z ;
   return os;
}

# endif
