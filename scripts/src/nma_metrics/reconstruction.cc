# include <string>
# include <cstdlib>
# include <iostream>
# include <ctime>
# include "amino.h"
# include "reconstruction.h"
using namespace std;


int main(int argc, char** argv){
	//srand(static_cast<unsigned int>(time(0)));
	
	
    Coordinate x1(24.969, 13.428, 30.692), 
			   x2(24.044, 12.661, 29.808), 
			   x3(22.785, 13.482, 29.543), 
			   x4(21.951, 13.670, 30.431);
	
	/*		   
	Coordinate x1(0.f, 0.f, 0.f), 
	           x2(1.f, 0.f, 0.f), 
			   x3(1.f, 1.f, 0.f), 
			   x4(1.f, 1.f, 1.f);
	*/
	map<string, AminoBlocks> lib = read_reconstruction_dictionary();
	
	cout << "Size: " << lib.size() << endl;
	
	float dist = euclidean_distance(x3,x4);
	float ang = calculate_angle(x2,x3,x4);
	float dihed = dihedral_angle(x1,x2,x3,x4);
	
	cout << dist << " " << ang  << " " << dihed << endl;
		 
	Coordinate x4c = reconstruct_coordinate(x1, x2, x3, dist, DEG2RAD(ang), DEG2RAD(dihed));
	
	cout << "Reconstruction Error: " <<  euclidean_distance(x4, x4c) << endl;
	
	cout << rotation_matrix(Vector3D(1,0,0), M_PI/2) << endl;
	
	return 0;
}