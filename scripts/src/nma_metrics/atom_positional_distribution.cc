# include <map>
# include <string>
# include <vector>
# include <unistd.h>
# include <cstdio>
# include "pdb.h"
# include "utility.h"
# include "reconstruction.h"
# include "partial_charge.h"
using namespace std;

void usage(char const* prog)
{
	fprintf(stderr, "Usage: %s -p <PDB file> -c <Chain ID> [-h] \n", prog );
	exit(1);
}


int main( int argc, char ** argv){
	int c;
	char chain = '\0';
	string pdbfile;
	string prog(argv[0]);
	while(( c = getopt(argc, argv, "p:c:")) != -1 ){
		switch(c)
		{
			case 'p':
			   pdbfile = optarg;
			   if( ! is_file(pdbfile) )
				   usage(prog.c_str());
			   break;
			case 'c':
			   chain = optarg[0];
			   break;
			default:
			   usage(prog.c_str());
			   break;
		}
	}
	
	if( pdbfile == "" || chain == '\0')
		usage(prog.c_str());
	
	vector<ResidueAtomPosition> data = residue_positions_from_pdb(pdbfile, chain, __degree);
	for( auto it = data.begin(); it != data.end(); ++it ){
		printf("%s %s %s %s %s %.2f %.2f %.2f\n", it->residue_pre.c_str(), 
		                                       it->residue_name.c_str(),
							 			       it->residue_post.c_str(),
											   it->atom_name.c_str(),
											   ::lookup_conversion(it->calpha_angle).c_str(),
											   it->atom_distance,
											   it->pseudo_angle,
											   it->pseudo_dihed);
	}
	return 0;
}