# include <map>
# include <string>
# include <vector>
# include <unistd.h>
# include <cstdlib>
# include <cstdio>
# include "config_reader.h"
using namespace std;

void usage(char const* prog)
{
	fprintf(stderr, "Usage: %s -c <base config file> [-h] \n", prog );
	exit(1);
}


int main(int argc, char** argv)
{
	int c;
	char chain = '\0';
	string config_file;
	string prog(argv[0]);
	while(( c = getopt(argc, argv, "c:h")) != -1 ){
		switch(c)
		{
			case 'c':
			   config_file = optarg;
			   if( ! is_file(config_file) )
				   usage(prog.c_str());
			   break;
			default:
			   usage(prog.c_str());
			   break;
		}
	}
	
	if( config_file.size() == 0  || ! is_file(config_file))
		usage(prog.c_str());
		
	BoundingBoxConfig config = ::read_boundingbox_configfile(config_file.c_str());
	cout << config << endl;
	return 0;
}