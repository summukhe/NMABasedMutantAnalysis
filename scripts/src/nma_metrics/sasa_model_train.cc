# include <string>
# include <unistd.h>
# include <iomanip>
# include <cstdlib>

# include "utility.h"
# include "sasa_calculation.h"
using namespace std;

void usage(char const* prog)
{
	fprintf(stderr, "Usage: %s -i <sasa library> -m <sasa model location> [-h] \n", prog );
	exit(1);
}


int main(int argc, char** argv)
{
   int c;
   string prog(argv[0]);
   string data_infile, model_outfile;
   while(( c = getopt(argc, argv, "i:m:h")) != -1 ){
	 switch(c){
	   case 'i':
	     data_infile = optarg;
	     if( ! is_file(data_infile)) 
			 usage(prog.c_str());
		 break;
	   case 'm':
	     model_outfile = optarg;
		 break;
	   default:
	     usage(prog.c_str());
		 break;
	 }
   }
   if( data_infile == "" || model_outfile == "" )
	   usage(prog.c_str());
   RegressorInput inp = ::read_dnn_regressor_input(data_infile);
   float mae = train_network(inp, model_outfile);
   cout << "Test set error: " << mae << endl;
   return 0;
}