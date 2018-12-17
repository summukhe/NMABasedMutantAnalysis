# include <map>
# include <vector>
# include <string>
# include <cassert>
# include <unistd.h>
# include <cstdlib>
# include "utility.h"
# include "trajectory.h"
# include "config_reader.h"
# include "sasa_calculation.h"
# include "trajectory_analysis.h"
using namespace std;

# ifndef GRID_SIZE
# define GRID_SIZE     1.5f
# endif

# ifndef VOL_DISCRETE
# define VOL_DISCRETE  0.3f
# endif

# ifndef MIN_GRID
# define MIN_GRID      1.0f
# endif


void usage(char const* prog)
{
	fprintf(stderr, "Usage: %s -s <source trajectory> -t <target trajectory> -b <box configuration> [-g <grid size>] [-h] \n", prog );
	exit(1);
}

int main(int argc, char** argv)
{
   int c;
   string prog(argv[0]);
   float grid_size = GRID_SIZE;
   string src_trjfile, tgt_trjfile, box_configfile;
   while(( c = getopt(argc, argv, "s:t:b:g:h")) != -1 ){
	 switch(c){
	   case 's':
	     src_trjfile = optarg;
	     if( ! is_file(src_trjfile)) 
			 usage(prog.c_str());
		 break;
	   case 't':
	     tgt_trjfile = optarg;
		 if( ! is_file(tgt_trjfile) )
			 usage(prog.c_str());
		 break;
	   case 'b':
	     box_configfile = optarg;
		 if( ! is_file(box_configfile) )
			 usage(prog.c_str());
		 break;
	   case 'g':
	     grid_size = atof(optarg);
		 if( grid_size < MIN_GRID )
			 usage(prog.c_str());
	     break;
	   default:
	     usage(prog.c_str());
		 break;
	 }
   }
   
   if( src_trjfile == "" || tgt_trjfile == "" || box_configfile == "" )
	   usage(prog.c_str());
   
   Trajectory trajectory_ref = ::read_trajectory(src_trjfile);
   Trajectory trajectory_tgt = ::read_trajectory(tgt_trjfile);
   
   Trajectory trajectory_rev = flip_trajectory(trajectory_tgt);
   
   if(trajectory_distance(trajectory_ref, trajectory_rev ) > 
         trajectory_distance(trajectory_ref, trajectory_tgt))
           trajectory_tgt = trajectory_rev;
   
   BoundingBoxConfig config = ::read_boundingbox_configfile(box_configfile.c_str());
   Grid3D grid(config.box_max, config.box_min, grid_size);
   std::vector<float>  vsig = ::sasa_score(trajectory_ref,trajectory_tgt,grid);
   cout << fixed << setprecision(5) << (::mean(vsig.begin(), vsig.end()) / grid.grid_volume()) * 1000.f <<endl;
   return 0;
}