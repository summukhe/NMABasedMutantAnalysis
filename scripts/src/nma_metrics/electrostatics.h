# ifndef __ELECTRO_STATICS_H
# define __ELECTRO_STATICS_H

# include <cmath>
# include <vector>
# include <string>
# include <cassert>
# include "grid.h"
# include "utility.h"
# include "geometry.h"
# include "trajectory.h"
# include "reconstruction.h"
# include "partial_charge.h"
using namespace std;


# ifndef M_PI
# define M_PI 3.14159f
# endif

# ifndef EPSILON
# define EPSILON 1.f
# endif

# ifndef EPS
# define EPS 1e-3
# endif


struct PointCharge{
	Coordinate c;
	float      q;
	float      p;
	
public:
	PointCharge( Coordinate const& pt, 
				 const float charge = 0.f, 
				 const float prob = 1.f):c(pt),q(charge),p(prob)
	{
		assert( this->p >= 0.f && this->p <= 1.f);
	}
};


float epot( const float, 
            const float,
            const float = EPSILON );

float epot( Coordinate const& ,
            std::vector<PointCharge> const&,
			const float = EPSILON );


float epot( const float q, const float r, const float epsilon)
{
	float R = ( r > EPS )?r:EPS;
	return q/(4*M_PI*epsilon*R);
}


float epot( Coordinate const& pt, std::vector<PointCharge> const& charges, const float epsilon )
{
	float val, pot = 0.f;
	int n = charges.size();
	#pragma omp parallel for private(val) reduction(+:pot)
	for( int i=0; i < n; ++i )
	{
		val = epot(charges[i].q, 
				   ::euclidean_distance(pt, charges[i].c), 
					epsilon) * charges[i].p;
		pot += val;
	}
	return pot;
}


GridSignature epotential_signature(Snapshot const& inst, Grid3D const& box, const float epsilon = EPSILON)
{
	int n = box.cell_count();
	GridSignature score;
	std::vector<Coordinate> coordinates;
	for(int i=0; i < n; ++i )
		coordinates.push_back(box.get_coordinate(i));
	CASequence seq = convert_snapshot2caseq(inst);
	std::vector<int> resids = seq.residue_ids(true);
	std::vector<PointCharge> vpts;
	for( int i=0; i < resids.size(); ++i )
	{
		CaReconstructionUnit const& config = seq.residue_config( resids[i] );
		std::vector<string> anames = config.atom_names();
		for( int j=0; j < anames.size(); ++j )
		{
			int m = config.nconfigurations(anames[j]);
			for( int k=0; k < m; ++k )
			{
				Configuration cf = config.configuration(anames[j], k);
				vpts.push_back( PointCharge(cf.coordinate, cf.charge, cf.prob) );
			}
		}
	}	
	for(int i=0; i < n; ++i )
		score[i] = epot(coordinates[i], vpts, epsilon);
	return score;
}


# endif