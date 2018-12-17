# ifndef __TRAJECTORY_ANALYSIS_H
# define __TRAJECTORY_ANALYSIS_H

# include <set>
# include <cmath>
# include <cassert>
# include <algorithm>
# include "grid.h"
# include "octree.h"
# include "geometry.h" 
# include "trajectory.h"
# include "electrostatics.h"
# include "residue_volume.h"
# include "sasa_calculation.h"
using namespace std;

# ifndef EPSILON
# define EPSILON 1.f
# endif

typedef std::vector<GridSignature> TrajectorySignature;

enum Approximation {  _grid = 1, _octree = 2};

void fill_sphere_approximation( Grid3D&,
                                Coordinate const&,
                                const float,
								const float);

void fill_sphere_in_grid( Grid3D& , 
                          Coordinate const& , 
                          const float , 
                          const float ,
                          const Approximation = _grid);

void fill_grid_fast(Snapshot const&, 
                    Grid3D&,
					const float);

void fill_grid( Snapshot const& ,
				Grid3D& ,
                const float ,
                const Approximation = _grid);

bool get_bounding_box( Trajectory const&, 
                       std::vector<int> const&, 
                       Coordinate&, 
                       Coordinate&);


bool is_box( Coordinate const& , 
             Coordinate const&  );
             
  
float box_volume( Coordinate const& , 
                  Coordinate const& );  


bool is_inside( Coordinate const& , 
                Coordinate const& , 
                Coordinate const& );
         
int box_intersect( Coordinate const& , 
                   Coordinate const& ,
                   Coordinate const& , 
                   Coordinate const& );
                   
bool inclusive( Coordinate const& , 
                Coordinate const& ,
                Coordinate const& , 
                Coordinate const& );
                
bool union_box( Coordinate const& , 
                Coordinate const& ,
                Coordinate const& , 
                Coordinate const& ,
                Coordinate& ,  
                Coordinate& );
                   
TrajectorySignature generate_volume_signature( Trajectory const& , 
											   Coordinate const& ,
                                               Coordinate const& , 
                                               const float , 
                                               const float );

TrajectorySignature generate_electrostatic_signature(Trajectory const&, 
													 Coordinate const&,
													 Coordinate const&,
													 const float,
													 const float = EPSILON);

GridSignature occlusion_map(TrajectorySignature const&);    

float occlusion_score( TrajectorySignature const& , 
					   TrajectorySignature const& );

bool is_intersecting( Coordinate const&, 
	                  Coordinate const&, 
                      Coordinate const&, 
					  Coordinate const&);
					  
std::vector<float> trajectory_signature_difference( TrajectorySignature const& ,  
                                              TrajectorySignature const& ,
                                              const int = 0);
					  
float grid_signature_difference( GridSignature const&, 
                                 GridSignature const& );
					

TrajectorySignature sasa_map(Trajectory const& , 
								 Grid3D const& );


std::vector<float> sasa_score( Trajectory const& , Trajectory const&, 
                               Grid3D const& , const int = 0);





void fill_sphere_approximation( Grid3D& box,
                                Coordinate const& center,
                                const float radius,
	   					        const float precision )
{
    assert( (precision < radius) && (precision > 0) );
    float xmx = center.x + radius, ymx = center.y + radius, zmx = center.z + radius;
    float xmn = center.x - radius, ymn = center.y - radius, zmn = center.z - radius;
    Grid3D sphere(xmx,ymx,zmx,
	              xmn,ymn,zmn,
				  precision, precision, precision);
				  
    if( is_intersecting( sphere.bounding_box(), box.bounding_box() ) )
	{
	    int n =  sphere.cell_count();
	    float v = sphere.cell_volume();
        for( int i=0; i < n; ++i)
		{
          Coordinate crd = sphere.get_coordinate(i);
		  if( euclidean_distance(center, crd) < (radius + 0.5 * precision) )
		  {
            sphere.fill_volume(crd, v);
		  }
		}
	    for(int i=0; i < sphere.cell_count(); ++i )
		{
          Coordinate crd = sphere.get_coordinate(i);
          box.fill_volume(crd, sphere.cell_volume(i));
		}
	}
}


void fill_sphere_in_grid( Grid3D& grid, 
                          Coordinate const& center, 
                          const float radius, 
                          const float precision,
                          const Approximation approx_type)
{
    fill_sphere_approximation( grid, center, radius, precision);
}


void fill_grid_fast(Snapshot const& snapshot, 
                    Grid3D& grid,
					const float precision)
{
	std::vector<Sphere3D> rep;
	int n = snapshot_size(snapshot);
	std::vector<AABB>     vbox;
	std::vector<Sphere3D> vsph;
	AABB box = grid.bounding_box();
	int m = grid.cell_count();
	for( int i=0; i < m; ++i){
		AABB cell = grid.grid_cell_by_index(i);
		vbox.push_back(cell);
		vsph.push_back(::inscribed_sphere(cell));
	}
	assert( n > 1 );
	Coordinate c1, c2, c3;
	for(int i=0; i < n; ++i){
		c3 = ::fetch_coordinate(snapshot, snapshot.residueIds[i]);
		Sphere3D sph(c3, CARadius);
		if( ::is_intersecting(box, sph) )
			rep.push_back(sph);
		if(i > 2){
		  SideChainModel sm = fix_sidechain(c1, c2, c3, snapshot.residueType[i-1]);
		  Sphere3D sc(sm.center(), sm.radius());
		  if( ::is_intersecting(box, sc))
			  rep.push_back(sc);
		}
		c1 = c2;
		c2 = c3;
	}
	
	std::map<int,float> vol_frac;
	for(int i=0; i < rep.size(); ++i)
	  for(int j=0; j < vbox.size(); ++j)
		 if( ::is_intersecting(vbox[j], rep[i]) )
		 {
           float f = (intersecting_volume(rep[i],vsph[j]) * vbox[j].volume())/vsph[j].volume();
		   if( ! ::is_in(vol_frac, j) )
			   vol_frac[j] = 0;
		   vol_frac[j] += f ;
		 }
	for(typename std::map<int,float>::const_iterator it=vol_frac.begin(); it != vol_frac.end(); ++it)
		grid.add_to_cell(it->first, it->second);
}					


void fill_grid( Snapshot const& snapshot,   /* Snapshot containing all calpha positions */ 
                Grid3D& grid,  
                const float precision,      /* finer volume calculation smaller grid size */
                const Approximation approx_type )
{
   float granularity = grid.granularity();
   assert( precision < granularity );
   assert( precision > 0.f );
   int n = snapshot_size(snapshot);   

   fill_sphere_in_grid( grid, 
                        ::fetch_coordinate(snapshot, snapshot.residueIds[0]), 
                        CARadius,
                        precision, 
                        approx_type );
   fill_sphere_in_grid( grid, 
                        ::fetch_coordinate(snapshot, snapshot.residueIds[n-1]), 
                        CARadius, 
                        precision, 
                        approx_type );

   for( int i=1; i < n-1; ++i ){
      int residue0 = snapshot.residueIds[i-1];
      int residue1 = snapshot.residueIds[i];
      int residue2 = snapshot.residueIds[i+1];

      Coordinate crd0 = ::fetch_coordinate(snapshot, residue0); 
      Coordinate crd1 = ::fetch_coordinate(snapshot, residue1);
      Coordinate crd2 = ::fetch_coordinate(snapshot, residue2);

      SideChainModel sm = fix_sidechain(crd0, crd1, crd2, snapshot.residueType[i]);

      fill_sphere_in_grid(grid, 
                          crd1, 
                          CARadius, 
                          precision, 
                          approx_type);
      fill_sphere_in_grid(grid, 
                          sm.center(), 
                          sm.radius(), 
                          precision, 
                          approx_type );
   }
}


bool get_bounding_box( Trajectory const& trajectory, 
                       std::vector<int> const& residues, 
                       Coordinate& max_point, 
                       Coordinate& min_point )
{
   float mx, my, mz;  mx = my = mz = -9999.f;
   float nx, ny, nz;  nx = ny = nz =  9999.f;
   bool status = true;

   for( int i=0; i < trajectory_size(trajectory); ++i ){
     for( int j=0; j < residues.size(); ++j ){
        if( residue_position(trajectory[i], residues[j]) != -1 ){
           Coordinate crd = fetch_coordinate(trajectory[i], residues[j] );
           mx = (mx < crd.x)?crd.x:mx;
           my = (my < crd.y)?crd.y:my;
           mz = (mz < crd.z)?crd.z:mz;

           nx = (nx > crd.x)?crd.x:nx;
           ny = (ny > crd.y)?crd.y:ny;
           nz = (nz > crd.z)?crd.z:nz;
        }else{
           status = false;
        }
     }
   }
   max_point.x = mx; max_point.y = my; max_point.z = mz;
   min_point.x = nx; min_point.y = ny; min_point.z = ny;
   return status;
}


bool is_box( Coordinate const& mx_pt, Coordinate const& mn_pt ) {
   return (mx_pt.x > mn_pt.x) && (mx_pt.y > mn_pt.y) && (mx_pt.z > mn_pt.z);
}


float box_volume( Coordinate const& max_point, 
                  Coordinate const& min_point ) 
{
   assert( is_box(max_point, min_point) );
   return (max_point.x - min_point.x)*(max_point.y - min_point.y)*(max_point.z - min_point.z);
}


bool is_inside( Coordinate const& mx_point, 
                Coordinate const& mn_point, 
                Coordinate const& crd )
{
   if( ! is_box(mx_point, mn_point) ) return false;
   return ( mx_point.x > crd.x && crd.x > mn_point.x ) && 
          ( mx_point.y > crd.y && crd.y > mn_point.y ) && 
          ( mx_point.z > crd.z && crd.z > mn_point.z );
}


int box_intersect( Coordinate const& max_point1, 
                   Coordinate const& min_point1,
                   Coordinate const& max_point2, 
                   Coordinate const& min_point2 )
{
   Coordinate max_point_ref = max_point1, 
              min_point_ref = min_point1,
              max_point_tgt = max_point2,
              min_point_tgt = min_point2;

   if( box_volume(max_point1, min_point1) < box_volume(max_point2, min_point2) )
   {
      max_point_ref = max_point2;
      min_point_ref = min_point2;
      max_point_tgt = max_point1;
      min_point_tgt = min_point1;
   }
   
   int s = 0;
   if( is_inside( max_point_ref, min_point_ref, min_point_tgt ) ) s++;
   if( is_inside(max_point_ref, min_point_ref, Coordinate(min_point_tgt.x, min_point_tgt.y, max_point_tgt.z)) ) s++;
   if( is_inside(max_point_ref, min_point_ref, Coordinate(min_point_tgt.x, max_point_tgt.y, min_point_tgt.z)) ) s++;
   if( is_inside(max_point_ref, min_point_ref, Coordinate(min_point_tgt.x, max_point_tgt.y, max_point_tgt.z)) ) s++;
   if( is_inside(max_point_ref, min_point_ref, Coordinate(max_point_tgt.x, min_point_tgt.y, min_point_tgt.z)) ) s++;
   if( is_inside(max_point_ref, min_point_ref, Coordinate(max_point_tgt.x, min_point_tgt.y, max_point_tgt.z)) ) s++;
   if( is_inside(max_point_ref, min_point_ref, Coordinate(max_point_tgt.x, max_point_tgt.y, min_point_tgt.z)) ) s++;
   if( is_inside(max_point_ref, min_point_ref, max_point_tgt ) ) s++;
   
   return s;
}


bool is_intersecting( Coordinate const& max_point1, 
                      Coordinate const& min_point1, 
                      Coordinate const& max_point2, 
                      Coordinate const& min_point2 )
{
  return (box_intersect(max_point1, min_point1, max_point2, min_point2) > 0);
}


bool inclusive( Coordinate const& max_point1, 
                Coordinate const& min_point1,
                Coordinate const& max_point2, 
                Coordinate const& min_point2)
{
   return (box_intersect(max_point1, min_point1, max_point2, min_point2) == 8);
}


bool union_box( Coordinate const& max_point1, 
                Coordinate const& min_point1,
                Coordinate const& max_point2, 
                Coordinate const& min_point2,
                Coordinate& max_point,  
                Coordinate& min_point )
{
  assert( is_box(max_point1, min_point1) && is_box(max_point2, min_point2) );
  bool status = ( ! is_intersecting( max_point1, min_point1, max_point2, min_point2) );

  max_point.x = max( max_point1.x, max_point2.x );
  max_point.y = max( max_point1.y, max_point2.y );
  max_point.z = max( max_point1.z, max_point2.z );
  
  min_point.x = min( min_point1.x, min_point2.x );
  min_point.y = min( min_point1.y, min_point2.y );
  min_point.z = min( min_point1.z, min_point2.z );
  return status;
}


TrajectorySignature generate_volume_signature( Trajectory const& trajectory, 
                                               Coordinate const& max_box,
                                               Coordinate const& min_box, 
                                               const float granularity, 
                                               const float precision )
{
   TrajectorySignature vsig;
   std::vector<Grid3D> vgrid;
   for(int i=0; i < trajectory_size(trajectory); ++i ){
	      vgrid.push_back( Grid3D( max_box.x, max_box.y, max_box.z, 
				                   min_box.x, min_box.y, min_box.z, 
				                   granularity, granularity, granularity ) );
				  
   }
		
   #pragma omp parallel for
   for(int i=0; i < trajectory_size(trajectory); ++i){
	 fill_grid_fast(trajectory[i], vgrid[i], precision);
   }
   for(int i=0; i < trajectory_size(trajectory); ++i)
      vsig.push_back( vgrid[i].volume_signature() );
   return vsig;
}


TrajectorySignature generate_electrostatic_signature(Trajectory const& trajectory, 
													 Coordinate const& max_box,
													 Coordinate const& min_box,
													 const float granularity,
													 const float epsilon )
{
	TrajectorySignature vsig;
	int n = ::trajectory_size(trajectory);
	Grid3D grid( max_box.x,   max_box.y,   max_box.z, 
	             min_box.x,   min_box.y,   min_box.z, 
				 granularity, granularity, granularity);
	for( int i=0; i < n; ++i ){
		vsig.push_back( epotential_signature(trajectory[i], grid, epsilon) );
	}
	return vsig;
}


float grid_signature_difference( GridSignature const& sig1, 
                                 GridSignature const& sig2 )
{
   float score, s1, s2;
   score = s1 = s2 = 0.f;
   std::set<int> grid_index; 
   typename GridSignature::const_iterator it;
   for( it = sig1.begin(); it != sig1.end(); ++it )
        grid_index.insert(it->first);
   for( it = sig2.begin(); it != sig2.end(); ++it )
        grid_index.insert(it->first);
   for( typename std::set<int>::const_iterator i = grid_index.begin(); i != grid_index.end(); ++i )
   {
      s1 = ( (it = sig1.find(*i)) == sig1.end() )?0.f:it->second;
      s2 = ( (it = sig2.find(*i)) == sig2.end() )?0.f:it->second;
      score += fabs(s1 - s2);
   }
   return score;
}


std::vector<float> trajectory_signature_difference( TrajectorySignature const& trj_sig1,  
                                              TrajectorySignature const& trj_sig2,
                                              const int theta )
{
  assert( trj_sig1.size() == trj_sig2.size() );
  std::vector<float> v;
  int n = trj_sig1.size();
  for( int i=0; i < n; ++i )
     v.push_back( grid_signature_difference(trj_sig1[i], trj_sig2[ (i + theta) % n ]) );
  return v;
}


GridSignature occlusion_map( TrajectorySignature const& trj_sig )
{
  std::set<int> grid_index;
  GridSignature omap;
  typename GridSignature::iterator iter;
  typename GridSignature::const_iterator const_iter;
  for( int i=0; i < trj_sig.size(); ++i )
  {
   for( const_iter = trj_sig[i].begin(); const_iter != trj_sig[i].end(); ++const_iter )
   {
      if( (iter = omap.find(const_iter->first) ) == omap.end() )
          omap[const_iter->first] = const_iter->second;
      else if( iter->second < const_iter->second )
          iter->second = const_iter->second;
   }
  }
  return omap;
}


TrajectorySignature sasa_map(Trajectory const& trj, Grid3D const& grid )
{
	TrajectorySignature trj_sig;
	int n = ::trajectory_size(trj);
	for( int i=0; i < n; ++i )
	{
		int nresidues = ::snapshot_size( trj[i] );
		std::vector<int> residue_ids;
		for( int j=0; j < nresidues; ++j )
		{
			if( grid.inside_grid( trj[i].coordinates[j] ) )
				residue_ids.push_back( trj[i].residueIds[j] );
		}
		GridSignature smap = ::sasa_values(trj[i], residue_ids);
		trj_sig.push_back(smap);
	}
	return trj_sig;
}


float occlusion_score( TrajectorySignature const& trj_ref, 
					   TrajectorySignature const& trj_tgt)
{
	assert( trj_ref.size() == trj_tgt.size() );
	return fabs( grid_signature_difference(occlusion_map(trj_ref),
	                                       occlusion_map(trj_tgt)));
}

std::vector<float> sasa_score( Trajectory const& trj_ref, 
							   Trajectory const& trj_tgt,
				               Grid3D const& box,
                               const int theta )
{
	assert( ::trajectory_size(trj_ref) == ::trajectory_size(trj_tgt) );
	return ::trajectory_signature_difference( sasa_map(trj_ref, box), 
	                                          sasa_map(trj_tgt, box), theta );
}

# endif
