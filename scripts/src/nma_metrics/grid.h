# ifndef __GRID_H
# define __GRID_H 1

# include <map>
# include <cmath>
# include <vector>
# include <cstring>
# include <cassert>
# include <iomanip>
# include <iostream>
# include <algorithm>
# include "geometry.h"
using namespace std;

# ifndef EPS
# define EPS 1e-5
# endif

# ifndef INVALID_INDEX
# define INVALID_INDEX -1
# endif



struct GridIndex
{
  int xi, yi, zi;
  public:
      GridIndex(const int x, 
                const int y, 
                const int z):xi(x),
                             yi(y),
                             zi(z) 
      {}
};


typedef std::map<int,float> GridSignature;

bool valid_grid_index( GridIndex const& );
bool check_bounding_coordinate(Coordinate const& , Coordinate const& );


bool valid_grid_index( GridIndex const& gidx )
{
   return (gidx.xi != INVALID_INDEX) && 
		  (gidx.yi != INVALID_INDEX) && 
		  (gidx.zi != INVALID_INDEX);
}


class Grid3D
{
  protected:
   std::map<int,float> grid_cells;
   float               x_max, y_max, z_max;
   float               x_min, y_min, z_min;
   float               dx, dy, dz;
   
  protected:
     Grid3D():grid_cells()
	 {}
   
  public:
     Grid3D(Coordinate const& max_pt, 
	        Coordinate const& min_pt, 
			const float hx):x_max(max_pt.x),y_max(max_pt.y),z_max(max_pt.z),
							x_min(min_pt.x),y_min(min_pt.y),z_min(min_pt.z),
							dx(hx),dy(hx),dz(hx),grid_cells()
     {}
	 
     Grid3D( const float xmx, 
             const float ymx, 
             const float zmx, 
             const float xmn, 
             const float ymn, 
             const float zmn,
             const float hx,
             const float hy,
             const float hz):x_max(xmx),
			                 y_max(ymx),
							 z_max(zmx),
							 x_min(xmn),
							 y_min(ymn),
							 z_min(zmn),
							 dx(hx),
							 dy(hy),
							 dz(hz),
							 grid_cells()
     {
        assert( (this->x_max > this->x_min) && (this->x_max - this->x_min > this->dx) );
        assert( (this->y_max > this->y_min) && (this->y_max - this->y_min > this->dy) );
        assert( (this->z_max > this->z_min) && (this->z_max - this->z_min > this->dz) );
        assert( this->dx > EPS && this->dy > EPS && this->dz > EPS );
		int n = this->cell_count();
		assert( n > 0);
	 }
	 
	 Grid3D(Grid3D const& grid):grid_cells(grid.grid_cells),
	                            x_max(grid.x_max),
								y_max(grid.y_max),
								z_max(grid.z_max),
								x_min(grid.x_min),
								y_min(grid.y_min),
								z_min(grid.z_min),
								dx(grid.dx),
								dy(grid.dy),
								dz(grid.dz)
	 {}
		
	 Grid3D& operator = (Grid3D const& rhs)
	 {
		 if(this == &rhs ) return *this;
		 this->x_max = rhs.x_max;
		 this->y_max = rhs.y_max;
		 this->z_max = rhs.z_max;
		 this->x_min = rhs.x_min;
		 this->y_min = rhs.y_min;
		 this->z_min = rhs.z_min;
		 this->dx = rhs.dx;
		 this->dy = rhs.dy;
		 this->dz = rhs.dz;
		 this->grid_cells = rhs.grid_cells;
	 }

     inline int nx()const 
	 { 
		 return static_cast<int>( ceil( (this->x_max - this->x_min)/ this->dx) ); 
	 }
	 
     inline int ny()const 
	 { 
		 return static_cast<int>( ceil( (this->y_max - this->y_min)/ this->dy) ); 
	 }
	 
     inline int nz()const 
	 { 
		 return static_cast<int>( ceil( (this->z_max - this->z_min)/ this->dz) ); 
	 }

     float grid_volume()const 
	 { 
		 
		 return static_cast<float>( (this->x_max - this->x_min) * 
                                    (this->y_max - this->y_min) * 
                                    (this->z_max - this->z_min) ); 
	 }
     
	 float granularity()const
	 {
		 return std::min( std::min(this->dx, this->dy), this->dz);
	 }
	 
	 float cell_volume(void)const 
	 { 
		 return static_cast<float>(this->dx * this->dy * this->dz); 
	 }
 
     float cell_volume(const int i)const 
	 {
        if(i < 0 || i > this->cell_count()) return 0;
        return static_cast<float>( this->cell_fraction(i) * this->cell_volume() );
     }
	 
	 float cell_fraction(const int i)const
	 {
		typename std::map<int,float>::const_iterator it;
	    if( (it = this->grid_cells.find(i)) != this->grid_cells.end()) 
			return it->second;
		return 0.f;
	 }
	 
	 bool add_to_cell(const int i, const float v)
	 {
		 if( v <= 0.f || i < 0 || i > this->cell_count() ) 
			 return false;
		 typename std::map<int,float>::iterator it;
		 float f = (v / this->cell_volume());
		 f = ( f > 1.f )?1.:f;
		 if( (it = this->grid_cells.find(i)) != this->grid_cells.end() ){
			 it->second += f;
			 if(it->second > 1.f) it->second = 1.f;
		 }else{
			 this->grid_cells[i] = f;
		 }
		 return true;
	 }
	 
	 bool free_from_cell(const int i, const float v)
	 {
		 if( v <= 0.f || i < 0 || i > this->cell_count() ) 
			 return false;
		 typename std::map<int,float>::iterator it;
		 float f = (v / this->cell_volume());
		 f = ( f > 1.f )?1.:f;
		 if( (it = this->grid_cells.find(i)) != this->grid_cells.end() ){
			 if( it->second <= f )
				 this->grid_cells.erase(it);
			 else
				 it->second -= f;
			 return true;
		 }
		 return false;
	 }
       
     float cell_volume(GridIndex const& gidx)const 
	 {
        return this->cell_volume( this->to_linear_index( gidx ) );
     }
    
     inline int  cell_count()const  
	 { 
		 return static_cast<int>(this->nx() * this->ny() * this->nz()); 
	 }

     bool inside_grid( Coordinate const& crd )const 
	 { 
         return ( (crd.x < this->x_max) && (crd.y < this->y_max) && (crd.z < this->z_max) && 
                  (crd.x > this->y_min) && (crd.y > this->y_min) && (crd.z > this->z_min) );
     }     
   
     bool inside_grid( const float x, const float y, const float z)const 
	 {
        return this->inside_grid( Coordinate(x,y,z) );
     }

     GridIndex grid_coordinate( Coordinate const& crd )const 
	 {
         GridIndex gidx(INVALID_INDEX,INVALID_INDEX,INVALID_INDEX);
         if( this->inside_grid(crd) ){
            gidx.xi = static_cast<int>( (crd.x - this->x_min)/this->dx );
            gidx.yi = static_cast<int>( (crd.y - this->y_min)/this->dy );
            gidx.zi = static_cast<int>( (crd.z - this->z_min)/this->dz ); 
         }
         return gidx;
     }
     
     GridIndex grid_coordinate( const float x, const float y, const float z)const 
	 {
       return this->grid_coordinate( Coordinate(x,y,z) ); 
     }
	 
	 AABB bounding_box()const
	 {
		 return AABB( Coordinate(this->x_max, this->y_max, this->z_max), 
		              Coordinate(this->x_min, this->y_min, this->z_min) );
	 }

     AABB grid_cell( Coordinate const& crd ) const 
	 {
        return  this->grid_cell_by_index( this->grid_coordinate(crd) );
     }
     
     AABB grid_cell_by_index( GridIndex const& gidx ) const 
	 {
        return AABB( Coordinate( (gidx.xi+1) * this->dx + this->x_min, 
                                 (gidx.yi+1) * this->dy + this->y_min, 
                                 (gidx.zi+1) * this->dz + this->z_min ),
                     Coordinate( gidx.xi * this->dx + this->x_min, 
                                 gidx.yi * this->dy + this->y_min, 
                                 gidx.zi * this->dz + this->z_min ) );
     }
     
     AABB grid_cell_by_index( const int i ) const 
	 {
          return this->grid_cell_by_index( this->to_grid_index(i) );
     }

     int to_linear_index( GridIndex const& gidx )const 
	 {
        int index = INVALID_INDEX;
        if( valid_grid_index(gidx) )
           index = gidx.xi * this->ny() * this->nz() + gidx.yi * this->nz() + gidx.zi;
        return index;
     }

     GridIndex to_grid_index( const int linear_index )const 
	 {
        GridIndex gidx(INVALID_INDEX, INVALID_INDEX, INVALID_INDEX);
        if( linear_index >= 0 && linear_index < this->cell_count() ){
           gidx.xi = static_cast<int>( linear_index / (this->nz() * this->ny() ) );
           gidx.yi = static_cast<int>( (linear_index % (this->ny() * this->nz())) / this->nz()  );
           gidx.zi = static_cast<int>( (linear_index %  (this->ny() * this->nz())) % this->nz() );
        }
        return gidx;
     }

     Coordinate get_coordinate( GridIndex const& gidx )const 
	 {
        assert( valid_grid_index(gidx) );
        Coordinate crd( (gidx.xi + 0.5) * this->dx + this->x_min, 
                        (gidx.yi + 0.5) * this->dy + this->y_min, 
                        (gidx.zi + 0.5) * this->dz + this->z_min );
        return crd;
     }

     Coordinate get_coordinate( const int linear_index )const
	 {
        return this->get_coordinate( this->to_grid_index(linear_index) );
     }

     bool fill_volume( Coordinate const& crd, const float vol )
	 {
         if( ! this->inside_grid( crd ) || vol < 0.f ) return false;
         int idx = this->to_linear_index( this->grid_coordinate( crd ) );
         return this->add_to_cell(idx,vol);
     }

     bool free_volume( Coordinate const& crd, const float vol )
	 {
        if( ! this->inside_grid(crd) || vol < 0.f ) return false;
        int idx = this->to_linear_index( this->grid_coordinate(crd) );
        return this->free_from_cell(idx,vol);
     }

     float filled_volume()const
	 {
         float f = 0;
         for(typename std::map<int,float>::const_iterator it = this->grid_cells.begin(); it != this->grid_cells.end(); ++it )
            f += it->second;
         return this->cell_volume() *  f;
     }

     Coordinate min_coordinate()const 
	 {
        return Coordinate( this->x_min, this->y_min, this->z_min );
     }

     Coordinate max_coordinate()const 
	 {
        return Coordinate( this->x_max, this->y_max, this->z_max );
     }
     
	 GridSignature volume_signature()const
	 {
		 return this->grid_cells;
	 }
	 
	 virtual ~Grid3D()
	 {
	 }
};
	 

bool check_bounding_coordinate(Coordinate const& mx_crd, Coordinate const& mn_crd)
{
	return ( (mx_crd.x > mn_crd.x) && (mx_crd.y > mn_crd.y) && (mx_crd.z > mn_crd.z) );
}


ostream& operator << (ostream& os, GridSignature const& sig)
{
   typename GridSignature::const_iterator it;
   for( it = sig.begin(); it != sig.end(); ++it )
       os << it->first << ":" << fixed << setprecision(3) << it->second << endl;
   return os;
}


# endif /* __GRID_H */
