#ifndef ET_INTEGRATION_AMR_CELL_TAG_K_H
#define ET_INTEGRATION_AMR_CELL_TAG_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
bool
state_is_tagged(int i, int j, int k,
        int lev,
        amrex::Array4<amrex::Real const> const& state_fab,
        amrex::Real error_threshold,
		const amrex::Real time,
		const amrex::Real time_lower, const amrex::Real time_upper,
		const int x0_lower, const int x0_upper,
		const int x1_lower, const int x1_upper,		
        const amrex::GeometryData& geom) noexcept 
{

        const auto domain_xlo = geom.ProbLo(); 

        amrex::Real x0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1]; 
        amrex::Real x2 = (k + 0.5)*geom.CellSize(2) + domain_xlo[2]; 
        
        amrex::Real x0c = 0;
        amrex::Real x1c = 0;
        
        amrex::Real x0c2 = 0;
        amrex::Real x1c2 = 0;
        
        
        //return !(((x0 - x0c >= -16 && x0 - x0c < 16) && (x1 - x1c >= -16 && x1 - x1c < 16)) || ((x0 - x0c2 >= -16 && x0 - x0c2 < 16) && (x1 - x1c2 >= -16 && x1 - x1c2 < 16))) || (time < 500.0 || time > 600.0);
        
        //return ((x0-x0c)*(x0-x0c) + (x1-x1c)*(x1-x1c) > 1024) || (time < 500.0 || time > 600.0);
        
        return (!((x0 >= x0_lower && x0 < x0_upper) && (x1 >= x1_lower && x1< x1_upper))) || (time < time_lower || time > time_upper);
        //return true;
        
        //(time < 100 || time > 200) || (x1 < -35 || x1 > 35 || x0 < -35 || x0 > 35);
        
        
        // || ((i-32)*(i-32) + (j-32)*(j-32) > 400); //(time < 100.0 || time > 200.0); // || ((i-32)*(i-32) + (j-32)*(j-32) > 400); //(i > 8 && i < 24);//(time < 100.0) || ((i < 5 || i > 27) || (j < 5 || j > 27));  //(i == 0); // && (((i == 8) && (j == 8)) || ((i == 24) && (j == 24))); //(time < 10.0) || ((time > 20.0) && (i == 32) && (j == 16));
}
#endif
