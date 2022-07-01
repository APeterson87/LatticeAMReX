#ifndef ET_INTEGRATION_INIT_K_H
#define ET_INTEGRATION_INIT_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
state_init(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& state_fab, 
        amrex::Real time, const amrex::GeometryData& geom) noexcept 
{
        const auto domain_xlo = geom.ProbLo(); 
        
        amrex::Real x0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1]; 

        
        amrex::Real theta_0 = 2*M_PI*(amrex::Random()-0.5);
        amrex::Real theta_1 = 2*M_PI*(amrex::Random()-0.5);
        
        amrex::Real P_0 = amrex::RandomNormal(0,1);
        amrex::Real P_1 = amrex::RandomNormal(0,1);
        
        state_fab(i, j, k, Idx::U_0_Real) = std::cos(theta_0); 
        state_fab(i, j, k, Idx::U_0_Imaginary) = std::sin(theta_0);
        state_fab(i, j, k, Idx::U_1_Real) = std::cos(theta_1);
        state_fab(i, j, k, Idx::U_1_Imaginary) = std::sin(theta_1);

        state_fab(i, j, k, Idx::P_0) = P_0;
        state_fab(i, j, k, Idx::P_1) = P_1;
        
        
}
#endif