#ifndef ET_INTEGRATION_DIAG_K_H
#define ET_INTEGRATION_DIAG_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
state_diagnostics(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& diag, 
        amrex::Array4<amrex::Real const> const& state_fab,
        const amrex::Real time_lev,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx, 
        const amrex::GeometryData& geom) noexcept 
{
        const auto domain_xlo = geom.ProbLo(); 
        
        amrex::Real x0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1];
    
        amrex::Real dx0 = geom.CellSize(0);
        amrex::Real dx1 = geom.CellSize(1);
        
        /*amrex::Real U1_r = std::cos(state_fab(i, j, k, Idx::Theta_0));
        amrex::Real U1_i = std::sin(state_fab(i, j, k, Idx::Theta_0));
        
        amrex::Real U2_r = std::cos(state_fab(i+1, j, k, Idx::Theta_1));
        amrex::Real U2_i = std::sin(state_fab(i+1, j, k, Idx::Theta_1));
        
        amrex::Real U3_r = std::cos(state_fab(i, j+1, k, Idx::Theta_0));
        amrex::Real U3_i = -std::sin(state_fab(i, j+1, k, Idx::Theta_0));
        
        amrex::Real U4_r = std::cos(state_fab(i, j, k, Idx::Theta_1));
        amrex::Real U4_i = -std::sin(state_fab(i, j, k, Idx::Theta_1));
    
        amrex::Real plaq_re = U1_i*U2_i*U3_i*U4_i - U1_i*U2_i*U3_r*U4_r - U1_i*U2_r*U3_i*U4_r - U1_i*U2_r*U3_r*U4_i - U1_r*U2_i*U3_i*U4_r - U1_r*U2_i*U3_r*U4_i - U1_r*U2_r*U3_i*U4_i + U1_r*U2_r*U3_r*U4_r;
        
        diag(i, j, k, Diag::ActionDensity) = 1-plaq_re;*/
        
        
}
#endif