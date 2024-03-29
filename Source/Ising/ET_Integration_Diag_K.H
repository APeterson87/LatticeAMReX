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
        amrex::Real x2 = (k + 0.5)*geom.CellSize(2) + domain_xlo[2]; 
        
        amrex::Real s = state_fab(i, j, k, Idx::s);
        amrex::Real sup = state_fab(i, j+1, k, Idx::s);
        amrex::Real sdown = state_fab(i, j-1, k, Idx::s);
        amrex::Real sleft = state_fab(i-1, j, k, Idx::s);
        amrex::Real sright = state_fab(i+1, j, k, Idx::s);
        
        diag(i, j, k, Diag::Edensity) = -s*(sup+sdown+sleft+sright);  
        

}
#endif