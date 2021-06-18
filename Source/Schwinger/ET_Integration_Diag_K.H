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
        
        amrex::Real theta_0_00 = std::acos(state_fab(i, j, k, Idx::U_0_Re));
        amrex::Real theta_1_00 = std::acos(state_fab(i, j, k, Idx::U_1_Re));
        amrex::Real theta_0_01 = std::acos(state_fab(i, j+1, k, Idx::U_0_Re));
        amrex::Real theta_1_10 = std::acos(state_fab(i+1, j, k, Idx::U_1_Re));
        
        amrex::Real theta_plaq = theta_0_00 + theta_1_10 - theta_0_01 - theta_1_00;
        
        diag(i, j, k, Diag::ActionDensity) = 1-std::cos(theta_plaq);
        
        amrex::Real AbsU_0_sqrd = state_fab(i, j, k, Idx::U_0_Re)*state_fab(i, j, k, Idx::U_0_Re);
        
        amrex::Real AbsU_1_sqrd = state_fab(i, j, k, Idx::U_1_Re)*state_fab(i, j, k, Idx::U_1_Re);
        
        diag(i, j, k, Diag::AbsU_0_Re) = std::sqrt(AbsU_0_sqrd);
        diag(i, j, k, Diag::AbsU_1_Re) = std::sqrt(AbsU_1_sqrd);

}
#endif