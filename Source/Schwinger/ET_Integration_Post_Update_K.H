#ifndef ET_INTEGRATION_POST_UPDATE_K_H
#define ET_INTEGRATION_POST_UPDATE_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
state_post_update(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& state_fab,
        const amrex::Real time,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
        const amrex::GeometryData& geom,
        const amrex::Real Temp_T,
        amrex::Real epsilon) noexcept 
{
        const auto domain_xlo = geom.ProbLo();
        
        amrex::Real U_0_Real = state_fab(i, j, k, Idx::U_0_Real);
        amrex::Real U_0_Imaginary = state_fab(i, j, k, Idx::U_0_Imaginary);
        
        amrex::Real U_1_Real = state_fab(i, j, k, Idx::U_1_Real);
        amrex::Real U_1_Imaginary = state_fab(i, j, k, Idx::U_1_Imaginary);
        
        amrex::Real P_0 = state_fab(i, j, k, Idx::P_0);
        amrex::Real P_1 = state_fab(i, j, k, Idx::P_1);
        
        
        state_fab(i, j, k, Idx::U_0_Real) = std::cos(epsilon*P_0)*U_0_Real - std::sin(epsilon*P_0)*U_0_Imaginary;
        state_fab(i, j, k, Idx::U_0_Imaginary) = std::sin(epsilon*P_0)*U_0_Real + std::cos(epsilon*P_0)*U_0_Imaginary;
        
        state_fab(i, j, k, Idx::U_1_Real) = std::cos(epsilon*P_1)*U_1_Real - std::sin(epsilon*P_1)*U_1_Imaginary;
        state_fab(i, j, k, Idx::U_1_Imaginary) = std::sin(epsilon*P_1)*U_1_Real + std::cos(epsilon*P_1)*U_1_Imaginary;

}
#endif