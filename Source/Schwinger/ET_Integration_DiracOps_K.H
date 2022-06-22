#ifndef ET_INTEGRATION_DIRACOPS_K_H
#define ET_INTEGRATION_DIRACOPS_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 


AMREX_GPU_DEVICE 
inline 
void 
state_set_g3p(int i, int j, int k,
        amrex::Array4<amrex::Real> const& g3p_fab,
        amrex::Array4<amrex::Real> const& p_fab,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
        const amrex::GeometryData& geom) noexcept 
{
        
        g3p_fab(i, j, k, cIdx::Real_0) = p_fab(i, j, k, cIdx::Real_0);
        g3p_fab(i, j, k, cIdx::Imaginary_0) = p_fab(i, j, k, cIdx::Imaginary_0);
        g3p_fab(i, j, k, cIdx::Real_1) = -p_fab(i, j, k, cIdx::Real_1);
        g3p_fab(i, j, k, cIdx::Imaginary_1) = -p_fab(i, j, k, cIdx::Imaginary_1);
        
}

AMREX_GPU_DEVICE 
inline 
void 
state_set_g1p(int i, int j, int k,
        amrex::Array4<amrex::Real> const& g1p_fab,
        amrex::Array4<amrex::Real> const& p_fab,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
        const amrex::GeometryData& geom) noexcept 
{
        
        g1p_fab(i, j, k, cIdx::Real_0) = p_fab(i, j, k, cIdx::Real_1);
        g1p_fab(i, j, k, cIdx::Imaginary_0) = p_fab(i, j, k, cIdx::Imaginary_1);
        g1p_fab(i, j, k, cIdx::Real_1) = p_fab(i, j, k, cIdx::Real_0);
        g1p_fab(i, j, k, cIdx::Imaginary_1) = p_fab(i, j, k, cIdx::Imaginary_0);
        
}

AMREX_GPU_DEVICE 
inline 
void 
state_set_g2p(int i, int j, int k,
        amrex::Array4<amrex::Real> const& g2p_fab,
        amrex::Array4<amrex::Real> const& p_fab,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
        const amrex::GeometryData& geom) noexcept 
{
        
        g2p_fab(i, j, k, cIdx::Real_0) = p_fab(i, j, k, cIdx::Imaginary_1);
        g2p_fab(i, j, k, cIdx::Imaginary_0) = -p_fab(i, j, k, cIdx::Real_1);
        g2p_fab(i, j, k, cIdx::Real_1) = -p_fab(i, j, k, cIdx::Imaginary_0);
        g2p_fab(i, j, k, cIdx::Imaginary_1) = p_fab(i, j, k, cIdx::Real_0);
        
}

#endif