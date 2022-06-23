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

AMREX_GPU_DEVICE 
inline 
void 
state_set_Dp(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& Dp_fab,
        amrex::Array4<amrex::Real> const& p_fab,
        amrex::Array4<amrex::Real> const& U_fab,
        amrex::Real m_0) noexcept 
{
        
        amrex::Real Dp_Real_0 = (m_0 + 2)*p_fab(i, j, k, cIdx::Real_0) + U_fab(i - 1, j, k, cIdx::Imaginary_0)*p_fab(i - 1, j, k, cIdx::Imaginary_0) + U_fab(i - 1, j, k, cIdx::Imaginary_0)*p_fab(i - 1, j, k, cIdx::Imaginary_1) + U_fab(i - 1, j, k, cIdx::Real_0)*p_fab(i - 1, j, k, cIdx::Real_0) + U_fab(i - 1, j, k, cIdx::Real_0)*p_fab(i - 1, j, k, cIdx::Real_1) + U_fab(i, j - 1, k, cIdx::Imaginary_1)*p_fab(i, j - 1, k, cIdx::Imaginary_0) - U_fab(i, j - 1, k, cIdx::Imaginary_1)*p_fab(i, j - 1, k, cIdx::Real_1) + U_fab(i, j - 1, k, cIdx::Real_1)*p_fab(i, j - 1, k, cIdx::Imaginary_1) + U_fab(i, j - 1, k, cIdx::Real_1)*p_fab(i, j - 1, k, cIdx::Real_0) - U_fab(i, j, k, cIdx::Imaginary_0)*p_fab(i + 1, j, k, cIdx::Imaginary_0) + U_fab(i, j, k, cIdx::Imaginary_0)*p_fab(i + 1, j, k, cIdx::Imaginary_1) - U_fab(i, j, k, cIdx::Imaginary_1)*p_fab(i, j + 1, k, cIdx::Imaginary_0) - U_fab(i, j, k, cIdx::Imaginary_1)*p_fab(i, j + 1, k, cIdx::Real_1) + U_fab(i, j, k, cIdx::Real_0)*p_fab(i + 1, j, k, cIdx::Real_0) - U_fab(i, j, k, cIdx::Real_0)*p_fab(i + 1, j, k, cIdx::Real_1) - U_fab(i, j, k, cIdx::Real_1)*p_fab(i, j + 1, k, cIdx::Imaginary_1) + U_fab(i, j, k, cIdx::Real_1)*p_fab(i, j + 1, k, cIdx::Real_0);
       
        amrex::Real Dp_Real_1 = (m_0 + 2)*p_fab(i, j, k, cIdx::Real_1) + U_fab(i - 1, j, k, cIdx::Imaginary_0)*p_fab(i - 1, j, k, cIdx::Imaginary_0) + U_fab(i - 1, j, k, cIdx::Imaginary_0)*p_fab(i - 1, j, k, cIdx::Imaginary_1) + U_fab(i - 1, j, k, cIdx::Real_0)*p_fab(i - 1, j, k, cIdx::Real_0) + U_fab(i - 1, j, k, cIdx::Real_0)*p_fab(i - 1, j, k, cIdx::Real_1) + U_fab(i, j - 1, k, cIdx::Imaginary_1)*p_fab(i, j - 1, k, cIdx::Imaginary_1) + U_fab(i, j - 1, k, cIdx::Imaginary_1)*p_fab(i, j - 1, k, cIdx::Real_0) - U_fab(i, j - 1, k, cIdx::Real_1)*p_fab(i, j - 1, k, cIdx::Imaginary_0) + U_fab(i, j - 1, k, cIdx::Real_1)*p_fab(i, j - 1, k, cIdx::Real_1) + U_fab(i, j, k, cIdx::Imaginary_0)*p_fab(i + 1, j, k, cIdx::Imaginary_0) - U_fab(i, j, k, cIdx::Imaginary_0)*p_fab(i + 1, j, k, cIdx::Imaginary_1) - U_fab(i, j, k, cIdx::Imaginary_1)*p_fab(i, j + 1, k, cIdx::Imaginary_1) + U_fab(i, j, k, cIdx::Imaginary_1)*p_fab(i, j + 1, k, cIdx::Real_0) - U_fab(i, j, k, cIdx::Real_0)*p_fab(i + 1, j, k, cIdx::Real_0) + U_fab(i, j, k, cIdx::Real_0)*p_fab(i + 1, j, k, cIdx::Real_1) + U_fab(i, j, k, cIdx::Real_1)*p_fab(i, j + 1, k, cIdx::Imaginary_0) + U_fab(i, j, k, cIdx::Real_1)*p_fab(i, j + 1, k, cIdx::Real_1);


        amrex::Real Dp_Imaginary_0 = (m_0 + 2)*p_fab(i, j, k, cIdx::Imaginary_0) - U_fab(i - 1, j, k, cIdx::Imaginary_0)*p_fab(i - 1, j, k, cIdx::Real_0) - U_fab(i - 1, j, k, cIdx::Imaginary_0)*p_fab(i - 1, j, k, cIdx::Real_1) + U_fab(i - 1, j, k, cIdx::Real_0)*p_fab(i - 1, j, k, cIdx::Imaginary_0) + U_fab(i - 1, j, k, cIdx::Real_0)*p_fab(i - 1, j, k, cIdx::Imaginary_1) - U_fab(i, j - 1, k, cIdx::Imaginary_1)*p_fab(i, j - 1, k, cIdx::Imaginary_1) - U_fab(i, j - 1, k, cIdx::Imaginary_1)*p_fab(i, j - 1, k, cIdx::Real_0) + U_fab(i, j - 1, k, cIdx::Real_1)*p_fab(i, j - 1, k, cIdx::Imaginary_0) - U_fab(i, j - 1, k, cIdx::Real_1)*p_fab(i, j - 1, k, cIdx::Real_1) + U_fab(i, j, k, cIdx::Imaginary_0)*p_fab(i + 1, j, k, cIdx::Real_0) - U_fab(i, j, k, cIdx::Imaginary_0)*p_fab(i + 1, j, k, cIdx::Real_1) - U_fab(i, j, k, cIdx::Imaginary_1)*p_fab(i, j + 1, k, cIdx::Imaginary_1) + U_fab(i, j, k, cIdx::Imaginary_1)*p_fab(i, j + 1, k, cIdx::Real_0) + U_fab(i, j, k, cIdx::Real_0)*p_fab(i + 1, j, k, cIdx::Imaginary_0) - U_fab(i, j, k, cIdx::Real_0)*p_fab(i + 1, j, k, cIdx::Imaginary_1) + U_fab(i, j, k, cIdx::Real_1)*p_fab(i, j + 1, k, cIdx::Imaginary_0) + U_fab(i, j, k, cIdx::Real_1)*p_fab(i, j + 1, k, cIdx::Real_1);
        
        amrex::Real Dp_Imaginary_1 = (m_0 + 2)*p_fab(i, j, k, cIdx::Imaginary_1) - U_fab(i - 1, j, k, cIdx::Imaginary_0)*p_fab(i - 1, j, k, cIdx::Real_0) - U_fab(i - 1, j, k, cIdx::Imaginary_0)*p_fab(i - 1, j, k, cIdx::Real_1) + U_fab(i - 1, j, k, cIdx::Real_0)*p_fab(i - 1, j, k, cIdx::Imaginary_0) + U_fab(i - 1, j, k, cIdx::Real_0)*p_fab(i - 1, j, k, cIdx::Imaginary_1) + U_fab(i, j - 1, k, cIdx::Imaginary_1)*p_fab(i, j - 1, k, cIdx::Imaginary_0) - U_fab(i, j - 1, k, cIdx::Imaginary_1)*p_fab(i, j - 1, k, cIdx::Real_1) + U_fab(i, j - 1, k, cIdx::Real_1)*p_fab(i, j - 1, k, cIdx::Imaginary_1) + U_fab(i, j - 1, k, cIdx::Real_1)*p_fab(i, j - 1, k, cIdx::Real_0) - U_fab(i, j, k, cIdx::Imaginary_0)*p_fab(i + 1, j, k, cIdx::Real_0) + U_fab(i, j, k, cIdx::Imaginary_0)*p_fab(i + 1, j, k, cIdx::Real_1) + U_fab(i, j, k, cIdx::Imaginary_1)*p_fab(i, j + 1, k, cIdx::Imaginary_0) + U_fab(i, j, k, cIdx::Imaginary_1)*p_fab(i, j + 1, k, cIdx::Real_1) - U_fab(i, j, k, cIdx::Real_0)*p_fab(i + 1, j, k, cIdx::Imaginary_0) + U_fab(i, j, k, cIdx::Real_0)*p_fab(i + 1, j, k, cIdx::Imaginary_1) + U_fab(i, j, k, cIdx::Real_1)*p_fab(i, j + 1, k, cIdx::Imaginary_1) - U_fab(i, j, k, cIdx::Real_1)*p_fab(i, j + 1, k, cIdx::Real_0);
        
        Dp_fab(i, j, k, cIdx::Real_0) = Dp_Real_0;
        Dp_fab(i, j, k, cIdx::Imaginary_0) = Dp_Imaginary_0;
        Dp_fab(i, j, k, cIdx::Real_1) = Dp_Real_1;
        Dp_fab(i, j, k, cIdx::Imaginary_1) = Dp_Imaginary_1;
        
}

AMREX_GPU_DEVICE 
inline 
void 
state_set_g3Dp(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& g3Dp_fab,
        amrex::Array4<amrex::Real> const& p_fab,
        amrex::Array4<amrex::Real> const& U_fab,
        amrex::Real m_0) noexcept 
{
        
        amrex::Real g3Dp_Real_0 = (m_0 + 2)*p_fab(i, j, k, cIdx::Real_0) + U_fab(i - 1, j, k, cIdx::Imaginary_0)*p_fab(i - 1, j, k, cIdx::Imaginary_0) + U_fab(i - 1, j, k, cIdx::Imaginary_0)*p_fab(i - 1, j, k, cIdx::Imaginary_1) + U_fab(i - 1, j, k, cIdx::Real_0)*p_fab(i - 1, j, k, cIdx::Real_0) + U_fab(i - 1, j, k, cIdx::Real_0)*p_fab(i - 1, j, k, cIdx::Real_1) + U_fab(i, j - 1, k, cIdx::Imaginary_1)*p_fab(i, j - 1, k, cIdx::Imaginary_0) - U_fab(i, j - 1, k, cIdx::Imaginary_1)*p_fab(i, j - 1, k, cIdx::Real_1) + U_fab(i, j - 1, k, cIdx::Real_1)*p_fab(i, j - 1, k, cIdx::Imaginary_1) + U_fab(i, j - 1, k, cIdx::Real_1)*p_fab(i, j - 1, k, cIdx::Real_0) - U_fab(i, j, k, cIdx::Imaginary_0)*p_fab(i + 1, j, k, cIdx::Imaginary_0) + U_fab(i, j, k, cIdx::Imaginary_0)*p_fab(i + 1, j, k, cIdx::Imaginary_1) - U_fab(i, j, k, cIdx::Imaginary_1)*p_fab(i, j + 1, k, cIdx::Imaginary_0) - U_fab(i, j, k, cIdx::Imaginary_1)*p_fab(i, j + 1, k, cIdx::Real_1) + U_fab(i, j, k, cIdx::Real_0)*p_fab(i + 1, j, k, cIdx::Real_0) - U_fab(i, j, k, cIdx::Real_0)*p_fab(i + 1, j, k, cIdx::Real_1) - U_fab(i, j, k, cIdx::Real_1)*p_fab(i, j + 1, k, cIdx::Imaginary_1) + U_fab(i, j, k, cIdx::Real_1)*p_fab(i, j + 1, k, cIdx::Real_0);
       
        amrex::Real g3Dp_Real_1 = (-m_0 - 2)*p_fab(i, j, k, cIdx::Real_1) - U_fab(i - 1, j, k, cIdx::Imaginary_0)*p_fab(i - 1, j, k, cIdx::Imaginary_0) - U_fab(i - 1, j, k, cIdx::Imaginary_0)*p_fab(i - 1, j, k, cIdx::Imaginary_1) - U_fab(i - 1, j, k, cIdx::Real_0)*p_fab(i - 1, j, k, cIdx::Real_0) - U_fab(i - 1, j, k, cIdx::Real_0)*p_fab(i - 1, j, k, cIdx::Real_1) - U_fab(i, j - 1, k, cIdx::Imaginary_1)*p_fab(i, j - 1, k, cIdx::Imaginary_1) - U_fab(i, j - 1, k, cIdx::Imaginary_1)*p_fab(i, j - 1, k, cIdx::Real_0) + U_fab(i, j - 1, k, cIdx::Real_1)*p_fab(i, j - 1, k, cIdx::Imaginary_0) - U_fab(i, j - 1, k, cIdx::Real_1)*p_fab(i, j - 1, k, cIdx::Real_1) - U_fab(i, j, k, cIdx::Imaginary_0)*p_fab(i + 1, j, k, cIdx::Imaginary_0) + U_fab(i, j, k, cIdx::Imaginary_0)*p_fab(i + 1, j, k, cIdx::Imaginary_1) + U_fab(i, j, k, cIdx::Imaginary_1)*p_fab(i, j + 1, k, cIdx::Imaginary_1) - U_fab(i, j, k, cIdx::Imaginary_1)*p_fab(i, j + 1, k, cIdx::Real_0) + U_fab(i, j, k, cIdx::Real_0)*p_fab(i + 1, j, k, cIdx::Real_0) - U_fab(i, j, k, cIdx::Real_0)*p_fab(i + 1, j, k, cIdx::Real_1) - U_fab(i, j, k, cIdx::Real_1)*p_fab(i, j + 1, k, cIdx::Imaginary_0) - U_fab(i, j, k, cIdx::Real_1)*p_fab(i, j + 1, k, cIdx::Real_1);


        amrex::Real g3Dp_Imaginary_0 = (m_0 + 2)*p_fab(i, j, k, cIdx::Imaginary_0) - U_fab(i - 1, j, k, cIdx::Imaginary_0)*p_fab(i - 1, j, k, cIdx::Real_0) - U_fab(i - 1, j, k, cIdx::Imaginary_0)*p_fab(i - 1, j, k, cIdx::Real_1) + U_fab(i - 1, j, k, cIdx::Real_0)*p_fab(i - 1, j, k, cIdx::Imaginary_0) + U_fab(i - 1, j, k, cIdx::Real_0)*p_fab(i - 1, j, k, cIdx::Imaginary_1) - U_fab(i, j - 1, k, cIdx::Imaginary_1)*p_fab(i, j - 1, k, cIdx::Imaginary_1) - U_fab(i, j - 1, k, cIdx::Imaginary_1)*p_fab(i, j - 1, k, cIdx::Real_0) + U_fab(i, j - 1, k, cIdx::Real_1)*p_fab(i, j - 1, k, cIdx::Imaginary_0) - U_fab(i, j - 1, k, cIdx::Real_1)*p_fab(i, j - 1, k, cIdx::Real_1) + U_fab(i, j, k, cIdx::Imaginary_0)*p_fab(i + 1, j, k, cIdx::Real_0) - U_fab(i, j, k, cIdx::Imaginary_0)*p_fab(i + 1, j, k, cIdx::Real_1) - U_fab(i, j, k, cIdx::Imaginary_1)*p_fab(i, j + 1, k, cIdx::Imaginary_1) + U_fab(i, j, k, cIdx::Imaginary_1)*p_fab(i, j + 1, k, cIdx::Real_0) + U_fab(i, j, k, cIdx::Real_0)*p_fab(i + 1, j, k, cIdx::Imaginary_0) - U_fab(i, j, k, cIdx::Real_0)*p_fab(i + 1, j, k, cIdx::Imaginary_1) + U_fab(i, j, k, cIdx::Real_1)*p_fab(i, j + 1, k, cIdx::Imaginary_0) + U_fab(i, j, k, cIdx::Real_1)*p_fab(i, j + 1, k, cIdx::Real_1);
        
        amrex::Real g3Dp_Imaginary_1 = (-m_0 - 2)*p_fab(i, j, k, cIdx::Imaginary_1) + U_fab(i - 1, j, k, cIdx::Imaginary_0)*p_fab(i - 1, j, k, cIdx::Real_0) + U_fab(i - 1, j, k, cIdx::Imaginary_0)*p_fab(i - 1, j, k, cIdx::Real_1) - U_fab(i - 1, j, k, cIdx::Real_0)*p_fab(i - 1, j, k, cIdx::Imaginary_0) - U_fab(i - 1, j, k, cIdx::Real_0)*p_fab(i - 1, j, k, cIdx::Imaginary_1) - U_fab(i, j - 1, k, cIdx::Imaginary_1)*p_fab(i, j - 1, k, cIdx::Imaginary_0) + U_fab(i, j - 1, k, cIdx::Imaginary_1)*p_fab(i, j - 1, k, cIdx::Real_1) - U_fab(i, j - 1, k, cIdx::Real_1)*p_fab(i, j - 1, k, cIdx::Imaginary_1) - U_fab(i, j - 1, k, cIdx::Real_1)*p_fab(i, j - 1, k, cIdx::Real_0) + U_fab(i, j, k, cIdx::Imaginary_0)*p_fab(i + 1, j, k, cIdx::Real_0) - U_fab(i, j, k, cIdx::Imaginary_0)*p_fab(i + 1, j, k, cIdx::Real_1) - U_fab(i, j, k, cIdx::Imaginary_1)*p_fab(i, j + 1, k, cIdx::Imaginary_0) - U_fab(i, j, k, cIdx::Imaginary_1)*p_fab(i, j + 1, k, cIdx::Real_1) + U_fab(i, j, k, cIdx::Real_0)*p_fab(i + 1, j, k, cIdx::Imaginary_0) - U_fab(i, j, k, cIdx::Real_0)*p_fab(i + 1, j, k, cIdx::Imaginary_1) - U_fab(i, j, k, cIdx::Real_1)*p_fab(i, j + 1, k, cIdx::Imaginary_1) + U_fab(i, j, k, cIdx::Real_1)*p_fab(i, j + 1, k, cIdx::Real_0);
        
        g3Dp_fab(i, j, k, cIdx::Real_0) = g3Dp_Real_0;
        g3Dp_fab(i, j, k, cIdx::Imaginary_0) = g3Dp_Imaginary_0;
        g3Dp_fab(i, j, k, cIdx::Real_1) = g3Dp_Real_1;
        g3Dp_fab(i, j, k, cIdx::Imaginary_1) = g3Dp_Imaginary_1;
        
}

#endif