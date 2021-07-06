#ifndef ET_INTEGRATION_RHS_K_H
#define ET_INTEGRATION_RHS_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
state_rhs(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& rhs_fab, 
        amrex::Array4<amrex::Real const> const& state_fab,
        const amrex::Real time,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx, 
        const amrex::GeometryData& geom,
        amrex::Real factor) noexcept 
{


        amrex::Real U_0_00_r = state_fab(i, j, k, Idx::U_0_Real);
        amrex::Real U_0_00_i = state_fab(i, j, k, Idx::U_0_Imaginary);
        
        amrex::Real U_1_p0_r = state_fab(i+1, j, k, Idx::U_1_Real);
        amrex::Real U_1_p0_i = state_fab(i+1, j, k, Idx::U_1_Imaginary);
        
        amrex::Real U_0_0p_r = state_fab(i, j+1, k, Idx::U_0_Real);
        amrex::Real U_0_0p_i = state_fab(i, j+1, k, Idx::U_0_Imaginary);
        
        amrex::Real U_1_00_r = state_fab(i, j, k, Idx::U_1_Real);
        amrex::Real U_1_00_i = state_fab(i, j, k, Idx::U_1_Imaginary);
        
        amrex::Real U_0_0m_r = state_fab(i, j-1, k, Idx::U_0_Real);
        amrex::Real U_0_0m_i = state_fab(i, j-1, k, Idx::U_0_Imaginary);
        
        amrex::Real U_1_pm_r = state_fab(i+1, j-1, k, Idx::U_1_Real);
        amrex::Real U_1_pm_i = state_fab(i+1, j-1, k, Idx::U_1_Imaginary);
        
        amrex::Real U_1_0m_r = state_fab(i, j-1, k, Idx::U_1_Real);
        amrex::Real U_1_0m_i = state_fab(i, j-1, k, Idx::U_1_Imaginary);
        
        amrex::Real U_0_mp_r = state_fab(i-1, j+1, k, Idx::U_0_Real);
        amrex::Real U_0_mp_i = state_fab(i-1, j+1, k, Idx::U_0_Imaginary);
        
        amrex::Real U_1_m0_r = state_fab(i-1, j, k, Idx::U_1_Real);
        amrex::Real U_1_m0_i = state_fab(i-1, j, k, Idx::U_1_Imaginary);
        
        amrex::Real U_0_m0_r = state_fab(i-1, j, k, Idx::U_0_Real);
        amrex::Real U_0_m0_i = state_fab(i-1, j, k, Idx::U_0_Imaginary);
        
        
        amrex::Real Plaq_0_im = -U_0_00_i*U_0_0p_i*U_1_00_i*U_1_p0_r + U_0_00_i*U_0_0p_i*U_1_00_r*U_1_p0_i + U_0_00_i*U_0_0p_r*U_1_00_i*U_1_p0_i + U_0_00_i*U_0_0p_r*U_1_00_r*U_1_p0_r - U_0_00_r*U_0_0p_i*U_1_00_i*U_1_p0_i - U_0_00_r*U_0_0p_i*U_1_00_r*U_1_p0_r - U_0_00_r*U_0_0p_r*U_1_00_i*U_1_p0_r + U_0_00_r*U_0_0p_r*U_1_00_r*U_1_p0_i;
        
        amrex::Real Plaq_1_im = -U_0_00_i*U_0_0m_i*U_1_0m_i*U_1_pm_r + U_0_00_i*U_0_0m_i*U_1_0m_r*U_1_pm_i - U_0_00_i*U_0_0m_r*U_1_0m_i*U_1_pm_i - U_0_00_i*U_0_0m_r*U_1_0m_r*U_1_pm_r + U_0_00_r*U_0_0m_i*U_1_0m_i*U_1_pm_i + U_0_00_r*U_0_0m_i*U_1_0m_r*U_1_pm_r - U_0_00_r*U_0_0m_r*U_1_0m_i*U_1_pm_r + U_0_00_r*U_0_0m_r*U_1_0m_r*U_1_pm_i;
        
        amrex::Real F_0 = (Plaq_0_im - Plaq_1_im);

        amrex::Real Plaq_2_im = U_0_m0_i*U_0_mp_i*U_1_00_i*U_1_m0_r - U_0_m0_i*U_0_mp_i*U_1_00_r*U_1_m0_i + U_0_m0_i*U_0_mp_r*U_1_00_i*U_1_m0_i + U_0_m0_i*U_0_mp_r*U_1_00_r*U_1_m0_r - U_0_m0_r*U_0_mp_i*U_1_00_i*U_1_m0_i - U_0_m0_r*U_0_mp_i*U_1_00_r*U_1_m0_r + U_0_m0_r*U_0_mp_r*U_1_00_i*U_1_m0_r - U_0_m0_r*U_0_mp_r*U_1_00_r*U_1_m0_i;
        
        amrex::Real F_1 = (Plaq_2_im - Plaq_0_im);
        
        rhs_fab(i, j, k, Idx::P_0) = -factor*F_0;
        
        rhs_fab(i, j, k, Idx::P_1) = -factor*F_1;
        
        
        /*
        amrex::Real theta_A_0_1 = state_fab(i+1, j, k, Idx::Q_1)-state_fab(i, j+1, k, Idx::Q_0)-state_fab(i, j, k, Idx::Q_1);
        
        amrex::Real theta_A_0_2 = -state_fab(i+1, j-1, k, Idx::Q_1)-state_fab(i, j-1, k, Idx::Q_0)+state_fab(i, j-1, k, Idx::Q_1);
        
        amrex::Real theta_A_1_1 = state_fab(i, j+1, k, Idx::Q_0)-state_fab(i+1, j, k, Idx::Q_1)-state_fab(i, j, k, Idx::Q_0);
        
        amrex::Real theta_A_1_2 = -state_fab(i-1, j+1, k, Idx::Q_0)-state_fab(i-1, j, k, Idx::Q_1)+state_fab(i-1, j, k, Idx::Q_0);
        
        amrex::Real Q_0 = state_fab(i, j, k, Idx::Q_0);
        
        amrex::Real Q_1 = state_fab(i, j, k, Idx::Q_1);
        
        amrex::Real F_0 = std::sin(Q_0 + theta_A_0_1) + std::sin(Q_0 + theta_A_0_2);
        
        amrex::Real F_1 = std::sin(Q_1 + theta_A_1_1) + std::sin(Q_1 + theta_A_1_2);

        rhs_fab(i, j, k, Idx::P_0) = -factor*F_0;
        
        rhs_fab(i, j, k, Idx::P_1) = -factor*F_1;
        */
        
}
#endif