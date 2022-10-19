#ifndef ET_INTEGRATION_RHS_K_H
#define ET_INTEGRATION_RHS_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 
#include <complex>

AMREX_GPU_DEVICE 
inline 
void 
state_update_momentum(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& state_fab,
        const amrex::Real time,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
        const amrex::GeometryData& geom,
        amrex::Real m_0,
        amrex::Real r,
        amrex::Real beta,
        bool is_fermi_dynamical,
        amrex::Real dtau) noexcept 
{
         
        /*amrex::Real U_0_00_r = state_fab(i, j, k, Idx::U_0_Real);
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
        amrex::Real U_0_m0_i = state_fab(i-1, j, k, Idx::U_0_Imaginary);*/
        
        
        amrex::Real U_0_00_r = std::cos(state_fab(i, j, k, Idx::theta_0));
        amrex::Real U_0_00_i = std::sin(state_fab(i, j, k, Idx::theta_0));
        
        amrex::Real U_1_p0_r = std::cos(state_fab(i+1, j, k, Idx::theta_1));
        amrex::Real U_1_p0_i = std::sin(state_fab(i+1, j, k, Idx::theta_1));
        
        amrex::Real U_0_0p_r = std::cos(state_fab(i, j+1, k, Idx::theta_0));
        amrex::Real U_0_0p_i = std::sin(state_fab(i, j+1, k, Idx::theta_0));
        
        amrex::Real U_1_00_r = std::cos(state_fab(i, j, k, Idx::theta_1));
        amrex::Real U_1_00_i = std::sin(state_fab(i, j, k, Idx::theta_1));
        
        amrex::Real U_0_0m_r = std::cos(state_fab(i, j-1, k, Idx::theta_0));
        amrex::Real U_0_0m_i = std::sin(state_fab(i, j-1, k, Idx::theta_0));
        
        amrex::Real U_1_pm_r = std::cos(state_fab(i+1, j-1, k, Idx::theta_1));
        amrex::Real U_1_pm_i = std::sin(state_fab(i+1, j-1, k, Idx::theta_1));
        
        amrex::Real U_1_0m_r = std::cos(state_fab(i, j-1, k, Idx::theta_1));
        amrex::Real U_1_0m_i = std::sin(state_fab(i, j-1, k, Idx::theta_1));
        
        amrex::Real U_0_mp_r = std::cos(state_fab(i-1, j+1, k, Idx::theta_0));
        amrex::Real U_0_mp_i = std::sin(state_fab(i-1, j+1, k, Idx::theta_0));
        
        amrex::Real U_1_m0_r = std::cos(state_fab(i-1, j, k, Idx::theta_1));
        amrex::Real U_1_m0_i = std::sin(state_fab(i-1, j, k, Idx::theta_1));
        
        amrex::Real U_0_m0_r = std::cos(state_fab(i-1, j, k, Idx::theta_0));
        amrex::Real U_0_m0_i = std::sin(state_fab(i-1, j, k, Idx::theta_0));
        
        amrex::Real Plaq_0_im = -U_0_00_i*U_0_0p_i*U_1_00_i*U_1_p0_r + U_0_00_i*U_0_0p_i*U_1_00_r*U_1_p0_i + U_0_00_i*U_0_0p_r*U_1_00_i*U_1_p0_i + U_0_00_i*U_0_0p_r*U_1_00_r*U_1_p0_r - U_0_00_r*U_0_0p_i*U_1_00_i*U_1_p0_i - U_0_00_r*U_0_0p_i*U_1_00_r*U_1_p0_r - U_0_00_r*U_0_0p_r*U_1_00_i*U_1_p0_r + U_0_00_r*U_0_0p_r*U_1_00_r*U_1_p0_i;
        
        amrex::Real Plaq_1_im = -U_0_00_i*U_0_0m_i*U_1_0m_i*U_1_pm_r + U_0_00_i*U_0_0m_i*U_1_0m_r*U_1_pm_i - U_0_00_i*U_0_0m_r*U_1_0m_i*U_1_pm_i - U_0_00_i*U_0_0m_r*U_1_0m_r*U_1_pm_r + U_0_00_r*U_0_0m_i*U_1_0m_i*U_1_pm_i + U_0_00_r*U_0_0m_i*U_1_0m_r*U_1_pm_r - U_0_00_r*U_0_0m_r*U_1_0m_i*U_1_pm_r + U_0_00_r*U_0_0m_r*U_1_0m_r*U_1_pm_i;
        
        amrex::Real F_0 = beta*(Plaq_0_im - Plaq_1_im);

        amrex::Real Plaq_2_im = U_0_m0_i*U_0_mp_i*U_1_00_i*U_1_m0_r - U_0_m0_i*U_0_mp_i*U_1_00_r*U_1_m0_i + U_0_m0_i*U_0_mp_r*U_1_00_i*U_1_m0_i + U_0_m0_i*U_0_mp_r*U_1_00_r*U_1_m0_r - U_0_m0_r*U_0_mp_i*U_1_00_i*U_1_m0_i - U_0_m0_r*U_0_mp_i*U_1_00_r*U_1_m0_r + U_0_m0_r*U_0_mp_r*U_1_00_i*U_1_m0_r - U_0_m0_r*U_0_mp_r*U_1_00_r*U_1_m0_i;
        
        amrex::Real F_1 = beta*(Plaq_2_im - Plaq_0_im);
        
        amrex::Real fD_0 = 0;
        amrex::Real fD_1 = 0;
        
        
        state_fab(i, j, k, Idx::P_0) += -(F_0-fD_0)*dtau;
        
        state_fab(i, j, k, Idx::P_1) += -(F_1-fD_1)*dtau;
        

}

AMREX_GPU_DEVICE 
inline 
void 
state_update_gauge(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& state_fab,
        const amrex::Real time,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
        const amrex::GeometryData& geom,
        amrex::Real dtau) noexcept 
{

        amrex::Real theta_0 = state_fab(i, j, k, Idx::theta_0);
        amrex::Real theta_1 = state_fab(i, j, k, Idx::theta_1);
        
        amrex::Real U_0_Real = state_fab(i, j, k, Idx::U_0_Real);
        amrex::Real U_0_Imaginary = state_fab(i, j, k, Idx::U_0_Imaginary);
        
        amrex::Real U_1_Real = state_fab(i, j, k, Idx::U_1_Real);
        amrex::Real U_1_Imaginary = state_fab(i, j, k, Idx::U_1_Imaginary);
        
        amrex::Real P_0 = state_fab(i, j, k, Idx::P_0);
        amrex::Real P_1 = state_fab(i, j, k, Idx::P_1);
        
        state_fab(i, j, k, Idx::theta_0) = theta_0 + dtau*P_0;
        state_fab(i, j, k, Idx::theta_1) = theta_1 + dtau*P_1;
        
        state_fab(i, j, k, Idx::U_0_Real) = std::cos(dtau*P_0)*U_0_Real - std::sin(dtau*P_0)*U_0_Imaginary;
        state_fab(i, j, k, Idx::U_0_Imaginary) = std::sin(dtau*P_0)*U_0_Real + std::cos(dtau*P_0)*U_0_Imaginary;
        
        state_fab(i, j, k, Idx::U_1_Real) = std::cos(dtau*P_1)*U_1_Real - std::sin(dtau*P_1)*U_1_Imaginary;
        state_fab(i, j, k, Idx::U_1_Imaginary) = std::sin(dtau*P_1)*U_1_Real + std::cos(dtau*P_1)*U_1_Imaginary;
        
        
}

AMREX_GPU_DEVICE 
inline 
void 
state_normalize_gauge(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& state_fab) noexcept 
{
        amrex::Real U_0_Real = state_fab(i, j, k, Idx::U_0_Real);
        amrex::Real U_0_Imaginary = state_fab(i, j, k, Idx::U_0_Imaginary);
        
        amrex::Real U_1_Real = state_fab(i, j, k, Idx::U_1_Real);
        amrex::Real U_1_Imaginary = state_fab(i, j, k, Idx::U_1_Imaginary);
        
        amrex::Real P_0 = state_fab(i, j, k, Idx::P_0);
        amrex::Real P_1 = state_fab(i, j, k, Idx::P_1);
        
        
        state_fab(i, j, k, Idx::U_0_Real) = U_0_Real*std::pow(U_0_Real*U_0_Real + U_0_Imaginary*U_0_Imaginary, -0.5);
        state_fab(i, j, k, Idx::U_0_Imaginary) = U_0_Imaginary*std::pow(U_0_Real*U_0_Real + U_0_Imaginary*U_0_Imaginary, -0.5);
        
        state_fab(i, j, k, Idx::U_1_Real) = U_1_Real*std::pow(U_1_Real*U_1_Real + U_1_Imaginary*U_1_Imaginary, -0.5);
        state_fab(i, j, k, Idx::U_1_Imaginary) = U_1_Imaginary*std::pow(U_1_Real*U_1_Real + U_1_Imaginary*U_1_Imaginary, -0.5);
}


#endif