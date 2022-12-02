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
         
        const auto domain_xlo = geom.ProbLo(); 
        
        amrex::Real x0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1];
        
        amrex::Real dx0 = geom.CellSize(0);
        amrex::Real dx1 = geom.CellSize(1);
        
         
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
        
        /*
        amrex::Real U_0_00_r = std::cos(state_fab(i, j, k, Idx::Theta_0));
        amrex::Real U_0_00_i = std::sin(state_fab(i, j, k, Idx::Theta_0));
        
        amrex::Real U_1_p0_r = std::cos(state_fab(i+1, j, k, Idx::Theta_1));
        amrex::Real U_1_p0_i = std::sin(state_fab(i+1, j, k, Idx::Theta_1));
        
        amrex::Real U_0_0p_r = std::cos(state_fab(i, j+1, k, Idx::Theta_0));
        amrex::Real U_0_0p_i = std::sin(state_fab(i, j+1, k, Idx::Theta_0));
        
        amrex::Real U_1_00_r = std::cos(state_fab(i, j, k, Idx::Theta_1));
        amrex::Real U_1_00_i = std::sin(state_fab(i, j, k, Idx::Theta_1));
        
        amrex::Real U_0_0m_r = std::cos(state_fab(i, j-1, k, Idx::Theta_0));
        amrex::Real U_0_0m_i = std::sin(state_fab(i, j-1, k, Idx::Theta_0));
        
        amrex::Real U_1_pm_r = std::cos(state_fab(i+1, j-1, k, Idx::Theta_1));
        amrex::Real U_1_pm_i = std::sin(state_fab(i+1, j-1, k, Idx::Theta_1));
        
        amrex::Real U_1_0m_r = std::cos(state_fab(i, j-1, k, Idx::Theta_1));
        amrex::Real U_1_0m_i = std::sin(state_fab(i, j-1, k, Idx::Theta_1));
        
        amrex::Real U_0_mp_r = std::cos(state_fab(i-1, j+1, k, Idx::Theta_0));
        amrex::Real U_0_mp_i = std::sin(state_fab(i-1, j+1, k, Idx::Theta_0));
        
        amrex::Real U_1_m0_r = std::cos(state_fab(i-1, j, k, Idx::Theta_1));
        amrex::Real U_1_m0_i = std::sin(state_fab(i-1, j, k, Idx::Theta_1));
        
        amrex::Real U_0_m0_r = std::cos(state_fab(i-1, j, k, Idx::Theta_0));
        amrex::Real U_0_m0_i = std::sin(state_fab(i-1, j, k, Idx::Theta_0));
        */
        
        amrex::Real Plaq_0_im = -U_0_00_i*U_0_0p_i*U_1_00_i*U_1_p0_r + U_0_00_i*U_0_0p_i*U_1_00_r*U_1_p0_i + U_0_00_i*U_0_0p_r*U_1_00_i*U_1_p0_i + U_0_00_i*U_0_0p_r*U_1_00_r*U_1_p0_r - U_0_00_r*U_0_0p_i*U_1_00_i*U_1_p0_i - U_0_00_r*U_0_0p_i*U_1_00_r*U_1_p0_r - U_0_00_r*U_0_0p_r*U_1_00_i*U_1_p0_r + U_0_00_r*U_0_0p_r*U_1_00_r*U_1_p0_i;
        
        amrex::Real Plaq_1_im = -U_0_00_i*U_0_0m_i*U_1_0m_i*U_1_pm_r + U_0_00_i*U_0_0m_i*U_1_0m_r*U_1_pm_i - U_0_00_i*U_0_0m_r*U_1_0m_i*U_1_pm_i - U_0_00_i*U_0_0m_r*U_1_0m_r*U_1_pm_r + U_0_00_r*U_0_0m_i*U_1_0m_i*U_1_pm_i + U_0_00_r*U_0_0m_i*U_1_0m_r*U_1_pm_r - U_0_00_r*U_0_0m_r*U_1_0m_i*U_1_pm_r + U_0_00_r*U_0_0m_r*U_1_0m_r*U_1_pm_i;
        
        amrex::Real F_0 = beta*(Plaq_0_im - Plaq_1_im);

        amrex::Real Plaq_2_im = U_0_m0_i*U_0_mp_i*U_1_00_i*U_1_m0_r - U_0_m0_i*U_0_mp_i*U_1_00_r*U_1_m0_i + U_0_m0_i*U_0_mp_r*U_1_00_i*U_1_m0_i + U_0_m0_i*U_0_mp_r*U_1_00_r*U_1_m0_r - U_0_m0_r*U_0_mp_i*U_1_00_i*U_1_m0_i - U_0_m0_r*U_0_mp_i*U_1_00_r*U_1_m0_r + U_0_m0_r*U_0_mp_r*U_1_00_i*U_1_m0_r - U_0_m0_r*U_0_mp_r*U_1_00_r*U_1_m0_i;
        
        amrex::Real F_1 = beta*(Plaq_2_im - Plaq_0_im);
        
        amrex::Real fD_0 = 0;
        amrex::Real fD_1 = 0;
        
        state_fab(i, j, k, Idx::Pi_0) += -(F_0-fD_0)*dtau/(dx0*dx1); //dx0; //(dx1*dx0);
        
        state_fab(i, j, k, Idx::Pi_1) += -(F_1-fD_1)*dtau/(dx0*dx1); //dx1; //(dx0*dx1);
        /*amrex::Real ddA_00 = (state_fab(i+1, j, k, Idx::A_0)-2*state_fab(i, j, k, Idx::A_0)+state_fab(i-1, j, k, Idx::A_0))/(dx0*dx0);
        amrex::Real ddA_01 = (state_fab(i, j+1, k, Idx::A_0)-2*state_fab(i, j, k, Idx::A_0)+state_fab(i, j-1, k, Idx::A_0))/(dx1*dx1);
        
        amrex::Real ddA_10 = (state_fab(i+1, j, k, Idx::A_1)-2*state_fab(i, j, k, Idx::A_1)+state_fab(i-1, j, k, Idx::A_1))/(dx0*dx0);
        amrex::Real ddA_11 = (state_fab(i, j+1, k, Idx::A_1)-2*state_fab(i, j, k, Idx::A_1)+state_fab(i, j-1, k, Idx::A_1))/(dx1*dx1);
        
        state_fab(i, j, k, Idx::Pi_0
        ) += beta*dtau*(ddA_00 + ddA_01);
        state_fab(i, j, k, Idx::Pi_1) += beta*dtau*(ddA_10 + ddA_11);*/
        
        
        
        //state_fab(i, j, k, Idx::P_0) += -(F_0-fD_0)*dtau;
        
        //state_fab(i, j, k, Idx::P_1) += -(F_1-fD_1)*dtau;
        

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
        const auto domain_xlo = geom.ProbLo();
        
        amrex::Real x0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1];
        
        amrex::Real dx0 = geom.CellSize(0);
        amrex::Real dx1 = geom.CellSize(1);

        //amrex::Real Theta_0 = state_fab(i, j, k, Idx::Theta_0);
        //amrex::Real Theta_1 = state_fab(i, j, k, Idx::Theta_1);
        
        amrex::Real Pi_0 = state_fab(i, j, k, Idx::Pi_0);
        amrex::Real Pi_1 = state_fab(i, j, k, Idx::Pi_1);
        /*
        state_fab(i, j, k, Idx::Theta_0) = Theta_0 + dtau*Pi_0;
        state_fab(i, j, k, Idx::Theta_1) = Theta_1 + dtau*Pi_1;
        */

        //amrex::Real theta_0 = state_fab(i, j, k, Idx::theta_0);
        //amrex::Real theta_1 = state_fab(i, j, k, Idx::theta_1);
        
        /*amrex::Real U_0_Real = std::cos(Theta_0);
        amrex::Real U_0_Imaginary = std::sin(Theta_0);
        
        amrex::Real U_1_Real = std::cos(Theta_1);
        amrex::Real U_1_Imaginary = std::sin(Theta_1);*/
        
        amrex::Real U_0_Real = state_fab(i, j, k, Idx::U_0_Real);
        amrex::Real U_0_Imaginary = state_fab(i, j, k, Idx::U_0_Imaginary);
        
        amrex::Real U_1_Real = state_fab(i, j, k, Idx::U_1_Real);
        amrex::Real U_1_Imaginary = state_fab(i, j, k, Idx::U_1_Imaginary);
        
        //amrex::Real P_0 = state_fab(i, j, k, Idx::P_0);
        //amrex::Real P_1 = state_fab(i, j, k, Idx::P_1);
        
        
        
        //state_fab(i, j, k, Idx::theta_0) = theta_0 + dtau*P_0;
        //state_fab(i, j, k, Idx::theta_1) = theta_1 + dtau*P_1;
        
        state_fab(i, j, k, Idx::U_0_Real) = std::cos(dtau*Pi_0*dx0)*U_0_Real - std::sin(dtau*Pi_0*dx0)*U_0_Imaginary;
        state_fab(i, j, k, Idx::U_0_Imaginary) = std::sin(dtau*Pi_0*dx0)*U_0_Real + std::cos(dtau*Pi_0*dx0)*U_0_Imaginary;
        
        state_fab(i, j, k, Idx::U_1_Real) = std::cos(dtau*Pi_1*dx1)*U_1_Real - std::sin(dtau*Pi_1*dx1)*U_1_Imaginary;
        state_fab(i, j, k, Idx::U_1_Imaginary) = std::sin(dtau*Pi_1*dx1)*U_1_Real + std::cos(dtau*Pi_1*dx1)*U_1_Imaginary;
        
        
}

AMREX_GPU_DEVICE 
inline 
void 
state_normalize_gauge(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& state_fab,
        const amrex::GeometryData& geom,
        amrex::Real beta) noexcept 
{

        const auto domain_xlo = geom.ProbLo();
        
        amrex::Real x0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1];
        
        amrex::Real dx0 = geom.CellSize(0);
        amrex::Real dx1 = geom.CellSize(1);
        
        amrex::Real U_0_Real = state_fab(i, j, k, Idx::U_0_Real);
        amrex::Real U_0_Imaginary = state_fab(i, j, k, Idx::U_0_Imaginary);
        
        amrex::Real U_1_Real = state_fab(i, j, k, Idx::U_1_Real);
        amrex::Real U_1_Imaginary = state_fab(i, j, k, Idx::U_1_Imaginary);
        
        amrex::Real P_0 = state_fab(i, j, k, Idx::Pi_0);
        amrex::Real P_1 = state_fab(i, j, k, Idx::Pi_1);
        
        
        /*state_fab(i, j, k, Idx::U_0_Real) = U_0_Real*std::pow(U_0_Real*U_0_Real + U_0_Imaginary*U_0_Imaginary, -0.5);
        state_fab(i, j, k, Idx::U_0_Imaginary) = U_0_Imaginary*std::pow(U_0_Real*U_0_Real + U_0_Imaginary*U_0_Imaginary, -0.5);
        
        state_fab(i, j, k, Idx::U_1_Real) = U_1_Real*std::pow(U_1_Real*U_1_Real + U_1_Imaginary*U_1_Imaginary, -0.5);
        state_fab(i, j, k, Idx::U_1_Imaginary) = U_1_Imaginary*std::pow(U_1_Real*U_1_Real + U_1_Imaginary*U_1_Imaginary, -0.5);*/
        
        amrex::Real U1_r = state_fab(i, j, k, Idx::U_0_Real);
        amrex::Real U1_i = state_fab(i, j, k, Idx::U_0_Imaginary);
        
        amrex::Real U2_r = state_fab(i+1, j, k, Idx::U_1_Real);
        amrex::Real U2_i = state_fab(i+1, j, k, Idx::U_1_Imaginary);
        
        amrex::Real U3_r = state_fab(i, j+1, k, Idx::U_0_Real);
        amrex::Real U3_i = -state_fab(i, j+1, k, Idx::U_0_Imaginary);
        
        amrex::Real U4_r = state_fab(i, j, k, Idx::U_1_Real);
        amrex::Real U4_i = -state_fab(i, j, k, Idx::U_1_Imaginary);
        
        /*
        amrex::Real U1_r = std::cos(state_fab(i, j, k, Idx::Theta_0));
        amrex::Real U1_i = std::sin(state_fab(i, j, k, Idx::Theta_0));
        
        amrex::Real U2_r = std::cos(state_fab(i+1, j, k, Idx::Theta_1));
        amrex::Real U2_i = std::sin(state_fab(i+1, j, k, Idx::Theta_1));
        
        amrex::Real U3_r = std::cos(state_fab(i, j+1, k, Idx::Theta_0));
        amrex::Real U3_i = -std::sin(state_fab(i, j+1, k, Idx::Theta_0));
        
        amrex::Real U4_r = std::cos(state_fab(i, j, k, Idx::Theta_1));
        amrex::Real U4_i = -std::sin(state_fab(i, j, k, Idx::Theta_1));
        */
        amrex::Real plaq_re = U1_i*U2_i*U3_i*U4_i - U1_i*U2_i*U3_r*U4_r - U1_i*U2_r*U3_i*U4_r - U1_i*U2_r*U3_r*U4_i - U1_r*U2_i*U3_i*U4_r - U1_r*U2_i*U3_r*U4_i - U1_r*U2_r*U3_i*U4_i + U1_r*U2_r*U3_r*U4_r;
        
        state_fab(i, j, k, Idx::ActionDensity) = 0.5*state_fab(i, j, k, Idx::Pi_0)*state_fab(i, j, k, Idx::Pi_0) + 0.5*state_fab(i, j, k, Idx::Pi_1)*state_fab(i, j, k, Idx::Pi_1) + beta*(1-plaq_re)/(dx0*dx1);
}


#endif