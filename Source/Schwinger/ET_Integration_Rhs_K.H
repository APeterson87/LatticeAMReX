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
        amrex::Array4<amrex::Real> const& aux_fab,
        const amrex::Real time,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
        const amrex::GeometryData& geom,
        amrex::Real m_0,
        amrex::Real r,
        amrex::Real beta,
        bool is_fermi_dynamical,
        amrex::Real dtau) noexcept 
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
        
        amrex::Real F_0 = beta*(Plaq_0_im - Plaq_1_im);

        amrex::Real Plaq_2_im = U_0_m0_i*U_0_mp_i*U_1_00_i*U_1_m0_r - U_0_m0_i*U_0_mp_i*U_1_00_r*U_1_m0_i + U_0_m0_i*U_0_mp_r*U_1_00_i*U_1_m0_i + U_0_m0_i*U_0_mp_r*U_1_00_r*U_1_m0_r - U_0_m0_r*U_0_mp_i*U_1_00_i*U_1_m0_i - U_0_m0_r*U_0_mp_i*U_1_00_r*U_1_m0_r + U_0_m0_r*U_0_mp_r*U_1_00_i*U_1_m0_r - U_0_m0_r*U_0_mp_r*U_1_00_r*U_1_m0_i;
        
        amrex::Real F_1 = beta*(Plaq_2_im - Plaq_0_im);
        
        //amrex::Real m_0 = 0.1;
        //amrex::Real r  = 1.0;
        amrex::Real fD_0 = 0;
        amrex::Real fD_1 = 0;
        
        if(is_fermi_dynamical)
        {
            fD_0 = (aux_fab(i + 1, j, k, auxIdx::DDinvPhi_0_Imaginary)*aux_fab(i, j, k, auxIdx::g3DinvPhi_0_Imaginary) + aux_fab(i + 1, j, k, auxIdx::DDinvPhi_0_Imaginary)*aux_fab(i, j, k, auxIdx::g3DinvPhi_1_Imaginary) + aux_fab(i + 1, j, k, auxIdx::DDinvPhi_0_Real)*aux_fab(i, j, k, auxIdx::g3DinvPhi_0_Real) + aux_fab(i + 1, j, k, auxIdx::DDinvPhi_0_Real)*aux_fab(i, j, k, auxIdx::g3DinvPhi_1_Real))*state_fab(i, j, k, Idx::U_0_Imaginary) + (aux_fab(i + 1, j, k, auxIdx::DDinvPhi_0_Imaginary)*aux_fab(i, j, k, auxIdx::g3DinvPhi_0_Real) + aux_fab(i + 1, j, k, auxIdx::DDinvPhi_0_Imaginary)*aux_fab(i, j, k, auxIdx::g3DinvPhi_1_Real) - aux_fab(i + 1, j, k, auxIdx::DDinvPhi_0_Real)*aux_fab(i, j, k, auxIdx::g3DinvPhi_0_Imaginary) - aux_fab(i + 1, j, k, auxIdx::DDinvPhi_0_Real)*aux_fab(i, j, k, auxIdx::g3DinvPhi_1_Imaginary))*state_fab(i, j, k, Idx::U_0_Real) - (aux_fab(i + 1, j, k, auxIdx::DDinvPhi_1_Imaginary)*aux_fab(i, j, k, auxIdx::g3DinvPhi_0_Imaginary) + aux_fab(i + 1, j, k, auxIdx::DDinvPhi_1_Imaginary)*aux_fab(i, j, k, auxIdx::g3DinvPhi_1_Imaginary) + aux_fab(i + 1, j, k, auxIdx::DDinvPhi_1_Real)*aux_fab(i, j, k, auxIdx::g3DinvPhi_0_Real) + aux_fab(i + 1, j, k, auxIdx::DDinvPhi_1_Real)*aux_fab(i, j, k, auxIdx::g3DinvPhi_1_Real))*state_fab(i, j, k, Idx::U_0_Imaginary) - (aux_fab(i + 1, j, k, auxIdx::DDinvPhi_1_Imaginary)*aux_fab(i, j, k, auxIdx::g3DinvPhi_0_Real) + aux_fab(i + 1, j, k, auxIdx::DDinvPhi_1_Imaginary)*aux_fab(i, j, k, auxIdx::g3DinvPhi_1_Real) - aux_fab(i + 1, j, k, auxIdx::DDinvPhi_1_Real)*aux_fab(i, j, k, auxIdx::g3DinvPhi_0_Imaginary) - aux_fab(i + 1, j, k, auxIdx::DDinvPhi_1_Real)*aux_fab(i, j, k, auxIdx::g3DinvPhi_1_Imaginary))*state_fab(i, j, k, Idx::U_0_Real) + (aux_fab(i + 1, j, k, auxIdx::g3DinvPhi_0_Imaginary)*aux_fab(i, j, k, auxIdx::DDinvPhi_0_Imaginary) + aux_fab(i + 1, j, k, auxIdx::g3DinvPhi_0_Real)*aux_fab(i, j, k, auxIdx::DDinvPhi_0_Real) - aux_fab(i + 1, j, k, auxIdx::g3DinvPhi_1_Imaginary)*aux_fab(i, j, k, auxIdx::DDinvPhi_0_Imaginary) - aux_fab(i + 1, j, k, auxIdx::g3DinvPhi_1_Real)*aux_fab(i, j, k, auxIdx::DDinvPhi_0_Real))*state_fab(i, j, k, Idx::U_0_Imaginary) + (aux_fab(i + 1, j, k, auxIdx::g3DinvPhi_0_Imaginary)*aux_fab(i, j, k, auxIdx::DDinvPhi_0_Real) - aux_fab(i + 1, j, k, auxIdx::g3DinvPhi_0_Real)*aux_fab(i, j, k, auxIdx::DDinvPhi_0_Imaginary) - aux_fab(i + 1, j, k, auxIdx::g3DinvPhi_1_Imaginary)*aux_fab(i, j, k, auxIdx::DDinvPhi_0_Real) + aux_fab(i + 1, j, k, auxIdx::g3DinvPhi_1_Real)*aux_fab(i, j, k, auxIdx::DDinvPhi_0_Imaginary))*state_fab(i, j, k, Idx::U_0_Real) + (aux_fab(i + 1, j, k, auxIdx::g3DinvPhi_0_Imaginary)*aux_fab(i, j, k, auxIdx::DDinvPhi_1_Imaginary) + aux_fab(i + 1, j, k, auxIdx::g3DinvPhi_0_Real)*aux_fab(i, j, k, auxIdx::DDinvPhi_1_Real) - aux_fab(i + 1, j, k, auxIdx::g3DinvPhi_1_Imaginary)*aux_fab(i, j, k, auxIdx::DDinvPhi_1_Imaginary) - aux_fab(i + 1, j, k, auxIdx::g3DinvPhi_1_Real)*aux_fab(i, j, k, auxIdx::DDinvPhi_1_Real))*state_fab(i, j, k, Idx::U_0_Imaginary) + (aux_fab(i + 1, j, k, auxIdx::g3DinvPhi_0_Imaginary)*aux_fab(i, j, k, auxIdx::DDinvPhi_1_Real) - aux_fab(i + 1, j, k, auxIdx::g3DinvPhi_0_Real)*aux_fab(i, j, k, auxIdx::DDinvPhi_1_Imaginary) - aux_fab(i + 1, j, k, auxIdx::g3DinvPhi_1_Imaginary)*aux_fab(i, j, k, auxIdx::DDinvPhi_1_Real) + aux_fab(i + 1, j, k, auxIdx::g3DinvPhi_1_Real)*aux_fab(i, j, k, auxIdx::DDinvPhi_1_Imaginary))*state_fab(i, j, k, Idx::U_0_Real);



            fD_1 = (aux_fab(i, j + 1, k, auxIdx::DDinvPhi_0_Imaginary)*aux_fab(i, j, k, auxIdx::g3DinvPhi_0_Imaginary) - aux_fab(i, j + 1, k, auxIdx::DDinvPhi_0_Imaginary)*aux_fab(i, j, k, auxIdx::g3DinvPhi_1_Real) + aux_fab(i, j + 1, k, auxIdx::DDinvPhi_0_Real)*aux_fab(i, j, k, auxIdx::g3DinvPhi_0_Real) + aux_fab(i, j + 1, k, auxIdx::DDinvPhi_0_Real)*aux_fab(i, j, k, auxIdx::g3DinvPhi_1_Imaginary))*state_fab(i, j, k, Idx::U_1_Imaginary) + (aux_fab(i, j + 1, k, auxIdx::DDinvPhi_0_Imaginary)*aux_fab(i, j, k, auxIdx::g3DinvPhi_0_Real) + aux_fab(i, j + 1, k, auxIdx::DDinvPhi_0_Imaginary)*aux_fab(i, j, k, auxIdx::g3DinvPhi_1_Imaginary) - aux_fab(i, j + 1, k, auxIdx::DDinvPhi_0_Real)*aux_fab(i, j, k, auxIdx::g3DinvPhi_0_Imaginary) + aux_fab(i, j + 1, k, auxIdx::DDinvPhi_0_Real)*aux_fab(i, j, k, auxIdx::g3DinvPhi_1_Real))*state_fab(i, j, k, Idx::U_1_Real) + (aux_fab(i, j + 1, k, auxIdx::DDinvPhi_1_Imaginary)*aux_fab(i, j, k, auxIdx::g3DinvPhi_0_Imaginary) - aux_fab(i, j + 1, k, auxIdx::DDinvPhi_1_Imaginary)*aux_fab(i, j, k, auxIdx::g3DinvPhi_1_Real) + aux_fab(i, j + 1, k, auxIdx::DDinvPhi_1_Real)*aux_fab(i, j, k, auxIdx::g3DinvPhi_0_Real) + aux_fab(i, j + 1, k, auxIdx::DDinvPhi_1_Real)*aux_fab(i, j, k, auxIdx::g3DinvPhi_1_Imaginary))*state_fab(i, j, k, Idx::U_1_Real) - (aux_fab(i, j + 1, k, auxIdx::DDinvPhi_1_Imaginary)*aux_fab(i, j, k, auxIdx::g3DinvPhi_0_Real) + aux_fab(i, j + 1, k, auxIdx::DDinvPhi_1_Imaginary)*aux_fab(i, j, k, auxIdx::g3DinvPhi_1_Imaginary) - aux_fab(i, j + 1, k, auxIdx::DDinvPhi_1_Real)*aux_fab(i, j, k, auxIdx::g3DinvPhi_0_Imaginary) + aux_fab(i, j + 1, k, auxIdx::DDinvPhi_1_Real)*aux_fab(i, j, k, auxIdx::g3DinvPhi_1_Real))*state_fab(i, j, k, Idx::U_1_Imaginary) + (aux_fab(i, j + 1, k, auxIdx::g3DinvPhi_0_Imaginary)*aux_fab(i, j, k, auxIdx::DDinvPhi_0_Imaginary) + aux_fab(i, j + 1, k, auxIdx::g3DinvPhi_0_Real)*aux_fab(i, j, k, auxIdx::DDinvPhi_0_Real) - aux_fab(i, j + 1, k, auxIdx::g3DinvPhi_1_Imaginary)*aux_fab(i, j, k, auxIdx::DDinvPhi_0_Real) + aux_fab(i, j + 1, k, auxIdx::g3DinvPhi_1_Real)*aux_fab(i, j, k, auxIdx::DDinvPhi_0_Imaginary))*state_fab(i, j, k, Idx::U_1_Imaginary) + (aux_fab(i, j + 1, k, auxIdx::g3DinvPhi_0_Imaginary)*aux_fab(i, j, k, auxIdx::DDinvPhi_0_Real) - aux_fab(i, j + 1, k, auxIdx::g3DinvPhi_0_Real)*aux_fab(i, j, k, auxIdx::DDinvPhi_0_Imaginary) + aux_fab(i, j + 1, k, auxIdx::g3DinvPhi_1_Imaginary)*aux_fab(i, j, k, auxIdx::DDinvPhi_0_Imaginary) + aux_fab(i, j + 1, k, auxIdx::g3DinvPhi_1_Real)*aux_fab(i, j, k, auxIdx::DDinvPhi_0_Real))*state_fab(i, j, k, Idx::U_1_Real) + (aux_fab(i, j + 1, k, auxIdx::g3DinvPhi_0_Imaginary)*aux_fab(i, j, k, auxIdx::DDinvPhi_1_Imaginary) + aux_fab(i, j + 1, k, auxIdx::g3DinvPhi_0_Real)*aux_fab(i, j, k, auxIdx::DDinvPhi_1_Real) - aux_fab(i, j + 1, k, auxIdx::g3DinvPhi_1_Imaginary)*aux_fab(i, j, k, auxIdx::DDinvPhi_1_Real) + aux_fab(i, j + 1, k, auxIdx::g3DinvPhi_1_Real)*aux_fab(i, j, k, auxIdx::DDinvPhi_1_Imaginary))*state_fab(i, j, k, Idx::U_1_Real) - (aux_fab(i, j + 1, k, auxIdx::g3DinvPhi_0_Imaginary)*aux_fab(i, j, k, auxIdx::DDinvPhi_1_Real) - aux_fab(i, j + 1, k, auxIdx::g3DinvPhi_0_Real)*aux_fab(i, j, k, auxIdx::DDinvPhi_1_Imaginary) + aux_fab(i, j + 1, k, auxIdx::g3DinvPhi_1_Imaginary)*aux_fab(i, j, k, auxIdx::DDinvPhi_1_Imaginary) + aux_fab(i, j + 1, k, auxIdx::g3DinvPhi_1_Real)*aux_fab(i, j, k, auxIdx::DDinvPhi_1_Real))*state_fab(i, j, k, Idx::U_1_Imaginary);
            
        }
        
        
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
        amrex::Real U_0_Real = state_fab(i, j, k, Idx::U_0_Real);
        amrex::Real U_0_Imaginary = state_fab(i, j, k, Idx::U_0_Imaginary);
        
        amrex::Real U_1_Real = state_fab(i, j, k, Idx::U_1_Real);
        amrex::Real U_1_Imaginary = state_fab(i, j, k, Idx::U_1_Imaginary);
        
        amrex::Real P_0 = state_fab(i, j, k, Idx::P_0);
        amrex::Real P_1 = state_fab(i, j, k, Idx::P_1);
        
        
        state_fab(i, j, k, Idx::U_0_Real) = std::cos(dtau*P_0)*U_0_Real - std::sin(dtau*P_0)*U_0_Imaginary;
        state_fab(i, j, k, Idx::U_0_Imaginary) = std::sin(dtau*P_0)*U_0_Real + std::cos(dtau*P_0)*U_0_Imaginary;
        
        state_fab(i, j, k, Idx::U_1_Real) = std::cos(dtau*P_1)*U_1_Real - std::sin(dtau*P_1)*U_1_Imaginary;
        state_fab(i, j, k, Idx::U_1_Imaginary) = std::sin(dtau*P_1)*U_1_Real + std::cos(dtau*P_1)*U_1_Imaginary;
}

AMREX_GPU_DEVICE 
inline 
void 
state_update_topChargeDensity(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& state_fab,
        amrex::Array4<amrex::Real> const& smeared_fab,
        const amrex::Real time,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
        const amrex::GeometryData& geom) noexcept 
{
        std::complex<double> Sm_0_00 = std::complex<double>(smeared_fab(i, j, k, cIdx::Real_0), smeared_fab(i, j, k, cIdx::Imaginary_0));
        std::complex<double> Sm_1_p0 = std::complex<double>(smeared_fab(i+1, j, k, cIdx::Real_1), smeared_fab(i+1, j, k, cIdx::Imaginary_1));
        std::complex<double> Sm_0_0p = std::complex<double>(smeared_fab(i, j+1, k, cIdx::Real_0), smeared_fab(i, j+1, k, cIdx::Imaginary_0));
        std::complex<double> Sm_1_00 = std::complex<double>(smeared_fab(i, j, k, cIdx::Real_1), smeared_fab(i, j, k, cIdx::Imaginary_1));
    
        std::complex<double> w = Sm_0_00*Sm_1_p0*std::conj(Sm_0_0p)*std::conj(Sm_1_00);
        
        amrex::Real topC = std::arg(w);
        
        state_fab(i, j, k, Idx::TopChargeDensity) = topC/(2*M_PI);
        
        amrex::Real U1_r = state_fab(i, j, k, Idx::U_0_Real);
        amrex::Real U1_i = state_fab(i, j, k, Idx::U_0_Imaginary);
        
        amrex::Real U2_r = state_fab(i+1, j, k, Idx::U_1_Real);
        amrex::Real U2_i = state_fab(i+1, j, k, Idx::U_1_Imaginary);
        
        amrex::Real U3_r = state_fab(i, j+1, k, Idx::U_0_Real);
        amrex::Real U3_i = -state_fab(i, j+1, k, Idx::U_0_Imaginary);
        
        amrex::Real U4_r = state_fab(i, j, k, Idx::U_1_Real);
        amrex::Real U4_i = -state_fab(i, j, k, Idx::U_1_Imaginary);
        
        amrex::Real plaq_re = U1_i*U2_i*U3_i*U4_i - U1_i*U2_i*U3_r*U4_r - U1_i*U2_r*U3_i*U4_r - U1_i*U2_r*U3_r*U4_i - U1_r*U2_i*U3_i*U4_r - U1_r*U2_i*U3_r*U4_i - U1_r*U2_r*U3_i*U4_i + U1_r*U2_r*U3_r*U4_r;
        
        state_fab(i, j, k, Idx::ActionDensity) = 1-plaq_re;
}

#endif