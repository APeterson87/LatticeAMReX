#ifndef ET_INTEGRATION_FERMION_UPDATE_K_H
#define ET_INTEGRATION_FERMION_UPDATE_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
state_set_g3DinvPhi(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& state_fab,
        amrex::Array4<amrex::Real> const& aux_fab,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
        const amrex::GeometryData& geom,
        amrex::Real m_0,
        amrex::Real r) noexcept 
{
        
        /*amrex::Real g3DinvPhi_0_Real = (m_0 + 2)*aux_fab(i, j, k, auxIdx::DDinvPhi_0_Real) + 0.5*aux_fab(i + 1, j, k, auxIdx::DDinvPhi_0_Imaginary)*state_fab(i, j, k, Idx::U_0_Imaginary) - 0.5*aux_fab(i + 1, j, k, auxIdx::DDinvPhi_0_Real)*state_fab(i, j, k, Idx::U_0_Real) - 0.5*aux_fab(i + 1, j, k, auxIdx::DDinvPhi_1_Imaginary)*state_fab(i, j, k, Idx::U_0_Imaginary) + 0.5*aux_fab(i + 1, j, k, auxIdx::DDinvPhi_1_Real)*state_fab(i, j, k, Idx::U_0_Real) - 0.5*aux_fab(i - 1, j, k, auxIdx::DDinvPhi_0_Imaginary)*state_fab(i - 1, j, k, Idx::U_0_Imaginary) - 0.5*aux_fab(i - 1, j, k, auxIdx::DDinvPhi_0_Real)*state_fab(i - 1, j, k, Idx::U_0_Real) - 0.5*aux_fab(i - 1, j, k, auxIdx::DDinvPhi_1_Imaginary)*state_fab(i - 1, j, k, Idx::U_0_Imaginary) - 0.5*aux_fab(i - 1, j, k, auxIdx::DDinvPhi_1_Real)*state_fab(i - 1, j, k, Idx::U_0_Real) + 0.5*aux_fab(i, j + 1, k, auxIdx::DDinvPhi_0_Imaginary)*state_fab(i, j, k, Idx::U_1_Imaginary) - 0.5*aux_fab(i, j + 1, k, auxIdx::DDinvPhi_0_Real)*state_fab(i, j, k, Idx::U_1_Real) + 0.5*aux_fab(i, j + 1, k, auxIdx::DDinvPhi_1_Imaginary)*state_fab(i, j, k, Idx::U_1_Real) + 0.5*aux_fab(i, j + 1, k, auxIdx::DDinvPhi_1_Real)*state_fab(i, j, k, Idx::U_1_Imaginary) - 0.5*aux_fab(i, j - 1, k, auxIdx::DDinvPhi_0_Imaginary)*state_fab(i, j - 1, k, Idx::U_1_Imaginary) - 0.5*aux_fab(i, j - 1, k, auxIdx::DDinvPhi_0_Real)*state_fab(i, j - 1, k, Idx::U_1_Real) - 0.5*aux_fab(i, j - 1, k, auxIdx::DDinvPhi_1_Imaginary)*state_fab(i, j - 1, k, Idx::U_1_Real) + 0.5*aux_fab(i, j - 1, k, auxIdx::DDinvPhi_1_Real)*state_fab(i, j - 1, k, Idx::U_1_Imaginary);



        amrex::Real g3DinvPhi_0_Imaginary = (m_0 + 2)*aux_fab(i, j, k, auxIdx::DDinvPhi_0_Imaginary) - 0.5*aux_fab(i + 1, j, k, auxIdx::DDinvPhi_0_Imaginary)*state_fab(i, j, k, Idx::U_0_Real) - 0.5*aux_fab(i + 1, j, k, auxIdx::DDinvPhi_0_Real)*state_fab(i, j, k, Idx::U_0_Imaginary) + 0.5*aux_fab(i + 1, j, k, auxIdx::DDinvPhi_1_Imaginary)*state_fab(i, j, k, Idx::U_0_Real) + 0.5*aux_fab(i + 1, j, k, auxIdx::DDinvPhi_1_Real)*state_fab(i, j, k, Idx::U_0_Imaginary) - 0.5*aux_fab(i - 1, j, k, auxIdx::DDinvPhi_0_Imaginary)*state_fab(i - 1, j, k, Idx::U_0_Real) + 0.5*aux_fab(i - 1, j, k, auxIdx::DDinvPhi_0_Real)*state_fab(i - 1, j, k, Idx::U_0_Imaginary) - 0.5*aux_fab(i - 1, j, k, auxIdx::DDinvPhi_1_Imaginary)*state_fab(i - 1, j, k, Idx::U_0_Real) + 0.5*aux_fab(i - 1, j, k, auxIdx::DDinvPhi_1_Real)*state_fab(i - 1, j, k, Idx::U_0_Imaginary) - 0.5*aux_fab(i, j + 1, k, auxIdx::DDinvPhi_0_Imaginary)*state_fab(i, j, k, Idx::U_1_Real) - 0.5*aux_fab(i, j + 1, k, auxIdx::DDinvPhi_0_Real)*state_fab(i, j, k, Idx::U_1_Imaginary) + 0.5*aux_fab(i, j + 1, k, auxIdx::DDinvPhi_1_Imaginary)*state_fab(i, j, k, Idx::U_1_Imaginary) - 0.5*aux_fab(i, j + 1, k, auxIdx::DDinvPhi_1_Real)*state_fab(i, j, k, Idx::U_1_Real) - 0.5*aux_fab(i, j - 1, k, auxIdx::DDinvPhi_0_Imaginary)*state_fab(i, j - 1, k, Idx::U_1_Real) + 0.5*aux_fab(i, j - 1, k, auxIdx::DDinvPhi_0_Real)*state_fab(i, j - 1, k, Idx::U_1_Imaginary) + 0.5*aux_fab(i, j - 1, k, auxIdx::DDinvPhi_1_Imaginary)*state_fab(i, j - 1, k, Idx::U_1_Imaginary) + 0.5*aux_fab(i, j - 1, k, auxIdx::DDinvPhi_1_Real)*state_fab(i, j - 1, k, Idx::U_1_Real);



        amrex::Real g3DinvPhi_1_Real = (-m_0 - 2)*aux_fab(i, j, k, auxIdx::DDinvPhi_1_Real) + 0.5*aux_fab(i + 1, j, k, auxIdx::DDinvPhi_0_Imaginary)*state_fab(i, j, k, Idx::U_0_Imaginary) - 0.5*aux_fab(i + 1, j, k, auxIdx::DDinvPhi_0_Real)*state_fab(i, j, k, Idx::U_0_Real) - 0.5*aux_fab(i + 1, j, k, auxIdx::DDinvPhi_1_Imaginary)*state_fab(i, j, k, Idx::U_0_Imaginary) + 0.5*aux_fab(i + 1, j, k, auxIdx::DDinvPhi_1_Real)*state_fab(i, j, k, Idx::U_0_Real) + 0.5*aux_fab(i - 1, j, k, auxIdx::DDinvPhi_0_Imaginary)*state_fab(i - 1, j, k, Idx::U_0_Imaginary) + 0.5*aux_fab(i - 1, j, k, auxIdx::DDinvPhi_0_Real)*state_fab(i - 1, j, k, Idx::U_0_Real) + 0.5*aux_fab(i - 1, j, k, auxIdx::DDinvPhi_1_Imaginary)*state_fab(i - 1, j, k, Idx::U_0_Imaginary) + 0.5*aux_fab(i - 1, j, k, auxIdx::DDinvPhi_1_Real)*state_fab(i - 1, j, k, Idx::U_0_Real) + 0.5*aux_fab(i, j + 1, k, auxIdx::DDinvPhi_0_Imaginary)*state_fab(i, j, k, Idx::U_1_Real) + 0.5*aux_fab(i, j + 1, k, auxIdx::DDinvPhi_0_Real)*state_fab(i, j, k, Idx::U_1_Imaginary) - 0.5*aux_fab(i, j + 1, k, auxIdx::DDinvPhi_1_Imaginary)*state_fab(i, j, k, Idx::U_1_Imaginary) + 0.5*aux_fab(i, j + 1, k, auxIdx::DDinvPhi_1_Real)*state_fab(i, j, k, Idx::U_1_Real) - 0.5*aux_fab(i, j - 1, k, auxIdx::DDinvPhi_0_Imaginary)*state_fab(i, j - 1, k, Idx::U_1_Real) + 0.5*aux_fab(i, j - 1, k, auxIdx::DDinvPhi_0_Real)*state_fab(i, j - 1, k, Idx::U_1_Imaginary) + 0.5*aux_fab(i, j - 1, k, auxIdx::DDinvPhi_1_Imaginary)*state_fab(i, j - 1, k, Idx::U_1_Imaginary) + 0.5*aux_fab(i, j - 1, k, auxIdx::DDinvPhi_1_Real)*state_fab(i, j - 1, k, Idx::U_1_Real);



        amrex::Real g3DinvPhi_1_Imaginary = (-m_0 - 2)*aux_fab(i, j, k, auxIdx::DDinvPhi_1_Imaginary) - 0.5*aux_fab(i + 1, j, k, auxIdx::DDinvPhi_0_Imaginary)*state_fab(i, j, k, Idx::U_0_Real) - 0.5*aux_fab(i + 1, j, k, auxIdx::DDinvPhi_0_Real)*state_fab(i, j, k, Idx::U_0_Imaginary) + 0.5*aux_fab(i + 1, j, k, auxIdx::DDinvPhi_1_Imaginary)*state_fab(i, j, k, Idx::U_0_Real) + 0.5*aux_fab(i + 1, j, k, auxIdx::DDinvPhi_1_Real)*state_fab(i, j, k, Idx::U_0_Imaginary) + 0.5*aux_fab(i - 1, j, k, auxIdx::DDinvPhi_0_Imaginary)*state_fab(i - 1, j, k, Idx::U_0_Real) - 0.5*aux_fab(i - 1, j, k, auxIdx::DDinvPhi_0_Real)*state_fab(i - 1, j, k, Idx::U_0_Imaginary) + 0.5*aux_fab(i - 1, j, k, auxIdx::DDinvPhi_1_Imaginary)*state_fab(i - 1, j, k, Idx::U_0_Real) - 0.5*aux_fab(i - 1, j, k, auxIdx::DDinvPhi_1_Real)*state_fab(i - 1, j, k, Idx::U_0_Imaginary) + 0.5*aux_fab(i, j + 1, k, auxIdx::DDinvPhi_0_Imaginary)*state_fab(i, j, k, Idx::U_1_Imaginary) - 0.5*aux_fab(i, j + 1, k, auxIdx::DDinvPhi_0_Real)*state_fab(i, j, k, Idx::U_1_Real) + 0.5*aux_fab(i, j + 1, k, auxIdx::DDinvPhi_1_Imaginary)*state_fab(i, j, k, Idx::U_1_Real) + 0.5*aux_fab(i, j + 1, k, auxIdx::DDinvPhi_1_Real)*state_fab(i, j, k, Idx::U_1_Imaginary) + 0.5*aux_fab(i, j - 1, k, auxIdx::DDinvPhi_0_Imaginary)*state_fab(i, j - 1, k, Idx::U_1_Imaginary) + 0.5*aux_fab(i, j - 1, k, auxIdx::DDinvPhi_0_Real)*state_fab(i, j - 1, k, Idx::U_1_Real) + 0.5*aux_fab(i, j - 1, k, auxIdx::DDinvPhi_1_Imaginary)*state_fab(i, j - 1, k, Idx::U_1_Real) - 0.5*aux_fab(i, j - 1, k, auxIdx::DDinvPhi_1_Real)*state_fab(i, j - 1, k, Idx::U_1_Imaginary);

        
        aux_fab(i, j, k, auxIdx::g3DinvPhi_0_Real) = g3DinvPhi_0_Real;
        aux_fab(i, j, k, auxIdx::g3DinvPhi_0_Imaginary) = g3DinvPhi_0_Imaginary;
        aux_fab(i, j, k, auxIdx::g3DinvPhi_1_Real) = g3DinvPhi_1_Real;
        aux_fab(i, j, k, auxIdx::g3DinvPhi_1_Imaginary) = g3DinvPhi_1_Imaginary;*/
        
}
#endif