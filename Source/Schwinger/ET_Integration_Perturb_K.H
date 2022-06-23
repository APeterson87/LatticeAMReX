#ifndef ET_INTEGRATION_PERTURB_K_H
#define ET_INTEGRATION_MCMC_PERTURB_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
state_Perturbation(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& state_fab,
        const amrex::Real time) noexcept 
{
        
        amrex::Real P_0 = amrex::RandomNormal(0,1);
        amrex::Real P_1 = amrex::RandomNormal(0,1);
        
        state_fab(i, j, k, Idx::P_0) = P_0;
        state_fab(i, j, k, Idx::P_1) = P_1;
        
}

AMREX_GPU_DEVICE 
inline 
void 
fermi_Perturbation(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& fermi_fab,
        const amrex::Real time) noexcept 
{
        
        fermi_fab(i, j, k, cIdx::Real_0) = amrex::RandomNormal(0,1)/std::sqrt(2);
        fermi_fab(i, j, k, cIdx::Imaginary_0) = amrex::RandomNormal(0,1)/std::sqrt(2);
        fermi_fab(i, j, k, cIdx::Real_1) = amrex::RandomNormal(0,1)/std::sqrt(2);
        fermi_fab(i, j, k, cIdx::Imaginary_1) = amrex::RandomNormal(0,1)/std::sqrt(2);
        
}

#endif