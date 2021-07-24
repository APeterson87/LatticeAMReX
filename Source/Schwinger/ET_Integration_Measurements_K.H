#ifndef ET_INTEGRATION_MEASUREMENTS_K_H
#define ET_INTEGRATION_MEASUREMENTS_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

#include<complex>


AMREX_GPU_DEVICE 
inline 
amrex::Real
sum_plaq(const int i, const int j, const int k, 
        amrex::Array4<amrex::Real> const& state_fab) 
{
    
        amrex::Real U1_r = state_fab(i, j, k, Idx::U_0_Real);
        amrex::Real U1_i = state_fab(i, j, k, Idx::U_0_Imaginary);
        
        amrex::Real U2_r = state_fab(i+1, j, k, Idx::U_1_Real);
        amrex::Real U2_i = state_fab(i+1, j, k, Idx::U_1_Imaginary);
        
        amrex::Real U3_r = state_fab(i, j+1, k, Idx::U_0_Real);
        amrex::Real U3_i = -state_fab(i, j+1, k, Idx::U_0_Imaginary);
        
        amrex::Real U4_r = state_fab(i, j, k, Idx::U_1_Real);
        amrex::Real U4_i = -state_fab(i, j, k, Idx::U_1_Imaginary);
        
        amrex::Real plaq_re = U1_i*U2_i*U3_i*U4_i - U1_i*U2_i*U3_r*U4_r - U1_i*U2_r*U3_i*U4_r - U1_i*U2_r*U3_r*U4_i - U1_r*U2_i*U3_i*U4_r - U1_r*U2_i*U3_r*U4_i - U1_r*U2_r*U3_i*U4_i + U1_r*U2_r*U3_r*U4_r;
        
            
        return plaq_re;
}

AMREX_GPU_DEVICE 
inline 
amrex::Real
state_TopCharge(const int i, const int j, const int k, 
        amrex::Array4<amrex::Real> const& smeared_fab) 
{
        //enum {Real_0, Imaginary_0, Real_1, Imaginary_1};
        std::complex<double> Sm_0_00 = std::complex<double>(smeared_fab(i, j, k, cIdx::Real_0), smeared_fab(i, j, k, cIdx::Imaginary_0));
        std::complex<double> Sm_1_p0 = std::complex<double>(smeared_fab(i+1, j, k, cIdx::Real_1), smeared_fab(i+1, j, k, cIdx::Imaginary_1));
        std::complex<double> Sm_0_0p = std::complex<double>(smeared_fab(i, j+1, k, cIdx::Real_0), smeared_fab(i, j+1, k, cIdx::Imaginary_0));
        std::complex<double> Sm_1_00 = std::complex<double>(smeared_fab(i, j, k, cIdx::Real_1), smeared_fab(i, j, k, cIdx::Imaginary_1));
    
        std::complex<double> w = Sm_0_00*Sm_1_p0*std::conj(Sm_0_0p)*std::conj(Sm_1_00);
        amrex::Real topC = std::arg(w);
    
        
        
        return topC;
}
#endif
