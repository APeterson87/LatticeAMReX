#ifndef ET_INTEGRATION_SUM_ACTION_K_H
#define ET_INTEGRATION_SUM_ACTION_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
amrex::Real
sum_action_gauge(const int i, const int j, const int k, 
        amrex::Array4<amrex::Real> const& state_fab,
        amrex::Real beta) 
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
        
        amrex::Real H_Gauge = beta*(1-plaq_re);
        
        amrex::Real P_0 = state_fab(i, j, k, Idx::P_0);
        amrex::Real P_1 = state_fab(i, j, k, Idx::P_1);

        amrex::Real delta_action = H_Gauge;
            
        return delta_action;
}

AMREX_GPU_DEVICE 
inline 
amrex::Real
sum_action_mom(const int i, const int j, const int k, 
        amrex::Array4<amrex::Real> const& state_fab) 
{
        
        amrex::Real P_0 = state_fab(i, j, k, Idx::P_0);
        amrex::Real P_1 = state_fab(i, j, k, Idx::P_1);

        amrex::Real delta_action = 0.5*P_0*P_0 + 0.5*P_1*P_1;
            
        return delta_action;
}

AMREX_GPU_DEVICE 
inline 
amrex::Real
sum_action_D(const int i, const int j, const int k, 
        amrex::Array4<amrex::Real> const& state_fab,
        amrex::Array4<amrex::Real> const& Ainvp_fab) 
{
        
        amrex::Real PhiDDinvPhi_0_r = state_fab(i, j, k, Idx::Phi_0_Real)*Ainvp_fab(i, j, k, cIdx::Real_0); 
        amrex::Real PhiDDinvPhi_0_i = state_fab(i, j, k, Idx::Phi_0_Imaginary)*Ainvp_fab(i, j, k, cIdx::Imaginary_0); 
        amrex::Real PhiDDinvPhi_1_r = state_fab(i, j, k, Idx::Phi_1_Real)*Ainvp_fab(i, j, k, cIdx::Real_1); 
        amrex::Real PhiDDinvPhi_1_i = state_fab(i, j, k, Idx::Phi_1_Imaginary)*Ainvp_fab(i, j, k, cIdx::Imaginary_1); 
        
        amrex::Real delta_action = (PhiDDinvPhi_0_r+ PhiDDinvPhi_0_i + PhiDDinvPhi_1_r + PhiDDinvPhi_1_i) ;
    
        return delta_action;
}

#endif