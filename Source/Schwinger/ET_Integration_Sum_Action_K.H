#ifndef ET_INTEGRATION_SUM_ACTION_K_H
#define ET_INTEGRATION_SUM_ACTION_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
amrex::Real
sum_action(const int i, const int j, const int k, 
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
        /*
        amrex::Real Uc1_r = state_fab(i, j, k, Idx::Ucurrent_0_Real);
        amrex::Real Uc1_i = state_fab(i, j, k, Idx::Ucurrent_0_Imaginary);
        
        amrex::Real Uc2_r = state_fab(i+1, j, k, Idx::Ucurrent_1_Real);
        amrex::Real Uc2_i = state_fab(i+1, j, k, Idx::Ucurrent_1_Imaginary);
        
        amrex::Real Uc3_r = state_fab(i, j+1, k, Idx::Ucurrent_0_Real);
        amrex::Real Uc3_i = -state_fab(i, j+1, k, Idx::Ucurrent_0_Imaginary);
        
        amrex::Real Uc4_r = state_fab(i, j, k, Idx::Ucurrent_1_Real);
        amrex::Real Uc4_i = -state_fab(i, j, k, Idx::Ucurrent_1_Imaginary);
        
        amrex::Real plaq_c_re = Uc1_i*Uc2_i*Uc3_i*Uc4_i - Uc1_i*Uc2_i*Uc3_r*Uc4_r - Uc1_i*Uc2_r*Uc3_i*Uc4_r - Uc1_i*Uc2_r*Uc3_r*Uc4_i - Uc1_r*Uc2_i*Uc3_i*Uc4_r - Uc1_r*Uc2_i*Uc3_r*Uc4_i - Uc1_r*Uc2_r*Uc3_i*Uc4_i + Uc1_r*Uc2_r*Uc3_r*Uc4_r;
        
        amrex::Real H_Gauge_cur = (1-plaq_c_re);
        */
        amrex::Real P_0 = state_fab(i, j, k, Idx::P_0);
        amrex::Real P_1 = state_fab(i, j, k, Idx::P_1);
        /*
        amrex::Real P_cur_0 = state_fab(i, j, k, Idx::Pcurrent_0);
        amrex::Real P_cur_1 = state_fab(i, j, k, Idx::Pcurrent_1);
        */
        amrex::Real delta_action = H_Gauge + 0.5*P_0*P_0 + 0.5*P_1*P_1;
            
        return delta_action;
}

AMREX_GPU_DEVICE 
inline 
amrex::Real
sum_action_D(const int i, const int j, const int k, 
        amrex::Array4<amrex::Real> const& state_fab) 
{
        
        amrex::Real PhiDDinvPhi_0_r = state_fab(i, j, k, Idx::Phi_0_Real)*state_fab(i, j, k, Idx::DDinvPhi_0_Real); 
        amrex::Real PhiDDinvPhi_0_i = state_fab(i, j, k, Idx::Phi_0_Imaginary)*state_fab(i, j, k, Idx::DDinvPhi_0_Imaginary); 
        amrex::Real PhiDDinvPhi_1_r = state_fab(i, j, k, Idx::Phi_1_Real)*state_fab(i, j, k, Idx::DDinvPhi_1_Real); 
        amrex::Real PhiDDinvPhi_1_i = state_fab(i, j, k, Idx::Phi_1_Imaginary)*state_fab(i, j, k, Idx::DDinvPhi_1_Imaginary); 
        
        amrex::Real delta_action = (PhiDDinvPhi_0_r+ PhiDDinvPhi_0_i + PhiDDinvPhi_1_r + PhiDDinvPhi_1_i) ;
    
        return delta_action;
}

AMREX_GPU_DEVICE 
inline 
amrex::Real
sum_action_chi(const int i, const int j, const int k, 
        amrex::Array4<amrex::Real> const& state_fab) 
{
        
        amrex::Real chi_0_r = state_fab(i, j, k, Idx::chi_0_Real)*state_fab(i, j, k, Idx::chi_0_Real); 
        amrex::Real chi_0_i = state_fab(i, j, k, Idx::chi_0_Imaginary)*state_fab(i, j, k, Idx::chi_0_Imaginary); 
        amrex::Real chi_1_r = state_fab(i, j, k, Idx::chi_1_Real)*state_fab(i, j, k, Idx::chi_1_Real); 
        amrex::Real chi_1_i = state_fab(i, j, k, Idx::chi_1_Imaginary)*state_fab(i, j, k, Idx::chi_1_Imaginary);  
        
        amrex::Real delta_action = (chi_0_r + chi_0_i + chi_1_r + chi_1_i);
        
    
        return delta_action;
}

AMREX_GPU_DEVICE 
inline 
amrex::Real
dotProd(const int i, const int j, const int k, 
        amrex::Array4<amrex::Real> const& state_fab, int comp_1, int comp_2) 
{
        
        amrex::Real dotProduct = state_fab(i, j, k, comp_1)*state_fab(i, j, k, comp_2); 
            
        return dotProduct;
}

AMREX_GPU_DEVICE 
inline 
void
ResDDPhi(const int i, const int j, const int k, 
        amrex::Array4<amrex::Real> const& state_fab) 
{
        
        state_fab(i, j, k, Idx::DDinvPhi_0_Real) = 0; 
        state_fab(i, j, k, Idx::DDinvPhi_0_Imaginary) = 0; 
        state_fab(i, j, k, Idx::DDinvPhi_1_Real) = 0; 
        state_fab(i, j, k, Idx::DDinvPhi_1_Imaginary) = 0; 
            
}

AMREX_GPU_DEVICE 
inline 
void
CopyComp(const int i, const int j, const int k, 
        amrex::Array4<amrex::Real> const& state_fab, int comp_1, int comp_2) 
{
        
        state_fab(i, j, k, comp_1) = state_fab(i, j, k, comp_2); 
            
}
#endif