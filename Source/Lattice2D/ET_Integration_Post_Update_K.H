#ifndef ET_INTEGRATION_POST_UPDATE_K_H
#define ET_INTEGRATION_POST_UPDATE_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
state_post_update(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& state_fab,
        const amrex::Real time,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
        const amrex::GeometryData& geom) noexcept 
{
        const auto domain_xlo = geom.ProbLo();
        
        amrex::Real x0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1]; 
        amrex::Real x2 = (k + 0.5)*geom.CellSize(2) + domain_xlo[2];
        
        amrex::Real theta_A_0_1 = std::acos(state_fab(i+1, j, k, Idx::U_1_Re))-std::acos(state_fab(i, j+1, k, Idx::U_0_Re))-std::acos(state_fab(i, j, k, Idx::U_1_Re));
        
        amrex::Real theta_A_0_2 = -std::acos(state_fab(i+1, j-1, k, Idx::U_1_Re))-std::acos(state_fab(i, j-1, k, Idx::U_0_Re))+std::acos(state_fab(i, j-1, k, Idx::U_1_Re));
        
        amrex::Real theta_A_1_1 = std::acos(state_fab(i, j+1, k, Idx::U_0_Re))-std::acos(state_fab(i+1, j, k, Idx::U_1_Re))-std::acos(state_fab(i, j, k, Idx::U_0_Re));
        
        amrex::Real theta_A_1_2 = -std::acos(state_fab(i-1, j+1, k, Idx::U_0_Re))-std::acos(state_fab(i-1, j, k, Idx::U_1_Re))+std::acos(state_fab(i-1, j, k, Idx::U_0_Re));
        
        amrex::Real theta_0 = std::acos(state_fab(i, j, k, Idx::U_0_Re));
        amrex::Real theta_1 = std::acos(state_fab(i, j, k, Idx::U_1_Re));
        
        amrex::Real theta_0_prime = theta_0 + 0.1*(amrex::Random()-0.5);
        amrex::Real theta_1_prime = theta_1 + 0.1*(amrex::Random()-0.5);
        
        amrex::Real deltaS_0 = -(std::cos(theta_0_prime + theta_A_0_1)+std::cos(theta_0_prime + theta_A_0_2)-std::cos(theta_0 + theta_A_0_1)-std::cos(theta_0 + theta_A_0_2));
        
        amrex::Real deltaS_1 = -(std::cos(theta_1_prime + theta_A_1_1)+std::cos(theta_1_prime + theta_A_1_2)-std::cos(theta_1 + theta_A_1_1)-std::cos(theta_1 + theta_A_1_2));
        
        amrex::Real r0 = amrex::Random();
        
        amrex::Real r1 = amrex::Random();
        
        amrex::Real T = 0.01;
        
        if (amrex::Random() < 0.1)
        {
            state_fab(i, j, k, Idx::U_0_Re) = (r0 < std::min(1.0,std::exp(-deltaS_0/T)))*std::cos(theta_0_prime) + (r0 >= std::min(1.0,std::exp(-deltaS_0/T)))*std::cos(theta_0);
        
            state_fab(i, j, k, Idx::U_0_Im) = (r0 < std::min(1.0,std::exp(-deltaS_0/T)))*std::sin(theta_0_prime) + (r0 >= std::min(1.0,std::exp(-deltaS_0/T)))*std::sin(theta_0);
        }
        
        if (amrex::Random() < 0.1)
        {
            state_fab(i, j, k, Idx::U_1_Re) = (r1 < std::min(1.0,std::exp(-deltaS_1/T)))*std::cos(theta_1_prime) + (r1 >= std::min(1.0,std::exp(-deltaS_1/T)))*std::cos(theta_1);
        
            state_fab(i, j, k, Idx::U_1_Im) = (r1 < std::min(1.0,std::exp(-deltaS_1/T)))*std::sin(theta_1_prime) + (r1 >= std::min(1.0,std::exp(-deltaS_1/T)))*std::sin(theta_1);
        }
        
        
        
        
        

}
#endif