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
        const amrex::GeometryData& geom,
        const amrex::Real Temp_T) noexcept 
{
        const auto domain_xlo = geom.ProbLo();  
        
        amrex::Real s = state_fab(i, j, k, Idx::s);
        amrex::Real sup = state_fab(i, j+1, k, Idx::s);
        amrex::Real sdown = state_fab(i, j-1, k, Idx::s);
        amrex::Real sleft = state_fab(i-1, j, k, Idx::s);
        amrex::Real sright = state_fab(i+1, j, k, Idx::s);
        
        //Get the energy difference from a flip of s(i, j)
        amrex::Real deltaU = 2*s*(sup+sdown+sleft+sright);  
        
        
        //Calculate the change in probability density for flipping s(i, j)
        amrex::Real Pdens = std::exp(-deltaU/Temp_T);  
        
        if (amrex::Random() < 0.1)  //Only about 10% of s(i,j) are considered candidates.  Otherwise strange effects occur.
        {
        
            //Metropolis MCMC step
            if (deltaU <= 0)
                state_fab(i, j, k, Idx::s) *= -1;
            else
            {
                if(amrex::Random() < Pdens)
                    state_fab(i, j, k, Idx::s) *= -1;
            }
        }

}
#endif