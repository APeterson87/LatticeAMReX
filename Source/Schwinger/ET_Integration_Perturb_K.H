#ifndef ET_INTEGRATION_PERTURB_K_H
#define ET_INTEGRATION_PERTURB_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

#define TWO_PI 6.283185307179586

AMREX_GPU_DEVICE 
inline 
void 
mom_refresh(int i, int j, int k, 
	    amrex::Array4<amrex::Real> const& state_fab,
	    const amrex::Real sigma_lev) noexcept 
{
  amrex::Real Pi_0 = amrex::RandomNormal(0, sigma_lev);
  amrex::Real Pi_1 = amrex::RandomNormal(0, sigma_lev);
  state_fab(i, j, k, Idx::Pi_0) = Pi_0;
  state_fab(i, j, k, Idx::Pi_1) = Pi_1;
  
  /*
    Box Muller
    amrex::Real r1 = sqrt(-2.0*log(drand48()));
    amrex::Real theta1 = TWO_PI*drand48();
    state_fab(i, j, k, Idx::Pi_0) = r1*cos(theta1);
    amrex::Real r2 = sqrt(-2.0*log(drand48()));
    amrex::Real theta2 = TWO_PI*drand48();
    state_fab(i, j, k, Idx::Pi_1) = r2*cos(theta2); 
  */
}

AMREX_GPU_DEVICE 
inline 
void 
fermion_Perturbation(int i, int j, int k, 
		     amrex::Array4<amrex::Real> const& fermion_fab) noexcept 
{  
  fermion_fab(i, j, k, cIdx::Real_0) = amrex::RandomNormal(0,1)/std::sqrt(2);
  fermion_fab(i, j, k, cIdx::Imag_0) = amrex::RandomNormal(0,1)/std::sqrt(2);
  fermion_fab(i, j, k, cIdx::Real_1) = amrex::RandomNormal(0,1)/std::sqrt(2);
  fermion_fab(i, j, k, cIdx::Imag_1) = amrex::RandomNormal(0,1)/std::sqrt(2); 
}



#endif
