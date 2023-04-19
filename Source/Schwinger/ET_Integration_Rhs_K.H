#ifndef ET_INTEGRATION_RHS_K_H
#define ET_INTEGRATION_RHS_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 
#include <complex>

AMREX_GPU_DEVICE 
inline 
void 
state_compute_gauge_force(int i, int j, int k, 
		    amrex::Array4<amrex::Real> const& state_fab,
		    const amrex::GeometryData& geom,
		    amrex::Real beta,
		    amrex::Real dtau) noexcept 
{
  
  const auto domain_xlo = geom.ProbLo(); 
        
  amrex::Real x0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
  amrex::Real x1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1];        
  amrex::Real dx0 = geom.CellSize(0);
  amrex::Real dx1 = geom.CellSize(1);
	
  // U(x,y)_mu
  amrex::Real U_0_00_r = state_fab(i, j, k, Idx::U_0_Real);
  amrex::Real U_0_00_i = state_fab(i, j, k, Idx::U_0_Imag);

  // U(x+1,y)_nu
  amrex::Real U_1_p0_r = state_fab(i+1, j, k, Idx::U_1_Real);
  amrex::Real U_1_p0_i = state_fab(i+1, j, k, Idx::U_1_Imag);

  // U(x,y+1)_mu
  amrex::Real U_0_0p_r = state_fab(i, j+1, k, Idx::U_0_Real);
  amrex::Real U_0_0p_i = state_fab(i, j+1, k, Idx::U_0_Imag);

  // U(x,y)_nu
  amrex::Real U_1_00_r = state_fab(i, j, k, Idx::U_1_Real);
  amrex::Real U_1_00_i = state_fab(i, j, k, Idx::U_1_Imag);


  // U(x,y-1)_mu
  amrex::Real U_0_0m_r = state_fab(i, j-1, k, Idx::U_0_Real);
  amrex::Real U_0_0m_i = state_fab(i, j-1, k, Idx::U_0_Imag);

  // U(x+1,y-1)_nu
  amrex::Real U_1_pm_r = state_fab(i+1, j-1, k, Idx::U_1_Real);
  amrex::Real U_1_pm_i = state_fab(i+1, j-1, k, Idx::U_1_Imag);

  // U(x,y-1)_nu
  amrex::Real U_1_0m_r = state_fab(i, j-1, k, Idx::U_1_Real);
  amrex::Real U_1_0m_i = state_fab(i, j-1, k, Idx::U_1_Imag);

  // U(x-1,y+1)_mu
  amrex::Real U_0_mp_r = state_fab(i-1, j+1, k, Idx::U_0_Real);
  amrex::Real U_0_mp_i = state_fab(i-1, j+1, k, Idx::U_0_Imag);

  // U(x-1,y)_nu
  amrex::Real U_1_m0_r = state_fab(i-1, j, k, Idx::U_1_Real);
  amrex::Real U_1_m0_i = state_fab(i-1, j, k, Idx::U_1_Imag);

  // U(x-1,y)_mu
  amrex::Real U_0_m0_r = state_fab(i-1, j, k, Idx::U_0_Real);
  amrex::Real U_0_m0_i = state_fab(i-1, j, k, Idx::U_0_Imag);
        
  // Compute the gauge force.
  // KEY: The U(x,y)_mu link is denoted at U_mu_xy
  //      where p = x/y + 1, m = x/y - 1
  //      For the imaginary componenets, there will either be 1 real 3 imag
  //      or 3 real 1 imag.
	
  // U(x,y)_mu * U(x+1,y)_nu * U^\dagger(x,y+1)_mu * U^\dagger(x,y)_nu 
  amrex::Real Plaq_0_im = -U_0_00_i*U_0_0p_i*U_1_00_i*U_1_p0_r + U_0_00_i*U_0_0p_i*U_1_00_r*U_1_p0_i +
    U_0_00_i*U_0_0p_r*U_1_00_i*U_1_p0_i + U_0_00_i*U_0_0p_r*U_1_00_r*U_1_p0_r -
    U_0_00_r*U_0_0p_i*U_1_00_i*U_1_p0_i - U_0_00_r*U_0_0p_i*U_1_00_r*U_1_p0_r -
    U_0_00_r*U_0_0p_r*U_1_00_i*U_1_p0_r + U_0_00_r*U_0_0p_r*U_1_00_r*U_1_p0_i;

  // U(x,y-1)_mu * U(x+1,y-1)_nu * U^\dagger(x,y)_mu * U^\dagger(x,y-1)_nu 
  amrex::Real Plaq_1_im = -U_0_00_i*U_0_0m_i*U_1_0m_i*U_1_pm_r + U_0_00_i*U_0_0m_i*U_1_0m_r*U_1_pm_i -
    U_0_00_i*U_0_0m_r*U_1_0m_i*U_1_pm_i - U_0_00_i*U_0_0m_r*U_1_0m_r*U_1_pm_r +
    U_0_00_r*U_0_0m_i*U_1_0m_i*U_1_pm_i + U_0_00_r*U_0_0m_i*U_1_0m_r*U_1_pm_r -
    U_0_00_r*U_0_0m_r*U_1_0m_i*U_1_pm_r + U_0_00_r*U_0_0m_r*U_1_0m_r*U_1_pm_i;
        
  amrex::Real F_0 = beta*(Plaq_0_im - Plaq_1_im);

  // U(x,y)_nu * U^\dagger(x+1,y-1)_mu * U^\dagger(x-1,y)_nu * U^\dagger(x-1,y)_mu 
  amrex::Real Plaq_2_im = U_0_m0_i*U_0_mp_i*U_1_00_i*U_1_m0_r - U_0_m0_i*U_0_mp_i*U_1_00_r*U_1_m0_i +
    U_0_m0_i*U_0_mp_r*U_1_00_i*U_1_m0_i + U_0_m0_i*U_0_mp_r*U_1_00_r*U_1_m0_r -
    U_0_m0_r*U_0_mp_i*U_1_00_i*U_1_m0_i - U_0_m0_r*U_0_mp_i*U_1_00_r*U_1_m0_r +
    U_0_m0_r*U_0_mp_r*U_1_00_i*U_1_m0_r - U_0_m0_r*U_0_mp_r*U_1_00_r*U_1_m0_i;
        
  amrex::Real F_1 = beta*(Plaq_2_im - Plaq_0_im);
        
  state_fab(i, j, k, Idx::F_U_0) = F_0;
  state_fab(i, j, k, Idx::F_U_1) = F_1;
}

AMREX_GPU_DEVICE 
inline 
void 
state_update_momentum(int i, int j, int k, 
		      amrex::Array4<amrex::Real> const& state_fab,
		      const amrex::GeometryData& geom,
		      amrex::Real beta,
		      amrex::Real dtau) noexcept 
{
         
  const auto domain_xlo = geom.ProbLo(); 
        
  amrex::Real x0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
  amrex::Real x1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1];        
  amrex::Real dx0 = geom.CellSize(0);
  amrex::Real dx1 = geom.CellSize(1);
  
  state_fab(i, j, k, Idx::Pi_0) += (state_fab(i, j, k, Idx::F_Phi_0) - state_fab(i, j, k, Idx::F_U_0))*dtau/(dx0*dx1);
  state_fab(i, j, k, Idx::Pi_1) += (state_fab(i, j, k, Idx::F_Phi_1) - state_fab(i, j, k, Idx::F_U_1))*dtau/(dx0*dx1);
}


// Update the gauge field U with the momentum kick P
//U_{k} = exp(i dtau P_{k-1/2}) * U_{k-1}
AMREX_GPU_DEVICE 
inline 
void 
state_update_gauge(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& state_fab,
        const amrex::GeometryData& geom,
        amrex::Real dtau) noexcept 
{
        const auto domain_xlo = geom.ProbLo();
        
        amrex::Real x0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1];
        
        amrex::Real dx0 = geom.CellSize(0);
        amrex::Real dx1 = geom.CellSize(1);
	
	// U(x,y)_mu
        amrex::Real U_0_Re = state_fab(i, j, k, Idx::U_0_Real);
        amrex::Real U_0_Im = state_fab(i, j, k, Idx::U_0_Imag);

	// U(x,y)_nu
        amrex::Real U_1_Re = state_fab(i, j, k, Idx::U_1_Real);
        amrex::Real U_1_Im = state_fab(i, j, k, Idx::U_1_Imag);

	// P(x,y)
	amrex::Real Pmu = dtau*state_fab(i, j, k, Idx::Pi_0)*dx0;
	amrex::Real Pnu = dtau*state_fab(i, j, k, Idx::Pi_1)*dx1; 

	// Update U(x,y)_mu
        state_fab(i, j, k, Idx::U_0_Real)      = std::cos(Pmu)*U_0_Re - std::sin(Pmu)*U_0_Im;
        state_fab(i, j, k, Idx::U_0_Imag) = std::sin(Pmu)*U_0_Re + std::cos(Pmu)*U_0_Im;

	// Update U(x,y)_nu
        state_fab(i, j, k, Idx::U_1_Real)      = std::cos(Pnu)*U_1_Re - std::sin(Pnu)*U_1_Im;
        state_fab(i, j, k, Idx::U_1_Imag) = std::sin(Pnu)*U_1_Re + std::cos(Pnu)*U_1_Im;
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
        amrex::Real U_0_Imag = state_fab(i, j, k, Idx::U_0_Imag);
        
        amrex::Real U_1_Real = state_fab(i, j, k, Idx::U_1_Real);
        amrex::Real U_1_Imag = state_fab(i, j, k, Idx::U_1_Imag);
        
        amrex::Real P_0 = state_fab(i, j, k, Idx::Pi_0);
        amrex::Real P_1 = state_fab(i, j, k, Idx::Pi_1);
        
        
        /*state_fab(i, j, k, Idx::U_0_Real) = U_0_Real*std::pow(U_0_Real*U_0_Real + U_0_Imag*U_0_Imag, -0.5);
        state_fab(i, j, k, Idx::U_0_Imag) = U_0_Imag*std::pow(U_0_Real*U_0_Real + U_0_Imag*U_0_Imag, -0.5);
        
        state_fab(i, j, k, Idx::U_1_Real) = U_1_Real*std::pow(U_1_Real*U_1_Real + U_1_Imag*U_1_Imag, -0.5);
        state_fab(i, j, k, Idx::U_1_Imag) = U_1_Imag*std::pow(U_1_Real*U_1_Real + U_1_Imag*U_1_Imag, -0.5);*/
        
        amrex::Real U1_r = state_fab(i, j, k, Idx::U_0_Real);
        amrex::Real U1_i = state_fab(i, j, k, Idx::U_0_Imag);
        
        amrex::Real U2_r = state_fab(i+1, j, k, Idx::U_1_Real);
        amrex::Real U2_i = state_fab(i+1, j, k, Idx::U_1_Imag);
        
        amrex::Real U3_r = state_fab(i, j+1, k, Idx::U_0_Real);
        amrex::Real U3_i = -state_fab(i, j+1, k, Idx::U_0_Imag);
        
        amrex::Real U4_r = state_fab(i, j, k, Idx::U_1_Real);
        amrex::Real U4_i = -state_fab(i, j, k, Idx::U_1_Imag);
        
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
        amrex::Real plaq_re = (U1_i*U2_i*U3_i*U4_i - U1_i*U2_i*U3_r*U4_r -
			       U1_i*U2_r*U3_i*U4_r - U1_i*U2_r*U3_r*U4_i -
			       U1_r*U2_i*U3_i*U4_r - U1_r*U2_i*U3_r*U4_i -
			       U1_r*U2_r*U3_i*U4_i + U1_r*U2_r*U3_r*U4_r);
        
        state_fab(i, j, k, Idx::ActionDensity) = 0.5*state_fab(i, j, k, Idx::Pi_0)*state_fab(i, j, k, Idx::Pi_0) + 0.5*state_fab(i, j, k, Idx::Pi_1)*state_fab(i, j, k, Idx::Pi_1) + beta*(1-plaq_re)/(dx0*dx1);
}

#endif
