#ifndef ET_INTEGRATION_UTILS_K_H
#define ET_INTEGRATION_UTILS_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

using namespace amrex;

struct Observables {
            int num_accepted;
            Real TimeSumAvePlaq;
            Real TimeSumDeltaH;
    
            std::vector<int> TopologicalCharge;
            std::vector<Real> deltaH;
    
            std::vector<std::vector<double>> sigma_meas;
            std::vector<double> sigma_ave;
            std::vector<double> sigma_stdev;
    
            std::vector<std::vector<double>> pion_correlation;
            std::vector<double> pion_corr_ave;
            std::vector<double> pion_corr_stdev;
    
            Observables() : num_accepted(0), 
                            TimeSumAvePlaq(0), 
                            TimeSumDeltaH(0) { }
        };

struct Parameters {
    
    //Physical Parameters
            Real beta = 3.0;
            Real mass = 0.1;
            Real r = 1.0;
            Real Temp = 1.0;
            
            std::vector<Real> beta_lev;
            std::vector<Real> sigma_lev;
            std::vector<Real> Temp_lev;
            std::vector<Real> hmc_tau_lev;
    
    //Problem Parameters
            int Nx = 8;
            int Ny = 8;
    
            int hmc_substeps = 15;
            int therm_steps = 100;
            int starting_meas_step = 200;
    
            Real tau = 1.0;
    
            bool use_dynamical_fermions = true;
    
            Real BiCG_Thresh = 1e-6;
            int BiCG_Max_Iters = 100;
            bool stabilized = false;
    
            int APE_smear_iter = 1.0;
            Real APE_smear_alpha = 0.5;
            int measWilsonLoops_Interval = 10;
            int checkrevtraj_Interval = 5;
    

};


// DMH: Use this as a template for Wilson Flow.
//      This is a cheap, gauge only physicial measurement

AMREX_GPU_DEVICE 
inline 
void 
state_smear_gauge(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& smeared_fab,
	const amrex::Real time,
        amrex::Real alpha) noexcept 
{

        amrex::Real U_0_00_r = smeared_fab(i, j, k, cIdx::Real_0);
        amrex::Real U_0_00_i = smeared_fab(i, j, k, cIdx::Imaginary_0);
        std::complex<double> U_0_00 = std::complex<double>(U_0_00_r, U_0_00_i);
        
        amrex::Real U_1_p0_r = smeared_fab(i+1, j, k, cIdx::Real_1);
        amrex::Real U_1_p0_i = smeared_fab(i+1, j, k, cIdx::Imaginary_1);
        std::complex<double> U_1_p0 = std::complex<double>(U_1_p0_r, U_1_p0_i);
        
        amrex::Real U_0_0p_r = smeared_fab(i, j+1, k, cIdx::Real_0);
        amrex::Real U_0_0p_i = smeared_fab(i, j+1, k, cIdx::Imaginary_0);
        std::complex<double> U_0_0p = std::complex<double>(U_0_0p_r, U_0_0p_i);
        
        amrex::Real U_1_00_r = smeared_fab(i, j, k, cIdx::Real_1);
        amrex::Real U_1_00_i = smeared_fab(i, j, k, cIdx::Imaginary_1);
        std::complex<double> U_1_00 = std::complex<double>(U_1_00_r, U_1_00_i);
        
        amrex::Real U_0_0m_r = smeared_fab(i, j-1, k, cIdx::Real_0);
        amrex::Real U_0_0m_i = smeared_fab(i, j-1, k, cIdx::Imaginary_0);
        std::complex<double> U_0_0m = std::complex<double>(U_0_0m_r, U_0_0m_i);
        
        amrex::Real U_1_pm_r = smeared_fab(i+1, j-1, k, cIdx::Real_1);
        amrex::Real U_1_pm_i = smeared_fab(i+1, j-1, k, cIdx::Imaginary_1);
        std::complex<double> U_1_pm = std::complex<double>(U_1_pm_r, U_1_pm_i);
        
        amrex::Real U_1_0m_r = smeared_fab(i, j-1, k, cIdx::Real_1);
        amrex::Real U_1_0m_i = smeared_fab(i, j-1, k, cIdx::Imaginary_1);
        std::complex<double> U_1_0m = std::complex<double>(U_1_0m_r, U_1_0m_i);
        
        amrex::Real U_0_mp_r = smeared_fab(i-1, j+1, k, cIdx::Real_0);
        amrex::Real U_0_mp_i = smeared_fab(i-1, j+1, k, cIdx::Imaginary_0);
        std::complex<double> U_0_mp = std::complex<double>(U_0_mp_r, U_0_mp_i);
        
        amrex::Real U_1_m0_r = smeared_fab(i-1, j, k, cIdx::Real_1);
        amrex::Real U_1_m0_i = smeared_fab(i-1, j, k, cIdx::Imaginary_1);
        std::complex<double> U_1_m0 = std::complex<double>(U_1_m0_r, U_1_m0_i);
        
        amrex::Real U_0_m0_r = smeared_fab(i-1, j, k, cIdx::Real_0);
        amrex::Real U_0_m0_i = smeared_fab(i-1, j, k, cIdx::Imaginary_0);
        std::complex<double> U_0_m0 = std::complex<double>(U_0_m0_r, U_0_m0_i);
        
        std::complex<double> staple_0 = U_1_00*U_0_0p*std::conj(U_1_p0) + std::conj(U_1_0m)*U_0_0m*U_1_pm;
        std::complex<double> staple_1 = U_0_00*U_1_p0*std::conj(U_0_0p) + std::conj(U_0_m0)*U_1_m0*U_0_mp;
        
        std::complex<double> SmU_0_00 = (1.0-alpha)*U_0_00 + alpha/2.0*staple_0;
        std::complex<double> SmU_1_00 = (1.0-alpha)*U_1_00 + alpha/2.0*staple_1;
                
	// Add the projection here
	amrex::Real SmU_tmp_0_00_Real = std::real(SmU_0_00);
        amrex::Real SmU_tmp_0_00_Imaginary = std::imag(SmU_0_00);
        std::complex<double> SmU_tmp_0_00 = std::complex<double>(SmU_tmp_0_00_Real, SmU_tmp_0_00_Imaginary);
        
        amrex::Real SmU_tmp_1_00_Real = std::real(SmU_1_00);
        amrex::Real SmU_tmp_1_00_Imaginary = std::imag(SmU_1_00);
        std::complex<double> SmU_tmp_1_00 = std::complex<double>(SmU_tmp_1_00_Real, SmU_tmp_1_00_Imaginary);
        
        smeared_fab(i, j, k, cIdx::Real_0) = std::real(std::polar(1.0,std::arg(SmU_tmp_0_00)));
        smeared_fab(i, j, k, cIdx::Imaginary_0) = std::imag(std::polar(1.0,std::arg(SmU_tmp_0_00)));
        smeared_fab(i, j, k, cIdx::Real_1) = std::real(std::polar(1.0,std::arg(SmU_tmp_1_00)));
        smeared_fab(i, j, k, cIdx::Imaginary_1) = std::imag(std::polar(1.0,std::arg(SmU_tmp_1_00)));

	
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void adam_avgdown_nodes (Box const& bx, Array4<T> const& crse,
                          Array4<T const> const& fine,
                          int ccomp, int fcomp, int ncomp,
                          IntVect const& ratio) noexcept
{
    const auto clo = lbound(bx);
    const auto chi = ubound(bx);
    const int facx = ratio[0];
    const int facy = ratio[1];
    Real alpha = 0.1;

        for     (int j = clo.y; j <= chi.y; ++j) {
            int jj = j*facy;
            AMREX_PRAGMA_SIMD
            for (int i = clo.x; i <= chi.x; ++i) {
                int ii = i*facx;
                
                int itake = amrex::Random_int(2);
                int jtake = amrex::Random_int(2);
                    /*
                    crse(i,j,0,ccomp) = fine(ii+itake,jj+jtake,0,fcomp);
                    crse(i,j,0,1+ccomp) = fine(ii+itake,jj+jtake,0,1+fcomp);
                    */
                    
                    Real Uf_0_R_00 = fine(ii, jj, 0, fcomp);
                    Real Uf_0_I_00 = fine(ii, jj, 0, 1+fcomp);
                    
                    Real Thetaf_0_00 = std::atan2(Uf_0_I_00,Uf_0_R_00);
                    
                    Real Uf_0_R_0p = fine(ii, jj+1, 0, fcomp);
                    Real Uf_0_I_0p = fine(ii, jj+1, 0, 1+fcomp);
                    
                    Real Thetaf_0_0p = std::atan2(Uf_0_I_0p,Uf_0_R_0p);
                    
                    Real Uf_1_R_00 = fine(ii, jj, 0, 2+fcomp);
                    Real Uf_1_I_00 = fine(ii, jj, 0, 3+fcomp);
                    
                    Real Thetaf_1_00 = std::atan2(Uf_1_I_00,Uf_1_R_00);
                    
                    Real Uf_1_R_p0 = fine(ii+1, jj, 0, 2+fcomp);
                    Real Uf_1_I_p0 = fine(ii+1, jj, 0, 3+fcomp);
                    
                    Real Thetaf_1_p0 = std::atan2(Uf_1_I_p0,Uf_1_R_p0);
                    
                    Real Theta_0 = 0.5*(Thetaf_0_00 + Thetaf_0_0p);
                    Real Theta_1 = 0.5*(Thetaf_1_00 + Thetaf_1_p0);
                    
                    crse(i,j,0,ccomp) = std::cos(Theta_0);
                    crse(i,j,0,1+ccomp) = std::sin(Theta_0);
                    crse(i,j,0,2+ccomp) = std::cos(Theta_1);
                    crse(i,j,0,3+ccomp) = std::sin(Theta_1);
                    
                    //crse(i,j,0,ccomp) = 0.5*(fine(ii,jj,0,fcomp)+fine(ii,jj+1,0,fcomp));
                    //crse(i,j,0,1+ccomp) = 0.5*(fine(ii,jj,0,1+fcomp)+fine(ii+1,jj,0,1+fcomp));
                    
                    /*
                    Real staple0 = 0;
                    staple0 += alpha/2*(fine(ii,jj,0,1+fcomp)+fine(ii,jj+1,0,fcomp)-fine(ii+1,jj,0,1+fcomp));
                    staple0 += alpha/2*(-fine(ii,jj-1,0,1+fcomp) + fine(ii,jj-1,0,fcomp) + fine(ii+1,jj-1,0,1+fcomp));
                    
                    Real staple1 = 0;
                    staple1 += alpha/2*(fine(ii,jj,0,fcomp) + fine(ii+1,jj,0,fcomp+1) - fine(ii,jj+1,0,fcomp)); 
                    staple1 += alpha/2*(-fine(ii-1,jj,0,fcomp) + fine(ii-1,jj,0,fcomp+1) + fine(ii-1,jj+1,0,fcomp));
                    
                    crse(i,j,0,ccomp) = (1-alpha)*fine(ii,jj,0,fcomp) + staple0;
                    crse(i,j,0,1+ccomp) = (1-alpha)*fine(ii,jj,0,1+fcomp) + staple1;
                    
                    
                    /*crse(i,j,0,2+ccomp) = 0.5*(fine(ii,jj,0,2+fcomp) + fine(ii+1,jj,0,2+fcomp));
                    crse(i,j,0,3+ccomp) = 0.5*(fine(ii,jj,0,3+fcomp) + fine(ii,jj+1,0,3+fcomp));
                    crse(i,j,0,8+ccomp) = 0*(fine(ii,jj,0,8+fcomp) + fine(ii,jj+1,0,8+fcomp)+fine(ii+1,jj,0,8+fcomp) + fine(ii+1,jj+1,0,8+fcomp));
                    */
                    
                    //crse(i,j,0,ccomp) = 0.25*(fine(ii,jj,0,0+fcomp) + fine(ii,jj+1,0,0+fcomp)+fine(ii+1,jj,0,0+fcomp) + fine(ii+1,jj+1,0,0+fcomp));
                    //crse(i,j,0,1+ccomp) = 0.25*(fine(ii,jj,0,1+fcomp) + fine(ii,jj+1,0,1+fcomp)+fine(ii+1,jj,0,1+fcomp) + fine(ii+1,jj+1,0,1+fcomp));
                    crse(i,j,0,4+ccomp) = fine(ii+itake,jj+jtake,0,4+fcomp); //0.25*(fine(ii,jj,0,2+fcomp) + fine(ii,jj+1,0,2+fcomp)+fine(ii+1,jj,0,2+fcomp) + fine(ii+1,jj+1,0,2+fcomp));
                    crse(i,j,0,5+ccomp) = fine(ii+itake,jj+jtake,0,5+fcomp); //0.25*(fine(ii,jj,0,3+fcomp) + fine(ii,jj+1,0,3+fcomp)+fine(ii+1,jj,0,3+fcomp) + fine(ii+1,jj+1,0,3+fcomp));
                    crse(i,j,0,10+ccomp) = 0*(fine(ii,jj,0,10+fcomp) + fine(ii,jj+1,0,10+fcomp)+fine(ii+1,jj,0,10+fcomp) + fine(ii+1,jj+1,0,8+fcomp));
                    
                    crse(i,j,0,4+ccomp) = 0;
                    crse(i,j,0,5+ccomp) = 0;
                    crse(i,j,0,6+ccomp) = 0;
                    crse(i,j,0,7+ccomp) = 0;
                
            }
        }    
}
/*
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void adam_avgdown_nodes (Box const& bx, Array4<Real> const& crse,
                    Array4<Real const> const& fine,
                    int ccomp, int fcomp, int ncomp,
                    IntVect const& ratio) noexcept
{
    const auto clo = lbound(bx);
    const auto chi = ubound(bx);
    const int facx = ratio[0];
    const int facy = ratio[1];
    //const Real volfrac = static_cast<Real>(1.0)/static_cast<Real>(facx*facy);

    for (int n = 0; n < ncomp; ++n) {
        for (int j = clo.y; j <= chi.y; ++j) {
        for (int i = clo.x; i <= chi.x; ++i) {
            int ii = i*facx;
            int jj = j*facy;
            Real c = 0.;
            for (int jref = 0; jref < facy; ++jref) {
            for (int iref = 0; iref < facx; ++iref) {
                c += fine(ii+iref,jj+jref,0,n+fcomp);
            }}
            crse(i,j,0,n+ccomp) = c;
        }}
    }
}
*/
#endif
