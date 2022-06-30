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


AMREX_GPU_DEVICE 
inline 
void 
state_smear_gauge(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& smeared_fab,
        amrex::Array4<amrex::Real> const& smeared_tmp_fab,
        const amrex::Real time,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
        const amrex::GeometryData& geom,
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
        
        
        smeared_tmp_fab(i, j, k, cIdx::Real_0) = std::real(SmU_0_00);
        smeared_tmp_fab(i, j, k, cIdx::Imaginary_0) = std::imag(SmU_0_00);
        smeared_tmp_fab(i, j, k, cIdx::Real_1) = std::real(SmU_1_00);
        smeared_tmp_fab(i, j, k, cIdx::Imaginary_1) = std::imag(SmU_1_00);

}


AMREX_GPU_DEVICE 
inline 
void 
state_project_smeared(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& smeared_fab,
        amrex::Array4<amrex::Real> const& smeared_tmp_fab,
        const amrex::Real time,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
        const amrex::GeometryData& geom) noexcept 
{
        amrex::Real SmU_tmp_0_00_Real = smeared_tmp_fab(i, j, k, cIdx::Real_0);
        amrex::Real SmU_tmp_0_00_Imaginary = smeared_tmp_fab(i, j, k, cIdx::Imaginary_0);
        std::complex<double> SmU_tmp_0_00 = std::complex<double>(SmU_tmp_0_00_Real, SmU_tmp_0_00_Imaginary);
        
        amrex::Real SmU_tmp_1_00_Real = smeared_tmp_fab(i, j, k, cIdx::Real_1);
        amrex::Real SmU_tmp_1_00_Imaginary = smeared_tmp_fab(i, j, k, cIdx::Imaginary_1);
        std::complex<double> SmU_tmp_1_00 = std::complex<double>(SmU_tmp_1_00_Real, SmU_tmp_1_00_Imaginary);
        
        smeared_fab(i, j, k, cIdx::Real_0) = std::real(std::polar(1.0,std::arg(SmU_tmp_0_00)));
        smeared_fab(i, j, k, cIdx::Imaginary_0) = std::imag(std::polar(1.0,std::arg(SmU_tmp_0_00)));
        smeared_fab(i, j, k, cIdx::Real_1) = std::real(std::polar(1.0,std::arg(SmU_tmp_1_00)));
        smeared_fab(i, j, k, cIdx::Imaginary_1) = std::imag(std::polar(1.0,std::arg(SmU_tmp_1_00)));

}
#endif
