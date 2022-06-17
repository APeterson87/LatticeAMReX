#ifndef ET_INTEGRATION_UTILS_K_H
#define ET_INTEGRATION_UTILS_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

using namespace amrex;

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

#endif