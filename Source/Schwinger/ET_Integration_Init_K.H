#ifndef ET_INTEGRATION_INIT_K_H
#define ET_INTEGRATION_INIT_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
state_init(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& state_fab,
        amrex::Real cur_time, int lev, const amrex::GeometryData& geom) noexcept 
{
        const auto domain_xlo = geom.ProbLo(); 
        
        amrex::Real x0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1]; 
        
        amrex::Real dx0 = geom.CellSize(0);
        amrex::Real dx1 = geom.CellSize(1);
        
        srand (static_cast <unsigned> (time(0)));
        double r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        
        amrex::Real Theta_0 =  M_PI/4*(r-0.5); //0.1*std::sin(2*M_PI*(x0)/(128)); //0.1*std::exp(-(x0*x0 + x1*x1)/100);
        amrex::Real Theta_1 =  M_PI/4*(r-0.5); //0.1*std::sin(2*M_PI*(x0)/(128)); //0.1*std::exp(-(x0*x0 + x1*x1)/100);
        
        amrex::Real Pi_0 = amrex::RandomNormal(0,1.0); //-0.1*2*M_PI/64*std::sin(2*M_PI*(x0)/(128)); // -0.1*std::exp(-(x0*x0 + x1*x1)/100); 
        amrex::Real Pi_1 = amrex::RandomNormal(0,1.0); //-0.1*2*M_PI/64*std::sin(2*M_PI*(x0)/(128)); //-0.1*std::exp(-(x0*x0 + x1*x1)/100);
        
        state_fab(i, j, k, Idx::U_0_Real) = std::cos(Theta_0); 
        state_fab(i, j, k, Idx::U_0_Imaginary) = std::sin(Theta_0);
        state_fab(i, j, k, Idx::U_1_Real) = std::cos(Theta_1);
        state_fab(i, j, k, Idx::U_1_Imaginary) = std::sin(Theta_1);
        
        //state_fab(i, j, k, Idx::Pi_0) = Pi_0;
        //state_fab(i, j, k, Idx::Pi_1) = Pi_1;
        
        state_fab(i, j, k, Idx::Phi_0_Real) = 0; 
        state_fab(i, j, k, Idx::Phi_0_Imaginary) = 0;
        state_fab(i, j, k, Idx::Phi_1_Real) = 0;
        state_fab(i, j, k, Idx::Phi_1_Imaginary) = 0;
        
        
}

AMREX_GPU_DEVICE 
inline 
void 
aux_init(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& aux_fab, 
        amrex::Real time, const amrex::GeometryData& geom) noexcept 
{
        const auto domain_xlo = geom.ProbLo(); 
        
        amrex::Real x0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1]; 
        
        aux_fab(i, j, k, auxIdx::DDinvPhi_0_Real) = 0; 
        aux_fab(i, j, k, auxIdx::DDinvPhi_0_Imaginary) = 0;
        aux_fab(i, j, k, auxIdx::DDinvPhi_1_Real) = 0;
        aux_fab(i, j, k, auxIdx::DDinvPhi_1_Imaginary) = 0;
        
        aux_fab(i, j, k, auxIdx::g3DinvPhi_0_Real) = 0; 
        aux_fab(i, j, k, auxIdx::g3DinvPhi_0_Imaginary) = 0;
        aux_fab(i, j, k, auxIdx::g3DinvPhi_1_Real) = 0;
        aux_fab(i, j, k, auxIdx::g3DinvPhi_1_Imaginary) = 0;


        
        
}
#endif