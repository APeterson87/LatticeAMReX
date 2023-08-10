#include <Adam_Interpolater.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_Interpolater.H>
#include <AMReX_Interp_C.H>
#include <AMReX_MFInterp_C.H>

#include <climits>

namespace amrex {
    
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
adaminterp_interp (Box const& bx,
                 Array4<Real> const& fine, const int fcomp, const int ncomp,
                 Array4<Real const> const& crse, const int ccomp,
                 IntVect const& ratio) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    //amrex::Print() << "fcomp = " << fcomp << std::endl;
    //amrex::Print() << "ccomp = " << ccomp << std::endl;
    for (int n = 0; n < ncomp; ++n) {
        for (int j = lo.y; j <= hi.y; ++j) {
            const int jc = amrex::coarsen(j,ratio[1]);
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                const int ic = amrex::coarsen(i,ratio[0]);
                
                int itakeA = Random_int(2);
                int jtakeA = Random_int(2);
                
                int itake = -2+Random_int(5);
                int jtake = -2+Random_int(5);
                
                
                
                if (n == 0 || n == 1 || n == 2 || n == 3)
                {
                    Real Uc_0_R_00 = crse(ic, jc, 0, 0+ccomp);
                    Real Uc_0_I_00 = crse(ic, jc, 0, 1+ccomp);
                    
                    Real Thetac_0_00 = std::atan2(Uc_0_I_00, Uc_0_R_00);
                    
                    Real Uc_0_R_0p = crse(ic, jc+1, 0, 0+ccomp);
                    Real Uc_0_I_0p = crse(ic, jc+1, 0, 1+ccomp);
                    
                    Real Thetac_0_0p = std::atan2(Uc_0_I_0p, Uc_0_R_0p);
                    
                    Real Uc_0_R_0m = crse(ic, jc-1, 0, 0+ccomp);
                    Real Uc_0_I_0m = crse(ic, jc-1, 0, 1+ccomp);
                    
                    Real Thetac_0_0m = std::atan2(Uc_0_I_0m, Uc_0_R_0m);
                    
                    Real Uc_1_R_p0 = crse(ic+1, jc, 0, 2+ccomp);
                    Real Uc_1_I_p0 = crse(ic+1, jc, 0, 3+ccomp);
                    
                    Real Thetac_1_p0 = std::atan2(Uc_1_I_p0, Uc_1_R_p0);
                    
                    Real Uc_1_R_00 = crse(ic, jc, 0, 2+ccomp);
                    Real Uc_1_I_00 = crse(ic, jc, 0, 3+ccomp);
                    
                    Real Thetac_1_00 = std::atan2(Uc_1_I_00, Uc_1_R_00);
                    
                    Real Uc_1_R_m0 = crse(ic-1, jc, 0, 2+ccomp);
                    Real Uc_1_I_m0 = crse(ic-1, jc, 0, 3+ccomp);
                    
                    Real Thetac_1_m0 = std::atan2(Uc_1_I_m0, Uc_1_R_m0);
                    
                    //Real Theta_0 = 1.0/6.0*(Thetac_0_0p + Thetac_0_00 + Thetac_0_0m);
                    //Real Theta_1 = 1.0/6.0*(Thetac_1_p0 + Thetac_1_00 + Thetac_1_m0);
                    
                    Real Theta_0 = Thetac_0_00;
                    Real Theta_1 = Thetac_1_00;
                    
                    if(i%2+j%2 == 0)
                    {
                        fine(i,j,0,0+fcomp) = std::cos(Theta_0);
                        fine(i,j,0,1+fcomp) = std::sin(Theta_0);
                        fine(i,j,0,2+fcomp) = std::cos(Theta_1);
                        fine(i,j,0,3+fcomp) = std::sin(Theta_1);
                    }
                    else
                    {
                        fine(i,j,0,0+fcomp) = 1; //std::cos(Theta_0);
                        fine(i,j,0,1+fcomp) = 0; //std::sin(Theta_0);
                        fine(i,j,0,2+fcomp) = 1; //std::cos(Theta_1);
                        fine(i,j,0,3+fcomp) = 0; //std::sin(Theta_1);
                    }
                    
                    
                    //fine(i,j,0,0+fcomp) = 1.0/6.0*(crse(ic,jc+1,0,0+ccomp)+crse(ic,jc,0,0+ccomp)+crse(ic,jc-1,0,0+ccomp));
                    //fine(i,j,0,1+fcomp) = 1.0/6.0*(crse(ic+1,jc,0,1+ccomp)+crse(ic,jc,0,1+ccomp)+crse(ic-1,jc,0,1+ccomp));
                    
                    /*if(i%2+j%2 == 0)
                    {
                        fine(i,j,0,0+fcomp) = crse(ic,jc,0,0+ccomp);
                        fine(i,j,0,1+fcomp) = crse(ic,jc,0,1+ccomp);
                    }
                    else
                    {
                        fine(i,j,0,0+fcomp) = 0;
                        fine(i,j,0,1+fcomp) = 0;
                    }*/
                        
                }
                
                else if (n == 4 || n == 5) fine(i,j,0,n+fcomp) = crse(ic+itake,jc+jtake,0,n+ccomp);
                
                //else
                  // fine(i,j,0,n+fcomp) = crse(ic,jc,0,n+ccomp);
                
                //if(n < 2)
                //{
                    /*
                    if(i%2+j%2 == 0)
                    {
                        Real dxC = 0.5;
                        Real dxF = 0.5*dxC;
                        Real thetaC_1 = dxC*crse(ic,jc,0,0+ccomp);
                        Real thetaC_2 = dxC*crse(ic+1,jc,0,1+ccomp);
                        Real thetaC_3 = dxC*crse(ic,jc+1,0,0+ccomp);
                        Real thetaC_4 = dxC*crse(ic,jc,0,1+ccomp);
                        
                        Real Plaq = std::cos(thetaC_1 + thetaC_2 - thetaC_3 - thetaC_4);
                        
                        fine(i,j,0,1+fcomp) = 0.5*thetaC_1/dxF;
                        fine(i+1,j,0,1+fcomp) = 0.5*thetaC_1/dxF;
                        fine(i,j,0,0+fcomp) = 0.5*thetaC_4/dxF;
                        fine(i,j+1,0,0+fcomp) = 0.5*thetaC_4/dxF;
                        
                        Real y = 0.25*(thetaC_1+thetaC_2 - thetaC_3 -thetaC_4);
                        Real thetaF = std::acos(0.5*std::sqrt(2+(Plaq/std::cos(y))));
                        
                        fine(i+1,j,0,0+fcomp) = thetaF/dxF;
                        fine(i+1,j+1,0,1+fcomp) = thetaF/dxF;
                        fine(i+1,j+1,0,0+fcomp) = thetaF/dxF;
                        fine(i,j+1,0,1+fcomp) = thetaF/dxF;
                    }
                   */     
                    
                    /*
                    int itake = -1+2*amrex::Random_int(2);
                    int jtake = -1+2*amrex::Random_int(2);
                    
                    if(i%2 == 0)
                    {
                        fine(i,j,0,0+fcomp) = crse(ic-1,jc+jtake,0,0+ccomp);
                        
                        //fine(i,j,0,1+fcomp) = crse(ic+itake,jc,0,1+ccomp);
                        //if (j != hi.y) fine(i,j+1,0,1+fcomp) = crse(ic+itake,jc,0,1+ccomp);
                    }
                    else
                        fine(i,j,0,0+fcomp) = crse(ic,jc+jtake,0,0+ccomp);
                    */
                   
               // }
               else
                   fine(i,j,0,n+fcomp) = crse(ic,jc,0,n+ccomp);
                    
            }
        }
    }
   Array4<Real> tmp = fine;
   amrex::Real alpha = 0.5;
   for(int iter = 0; iter < 0; iter++)
   {
       for (int j = lo.y; j <= hi.y; ++j) {
           AMREX_PRAGMA_SIMD
           for (int i = lo.x; i <= hi.x; ++i) {
               tmp(i,j,0,fcomp) *= (1-alpha);
               tmp(i,j,0,1+fcomp) *= (1-alpha);
               
               tmp(i,j,0,fcomp) += alpha/2*(fine(i,j,0,1+fcomp)+fine(i,j+1,0,fcomp)-fine(i+1,j,0,1+fcomp));
               tmp(i,j,0,fcomp) += alpha/2*(-fine(i,j-1,0,1+fcomp) + fine(i,j-1,0,fcomp) + fine(i+1,j-1,0,1+fcomp));

               tmp(i,j,0,1+fcomp) += alpha/2*(fine(i,j,0,fcomp) + fine(i+1,j,0,fcomp+1) - fine(i,j+1,0,fcomp)); 
               tmp(i,j,0,1+fcomp) += alpha/2*(-fine(i-1,j,0,fcomp) + fine(i-1,j,0,fcomp+1) + fine(i-1,j+1,0,fcomp));//0.05*crse(ic,jc,0,n+ccomp);


           }

       }
       
       for (int j = lo.y; j <= hi.y; ++j) {
           AMREX_PRAGMA_SIMD
           for (int i = lo.x; i <= hi.x; ++i) {
               fine(i, j, 0, fcomp) = tmp(i, j, 0, fcomp);
               fine(i, j, 0, 1+fcomp) = tmp(i, j, 0, 1+fcomp);
           }
       }
               
   }
 
}
/*
 * PCInterp, NodeBilinear, FaceLinear, CellConservativeLinear, and
 * CellBilinear are supported for all dimensions on cpu and gpu.
 *
 * CellConservativeProtected only works in 2D and 3D on cpu and gpu.
 *
 * CellQuadratic only works in 2D and 3D on cpu and gpu.
 *
 * CellConservativeQuartic only works with ref ratio of 2 on cpu and gpu.
 *
 * FaceDivFree works in 2D and 3D on cpu and gpu.
 * The algorithm is restricted to ref ratio of 2.
 */

//
// CONSTRUCT A GLOBAL OBJECT OF EACH VERSION.
//
AdamInterp                adam_interp;


    
AdamInterp::~AdamInterp () {}

Box
AdamInterp::CoarseBox (const Box& fine,
                     int        ratio)
{
    Box crse = amrex::coarsen(fine,ratio);
    crse.grow(2);
    return crse;
}

Box
AdamInterp::CoarseBox (const Box&     fine,
                     const IntVect& ratio)
{
    Box crse = amrex::coarsen(fine,ratio);
    crse.grow(2);
    return crse;
}

void
AdamInterp::interp (const FArrayBox& crse,
                  int              crse_comp,
                  FArrayBox&       fine,
                  int              fine_comp,
                  int              ncomp,
                  const Box&       fine_region,
                  const IntVect&   ratio,
                  const Geometry& /*crse_geom*/,
                  const Geometry& /*fine_geom*/,
                  Vector<BCRec> const& /*bcr*/,
                  int               /*actual_comp*/,
                  int               /*actual_state*/,
                  RunOn             runon)
{
    BL_PROFILE("AdamInterp::interp()");

    Array4<Real const> const& crsearr = crse.const_array();
    Array4<Real> const& finearr = fine.array();;

    AMREX_LAUNCH_HOST_DEVICE_LAMBDA_FLAG ( runon, fine_region, tbx,
    {
        amrex::adaminterp_interp(tbx,finearr,fine_comp,ncomp,crsearr,crse_comp,ratio);
    });
}



}
