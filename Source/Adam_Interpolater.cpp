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
                //fine(i,j,0,n+fcomp) = crse(ic,jc,0,n+ccomp);
                if(n < 4)
                {
                    /*if(i%2+j%2 == 0)
                        fine(i,j,0,n+fcomp) = crse(ic,jc,0,n+ccomp);
                    else
                    {
                        if(n%2 == 0)
                            fine(i,j,0,n+fcomp) = 1;
                        else
                            fine(i,j,0,n+fcomp) = 0;
                    }*/
                    
                
                }
                else if(n == 4 || n == 5)
                {
                    fine(i,j,0,n+fcomp) = 0.5*crse(ic,jc,0,n+ccomp);
                    /*if(i%2+j%2 == 0)
                        fine(i,j,0,n+fcomp) = 0*crse(ic,jc,0,n+ccomp);
                    else
                        fine(i,j,0,n+fcomp) = 0;*/
                }
                else
                    fine(i,j,0,n+fcomp) = 0.5*crse(ic,jc,0,n+ccomp);
                    
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
    crse.grow(1);
    return crse;
}

Box
AdamInterp::CoarseBox (const Box&     fine,
                     const IntVect& ratio)
{
    Box crse = amrex::coarsen(fine,ratio);
    crse.grow(1);
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
