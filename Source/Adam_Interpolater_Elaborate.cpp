#include <Adam_Interpolater.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_Interpolater.H>
#include <AMReX_Interp_C.H>
#include <AMReX_MFInterp_C.H>

#include <climits>

namespace amrex {
    
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void adaminterp_interp (int i, int j, int, int n, Array4<Real> const& fine, int fcomp,
                          Array4<Real const> const& crse, int ccomp, IntVect const& ratio) noexcept
{
    int ic = amrex::coarsen(i,ratio[0]);
    int jc = amrex::coarsen(j,ratio[1]);
    
    int ioff = i - ic*ratio[0];
    int joff = j - jc*ratio[1];
    Real rxinv = Real(1.0) / Real(ratio[0]);
    Real ryinv = Real(1.0) / Real(ratio[1]);
    if (ioff != 0 && joff != 0) {
        // Node on a X-Y face
        fine(i,j,0,n+fcomp) = rxinv * ryinv *
            (crse(ic  ,jc  ,0,n+ccomp) * (ratio[0]-ioff) * (ratio[1]-joff) +
             crse(ic+1,jc  ,0,n+ccomp) * (         ioff) * (ratio[1]-joff) +
             crse(ic  ,jc+1,0,n+ccomp) * (ratio[0]-ioff) * (         joff) +
             crse(ic+1,jc+1,0,n+ccomp) * (         ioff) * (         joff));
    } else if (ioff != 0) {
        // Node on X line
        fine(i,j,0,n+fcomp) = rxinv*((ratio[0]-ioff)*crse(ic  ,jc,0,n+ccomp) +
                                     (         ioff)*crse(ic+1,jc,0,n+ccomp));
    } else if (joff != 0) {
        // Node on Y line
        fine(i,j,0,n+fcomp) = ryinv*((ratio[1]-joff)*crse(ic,jc  ,0,n+ccomp) +
                                     (         joff)*crse(ic,jc+1,0,n+ccomp));
    } else {
        // Node coincident with coarse node
        fine(i,j,0,n+fcomp) = crse(ic,jc,0,n+ccomp);
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
    Box b = amrex::coarsen(fine,ratio);

    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        if (b.length(i) < 2)
        {
            //
            // Don't want degenerate boxes.
            //
            b.growHi(i,1);
        }
    }

    return b;
}

Box
AdamInterp::CoarseBox (const Box&     fine,
                         const IntVect& ratio)
{
    Box b = amrex::coarsen(fine,ratio);

    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        if (b.length(i) < 2)
        {
            //
            // Don't want degenerate boxes.
            //
            b.growHi(i,1);
        }
    }

    return b;
}

void
AdamInterp::interp (const FArrayBox&  crse,
                      int               crse_comp,
                      FArrayBox&        fine,
                      int               fine_comp,
                      int               ncomp,
                      const Box&        fine_region,
                      const IntVect&    ratio,
                      const Geometry& /*crse_geom */,
                      const Geometry& /*fine_geom */,
                      Vector<BCRec> const& /*bcr*/,
                      int               /*actual_comp*/,
                      int               /*actual_state*/,
                      RunOn             runon)
{
    BL_PROFILE("NodeBilinear::interp()");

    Array4<Real const> const& crsearr = crse.const_array();
    Array4<Real> const& finearr = fine.array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon,fine_region,ncomp,i,j,k,n,
    {
        adaminterp_interp(i,j,k,n, finearr, fine_comp, crsearr, crse_comp, ratio);
    });
}



}
