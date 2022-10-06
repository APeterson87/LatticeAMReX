#ifndef ADAM_INTERPOLATER_H_
#define ADAM_INTERPOLATER_H_
#include <AMReX_Config.H>
#include <AMReX_Interpolater.H>

#include <AMReX_Extension.H>
#include <AMReX_GpuControl.H>
#include <AMReX_InterpBase.H>

namespace amrex {

class Geometry;
class FArrayBox;
class IArrayBox;

/**
* \brief Virtual base class for interpolaters.
*
* Specifies interpolater interface for coarse-to-fine interpolation in space.
*/


    
class AdamInterp
    :
    public Interpolater
{
public:

    /**
    * \brief The destructor.
    */
    virtual ~AdamInterp () override;

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    virtual Box CoarseBox (const Box& fine, int ratio) override;

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    virtual Box CoarseBox (const Box& fine, const IntVect& ratio) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    */
    virtual void interp (const FArrayBox& crse,
                         int              crse_comp,
                         FArrayBox&       fine,
                         int              fine_comp,
                         int              ncomp,
                         const Box&       fine_region,
                         const IntVect&   ratio,
                         const Geometry&  crse_geom,
                         const Geometry&  fine_geom,
                         Vector<BCRec> const&  bcr,
                         int              actual_comp,
                         int              actual_state,
                         RunOn            gpu_or_cpu) override;
};





//! CONSTRUCT A GLOBAL OBJECT OF EACH VERSION.
extern AMREX_EXPORT AdamInterp                adam_interp;    

}

#endif /*_INTERPOLATER_H_*/
