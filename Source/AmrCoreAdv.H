#ifndef AmrCoreAdv_H_
#define AmrCoreAdv_H_

#include <string>
#include <limits>
#include <memory>
#include <complex>
#include <thread>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_AmrCore.H>
#include <AMReX_BCRec.H>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#include "bc_fill.H"

#include "ET_Integration.H"
#include "AMReX_TimeIntegrator.H"

using namespace amrex;

namespace InterpType {
    enum {
        PCInterp = 0,
        NodeBilinear,
        CellConservativeLinear,
        CellBilinear,
        CellQuadratic,
        CellConservativeProtected,
        CellConservativeQuartic
    };
}

class AmrCoreAdv
    : public amrex::AmrCore
{
public:

    ////////////////
    // public member functions

    // constructor - reads in parameters from inputs file
    //             - sizes multilevel arrays and data structures
    AmrCoreAdv ();
    virtual ~AmrCoreAdv();

    // advance solution to final time
    void Evolve ();

    // initializes multilevel data
    void InitData ();

    // Make a new level using provided BoxArray and DistributionMapping and
    // fill with interpolated coarse level data.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
					 const amrex::DistributionMapping& dm) override;

    // Remake an existing level using provided BoxArray and DistributionMapping and
    // fill with existing fine and coarse data.
    // overrides the pure virtual function in AmrCore
    virtual void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
			      const amrex::DistributionMapping& dm) override;
    
    void ResetLevel (int lev, Real time, amrex::Vector<amrex::MultiFab>& hold_grid, bool grid_fractured);
    
    
    // Delete level data
    // overrides the pure virtual function in AmrCore
    virtual void ClearLevel (int lev) override;

    // Make a new level from scratch using provided BoxArray and DistributionMapping.
    // Only used during initialization.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
					  const amrex::DistributionMapping& dm) override;

    // tag all cells for refinement
    // overrides the pure virtual function in AmrCore
    virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;
    

    void BiCG_Solve(MultiFab& x_mf, MultiFab& b_mf, MultiFab& state_mf, bool useInitGuess, int lev, Real time, const amrex::Geometry& geom_lev, Parameters Param, bool grid_fractured);

private:

    ////////////////
    // private member functions

    // read in some parameters from inputs file
    void ReadParameters();

    // set covered coarse cells to be the average of overlying fine cells
    void AverageDown ();

    // more flexible version of AverageDown() that lets you average down across multiple levels
    void AverageDownTo (int crse_lev);
    
    void AverageDownNodal();
    
    void AverageDownToNodal(int crse_lev);
    
    template <typename FAB>
    void adam_down_nodal (const FabArray<FAB>& S_fine,
                             FabArray<FAB>& S_crse,
                             const IntVect& ratio,
                             int ngcrse = 0,
                             bool mfiter_is_definitely_safe=false);
    
    void FlipSigns(int lev, MultiFab& mf, int icomp, int ncomp);

    // compute a new multifab by coping in phi from valid region and filling ghost cells
    // works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
    void FillPatch (int lev, amrex::Real time, amrex::MultiFab& mf, int icomp, int ncomp,  bool grid_fractured);

    // compute a new multifab by coping in data from valid region and filling ghost cells
    // works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
    // unlike FillPatch, FillIntermediatePatch will use the supplied multifab instead of fine level data.
    // This is to support filling boundary cells at an intermediate time between old/new times
    // on the fine level when valid data at a specific time is already available (such as
    // at each RK stage when integrating between initial and final times at a given level).
    void FillIntermediatePatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp);

    // fill an entire multifab by interpolating from the coarser level
    // this comes into play when a new level of refinement appears
    void FillCoarsePatch (int lev, amrex::Real time, amrex::MultiFab& mf, int icomp, int ncomp);

    // utility to copy in data from old and/or new state into another multifab
    void GetData (int lev, amrex::Real time, amrex::Vector<amrex::MultiFab*>& data,
                  amrex::Vector<amrex::Real>& datatime);

    // advance a level by dt
    // includes a recursive call for finer levels
    void timeStepRegrid (int lev, amrex::Real time, int iteration, Parameters Param);
    void timeStepAdvance (int lev, amrex::Real time, Parameters Param);

    // advance a single level for a single time step
    void Advance (int lev, amrex::Real time, amrex::Real dt_lev, int ncycle, Parameters Param);

    // a wrapper for EstTimeStep(0
    void ComputeDt ();

    // compute dt from CFL considerations
    Real EstTimeStep (int lev, bool local=false) const;

    // get plotfile name
    std::string PlotFileName (int lev) const;

    // put together an array of multifabs for writing
    amrex::Vector<const amrex::MultiFab*> PlotFileMF () const;

    // set plotfile variables names
    amrex::Vector<std::string> PlotFileVarNames () const;

    // write plotfile to disk
    void WritePlotFile () const;
    
    // get plotfile name
    std::string DiagFileName (int lev) const;

    // set plotfile variables names
    amrex::Vector<std::string> DiagFileVarNames () const;

    // write plotfile to disk
    void ComputeAndWriteDiagnosticFile () const;

    // write checkpoint file to disk
    void WriteCheckpointFile () const;

    // read checkpoint file from disk
    void ReadCheckpointFile ();

    // Read the file passed to amr.restart and use it as an initial condition for
    // the current simulation. Supports a different number of components and
    // ghost cells.
    void InitializeFromFile ();

    // Initialize the new-time data at a level from the initial_data MultiFab
    void InitializeLevelFromData (int lev, const amrex::MultiFab& initial_data);

    // utility to skip to next line in Header
    static void GotoNextLine (std::istream& is);

    // Single level functions called by advance()
    void post_update (amrex::MultiFab& state_mf, const amrex::Real time, const amrex::Geometry& geom);
    //void MCMC_update (amrex::MultiFab& state_mf, const amrex::Real time, const amrex::Geometry& geom, bool state_is_accepted);
    void fill_rhs (amrex::MultiFab& rhs_mf, const amrex::MultiFab& state_mf, const amrex::Real time, const amrex::Geometry& geom, amrex::Real factor);
    void fill_state_diagnostics (amrex::MultiFab& diag_mf, const amrex::MultiFab& state_mf, const amrex::Real time_lev, const amrex::Geometry& geom) const;
    
    /* My Functions */    
    
    void Perturb (amrex::MultiFab& state, int lev, amrex::Real time, const amrex::Geometry& geom, Parameters Param);
    void StatePerturb (amrex::Vector<amrex::MultiFab>& state, Real time, const amrex::Vector<amrex::Geometry>& geom, Parameters Param, bool grid_fractured);
    void StateAction (amrex::Vector<amrex::MultiFab>& state, Real time, const amrex::Vector<amrex::Geometry>& geom, Parameters Param);
    
    void Trajectory(MultiFab& state_mf, int lev, const amrex::Real time, const Geometry& geom, Parameters Param);
    void StateTrajectory(amrex::Vector<amrex::MultiFab>& state, const amrex::Real time, const amrex::Vector<amrex::Geometry>& geom, bool enable_smearing, Parameters Param);
    
    void update_momentum (MultiFab& state_mf, int lev, const amrex::Real time, const amrex::Geometry& geom, Parameters Param, amrex::Real dtau);
    
    void update_gauge (MultiFab& state_mf, int lev, const amrex::Real time, const Geometry& geom, Parameters Param, amrex::Real dtau);
    void update_action (MultiFab& state_mf, int lev, const amrex::Real time, const Geometry& geom, Parameters Param);
    
    void update_topChargeDensity (MultiFab& state_mf, int lev, const amrex::Real time, const Geometry& geom, Parameters Param);
    void update_State_topChargeDensity (amrex::Vector<amrex::MultiFab>& state, Real time, const amrex::Vector<amrex::Geometry>& geom, Parameters Param);
    
    void set_gauge_from_theta (MultiFab& state_mf, int lev);
    void ResetDDinvPhi(amrex::MultiFab& state_mf);

    void Set_DDinvPhi(MultiFab& state_mf, MultiFab& p_mf, const Geometry& geom, amrex::Real alpha);

    void Set_res(MultiFab& res_mf, MultiFab& state_mf, const Geometry& geom, amrex::Real alpha);
    void Set_p(MultiFab& state_mf, MultiFab& res_mf, const Geometry& geom, amrex::Real beta);
    void SetPhiFromChi(amrex::MultiFab& state_mf, int lev, amrex::Real time, amrex::Real m_0, amrex::Real r);
    void Set_g3DinvPhi (MultiFab& state_mf, MultiFab& aux_mf, const Geometry& geom, int lev, amrex::Real time, amrex::Real m_0, amrex::Real r);
    void Set_DDp (MultiFab& DDp_mf, MultiFab& p_mf, MultiFab& state_mf, const Geometry& geom, amrex::Real m_0, amrex::Real r);
    
    void Set_g3p (MultiFab& g3p_mf, MultiFab& p_mf, int lev, Real time);
    void Set_g3Dp (MultiFab& g3Dp_mf, MultiFab& p_mf, MultiFab& state_mf, Parameters Param);
    
    void Set_g3Dp_new (MultiFab& g3Dp_mf, MultiFab& p_mf, MultiFab& U_mf, Parameters Param);
    
    void Set_Dp (MultiFab& Dp_mf, MultiFab& p_mf, MultiFab& U_mf, int lev, Real time, Parameters Param);
    
    void smear_gauge (MultiFab& smearedU_mf, MultiFab& U_mf, int lev, const amrex::Real time, Parameters Param);
    
    amrex::Real Total_Action(MultiFab& state_mf, int lev, const amrex::Geometry& geom, Parameters Param);
    amrex::Real Action_Gauge(amrex::MultiFab& state_mf, int lev, const amrex::Geometry& geom, amrex::Real beta);
    amrex::Real Action_Momentum(amrex::MultiFab& state_mf, int lev, const amrex::Geometry& geom);
    amrex::Real Action_Fermi(amrex::MultiFab& state_mf, amrex::MultiFab& aux_mf, int lev);
    //amrex::Real Action_Fermi_From_Chi (MultiFab& state_mf);
    
    
    Real Total_Action_Lev (amrex::Vector<amrex::MultiFab>& state, const amrex::Vector<amrex::Geometry>& geom, Parameters Param, int high_lev);
    Real Total_Gauge_Action_Lev(amrex::Vector<amrex::MultiFab>& state, const amrex::Vector<amrex::Geometry>& geom, Parameters Param, int high_lev);
    Real Total_Momentum_Action_Lev(amrex::Vector<amrex::MultiFab>& state, const amrex::Vector<amrex::Geometry>& geom, Parameters Param, int high_lev);
    Real Action_Gauge_Lev (amrex::Vector<amrex::MultiFab>& state, int lev, const amrex::Geometry& geom, Parameters Param);
    Real Action_Momentum_Lev (amrex::Vector<amrex::MultiFab>& state, int lev, const amrex::Geometry& geom);
    
    Real meas_TopCharge (MultiFab& state_mf, int lev, const amrex::Real time, const Geometry& geom, Parameters Param);
    
    Real meas_TopCharge_Lev (amrex::Vector<amrex::MultiFab>& state, int lev, const amrex::Real time, const amrex::Geometry& geom, Parameters Param);
    
    Real meas_Total_TopCharge_Lev(amrex::Vector<amrex::MultiFab>& state, const amrex::Real time, const amrex::Vector<amrex::Geometry>& geom, Parameters Param, int high_lev);
    
    void meas_SigmaLoops(MultiFab& U_mf, int lev, amrex::Real time, const Geometry& geom, Parameters Param, amrex::Vector<double> &sigma, std::ofstream& Sigma);
    void meas_WilsonLoops(MultiFab& U_mf, int lev, amrex::Real time, int step, const Geometry& geom, Parameters Param, amrex::Vector<double> &sigma, std::ofstream& Sigma);
    Real Measure_Plaq (MultiFab& U_mf, int lev);
    
    void PionCorrelation (MultiFab& state_mf, int lev, const Real time_lev, int step, const amrex::Geometry& geom, Parameters Param, amrex::Vector<double>& picorr, std::ofstream& PiCorr, bool grid_fractured);


    /* End of My Functions*/
    
    
    ////////////////
    
    //amrex::Real SumAction(amrex::MultiFab& state_mf);
    // private data members

    amrex::Vector<int> istep;      // which step?
    amrex::Vector<int> nsubsteps;  // how many substeps on each level?

    // keep track of old time, new time, and time step at each level
    amrex::Vector<amrex::Real> t_new;
    amrex::Vector<amrex::Real> t_old;
    amrex::Vector<amrex::Real> dt;

    // array of multifabs to store the gauge at each level of refinement
    // after advancing a level we use "swap".
    amrex::Vector<amrex::MultiFab> grid_new;
    amrex::Vector<amrex::MultiFab> grid_hold;

    // array of time integrator objects at each level of refinement
    //amrex::Vector<std::unique_ptr<TimeIntegrator<MultiFab> > > integrator;

    // BCRec is essentially a 2*DIM integer array storing the physical boundary
    // condition types at the lo/hi walls in each direction
    amrex::Vector<BCRec> bcs;

    // Interpolater to use between levels
    amrex::Interpolater* mapper = nullptr;

    ////////////////
    // runtime parameters

    // maximum number of steps and stop time
    int max_step = std::numeric_limits<int>::max();
    int num_therm_steps = 100;
    int start_meas_step = 200;
    
    amrex::Real stop_time = std::numeric_limits<amrex::Real>::max();
    
    int num_hmc_substeps = 10;
    Real hmc_tau = 1.0;
    
    amrex::Real coupling_beta = 1.0;
    amrex::Real mass_0 = 0.1;
    amrex::Real wilson_r = 1.0;
    amrex::Real Temp_T = 1.0;
    
    
    bool use_dynamical_fermi = true;
    amrex::Real BiCGThresh = 1e-5;
    int BiCG_Max_Iter = 100;
    bool Use_BiCG_Stab = false;
    
    amrex::Real APE_alpha = 0.5;
    int APE_iter = 2;
    
    int measSigma_Interval = 10;
    int Check_revTraj_Int = 5;
    
    Vector<int> numcells;
    
    Vector<Real> coupling_beta_vec;
    Vector<Real> mom_sigma_vec;
    Vector<Real> temp_vec;
    Vector<Real> hmc_tau_vec;

    // if >= 0 we restart from a checkpoint
    std::string restart_chkfile = "";

    // if this flag is set, then the amr.restart file is
    // used to set the solver's initial conditions.
    int restart_is_initial_data = 0;

    // advective cfl number - dt = cfl*dx/umax
    amrex::Real cfl = 0.7;
    amrex::Real epsilon = 0.1;

    // how often each level regrids the higher levels of refinement
    // (after a level advances that many time steps)
    int regrid_int = 2;

    // plotfile prefix and frequency
    std::string plot_file {"plt"};
    int plot_int = -1;

    // diagnostic prefix and frequency
    std::string diag_file {"diag"};
    int diag_int = -1;

    // checkpoint prefix and frequency
    std::string chk_file {"chk"};
    int chk_int = -1;

    // domain and boundary conditions
    amrex::Vector<int> domain_lo_bc_types;
    amrex::Vector<int> domain_hi_bc_types;

    // elliptic?
    int elliptic = 0;

    // interpolation type between levels
    int interpolation_type = 0;

    // If the grid has been regridded
    bool is_regridded = false;

    // tagging criteria
    amrex::Real time_lower = 500;
    amrex::Real time_upper = 600;
    int x0_lower = -8;
    int x0_upper = 8;
    int x1_lower = -8;
    int x1_upper = 8;
  
};

#endif
