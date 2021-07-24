#include <AmrCoreAdv.H>

using namespace amrex;

// constructor - reads in parameters from inputs file
//             - sizes multilevel arrays and data structures
//             - initializes BCRec boundary condition object
AmrCoreAdv::AmrCoreAdv ()
{
    ReadParameters();

    // Geometry on all levels has been defined already.

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    int nlevs_max = max_level + 1;

    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    for (int lev = 1; lev <= max_level; ++lev) {
        if (elliptic)
            nsubsteps[lev] = MaxRefRatio(lev-1) * MaxRefRatio(lev-1);
        else
            nsubsteps[lev] = MaxRefRatio(lev-1);
    }

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, 1.e100);

    grid_new.resize(nlevs_max);
    grid_old.resize(nlevs_max);

    integrator.resize(nlevs_max);

    bcs.resize(Idx::NumScalars);

    // boundary conditions
    for (int n = 0; n < Idx::NumScalars; ++n)
    {
        for (int i = 0; i < AMREX_SPACEDIM; ++i)
        {
            // is_periodic overrides inputs in domain_(lo/hi)_bc_type
            if (Geom(0).isPeriodic(i))
            {
                bcs[n].setLo(i, BCType::int_dir);
                bcs[n].setHi(i, BCType::int_dir);
            }
            else
            {
                bcs[n].setLo(i, domain_lo_bc_types[i]);
                bcs[n].setHi(i, domain_hi_bc_types[i]);
                
            }
        }
    }

    // set interpolation type between levels
    if (interpolation_type == InterpType::CellConservativeLinear)
        mapper = &cell_cons_interp;
    else if (interpolation_type == InterpType::CellConservativeQuartic)
        mapper = &quartic_interp;
    else {
        amrex::Error("Unsupported interpolation type");
    }
}

AmrCoreAdv::~AmrCoreAdv ()
{
}

// advance solution to final time
void
AmrCoreAdv::Evolve ()
{
    Real cur_time = t_new[0];
    int last_plot_file_step = 0;
    int last_chk_file_step = 0;
    int last_diag_file_step = 0;
    srand (static_cast <unsigned> (time(0)));

    Observables Obs;
    
    Parameters Param{.beta = coupling_beta, 
                     .mass = mass_0, 
                     .r = wilson_r, 
                     .Temp = Temp_T, 
                     
                     .Nx = numcells[0],
                     .Ny = numcells[1],
                     .hmc_substeps = num_hmc_substeps,
                     .therm_steps = num_therm_steps,
                     .starting_meas_step = start_meas_step,
                     .tau = hmc_tau,
                     .BiCG_Thresh = BiCGThresh,
                     .BiCG_Max_Iters = BiCG_Max_Iter,
                     .APE_smear_iter = APE_iter, 
                     .APE_smear_alpha = APE_alpha};
    
    amrex::Print() << "\nbeta = " << Param.beta << std::endl;
    amrex::Print() << "Dirac Mass = " << Param.mass << std::endl;
    amrex::Print() << "Wilson Factor = " << Param.r << std::endl;
    amrex::Print() << "Temperature = " << Param.Temp << std::endl;
    
    amrex::Print() << "\nLattice Number x = " << Param.Nx << std::endl;
    amrex::Print() << "Lattice Number y = " << Param.Ny << std::endl;
    
    amrex::Print() << "\nHMC substeps = " << Param.hmc_substeps << std::endl;
    amrex::Print() << "HMC tau = " << Param.tau << std::endl;
    amrex::Print() << "HMC dtau = " << Param.tau/Param.hmc_substeps << std::endl;
    amrex::Print() << "Thermalization steps = " << Param.therm_steps << std::endl;
    amrex::Print() << "Starting measurements at step " << Param.starting_meas_step << std::endl; 
    
    amrex::Print() << "\nBiconjugate tolerance  = " << Param.BiCG_Thresh << std::endl;
    amrex::Print() << "Max Biconjugate steps = " << Param.BiCG_Max_Iters << std::endl;
    
    
    std::this_thread::sleep_for(std::chrono::seconds(3));

    for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step)
    {
        amrex::Print() << "\nCoarse STEP " << step+1 << " starts ..." << std::endl;

        ComputeDt();

        int lev = 0;
        int iteration = 1;
        


        timeStep(lev, cur_time, iteration, Param, Obs);

        cur_time += dt[0];

        amrex::Print() << "Coarse STEP " << step+1 << " ends." << " TIME = " << cur_time
                       << " DT = " << dt[0]  << std::endl;

        // sync up time
        for (lev = 0; lev <= finest_level; ++lev) {
            t_new[lev] = cur_time;
        }

        if (plot_int > 0 && (step+1) % plot_int == 0) {
            last_plot_file_step = step+1;
            WritePlotFile();
        }

        if (chk_int > 0 && (step+1) % chk_int == 0) {
            last_chk_file_step = step+1;
            WriteCheckpointFile();
        }

        // Write a diagnostic file?
        if (diag_int > 0 && (step+1) % diag_int == 0) {
            last_diag_file_step = step + 1;
            ComputeAndWriteDiagnosticFile();
        }

#ifdef AMREX_MEM_PROFILING
        {
            std::ostringstream ss;
            ss << "[STEP " << step+1 << "]";
            MemProfiler::report(ss.str());
        }
#endif

        if (cur_time >= stop_time - 1.e-6*dt[0]) break;
    }

    if (plot_int > 0 && istep[0] > last_plot_file_step) {
        WritePlotFile();
    }

    if (chk_int > 0 && istep[0] > last_chk_file_step) {
        WriteCheckpointFile();
    }
    
    if (diag_int > 0 && istep[0] > last_diag_file_step) {
        ComputeAndWriteDiagnosticFile();
    }
    
    amrex::Print() << "\nFinal Acceptance Rate = " << static_cast<float>(Obs.num_accepted)/static_cast<float>(max_step-start_meas_step) << std::endl;
    
    amrex::Print() << "\nTime average Plaq = " << Obs.TimeSumAvePlaq/static_cast<float>(max_step-start_meas_step) << std::endl;
    
    amrex::Print() << "\nTime average deltaH = " << Obs.TimeSumDeltaH/static_cast<float>(max_step-start_meas_step) << std::endl;
    
    std::ofstream TopCharge("TopologicalCharges.dat");
    for (const auto &e : Obs.TopologicalCharge) TopCharge << e << "\n";
    TopCharge.close();
    
    std::ofstream DeltaH("DeltaH.dat");
    for (const auto &e : Obs.deltaH) DeltaH << e << "\n";
    DeltaH.close();
}

// initializes multilevel data
void
AmrCoreAdv::InitData ()
{
    if (restart_chkfile == "") {
        // start simulation from the beginning
        const Real time = 0.0;
        InitFromScratch(time);
        AverageDown();

        if (chk_int > 0) {
            WriteCheckpointFile();
        }
    }
    else if (restart_is_initial_data) {
        // Load the initial data from a file
        InitializeFromFile();
    }
    else {
        // restart from a checkpoint
        ReadCheckpointFile();
    }

    if (plot_int > 0) {
        WritePlotFile();
    }
    
    if (diag_int > 0) {
        ComputeAndWriteDiagnosticFile();
    }
}

// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data.
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
				                    const DistributionMapping& dm)
{
    const int ncomp = grid_new[lev-1].nComp();
    const int nghost = grid_new[lev-1].nGrow();

    grid_new[lev].define(ba, dm, ncomp, nghost);
    grid_old[lev].define(ba, dm, ncomp, nghost);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    FillCoarsePatch(lev, time, grid_new[lev], 0, ncomp);

    // also create the time integrator for this level
    integrator[lev] = std::make_unique<TimeIntegrator<MultiFab> >(grid_old[lev]);
}

// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::RemakeLevel (int lev, Real time, const BoxArray& ba,
			             const DistributionMapping& dm)
{
    const int ncomp = grid_new[lev].nComp();
    const int nghost = grid_new[lev].nGrow();

    MultiFab new_state(ba, dm, ncomp, nghost);
    MultiFab old_state(ba, dm, ncomp, nghost);

    FillPatch(lev, time, new_state, 0, ncomp);

    std::swap(new_state, grid_new[lev]);
    std::swap(old_state, grid_old[lev]);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    // also recreate the time integrator for this level
    integrator[lev] = std::make_unique<TimeIntegrator<MultiFab> >(grid_old[lev]);
}

// Delete level data
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::ClearLevel (int lev)
{
    grid_new[lev].clear();
    grid_old[lev].clear();
}

// Make a new level from scratch using provided BoxArray and DistributionMapping.
// Only used during initialization.
// overrides the pure virtual function in AmrCore
void AmrCoreAdv::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
					                      const DistributionMapping& dm)
{
    const int ncomp = Idx::NumScalars;
    const int nghost = NUM_GHOST_CELLS;

    grid_new[lev].define(ba, dm, ncomp, nghost);
    grid_old[lev].define(ba, dm, ncomp, nghost);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    Real cur_time = t_new[lev];
    MultiFab& state = grid_new[lev];

    const auto& geom = Geom(lev);
    const auto dx = geom.CellSizeArray();

    // Loop over grids at this level
#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( MFIter mfi(state, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      const Box& bx = mfi.tilebox();

      const auto& state_arr = state.array(mfi);

      // For each grid, loop over all the valid points
      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        state_init(i, j, k, state_arr, cur_time, geom.data());
      });
    }

    // also create the time integrator for this level
    integrator[lev] = std::make_unique<TimeIntegrator<MultiFab> >(grid_old[lev]);
}

// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow)
{
    static bool first = true;
    static Vector<Real> s_error;
    static int error_comp;
    
    const amrex::Geometry& geom = Geom(lev);
    const auto dx = geom.CellSizeArray();

    // only do this during the first call to ErrorEst
    if (first)
    {
        first = false;

        // in subroutine state_error, you could use more elaborate tagging, such
        // as more advanced logical expressions, or gradients, etc.
        ParmParse pp("problem");
        int n = pp.countval("s_error");
        if (n > 0) {
            pp.getarr("s_error", s_error, 0, n);
        }

        pp.get("error_comp", error_comp);
    }

    if (lev >= s_error.size()) return;

    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    const MultiFab& state = grid_new[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        for (MFIter mfi(state, true); mfi.isValid(); ++mfi)
        {
            const Box& tilebox  = mfi.tilebox();

            auto tagarr = tags.array(mfi);
            const auto& state_fab = state.array(mfi);
            const Real error_threshold = s_error[lev];

            // For each grid, loop over all the valid points
            amrex::ParallelFor(tilebox,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {          
                tagarr(i, j, k) = state_is_tagged(i, j, k, lev, state_fab, error_threshold, time, dx, geom.data()) ? tagval : clearval;
            });
        }
    }
}

// read in some parameters from inputs file
void
AmrCoreAdv::ReadParameters ()
{
    {
        ParmParse pp;  // Traditionally, max_step and stop_time do not have prefix.
        pp.query("max_step", max_step);
        
        pp.get("num_therm_steps", num_therm_steps);
        pp.get("starting_measurements_step", start_meas_step); 
        
        pp.query("stop_time", stop_time);
        pp.query("num_hmc_substeps", num_hmc_substeps);
        pp.query("hmc_tau", hmc_tau);
        pp.get("coupling_beta", coupling_beta);
        pp.get("mass_0", mass_0);
        pp.get("wilson_r", wilson_r);
        pp.get("Temp_T", Temp_T);
        pp.get("BiCGThresh", BiCGThresh);
        pp.get("BiCG_Max_Iter", BiCG_Max_Iter);
        
        
        pp.get("APE_alpha", APE_alpha);
        pp.get("APE_iter", APE_iter);
        

        // Query domain periodicity
        pp.getarr("domain_lo_bc_types", domain_lo_bc_types);
        pp.getarr("domain_hi_bc_types", domain_hi_bc_types);
    }

    {
        ParmParse pp("amr"); // Traditionally, these have prefix, amr.
        
        pp.getarr("n_cell", numcells);

        pp.query("interpolation_type", interpolation_type);

        pp.query("regrid_int", regrid_int);
        pp.query("plot_file", plot_file);
        pp.query("plot_int", plot_int);
        pp.query("chk_file", chk_file);
        pp.query("chk_int", chk_int);

        pp.query("restart", restart_chkfile);
        pp.query("restart_is_initial_data", restart_is_initial_data);

        // Diagnostics
        // Default diag_int to -1, allow us to set it to something else in the inputs file
        // If diag_int < 0 then no diagnostic files will be written
        diag_file = "diag";
        diag_int = -1;
        pp.query("diag_file", diag_file);
        pp.query("diag_int", diag_int);
    }

    {
        ParmParse pp("problem");

        // CFL
        pp.get("cfl", cfl);

        // Elliptic?
        pp.get("elliptic", elliptic);
        pp.get("epsilon", epsilon);
    }
}

// set covered coarse cells to be the average of overlying fine cells
void
AmrCoreAdv::AverageDown ()
{
    for (int lev = finest_level-1; lev >= 0; --lev)
    {
        amrex::average_down(grid_new[lev+1], grid_new[lev],
                            geom[lev+1], geom[lev],
                            0, grid_new[lev].nComp(), refRatio(lev));
    }
}

// more flexible version of AverageDown() that lets you average down across multiple levels
void
AmrCoreAdv::AverageDownTo (int crse_lev)
{
    amrex::average_down(grid_new[crse_lev+1], grid_new[crse_lev],
                        geom[crse_lev+1], geom[crse_lev],
                        0, grid_new[crse_lev].nComp(), refRatio(crse_lev));
}

void 
AmrCoreAdv::FlipSigns(int lev, MultiFab& mf, int icomp, int ncomp)
{
  const auto& geom = Geom(lev);  
  const Box& bx = geom.Domain();
    
#ifdef _OPENMP
#pragma omp parallel
#endif
  for ( MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
  {
    const Box& gbx = mfi.growntilebox();

    
    const auto Total_comps = mf.nComp();
    AMREX_ASSERT(icomp + ncomp - 1 <= Total_comps);

    const auto& fab = mf.array(mfi);

    // For each grid, loop over all the points including ghost cells
    amrex::ParallelFor(gbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        //Only flip sign if j is outside bx
        if(j < bx.smallEnd(1) || j > bx.bigEnd(1))
        {
            for(int n = icomp; n <= icomp+ncomp - 1; n++)
            {
                fab(i, j, k, n) *= -1;
            }
        }
    });
  }
}

// compute a new multifab by coping in data from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void
AmrCoreAdv::FillPatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
    if (lev == 0)
    {
        Vector<MultiFab*> smf;
        Vector<Real> stime;
        GetData(0, time, smf, stime);

        if(Gpu::inLaunchRegion())
        {
            GpuBndryFuncFab<AmrCoreFillGpu> gpu_bndry_func(AmrCoreFillGpu{});
            PhysBCFunct<GpuBndryFuncFab<AmrCoreFillGpu> > physbc(geom[lev],bcs,gpu_bndry_func);
            //! make sure this is passing the right bc
            amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp, 
                                        geom[lev], physbc, 0);
        }
        else
        {
            CpuBndryFuncFab bndry_func(&AmrCoreFillCpu);  // Without EXT_DIR, we can pass a nullptr.
            PhysBCFunct<CpuBndryFuncFab> physbc(geom[lev],bcs,bndry_func);
            //! make sure this is passing the right bc
            amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp, 
                                        geom[lev], physbc, 0);
        }
    }
    else
    {
        Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetData(lev-1, time, cmf, ctime);
        GetData(lev  , time, fmf, ftime);

        if(Gpu::inLaunchRegion())
        {
            GpuBndryFuncFab<AmrCoreFillGpu> gpu_bndry_func(AmrCoreFillGpu{});
            PhysBCFunct<GpuBndryFuncFab<AmrCoreFillGpu> > cphysbc(geom[lev-1],bcs,gpu_bndry_func);
            PhysBCFunct<GpuBndryFuncFab<AmrCoreFillGpu> > fphysbc(geom[lev],bcs,gpu_bndry_func);

            //! make sure this is passing the right bc
            amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                      0, icomp, ncomp, geom[lev-1], geom[lev],
                                      cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                      mapper, bcs, 0);
        }
        else
        {
            CpuBndryFuncFab bndry_func(&AmrCoreFillCpu);  // Without EXT_DIR, we can pass a nullptr.
            PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bndry_func);
            PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev],bcs,bndry_func);

            //! make sure this is passing the right bc
            amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                      0, icomp, ncomp, geom[lev-1], geom[lev],
                                      cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                      mapper, bcs, 0);
        }
    }

}

// compute a new multifab by coping in data from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
// unlike FillPatch, FillIntermediatePatch will use the supplied multifab instead of fine level data.
// This is to support filling boundary cells at an intermediate time between old/new times
// on the fine level when valid data at a specific time is already available (such as
// at each RK stage when integrating between initial and final times at a given level).
void
AmrCoreAdv::FillIntermediatePatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
    if (lev == 0)
    {
        // on lev, use the mf data and time passed to FillIntermediatePatch().
        Vector<MultiFab*> smf { &mf };
        Vector<Real> stime { time };

        if(Gpu::inLaunchRegion())
        {
            GpuBndryFuncFab<AmrCoreFillGpu> gpu_bndry_func(AmrCoreFillGpu{});
            PhysBCFunct<GpuBndryFuncFab<AmrCoreFillGpu> > physbc(geom[lev],bcs,gpu_bndry_func);
            //! make sure this is passing the right bc
            amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
                                        geom[lev], physbc, 0);
        }
        else
        {
            CpuBndryFuncFab bndry_func(&AmrCoreFillCpu);  // Without EXT_DIR, we can pass a nullptr.
            PhysBCFunct<CpuBndryFuncFab> physbc(geom[lev],bcs,bndry_func);
            //! make sure this is passing the right bc
            amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
                                        geom[lev], physbc, 0);
        }
    }
    else
    {
        MultiFab mf_temp(mf.boxArray(), mf.DistributionMap(), mf.nComp(), NUM_GHOST_CELLS);

        Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;

        // interpolate coarse data in time
        MultiFab mf_coarse_temp(grid_old[lev-1].boxArray(), grid_old[lev-1].DistributionMap(),
                                grid_old[lev-1].nComp(), NUM_GHOST_CELLS);

        // fill mf_coarse_temp using MC Equation 39 at "time"
        Real timestep_fraction = (time - t_old[lev-1])/dt[lev-1];
        integrator[lev-1]->time_interpolate(grid_new[lev-1], grid_old[lev-1], timestep_fraction, mf_coarse_temp);

        // we'll want to interpolate from mf_coarse_temp at "time"
        cmf.push_back(&mf_coarse_temp);
        ctime.push_back(time);

        // on lev, use the mf data and time passed to FillIntermediatePatch().
        fmf.push_back(&mf);
        ftime.push_back(time);

        if(Gpu::inLaunchRegion())
        {
            GpuBndryFuncFab<AmrCoreFillGpu> gpu_bndry_func(AmrCoreFillGpu{});
            PhysBCFunct<GpuBndryFuncFab<AmrCoreFillGpu> > cphysbc(geom[lev-1],bcs,gpu_bndry_func);
            PhysBCFunct<GpuBndryFuncFab<AmrCoreFillGpu> > fphysbc(geom[lev],bcs,gpu_bndry_func);

            //! make sure this is passing the right bc
            amrex::FillPatchTwoLevels(mf_temp, time, cmf, ctime, fmf, ftime,
                                      0, icomp, ncomp, geom[lev-1], geom[lev],
                                      cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                      mapper, bcs, 0);
        }
        else
        {
            CpuBndryFuncFab bndry_func(&AmrCoreFillCpu);  // Without EXT_DIR, we can pass a nullptr.
            PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bndry_func);
            PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev],bcs,bndry_func);

            //! make sure this is passing the right bc
            amrex::FillPatchTwoLevels(mf_temp, time, cmf, ctime, fmf, ftime,
                                      0, icomp, ncomp, geom[lev-1], geom[lev],
                                      cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                      mapper, bcs, 0);
        }

        // Replace mf with mf_temp
        std::swap(mf_temp, mf);
    }

}

// fill an entire multifab by interpolating from the coarser level
// this comes into play when a new level of refinement appears
void
AmrCoreAdv::FillCoarsePatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
    BL_ASSERT(lev > 0);

    Vector<MultiFab*> cmf;
    Vector<Real> ctime;
    GetData(lev-1, time, cmf, ctime);

    if (cmf.size() != 1) {
	amrex::Abort("FillCoarsePatch: how did this happen?");
    }

    if(Gpu::inLaunchRegion())
    {
        GpuBndryFuncFab<AmrCoreFillGpu> gpu_bndry_func(AmrCoreFillGpu{});
        PhysBCFunct<GpuBndryFuncFab<AmrCoreFillGpu> > cphysbc(geom[lev-1],bcs,gpu_bndry_func);
        PhysBCFunct<GpuBndryFuncFab<AmrCoreFillGpu> > fphysbc(geom[lev],bcs,gpu_bndry_func);

        //! make sure this is passing the right bc
        amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
                                     cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                     mapper, bcs, 0);
    }
    else
    {
        CpuBndryFuncFab bndry_func(&AmrCoreFillCpu);  // Without EXT_DIR, we can pass a nullptr.
        PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bndry_func);
        PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev],bcs,bndry_func);

        //! make sure this is passing the right bc
        amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
                                     cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                     mapper, bcs, 0);
    }

}

// utility to copy in data from old/new data into another multifab
void
AmrCoreAdv::GetData (int lev, Real time, Vector<MultiFab*>& data, Vector<Real>& datatime)
{
    data.clear();
    datatime.clear();

    const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

    if (time > t_new[lev] - teps && time < t_new[lev] + teps)
    {
        data.push_back(&grid_new[lev]);
        datatime.push_back(t_new[lev]);
    }
    else if (time > t_old[lev] - teps && time < t_old[lev] + teps)
    {
        data.push_back(&grid_old[lev]);
        datatime.push_back(t_old[lev]);
    }
    else
    {
        data.push_back(&grid_old[lev]);
        data.push_back(&grid_new[lev]);
        datatime.push_back(t_old[lev]);
        datatime.push_back(t_new[lev]);
    }
}


// advance a level by dt
// includes a recursive call for finer levels
void
AmrCoreAdv::timeStep (int lev, Real time, int iteration, Parameters Param, Observables& Obs)
{
    if (regrid_int > 0)  // We may need to regrid
    {
        // help keep track of whether a level was already regridded
        // from a coarser level call to regrid
        static Vector<int> last_regrid_step(max_level+1, 0);

        // regrid changes level "lev+1" so we don't regrid on max_level
        // also make sure we don't regrid fine levels again if
        // it was taken care of during a coarser regrid
        if (lev < max_level && istep[lev] > last_regrid_step[lev])
        {
            if (istep[lev] % regrid_int == 0)
            {
                // regrid could add newly refine levels (if finest_level < max_level)
                // so we save the previous finest level index
                int old_finest = finest_level;
                regrid(lev, time);

                // mark that we have regridded this level already
                for (int k = lev; k <= finest_level; ++k) {
                    last_regrid_step[k] = istep[k];
                }

                // if there are newly created levels, set the time step
                for (int k = old_finest+1; k <= finest_level; ++k) {
                    dt[k] = dt[k-1] / MaxRefRatio(k-1);
                }
            }
        }
    }

    if (Verbose()) {
        amrex::Print() << "[Level " << lev << " step " << istep[lev]+1 << "] ";
        amrex::Print() << "ADVANCE with time = " << t_new[lev]
                       << " dt = " << dt[lev] << std::endl;
    }

    // advance a single level for a single time step, updates flux registers
    Advance(lev, time, dt[lev], iteration, nsubsteps[lev], Param, Obs);

    ++istep[lev];

    if (Verbose())
    {
        amrex::Print() << "[Level " << lev << " step " << istep[lev] << "] ";
        amrex::Print() << "Advanced " << CountCells(lev) << " cells" << std::endl;
    }

    if (lev < finest_level)
    {
        // recursive call for next-finer level
        for (int i = 1; i <= nsubsteps[lev+1]; ++i)
        {
            timeStep(lev+1, time+(i-1)*dt[lev+1], i, Param, Obs);
        }

        AverageDownTo(lev); // average lev+1 down to lev
    }
}

// advance a single level for a single time step, updates flux registers
void
AmrCoreAdv::Advance (int lev, Real time, Real dt_lev, int iteration, int ncycle, Parameters Param, Observables& Obs)
{
    constexpr int num_grow = NUM_GHOST_CELLS;

    t_old[lev] = t_new[lev];
    t_new[lev] += dt_lev;

    const Real old_time = t_old[lev];
    const Real new_time = t_new[lev];
    
    Real dtau = 1.0/num_hmc_substeps;
        
    Real HGaugecurrent;
    Real HMomcurrent;
    Real HFermicurrent;
    Real HTotalcurrent;
    Real HGauge;
    Real HMom;
    Real HFermi;
    Real HTotal;
    
    
    MultiFab::Copy(grid_old[lev],grid_new[lev], 0, 0, Idx::NumScalars, NUM_GHOST_CELLS);
    MultiFab& S_new = grid_new[lev];
    

    // State with ghost cells as the integrator initial condition
    //MultiFab Sborder(grids[lev], dmap[lev], S_new.nComp(), num_grow);
    //FillPatch(lev, time, Sborder, 0, Sborder.nComp());
    
    //MultiFab& Sborder = grid_old[lev];
    //FillPatch(lev, time, Sborder, 0, Sborder.nComp());

    /* Integrate from (y,t) = (Sborder, time) by dt_lev to set S_new. */
    const auto geom_lev = geom[lev];
    
    FillIntermediatePatch(lev, time, S_new, 0, S_new.nComp());
    FlipSigns(lev, S_new, Idx::Phi_0_Real, 4);
    FlipSigns(lev, S_new, Idx::DDinvPhi_0_Real, 4);
    FlipSigns(lev, S_new, Idx::chi_0_Real, 4);
    FlipSigns(lev, S_new, Idx::g3DinvPhi_0_Real, 4);
    
    SetPhiFromChi(S_new, lev, time, mass_0, wilson_r);
    
    Fermi_BiCG_solve(S_new, lev, time, geom_lev, Param);
    
    HGaugecurrent = Action_Gauge(S_new, Param.beta);
    HMomcurrent = Action_Momentum(S_new);
    HFermicurrent = Action_Fermi_From_Chi(S_new);
    HTotalcurrent = HGaugecurrent + HMomcurrent + HFermicurrent;

    Trajectory(S_new, lev,  time, geom_lev, Param);
    
    if((istep[0]+1)%5 == 0  && istep[0] <= Param.therm_steps)
    {
        Real L_2 = Check_Reversed_Trajectory(S_new, lev, time, geom_lev, Param);
        amrex::Print() << "The L2 norm is " << L_2/(2*numcells[0]*numcells[1]) << std::endl;
    }
    
    
    HGauge = Action_Gauge(S_new, Param.beta);
    HMom = Action_Momentum(S_new);
    HFermi = Action_Fermi(S_new);
    HTotal = HGauge + HMom + HFermi;
    
    Real r_loc = std::rand()/(static_cast <float> (RAND_MAX));
    
    amrex::Print() << "MCMC random number = " << r_loc << std::endl;
    amrex::Print() << "Exp(-deltaH/T) = " << std::exp(-(HTotal-HTotalcurrent)/Temp_T) << std::endl;
    amrex::Print() << "deltaH = " << HTotal - HTotalcurrent << std::endl;
    amrex::Print() << "Total H = " << HTotalcurrent << std::endl;

    if(r_loc > std::exp(-(HTotal-HTotalcurrent)/Temp_T) && istep[0] >= Param.therm_steps)
    {
        MultiFab::Copy(S_new, S_new,  Idx::Ucurrent_0_Real, Idx::U_0_Real, 4, NUM_GHOST_CELLS);
        amrex::Print() << "NEW STATE REJECTED " << std::endl;
    }
    else
    {
        MultiFab::Copy(S_new, S_new, Idx::U_0_Real, Idx::Ucurrent_0_Real, 4, NUM_GHOST_CELLS);
        amrex::Print() << "NEW STATE ACCEPTED " << std::endl;
        if(istep[0] > Param.starting_meas_step)
            Obs.num_accepted += 1;
    }
    
    if(istep[0] >= Param.starting_meas_step)
    {
        Obs.TimeSumAvePlaq += Measure_Plaq(S_new);
        Obs.TimeSumDeltaH += HTotal-HTotalcurrent;
        Real InstantonNumber = meas_TopCharge(S_new, lev, time, geom_lev, Param);
        Obs.TopologicalCharge.push_back((int)std::round(InstantonNumber));
        amrex::Print() << "Topological Charge = " << (int)std::round(InstantonNumber) << std::endl;
        
        Obs.deltaH.push_back(HTotal-HTotalcurrent);
        
        
    }
    
    Perturb(S_new, lev, time, geom_lev);
    
}


// a wrapper for EstTimeStep
void
AmrCoreAdv::ComputeDt ()
{
    Vector<Real> dt_tmp(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        dt_tmp[lev] = EstTimeStep(lev, true);
    }

    ParallelDescriptor::ReduceRealMin(&dt_tmp[0], dt_tmp.size());

    constexpr Real change_max = 1.1;
    Real dt_0 = dt_tmp[0];
    int n_factor = 1;
    for (int lev = 0; lev <= finest_level; ++lev) {
        dt_tmp[lev] = std::min(dt_tmp[lev], change_max*dt[lev]);
        n_factor *= nsubsteps[lev];
        dt_0 = std::min(dt_0, n_factor*dt_tmp[lev]);
    }

    // Limit dt's by the value of stop_time.
    const Real eps = 1.e-3*dt_0;
    if (t_new[0] + dt_0 > stop_time - eps) {
        dt_0 = stop_time - t_new[0];
    }

    dt[0] = dt_0;
    for (int lev = 1; lev <= finest_level; ++lev) {
        dt[lev] = dt[lev-1] / nsubsteps[lev];
    }
}

// compute dt from CFL considerations
Real
AmrCoreAdv::EstTimeStep (int lev, bool local) const
{
    BL_PROFILE("AmrCoreAdv::EstTimeStep()");

    Real dt = std::numeric_limits<Real>::max();

    const auto dx = geom[lev].CellSizeArray();

    if (elliptic)
    {
        Print() << "using elliptic timestep\n";
        dt = dx[0]*dx[0];
#if AMREX_SPACEDIM > 1
        dt = amrex::min(dt, dx[1]*dx[1]);
#endif
#if AMREX_SPACEDIM > 2
        dt = amrex::min(dt, dx[2]*dx[2]);
#endif
    }
    else
    {
        dt = dx[0];
#if AMREX_SPACEDIM > 1
        dt = amrex::min(dt, dx[1]);
#endif
#if AMREX_SPACEDIM > 2
        dt = amrex::min(dt, dx[2]);
#endif
    }

    dt *= cfl;

    return dt;
}

// get plotfile name
std::string
AmrCoreAdv::PlotFileName (int lev) const
{
    return amrex::Concatenate(plot_file, lev, 7);
}

// put together an array of multifabs for writing
Vector<const MultiFab*>
AmrCoreAdv::PlotFileMF () const
{
    Vector<const MultiFab*> r;
    for (int i = 0; i <= finest_level; ++i) {
        r.push_back(&grid_new[i]);
    }
    return r;
}

// set plotfile variable names
Vector<std::string>
AmrCoreAdv::PlotFileVarNames () const
{
    return Variable::names;
}

// write plotfile to disk
void
AmrCoreAdv::WritePlotFile () const
{
    const std::string& plotfilename = PlotFileName(istep[0]);
    const auto& mf = PlotFileMF();
    const auto& varnames = PlotFileVarNames();

    amrex::Print() << "Writing plotfile " << plotfilename << "\n";

    amrex::WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf, varnames,
                                   Geom(), t_new[0], istep, refRatio());
}

std::string
AmrCoreAdv::DiagFileName (int lev) const
{
    return amrex::Concatenate(diag_file, lev, 7);
}

Vector<std::string>
AmrCoreAdv::DiagFileVarNames () const
{
    return Diagnostics::names;
}

void
AmrCoreAdv::ComputeAndWriteDiagnosticFile () const
{
    const std::string& diagfilename = DiagFileName(istep[0]);
    //const auto& mf = PlotFileMF();
    const int nDiagComp = Diag::NumScalars;
    Vector<MultiFab>  grid_diag(finest_level+1);
    Vector<const MultiFab*> grid_diag_mf;
    
    // sets the size and distribution of diagnostic data
    for (int lev = 0; lev <= finest_level; ++lev) {
        grid_diag[lev].define(boxArray(lev), DistributionMap(lev), nDiagComp, 0);
        grid_diag_mf.push_back(&grid_diag[lev]);
    }
    
    for (int lev = 0; lev <= finest_level; ++lev) {  
        const amrex::Geometry& geom = Geom(lev);
        fill_state_diagnostics(grid_diag[lev], grid_new[lev], t_new[lev], geom);
        
    }
    
    const auto& diagnames = DiagFileVarNames();

    amrex::Print() << "Writing diagfile " << diagfilename << "\n";

    amrex::WriteMultiLevelPlotfile(diagfilename, finest_level+1, grid_diag_mf, diagnames,
                                   Geom(), t_new[0], istep, refRatio());
}

void
AmrCoreAdv::WriteCheckpointFile () const
{

    // chk00010            write a checkpoint file with this root directory
    // chk00010/Header     this contains information you need to save (e.g., finest_level, t_new, etc.) and also
    //                     the BoxArrays at each level
    // chk00010/Level_0/
    // chk00010/Level_1/
    // etc.                these subdirectories will hold the MultiFab data at each level of refinement

    // checkpoint file name, e.g., chk00010
    const std::string& checkpointname = amrex::Concatenate(chk_file,istep[0],7);

    amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

    const int nlevels = finest_level+1;

    // ---- prebuild a hierarchy of directories
    // ---- dirName is built first.  if dirName exists, it is renamed.  then build
    // ---- dirName/subDirPrefix_0 .. dirName/subDirPrefix_nlevels-1
    // ---- if callBarrier is true, call ParallelDescriptor::Barrier()
    // ---- after all directories are built
    // ---- ParallelDescriptor::IOProcessor() creates the directories
    amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);

    // write Header file
   if (ParallelDescriptor::IOProcessor()) {

       std::string HeaderFileName(checkpointname + "/Header");
       VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
       std::ofstream HeaderFile;
       HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
       HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                               std::ofstream::trunc |
                                               std::ofstream::binary);
       if( ! HeaderFile.good()) {
           amrex::FileOpenFailed(HeaderFileName);
       }

       HeaderFile.precision(17);

       // write out title line
       HeaderFile << "Checkpoint file for AmrCoreAdv\n";

       // write out finest_level
       HeaderFile << finest_level << "\n";

       // write the number of components
       HeaderFile << Idx::NumScalars << "\n";

       // write the number of ghost cells
       HeaderFile << NUM_GHOST_CELLS << "\n";

       // write out array of istep
       for (int i = 0; i < istep.size(); ++i) {
           HeaderFile << istep[i] << " ";
       }
       HeaderFile << "\n";

       // write out array of dt
       for (int i = 0; i < dt.size(); ++i) {
           HeaderFile << dt[i] << " ";
       }
       HeaderFile << "\n";

       // write out array of t_new
       for (int i = 0; i < t_new.size(); ++i) {
           HeaderFile << t_new[i] << " ";
       }
       HeaderFile << "\n";

       // write the BoxArray at each level
       for (int lev = 0; lev <= finest_level; ++lev) {
           boxArray(lev).writeOn(HeaderFile);
           HeaderFile << '\n';
       }
   }

   // write the MultiFab data to, e.g., chk00010/Level_0/
   for (int lev = 0; lev <= finest_level; ++lev) {
       VisMF::Write(grid_new[lev],
                    amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "Cell"));
   }

}


void
AmrCoreAdv::ReadCheckpointFile ()
{

    amrex::Print() << "Restart from checkpoint " << restart_chkfile << "\n";

    // Header
    std::string File(restart_chkfile + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    int chk_ncomp, chk_nghost;

    // read in title line
    std::getline(is, line);

    // read in finest_level
    is >> finest_level;
    GotoNextLine(is);

    // read in number of components & assert they are the same as here
    is >> chk_ncomp;
    GotoNextLine(is);
    AMREX_ASSERT(chk_ncomp == Idx::NumScalars);

    // read in number of ghost cells & assert they are the same as here
    is >> chk_nghost;
    GotoNextLine(is);
    AMREX_ASSERT(chk_nghost == NUM_GHOST_CELLS);

    // read in array of istep
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            istep[i++] = std::stoi(word);
        }
    }

    // read in array of dt
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            dt[i++] = std::stod(word);
        }
    }

    // read in array of t_new
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            t_new[i++] = std::stod(word);
        }
    }

    for (int lev = 0; lev <= finest_level; ++lev) {

        // read in level 'lev' BoxArray from Header
        BoxArray ba;
        ba.readFrom(is);
        GotoNextLine(is);

        // create a distribution mapping
        DistributionMapping dm { ba, ParallelDescriptor::NProcs() };

        // set BoxArray grids and DistributionMapping dmap in AMReX_AmrMesh.H class
        SetBoxArray(lev, ba);
        SetDistributionMap(lev, dm);

        // build MultiFab data
        int ncomp = Idx::NumScalars;
        int nghost = NUM_GHOST_CELLS;

        grid_old[lev].define(grids[lev], dmap[lev], ncomp, nghost);
        grid_new[lev].define(grids[lev], dmap[lev], ncomp, nghost);

        // also create the time integrator for this level
        integrator[lev] = std::make_unique<TimeIntegrator<MultiFab> >(grid_old[lev]);
    }

    // read in the MultiFab data
    for (int lev = 0; lev <= finest_level; ++lev) {
        VisMF::Read(grid_new[lev],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Cell"));
    }

}

// Read the file in restart_chkfile and use it as an initial condition for
// the current simulation. Supports a different number of components and
// ghost cells.
void
AmrCoreAdv::InitializeFromFile ()
{
    amrex::Print() << "Initializing from file " << restart_chkfile << "\n";

    const Real time = 0.0;

    // Header
    std::string File(restart_chkfile + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;
    int chk_ncomp, chk_nghost;

    // read in title line
    std::getline(is, line);

    // read in finest_level
    is >> finest_level;
    GotoNextLine(is);

    // read in number of components and ghost cells
    // we don't assert those are the same as in the current solver
    // but we need this to build the temporary MultiFab for reading
    // the data we use to initialize our grid.

    // read in number of components
    is >> chk_ncomp;
    GotoNextLine(is);

    // read in number of ghost cells
    is >> chk_nghost;
    GotoNextLine(is);

    // read in array of istep & ignore
    std::getline(is, line);

    // read in array of dt & ignore
    std::getline(is, line);

    // read in array of t_new & ignore
    std::getline(is, line);

    for (int lev = 0; lev <= max_level; ++lev) {

        // read in level 'lev' BoxArray from Header
        BoxArray ba;
        ba.readFrom(is);
        GotoNextLine(is);

        // create a distribution mapping
        DistributionMapping dm { ba, ParallelDescriptor::NProcs() };

        // set BoxArray grids and DistributionMapping dmap in AMReX_AmrMesh.H class
        SetBoxArray(lev, ba);
        SetDistributionMap(lev, dm);

        // build MultiFab data
        int ncomp = Idx::NumScalars;
        int nghost = NUM_GHOST_CELLS;
        grid_old[lev].define(grids[lev], dmap[lev], ncomp, nghost);
        grid_new[lev].define(grids[lev], dmap[lev], ncomp, nghost);

    }

    // read in the MultiFab data & initialize from it
    for (int lev = 0; lev <= max_level; ++lev) {
        MultiFab initial_data_lev(grids[lev], dmap[lev], chk_ncomp, chk_nghost);
        VisMF::Read(initial_data_lev,
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Cell"));

        // This fills grid_new[lev] from initial_data_lev
        InitializeLevelFromData(lev, initial_data_lev);

        //! Calling RemakeLevel makes sure grid_new[lev] ghost cells
        // are properly filled (can make this more efficient by just
        // interpolating to fill boundary cells in grid_new[lev] instead.
        RemakeLevel(lev, time, grids[lev], dmap[lev]);
    }
}

// Initialize the new-time data at a level from the initial_data MultiFab
void
AmrCoreAdv::InitializeLevelFromData(int lev, const MultiFab& initial_data)
{
    const auto& geom = Geom(lev);
    const auto dx = geom.CellSizeArray();
#ifdef PROBLEM_LOADS_INITDATA
    auto& state_mf = grid_new[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( MFIter mfi(initial_data, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.tilebox();
        auto state_fab = state_mf.array(mfi);
        const auto idata = initial_data.array(mfi);

        // Call a user-supplied function to initialize the state data
        // from the input data file.
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            initialize_from_data(i, j, k, state_fab, idata, dx, geom.data());
        });
    }
#else
    amrex::Error("Custom initialization from external data is not implemented for this problem.");
#endif
}

// utility to skip to next line in Header
void
AmrCoreAdv::GotoNextLine (std::istream& is)
{
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}

void AmrCoreAdv::Trajectory(MultiFab& state_mf, int lev, const amrex::Real time, const Geometry& geom, Parameters Param)
{
    Real dtau = Param.tau/Param.hmc_substeps;
    update_momentum(state_mf, time, geom, Param, dtau/2.0);
    
    
    for(int i = 1; i < Param.hmc_substeps; i++)
    {
        update_gauge(state_mf, lev, time, geom, dtau);
        
        Fermi_BiCG_solve(state_mf, lev, time, geom, Param);
        
        update_momentum(state_mf, time, geom, Param, dtau);
        
    }
    
    update_gauge(state_mf, lev, time, geom, dtau);
    
    Fermi_BiCG_solve(state_mf, lev, time, geom, Param);
    
    update_momentum(state_mf, time, geom, Param, dtau/2.0);
}

Real AmrCoreAdv::Check_Reversed_Trajectory(const MultiFab& state_mf, int lev, const amrex::Real time, const Geometry& geom, Parameters Param)
{
    Real dtau = Param.tau/Param.hmc_substeps;

    int ncomp = state_mf.nComp();
    const BoxArray& ba = state_mf.boxArray();
    const DistributionMapping& dm = state_mf.DistributionMap();
    
    MultiFab* initial_mf = new MultiFab; 
    initial_mf -> define(ba, dm, ncomp, NUM_GHOST_CELLS);
    
    MultiFab::Copy(*initial_mf, state_mf, 0, 0, ncomp, NUM_GHOST_CELLS);
    
    MultiFab* final_mf = new MultiFab; 
    final_mf -> define(ba, dm, ncomp, NUM_GHOST_CELLS);
    
    MultiFab::Copy(*final_mf, state_mf, 0, 0, ncomp, NUM_GHOST_CELLS);
    
    Trajectory(*final_mf, lev, time, geom, Param);
    
    final_mf -> MultiFab::mult(-1, Idx::P_0, 2, NUM_GHOST_CELLS);
    
    Trajectory(*final_mf, lev, time, geom, Param);
    
    MultiFab* diff_mf = new MultiFab;
    diff_mf -> define(ba, dm, 4, NUM_GHOST_CELLS);
    
    MultiFab::LinComb(*diff_mf, 1.0, *final_mf, Idx::U_0_Real, -1.0, *initial_mf, Idx::U_0_Real, 0, 4, NUM_GHOST_CELLS);
    
    Real L_2 = (diff_mf -> MultiFab::norm2(0)) + (diff_mf -> MultiFab::norm2(1)) + (diff_mf -> MultiFab::norm2(2)) + (diff_mf -> MultiFab::norm2(3));
    
    delete initial_mf;
    delete final_mf;
    delete diff_mf;
    
    return L_2;
    
}

void AmrCoreAdv::update_gauge (MultiFab& state_mf, int lev, const amrex::Real time, const Geometry& geom, amrex::Real dtau)
{
    const auto dx = geom.CellSizeArray();
    int ncomp = state_mf.nComp();
    
#ifdef _OPENMP
#pragma omp parallel
#endif
  for ( MFIter mfi(state_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
  {
    const Box& bx = mfi.tilebox();
    const auto ncomp = state_mf.nComp();

    const auto& state_fab = state_mf.array(mfi);
    
    // For each grid, loop over all the valid points
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        state_update_gauge(i, j, k, state_fab, time, dx, geom.data(), dtau);
    });
  }
    
    FillIntermediatePatch(lev, time, state_mf, 0, state_mf.nComp());
    FlipSigns(lev, state_mf, Idx::Phi_0_Real, 4);
    FlipSigns(lev, state_mf, Idx::DDinvPhi_0_Real, 4);
    FlipSigns(lev, state_mf, Idx::chi_0_Real, 4);
    FlipSigns(lev, state_mf, Idx::g3DinvPhi_0_Real, 4);
}

void AmrCoreAdv::update_momentum (MultiFab& state_mf, const amrex::Real time, const amrex::Geometry& geom, Parameters Param, amrex::Real dtau)
{
  const auto dx = geom.CellSizeArray();

#ifdef _OPENMP
#pragma omp parallel
#endif
  for ( MFIter mfi(state_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
  {
    const Box& bx = mfi.tilebox();
    const auto ncomp = state_mf.nComp();

    const auto& state_fab = state_mf.array(mfi);

    // For each grid, loop over all the valid points
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      state_update_momentum(i, j, k, state_fab, time, dx, geom.data(), Param.mass, Param.r, Param.beta, dtau);
    });
  }
}

void
AmrCoreAdv::Fermi_BiCG_solve(MultiFab& S_new, int lev, Real time, const amrex::Geometry& geom_lev, Parameters Param)
{
    amrex::Real alpha, beta, denom, rsq, rsq_new, bsqrt, bnorm;
    
    //amrex::Print() << "Initial Action = " <<  SumActionD(S_new) << std::endl;
    
    //S_new.MultiFab::negate(Idx::DDinvPhi_0_Real, 4, 0);
    
    S_new.MultiFab::setVal(0, Idx::DDinvPhi_0_Real, 4, 0);
    
    //ResetDDinvPhi(S_new);
    
    const BoxArray& ba = S_new.boxArray();
    const DistributionMapping& dm = S_new.DistributionMap();
    
    MultiFab* mfres = new MultiFab;
    mfres -> define(ba,dm,4,NUM_GHOST_CELLS);
    
    MultiFab::Copy(*mfres, S_new, Idx::Phi_0_Real, 0, 4, NUM_GHOST_CELLS);
    
    //FillIntermediatePatch(lev, time, S_new, 0, S_new.nComp());
    
    MultiFab* mfp = new MultiFab;
    mfp -> define(ba,dm,4,NUM_GHOST_CELLS);
    
    MultiFab::Copy(*mfp, *mfres, resIdx::res_0_Real, pIdx::p_0_Real, 4, NUM_GHOST_CELLS);
    
    //FillIntermediatePatch(lev, time, S_new, 0, S_new.nComp());
    
    MultiFab* mfDDp = new MultiFab;
    mfDDp -> define(ba,dm,4,NUM_GHOST_CELLS);
    
    rsq = MultiFab::Dot(*mfres, resIdx::res_0_Real, 4, 0);
    
    //amrex::Print() << rsq << "\n";
    
    for(int i = 0; i < Param.BiCG_Max_Iters; i++)
    {
        
        Set_DDp(*mfDDp, *mfp, S_new, geom_lev, mass_0, wilson_r);
        
        //FillIntermediatePatch(lev, time, S_new, 0, S_new.nComp());
        
        //amrex::Print() << "Made it here \n";
        
        denom = MultiFab::Dot(*mfp, pIdx::p_0_Real, *mfDDp, DDpIdx::DDp_0_Real, 4, 0);
        
        alpha = rsq/denom;
        
        MultiFab::Saxpy(S_new, alpha, *mfp, pIdx::p_0_Real, Idx::DDinvPhi_0_Real, 4, 0); 
        
        //Set_DDinvPhi(S_new, *mfp, geom_lev, alpha);
        
        //FillIntermediatePatch(lev, time, S_new, 0, S_new.nComp());
        
        MultiFab::Saxpy(*mfres, -alpha, *mfDDp, DDpIdx::DDp_0_Real, resIdx::res_0_Real, 4, 0);
        //Set_res(*mfres, *mfDDp, geom_lev, alpha);
        
        //FillIntermediatePatch(lev, time, S_new, 0, S_new.nComp());
        
        rsq_new = MultiFab::Dot(*mfres, resIdx::res_0_Real, 4, 0);
        if (rsq_new < Param.BiCG_Thresh)
            break;
        
        beta = rsq_new/rsq;
        rsq = rsq_new;
        
        MultiFab::LinComb(*mfp, beta, *mfp, pIdx::p_0_Real, 1.0, *mfres, resIdx::res_0_Real, pIdx::p_0_Real, 4, 0); 
        
        //Set_p(*mfp, *mfres, geom_lev, beta);
        
        //FillIntermediatePatch(lev, time, S_new, Idx::chi_0_Real, 16);
        FillIntermediatePatch(lev, time, S_new, 0, S_new.nComp());
        //FlipSigns(lev, S_new, Idx::chi_0_Real, 16);        
        FlipSigns(lev, S_new, Idx::Phi_0_Real, 4);
        FlipSigns(lev, S_new, Idx::DDinvPhi_0_Real, 4);
        FlipSigns(lev, S_new, Idx::chi_0_Real, 4);
        FlipSigns(lev, S_new, Idx::g3DinvPhi_0_Real, 4);
        
        FillIntermediatePatch(lev, time, *mfres, 0, mfres -> nComp());
        FlipSigns(lev, *mfres, resIdx::res_0_Real, 4);
        
        FillIntermediatePatch(lev, time, *mfp, 0, mfp -> nComp());
        FillIntermediatePatch(lev, time, *mfDDp, 0, mfDDp -> nComp());
        FlipSigns(lev, *mfp, pIdx::p_0_Real, 4);
        FlipSigns(lev, *mfDDp, DDpIdx::DDp_0_Real, 4);
            

        
        if (i == (Param.BiCG_Max_Iters - 1) && rsq_new > Param.BiCG_Thresh)
        {
            amrex::Print() << "Failed to converge after " << Param.BiCG_Max_Iters << " steps!" << std::endl;
            amrex::Print() << "Final rsq = " << rsq_new << std::endl;
            
        }
        
        
    }
    
    delete mfres;
    delete mfp;
    delete mfDDp;
    
    Set_g3DinvPhi(S_new, geom_lev, mass_0, wilson_r);
    //FillIntermediatePatch(lev, time, S_new, 0, S_new.nComp());
    
    FillIntermediatePatch(lev, time, S_new, 0, S_new.nComp());
    FlipSigns(lev, S_new, Idx::Phi_0_Real, 4);
    FlipSigns(lev, S_new, Idx::DDinvPhi_0_Real, 4);
    FlipSigns(lev, S_new, Idx::chi_0_Real, 4);
    FlipSigns(lev, S_new, Idx::g3DinvPhi_0_Real, 4);
    
}

void AmrCoreAdv::Set_DDp (MultiFab& DDp_mf, MultiFab& p_mf, MultiFab& state_mf, const Geometry& geom, amrex::Real m_0, amrex::Real r)
{
    const auto dx = geom.CellSizeArray();
    int ncomp = state_mf.nComp();
    
#ifdef _OPENMP
#pragma omp parallel
#endif
  for ( MFIter mfi(state_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
  {
    const Box& bx = mfi.tilebox();
    const auto ncomp = state_mf.nComp();
    
      
    const auto& DDp_fab = DDp_mf.array(mfi);
    const auto& p_fab = p_mf.array(mfi);  
    const auto& state_fab = state_mf.array(mfi);
    
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        state_set_DDp(i, j, k, DDp_fab, p_fab, state_fab, dx, geom.data(), m_0, r);
    });
      
  }
}

void AmrCoreAdv::Set_g3DinvPhi (MultiFab& state_mf, const Geometry& geom, amrex::Real m_0, amrex::Real r)
{
    const auto dx = geom.CellSizeArray();
    int ncomp = state_mf.nComp();
    
#ifdef _OPENMP
#pragma omp parallel
#endif
  for ( MFIter mfi(state_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
  {
    const Box& bx = mfi.tilebox();
    const auto ncomp = state_mf.nComp();

    const auto& state_fab = state_mf.array(mfi);
    
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        state_set_g3DinvPhi(i, j, k, state_fab, dx, geom.data(), m_0, r);
    });
      
  }
}

void AmrCoreAdv::SetPhiFromChi (MultiFab& state_mf, int lev, amrex::Real time, amrex::Real m_0, amrex::Real r)
{

#ifdef _OPENMP
#pragma omp parallel
#endif

  for ( MFIter mfi(state_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
  {
    const Box& bx = mfi.tilebox();
    const auto ncomp = state_mf.nComp();

    const auto& state_fab = state_mf.array(mfi); 

    // For each grid, loop over all the valid points
      
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      state_set_phi_from_chi(i,j,k,state_fab, m_0, r);
    });
      
  }
    
  FillIntermediatePatch(lev, time, state_mf, 0, state_mf.nComp());
  FlipSigns(lev, state_mf, Idx::Phi_0_Real, 4);
  FlipSigns(lev, state_mf, Idx::DDinvPhi_0_Real, 4);
  FlipSigns(lev, state_mf, Idx::chi_0_Real, 4);
  FlipSigns(lev, state_mf, Idx::g3DinvPhi_0_Real, 4);
   
}



void AmrCoreAdv::Perturb (MultiFab& state_mf, int lev, const amrex::Real time, const Geometry& geom)
{
    const auto dx = geom.CellSizeArray();
    int ncomp = state_mf.nComp();

#ifdef _OPENMP
#pragma omp parallel
#endif
  for ( MFIter mfi(state_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
  {
    const Box& bx = mfi.tilebox();
    const auto ncomp = state_mf.nComp();

    const auto& state_fab = state_mf.array(mfi);

    // For each grid, loop over all the valid points

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        state_Perturbation(i, j, k, state_fab, time, dx, geom.data());
    });
  }
    
    FillIntermediatePatch(lev, time, state_mf, 0, state_mf.nComp());
    FlipSigns(lev, state_mf, Idx::Phi_0_Real, 4);
    FlipSigns(lev, state_mf, Idx::DDinvPhi_0_Real, 4);
    FlipSigns(lev, state_mf, Idx::chi_0_Real, 4);
    FlipSigns(lev, state_mf, Idx::g3DinvPhi_0_Real, 4);
}


Real AmrCoreAdv::Action_Gauge (MultiFab& state_mf, Real beta)
{
    
    ReduceOps<ReduceOpSum> reduce_operations;
    ReduceData<Real> reduce_data(reduce_operations);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif

  for ( MFIter mfi(state_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
  {
    const Box& bx = mfi.tilebox();
    const auto ncomp = state_mf.nComp();

    const auto& state_fab = state_mf.array(mfi); 

    // For each grid, loop over all the valid points
    reduce_operations.eval(bx, reduce_data,
    [=] AMREX_GPU_DEVICE (const int i, const int j, const int k) -> ReduceTuple
    {
        return {sum_action_gauge(i,j,k,state_fab, beta)};
    });
  }
    ReduceTuple reduced_values = reduce_data.value();
    // MPI reduction
    ParallelDescriptor::ReduceRealSum(amrex::get<0>(reduced_values));
    Real action = amrex::get<0>(reduced_values);
    
    return action;
}

Real AmrCoreAdv::Action_Momentum (MultiFab& state_mf)
{
    
    ReduceOps<ReduceOpSum> reduce_operations;
    ReduceData<Real> reduce_data(reduce_operations);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif

  for ( MFIter mfi(state_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
  {
    const Box& bx = mfi.tilebox();
    const auto ncomp = state_mf.nComp();

    const auto& state_fab = state_mf.array(mfi); 

    // For each grid, loop over all the valid points
    reduce_operations.eval(bx, reduce_data,
    [=] AMREX_GPU_DEVICE (const int i, const int j, const int k) -> ReduceTuple
    {
        return {sum_action_mom(i,j,k,state_fab)};
    });
  }
    ReduceTuple reduced_values = reduce_data.value();
    // MPI reduction
    ParallelDescriptor::ReduceRealSum(amrex::get<0>(reduced_values));
    Real action = amrex::get<0>(reduced_values);
    
    return action;
}

Real AmrCoreAdv::Action_Fermi (MultiFab& state_mf)
{
    
    ReduceOps<ReduceOpSum> reduce_operations;
    ReduceData<Real> reduce_data(reduce_operations);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif

  for ( MFIter mfi(state_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
  {
    const Box& bx = mfi.tilebox();
    const auto ncomp = state_mf.nComp();

    const auto& state_fab = state_mf.array(mfi); 

    // For each grid, loop over all the valid points
    reduce_operations.eval(bx, reduce_data,
    [=] AMREX_GPU_DEVICE (const int i, const int j, const int k) -> ReduceTuple
    {
        return {sum_action_D(i,j,k,state_fab)};
    });
  }
    ReduceTuple reduced_values = reduce_data.value();
    // MPI reduction
    ParallelDescriptor::ReduceRealSum(amrex::get<0>(reduced_values));
    Real action = amrex::get<0>(reduced_values);
    
    return action;
}

Real AmrCoreAdv::Action_Fermi_From_Chi (MultiFab& state_mf)
{
    
    ReduceOps<ReduceOpSum> reduce_operations;
    ReduceData<Real> reduce_data(reduce_operations);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif

  for ( MFIter mfi(state_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
  {
    const Box& bx = mfi.tilebox();
    const auto ncomp = state_mf.nComp();

    const auto& state_fab = state_mf.array(mfi); 

    // For each grid, loop over all the valid points
    reduce_operations.eval(bx, reduce_data,
    [=] AMREX_GPU_DEVICE (const int i, const int j, const int k) -> ReduceTuple
    {
        return {sum_action_chi(i,j,k,state_fab)};
    });
  }
    ReduceTuple reduced_values = reduce_data.value();
    // MPI reduction
    ParallelDescriptor::ReduceRealSum(amrex::get<0>(reduced_values));
    Real action = amrex::get<0>(reduced_values);
    
    return action;
}


Real AmrCoreAdv::Measure_Plaq (MultiFab& state_mf)
{
    
    ReduceOps<ReduceOpSum> reduce_operations;
    ReduceData<Real> reduce_data(reduce_operations);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif

  for ( MFIter mfi(state_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
  {
    const Box& bx = mfi.tilebox();
    const auto ncomp = state_mf.nComp();

    const auto& state_fab = state_mf.array(mfi); 

    // For each grid, loop over all the valid points
    reduce_operations.eval(bx, reduce_data,
    [=] AMREX_GPU_DEVICE (const int i, const int j, const int k) -> ReduceTuple
    {
        return {sum_plaq(i,j,k,state_fab)/(numcells[0]*numcells[1])};
    });
  }
    ReduceTuple reduced_values = reduce_data.value();
    // MPI reduction
    ParallelDescriptor::ReduceRealSum(amrex::get<0>(reduced_values));
    Real action = amrex::get<0>(reduced_values);
    
    return action;
}


void AmrCoreAdv::smear_gauge (MultiFab& smeared_mf, MultiFab& state_mf, int lev, const amrex::Real time, const Geometry& geom, Parameters Param)
{
    const auto dx = geom.CellSizeArray();
    int ncomp = state_mf.nComp();
    
    const BoxArray& ba = state_mf.boxArray();
    const DistributionMapping& dm = state_mf.DistributionMap();
    
    MultiFab* smeared_tmp_mf = new MultiFab;
    smeared_tmp_mf -> define(ba,dm,4,NUM_GHOST_CELLS);
    
    MultiFab::Copy(smeared_mf, state_mf, Idx::U_0_Real, 0, 4, NUM_GHOST_CELLS);
    
    MultiFab::Copy(*smeared_tmp_mf, smeared_mf, cIdx::Real_0, cIdx::Real_0, 4, NUM_GHOST_CELLS);
    
    for(int iter = 0; iter < Param.APE_smear_iter; iter++)
    {
#ifdef _OPENMP
#pragma omp parallel
#endif

      for ( MFIter mfi(state_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
      {
        const Box& bx = mfi.tilebox();
        const auto ncomp = state_mf.nComp();
    
        const auto& state_fab = state_mf.array(mfi);
        const auto& smeared_fab = smeared_mf.array(mfi);
        const auto& smeared_tmp_fab = smeared_tmp_mf -> array(mfi);
    
        // For each grid, loop over all the valid points
        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            state_smear_gauge(i, j, k, smeared_fab, smeared_tmp_fab, time, dx, geom.data(), Param.APE_smear_alpha);
        });
      }
      
      FillIntermediatePatch(lev, time, *smeared_tmp_mf, 0, 4);
        
      for ( MFIter mfi(state_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
      {
         const Box& bx = mfi.tilebox();
         const auto ncomp = state_mf.nComp();
    
         const auto& state_fab = state_mf.array(mfi);
         const auto& smeared_fab = smeared_mf.array(mfi);
         const auto& smeared_tmp_fab = smeared_tmp_mf -> array(mfi);
    
          // For each grid, loop over all the valid points
         amrex::ParallelFor(bx,
         [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
             state_project_smeared(i, j, k, smeared_fab, smeared_tmp_fab, time, dx, geom.data());
         });
       }
        
    }
    
    delete smeared_tmp_mf;

    FillIntermediatePatch(lev, time, smeared_mf, 0, 4);
  
}

Real AmrCoreAdv::meas_TopCharge (MultiFab& state_mf, int lev, const amrex::Real time, const Geometry& geom, Parameters Param)
{
    const BoxArray& ba = state_mf.boxArray();
    const DistributionMapping& dm = state_mf.DistributionMap();
    
    MultiFab* smeared_mf = new MultiFab;
    smeared_mf -> define(ba, dm, 4, NUM_GHOST_CELLS);
    
    smear_gauge(*smeared_mf, state_mf, lev, time, geom, Param);
    
    ReduceOps<ReduceOpSum> reduce_operations;
    ReduceData<Real> reduce_data(reduce_operations);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif

  for ( MFIter mfi(state_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
  {
    const Box& bx = mfi.tilebox();
    const auto ncomp = state_mf.nComp();

    const auto& state_fab = state_mf.array(mfi);
    const auto& smeared_fab = smeared_mf -> array(mfi);

    // For each grid, loop over all the valid points
    reduce_operations.eval(bx, reduce_data,
    [=] AMREX_GPU_DEVICE (const int i, const int j, const int k) -> ReduceTuple
    {
        return {state_TopCharge(i,j,k,smeared_fab)};
    });
  }
    ReduceTuple reduced_values = reduce_data.value();
    // MPI reduction
    ParallelDescriptor::ReduceRealSum(amrex::get<0>(reduced_values));
    Real TopCharge = amrex::get<0>(reduced_values);
    
    delete smeared_mf;
    
    return TopCharge/(2.0*M_PI);
}



void AmrCoreAdv::fill_state_diagnostics (MultiFab& diag_mf, const MultiFab& state_mf, const Real time_lev, const amrex::Geometry& geom) const
{
  const auto dx = geom.CellSizeArray();

#ifdef _OPENMP
#pragma omp parallel
#endif
  for ( MFIter mfi(diag_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
  {
    const Box& bx = mfi.tilebox();
    const auto ncomp = state_mf.nComp();

    const auto& diag_fab = diag_mf.array(mfi);
    const auto& state_fab = state_mf.array(mfi);

    // For each grid, loop over all the valid points
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      state_diagnostics(i, j, k, diag_fab, state_fab, time_lev, dx, geom.data());
    });
  }
}
