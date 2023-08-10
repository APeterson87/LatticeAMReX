#include <AmrCoreAdv.H>
#include <Adam_Interpolater.H>

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
    grid_hold.resize(nlevs_max);
    
    //integrator.resize(nlevs_max);

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
    if (interpolation_type == InterpType::PCInterp)
        mapper = &pc_interp;
    else if (interpolation_type == InterpType::NodeBilinear)
        mapper = &node_bilinear_interp;
    else if (interpolation_type == InterpType::CellConservativeLinear)
        mapper = &cell_cons_interp;
    else if (interpolation_type == InterpType::CellBilinear)
        mapper = &cell_bilinear_interp;
    else if (interpolation_type == InterpType::CellConservativeQuartic)
        mapper = &quartic_interp;
    else {
        amrex::Error("Unsupported interpolation type");
    }
    
    
   mapper = &adam_interp;

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
                     .beta_lev = coupling_beta_vec,
                     .sigma_lev = mom_sigma_vec,
                     .Temp_lev = temp_vec,
                     .hmc_tau_lev = hmc_tau_vec,
                     
                     .Nx = numcells[0],
                     .Ny = numcells[1],
                     .hmc_substeps = num_hmc_substeps,
                     .therm_steps = num_therm_steps,
                     .starting_meas_step = start_meas_step,
                     .tau = hmc_tau,
                     .use_dynamical_fermions = use_dynamical_fermi,
                     .BiCG_Thresh = BiCGThresh,
                     .BiCG_Max_Iters = BiCG_Max_Iter,
                     .stabilized = Use_BiCG_Stab,
                     .APE_smear_iter = APE_iter, 
                     .APE_smear_alpha = APE_alpha,
                     .measSigmaLoops_Interval = measSigma_Interval,
                     .checkrevtraj_Interval = Check_revTraj_Int};

    amrex::Print() << std::endl << "***************** PARAMETERS *****************" << std::endl;

    amrex::Print() << "\nbeta = " << Param.beta << std::endl;
    amrex::Print() << "Dirac Mass = " << Param.mass << std::endl;
    amrex::Print() << "Wilson Factor = " << Param.r << std::endl;
    amrex::Print() << "Temperature = " << Param.Temp << std::endl;
    amrex::Print() << "\beta at lev 1 " << Param.beta_lev[1] << std::endl;
    
    amrex::Print() << "\nLattice Number x = " << Param.Nx << std::endl;
    amrex::Print() << "Lattice Number y = " << Param.Ny << std::endl;
    
    amrex::Print() << "\nHMC substeps = " << Param.hmc_substeps << std::endl;
    amrex::Print() << "HMC tau = " << Param.tau << std::endl;
    amrex::Print() << "HMC dtau = " << Param.tau/Param.hmc_substeps << std::endl;
    amrex::Print() << "Thermalization steps = " << Param.therm_steps << std::endl;
    amrex::Print() << "Starting measurements at step " << Param.starting_meas_step << std::endl; 
    
    if(Param.use_dynamical_fermions)
        amrex::Print() << "\nUsing dynamical fermions" << std::endl;
    
    amrex::Print() << "\nBiconjugate tolerance  = " << Param.BiCG_Thresh << std::endl;
    amrex::Print() << "Max Biconjugate steps = " << Param.BiCG_Max_Iters << std::endl;    
    
    amrex::Print() << std:: endl << "**********************************************" << std::endl << std::endl;
    
    std::this_thread::sleep_for(std::chrono::seconds(2));
    
    std::ofstream ActionLev0("ActionLev0.dat");
    std::ofstream ActionLev1("ActionLev1.dat");
    
    std::ofstream Sigma("Sigma.dat");
    std::ofstream WilsonLoop("WilsonLoop.dat");
    std::ofstream TopCharge("TopologicalCharges.dat");
    std::ofstream PiCorr("PionCorr.dat");
    
    amrex::Vector<double> pion_corr_total(2*Param.Ny/2+1, 0.0); 
    amrex::Vector<double> sigma_total(2*Param.Nx/2-1, 0.0);
    amrex::Vector<double> wilsonLoop_total(2*Param.Nx/2-1, 0.0);
    
    int histL = 101;
    std::vector<int> histQ(histL);
    for(int i = 0; i < histL; i++) histQ[i] = 0;
    
    int WLcount = 0;
    int num_accepted = 0;
    int num_checked = 0;
    
    for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step)
    {
        
        bool grid_fractured = (((int)std::round(cur_time)%100 <= 2) && (cur_time > 500));
        bool enable_smearing = (((int)std::round(cur_time)%100 <= 4) && (cur_time > 500));
        bool thermalizing = (((int)std::round(cur_time)%100 <= 20) && (cur_time > 500));
            
        amrex::Print() << "\nCoarse STEP " << step+1 << " starts ..." << std::endl;
        
        ComputeDt();  //This is needed to keep track of the current time on each level.
        cur_time += dt[0];

        int lev = 0;
        int iteration = 1;
        
        timeStepRegrid(lev, cur_time, iteration, Param);
        /*
        if(enable_smearing)
        {
            //amrex::Print() << "MADE IT HERE" << std::endl;

                MultiFab U_mf(grid_new[1].boxArray(), grid_new[1].DistributionMap(), 4, NUM_GHOST_CELLS);
                
                //FillPatch(lev, cur_time, grid_new[1], 0, grid_new[1].nComp());
                MultiFab::Copy(U_mf, grid_new[1], Idx::U_0_Real, cIdx::Real_0, 4, NUM_GHOST_CELLS);
                
                MultiFab SmearedU_mf(grid_new[1].boxArray(), grid_new[1].DistributionMap(), 4, NUM_GHOST_CELLS);
                MultiFab::Copy(SmearedU_mf, U_mf, cIdx::Real_0, Idx::U_0_Real, 4, NUM_GHOST_CELLS);
                
                smear_gauge(SmearedU_mf, U_mf, lev, cur_time, Param);
                
                MultiFab::Copy(grid_new[1], SmearedU_mf, cIdx::Real_0, Idx::U_0_Real, 4, NUM_GHOST_CELLS);
                
        }
        */
        //StateAction(grid_new, cur_time, geom, Param);
        
        for (int level = 0; level <= finest_level; ++level) {
            //Copy U and P into grid_hold
            MultiFab::Copy(grid_hold[level], grid_new[level], Idx::U_0_Real, Idx::U_0_Real, 6, grid_new[level].nGrow());    
        }

	// Fill `grid_new` with random values.
        
        StatePerturb(grid_new, cur_time, geom, Param, enable_smearing);

        amrex::Print() << std::endl << "*************** OLD STATE DATA ***************" << std::endl;
        Real HTotalcurrentLev = Total_Action_Lev(grid_new, geom, Param, finest_level); 
        amrex::Print() << "**********************************************" << std::endl << std::endl;

        amrex::Print() << std::endl << "****************** DO MCMC *******************" << std::endl;

        StateTrajectory(grid_new, cur_time, geom, enable_smearing, Param);
        
        /*
        if(istep[0] == 30)
        {
            int level = 1;
            
            const auto geom_lev = geom[level];
        
            const BoxArray& ba = grid_new[level].boxArray();
            
            const DistributionMapping& dm = grid_new[level].DistributionMap();
            
            //FillIntermediatePatch(level, cur_time, grid_aux[level], 0, grid_aux[level].nComp());
                
            MultiFab x_mf(ba, dm, 4, grid_new[level].nGrow());
            MultiFab::Swap(x_mf, grid_new[level], Idx::propUp_0_Real, cIdx::Real_0, 4, grid_new[level].nGrow()); 
            
            MultiFab b_mf(ba, dm, 4, grid_new[level].nGrow());
            MultiFab::Swap(b_mf, grid_new[level], Idx::Phi_0_Real, cIdx::Real_0, 4, grid_new[level].nGrow());
                
            //FillIntermediatePatch(level, cur_time, grid_new[level], 0, grid_new[level].nComp());
            //MultiFab U_mf(ba, dm, 4, grid_new[level].nGrow());
            //MultiFab::Swap(U_mf, grid_new[level], Idx::U_0_Real, cIdx::Real_0, 4, grid_new[level].nGrow());
            
            BiCG_Solve(x_mf, b_mf, grid_new[level], false, level, cur_time, geom_lev, Param);
                
            MultiFab::Swap(x_mf, grid_new[level], Idx::propUp_0_Real, cIdx::Real_0, 4, grid_new[level].nGrow()); 
            MultiFab::Swap(b_mf, grid_new[level], Idx::Phi_0_Real, cIdx::Real_0, 4, grid_new[level].nGrow());
            //MultiFab::Swap(U_mf, grid_new[level], Idx::U_0_Real, cIdx::Real_0, 4, grid_new[level].nGrow());
            
        }
        */
        
        amrex::Print() << std::endl;
        
        amrex::Print() << "*************** NEW STATE DATA ***************" << std::endl;
        Real HTotalLev = Total_Action_Lev(grid_new, geom, Param, finest_level);
        amrex::Print() << "**********************************************" << std::endl;
        
        Real r_loc = std::rand()/(static_cast <float> (RAND_MAX));        
        amrex::Print() << "Metropolis random number = " << r_loc << std::endl;
        amrex::Print() << "Exp(-deltaH/T) = " << std::exp(-(HTotalLev-HTotalcurrentLev)/Temp_T) << std::endl;
        amrex::Print() << "deltaHLev = " << HTotalLev - HTotalcurrentLev << std::endl;
        
        //bool is_fractured = (cur_time >= time_lower && cur_time <= time_upper);//Check if there are uncovered cells on the lowest level?
        
        if(/*!thermalizing &&*/ step >= Param.starting_meas_step)
            num_checked++;
        
        if(r_loc > std::exp(-(HTotalLev-HTotalcurrentLev)/Temp_T) && step >= Param.therm_steps /*&& !thermalizing*/)
        {
            for (int level = 0; level <= finest_level; ++level) {
                 MultiFab::Copy(grid_new[level], grid_hold[level], Idx::U_0_Real, Idx::U_0_Real, 6, grid_new[level].nGrow()); 
            }
            amrex::Print() << "NEW STATE REJECTED " << std::endl;
        }
        else
        {
            if(step >= Param.starting_meas_step /*&& !thermalizing*/) 
            {
                num_accepted++;
            }
            
            amrex::Print() << "NEW STATE ACCEPTED " << std::endl;
        }
        amrex::Print() << "**********************************************" << std::endl << std::endl;
        
        
        Real originalGAction = Total_Gauge_Action_Lev(grid_new, geom, Param, finest_level);
        Real originalMAction = Total_Momentum_Action_Lev(grid_new, geom, Param, finest_level);
        Real originalTotAction = originalGAction + originalMAction;
        
        //if (step >= Param.starting_meas_step && !enable_smearing)
            ActionLev1 << step << ' ' << originalGAction << ' ' << originalMAction << ' ' << originalTotAction << std::endl;
        
        if(step >= Param.starting_meas_step)
        {
            amrex::Print() << "**************** MEASUREMENTS ****************" << std::endl;

            
            const BoxArray& ba_lev = grid_new[1].boxArray();
            
            const DistributionMapping& dm_lev = grid_new[1].DistributionMap();
            
            MultiFab U_mf(ba_lev, dm_lev, 4, grid_new[1].nGrow());
            
            
            
            if (step%measSigma_Interval == 0 && !enable_smearing)
            {
                WLcount++;
                 
                /*
                amrex::Print() << "Measuring String Tension" << std::endl;
            
                amrex::Vector<double> sigma(2*Param.Nx/2-1, 0.0);
                meas_SigmaLoops(U_mf, 1, cur_time, geom[1], Param, sigma, Sigma);
                
                std::transform(sigma_total.begin(), sigma_total.end(), sigma.begin(), sigma_total.begin(), std::plus<Real>());
                for (Real e: sigma_total) amrex::Print() << std::fixed << std::setprecision(9) << e/WLcount << ' ';
                
                amrex::Print() << std::endl;
                */
                /*
                amrex::Print() << "Measuring Wilson Loops" << std::endl;
                
                amrex::Vector<double> wilsonLoop(2*Param.Nx/2-1, 0.0);
                MultiFab::Swap(U_mf, grid_new[1], Idx::U_0_Real, cIdx::Real_0, 4, grid_new[1].nGrow());
                meas_WilsonLoops(U_mf, 1, cur_time, step, geom[1], Param, wilsonLoop, WilsonLoop);
                MultiFab::Swap(U_mf, grid_new[1], Idx::U_0_Real, cIdx::Real_0, 4, grid_new[1].nGrow());
                


                std::transform(wilsonLoop_total.begin(), wilsonLoop_total.end(), wilsonLoop.begin(), wilsonLoop_total.begin(), std::plus<Real>());
                for (Real e: wilsonLoop_total) amrex::Print() << std::fixed << std::setprecision(9) << e/WLcount << ' ';
                

                amrex::Print() << std::endl;
                
                
                amrex::Print() << "Measuring Pion Correlation" << std::endl;
            
                amrex::Vector<double> pi_corr(2*Param.Ny/2+1, 0.0);
                PionCorrelation(grid_new[1], 1, cur_time, step, geom[1], Param, pi_corr, PiCorr, grid_fractured);
                
                std::transform(pion_corr_total.begin(), pion_corr_total.end(), pi_corr.begin(), pion_corr_total.begin(), std::plus<Real>());
                for (Real e: pion_corr_total) amrex::Print() << std::fixed << std::setprecision(9) << e/WLcount << ' ';
                
                amrex::Print() << std::endl << std::endl;
                */
                
                
            }
            
            
            
            
        }
        
        int TC_lev = 1;
        const BoxArray& ba_lev = grid_new[TC_lev].boxArray();    
        const DistributionMapping& dm_lev = grid_new[TC_lev].DistributionMap();
        MultiFab U_mf(ba_lev, dm_lev, 4, grid_new[TC_lev].nGrow());  //Copy gauge part into its own multifab U_mf.
        MultiFab::Copy(U_mf, grid_new[TC_lev], Idx::U_0_Real, cIdx::Real_0, 4, grid_new[TC_lev].nGrow());

        // Measure the TC integer

        Real InstantonNumberLev = meas_TopCharge(U_mf, TC_lev, cur_time, geom[TC_lev], Param);
        Real InstantonNumber = meas_Total_TopCharge_Lev(grid_new, cur_time, geom, Param, 0);
        
        if(step >= Param.starting_meas_step && !enable_smearing)
        {
            TopCharge << step << " " << (int)std::round(InstantonNumberLev) << " " << (int)std::round(InstantonNumber) << std::endl;
            amrex::Print() << "Top Charge " << (int)std::round(InstantonNumberLev) << std::endl;
        }
        


        //Real InstantonNumber = meas_TopCharge(U_mf, 1, cur_time, geom[1], Param);
        //TopCharge << step << " " << (int)std::round(InstantonNumber) << std::endl;
	int top_index = (int)std::round(InstantonNumberLev) + (histL-1)/2;
	histQ[top_index]++;

	//Update topological charge histogram
	std::string name = "top_histogram.dat";
	char fname[256];
	snprintf(fname, 100, "%s", name.c_str());
	FILE *fp = fopen(fname, "w");
	for(int i=0; i<histL; i++) fprintf(fp, "%d %d\n", i - (histL-1)/2, histQ[i]);
	fclose(fp);
	
        // Collect TC density data
        // Smear here
        update_State_topChargeDensity(grid_new, cur_time, geom, Param);
        
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
    
    amrex::Print() << "Total Accepted = " << num_accepted << std::endl;
    amrex::Print() << "Total Checked = " << num_checked << std::endl;
    amrex::Print() << "Final Acceptance Rate = " << static_cast<float>(num_accepted)/static_cast<float>(num_checked) << std::endl;
    
    ActionLev0.close();
    ActionLev1.close();
    
    Sigma.close();
    WilsonLoop.close();
    TopCharge.close();
    PiCorr.close();

    if (plot_int > 0 && istep[0] > last_plot_file_step) {
        WritePlotFile();
    }

    if (chk_int > 0 && istep[0] > last_chk_file_step) {
        WriteCheckpointFile();
    }
    
    if (diag_int > 0 && istep[0] > last_diag_file_step) {
        ComputeAndWriteDiagnosticFile();
    }
}

// initializes multilevel data
void
AmrCoreAdv::InitData ()
{
    if (restart_chkfile == "") {
        // start simulation from the beginning
        const Real time = 0.0;
        InitFromScratch(time);
        AverageDownNodal();

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
    grid_hold[lev].define(ba, dm, ncomp, nghost);
    

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    FillCoarsePatch(lev, time, grid_new[lev], 0, ncomp);
            
    //Copy U and P into grid_hold
    MultiFab::Copy(grid_hold[lev], grid_new[lev], Idx::U_0_Real, Idx::U_0_Real, 6, grid_new[lev].nGrow()); 
    
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
    MultiFab hold_state(ba, dm, ncomp, nghost);
    
    FillPatch(lev, time, new_state, 0, ncomp, 0);
    
    std::swap(new_state, grid_new[lev]);
    std::swap(hold_state, grid_hold[lev]);


    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

}

void
AmrCoreAdv::ResetLevel (int lev, Real time, amrex::Vector<amrex::MultiFab>& hold_grid, bool grid_fractured)
{
    const int ncomp = grid_new[lev].nComp();
    const int nghost = grid_new[lev].nGrow();
    
    auto ba = hold_grid[lev].boxArray();
    auto dm = hold_grid[lev].DistributionMap();
    
    
    MultiFab new_state(ba, dm, ncomp, nghost);
    MultiFab old_state(ba, dm, ncomp, nghost);

    FillPatch(lev, time, new_state, 0, ncomp, grid_fractured);

    std::swap(new_state, grid_new[lev]);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;
    
    MultiFab::Copy(grid_new[lev], hold_grid[lev], Idx::U_0_Real, Idx::U_0_Real, 6, nghost);
}

// Delete level data
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::ClearLevel (int lev)
{
    grid_new[lev].clear();
    grid_hold[lev].clear();

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
    grid_hold[lev].define(ba, dm, ncomp, nghost);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    Real cur_time = t_new[lev];
    MultiFab& state = grid_new[lev];

    const auto& geom = Geom(lev);
    const auto dx = geom.CellSizeArray();

    //state.setVal(1.1e30);
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
        state_init(i, j, k, state_arr, cur_time, lev, geom.data());
      });
    }
    
    //FillPatch(lev, time, state, 0, state.nComp());  //Perhaps find a new place for this.
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
    
    bool grid_fractured = (((int)std::round(time)%100 <= 2) && (time > 500));

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
	      tagarr(i, j, k) = state_is_tagged(i, j, k, lev, state_fab, error_threshold, grid_fractured, time, time_lower, time_upper, x0_lower, x0_upper, x1_lower, x1_upper, geom.data()) ? tagval : clearval;
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
        pp.getarr("coupling_beta_lev", coupling_beta_vec);
        pp.getarr("sigma_lev", mom_sigma_vec);
        pp.getarr("temp_lev", temp_vec);
        pp.getarr("hmc_tau_lev", hmc_tau_vec);

        
        pp.get("Use_Dynamical_Fermions", use_dynamical_fermi);
        pp.get("BiCGThresh", BiCGThresh);
        pp.get("BiCG_Max_Iter", BiCG_Max_Iter);
        pp.get("Use_BiCG_Stab", Use_BiCG_Stab);
        
        pp.get("time_lower", time_lower);
        pp.get("time_upper", time_upper);
        pp.get("x0_lower", x0_lower);
        pp.get("x0_upper", x0_upper);
        pp.get("x1_lower", x1_lower);
        pp.get("x1_upper", x1_upper);

	
        pp.get("APE_alpha", APE_alpha);
        pp.get("APE_iter", APE_iter);
        
        pp.get("measSigma_Interval", measSigma_Interval);
        pp.get("CheckReversedTrajectoryInterval", Check_revTraj_Int);
        

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
/*
void
AmrCoreAdv::AverageDownNodal ()
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
AmrCoreAdv::AverageDownToNodal (int crse_lev)
{
    amrex::average_down(grid_new[crse_lev+1], grid_new[crse_lev],
                        geom[crse_lev+1], geom[crse_lev],
                        0, grid_new[crse_lev].nComp(), refRatio(crse_lev));
}
*/
/*
void
AmrCoreAdv::AdamDown ()
{
    for (int lev = finest_level-1; lev >= 0; --lev)
    {
        adam_down(grid_new[lev+1], grid_new[lev],
                            geom[lev+1], geom[lev],
                            0, grid_new[lev].nComp(), refRatio(lev));
    }
}

// more flexible version of AverageDown() that lets you average down across multiple levels
void
AmrCoreAdv::AdamDownTo (int crse_lev)
{
    adam_down(grid_new[crse_lev+1], grid_new[crse_lev],
                        geom[crse_lev+1], geom[crse_lev],
                        0, grid_new[crse_lev].nComp(), refRatio(crse_lev));
}
*/

void AmrCoreAdv::AverageDownNodal()
{
    for (int lev = finest_level-1; lev >= 0; --lev)
    {
        adam_down_nodal(grid_new[lev+1], grid_new[lev], refRatio(lev));
    }
}

void
AmrCoreAdv::AverageDownToNodal (int crse_lev)
{
    adam_down_nodal(grid_new[crse_lev+1], grid_new[crse_lev], refRatio(crse_lev));
}

template <typename FAB>
void AmrCoreAdv::adam_down_nodal (const FabArray<FAB>& fine, FabArray<FAB>& crse,
                         const IntVect& ratio, int ngcrse, bool mfiter_is_definitely_safe)
{
    //AMREX_ASSERT(fine.is_nodal());
    //AMREX_ASSERT(crse.is_nodal());
    AMREX_ASSERT(crse.nComp() == fine.nComp());

    int ncomp = crse.nComp();
    using value_type = typename FAB::value_type;

    if (mfiter_is_definitely_safe || isMFIterSafe(fine, crse))
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(crse,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(ngcrse);
            Array4<value_type> const& crsearr = crse.array(mfi);
            Array4<value_type const> const& finearr = fine.const_array(mfi);

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                adam_avgdown_nodes(tbx,crsearr,finearr,0,0,ncomp,ratio);
            });
        }
    }
    else
    {
        FabArray<FAB> ctmp(amrex::coarsen(fine.boxArray(),ratio), fine.DistributionMap(),
                           ncomp, ngcrse);
        adam_down_nodal(fine, ctmp, ratio, ngcrse);
        crse.ParallelCopy(ctmp,0,0,ncomp,ngcrse,ngcrse);
    }
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
                fab(i, j, k, n) *= 1;
            }
        }
    });
  }
}

// compute a new multifab by coping in data from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void
AmrCoreAdv::FillPatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp, bool grid_fractured)
{
    if (lev == 0)
    {
        Vector<MultiFab*> smf;
        Vector<Real> stime;
        GetData(lev, time, smf, stime);

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
        MultiFab mf_coarse_temp(grid_hold[lev-1].boxArray(), grid_hold[lev-1].DistributionMap(),
                                grid_hold[lev-1].nComp(), NUM_GHOST_CELLS);

        // fill mf_coarse_temp using MC Equation 39 at "time"
        Real timestep_fraction = (time - t_old[lev-1])/dt[lev-1];
        //integrator[lev-1]->time_interpolate(grid_new[lev-1], grid_old[lev-1], timestep_fraction, mf_coarse_temp);

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
    //mf.FillBoundary(geom[lev].periodicity());
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
        data.push_back(&grid_hold[lev]);
        datatime.push_back(t_old[lev]);
    }
    else
    {
        data.push_back(&grid_hold[lev]);
        data.push_back(&grid_new[lev]);
        datatime.push_back(t_old[lev]);
        datatime.push_back(t_new[lev]);
    }
}


// advance a level by dt
// includes a recursive call for finer levels
void
AmrCoreAdv::timeStepRegrid (int lev, Real time, int iteration, Parameters Param)
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
	      is_regridded = true;
                // regrid could add newly refine levels (if finest_level < max_level)
                // so we save the previous finest level index
                int old_finest = finest_level;
                regrid(lev, time);
                
                //amrex::Print() << "Finest level " << finest_level << std::endl;

                // mark that we have regridded this level already
                for (int k = lev; k <= finest_level; ++k) {
                    last_regrid_step[k] = istep[k];
                }

                // if there are newly created levels, set the time step
                /*for (int k = old_finest+1; k <= finest_level; ++k) {
                    dt[k] = dt[k-1] / MaxRefRatio(k-1);
                }*/
            } else {
	      is_regridded = false;
	    }
        }
    }
    
    ++istep[lev];
    if (lev < finest_level)
    {
        // recursive call for next-finer level
        for (int i = 1; i <= nsubsteps[lev+1]; ++i)
        {
            timeStepRegrid(lev+1, time+(i-1)*dt[lev+1], i, Param);
        }
    }
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

        grid_hold[lev].define(grids[lev], dmap[lev], ncomp, nghost);
        grid_new[lev].define(grids[lev], dmap[lev], ncomp, nghost);

        // also create the time integrator for this level
        //integrator[lev] = std::make_unique<TimeIntegrator<MultiFab> >(grid_old[lev]);
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
        grid_hold[lev].define(grids[lev], dmap[lev], ncomp, nghost);
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

void AmrCoreAdv::Perturb (MultiFab& state_mf, int lev, Real time, const Geometry& geom, Parameters Param)
{
    const BoxArray& ba_lev = state_mf.boxArray();
    const DistributionMapping& dm_lev = state_mf.DistributionMap();
    
    amrex::Real sigma_lev = Param.sigma_lev[lev];
    
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
      state_Perturbation(i, j, k, state_fab, sigma_lev, time, geom.data());
    });
  }  
}

void AmrCoreAdv::StatePerturb (amrex::Vector<amrex::MultiFab>& state, Real time, const amrex::Vector<amrex::Geometry>& geom, Parameters Param, bool grid_fractured)
{
    for(int lev = finest_level; lev >= 0; lev--) {
      Perturb(state[lev], lev, time, geom[lev], Param);
    }
    //if(grid_fractured)
        AverageDownNodal();
}

void AmrCoreAdv::StateAction (amrex::Vector<amrex::MultiFab>& state, Real time, const amrex::Vector<amrex::Geometry>& geom, Parameters Param)
{
    for(int lev = finest_level; lev >= 0; lev--)
    {
        update_action(state[lev], lev, time, geom[lev], Param);
    }
}

/* Evolve the `state` MF with one Leapfrog trajectory */
void
AmrCoreAdv::StateTrajectory(amrex::Vector<amrex::MultiFab>& state, const amrex::Real time, const amrex::Vector<amrex::Geometry>& geom, bool enable_smearing, Parameters Param) {
  
  int expansion_factor = 1;
  if (enable_smearing) 
      expansion_factor = 1;
  
  for(int lev = finest_level; lev >= 0; lev--) {
    Real dtau = Param.hmc_tau_lev[lev]/(Param.hmc_substeps*expansion_factor);
    //Preping for the next step, FillPatches at each level. 
    FillPatch(lev, time, state[lev], 0, state[lev].nComp(), enable_smearing);
    update_momentum(state[lev], lev, time, geom[lev], Param, dtau/2.0);       
  }
  //if(enable_smearing)
      AverageDownNodal();  //Now do all the AveragingDown after time stepping is complete...?
      
  for(int i = 1; i < Param.hmc_substeps*expansion_factor; i++) {
    for(int lev = finest_level; lev >= 0; lev--) { //Update finer levels first.
      
      Real dtau = Param.hmc_tau_lev[lev]/(Param.hmc_substeps*expansion_factor);
      
      //Fill Ghost Cells right before I need them.  Less confusing/ambiguous.
      FillPatch(lev, time, state[lev], 0, state[lev].nComp(), enable_smearing);
      update_gauge(state[lev], lev, time, geom[lev], Param, dtau); 
      
      FillPatch(lev, time, state[lev], 0, state[lev].nComp(), enable_smearing);
      update_momentum(state[lev], lev, time, geom[lev], Param, dtau);
    }
    //if(enable_smearing)
        AverageDownNodal();  //Now do all the AveragingDown after time stepping is complete...?
  }
    
  /////////////////////////////////////////////////////////////////////////////////////////
  for(int lev = finest_level; lev >= 0; lev--) {
    
    Real dtau = Param.hmc_tau_lev[lev]/(Param.hmc_substeps*expansion_factor);   
    
    FillPatch(lev, time, state[lev], 0, state[lev].nComp(), enable_smearing);  //Preping for the next step, FillPatches at each level.
    update_gauge(state[lev], lev, time, geom[lev], Param, dtau);
    
    FillPatch(lev, time, state[lev], 0, state[lev].nComp(), enable_smearing);
    update_momentum(state[lev], lev, time, geom[lev], Param, dtau/2.0);
    
  }
  //if(enable_smearing)
      AverageDownNodal();  //Now do all the AveragingDown after time stepping is complete...?
}

void AmrCoreAdv::update_gauge (MultiFab& state_mf, int lev, const amrex::Real time, const Geometry& geom, Parameters Param, amrex::Real dtau)
{
    //auto& nodal_mask = grid_msk[lev];
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
    //const auto& mask_arr = nodal_mask -> array(mfi);
    
    // For each grid, loop over all the valid points
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
            state_update_gauge(i, j, k, state_fab, time, dx, geom.data(), dtau);
            state_normalize_gauge(i, j, k, state_fab, geom.data(), Param.beta_lev[0]);
    });
  }
   
  
}

void AmrCoreAdv::update_momentum (MultiFab& state_mf, int lev, const amrex::Real time, const amrex::Geometry& geom, Parameters Param, amrex::Real dtau)
{
    
  const auto dx = geom.CellSizeArray();
    
  const BoxArray& ba_lev = state_mf.boxArray();
  const DistributionMapping& dm_lev = state_mf.DistributionMap(); 

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
          state_update_momentum(i, j, k, state_fab, time, dx, geom.data(), Param.mass, Param.r, Param.beta_lev[lev], Param.use_dynamical_fermions, dtau);
    });
  }

    
}

void AmrCoreAdv::update_action (MultiFab& state_mf, int lev, const amrex::Real time, const Geometry& geom, Parameters Param)
{
    //auto& nodal_mask = grid_msk[lev];
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
        state_normalize_gauge(i, j, k, state_fab, geom.data(), Param.beta_lev[0]);
    });
  }
   
  
}

void AmrCoreAdv::update_topChargeDensity (MultiFab& state_mf, int lev, const amrex::Real time, const Geometry& geom, Parameters Param)
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
    //const auto& mask_arr = nodal_mask -> array(mfi);
    
    // For each grid, loop over all the valid points
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        state_topChargeDensity(i, j, k, state_fab, geom.data());
    });
  }
   
}

void AmrCoreAdv::update_State_topChargeDensity (amrex::Vector<amrex::MultiFab>& state, Real time, const amrex::Vector<amrex::Geometry>& geom, Parameters Param)
{
    for(int lev = finest_level; lev >= 0; lev--)
    {
        update_topChargeDensity (state[lev], lev, time, geom[lev], Param);
    }
}


Real AmrCoreAdv::Total_Action(MultiFab& state_mf, int lev, const amrex::Geometry& geom, Parameters Param)
{
    Real GaugeAction = Action_Gauge(state_mf, lev, geom, Param.beta_lev[lev]);
    Real MomAction = Action_Momentum(state_mf, lev, geom);
    Real TotalAction = GaugeAction + MomAction;
    
    amrex::Print() << "Gauge action = " << GaugeAction << std::endl;
    amrex::Print() << "Momentum action = " << MomAction << std::endl;
    amrex::Print() << "Total action = " << TotalAction << std::endl;
    
    return TotalAction;
}

Real AmrCoreAdv::Total_Action_Lev(amrex::Vector<amrex::MultiFab>& state, const amrex::Vector<amrex::Geometry>& geom, Parameters Param, int high_lev)
{
    Real GaugeAction = 0;
    Real MomAction = 0;
    for(int i = 0; i <= high_lev; i++)
    {
        GaugeAction += Action_Gauge_Lev(state, i, geom[i], Param)/Param.Temp_lev[i];
        MomAction += Action_Momentum_Lev(state, i, geom[i])/Param.Temp_lev[i];
        
    }
    Real TotalAction = GaugeAction + MomAction;
    
    amrex::Print() << "Gauge action = " << GaugeAction << std::endl;
    amrex::Print() << "Momentum action = " << MomAction << std::endl;
    amrex::Print() << "Total action = " << TotalAction << std::endl;
    
    return TotalAction;
}

Real AmrCoreAdv::Total_Gauge_Action_Lev(amrex::Vector<amrex::MultiFab>& state, const amrex::Vector<amrex::Geometry>& geom, Parameters Param, int high_lev)
{
    Real GaugeAction = 0;
    for(int i = 0; i <= high_lev; i++)
    {
        GaugeAction += Action_Gauge_Lev(state, i, geom[i], Param)/Param.Temp_lev[i];   
    }

    return GaugeAction;
}

Real AmrCoreAdv::Total_Momentum_Action_Lev(amrex::Vector<amrex::MultiFab>& state, const amrex::Vector<amrex::Geometry>& geom, Parameters Param, int high_lev)
{
    Real MomAction = 0;
    for(int i = 0; i <= high_lev; i++)
    {
        MomAction += Action_Momentum_Lev(state, i, geom[i])/Param.Temp_lev[i];  
    }

    return MomAction;
}

Real AmrCoreAdv::Action_Gauge (MultiFab& state_mf, int lev, const amrex::Geometry& geom, Real beta)
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
        return {sum_action_gauge(i,j,k,state_fab, geom.data(), beta)};
    });
  }
    ReduceTuple reduced_values = reduce_data.value();
    // MPI reduction
    ParallelDescriptor::ReduceRealSum(amrex::get<0>(reduced_values));
    Real action = amrex::get<0>(reduced_values);
    
    return action;
}

Real AmrCoreAdv::Action_Momentum (MultiFab& state_mf, int lev, const amrex::Geometry& geom)
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
        return {sum_action_mom(i,j,k,state_fab, geom.data())};
    });
  }
    ReduceTuple reduced_values = reduce_data.value();
    // MPI reduction
    ParallelDescriptor::ReduceRealSum(amrex::get<0>(reduced_values));
    Real action = amrex::get<0>(reduced_values);
    
    return action;
}

Real AmrCoreAdv::Action_Gauge_Lev (amrex::Vector<amrex::MultiFab>& state, int lev, const amrex::Geometry& geom, Parameters Param)
{
    if(lev != finest_level)
    {
        auto& state_mf = state[lev];


        BoxArray fba = state[lev+1].boxArray();
        const iMultiFab& mask = makeFineMask(state[lev], fba, IntVect(2));  // the last argument is the refinement ratio

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
        const auto& fine_mask_arr = mask.array(mfi);

        // For each grid, loop over all the valid points
        reduce_operations.eval(bx, reduce_data,
        [=] AMREX_GPU_DEVICE (const int i, const int j, const int k) -> ReduceTuple
        {
            return (!fine_mask_arr(i, j, k)*sum_action_gauge(i,j,k,state_fab, geom.data(), Param.beta_lev[lev]));
        });
      }
        ReduceTuple reduced_values = reduce_data.value();
        // MPI reduction
        ParallelDescriptor::ReduceRealSum(amrex::get<0>(reduced_values));
        Real action = amrex::get<0>(reduced_values);
    
    return action;
    }
    else
    {
        auto& state_mf = state[lev];

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
            return {sum_action_gauge(i,j,k,state_fab, geom.data(), Param.beta_lev[lev])};
        });
      }
        ReduceTuple reduced_values = reduce_data.value();
        // MPI reduction
        ParallelDescriptor::ReduceRealSum(amrex::get<0>(reduced_values));
        Real action = amrex::get<0>(reduced_values);
    
    return action;
    }
        
}

Real AmrCoreAdv::Action_Momentum_Lev (amrex::Vector<amrex::MultiFab>& state, int lev, const amrex::Geometry& geom)
{
    if(lev != finest_level)
    {
        auto& state_mf = state[lev];


        BoxArray fba = state[lev+1].boxArray();
        const iMultiFab& mask = makeFineMask(state[lev], fba, IntVect(2));  // the last argument is the refinement ratio

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
        const auto& fine_mask_arr = mask.array(mfi);

        // For each grid, loop over all the valid points
        reduce_operations.eval(bx, reduce_data,
        [=] AMREX_GPU_DEVICE (const int i, const int j, const int k) -> ReduceTuple
        {
            return (!fine_mask_arr(i, j, k)*sum_action_mom(i,j,k, state_fab,geom.data()));
        });
      }
        ReduceTuple reduced_values = reduce_data.value();
        // MPI reduction
        ParallelDescriptor::ReduceRealSum(amrex::get<0>(reduced_values));
        Real action = amrex::get<0>(reduced_values);
    
    return action;
    }
    else
    {
        auto& state_mf = state[lev];

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
            return {sum_action_mom(i,j,k,state_fab,geom.data())};
        });
      }
        ReduceTuple reduced_values = reduce_data.value();
        // MPI reduction
        ParallelDescriptor::ReduceRealSum(amrex::get<0>(reduced_values));
        Real action = amrex::get<0>(reduced_values);
    
    return action;
    }
        
}

void AmrCoreAdv::smear_gauge (MultiFab& smearedU_mf, MultiFab& U_mf, int lev, const amrex::Real time, Parameters Param)
{
    
    //amrex::Print() << "DOING SMEARING" << std::endl;
    
    const Box& domain_bx = geom[lev].Domain();
    
    const BoxArray& ba = U_mf.boxArray();
    const DistributionMapping& dm = U_mf.DistributionMap();
    
    
    MultiFab* smearedU_tmp_mf = new MultiFab;
    smearedU_tmp_mf -> define(ba,dm,4,NUM_GHOST_CELLS);
    
    //Set smearedU_mf to U_mf
    MultiFab::Copy(smearedU_mf, U_mf, cIdx::Real_0, cIdx::Real_0, 4, NUM_GHOST_CELLS);
    
    for(int iter = 0; iter < Param.APE_smear_iter; iter++)
    {
#ifdef _OPENMP
#pragma omp parallel
#endif
 
      for ( MFIter mfi(smearedU_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
      {
        const Box& bx = mfi.tilebox();
    
        const auto& U_fab = U_mf.array(mfi);
        const auto& smeared_fab = smearedU_mf.array(mfi);
        const auto& smeared_tmp_fab = smearedU_tmp_mf -> array(mfi);
    
        // For each grid, loop over all the valid points
        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            //if(i > 6 && i < 10 && j > 6 && j < 10)
                 state_smear_gauge(i, j, k, smeared_fab, smeared_tmp_fab, time, Param.APE_smear_alpha);
        });
      }

      
      for ( MFIter mfi(smearedU_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
      {
         const Box& bx = mfi.tilebox();
    
         const auto& smeared_fab = smearedU_mf.array(mfi);
         const auto& smeared_tmp_fab = smearedU_tmp_mf -> array(mfi);
    
          // For each grid, loop over all the valid points
         amrex::ParallelFor(bx,
         [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
             //if(i > 6 && i < 10 && j > 6 && j < 10)
                 state_project_smeared(i, j, k, smeared_fab, smeared_tmp_fab, time);
         });
       }
       
      

      // Copy smeared_tmp_fab into smeared_tmp_fab
      //MultiFab::Copy(smearedU_mf, U_mf, cIdx::Real_0, cIdx::Real_0, 4, NUM_GHOST_CELLS);
    }   
    
    
    delete smearedU_tmp_mf;
}

Real AmrCoreAdv::meas_TopCharge (MultiFab& U_mf, int lev, const amrex::Real time, const Geometry& geom, Parameters Param)
{
       
    const BoxArray& ba = U_mf.boxArray(); 
    const DistributionMapping& dm = U_mf.DistributionMap();
    
    MultiFab* smearedU_mf = new MultiFab;
    smearedU_mf -> define(ba, dm, 4, NUM_GHOST_CELLS);
    MultiFab::Copy(*smearedU_mf, U_mf, cIdx::Real_0, cIdx::Real_0, 4, NUM_GHOST_CELLS);
    //smear_gauge(*smearedU_mf, U_mf, lev, time, Param);
    
    ReduceOps<ReduceOpSum> reduce_operations;
    ReduceData<Real> reduce_data(reduce_operations);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif

  for ( MFIter mfi(U_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
  {
    const Box& bx = mfi.tilebox();
    Dim3 hi = ubound(bx);

    const auto& smeared_fab = smearedU_mf->array(mfi);

    // For each grid, loop over all the valid points
    reduce_operations.eval(bx, reduce_data,
    [=] AMREX_GPU_DEVICE (const int i, const int j, const int k) -> ReduceTuple
    {
        const auto domain_xlo = geom.ProbLo();
        
        return {state_TopCharge(i,j,k,smeared_fab)};
    });
  }
    ReduceTuple reduced_values = reduce_data.value();
    // MPI reduction
    ParallelDescriptor::ReduceRealSum(amrex::get<0>(reduced_values));
    Real TopCharge = amrex::get<0>(reduced_values);
    
    delete smearedU_mf;
    
    return TopCharge/(2.0*M_PI);
}

Real AmrCoreAdv::meas_TopCharge_Lev (amrex::Vector<amrex::MultiFab>& state, int lev, const amrex::Real time, const amrex::Geometry& geom, Parameters Param)
{
    if(lev != finest_level)
    {
        auto& state_mf = state[lev];
        
        const BoxArray& ba = state_mf.boxArray(); 
        const DistributionMapping& dm = state_mf.DistributionMap();

        BoxArray fba = state[lev+1].boxArray();
        const iMultiFab& mask = makeFineMask(state[lev], fba, IntVect(2));  // the last argument is the refinement ratio
        
        MultiFab* U_mf = new MultiFab;
        U_mf -> define(ba, dm, 4, NUM_GHOST_CELLS); 
        MultiFab::Copy(*U_mf, state_mf, Idx::U_0_Real, cIdx::Real_0, 4, NUM_GHOST_CELLS);
        
        MultiFab* smearedU_mf = new MultiFab;
        smearedU_mf -> define(ba, dm, 4, NUM_GHOST_CELLS); 
        
        smear_gauge(*smearedU_mf, *U_mf, lev, time, Param);

        ReduceOps<ReduceOpSum> reduce_operations;
        ReduceData<Real> reduce_data(reduce_operations);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        
        
#ifdef _OPENMP
#pragma omp parallel
#endif
      for ( MFIter mfi(*smearedU_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
      {
        const Box& bx = mfi.tilebox();
        const auto ncomp = state_mf.nComp();

        const auto& smearedU_fab = smearedU_mf -> array(mfi);
        const auto& fine_mask_arr = mask.array(mfi);

        // For each grid, loop over all the valid points
        reduce_operations.eval(bx, reduce_data,
        [=] AMREX_GPU_DEVICE (const int i, const int j, const int k) -> ReduceTuple
        {
            return (!fine_mask_arr(i, j, k)*state_TopCharge(i,j,k, smearedU_fab));
        });
      }
        ReduceTuple reduced_values = reduce_data.value();
        // MPI reduction
        ParallelDescriptor::ReduceRealSum(amrex::get<0>(reduced_values));
        Real action = amrex::get<0>(reduced_values);
    
        
        delete U_mf;
        delete smearedU_mf;
        return action/(2.0*M_PI);
    }
    
    else
    {
        auto& state_mf = state[lev];
        
        const BoxArray& ba = state_mf.boxArray(); 
        const DistributionMapping& dm = state_mf.DistributionMap();
        
        MultiFab* U_mf = new MultiFab;
        U_mf -> define(ba, dm, 4, NUM_GHOST_CELLS); 
        MultiFab::Copy(*U_mf, state_mf, Idx::U_0_Real, cIdx::Real_0, 4, NUM_GHOST_CELLS);
        
        MultiFab* smearedU_mf = new MultiFab;
        smearedU_mf -> define(ba, dm, 4, NUM_GHOST_CELLS); 
        
        smear_gauge(*smearedU_mf, *U_mf, lev, time, Param);

        ReduceOps<ReduceOpSum> reduce_operations;
        ReduceData<Real> reduce_data(reduce_operations);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        
        
#ifdef _OPENMP
#pragma omp parallel
#endif
      for ( MFIter mfi(*smearedU_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
      {
        const Box& bx = mfi.tilebox();
        const auto ncomp = state_mf.nComp();

        const auto& smearedU_fab = smearedU_mf -> array(mfi);

        // For each grid, loop over all the valid points
        reduce_operations.eval(bx, reduce_data,
        [=] AMREX_GPU_DEVICE (const int i, const int j, const int k) -> ReduceTuple
        {
            return (state_TopCharge(i,j,k, smearedU_fab));
        });
      }
        ReduceTuple reduced_values = reduce_data.value();
        // MPI reduction
        ParallelDescriptor::ReduceRealSum(amrex::get<0>(reduced_values));
        Real action = amrex::get<0>(reduced_values);
    
        
        delete U_mf;
        delete smearedU_mf;
        return action/(2.0*M_PI);
    }
    
    
}

Real AmrCoreAdv::meas_Total_TopCharge_Lev(amrex::Vector<amrex::MultiFab>& state, const amrex::Real time, const amrex::Vector<amrex::Geometry>& geom, Parameters Param, int high_lev)
{
    Real TopCharge = 0;
    for(int i = 0; i <= high_lev; i++)
    {
        TopCharge += meas_TopCharge_Lev(state, i, time, geom[i], Param);   
    }

    return TopCharge;
}

void AmrCoreAdv::meas_SigmaLoops(MultiFab& U_mf, int lev, amrex::Real time, const Geometry& geom, Parameters Param, amrex::Vector<double> &sigma, std::ofstream& Sigma)
{
    BL_PROFILE("AmrCoreAdv::meas_SigmaLoops"); //Compile with TINY_PROFILE = TRUE
    
    const BoxArray& ba = U_mf.boxArray();
    BoxArray nba = ba;
    //nba.surroundingNodes();
    

    const DistributionMapping& dm = U_mf.DistributionMap();
    
    MultiFab smearedU_mf;
    smearedU_mf.define(nba, dm, 4, NUM_GHOST_CELLS);
    
    smear_gauge(smearedU_mf, U_mf, lev, time, Param);
     
    Box domain_bx = geom.Domain();
    //domain_bx.surroundingNodes();
    
    const BoxArray domain_ba(domain_bx);
    
    BoxArray domain_nba = domain_ba;
    //domain_nba.surroundingNodes();
    

    DistributionMapping domain_dm{domain_nba};
    
    MultiFab domain_mf(domain_nba, domain_dm, 4, NUM_GHOST_CELLS);

    
    domain_mf.ParallelCopy(smearedU_mf, cIdx::Real_0, cIdx::Real_0, 4, NUM_GHOST_CELLS, NUM_GHOST_CELLS, geom.periodicity());

    //MultiFab::Copy(domain_mf, *smearedU_mf, cIdx::Real_0, cIdx::Real_0, 4, NUM_GHOST_CELLS);
    //domain_mf.ParallelCopy(U_mf, cIdx::Real_0, cIdx::Real_0, 4, geom.periodicity());
    //Real plaq_ave = Measure_Smeared_Plaq(domain_mf, lev, time, geom, Param);
    
    int Nx = std::pow(2,lev)*Param.Nx;
    int Ny = std::pow(2,lev)*Param.Ny;
    
    amrex::Vector<std::vector<std::complex<double>>> wLoops(Nx/2, std::vector<std::complex<double>> (Ny/2, 0.0));  
    //amrex::Vector<double> sigma(Nx/2-1, 0.0);
    
    Real plaq_ave_alias = Measure_Plaq(U_mf, lev);
    //Real plaq_ave_alias = Measure_Smeared_Plaq(domain_mf, lev, time, geom, Param);
    //amrex::Print() << plaq_ave;
        
    sigma[0] = -std::log(std::abs(plaq_ave_alias));
    
    int p1, p2;

    double inv_Lsq = 1.0/(Nx*Ny);
    
    int loop_max = Nx/2;
    
    for ( MFIter mfi(domain_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        const auto& domain_fab = domain_mf.array(mfi);

        #ifdef _OMP
        #pragma omp parallel for collapse(2)
        int num_threads;
        std::sscanf(std::getenv("OMP_NUM_THREADS"), "%d", &num_threads);
        omp_set_num_threads(4);
        #endif
        
        
        for(int Xrect=1; Xrect < loop_max; Xrect++) 
        {
      
        //Loop over all Y side sizes of rectangle
            for(int Yrect=1; Yrect< loop_max; Yrect++) 
            {
              //Loop over all x,y starting points
              for(int x=0; x<Nx; x++)
              {
                for(int y=0; y<Ny; y++)
                {
                    std::complex<double> w = std::complex<double>(1.0,0.0);
                    std::complex<double> dw_tmp;
                    
                    for(int dx = 0; dx < Xrect; dx++) 
                    {
                        w *= std::complex<double>(domain_fab((x+dx)%Nx,y,0,cIdx::Real_0),
                                                  domain_fab((x+dx)%Nx,y,0,cIdx::Imaginary_0));
                    }
                    
                    p1 = (x + Xrect)%Nx;
                    
                    for(int dy = 0; dy < Yrect; dy++) 
                    {
                        w *= std::complex<double>(domain_fab(p1,(y+dy)%Ny,0,cIdx::Real_1), 
                                                  domain_fab(p1,(y+dy)%Ny,0,cIdx::Imaginary_1));
                    }
                    
                    p2 = (y + Yrect)%Ny;
                    
                    for(int dx = Xrect - 1; dx >= 0; dx--)
                    {
                        dw_tmp = std::complex<double>(domain_fab((x+dx)%Nx,p2,0,cIdx::Real_0),
                                                      domain_fab((x+dx)%Nx,p2,0,cIdx::Imaginary_0));
                        dw_tmp = std::conj(dw_tmp);
                        w *= dw_tmp;
                    }
                    
                    for(int dy = Yrect - 1; dy >= 0; dy--)
                    {
                        dw_tmp = std::complex<double>(domain_fab(x,(y+dy)%Ny,0,cIdx::Real_1),
                                                      domain_fab(x,(y+dy)%Ny,0,cIdx::Imaginary_1));
                        dw_tmp = std::conj(dw_tmp);
                        w *= dw_tmp;
                    }
                    
                    wLoops[Xrect][Yrect] += w*inv_Lsq; 
                    
                }
              }
                
            }

        }
        
        #ifdef _OMP
        omp_set_num_threads(num_threads);
        #endif        
    }


#ifdef _OMP
#pragma omp parallel for
#endif

        for(int size = 2; size<loop_max; size++) 
        {
            sigma[size-1]  = -std::log(std::abs((std::real(wLoops[size][size])/std::real(wLoops[size-1][size]))*(std::real(wLoops[size-1][size-1])/std::real(wLoops[size][size-1]))));
    
            sigma[size-1] += -std::log(std::abs((std::real(wLoops[size][size])/std::real(wLoops[size][size-1]))*(std::real(wLoops[size-1][size-1])/std::real(wLoops[size-1][size]))));
    
            sigma[size-1] *= 0.5;
            //amrex::Print() << sigma[size - 1] << ", ";
        }
        
        for(int i = 0; i<sigma.size(); i++)
        {
            amrex::Print() << std::fixed << std::setprecision(9) << sigma[i] << ' ';
            Sigma << std::fixed << std::setprecision(9) << sigma[i] << ' ';
        }
        amrex::Print() << std::endl;
        Sigma << std::endl;
        
    
    //delete smearedU_mf;
}

void AmrCoreAdv::meas_WilsonLoops(MultiFab& U_mf, int lev, amrex::Real time, int step, const Geometry& geom, Parameters Param, amrex::Vector<double> &wilsonLoop, std::ofstream& WilsonLoop)
{
    BL_PROFILE("AmrCoreAdv::meas_WilsonLoops"); //Compile with TINY_PROFILE = TRUE
    
    const BoxArray& ba = U_mf.boxArray();
    BoxArray nba = ba;
    //nba.surroundingNodes();
    

    const DistributionMapping& dm = U_mf.DistributionMap();
    
    MultiFab smearedU_mf;
    smearedU_mf.define(nba, dm, 4, NUM_GHOST_CELLS);
    
    smear_gauge(smearedU_mf, U_mf, lev, time, Param);
     
    Box domain_bx = geom.Domain();
    //domain_bx.surroundingNodes();
    
    const BoxArray domain_ba(domain_bx);
    
    BoxArray domain_nba = domain_ba;
    //domain_nba.surroundingNodes();
    

    DistributionMapping domain_dm{domain_nba};
    
    MultiFab domain_mf(domain_nba, domain_dm, 4, NUM_GHOST_CELLS);

    
    domain_mf.ParallelCopy(smearedU_mf, cIdx::Real_0, cIdx::Real_0, 4, NUM_GHOST_CELLS, NUM_GHOST_CELLS, geom.periodicity());

    //MultiFab::Copy(domain_mf, *smearedU_mf, cIdx::Real_0, cIdx::Real_0, 4, NUM_GHOST_CELLS);
    //domain_mf.ParallelCopy(U_mf, cIdx::Real_0, cIdx::Real_0, 4, geom.periodicity());
    //Real plaq_ave = Measure_Smeared_Plaq(domain_mf, lev, time, geom, Param);
    
    int Nx = std::pow(2,lev)*Param.Nx;
    int Ny = std::pow(2,lev)*Param.Ny;
    
    amrex::Vector<std::vector<std::complex<double>>> wLoops(Nx/2, std::vector<std::complex<double>> (Ny/2, 0.0));  
    //amrex::Vector<double> sigma(Nx/2-1, 0.0);
    
    Real plaq_ave_alias = Measure_Plaq(U_mf, lev);
    //Real plaq_ave_alias = Measure_Smeared_Plaq(domain_mf, lev, time, geom, Param);
    //amrex::Print() << plaq_ave;
        
    wilsonLoop[0] = plaq_ave_alias;
    
    int p1, p2;

    double inv_Lsq = 1.0/(Nx*Ny);
    
    int loop_max = Nx/2;
    
    for ( MFIter mfi(domain_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        const auto& domain_fab = domain_mf.array(mfi);

        #ifdef _OMP
        #pragma omp parallel for collapse(2)
        int num_threads;
        std::sscanf(std::getenv("OMP_NUM_THREADS"), "%d", &num_threads);
        omp_set_num_threads(4);
        #endif
        
        
        for(int Xrect=1; Xrect < loop_max; Xrect++) 
        {
      
        //Loop over all Y side sizes of rectangle
            for(int Yrect=1; Yrect< loop_max; Yrect++) 
            {
              //Loop over all x,y starting points
              for(int x=0; x<Nx; x++)
              {
                for(int y=0; y<Ny; y++)
                {
                    std::complex<double> w = std::complex<double>(1.0,0.0);
                    std::complex<double> dw_tmp;
                    
                    for(int dx = 0; dx < Xrect; dx++) 
                    {
                        w *= std::complex<double>(domain_fab((x+dx)%Nx,y,0,cIdx::Real_0),
                                                  domain_fab((x+dx)%Nx,y,0,cIdx::Imaginary_0));
                    }
                    
                    p1 = (x + Xrect)%Nx;
                    
                    for(int dy = 0; dy < Yrect; dy++) 
                    {
                        w *= std::complex<double>(domain_fab(p1,(y+dy)%Ny,0,cIdx::Real_1), 
                                                  domain_fab(p1,(y+dy)%Ny,0,cIdx::Imaginary_1));
                    }
                    
                    p2 = (y + Yrect)%Ny;
                    
                    for(int dx = Xrect - 1; dx >= 0; dx--)
                    {
                        dw_tmp = std::complex<double>(domain_fab((x+dx)%Nx,p2,0,cIdx::Real_0),
                                                      domain_fab((x+dx)%Nx,p2,0,cIdx::Imaginary_0));
                        dw_tmp = std::conj(dw_tmp);
                        w *= dw_tmp;
                    }
                    
                    for(int dy = Yrect - 1; dy >= 0; dy--)
                    {
                        dw_tmp = std::complex<double>(domain_fab(x,(y+dy)%Ny,0,cIdx::Real_1),
                                                      domain_fab(x,(y+dy)%Ny,0,cIdx::Imaginary_1));
                        dw_tmp = std::conj(dw_tmp);
                        w *= dw_tmp;
                    }
                    
                    wLoops[Xrect][Yrect] += w*inv_Lsq; 
                    
                }
              }
                
            }

        }
        
        #ifdef _OMP
        omp_set_num_threads(num_threads);
        #endif        
        
    }

#ifdef _OMP
#pragma omp parallel for
#endif

        for(int size = 2; size<loop_max; size++) 
        {
            
            wilsonLoop[size-1] = std::abs(wLoops[size-1][size-1]);
           
        }
        
        
        WilsonLoop << step << " ";
        for(int i = 0; i<wilsonLoop.size(); i++)
        {
            amrex::Print() << std::fixed << std::setprecision(9) << wilsonLoop[i] << ' ';
            WilsonLoop << std::fixed << std::setprecision(9) << wilsonLoop[i] << ' ';
        }
        amrex::Print() << std::endl;
        WilsonLoop << std::endl;
        
    
    
    //delete smearedU_mf;
}

Real AmrCoreAdv::Measure_Plaq (MultiFab& U_mf, int lev)
{
    //auto& nodal_mask = grid_msk[lev];
    
    ReduceOps<ReduceOpSum> reduce_operations;
    ReduceData<Real> reduce_data(reduce_operations);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif

  for ( MFIter mfi(U_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
  {
    const Box& bx = mfi.tilebox();
    const auto ncomp = U_mf.nComp();

    const auto& U_fab = U_mf.array(mfi); 
    // auto& mask_arr = nodal_mask -> array(mfi);

    // For each grid, loop over all the valid points
    reduce_operations.eval(bx, reduce_data,
    [=] AMREX_GPU_DEVICE (const int i, const int j, const int k) -> ReduceTuple
    {
        return {/*mask_arr(i, j, k)*/sum_plaq(i,j,k,U_fab)/(std::pow(2,2*lev)*numcells[0]*numcells[1])};
    });
  }
    ReduceTuple reduced_values = reduce_data.value();
    // MPI reduction
    ParallelDescriptor::ReduceRealSum(amrex::get<0>(reduced_values));
    Real action = amrex::get<0>(reduced_values);
    
    return action;
}

void AmrCoreAdv::PionCorrelation (MultiFab& state, int lev, const Real time_lev, int step, const amrex::Geometry& geom, Parameters Param, amrex::Vector<double>& picorr, std::ofstream& PiCorr, bool grid_fractured)
{
    
    //auto& nodal_mask = grid_msk[lev];
    
    const BoxArray& ba = state.boxArray();
    BoxArray nba = ba;
    //nba.surroundingNodes();
    const DistributionMapping& dm = state.DistributionMap();
    MultiFab U_mf(nba, dm, 4, NUM_GHOST_CELLS);
    
    MultiFab propUp_mf(nba, dm, 4, NUM_GHOST_CELLS);
    
    MultiFab propDn_mf(nba, dm, 4, NUM_GHOST_CELLS);
    
    MultiFab source_mf(nba, dm, 4, NUM_GHOST_CELLS);
    
    MultiFab Dsource_mf(nba, dm, 4, NUM_GHOST_CELLS);
    
    MultiFab source_tmp_mf(nba, dm, 4, NUM_GHOST_CELLS);

    state.MultiFab::setVal(0, Idx::Source_0_Real, 4, NUM_GHOST_CELLS);
    Dsource_mf.MultiFab::setVal(0, cIdx::Real_0, 4, NUM_GHOST_CELLS);
    propUp_mf.MultiFab::setVal(0, cIdx::Real_0, 4, NUM_GHOST_CELLS);
    
    for ( MFIter mfi(state, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.tilebox();
        const auto ncomp = state.nComp();

        const auto& source_fab = state.array(mfi); 
        //const auto& mask_arr = nodal_mask -> array(mfi);
        
        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if(i == 0 && j == 0)
              source_fab(i, j, k, Idx::Source_0_Real) = 1;
        });
        

    }
    
    
    MultiFab::Swap(source_mf, state, Idx::Source_0_Real, cIdx::Real_0, 4, NUM_GHOST_CELLS);
    Set_g3p(source_tmp_mf, source_mf, lev, time_lev);
    
    
    
    MultiFab::Swap(state, source_tmp_mf, cIdx::Real_0, Idx::Source_0_Real, 4, NUM_GHOST_CELLS);
    
    FillPatch(lev, time_lev, state, 0, state.nComp(), grid_fractured); //new
    //FlipSigns(lev, state, Idx::Source_0_Real, 4);
    
    
    MultiFab::Swap(source_tmp_mf, state,  Idx::Source_0_Real, cIdx::Real_0, 4, NUM_GHOST_CELLS);
    
    
    
    
    MultiFab::Swap(U_mf, state, Idx::U_0_Real, cIdx::Real_0, 4, NUM_GHOST_CELLS);
    
    Set_Dp(Dsource_mf, source_tmp_mf, U_mf, lev, time_lev, Param);
    
    MultiFab::Swap(U_mf, state, Idx::U_0_Real, cIdx::Real_0, 4, NUM_GHOST_CELLS);
    
    Set_g3p(source_mf, Dsource_mf, lev, time_lev);
    
    
    MultiFab::Swap(state, source_mf, cIdx::Real_0, Idx::Source_0_Real, 4, NUM_GHOST_CELLS);
    FillPatch(lev, time_lev, state, 0, state.nComp(), grid_fractured);
    //FlipSigns(lev, state, Idx::Source_0_Real, 4);
    
    
    
    MultiFab::Swap(source_mf, state, Idx::Source_0_Real, cIdx::Real_0, 4, NUM_GHOST_CELLS);
    
    BiCG_Solve(propUp_mf, source_mf, state, 0, lev, time_lev, geom, Param, grid_fractured);
    
    MultiFab::Swap(state, propUp_mf, cIdx::Real_0, Idx::propUp_0_Real, 4, NUM_GHOST_CELLS);
    FillPatch(lev, time_lev, state, 0, state.nComp(), grid_fractured);
    //FlipSigns(lev, state, Idx::propUp_0_Real, 4);
    MultiFab::Swap(propUp_mf, state, Idx::propUp_0_Real,  cIdx::Real_0, 4, NUM_GHOST_CELLS);
        

    state.MultiFab::setVal(0, Idx::Source_0_Real, 4, NUM_GHOST_CELLS);
    Dsource_mf.MultiFab::setVal(0, cIdx::Real_0, 4, NUM_GHOST_CELLS);
    propDn_mf.MultiFab::setVal(0, cIdx::Real_0, 4, NUM_GHOST_CELLS);
    
    for ( MFIter mfi(state, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.tilebox();
        const auto ncomp = state.nComp();

        const auto& source_fab = state.array(mfi); 
        //const auto& mask_arr = nodal_mask -> array(mfi);
        
        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if(i == 0 && j == 0)
              source_fab(i, j, k, Idx::Source_1_Real) = 1;
        });
        

    }
    
    
    
    MultiFab::Swap(source_mf, state, Idx::Source_0_Real, cIdx::Real_0, 4, NUM_GHOST_CELLS);
    Set_g3p(source_tmp_mf, source_mf, lev, time_lev);
    
    
    
    MultiFab::Swap(state, source_tmp_mf, cIdx::Real_0, Idx::Source_0_Real, 4, NUM_GHOST_CELLS);
    
    FillPatch(lev, time_lev, state, 0, state.nComp(), grid_fractured); //new
    //FlipSigns(lev, state, Idx::Source_0_Real, 4);
    
    
    MultiFab::Swap(source_tmp_mf, state,  Idx::Source_0_Real, cIdx::Real_0, 4, NUM_GHOST_CELLS);
    
    
    
    MultiFab::Swap(U_mf, state, Idx::U_0_Real, cIdx::Real_0, 4, NUM_GHOST_CELLS);
    
    Set_Dp(Dsource_mf, source_tmp_mf, U_mf, lev, time_lev, Param);
    
    MultiFab::Swap(U_mf, state, Idx::U_0_Real, cIdx::Real_0, 4, NUM_GHOST_CELLS);
    
    Set_g3p(source_mf, Dsource_mf, lev, time_lev);
    
    
    MultiFab::Swap(state, source_mf, cIdx::Real_0, Idx::Source_0_Real, 4, NUM_GHOST_CELLS);
    FillPatch(lev, time_lev, state, 0, state.nComp(), grid_fractured);
    //FlipSigns(lev, state, Idx::Source_0_Real, 4);
    
    
    
    
    MultiFab::Swap(source_mf, state, Idx::Source_0_Real, cIdx::Real_0, 4, NUM_GHOST_CELLS);
    
    BiCG_Solve(propDn_mf, source_mf, state, 0, lev, time_lev, geom, Param, grid_fractured);
    
    MultiFab::Swap(state, propDn_mf, cIdx::Real_0, Idx::propDn_0_Real, 4, NUM_GHOST_CELLS);
    FillPatch(lev, time_lev, state, 0, state.nComp(), grid_fractured);
    //FlipSigns(lev, state, Idx::propDn_0_Real, 4);
    MultiFab::Swap(propDn_mf, state, Idx::propDn_0_Real,  cIdx::Real_0, 4, NUM_GHOST_CELLS);
    
    
    const Box& domain_bx = geom.Domain();
    const BoxArray domain_ba(&domain_bx, 1);
    BoxArray domain_nba = domain_ba;
    
    
    DistributionMapping domain_dm {domain_nba};
    MultiFab propUp_domain_mf(domain_nba, domain_dm, 4, NUM_GHOST_CELLS);
    
    MultiFab propDn_domain_mf(domain_nba, domain_dm, 4, NUM_GHOST_CELLS);
    
    propUp_domain_mf.MultiFab::setVal(0, cIdx::Real_0, 4, NUM_GHOST_CELLS);
    propDn_domain_mf.MultiFab::setVal(0, cIdx::Real_0, 4, NUM_GHOST_CELLS);
    
    //FillIntermediatePatch(lev, time_lev, propUp_mf, 0, 4);
    //FillIntermediatePatch(lev, time_lev, propDn_mf, 0, 4);
    propUp_domain_mf.ParallelCopy(propUp_mf, cIdx::Real_0, cIdx::Real_0, 4, NUM_GHOST_CELLS, NUM_GHOST_CELLS, geom.periodicity());
    propDn_domain_mf.ParallelCopy(propDn_mf, cIdx::Real_0, cIdx::Real_0, 4, NUM_GHOST_CELLS, NUM_GHOST_CELLS, geom.periodicity());
    
    int Nx = std::pow(2,lev)*Param.Nx;
    int Ny = std::pow(2,lev)*Param.Ny;
    
    double corr = 0.0; double tmp = 0.0;
    
    for ( MFIter mfi(propUp_domain_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        const auto& propUp_domain_fab = propUp_domain_mf.array(mfi);
        const auto& propDn_domain_fab = propDn_domain_mf.array(mfi);
        
        for(int y = 0; y < Ny; y++) 
        {
            //picorr[y] = 0.0;
            corr = 0.0;
            for(int x = 0; x < Nx; x++)
            {
                tmp = propUp_domain_fab(x, y, 0, cIdx::Real_0)*propUp_domain_fab(x, y, 0, cIdx::Real_0) +
                      propUp_domain_fab(x, y, 0, cIdx::Imaginary_0)*propUp_domain_fab(x, y, 0, cIdx::Imaginary_0) +
                    
                      propUp_domain_fab(x, y, 0, cIdx::Real_1)*propUp_domain_fab(x, y, 0, cIdx::Real_1) +
                      propUp_domain_fab(x, y, 0, cIdx::Imaginary_1)*propUp_domain_fab(x, y, 0, cIdx::Imaginary_1) +
                    
                      propDn_domain_fab(x, y, 0, cIdx::Real_0)*propDn_domain_fab(x, y, 0, cIdx::Real_0) +
                      propDn_domain_fab(x, y, 0, cIdx::Imaginary_0)*propDn_domain_fab(x, y, 0, cIdx::Imaginary_0) +
                    
                      propDn_domain_fab(x, y, 0, cIdx::Real_1)*propDn_domain_fab(x, y, 0, cIdx::Real_1) +
                      propDn_domain_fab(x, y, 0, cIdx::Imaginary_1)*propDn_domain_fab(x, y, 0, cIdx::Imaginary_1);
                
                corr += tmp;
                      
            }
                
            if ( y < ((Ny/2)+1) ) picorr[y] += corr;
            else {
              picorr[Ny-y] += corr;
              picorr[Ny-y] /= 2.0;
            }
            
        }
    
    }
    
    PiCorr << step << ' ';
    for(int i = 0; i < picorr.size(); i++) 
    {
        amrex::Print() << std::fixed << std::setprecision(9) << picorr[i] << ' ';
        PiCorr << std::fixed << std::setprecision(9) << picorr[i] << ' ';
    }
    amrex::Print() << std::endl;
    PiCorr << std::endl;
    
    //pion_corr_out = pion_corr;
    

    
}

void AmrCoreAdv::Set_g3p (MultiFab& g3p_mf, MultiFab& p_mf, int lev, Real time)
{
    //auto nodal_mask = p_mf.OwnerMask(geom[lev].periodicity());
    
#ifdef _OPENMP
#pragma omp parallel
#endif
  for ( MFIter mfi(p_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
  {
    const Box& bx = mfi.tilebox();

    const auto& g3p_fab = g3p_mf.array(mfi);
    const auto& p_fab = p_mf.array(mfi);
    //const auto& mask_arr = nodal_mask->array(mfi);
    
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        //if(mask_arr(i, j, k))
            state_set_g3p(i, j, k, g3p_fab, p_fab);
    });
      
  }
    
    
  //g3p_mf.OverrideSync(*nodal_mask, geom[lev].periodicity());
  //FillIntermediatePatch(lev, time, g3p_mf, 0, g3p_mf.nComp());
  //FlipSigns(lev, g3p_mf, cIdx::Real_0, 4);
}

void AmrCoreAdv::Set_Dp (MultiFab& Dp_mf, MultiFab& p_mf, MultiFab& U_mf, int lev, Real time, Parameters Param)
{
    //const auto dx = geom.CellSizeArray();
    //int ncomp = U_mf.nComp();
    
    //auto nodal_mask = p_mf.OwnerMask(geom[lev].periodicity());
    
#ifdef _OPENMP
#pragma omp parallel
#endif
  for ( MFIter mfi(U_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
  {
    const Box& bx = mfi.tilebox();
    const auto ncomp = U_mf.nComp();

    const auto& Dp_fab = Dp_mf.array(mfi);
    const auto& p_fab = p_mf.array(mfi); 
    const auto& U_fab = U_mf.array(mfi);
    //const auto& mask_arr = nodal_mask->array(mfi);
    
    
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        //if(mask_arr(i, j, k))
            state_set_Dp(i, j, k, Dp_fab, p_fab, U_fab, Param.mass);
    });
      
  }
  
  //Dp_mf.OverrideSync(*nodal_mask, geom[lev].periodicity());
  //FillIntermediatePatch(lev, time, Dp_mf, 0, Dp_mf.nComp());
  //FlipSigns(lev, Dp_mf, cIdx::Real_0, 4);
    
}

void
AmrCoreAdv::BiCG_Solve(MultiFab& x_mf, MultiFab& b_mf, MultiFab& state, bool useInitGuess, int lev, Real time, const amrex::Geometry& geom_lev, Parameters Param, bool grid_fractured)
{
    BL_PROFILE("BiCG_solve");
    amrex::Real alpha, beta, denom, rsq, rsq_new, bsqrt, bnorm;
    
    //MultiFab& state = grid_new[lev];
    
    
    const BoxArray& ba = x_mf.boxArray();
    
    const DistributionMapping& dm = x_mf.DistributionMap();
    
    MultiFab U_mf;
    U_mf.define(ba, dm, 4, NUM_GHOST_CELLS);
                
    MultiFab::Copy(U_mf, state, Idx::U_0_Real, cIdx::Real_0, 4, NUM_GHOST_CELLS);
    
    //MultiFab* mfres = new MultiFab;
    //mfres -> define(ba,dm,4,NUM_GHOST_CELLS);
    
    //FillIntermediatePatch(lev, time, b_mf, 0, b_mf.nComp()); //new
    MultiFab::Copy(state, b_mf, cIdx::Real_0, Idx::res_0_Real, 4, NUM_GHOST_CELLS);
    
    x_mf.MultiFab::setVal(0, cIdx::Real_0, 4, NUM_GHOST_CELLS);  
    
    
    MultiFab mfp;
    mfp.define(ba,dm,4,NUM_GHOST_CELLS);
    
    FillPatch(lev, time, state, 0, state.nComp(), grid_fractured); //new  
    MultiFab::Copy(state, state, Idx::res_0_Real, Idx::p_0_Real, 4, NUM_GHOST_CELLS);
    
    MultiFab mfDp;
    mfDp.define(ba,dm,4,NUM_GHOST_CELLS);
    
    MultiFab mfDDp;
    mfDDp.define(ba,dm,4,NUM_GHOST_CELLS);
    
    rsq = MultiFab::Dot(state, Idx::res_0_Real,  state, Idx::res_0_Real, 4, 0);
    
    //amrex::Print() << rsq << std::endl;
    
    for(int i = 0; i < Param.BiCG_Max_Iters; i++)
    {
       //MultiFab::Copy(mfp, state, Idx::p_0_Real, cIdx::Real_0, 4, NUM_GHOST_CELLS);

        
        FillPatch(lev, time, state, 0, state.nComp(), grid_fractured); //new
        //FlipSigns(lev, state, Idx::p_0_Real, 4);
        //FlipSigns(lev, state, Idx::Dp_0_Real, 4);
        //FlipSigns(lev, state, Idx::DDp_0_Real, 4);
        
        MultiFab::Swap(mfp, state, Idx::p_0_Real, cIdx::Real_0, 4, NUM_GHOST_CELLS);
        MultiFab::Swap(mfDp, state, Idx::Dp_0_Real, cIdx::Real_0, 4, NUM_GHOST_CELLS);
        
        Set_Dp(mfDp, mfp, U_mf, lev, time, Param);
        
        
        
        //FillIntermediatePatch(lev, time, *mfDp, 0, 4); //new
        Set_g3p(mfDp, mfDp, lev, time);
        
        MultiFab::Swap(mfp, state, Idx::p_0_Real, cIdx::Real_0, 4, NUM_GHOST_CELLS);
        MultiFab::Swap(mfDp, state, Idx::Dp_0_Real, cIdx::Real_0, 4, NUM_GHOST_CELLS);
        
        FillPatch(lev, time, state, 0, state.nComp(), grid_fractured); //new
        //FlipSigns(lev, state, Idx::p_0_Real, 4);
        //FlipSigns(lev, state, Idx::Dp_0_Real, 4);
        //FlipSigns(lev, state, Idx::DDp_0_Real, 4);
        
        //amrex::Print() << MultiFab::Dot(state, Idx::p_0_Real, state, Idx::p_0_Real, 4, NUM_GHOST_CELLS) << " ";
        //amrex::Print() << MultiFab::Dot(state, Idx::Dp_0_Real, state, Idx::Dp_0_Real, 4, NUM_GHOST_CELLS) << std::endl;
        
        MultiFab::Swap(mfDp, state, Idx::Dp_0_Real, cIdx::Real_0, 4, NUM_GHOST_CELLS);
        MultiFab::Swap(mfDDp, state, Idx::DDp_0_Real, cIdx::Real_0, 4, NUM_GHOST_CELLS);
        
        Set_Dp(mfDDp, mfDp, U_mf, lev, time, Param);
        
        //FillIntermediatePatch(lev, time, *mfDDp, 0, 4); //new
        Set_g3p(mfDDp, mfDDp, lev, time);
        
        MultiFab::Swap(mfDp, state, Idx::Dp_0_Real, cIdx::Real_0, 4, NUM_GHOST_CELLS);
        MultiFab::Swap(mfDDp, state, Idx::DDp_0_Real, cIdx::Real_0, 4, NUM_GHOST_CELLS);
        
        
        
        //FillIntermediatePatch(lev, time, *mfDDp, 0, 4);
        
        
        denom = MultiFab::Dot(state, Idx::p_0_Real, state, Idx::DDp_0_Real, 4,  0);
        
        alpha = rsq/denom;
        
        FillPatch(lev, time, state, 0, state.nComp(), grid_fractured);
        //FlipSigns(lev, state, Idx::p_0_Real, 4);
        MultiFab::Saxpy(x_mf, alpha, state, Idx::p_0_Real, cIdx::Real_0, 4, 0);
        

        //x_mf.OverrideSync(*nodal_mask, geom_lev.periodicity());
        //FillIntermediatePatch(lev, time, x_mf, 0, x_mf.nComp());
        //FlipSigns(lev, x_mf, cIdx::Real_0, 4);
        
        FillPatch(lev, time, state, 0, state.nComp(), grid_fractured);
        //FlipSigns(lev, state, Idx::DDp_0_Real, 4);
        MultiFab::Saxpy(state, -alpha, state, Idx::DDp_0_Real, Idx::res_0_Real, 4, 0);
              
        //mfres -> OverrideSync(*nodal_mask, geom_lev.periodicity());
        //FillIntermediatePatch(lev, time, *mfres, 0, mfres -> nComp());
        //FlipSigns(lev, *mfres, cIdx::Real_0, 4);

        rsq_new = MultiFab::Dot(state, Idx::res_0_Real, state, Idx::res_0_Real, 4, 0);
        //amrex::Print() << rsq_new << std::endl;
        if (rsq_new < Param.BiCG_Thresh)
            break;
        
        beta = rsq_new/rsq;
        
        //amrex::Print() << rsq << ' ' << rsq_new << std::endl;
        
        rsq = rsq_new;
        
        //amrex::Print() << rsq << std::endl;
        
        MultiFab::LinComb(state, beta, state, Idx::p_0_Real, 1.0, state, Idx::res_0_Real, Idx::p_0_Real, 4, 0);
        //amrex::Print() << MultiFab::Dot(*mfp, cIdx::Real_0, *mfp, cIdx::Real_0, 4, NUM_GHOST_CELLS) << std::endl;  
        //mfp -> OverrideSync(*nodal_mask, geom_lev.periodicity());
        //FillIntermediatePatch(lev, time, *mfp, 0, mfp -> nComp());
        //FlipSigns(lev, *mfp, cIdx::Real_0, 4);
        
        if (i == (Param.BiCG_Max_Iters - 1) && rsq_new > Param.BiCG_Thresh)
        {
            amrex::Print() << "Failed to converge after " << Param.BiCG_Max_Iters << " steps!" << std::endl;
            amrex::Print() << "Final rsq = " << rsq_new << std::endl;
            
        }
        
    }
    
    //amrex::Print() << MultiFab::Dot(x_mf,cIdx::Real_0,x_mf,cIdx::Real_0,4,0) << std::endl;
    
    //MultiFab::Swap(*mfU, S_new, Idx::U_0_Real, cIdx::Real_0, 4, NUM_GHOST_CELLS);
    
    //delete mfres;
    //delete mfU;
    //delete mfp;
    //delete mfDp;
    //delete mfDDp;
    
    //x_mf.OverrideSync(*nodal_mask, geom_lev.periodicity());
    //FillIntermediatePatch(lev, time, x_mf, 0, x_mf.nComp());
    //FlipSigns(lev, x_mf, cIdx::Real_0, 4);
    
}

