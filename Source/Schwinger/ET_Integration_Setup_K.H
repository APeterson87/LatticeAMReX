#ifndef ET_INTEGRATION_SETUP_K_H
#define ET_INTEGRATION_SETUP_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 

namespace Idx { 
  enum ETIndexes {// Gauge
		  U_0_Real,
		  U_0_Imag,
		  U_1_Real,
		  U_1_Imag,
		  // mom
		  Pi_0,
		  Pi_1,
		  // Fermion spin 0
		  Phi_0_Real, 
		  Phi_0_Imag,
		  // Fermion spin 1
		  Phi_1_Real,
		  Phi_1_Imag,
		  // Gauge force
		  F_U_0,
		  F_U_1,
		  // Fermion force
		  F_Phi_0,
		  F_Phi_1,
		  // measurements
		  ActionDensity,
		  TopChargeDensity,
		  NumScalars}; 
};

namespace auxIdx { 
         enum auxIndexes {DDinvPhi_0_Real, 
			  DDinvPhi_0_Imag, 
			  DDinvPhi_1_Real, 
			  DDinvPhi_1_Imag, 
			  g3DinvPhi_0_Real, 
			  g3DinvPhi_0_Imag, 
			  g3DinvPhi_1_Real, 
			  g3DinvPhi_1_Imag,
			  NumScalars}; 
};

namespace cIdx {
         enum cIndexes {Real_0,
                        Imag_0,
                        Real_1,
                        Imag_1,
                        NumIndexes};
};

namespace InitIdx { 
         enum ETInitIndexes {NumScalars}; 
};

namespace Diag { 
         enum DiagnosticIndexes {NumScalars}; 
};

#define NUM_GHOST_CELLS 1

#endif
