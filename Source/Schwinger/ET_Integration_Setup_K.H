#ifndef ET_INTEGRATION_SETUP_K_H
#define ET_INTEGRATION_SETUP_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 

namespace Idx { 
  enum ETIndexes {// Gauge
			 U_0_Real,
                         U_0_Imaginary,
                         U_1_Real,
                         U_1_Imaginary,
			 // mom
                         Pi_0,
                         Pi_1,
			 // fermion
                         Phi_0_Real, 
                         Phi_0_Imaginary, 
                         Phi_1_Real, 
                         Phi_1_Imaginary,
			 // measurements
                         ActionDensity,
                         TopChargeDensity,
                         NumScalars}; 
};

namespace auxIdx { 
         enum auxIndexes {DDinvPhi_0_Real, 
			  DDinvPhi_0_Imaginary, 
			  DDinvPhi_1_Real, 
			  DDinvPhi_1_Imaginary, 
			  g3DinvPhi_0_Real, 
			  g3DinvPhi_0_Imaginary, 
			  g3DinvPhi_1_Real, 
			  g3DinvPhi_1_Imaginary,
			  NumScalars}; 
};

namespace cIdx {
         enum cIndexes {Real_0,
                        Imaginary_0,
                        Real_1,
                        Imaginary_1,
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
