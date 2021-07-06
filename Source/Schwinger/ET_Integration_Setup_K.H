#ifndef ET_INTEGRATION_SETUP_K_H
#define ET_INTEGRATION_SETUP_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 

namespace Idx { 
         enum ETIndexes {U_0_Real, U_0_Imaginary, U_1_Real, U_1_Imaginary, Ucurrent_0_Real, Ucurrent_0_Imaginary, Ucurrent_1_Real, Ucurrent_1_Imaginary, P_0, Pcurrent_0, P_1, Pcurrent_1, NumScalars}; 
};

namespace InitIdx { 
         enum ETInitIndexes {NumScalars}; 
};

namespace Diag { 
         enum DiagnosticIndexes {ActionDensity, AbsU_0_Re, AbsU_1_Re, NumScalars}; 
};

#define NUM_GHOST_CELLS 1

#endif