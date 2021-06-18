#ifndef ET_INTEGRATION_SETUP_K_H
#define ET_INTEGRATION_SETUP_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 

namespace Idx { 
         enum ETIndexes {U_0_Re, U_0_Im, U_1_Re, U_1_Im, NumScalars}; 
};

namespace InitIdx { 
         enum ETInitIndexes {NumScalars}; 
};

namespace Diag { 
         enum DiagnosticIndexes {ActionDensity, AbsU_0_Re, AbsU_1_Re, NumScalars}; 
};

#define NUM_GHOST_CELLS 1

#endif