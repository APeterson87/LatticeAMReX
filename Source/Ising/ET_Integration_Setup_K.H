#ifndef ET_INTEGRATION_SETUP_K_H
#define ET_INTEGRATION_SETUP_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 

namespace Idx { 
         enum ETIndexes {s, NumScalars}; 
};

namespace InitIdx { 
         enum ETInitIndexes {NumScalars}; 
};

namespace Diag { 
         enum DiagnosticIndexes {Edensity, NumScalars}; 
};

#define NUM_GHOST_CELLS 1

#endif