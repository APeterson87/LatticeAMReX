#ifndef ET_INTEGRATION_SETUP_K_H
#define ET_INTEGRATION_SETUP_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 

namespace Idx { 
         enum ETIndexes {U_0_Real, 
                         U_0_Imaginary, 
                         U_1_Real, 
                         U_1_Imaginary,
                         theta_0,
                         theta_1,
                         P_0, 
                         P_1,
                         Phi_0_Real, 
                         Phi_0_Imaginary, 
                         Phi_1_Real, 
                         Phi_1_Imaginary,
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
         enum DiagnosticIndexes {ActionDensity, AbsU_0_Re, AbsU_1_Re, NumScalars}; 
};

#define NUM_GHOST_CELLS 1

#endif