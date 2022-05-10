#ifndef ET_INTEGRATION_SETUP_K_H
#define ET_INTEGRATION_SETUP_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 

namespace Idx { 
         enum ETIndexes {U_0_Real, 
                         U_0_Imaginary, 
                         U_1_Real, 
                         U_1_Imaginary, 
                         Ucurrent_0_Real, 
                         Ucurrent_0_Imaginary, 
                         Ucurrent_1_Real, 
                         Ucurrent_1_Imaginary, 
                         P_0, 
                         P_1,
                         chi_0_Real, 
                         chi_0_Imaginary, 
                         chi_1_Real, 
                         chi_1_Imaginary,
                         Phi_0_Real, 
                         Phi_0_Imaginary, 
                         Phi_1_Real, 
                         Phi_1_Imaginary, 
                         DDinvPhi_0_Real, 
                         DDinvPhi_0_Imaginary, 
                         DDinvPhi_1_Real, 
                         DDinvPhi_1_Imaginary, 
                         g3DinvPhi_0_Real, 
                         g3DinvPhi_0_Imaginary, 
                         g3DinvPhi_1_Real, 
                         g3DinvPhi_1_Imaginary,
                         TopChargeDensity,
                         ActionDensity,
                         NumScalars}; 
};

namespace cIdx {
         enum cIndexes {Real_0,
                        Imaginary_0,
                        Real_1,
                        Imaginary_1,
                        NumIndexes};
};

namespace pIdx { 
         enum pIndexes {p_0_Real, 
                        p_0_Imaginary, 
                        p_1_Real, 
                        p_1_Imaginary, 
                        p_NumScalars}; 
};

namespace DDpIdx { 
         enum DDpIndexes {DDp_0_Real, 
                        DDp_0_Imaginary, 
                        DDp_1_Real, 
                        DDp_1_Imaginary, 
                        DDp_NumScalars}; 
};

namespace resIdx { 
         enum resIndexes {res_0_Real, 
                          res_0_Imaginary, 
                          res_1_Real, 
                          res_1_Imaginary, 
                          res_NumScalars}; 
};

namespace InitIdx { 
         enum ETInitIndexes {NumScalars}; 
};

namespace Diag { 
         enum DiagnosticIndexes {ActionDensity, AbsU_0_Re, AbsU_1_Re, NumScalars}; 
};

#define NUM_GHOST_CELLS 2

#endif