#ifndef ET_INTEGRATION_UTILS_K_H
#define ET_INTEGRATION_UTILS_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

using namespace amrex;

struct Observables {
            int num_accepted;
            Real TimeSumAvePlaq;
            Real TimeSumDeltaH;
            std::vector<int> TopologicalCharge;
    
            Observables() : num_accepted(0), 
                            TimeSumAvePlaq(0), 
                            TimeSumDeltaH(0) { }
        };
#endif