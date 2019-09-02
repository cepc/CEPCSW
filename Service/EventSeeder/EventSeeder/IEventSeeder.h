#ifndef I_EVENT_SEEDER_H
#define I_EVENT_SEEDER_H

#include "GaudiKernel/IService.h"

class Algorithm;

class IEventSeeder: virtual public IService {
public:
    DeclareInterfaceID(IEventSeeder, 0, 1); // major/minor version
    
    virtual ~IEventSeeder() = default;

    // Register an algorithm for the seeding service
    virtual void registerAlg( Algorithm* alg ) = 0;

    // Get the seed corelated to current event and algorithm
    virtual unsigned int getSeed(Algorithm* alg, int eventNumber, int runNumber) = 0;
};

#endif
