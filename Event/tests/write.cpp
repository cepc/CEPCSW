#include "EventHeaderCollection.h"
#include "MCParticleCollection.h"

#include "podio/EventStore.h"
#include "podio/ROOTWriter.h"

#include <vector>
#include <iostream>

int main(int argc, char* argv[])
{
    std::cout << "start writing test" << std::endl;

    auto store = podio::EventStore();
    auto writer = podio::ROOTWriter("test.root", &store);

    auto& ehc = store.create<plcio::EventHeaderCollection>("EvtHeaders");
    auto& mcc = store.create<plcio::MCParticleCollection>("MCParticles");

    writer.registerForWrite("EvtHeaders");
    writer.registerForWrite("MCParticles");

    const unsigned int nEvt = 100;

    for ( unsigned int i = 0; i < nEvt; ++i ) {
        if ( i%10 == 0 ) {
            std::cout << "processing event " << i << std::endl;
        }

        auto header = plcio::EventHeader(i, -99, 9999, "SimDet");
        ehc.push_back(header);

        //////////
        writer.writeEvent();
        store.clearCollections();
    }

    writer.finish();

    return 0;
}
