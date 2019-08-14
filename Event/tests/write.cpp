#include "EventHeaderCollection.h"
#include "MCParticleCollection.h"

#include "podio/EventStore.h"
#include "podio/ROOTWriter.h"

#include <vector>
#include <iostream>

int main(int argc, char* argv[])
{
    //std::cout << "start writing test" << std::endl;

    auto store = podio::EventStore();
    auto writer = podio::ROOTWriter("test.root", &store);

    auto& ehc = store.create<plcio::EventHeaderCollection>("EvtHeaders");
    auto& mcc = store.create<plcio::MCParticleCollection>("MCParticles");

    writer.registerForWrite("EvtHeaders");
    writer.registerForWrite("MCParticles");

    const unsigned int nEvt = 100;

    for ( unsigned int i = 0; i < nEvt; ++i ) {
        //if ( i%10 == 0 ) {
            std::cout << "processing event " << i << std::endl;
        //}

        auto header = plcio::EventHeader(i, -99, 9999, "SimDet");
        ehc.push_back(header);

        for ( unsigned int i = 0; i < 10; ++i ) {
            mcc.create();
        }

        for ( unsigned int i = 0; i < 4; ++i ) {
            auto p = mcc[i];
            for ( unsigned int j = 0; j < 4; ++j ) {
                unsigned int idx = (2+j) + (i/2)*4;
                p.addDaughter( mcc[idx] );
                mcc[idx].addParent( p );
            }
        }

        for ( auto p : mcc ) {
            std::cout << " particle " << p.getObjectID().index << " has daughters: ";
            for ( auto it = p.daughters_begin(), end = p.daughters_end(); it != end; ++it ) {
                std::cout << " " << it->getObjectID().index;
            }
            std::cout << "  and parents: ";
            for ( auto it = p.parents_begin(), end = p.parents_end(); it != end ; ++it ){
                std::cout << " " << it->getObjectID().index;
            }
            std::cout << std::endl;
        }

        //////////
        writer.writeEvent();
        store.clearCollections();
    }

    writer.finish();

    return 0;
}
