#include "EventHeaderCollection.h"
#include "MCParticleCollection.h"

#include "podio/EventStore.h"
#include "podio/ROOTReader.h"

int main(int argc, char* argv[])
{
    auto reader = podio::ROOTReader();
    auto store = podio::EventStore();

    reader.openFile("test.root");
    store.setReader(&reader);

    unsigned int nEvt = reader.getEntries();

    for (unsigned int i = 0; i < nEvt; ++i ) {
        //if ( i%10 == 0 ) {
        //  std::cout << "processing event " << i << std::endl;
        //}

        auto& hds = store.get<plcio::EventHeaderCollection>("EvtHeaders");
        if ( hds.isValid() ) {
            auto header = hds[0];
            std::cout << "processing event " << header.getEventNumber() << std::endl;
        }
        else {
            return -1;
        }

        auto& mcc = store.get<plcio::MCParticleCollection>("MCParticles");
        if ( ! mcc.isValid() ) {
            return -1;
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

        store.clear();
        reader.endOfEvent();
    }

    reader.closeFile();

    return 0;
}
