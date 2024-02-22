#ifndef TBDATACLUSTERREADER_H
#define TBDATACLUSTERREADER_H

#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "DD4hep/Detector.h"
#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h> 

namespace edm4hep {
    class TrackerHitCollection;
}

namespace dd4hep {
    namespace rec {
        class CellIDPositionConverter;
    }
    namespace DDSegmentation {
        class BitFieldCoder;
    }
}

class TBDataClusterReader : public GaudiAlgorithm
{

    public :

        TBDataClusterReader(const std::string& name, ISvcLocator* svcLoc);

        virtual StatusCode initialize();
        virtual StatusCode execute();
        virtual StatusCode finalize();

    private :

        DataHandle<edm4hep::TrackerHitCollection> m_hitCol{"TrackerHit", Gaudi::DataHandle::Writer, this};
        dd4hep::DDSegmentation::BitFieldCoder m_decoder{"planeID:8,hitsNum:8,hitsNumX:8,hitsNumY:8"};
        // dd4hep::DDSegmentation::BitFieldCoder* m_decoder;

        int eventid = 0;
        TFile *tree_file = NULL;
        TTree *tree = NULL;
        std::vector<int>* planeID = NULL;
        std::vector<int>* hitsNum = NULL;
        std::vector<int>* hitsNumX = NULL;
        std::vector<int>* hitsNumY = NULL;
        std::vector<float>* row = NULL;
        std::vector<float>* col = NULL;
};

#endif  // TBDATACLUSTERREADER_H
