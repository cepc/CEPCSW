#ifndef TBDATAREADER_H
#define TBDATAREADER_H

#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "DD4hep/Detector.h"
#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h> 

namespace edm4hep {
    class RawTimeSeriesCollection;
}

namespace dd4hep {
    namespace rec {
        class CellIDPositionConverter;
    }
    namespace DDSegmentation {
        class BitFieldCoder;
    }
}

class TBDataReader : public GaudiAlgorithm
{

    public :

        TBDataReader(const std::string& name, ISvcLocator* svcLoc);

        virtual StatusCode initialize();
        virtual StatusCode execute();
        virtual StatusCode finalize();

    private :

        DataHandle<edm4hep::RawTimeSeriesCollection> m_hitCol{"RawTimeSeries", Gaudi::DataHandle::Writer, this};
        dd4hep::DDSegmentation::BitFieldCoder m_decoder{"planeID:8,row:12,col:12"};
        // dd4hep::DDSegmentation::BitFieldCoder* m_decoder;

        int eventid = 0;
        TFile *tree_file = NULL;
        TTree *tree = NULL;
        std::vector<int>* timeChip = NULL;
        std::vector<int>* planeID = NULL;
        std::vector<int>* row = NULL;
        std::vector<int>* col = NULL;
};

#endif  // TBDATAREADER_H
