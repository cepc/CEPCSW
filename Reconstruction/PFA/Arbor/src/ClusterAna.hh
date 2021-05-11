#ifndef _ClusterAna_hh_
#define _ClusterAna_hh_

#include "GaudiAlg/GaudiAlgorithm.h"
#include "Gaudi/Property.h"

#include "k4FWCore/DataHandle.h"

#include <string>
#include <iostream>
#include <fstream>
#include <TNtuple.h>
#include <TObject.h>
#include <TTree.h>
#include <TFile.h>

#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"

class ClusterAna  : public GaudiAlgorithm
{
public:

    ClusterAna(const std::string& name, ISvcLocator* svcLoc);

    ~ClusterAna() {};

    StatusCode initialize() override;

    StatusCode execute() override;

    StatusCode finalize() override;

protected:
    typedef DataHandle<edm4hep::MCParticleCollection> MCParticleColHandler;
    MCParticleColHandler m_mcParticle{"MCParticle", Gaudi::DataHandle::Reader, this};

    typedef DataHandle<edm4hep::CalorimeterHitCollection> CaloHitColHandler;
    CaloHitColHandler m_ecalbarrelhitcol{"ECALBarrel", Gaudi::DataHandle::Reader, this};
    CaloHitColHandler m_ecalendcaphitcol{"ECALEndcap", Gaudi::DataHandle::Reader, this};

    CaloHitColHandler m_hcalbarrelhitcol{"HCALBarrel", Gaudi::DataHandle::Reader, this};
    CaloHitColHandler m_hcalendcaphitcol{"HCALEndcap", Gaudi::DataHandle::Reader, this};
    CaloHitColHandler m_hcalotherhitcol {"HCALOther", Gaudi::DataHandle::Reader, this};

    typedef DataHandle<edm4hep::ClusterCollection> CluColHandler;
    CluColHandler m_Clu{"EHBushes", Gaudi::DataHandle::Reader, this};

    typedef DataHandle<edm4hep::ReconstructedParticleCollection> PFOColHandler;
    PFOColHandler m_PFO{"ArborPFO", Gaudi::DataHandle::Reader, this};

Gaudi::Property<std::string> _treeFileName{this,
            "TreeOutputFile", "Ana.root",
            "The name of the file to which the ROOT tree will be written"};
    Gaudi::Property<std::string> _treeName{this,
            "TreeName", "Evt",
            "The name of the ROOT tree"};
    std::string _colName;
    std::string _colAdcVals;

    Gaudi::Property<int> _overwrite{this,
            "OverwriteFile", 0,
            "If zero an already existing file will not be overwritten."};
    TTree *_outputTree, *_outputClu, *_outputPFO;
    
    float _EcalTotalE,_HcalTotalE,_CluE,_CluEn,_PFOEn;
    int _eventNr,_Num,_NClu,_Nhit,_NPFO,_PID;
    float CluPos[3];

    std::string _fileName;
    std::ostream *_output;
    std::string _histFileName;
};

#endif


