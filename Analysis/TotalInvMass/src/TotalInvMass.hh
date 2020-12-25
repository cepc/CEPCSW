#ifndef _TotalInvMass_hh_
#define _TotalInvMass_hh_

#include "GaudiAlg/GaudiAlgorithm.h"
#include <Gaudi/Property.h>

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
#include "edm4hep/ReconstructedParticleCollection.h"

class TotalInvMass  : public GaudiAlgorithm
{
public:

    TotalInvMass(const std::string& name, ISvcLocator* svcLoc);

    ~TotalInvMass() {};

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

    typedef DataHandle<edm4hep::ReconstructedParticleCollection> RecParticleColHandler;
    RecParticleColHandler m_reconep{"AncientPFOs", Gaudi::DataHandle::Reader, this};
    RecParticleColHandler m_arbopfo{"ArborLICHPFOs", Gaudi::DataHandle::Reader, this};

    Gaudi::Property<std::string> _treeFileName{this,
            "TreeOutputFile", "MCTruth.root",
            "The name of the file to which the ROOT tree will be written"};
    Gaudi::Property<std::string> _treeName{this,
            "TreeName", "MCPart",
            "The name of the ROOT tree"};
    std::string _colName;
    std::string _colAdcVals;

    Gaudi::Property<int> _overwrite{this,
            "OverwriteFile", 0,
            "If zero an already existing file will not be overwritten."};
    int _leptonID;
    float _cmsE; 
    TTree *_outputTree, *_outputPFO;
    float _ISRP[3];
    float _ISREn, _ISRPt; 
    float _N3En, _N3Pt; 
    float _NEn, _NPt; 

    int _HDPID, _OriQuarkID; 
    float _OQDir, _HDir;
    int _NMuP, _NMuM, _NChP, _NChM;
    float _P_MuP[4], _P_MuM[4], _P_DL[4];
    int _EventType; 
    float _InvMass, _RecoilMass; 
    float _J1CosTheta, _J2CosTheta, _MaxJetCosTheta; 
    float _OriJ1CosTheta, _OriJ2CosTheta, _MaxOriJetCosTheta; 
    float _Mass_a, _Mass_p; 

    int _PID1, _PID2;
    float _PL1[4], _PL2[4], _RPL1[4], _RPL2[4], _SM[4], _P_allCharged[4], _P_allNeutral[4], _P_Higgs[4], _P_allReco[4];
    float _Hmass; 
    int _Num;
    int _NHDaug; 
    int _HdaughterPID; 
    int _ZdaughterPID; 
    float _Pz[4], _Ph[4], _PzD1[4], _PzD2[4], _PhD1[4], _PhD2[4], _RPzD1[4], _RPzD2[4], _RPhD1[4], _RPhD2[4];
    float _P[4], _SumP[4], _VisP[4], _MissP[4];
    int _PID, _NFMCP, _MotherFlag, _NNeutrino; 
    float _ENeutrino, _DiPhMass, _DiPhMassCorr; 
    //		float _CosTheta, _Phi, _Charge;
    float _Mz, _Mrecoil, _MzReco, _MhReco, _MrecoilReco; 
    float KthEn[7][9];
    unsigned int _eventNr;

    float _Mass_p_Pisr, _Mass_a_Pisr, _Mass_a_Plcal;

    float TotalP_a[4], TotalP_p[4];
    int nCHPFO_a, nCHPFO_p, nNEPFO_a, nNEPFO_p;
    float NeCaloE_a[2], NeCaloE_p[2];
    float ElargeP[2], EequP[2], EsmallP[2];

    float ChP[4], FrP[4], PhP[4], NeP[4], UdP[4], FrPh[4], FrNe[4], KPF[4];
    float _EcalTotalE, _HcalTotalE, _EcalCluE, _HcalCluE, _EcalCluE_p, _HcalCluE_p; 
    float _HcalEn1, _HcalEn2,_HcalEn3,_HcalEn4,_HcalEn5;
    float _EcalEn1, _EcalEn2,_EcalEn3,_EcalEn4,_EcalEn5;

    int Type, Charge;
    float Energy, TrkSumEn; 
    float P[3], CluEnCom[2];
    float CluEn;

    int TrackHit;
    float StartPos[3], EndPos[3];
    float _visE;

    std::string _fileName;
    std::ostream *_output;
    std::string _histFileName;
};

#endif


