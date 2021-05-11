#ifndef SimHitMergeAlg_H
#define SimHitMergeAlg_H

#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "Gaudi/Property.h"
#include "edm4hep/EventHeader.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/SimCalorimeterHitConst.h"
#include "edm4hep/SimCalorimeterHit.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/MCRecoCaloAssociationCollection.h"
#include "edm4hep/MCParticleCollection.h"

#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h>
#include "DetInterface/IGeomSvc.h"

#include <string>
#include <iostream>
#include <fstream>


class SimHitMergeAlg  : public GaudiAlgorithm
{
public:

     SimHitMergeAlg(const std::string& name, ISvcLocator* svcLoc);

     virtual StatusCode initialize() ;

     virtual StatusCode execute() ;

     virtual StatusCode finalize() ;

protected:
     std::string GetLayerCoding(const std::string &encodingString) const;



     typedef DataHandle<edm4hep::SimCalorimeterHitCollection>  SimCaloType;

     Gaudi::Property<std::vector<std::string>> m_inputColNames{this, "InputCollections", {"EcalBarrelSiliconCollection", "EcalEndcapSiliconCollection", "EcalEndcapRingCollection"}, "Input Hit collection Names"};
     Gaudi::Property<std::vector<std::string>> m_outputColNames{this, "OutputCollections", {"ECALBarrel", "ECALEndcap", "ECALOther"}, "Name of merged Hit Collections"};
     std::vector<SimCaloType*> m_InputCollections;
     std::vector<SimCaloType*> m_OutputCollections;
     Gaudi::Property< bool > m_sanity_check{this, "sanity_check", false, "sanity check"};
     Gaudi::Property< float > m_cell_x{this, "cell_x", 10, ""};//mm
     Gaudi::Property< float > m_cell_y{this, "cell_y", 10, ""};//mm
     Gaudi::Property< float > m_cell_z{this, "cell_z", 10, ""};//mm

     SmartIF<IGeomSvc> m_geosvc;
     dd4hep::Detector* m_dd4hep_geo;
     dd4hep::rec::CellIDPositionConverter* m_cellIDConverter;
};

#endif


