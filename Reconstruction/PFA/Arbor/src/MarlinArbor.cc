#include "MarlinArbor.hh"
#include "ArborTool.h"
#include "ArborToolLCIO.hh"
#include "ArborHit.h"

#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "Gaudi/Property.h"
#include "edm4hep/EventHeader.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/SimCalorimeterHitConst.h"
#include "edm4hep/SimCalorimeterHit.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/Cluster.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/MCRecoCaloAssociationCollection.h"
#include "edm4hep/MCParticleCollection.h"

#include "cellIDDecoder.h"

#include "DD4hep/Detector.h"
#include "DD4hep/IDDescriptor.h"
#include "DD4hep/Plugins.h"

#include <values.h>
#include <string>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <Rtypes.h>
#include <sstream>
#include <set>
#include <TVector3.h>
#include <vector>
#include <algorithm>

#include "DetectorPos.hh"

using namespace std;

extern linkcoll InitLinks; 
extern linkcoll IterLinks_1; 
extern linkcoll IterLinks; 
extern linkcoll links_debug; 
extern branchcoll Trees; 
extern std::vector<int> IsoHitsIndex;

//std::vector<std::string> CaloHitCollections;

DECLARE_COMPONENT(MarlinArbor)

MarlinArbor::MarlinArbor(const std::string& name, ISvcLocator* svcLoc)
     : GaudiAlgorithm(name, svcLoc), 
          _eventNr(0),_output(0)
{
}

StatusCode MarlinArbor::initialize() {

	
	_cepc_thresholds.push_back(10);
	_cepc_thresholds.push_back(90);
	_cepc_thresholds.push_back(50);
	_cepc_thresholds.push_back(7.5);

     m_geosvc = service<IGeomSvc>("GeomSvc");

      ISvcLocator* svcloc = serviceLocator();
      m_ArborToolLCIO=new ArborToolLCIO("arborTools",svcloc);
         for(unsigned int i = 0; i < m_ecalReadoutNames.value().size(); i++){
             m_col_readout_map[m_ecalColNames.value().at(i)] = m_ecalReadoutNames.value().at(i);
         }
         for(unsigned int i = 0; i < m_hcalReadoutNames.value().size(); i++){
             m_col_readout_map[m_hcalColNames.value().at(i)] = m_hcalReadoutNames.value().at(i);
         }

     for (auto& ecal : m_ecalColNames) {
	  _ecalCollections.push_back( new CaloType(ecal, Gaudi::DataHandle::Reader, this) );
	  _calCollections.push_back( new CaloType(ecal, Gaudi::DataHandle::Reader, this) );
     }
     for (auto& hcal : m_hcalColNames) {
	  _hcalCollections.push_back( new CaloType(hcal, Gaudi::DataHandle::Reader, this) );
	  _calCollections.push_back( new CaloType(hcal, Gaudi::DataHandle::Reader, this) );
     }
	return GaudiAlgorithm::initialize();
}

void MarlinArbor::HitsPreparation()
{
	cout<<"Start to prepare Hits"<<endl;
}

void MarlinArbor::MakeIsoHits( std::vector<edm4hep::ConstCalorimeterHit> inputCaloHits, DataHandle<edm4hep::CalorimeterHitCollection>& m_hitcol)
{
	edm4hep::CalorimeterHitCollection* isohitcoll = m_hitcol.createAndPut();

	int nhit = inputCaloHits.size();

	for(int i = 0; i < nhit; i++)
	{
		auto a_hit = inputCaloHits[i];
		auto IsoHit = isohitcoll->create();
		IsoHit.setPosition(a_hit.getPosition());
		IsoHit.setCellID(a_hit.getCellID());
		IsoHit.setEnergy(a_hit.getEnergy());
		//isohitcoll->addElement(collhit);
	}

}

StatusCode MarlinArbor::execute()
{
     //if(_eventNr % m_reportEvery == 0) cout<<"eventNr: "<<_eventNr<<endl;
     _eventNr++;

	MarlinArbor::HitsPreparation();	//Absorb isolated hits; 

	TVector3 currHitPos;

	std::vector< TVector3 > inputHitsPos;
	std::vector< ArborHit > inputABHit; 
	std::vector< edm4hep::ConstCalorimeterHit > inputHits;  
	std::vector< edm4hep::ConstCalorimeterHit > inputECALHits;  
	std::vector< edm4hep::ConstCalorimeterHit > inputHCALHits;  
	std::vector< std::vector<int> > Sequence; 
	int LayerNum = 0; 
	int StaveNum = 0; 
	int SubDId = -10; 
	float Depth = 0; 
	int KShift = 0; 
	TVector3 TrkEndPointPos; 
	std::vector<edm4hep::ConstCalorimeterHit> IsoHits;

	for(unsigned int i1 = 0; i1 < _calCollections.size(); i1++)
	{

		std::cout<<i1<<"th collection"<<m_col_readout_map[m_ecalColNames.value().at(i1)]<<std::endl;
		std::string tmp_readout;
		
	      if(i1<2)tmp_readout = m_col_readout_map[m_ecalColNames.value().at(i1)];
	      else
		      tmp_readout = m_col_readout_map[m_hcalColNames.value().at(i1-2)];

	      std::cout<<tmp_readout<<std::endl;
              // get the DD4hep readout
              m_decoder = m_geosvc->getDecoder(tmp_readout);
			KShift = 0;
			SubDId = -1; 

			if( i1 < _EcalCalCollections.size() )
				SubDId = 1; 
			else if( i1 < _EcalCalCollections.size() + _HcalCalCollections.size() )
				SubDId = 2;
			else
				SubDId = 3; 

			if(i1 >  _EcalCalCollections.size() - 1)
				KShift = 100; 
			else if( i1 == _calCollections.size() - 2)	//HCAL Ring
				KShift = 50;

			auto CaloHitColl = _calCollections[i1]->get();
		
			//int NHitsCurrCol = CaloHitColl->getNumberOfElements();
			//CellIDDecoder<CalorimeterHit> idDecoder(CaloHitColl);
			for (auto a_hit: *CaloHitColl){
		       		currHitPos =  TVector3(a_hit.getPosition().x, a_hit.getPosition().y, a_hit.getPosition().z);
		       		Depth = DisSeedSurface(currHitPos);

				auto cellid = a_hit.getCellID();
			       LayerNum = m_decoder->get(cellid, "layer")+ KShift;
			       StaveNum=m_decoder->get(cellid, "stave");
		       		
		       		if(SubDId!=2 ){

		       			inputECALHits.push_back(a_hit);
		       		}
		       		else{
		       			inputHCALHits.push_back(a_hit);
		       		}
		       		ArborHit a_abhit(currHitPos, LayerNum, 0, Depth, StaveNum, SubDId);
		       		inputABHit.push_back(a_abhit);
		       		inputHits.push_back(a_hit);
		       }
		

			// cout<<i1<<"  Stat  "<<SubDId<<" ~~~ "<<inputABHit.size()<<endl; 

	}
	//cout<<"hit size"<<inputHits.size()<<endl;

	Sequence = Arbor(inputABHit, _cepc_thresholds);   

	m_ArborToolLCIO->ClusterBuilding( branchCol, inputHits, Trees, 0 );

	for(unsigned int i2 = 0; i2 < IsoHitsIndex.size(); i2++)
	{
		auto a_Isohit = inputHits[ IsoHitsIndex[i2] ];
		if(a_Isohit.getEnergy() > 0)	//Veto Trk End Hits
		{
			IsoHits.push_back(a_Isohit);
		}
	}

	MakeIsoHits(IsoHits, m_isohitcol);
	return StatusCode::SUCCESS;
}


StatusCode MarlinArbor::finalize()
{
	std::cout<<"Arbor Ends. Good luck"<<std::endl;
	return GaudiAlgorithm::finalize();
}
