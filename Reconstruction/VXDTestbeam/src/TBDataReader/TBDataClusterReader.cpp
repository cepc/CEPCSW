// read csv
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>

// dependence
#include "TBDataClusterReader.h"
#include "edm4hep/TrackerHitCollection.h"

// root
#include <TFile.h> 
#include <TTree.h>

DECLARE_COMPONENT(TBDataClusterReader)

TBDataClusterReader::TBDataClusterReader(const std::string& name, ISvcLocator* svcLoc)
    : GaudiAlgorithm(name, svcLoc)
{
    declareProperty("TrackerHitOut", m_hitCol);
}

StatusCode TBDataClusterReader::initialize()
{
    info() << "begin initialize TBDataClusterReader" << endmsg;

    // get file path
    tree_file = new TFile("/afs/ihep.ac.cn/users/l/lishuqi/scratchfs_2Weeks/pub/Clusterxy_all2_9.root", "READ");

    if (!tree_file->IsOpen()) {
        info() << "failed open file" << endmsg;
    }
    tree = (TTree*)tree_file->Get("Event");

    tree->SetBranchAddress("planeID", &planeID);
    tree->SetBranchAddress("col", &col);
    tree->SetBranchAddress("row", &row);
    tree->SetBranchAddress("hitsNum", &hitsNum);
    tree->SetBranchAddress("hitsNumX", &hitsNumX);
    tree->SetBranchAddress("hitsNumY", &hitsNumY);

    info() << "Total Number of Events = " << tree->GetEntries() << endmsg;
    return GaudiAlgorithm::initialize();
}

StatusCode TBDataClusterReader::execute()
{
    info() << "Executing Event " << eventid << endmsg;

    tree->GetEntry(eventid);
    auto hits = m_hitCol.createAndPut();
    hits->setID(eventid);
    int event_size = planeID->size();
    for (int i=0; i < event_size; i++){
        auto hit = hits->create();
        // use float type storage
        hit.setTime(col->at(i));
        hit.setEDep(row->at(i));
        // planeID:8,hitsNum:8,hitsNumX:8,hitsNumY:8
        dd4hep::rec::CellID ncellid;
        m_decoder.set(ncellid, "planeID", planeID->at(i));
        m_decoder.set(ncellid, "hitsNum", hitsNum->at(i));
        m_decoder.set(ncellid, "hitsNumX", hitsNumX->at(i));
        m_decoder.set(ncellid, "hitsNumY", hitsNumY->at(i));
        hit.setCellID(ncellid);
    }

    eventid += 1;
    return StatusCode::SUCCESS;
}

StatusCode TBDataClusterReader::finalize()
{
    info() << "finalize TBDataClusterReader" << endmsg;
    return GaudiAlgorithm::finalize();
}
