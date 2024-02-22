// read csv
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>

// dependence
#include "TBDataReader.h"
#include "edm4hep/RawTimeSeriesCollection.h"

// root
#include <TFile.h> 
#include <TTree.h>

DECLARE_COMPONENT(TBDataReader)

TBDataReader::TBDataReader(const std::string& name, ISvcLocator* svcLoc)
    : GaudiAlgorithm(name, svcLoc)
{
    declareProperty("RawTimeSeriesOut", m_hitCol);
}

StatusCode TBDataReader::initialize()
{
    info() << "begin initialize TBDataReader" << endmsg;

    // get file path
    tree_file = new TFile("/afs/ihep.ac.cn/users/l/lishuqi/public/TB/AllSameTimeHitsLoose2_9.root", "READ");

    if (!tree_file->IsOpen()) {
        info() << "failed open file" << endmsg;
    }
    tree = (TTree*)tree_file->Get("TimeInfo");

    tree->SetBranchAddress("timeChip", &timeChip);
    tree->SetBranchAddress("planeID", &planeID);
    tree->SetBranchAddress("row", &row);
    tree->SetBranchAddress("col", &col);

    info() << "Total Number of Events = " << tree->GetEntries() << endmsg;
    return GaudiAlgorithm::initialize();
}

StatusCode TBDataReader::execute()
{
    info() << "Executing Event " << eventid << endmsg;

    tree->GetEntry(eventid);
    auto hits = m_hitCol.createAndPut();
    hits->setID(eventid);
    int event_size = timeChip->size();
    for (int i=0; i < event_size; i++){
        auto hit = hits->create();
        hit.setTime(timeChip->at(i));
        // planeID:8,row:12,col:12
        dd4hep::rec::CellID ncellid;
        m_decoder.set(ncellid, "planeID", planeID->at(i));
        m_decoder.set(ncellid, "row", row->at(i));
        m_decoder.set(ncellid, "col", col->at(i));
        hit.setCellID(ncellid);
    }

    eventid += 1;
    return StatusCode::SUCCESS;
}

StatusCode TBDataReader::finalize()
{
    info() << "finalize TBDataReader" << endmsg;
    return GaudiAlgorithm::finalize();
}
